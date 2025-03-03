/*!

  @mainpage build_signature_nudb

  # Building Signature KMers

  The build_signature_nudb tool is used to construct the set of signature
  kmers from a set of fasta files of protein sequences and their functional assignments.

  The overall flow of the computation is as follows:

  -# Read command line parameters.  Use `build_sigature_nudb --help` for a summary of the
  supported parameters.
  -# Create a FunctionMap instance.
  -# Initialize `deleted_fids` with the contents of the parameter `=-deleted-features-file`
  -# Create output directory
  -# Initialize TBB global limit for parallelism to the value provided via the `--n-threads` parameter
  -# Create a KmerProcessor instance.
  -# Load the FunctionMap with the good roles specified with the parameter `--good-roles`
  -# Load the FunctionMap with the good functions specified with the parameter `--good-functions`
  -# For each of the function definition files found in the directories specfied with the parameters `--definition-dir` (these were loaded using populate_path_list()),
   - Load assignments using FunctionMap::load_id_assignments()
  -# Create the ref tbb::concurrent_vector all_fasta_data to store locations of the sequence data files
  -# For each file fasta data in `fasta_data` (loaded with populate_path_list() from the directories from parameters `--fasta-dir`)
   - Load the fasta into the FunctionMap using FunctionMap::load_fasta_file()
   - Add the file to the concurrent vector `all_fasta_data`
  -# For each file fasta data in `fasta_data_kept_functions` (loaded with populate_path_list() from the directories from parameters `--fasta-keep-functions-dir`)
   - Load the fasta into the FunctionMap using FunctionMap::load_fasta_file() with the `keep_function_flag` set to true.
   - Add the file to the concurrent vector `all_fasta_data`
  -# Use FunctionMap::process_kept_functions() to compute the set of functions to generate signatures kmers for.
  -# Write the function index file to the output directory. This maps the numeric function index as used in
  the internal data structures to the function string.
  -# In parallel, use load_fasta() to load the sequence data. This loads data into the large concurrent multimap
  that maps from kmer to information about all occurrences of this kmer in the sequence data. We log the function
  for each occurrence along with the kmer's offset from the end of the sequence and the size of the sequence.
  -# In parallel, process the kmers in the @ref KmerAttributeMap.
   - We group map entries having the same kmer into an instance of KmerSet.
   - When the set is complete, use process_set() to compute statistics on the set, determine
   if the kmer will be kept as a signature, and if so insert it into the KeptKmers map.
  -# In parallel across the KeptKmers map, compute signature weights.

  At this point the signature kmers are computed. We can calcuate a diagnostic recall of the input proteins, and save
  the output as a compact NuDB database.
   
*/


#include <cmath>
#include <array>
#include <algorithm>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <iostream>
#include <fstream>
#include <functional>
#include <memory>
#include <experimental/optional>
#include <unistd.h>
#include <mutex>

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>
#include <boost/thread/tss.hpp>
#include <boost/thread/thread.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "nudb_kmer_db.h"
#include "seed_utils.h"
#include "operators.h"
#include "fasta_parser.h"
#include "kmer_types_generic.h"
#include "kmer_generic.h"

#define TBB_PREVIEW_NUMA_SUPPORT 1
#define TBB_PREVIEW_CONCURRENT_ORDERED_CONTAINERS 1

#include "tbb/global_control.h"
#include "tbb/info.h"
#include "tbb/parallel_for.h"
#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_map.h"
#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_hash_map.h"
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_queue.h"

#include "welford.h"

using namespace seed_utils;

#define DEFINE_GLOBALS
#include "global.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace acc = boost::accumulators;
using namespace KmerGenericNS;

#include "function_map.h"

const int K = KMER_SIZE;
const int MaxSequencesPerFile = 100000;

// Kmer data type
typedef std::array<char, K> Kmer;

// Hash function on kmers
namespace std {
template<class T, size_t N>
struct hash<std::array<T, N>> {
    auto operator() (const array<T, N>& key) const {
	hash<T> hasher;
	size_t result = 0;
	for(size_t i = 0; i < N; ++i) {
	    result = result * 31 + hasher(key[i]); // ??
	}
	return result;
    }
};
}

struct KmerAttributes
{
    FunctionIndex func_index;
    OTUIndex otu_index;
    unsigned short offset;
    unsigned int seq_id;
    unsigned int protein_length;
};

namespace tbb {
    // Test whether tbb_hash_compare can be partially specialized as stated in Reference manual.
    template<> struct tbb_hash_compare<Kmer> {
	size_t hash( Kmer ) const {return 0;}
	bool equal( Kmer /*x*/, Kmer /*y*/ ) {return true;}
    };
}

struct tbb_hash {
    tbb_hash() {}
    size_t operator()(const Kmer& k) const {
        size_t h = 0;
	for (auto s: k)
	    h = (h * 17) ^ (unsigned int) s;
	return h;
    }
};

template<class T, size_t N>
struct MyHashCompare {
    static size_t hash( const std::array<T,N> &x ) {
        size_t h = 0;
	for (auto s: x)
	    h = (h*17)^s;
	return h;
    }
    //! True if strings are equal
    static bool equal( const std::array<T,N>&  x, const std::array<T,N> & y ) {
	return x==y;
    }
};

//typedef std::unordered_multimap<Kmer, KmerAttributes> KmerAttributeMap;
typedef tbb::concurrent_unordered_multimap<Kmer, KmerAttributes, tbb_hash> KmerAttributeMap;

inline std::ostream &operator<<(std::ostream &os, const Kmer &k)
{
    std::copy(k.begin(), k.end(),
	      std::ostream_iterator<char>(os, ""));
    return os;
}

struct KmerStatistics
{
    tbb::atomic<int> distinct_signatures = 0;
    tbb::concurrent_unordered_map<int, int> distinct_functions;
    tbb::concurrent_unordered_map<FunctionIndex, int> seqs_with_func;
    tbb::concurrent_unordered_set<unsigned int> seqs_with_a_signature;
};

//
// Whack some globals in here. We're maintaining enough state this
// needs to be wrapped up in an object.
//

KmerStatistics kmer_stats;

std::mutex io_mutex;
std::ofstream rejected_stream;
std::ofstream kept_stream;
std::ofstream kept_function_stream;

/*
 * Object to keep state for kmers we are keeping. We hang onto
 * some statistics in order to compute weights later on.
 */
struct KeptKmer
{
    Kmer kmer;
    unsigned short avg_from_end;		// Median offset from the end of the protein
    FunctionIndex function_index;
//    FunctionIndex function_index2;	// Secondary function in the case of ambiguity
    OTUIndex  otu_index;
    unsigned int seqs_containing_sig;	// Count of sequences containing this kmer
    unsigned int seqs_containing_function; // Count of sequences with the kmer that have the function
    float weight;
    unsigned short mean;
    unsigned short median;
    unsigned short var;
};

inline std::ostream &operator<<(std::ostream &os, const KeptKmer &k)
{
    os << k.kmer << " " << k.function_index << " " <<  k.seqs_containing_sig << " " << k.seqs_containing_function;
    return os;
}


// tbb::concurrent_vector<KeptKmer> kept_kmers;
typedef tbb::concurrent_unordered_map<Kmer, KeptKmer, tbb_hash> KeptKmers;
//tbb::concurrent_unordered_map<Kmer, KeptKmer, tbb_hash> kept_kmers;
KeptKmers kept_kmers;

/*
 * This is a little wrapper so we can use the KeptKmers map as
 * a database for kmer calling.
 */

class KeptKmerDB
{
public:
    static const int KmerSize = K;
    using KData = KeptKmer;

    KeptKmerDB(KeptKmers &kk) :
	kept_kmers_(kk) {
    }

    template <typename CB>
    void fetch(const Kmer &k, CB cb, int &ec) {
	auto iter = kept_kmers_.find(k);
	if (iter != kept_kmers_.end())
	{
	    cb(&(iter->second));
	}
	ec = 0;
    };

private:
    KeptKmers &kept_kmers_;
};

struct KmerSet
{
    KmerSet() : count(0) {}
    void reset() {
	count = 0;
	func_count.clear();
	set.clear();
    }
    // ~KmerSet() { std::cerr << "destroy " << this << "\n"; }
    Kmer kmer;
    std::map<FunctionIndex, int> func_count;
    int count;
    std::vector<KmerAttributes> set;
};

inline std::ostream &operator<<(std::ostream &os, const KmerSet &ks)
{
    os << "KmerSet: kmer=" << ks.kmer << " count=" << ks.count << "\n";
    for (auto iter = ks.func_count.begin(); iter != ks.func_count.end(); iter++)
    {
	os << iter->first << "\t" << iter->second << "\n";
    }
    return os;
}

void process_set(KmerSet &set);
typedef std::vector<KmerSet> KmerSetList;
typedef std::shared_ptr<KmerSetList> KmerSetListPtr;

struct KmerProcessor
{
    KmerProcessor(int n) : n_threads(n) {}
    boost::thread_group thread_pool;
    tbb::concurrent_bounded_queue<KmerSetListPtr> queue;
    int n_threads;
    void start() {
	for (int i = 0; i < n_threads; i++)
	{
	    std::cerr << "starting " << i << "\n";
	    thread_pool.create_thread([this, i]() {
		    thread_main(i);
		});
	}
    }
    void stop() {
	for (int i = 0; i < n_threads; i++)
	{
	    std::cerr << "stopping " << i << "\n";
	    KmerSetListPtr work = std::make_shared<KmerSetList>();
	    enqueue_work(work);
	}
	std::cerr << "Awaiting threads\n";
	thread_pool.join_all();
    }
	
    void enqueue_work(KmerSetListPtr work) {
	queue.push(work);
    }
    void thread_main(int i) {
	std::cerr << "running " << i << "\n";
	int sp = 0;
	while (1)
	{
	    KmerSetListPtr work;
	    queue.pop(work);
	    if (work->size() == 0)
	    {
		std::cerr << "shutting down " << i << "\n";
		break;
	    }

	    for (auto entry: *work)
	    {
		// std::cerr << i << " process " << entry << "\n";
		sp++;
		process_set(entry);
	    }
	}
	std::cerr << "thread " << i << " processed " << sp << " sets\n";
    }
};
KmerProcessor *g_kmer_processor;

std::set<unsigned char> ok_prot = { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
			   'a', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n', 'p', 'q', 'r', 's', 't', 'v', 'w', 'y'};

/*!
  @brief Load a single sequence.

  Use the provided FunctionMap to look up the function for the given sequence. If it is not present,
  skip this sequence.

  Update KmerStatistics::seqs_with_func to reflect the additional sequence having this function.

  For each kmer in the sequence,

  - If kmer has no invalid characters, insert it into the @ref KmerAttributeMap. This logs the existence
  of the kmer with the given function, offset from the end of its protein, and length of the source protein.

  
*/

void load_sequence(FunctionMap &fm, KmerAttributeMap &m, unsigned int &next_sequence_id,
		   const std::string &id, const std::string &def, const std::string &seq)
{
    if (id.empty())
	return;

    std::string func = fm.lookup_function(id);
    
    /*
     * Empty means empty (and perhaps deleted feature).
     */
    
    if (func.empty())
    {
	return;
    }
    
    unsigned int seq_id = next_sequence_id++;

    FunctionIndex function_index = fm.lookup_index(func);

    if (false)
    {
	if (function_index == UndefinedFunction)
	{
	    function_index = fm.lookup_index("hypothetical protein");
	    if (function_index == UndefinedFunction)
	    {
		std::cerr << "No function defined for hypothetical protein\n";
		exit(1);
	    }
	}
    }

    // std::cout << "Load seq " << id << " '" << func << "' " << function_index << " id=" << seq_id <<  "\n";
//    if (function_index == 19345)
//    {
//	std::cout << "FN 19345 " << id << " '" << func << "' " << function_index << " id=" << seq_id <<  "\n";
//    }
    // std::cout << "Load seq " << id << " '" << func << "' " << function_index << " id=" << seq_id <<  "\n";
    if (function_index == UndefinedFunction)
    {
//	std::cout << "Skipping undef seq " << id << " '" << func << "' " << function_index << " id=" << seq_id <<  "\n";
    	return;
    }

    kmer_stats.seqs_with_func[function_index]++;

    for (auto it = seq.begin(); it < seq.end() - K + 1; it++)
    {
        unsigned short n = (unsigned short) std::distance(it, seq.end());
	Kmer kmer;
	// Kmer::iterator kiter = kmer.begin();
	bool ok = true;
	std::copy_n(it, K, kmer.begin());
	for (auto x: kmer)
	{
	    if (ok_prot.find(x) == ok_prot.end())
	    {
		// std::cerr << "not ok prot "  << x << "\n";
		ok = false;
		break;
	    }
	}
	if (ok)
	{
	    /*
	    if ("PHSQRWTN" == kmer || function_index == 19345)
	    {
		std::cerr << kmer << " in " << id << " " << function_index << " " << seq_id << "\n";
	    }
	    */
	    m.insert({kmer, { function_index, UndefinedOTU, n, seq_id, static_cast<unsigned int>(seq.length())}});
	    // std::cout << kmer << " " << n << "\n";
	}
	else
	{
	    /*
	    if ("PHSQRWTN" == kmer || function_index == 19345)
	    {
		std::cerr << kmer << " rejected from " << id << " " << function_index << " " << seq_id << "\n";
	    }
	    */
	    // std::cout << "reject " << kmer << " " << n << "\n";
	}
    }
}


/*!
  @brief Load sequence data.

  Parse the fasta data using FastaParser.

  For each sequence, invoke load_sequence() to create and store kmers.
  
  If a function index is not defined, this is not a function we wish
  to process so skip this sequence.
  
 */
void load_fasta(FunctionMap &fm, KmerAttributeMap &m, unsigned file_number, const fs::path &file,
		const std::set<std::string> &deleted_fids)
{
    fs::ifstream ifstr(file);

    FastaParser parser;
    
    unsigned next_sequence_id = file_number * MaxSequencesPerFile;

    parser.set_def_callback([&fm, &m, &next_sequence_id, &deleted_fids](const std::string &id, const std::string &def, const std::string &seq) {
	if (deleted_fids.find(id) == deleted_fids.end())
	{
	    load_sequence(fm, m, next_sequence_id, id, def, seq);
	}
	return 0;
    });
    parser.parse(ifstr);
    parser.parse_complete();
}

void enqueue_set(KmerSetListPtr set_list)
{
    g_kmer_processor->enqueue_work(set_list);
}

//void process_set(Kmer &kmer, std::map<FunctionIndex, int> &func_count, int count, std::vector<KmerAttributes> &set)

/*! @brief Process a set of instances of a given kmer.

 */
void process_set(KmerSet &set)
{
    FunctionIndex best_func_1 = UndefinedFunction, best_func_2 = UndefinedFunction;
    int best_count_1 = -1, best_count_2 = -1;
    
    /*
    auto elt = std::max_element(set.func_count.begin(), set.func_count.end(),
				[](auto a, auto b) { return a.second < b.second; });

    FunctionIndex best_func = elt->first;
    int best_count = elt->second;
    */

    // we want the top two elements by value in the map; don't know how
    // with standard STL without copying to a vector, but if we're copying
    // anyway we can just search for them.
    for (auto x: set.func_count)
    {
	if (best_func_1 == UndefinedFunction)
	{
	    best_func_1 = x.first;
	    best_count_1 = x.second;
	}
	else if (x.second > best_count_1)
	{
	    best_func_2 = best_func_1;
	    best_count_2 = best_count_1;

	    best_func_1 = x.first;
	    best_count_1 = x.second;
	}
	else if (x.second > best_count_2)
	{
	    best_func_2 = x.first;
	    best_count_2 = x.second;
	}
    }

    if ("PHSQRWTN" == set.kmer)
    {
	std::cerr << set;
	std::cerr << "best=" << best_func_1 << " best2=" << best_func_2 << "\n";
    }


    float thresh = float(set.count) * 0.8f;
    int best_count = best_count_1;
    FunctionIndex best_func = best_func_1;

    if (rejected_stream.is_open() || kept_stream.is_open())
    {
	std::lock_guard<std::mutex> lock(io_mutex);

	/*
	kept_stream << "Process set for " << set.kmer << " best1=" << best_func_1 << " best_count_1=" << best_count_1
		    << " best2=" << best_func_2 << " best_count_2=" << best_count_2 << " count=" << set.count << " thresh=" << thresh << "\n";
	*/

/*
	for (auto x: set.func_count)
	{
	    kept_stream << x.first << " " << x.second << "\n";
	}
*/
	   
	/*
	std::sort(set.set.begin(), set.set.end(), [&set](auto a, auto b) {
		return set.func_count[b.func_index] < set.func_count[a.func_index];
	    });

	for (auto x: set.set)
	{
	    kept_stream << "  " << x.func_index << " " << x.seq_id << "\n";
	}
	*/
	
	if ((float) best_count < thresh)
	{
	    // std::cout << "discard for " << best_count << " < " << thresh << "\n";
	    if (rejected_stream.is_open())
		rejected_stream << set.kmer <<  " best_count=" << best_count << " thresh=" << thresh << "\n";
	    if ((float) (best_count_1 + best_count_2) >= thresh)
	    {
		if (kept_stream.is_open())
		    kept_stream << "AMBIG\t" << best_func_1 << "\t"
				<< best_func_2 << "\t"
				<< best_count_1 << "\t"
				<< best_count_2 << "\n";
	    }
	    return;
	}
    }
    else
    {
	if ((float) best_count < thresh)
	{
	    if ("PHSQRWTN" == set.kmer)
	    {
		std::cerr << "rejecting " << set.kmer << " best_count=" << best_count << " thresh=" << thresh << "\n";
	    }
	    return;
	}
    }

    unsigned int seqs_containing_func = 0;
    std::vector<unsigned short> offsets;

    acc::accumulator_set<unsigned short, acc::stats<acc::tag::mean,
						  acc::tag::median,
						  acc::tag::variance> > acc;

    for (auto item: set.set)
    {
	if (item.func_index == best_func)
	{
	    seqs_containing_func++;
	    acc(item.protein_length);
	}
	offsets.push_back(item.offset);
	kmer_stats.seqs_with_a_signature.insert(item.seq_id);
    }

    unsigned short mean = acc::mean(acc);
    unsigned short median = acc::median(acc);
    unsigned short var = acc::variance(acc);

    std::sort(offsets.begin(), offsets.end());
    unsigned short avg_from_end = offsets[offsets.size() / 2];
    // std::cout << seqs_containing_func << " " << avg_from_end<< "\n";

    kmer_stats.distinct_signatures++;
    kmer_stats.distinct_functions[best_func]++;

    kept_kmers.emplace(set.kmer, KeptKmer { set.kmer, avg_from_end, best_func, UndefinedOTU,
	    (unsigned int) set.set.size(), seqs_containing_func, 0.0, mean, median, var });
}


#if 0
void process_kmers(KmerAttributeMap &m)
{
    //std::map<int, int> func_count;
    //std::vector<KmerAttributes> set;
    //int count = 0;
    Kmer cur { 0 } ;
    
    KmerSetListPtr cur_set_list = std::make_shared<KmerSetList>();
    cur_set_list->emplace_back(KmerSet());
//    KmerSet cur_set;

    int sets_pushed = 0;
	
    for (auto ent: m)
    {
	const Kmer &kmer = ent.first;
	KmerAttributes &attr = ent.second;

	if (kmer != cur)
	{
	    if (cur_set_list->back().count > 0)
	    {
		sets_pushed++;
		if (cur_set_list->size() > 200000)
		{
		    /*
		    std::cerr << "push setlist size " << cur_set_list->size() << "\n";
		    for (auto x: *cur_set_list)
		    {
			std::cerr << "  " << x;
		    }
		    */
		    enqueue_set(cur_set_list);
		    cur_set_list = std::make_shared<KmerSetList>();
		}
		cur_set_list->emplace_back(KmerSet());
	    }
	    //cur_set.reset();
	    cur_set_list->back().kmer = kmer;
	    cur = kmer;
	}
	KmerSet &s = cur_set_list->back();
	s.func_count[attr.func_index]++;
	s.count++;
	s.set.emplace_back(attr);
    }
    sets_pushed++;
//    cur_set_list->emplace_back(cur_set);
    enqueue_set(cur_set_list);
    std::cout << "done, pushed " << sets_pushed << " sets\n";
}
#endif

void par_process_kmers(KmerAttributeMap &m)
{
    tbb::parallel_for(m.range(), [](auto r) {
	    KmerSet cur_set;
	    Kmer cur { 0 };
	    for (auto ent = r.begin(); ent != r.end(); ent++)
	    {
		const Kmer &kmer = ent->first;
		KmerAttributes &attr = ent->second;

		if (kmer != cur)
		{
		    if (cur_set.count > 0)
			process_set(cur_set);
		    
		    cur_set.reset();
		    cur_set.kmer = kmer;
		    cur = kmer;
		}
		cur_set.func_count[attr.func_index]++;
		cur_set.count++;
		cur_set.set.emplace_back(attr);
	    }
	    process_set(cur_set);
	});
}

#if 0
// Not used apparently
void process_kmer_block(KmerAttributeMap &m)
{
    Kmer cur { 0 } ;
    
    KmerSetListPtr cur_set_list = std::make_shared<KmerSetList>();
    cur_set_list->emplace_back(KmerSet());
//    KmerSet cur_set;

    int sets_pushed = 0;
	
    for (auto ent: m)
    {
	const Kmer &kmer = ent.first;
	KmerAttributes &attr = ent.second;

	if (kmer != cur)
	{
	    if (cur_set_list->back().count > 0)
	    {
		sets_pushed++;
		if (cur_set_list->size() > 200000)
		{
		    /*
		    std::cerr << "push setlist size " << cur_set_list->size() << "\n";
		    for (auto x: *cur_set_list)
		    {
			std::cerr << "  " << x;
		    }
		    */
		    enqueue_set(cur_set_list);
		    cur_set_list = std::make_shared<KmerSetList>();
		}
		cur_set_list->emplace_back(KmerSet());
	    }
	    //cur_set.reset();
	    cur_set_list->back().kmer = kmer;
	    cur = kmer;
	}
	KmerSet &s = cur_set_list->back();
	s.func_count[attr.func_index]++;
	s.count++;
	s.set.emplace_back(attr);
    }
    sets_pushed++;
//    cur_set_list->emplace_back(cur_set);
    enqueue_set(cur_set_list);
    std::cout << "done, pushed " << sets_pushed << " sets\n";
}
#endif

void compute_weight_of_signature(KeptKmer &kk)
{
    float NSF = (float) kmer_stats.seqs_with_a_signature.size();
    float KS  = (float) kmer_stats.distinct_signatures;
    // float KF  = kmer_stats.seqs_with_func.size();
    float NSi = (float) kk.seqs_containing_sig;
    float NFj = (float) kmer_stats.seqs_with_func[kk.function_index];
    float NSiFj = (float) kk.seqs_containing_function;

    kk.weight = std::log((NSiFj + 1.0f) / (NSi - NSiFj + 1.0f)) +
	std::log((NSF - NFj + KS) / (NFj + KS));

}

void write_function_index(const fs::path &dir, FunctionMap &fm)
{
    fm.write_function_index(dir);
}

void write_nudb_data(const std::string &nudb_file, tbb::concurrent_vector<KeptKmer> &kmers)
{
    typedef NuDBKmerDb<8> KDB;

    KDB db(nudb_file);

    if (!db.exists())
    {
	std::cerr << "creating new db\n";
	db.create();
    }
    db.open();
    
    for (auto k: kmers)
    {
	// std::cout << k << "\n";
	nudb::error_code ec;
//	db.insert(k.kmer, { k.otu_index, k.avg_from_end, k.function_index, k.weight }, ec);
	db.insert(k.kmer, { k.otu_index, k.avg_from_end, k.function_index, k.weight,
		k.mean, k.median, k.var, static_cast<unsigned short>(k.seqs_containing_function) }, ec);
//	if (!ec)
//	    std::cerr << "insert error: " << ec.message() << "\n";
    }
}

void write_nudb_data(const std::string &nudb_file, const KeptKmers &kmers)
{
    typedef NuDBKmerDb<8> KDB;

    KDB db(nudb_file);

    if (!db.exists())
    {
	std::cerr << "creating new db\n";
	db.create();
    }
    db.open();
    
    for (auto ent: kmers)
    {
	auto k = ent.second;
	// std::cout << k << "\n";
	nudb::error_code ec;
//	db.insert(k.kmer, { k.otu_index, k.avg_from_end, k.function_index, k.weight }, ec);
	db.insert(k.kmer, { k.otu_index, k.avg_from_end, k.function_index, k.weight,
		k.mean, k.median, k.var, static_cast<unsigned short>(k.seqs_containing_function) }, ec);
//	if (!ec)
//	    std::cerr << "insert error: " << ec.message() << "\n";
    }
}


void show_ps()
{
    std::string cmd = "ps uwww" + std::to_string(getpid());
    int rc = system(cmd.c_str());
}

/*!
  @brief Populate a list of paths from a list of directories.

  For each path in dirs, add any regular file found in that path to paths.

  @param dirs List of directories to load
  @param paths List of paths to populate
*/
void populate_path_list(const std::vector<std::string> &dirs, std::vector<fs::path> &paths)
{
    for (auto dir: dirs)
    {
	fs::path p(dir);
	for (auto dit: fs::directory_iterator(dir))
	{
	    if (fs::is_regular_file(dit.path()))
	    {
		paths.emplace_back(dit);
	    }
	}
    }
}    

void load_strings(const std::vector<std::string> &files, std::vector<std::string> &strings)
{
    for (auto f: files)
    {
	std::ifstream ifstr(f);
	if (ifstr.good())
	{
	    // std::cout << "load " << f << "\n";
	    std::string line;
	    while (std::getline(ifstr, line, '\n'))
	    {
		strings.emplace_back(line);
	    }
	}
	else
	{
	    std::cerr << "could not open " << f << "\n";
	}
    }
}

bool process_command_line_options(int argc, char *argv[],
				  std::vector<fs::path> &function_definitions,
				  std::vector<fs::path> &fasta_data,
				  std::vector<fs::path> &fasta_data_kept_functions,
				  std::vector<std::string> &good_functions,
				  std::vector<std::string> &good_roles,
				  fs::path &deleted_fids_file,
				  fs::path &ignored_functions_file,
				  int &min_reps_required,
				  fs::path &kmer_data_dir,
				  fs::path &final_kmers,
				  std::string &nudb_file,
				  int &n_threads)
{
    std::ostringstream x;
    x << "Usage: " << argv[0] << " [options]\nAllowed options";
    po::options_description desc(x.str());

    std::vector<std::string> definition_dirs;
    std::vector<std::string> fasta_dirs;
    std::vector<std::string> fasta_keep_dirs;
    std::vector<std::string> good_function_files;
    std::vector<std::string> good_role_files;

    n_threads = 1;

    desc.add_options()
	("definition-dir,D", po::value<std::vector<std::string>>(&definition_dirs)->multitoken(), "Directory of function definition files")
	("fasta-dir,F", po::value<std::vector<std::string>>(&fasta_dirs)->multitoken(), "Directory of fasta files of protein data")
	("fasta-keep-functions-dir,K", po::value<std::vector<std::string>>(&fasta_keep_dirs), "Directory of fasta files of protein data (keep functions defined here)")
	("good-functions", po::value<std::vector<std::string>>(&good_function_files), "File containing list of functions to be kept")
	("good-roles", po::value<std::vector<std::string>>(&good_role_files), "File containing list of roles to be kept")
	("deleted-features-file", po::value<fs::path>(&deleted_fids_file), "File containing list of deleted feature IDs")
	("ignored-functions-file", po::value<fs::path>(&ignored_functions_file), "File containing list of functions for which we do not create signatures")
	("kmer-data-dir", po::value<fs::path>(&kmer_data_dir), "Write kmer data files to this directory")
	("nudb-file", po::value<std::string>(&nudb_file), "Write saved kmers to this NuDB file base. Should be on a SSD drive.")
	("min-reps-required", po::value<int>(&min_reps_required), "Minimum number of genomes a function must be seen in to be considered for kmers")
	("final-kmers", po::value<fs::path>(&final_kmers), "Write final.kmers file to be consistent with km_build_Data")
	("n-threads", po::value<int>(&n_threads), "Number of threads to use")
	("help,h", "show this help message");

    po::variables_map vm;

    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help"))
    {
	std::cout << desc << "\n";
	return false;
    }

    po::notify(vm);

    /*
     * Read definition and fasta dirs to populate the path lists.
     */
    populate_path_list(definition_dirs, function_definitions);
    populate_path_list(fasta_dirs, fasta_data);
    populate_path_list(fasta_keep_dirs, fasta_data_kept_functions);

    std::cout << "definitions: ";
    for (auto x: definition_dirs)
	std::cout << x << " ";
    std::cout << std::endl;
    std::cout << "fasta: ";
    for (auto x: fasta_dirs)
	std::cout << x << " ";
    std::cout << std::endl;
    std::cout << "keep: ";
    for (auto x: fasta_keep_dirs)
	std::cout << x << " ";
    std::cout << std::endl;

    load_strings(good_function_files, good_functions);
    load_strings(good_role_files, good_roles);

    return true;
}

int main(int argc, char *argv[])
{
    KmerAttributeMap m;
    unsigned next_sequence_id = 0;

    // rejected_stream.open("/dev/shm/rejected.by.build_signature_kmers");
    // kept_stream.open("kept.by.build_signature_kmers");
    // kept_function_stream.open("function.kept.log");

    std::vector<fs::path> function_definitions;
    std::vector<fs::path> fasta_data;
    std::vector<fs::path> fasta_data_kept_functions;

    std::vector<std::string> good_functions;
    std::vector<std::string> good_roles;

    fs::path final_kmers;
    fs::path deleted_fids_file;
    fs::path kmer_data_dir;
    fs::path ignored_functions_file;
    
    int min_reps_required = 3;
    
    int n_threads;

    std::string nudb_file;

    if (!process_command_line_options(argc, argv,
				      function_definitions,
				      fasta_data,
				      fasta_data_kept_functions,
				      good_functions,
				      good_roles,
				      deleted_fids_file,
				      ignored_functions_file,
				      min_reps_required,
				      kmer_data_dir,
				      final_kmers,
				      nudb_file,
				      n_threads))
    {
	return 1;
    }

    FunctionMap fm((kmer_data_dir / "kept_functions.log").string());

    std::set<std::string> ignored_functions;
    /*
     * Read ignored functions if present
     */
    if (!ignored_functions_file.empty())
    {
	fs::ifstream ifstr(ignored_functions_file);
	std::string line;
	while (std::getline(ifstr, line, '\n'))
	{
	    ignored_functions.emplace(line);
	}
	
    }
    

    std::set<std::string> deleted_fids;
    /*
     * Read deleted fids if present.
     */
    if (!deleted_fids_file.empty())
    {
	fs::ifstream ifstr(deleted_fids_file);
	std::string line;
	while (std::getline(ifstr, line, '\n'))
	{
	    deleted_fids.emplace(line);
	}
	
    }


    /*
     * validate our kmer_data_dir
     */
    if (!kmer_data_dir.empty())
    {
	if (!fs::is_directory(kmer_data_dir))
	{
	    if (!fs::create_directory(kmer_data_dir))
	    {
		std::cerr << "Error creating " << kmer_data_dir << "\n";
		exit(1);
	    }
	}
    }

    tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, n_threads);

    g_kmer_processor = new KmerProcessor(n_threads);

    fm.add_good_roles(good_roles);
    fm.add_good_functions(good_functions);

    for (auto def: function_definitions)
    {
	fm.load_id_assignments(def);
    }

    tbb::concurrent_vector<fs::path> all_fasta_data;

    std::cerr << "load fasta\n";

    for (auto fasta: fasta_data)
    {
	fm.load_fasta_file(fasta, false, deleted_fids);
	all_fasta_data.emplace_back(fasta);
    }

    for (auto fasta: fasta_data_kept_functions)
    {
	fm.load_fasta_file(fasta, true, deleted_fids);
	all_fasta_data.emplace_back(fasta);
    }

    /*
     * Process the list of functions and
     * manage the set that we want to keep.
     */
    fm.process_kept_functions(min_reps_required, ignored_functions);

    if (!kmer_data_dir.empty())
    {
	write_function_index(kmer_data_dir, fm);
	/* Write an empty otu index and genomes file (genomes file
	 * can't be empty because kmer_search uses "-s genomes" to test
	 * for its existence). */
	{
	    fs::ofstream otu(kmer_data_dir / "otu.index");
	    otu.close();
	    fs::ofstream genomes(kmer_data_dir / "genomes");
	    genomes << "empty genomes\n";
	    genomes.close();
	}
    }
    
    if (0)
    {
	fm.dump();
	exit(0);
    }
    /*
     * With that done, go ahead and extract kmers.
     */

    if (n_threads < 0)
    {
	for (unsigned i = 0; i < (unsigned) all_fasta_data.size(); i++)
	{
	    load_fasta(fm, m, i, all_fasta_data[i], deleted_fids);
	}
    }
    else
    {
	size_t n = all_fasta_data.size();
	tbb::parallel_for(tbb::blocked_range<size_t>(0, n),
			  [&fm, &m, &next_sequence_id, &all_fasta_data, &deleted_fids](const tbb::blocked_range<size_t> &r) {
			      for (size_t i = r.begin(); i != r.end(); ++i)
			      {
				  auto fasta = all_fasta_data[i];
				  // std::cout << "load file " << i << " " << fasta << "\n";
				  load_fasta(fm, m, (unsigned) i, fasta, deleted_fids);
			      }
			  });
    }

    /*
    g_kmer_processor->start();
    process_kmers(m);
    g_kmer_processor->stop();
    */
    std::cerr << "processing kmers\n";
    par_process_kmers(m);
    std::cout << "Kept " << kept_kmers.size() << " kmers\n";
    std::cout << "distinct_signatures=" << kmer_stats.distinct_signatures << "\n";
    std::cout << "num_seqs_with_a_signature=" << kmer_stats.seqs_with_a_signature.size() << "\n";
    std::cerr << "computing weights\n";
//    std::for_each(kept_kmers.begin(), kept_kmers.end(), compute_weight_of_signature);

    tbb::parallel_for(kept_kmers.range(), [](auto r) {
	for (auto ent = r.begin(); ent != r.end(); ent++)
	{
	    compute_weight_of_signature(ent->second);
	}
    });
/*
    size_t nk = kept_kmers.size();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nk),
			  [](const tbb::blocked_range<size_t> &r) {
			      for (size_t i = r.begin(); i != r.end(); ++i)
			      {
				  compute_weight_of_signature(kept_kmers[i]);
			      }
			  });
*/
    
    std::cerr << "done\n";

#if 0
    if (!final_kmers.empty())
    {
	std::cerr << "writing kmers to " << final_kmers << "\n";
	fs::ofstream kf(final_kmers);
	std::for_each(kept_kmers.begin(), kept_kmers.end(), [&kf](const KeptKmer &k) {
		kf << k.kmer << "\t" << k.avg_from_end << "\t" << k.function_index << "\t" << k.weight << "\t" << k.otu_index << "\n";
		// kf << "\t" << k.seqs_containing_sig << "\t" << kmer_stats.seqs_with_func[k.function_index] << "\n";
	    });
	std::cerr << "done\n";
    }
#endif

    /*
     * Write kmer_stats.distinct_functions table ; this allows us
     * to reason about the coverage of functions in the database, as well as
     * which fusions have data for the component functions.
     */

    {
	fs::ofstream dfstr(kmer_data_dir / "distinct_functions");
	for (auto ent: kmer_stats.distinct_functions)
	{
	    dfstr << ent.first << "\t" << fm.lookup_function(ent.first) << "\t" << ent.second << "\n";
	}
    }

    /*    
    std::cerr << "writing hashtable to " << kmer_data_dir << " ...\n";
    KmerGuts *kguts = write_hashtable(kmer_data_dir, kept_kmers);
    std::cerr << "writing hashtable to " << kmer_data_dir << " done\n";
    */

    /*
     * Recall original data
     */

    auto hit_cb = [&fm](const hit_in_sequence_t<KeptKmerDB> &hit) {
	auto kd = hit.kdata;
	if (false)
	{
	    std::cout << hit.kmer << "\t" << fm.lookup_function(kd->function_index) << "\t" << kd->median << "\t" << kd->mean << "\t" << kd->var << "\t" << sqrt(kd->var) << "\t" << "\n";
	}
    };

    struct call_data
    {
	std::string id;
	std::string old_func;
	std::string old_func_stripped;
	std::string new_func;
	int func_index;
	float score;
    };
    tbb::concurrent_map<std::string, call_data> recall_report;

    auto call_cb = [&fm, &recall_report](const std::string &id, const std::string &func, int func_index, float score) {

	std::string orig, orig_stripped;
	fm.lookup_original_assignment(id, orig, orig_stripped);
	if (orig_stripped != func)
	{
	    recall_report.emplace(id, call_data { id, orig, orig_stripped, func, func_index, score});
	    // std::cout << "CALL "  << id << "\t" << orig_call << "\t" << func << "\t" << func_index << "\t" << score << "\n";
	}
	
    };

    struct saver
    {
	FunctionMap &fm;
	std::map<std::string, call_data> data;
	
	void operator()(const std::string &id, const std::string &func, int func_index, float score) {

	    std::string orig, orig_stripped;
	    fm.lookup_original_assignment(id, orig, orig_stripped);
	    if (orig_stripped != func)
	    {
		data.emplace(id, call_data { id, orig, orig_stripped, func, func_index, score});
		// std::cout << "CALL "  << id << "\t" << orig_call << "\t" << func << "\t" << func_index << "\t" << score << "\n";
	    }
	}
    };

    fs::path report_dir = kmer_data_dir / "recall.report.d";
    if (!fs::create_directory(report_dir))
    {
	std::cerr << "mkdir " << report_dir << " failed\n";
    }

    std::string fi_file = (kmer_data_dir / "function.index").string();
    KeptKmerDB  kdb(kept_kmers);
    KmerGeneric<KeptKmerDB> kmer_caller(kdb, fi_file);
    tbb::parallel_for(all_fasta_data.range(), [&kmer_caller, &hit_cb, &call_cb, &fm, &report_dir](auto r) {
	for (auto file: r)
	{
	    fs::path outfile(report_dir / file.filename());

	    saver s { fm } ;
	    fs::ifstream ifstr(file);

	    kmer_caller.process_fasta_stream(ifstr, hit_cb, s);
	    ifstr.close();

	    fs::ofstream ofstr(outfile);
	    for (auto ent: s.data)
	    {
		call_data &c = ent.second;
		ofstr << ent.first << "\t" << c.old_func << "\t" << c.old_func_stripped << "\t" << c.new_func << "\t" << c.func_index << "\t" << c.score << "\n";
	    }
	    ofstr.close();
	    
	}
    });

    /*
    {
	std::ofstream out("recall.report");
	for (auto ent: recall_report)
	{
	    call_data &c = ent.second;
	    out << ent.first << "\t" << c.old_func << "\t" << c.new_func << "\t" << c.func_index << "\t" << c.score << "\n";
	}
    }
    */
    
    if (!nudb_file.empty())
    {
	std::cerr << "write nudb data " << nudb_file << "\n";
	write_nudb_data(nudb_file, kept_kmers);
    }

    std::cerr << "all done\n";

    show_ps();
    rejected_stream.close();
    kept_stream.close();
    kept_function_stream.close();
    return 0;
}
