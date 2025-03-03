#define DEBUG_SCORING 0

#include "nudb_kmer_db.h"
#include "kmer_generic.h"
#include "fasta_parser.h"
#include "clamp.h"
#include <memory>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
//#include <boost/math/statistics/linear_regression.hpp>

#include <stdexcept>

namespace fs = boost::filesystem;
namespace po = boost::program_options;

using namespace KmerGenericNS;


bool process_command_line_options(int argc, char *argv[],
				  fs::path &kmer_data_dir,
				  fs::path &nudb_file,
				  bool &show_kmers,
				  bool &call_ambigs,
				  std::vector<fs::path> &input_files,
				  fs::path &output_file)
{

    std::ostringstream x;
    x << "Usage: " << argv[0] << " [options]\nAllowed options";
    po::options_description desc(x.str());
    po::options_description hidden("hidden options");
    po::options_description all("all options");

    hidden.add_options()
	("kmer-data-dir", po::value<fs::path>(&kmer_data_dir), "Read kmer data files from this directory.")
	("nudb-file", po::value<fs::path>(&nudb_file), "NuDB kmers database.")
	("input-files", po::value<std::vector<fs::path>>(&input_files), "NuDB kmers database.");

    desc.add_options()
	("output-file", po::value<fs::path>(&output_file), "Write output to this path")
	("show-kmers", po::bool_switch(&show_kmers), "Show all kmer hits.")
	("call-ambigs", po::bool_switch(&call_ambigs), "Call ambiguous calls with ?? notation.")
	("help,h", "show this help message");

    all.add(desc).add(hidden);

    po::positional_options_description pos;
    pos.add("kmer-data-dir", 1);
    pos.add("nudb-file", 1);
    pos.add("input-files", -1);

    po::variables_map vm;

    po::store(po::command_line_parser(argc, argv).options(all).positional(pos).run(), vm);

    if (vm.count("help"))
    {
	std::cout << all << "\n";
	return false;
    }

    po::notify(vm);

    return true;
}


int main(int argc, char **argv)
{
    fs::path kmer_data_dir;
    fs::path nudb_file;
    bool show_kmers;
    bool call_ambigs;
    std::vector<fs::path> input_files;
    fs::path output_file;

    process_command_line_options(argc, argv, kmer_data_dir, nudb_file, show_kmers, call_ambigs, input_files, output_file);

    const int K = 8;
    
    fs::path function_index_file = kmer_data_dir / "function.index";
    typedef NuDBKmerDb<K> KDB;
    KDB db(nudb_file.string());

    if (!db.exists())
    {
	std::cerr << "database does not exist\n";
	exit(1);
	// db.create();
    }
    db.open();

//    KmerNudb kmer_nudb(db, function_index_file);
    KmerGeneric<NuDBKmerDb<K>> kmer_nudb(db, function_index_file.string());

    std::ostream* out_fp = &std::cout;
    std::ofstream fout;
    if (!output_file.empty())
    {
	fout.open(output_file.string());
	out_fp = &fout;
    }

    try {
	FastaParser parser;
    
	parser.set_callback([&kmer_nudb, show_kmers, call_ambigs, out_fp](const std::string &id, const std::string &seq) {
	
	    double slen = static_cast<double>(seq.length());
	    auto calls = std::make_shared<std::vector<KmerCall>>();
	    // auto hits = std::make_shared<std::vector<hit_in_sequence_t>>();
	    auto ostats = std::make_shared<KmerOtuStats>();
	
	    auto hit_cb = [ slen, &kmer_nudb, show_kmers, out_fp](const hit_in_sequence_t<NuDBKmerDb<K>> &hit) {
		auto kd = hit.kdata;
		if (show_kmers)
		{
		    *out_fp <<
			hit.kmer << "\t" <<
			kmer_nudb.function_at_index(kd->function_index) << "\t" <<
			kd->median << "\t" <<
			kd->mean << "\t" <<
			kd->var << "\t" <<
			sqrt(kd->var) << "\t" <<
			kd->weight << "\t" <<
			+clamp_float<uint8_t>(30.0f * kd->weight) << "\t" << 
			kd->n_proteins <<
			"\n";
		}
	    };
	    kmer_nudb.process_aa_seq(id, seq, calls, hit_cb, ostats);
	    #if DEBUG_SCORING
	    std::cout << id << " calls:\n";
	    for (auto c: *calls)
	    {
		*out_fp << c.start << "-" << c.end << "\t" << c.count << "\t" << kmer_nudb.function_at_index(c.function_index)
			  << "\t" << c.function_index
			  << "\t" << c.weighted_hits
		    //R2 << "\t" << c.r2
			  << "\n";

	    }
	    #endif
	    FunctionIndex fi;
	    std::string func;
	    float score;
	    float wt;
	    float offset;
	    kmer_nudb.find_best_call(id, *calls, fi, func, score, wt, offset, call_ambigs);
	    *out_fp << id << "\t" << func << "\t" << fi << "\t" << score << "\n";
	    return 0;
	});

	for (auto file: input_files)
	{
	    std::cerr << "Process " << file << "\n";
	    fs::ifstream ifstr(file);
	    parser.parse(ifstr);
	    parser.parse_complete();
	}
    }
    catch (std::runtime_error &x)
    {
	std::cerr<< "caught " << x.what() << "\n";
	exit(1);
    }
}

