#ifndef _kmer_nudb_h
#define _kmer_nudb_h

#define _FILE_OFFSET_BITS 64
#include <stddef.h>
#include <netinet/in.h>
#include <string>
#include <memory>
#include <map>
#include <vector>
#include <functional>
#include <algorithm>
#include <ostream>


#include "kmer_params.h"
#include "kmer_types_generic.h"
#include "nudb_kmer_db.h"

class KmerNudb
{
public:

    NuDBKmerDb<KMER_SIZE> &db_;
    KmerParameters parameters_;

    KmerNudb(NuDBKmerDb<KMER_SIZE> &db, const std::string &function_index_file,
	     int min_hits = 5, float min_weighted_hits = 0.0, int max_gap = 200);
    ~KmerNudb(); 

    void gather_hits(const std::string &seqstr,
		     std::shared_ptr<std::vector<KmerCall>> calls,
		     std::function<void(const hit_in_sequence_t &)> hit_cb,
		     std::shared_ptr<KmerOtuStats> otu_stats);
    void process_seq(const char *id,const char *data,
		     std::shared_ptr<std::vector<KmerCall>> calls,
		     std::function<void(const hit_in_sequence_t &)> hit_cb,
		     std::shared_ptr<KmerOtuStats> otu_stats);
	
    void process_aa_seq(const std::string &id, const std::string &seq,
			std::shared_ptr<std::vector<KmerCall>> calls,
			std::function<void(const hit_in_sequence_t &)> hit_cb,
			std::shared_ptr<KmerOtuStats> otu_stats);

    void process_aa_seq_hits(const std::string &id, const std::string &seq,
			std::shared_ptr<std::vector<KmerCall>> calls,
			std::shared_ptr<std::vector<hit_in_sequence_t>> hits,
			std::shared_ptr<KmerOtuStats> otu_stats);

    std::string format_call(const KmerCall &c);
    std::string format_hit(const hit_in_sequence_t &h);
    std::string format_otu_stats(const std::string &id, size_t size, KmerOtuStats &otu_stats);

    void find_best_call(std::vector<KmerCall> &calls, FunctionIndex &function_index,
			std::string &function, float &score, float &weighted_score, float &score_offset);

    void read_function_index(const std::string &function_index_file);

    const std::string &function_at_index(int idx) {
	if (idx == UndefinedFunction)
	    return undefined_function_;
	else
	    return function_index_[idx];
    }

private:

    bool order_constraint_;
    int min_hits_;
    float min_weighted_hits_;
    int max_gap_;

    std::vector<std::string> function_index_;
    std::string undefined_function_;
};

template <int N, typename F>
void for_each_kmer(const std::string &str, F cb) {
    const char *ptr = str.c_str();
    const char *end = ptr  + str.length();
    const char *last_kmer= end - N;

    const char *next_ambig = std::find_if(ptr, end, [](char c) -> bool { return c == '*' || c == 'X'; });
    // std::cerr << "next_ambig=" << next_ambig << "\n";
    std::array<char, N> kmer;
    while (ptr <= last_kmer)
    {
	// std::cerr << "ptr=" << ptr << "\n";
	const char *kend = ptr + N;
	// std::cerr << "kend=" << kend << "\n";
	if (next_ambig != end && kend >= next_ambig)
	{
	    // std::cerr << "hit abmig\n";
	    ptr = next_ambig + 1;
	    next_ambig = std::find_if(ptr, end, [](char c) { return c == '*' || c == 'X'; });
	    continue;
	}
	std::copy(ptr, kend, kmer.data());
	// std::cerr << "cb " << kmer << "\n";
	cb(kmer, ptr - str.c_str());
	ptr++;
    }
}

inline std::ostream &operator<<(std::ostream &os, const KmerCall &c)
{
    os << "KmerCall(" << c.start << "-" << c.end << ": " << c.count << ", " << c.function_index
       << ", " << c.weighted_hits
       << ", " << c.r2
       << ")";
    return os;
}

#endif /* _kmer_guts_h */
