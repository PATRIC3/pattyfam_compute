#ifndef _kmer_generic_h
#define _kmer_generic_h

#define _FILE_OFFSET_BITS 64
#include <string>
#include <memory>
#include <map>
#include <vector>
#include <functional>
#include <algorithm>
#include <ostream>

#include <cstring>
#include <algorithm>
#include <memory>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <ctype.h>
#include <fcntl.h>
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>

#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>
#include <boost/math/statistics/univariate_statistics.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "kmer_params.h"
#include "kmer_types_generic.h"
#include "operators.h"

#include "fasta_parser.h"

#include "global.h"
#include "codet.h"

namespace KmerGenericNS
{

template <class Caller>
class KmerGeneric
{
public:
    
    Caller &caller_;
    KmerParameters parameters_;
    using KData = typename Caller::KData;

    KmerGeneric(Caller &caller, const std::string &function_index_file,
		int min_hits = 5, float min_weighted_hits = 0.0, int max_gap = 200);
    ~KmerGeneric(); 

    void gather_hits(const std::string &seqstr,
		     std::shared_ptr<std::vector<KmerCall>> calls,
		     std::function<void(const hit_in_sequence_t<Caller> &)> hit_cb,
		     std::shared_ptr<KmerOtuStats> otu_stats);
    void process_seq(const char *id,const char *data,
		     std::shared_ptr<std::vector<KmerCall>> calls,
		     std::function<void(const hit_in_sequence_t<Caller> &)> hit_cb,
		     std::shared_ptr<KmerOtuStats> otu_stats);
	
    void process_aa_seq(const std::string &id, const std::string &seq,
			std::shared_ptr<std::vector<KmerCall>> calls,
			std::function<void(const hit_in_sequence_t<Caller> &)> hit_cb,
			std::shared_ptr<KmerOtuStats> otu_stats);

    void process_aa_seq_hits(const std::string &id, const std::string &seq,
			std::shared_ptr<std::vector<KmerCall>> calls,
			     std::shared_ptr<std::vector<hit_in_sequence_t<Caller>>> hits,
			std::shared_ptr<KmerOtuStats> otu_stats);

    std::string format_call(const KmerCall &c);
    std::string format_hit(const hit_in_sequence_t<Caller> &h);
    std::string format_otu_stats(const std::string &id, size_t size, KmerOtuStats &otu_stats);

    void find_best_call(const std::string &id, std::vector<KmerCall> &calls, FunctionIndex &function_index,
			std::string &function, float &score, float &weighted_score, float &score_offset,
			bool call_ambigs = false);

    void read_function_index(const std::string &function_index_file);

    const std::string &function_at_index(int idx) {
	if (idx == UndefinedFunction)
	    return undefined_function_;
	else
	    return function_index_[idx];
    }

    /* Fasta processing methods */

    template <typename HitCB, typename CallCB>
    void process_fasta_stream(std::istream &istr, HitCB &hit_cb, CallCB &call_cb);

private:

    bool order_constraint_;
    int min_hits_;
    float min_weighted_hits_;
    int max_gap_;

    std::vector<std::string> function_index_;
    std::string undefined_function_;
};

#include "kmer_generic.tcc"

}

#endif /* _kmer_generic_h */
