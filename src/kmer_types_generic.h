#ifndef _kmer_types_generic_h
#define _kmer_types_generic_h

#include <climits>
#include <map>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <sstream>
#include <iostream>
#include <iterator>

#include "kmer_value_types.h"

class KmerResult
{
public:
    std::string id;
    size_t len;
    
};


/*
 * For compatibility with current kmer builds,
 * we use the original kguts types here rather than
 * the shorter OTUIndex and FunctionIndex.
 */
typedef struct sig_kmer {
    unsigned long long  which_kmer;
    int  otu_index;
    unsigned short  avg_from_end;
    int  function_index;
    float weight;
} sig_kmer_t;

// NuDB data struct

struct KData {
    OTUIndex otu_index;
    unsigned short  avg_from_end;
    FunctionIndex function_index;
    float weight;
};

// V690 The 'hit_in_sequence_t' class implements a copy constructor, but lacks the copy assignment operator. It is dangerous to use such a class.
template <class Caller>
struct hit_in_sequence_t {
    sig_kmer_t hit;
    unsigned int offset;

    typedef std::array<char, Caller::KmerSize> kmer_type;
    typedef typename Caller::KData kdata_type;
    kmer_type kmer;
    const kdata_type *kdata;
    
    hit_in_sequence_t(const sig_kmer_t &h, unsigned int o, const kmer_type &k, const kdata_type *kd)
	: hit(h), offset(o), kmer(k), kdata(kd) {}
    hit_in_sequence_t(const sig_kmer_t &h, unsigned int o) : hit(h), offset(o), kdata(0) {}
    hit_in_sequence_t(const hit_in_sequence_t &h) : hit(h.hit), offset(h.offset), kmer(h.kmer), kdata(h.kdata) {}

    hit_in_sequence_t(unsigned long long enc, const KData &k, unsigned int o) :
	hit {enc, k.otu_index, k.avg_from_end, k.function_index, k.weight}, offset{o}, kdata(0) {}

};

class KmerHit
{
public:
    OTUIndex   oI;
    unsigned int from0_in_prot;      /* offset from start of protein sequence */
    unsigned short avg_off_from_end;  /* average offset from the end */
    FunctionIndex fI;
    float weight;
    unsigned long long encoded_kmer;
};
typedef KmerHit hit_t;

class KmerCall
{
public:
    unsigned int start;
    unsigned int end;
    int count;
    FunctionIndex function_index;
    float weighted_hits;
    //R2 double r2;
    unsigned int protein_length_median;
    float protein_length_med_avg_dev;

    KmerCall() : start(0), end(0), count(0), function_index(UndefinedFunction), weighted_hits(0.0),
	//R2 r2(0.0),
	protein_length_median(0), protein_length_med_avg_dev(0.0) {
    }
    KmerCall(unsigned int s, unsigned int e, int c, FunctionIndex f, float w,
	     //R2 double r = 0.0,
	     unsigned int med = 0, float mad = 0.0)  :
    start(s), end(e), count(c), function_index(f), weighted_hits(w),
	//R2 r2(r),
	protein_length_median(med), protein_length_med_avg_dev(mad) { }
    KmerCall(KmerCall &&k) :
        start(k.start), end(k.end), count(k.count),
	function_index(k.function_index), weighted_hits(k.weighted_hits),
	    //R2 r2(k.r2),
	    protein_length_median(k.protein_length_median), protein_length_med_avg_dev(k.protein_length_med_avg_dev) { }
    KmerCall(const KmerCall &k) :
        start(k.start), end(k.end), count(k.count), function_index(k.function_index),
	    weighted_hits(k.weighted_hits),
	    //R2 r2(k.r2),
	    protein_length_median(k.protein_length_median), protein_length_med_avg_dev(k.protein_length_med_avg_dev)  { }

};

inline std::ostream &operator<<(std::ostream &os, const KmerCall &c)
{
    os << "KmerCall(" << c.start << "-" << c.end << ": " << c.count << ", " << c.function_index
       << ", " << c.weighted_hits
	//R2 << ", " << c.r2
       << ", " << c.protein_length_median
       << ", " << c.protein_length_med_avg_dev
       << ")";
    return os;
}

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

class KmerOtuStats
{
public:
    std::string contig_id;
    int contig_len;

    std::map<int, int> otu_map;
    std::vector<std::pair<int, int>> otus_by_count;
		    

    template<class T>
	struct less_second : std::binary_function<T,T,bool>
    {
	inline bool operator()( const T& lhs, const T& rhs )
	{
	    return rhs.second < lhs.second;
	}
    };

    void write(FILE *fh)
    {
	fprintf(fh, "OTU-COUNTS\t%s[%d]", contig_id.c_str(), contig_len);
	for (auto it = otus_by_count.begin(); it != otus_by_count.end(); it++)
	{
	    fprintf(fh, "\t%d-%d",it->second, it->first);
	}
	fprintf(fh, "\n");
    }

    void finalize()
    {
	otus_by_count.insert(otus_by_count.begin(), otu_map.begin(), otu_map.end());
	std::sort(otus_by_count.begin(), otus_by_count.end(), less_second<std::pair< int, int> > ());
    }
};

class KmerParameters
{
public:
KmerParameters() :
    param_map_(
	{
	    { "order_constraint", order_constraint },
	    { "min_hits", min_hits },
	    { "min_weighted_hits", min_weighted_hits },
	    { "max_gap", max_gap }
	})

    {
	set_default_parameters();
    }
    
    void set_default_parameters()
    {
	order_constraint = 0;
	min_hits = 5;
	min_weighted_hits = 0;
	max_gap  = 200;
    }

    std::map<std::string, int &> param_map_;
    void set_parameters(const std::map<std::string, std::string> &params)
    {
	set_default_parameters();
	for (auto p: params)
	{
	    auto myp = param_map_.find(p.first);
	    if (myp != param_map_.end())
	    {
		try {
		    int val = std::stoi(p.second);
		    myp->second = val;
		} catch (const std::invalid_argument& ia)
		  {
		      std::cerr << "Warning: invalid integer '" << p.second << "' passed for parameter " << p.first << "\n";
		  }
	    }
	}
    }

    int   order_constraint;
    int   min_hits;
    int   min_weighted_hits;
    int   max_gap;
};

template <class T, std::size_t N>
inline std::ostream& operator<<(std::ostream& o, const std::array<T, N>& arr)
{
    std::copy(arr.cbegin(), arr.cend(), std::ostream_iterator<T>(o, ""));
    return o;
}

std::string format_otu_stats(const std::string &id, size_t size, KmerOtuStats &otu_stats)
{
    std::ostringstream oss;
    oss << "OTU-COUNTS\t" << id << "[" << size << "]";

    int max_to_print = 5;
    for (auto x = otu_stats.otus_by_count.begin();
	 x != otu_stats.otus_by_count.end() && max_to_print > 0; x++, max_to_print--)
    {
	oss << "\t" << x->second << "-" << x->first;
    }
    oss << "\n";
    return oss.str();
}


template <std::size_t N>
inline bool operator==(const char *str, const std::array<char, N> &arr)
{
    if (strlen(str) != N)
	return false;
    return strcmp(str, arr.data()) == 0;
}

template <std::size_t N>
inline bool operator==(const std::array<char, N> &arr, const char *str)
{
    return str == arr;
}

#endif /* _kmer_types_generic_h */
