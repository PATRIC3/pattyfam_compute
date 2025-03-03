#ifndef _kmer_guts_h
#define _kmer_guts_h

/*

kmer_guts.c can be compiled into either a 5-mer or an 8-mer version.  I have 
labeled the critical changes with 

     ### CHANGE THIS FOR DIFFERING Ks ###

Note: to run use

   build_data_directory KmerData 8

Then use

   kmer_guts -w -D KmerData < input.contigs

to build an memory map -- it can take 10-15 minutes.

Then use

   kmer_guts -D KmerData < input.contigs

to run.  Eventually, the memory-mapped hash table, which is large,
will become resident, and these should go fairly quickly.

   kmer_guts -D KmerData -a  < amino_acid_sequences.fasta > hits

can be used to call translated PEGs.

where KmersData is a directory that must contain

       final.kmers     [a file of [kmer,avg-offset-from-end,function-index,otu-index]]
       function.index  [a file of [index,function] pairs]
       otu.index       [a file of [index,otu] pairs]
----------------

Conceptually, the data associated with each signature Kmer is

        1. the protein kmer
        2. the average offset from the end of the protein
        3. a set of numeric values that include

                function index
                OTU index

The program takes as input a  file that should be thought of as a sequence
of fasta files in which a single line

          >FLUSH

occurs after each infut file.  That is,

          >contig1
	  acgtacgt
	  >FLUSH
	  >contig2
	  acgtacgt

is an input file containing two input fasta files (each containing a single, short
contig).  The '>FLUSH' causes all files to be flushed.

In effect, the program processes a sequence of requests.  The output for each request
is a piece of stadout terminated by a line

          //

The lines in the output for a request are of four forms.  The first type of line begins
the processing of a contig and looks like

         processing contig[length]

Then, there will be a message like

         TRANSLATION contig length-contig frame offset-of-start  

before beginning output for each of the six frames.

Thus, there will be six such lines printed for each input contig.  Following
each of these "TRANSLATION..." lines, there will be a set of lines like

         CALL Start-in-protein-seq end-in-protein-seq number-hits function-index function

Finally, after all six frames have been processed (for DNA - just one sequence for aa input),
you get a line of the form

         OTU-COUNTS cnt1-otu1 cnt2-otu2 ...

This code uses a table indicating which K-mers are signatures. I call this the 
"kmer_bits" table.  You can run this code for K ranging from 5 to 8.

#################
The line

    #define K 8

defines the kmer size.  
The code is intended to work with only 5-mers or 8-mers.  That is:

YOU MUST SET K TO EITHER 5 OR 8 (sorry).
############################################

COMMAND LINE ARGUMENTS:

    -a         means amino acid input sequence (defaults to DNA)

    -d Level   sets debugging level (1 shows hits; after that it gets intense

    -m MinHits minimum number of hits to get CALLed

    -M MinWeighted  Kmers are now weighted; this is the min of the sum of the weights,

    -O         use order contraint

    -g MaxGap  sets maximum allowed gap between HITS

    -D Data    sets the Data directory where the memory map lives

    -s HashSize make sure that the value is the same when you save the memory map
                and when you use it to search

    -w          write the memory map (means Data must contain final.kmers and the indexes

    -l port	Run in server mode, listening on the given port. If port = 0, pick a port

    -L pfile	When running in server mode, write the port number into the given file
*/


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

#include "kmer_image.h"
#include "kmer_params.h"
#include "kmer_encoder.h"
#include "kmer_types.h"

// V690 The 'KmerGuts' class implements a copy constructor, but lacks the copy assignment operator. It is dangerous to use such a class.
class KmerGuts
{
public:

    KmerEncoder encoder_;
    typedef struct sig_kmer sig_kmer_t;
    using hit_in_sequence_t = ::hit_in_sequence_t;

    typedef struct kmer_handle {
	sig_kmer_t *kmer_table;
	unsigned long long num_sigs;
	char **function_array;   /* indexed by fI */
	char **otu_array;        /* OTU indexes point at a representation of multiple OTUs */
	int function_count;
	int otu_count;
    kmer_handle() : kmer_table(0), num_sigs(0), function_array(0), otu_array(0), function_count(0), otu_count(0) {}
    } kmer_handle_t;

/* the following stuff was added to condense sets of hits to specific calls.
   The basic strategy is to use a set of global variables to retain state and flush
   calls (bad, bad, bad...).
*/

    struct otu_count {
	OTUIndex oI;
	int count;
    };


/* parameters to main -- accessed globally */
    int debug;
    int aa;
    int hits_only;
    unsigned long long size_hash; /* 1400303159  tot_lookups=13474100 retry=2981020 for 5.contigs 4.684 sec */
    /* 2147483648  tot_lookups=13474100 retry=1736650  */
    /* 1073741824  tot_lookups=13474100 retry=4728020  */
    
    hit_t hits[MAX_HITS_PER_SEQ]; 
    int   num_hits;
    struct otu_count oI_counts[OI_BUFSZ];
    int num_oI;

    /*
     * Count used for loading.
     */
    unsigned long long kmers_loaded_;
    kmer_memory_image_t *kmer_image_for_loading_;
    unsigned long long image_size_for_loading_;
    sig_kmer_t *sig_kmers_for_loading_;

    void insert_kmer(const std::string &kmer,
		     FunctionIndex function_index, OTUIndex otu_index, unsigned short avg_offset,
		     float function_weight);
    void insert_kmer(unsigned long long encodedK,
		     FunctionIndex function_index, OTUIndex otu_index, unsigned short avg_offset,
		     float function_weight);
    void save_kmer_hash_table(const std::string &file);
     
    
    FunctionIndex   current_fI;
    char  current_id[300];
    size_t   current_length_contig;
    char  current_strand;
    unsigned short current_prot_off;
    int   order_constraint;
    int   min_hits;
    int   min_weighted_hits;
    int   max_gap;

    void set_default_parameters();
    std::map<std::string, int &> param_map_;
    void set_parameters(const std::map<std::string, std::string> &params);

    unsigned char *pIseq;
    char *data;
    char *cdata;
    char *pseq;
    // Work around lost pointers to internal allocation from the original kguts code.
    // This is a list of pointers to be freed when we're done.
    std::vector<char *> to_be_freed;

    int tot_lookups;
    int retry;

    std::shared_ptr<kmer_handle_t> kmersH;

    KmerGuts(const std::string &kmer_dir, std::shared_ptr<KmerImage> image);
    KmerGuts(const std::string &kmer_dir, unsigned long long num_buckets);
    KmerGuts(KmerGuts &);
    ~KmerGuts(); 
    void do_init();

    static unsigned char to_amino_acid_off(char c);
    char comp(char c);
    void rev_comp(const char *data,char *cdata);
    static unsigned long long encoded_kmer(unsigned char *p);
    static unsigned long long encoded_aa_kmer(const char *p);
    static void decoded_kmer(unsigned long long encodedK,char *decoded);
    int dna_char(char c);
    void translate(const char *seq,int off,char *pseq, unsigned char *pIseq);
    char **load_indexed_ar(const char *filename,int *sz);
    char **load_functions(const char *file, int *sz);
    char **load_otus(const char *file, int *sz);
    unsigned long long find_empty_hash_entry(sig_kmer_t sig_kmers[], unsigned long long encodedK);
    unsigned long long lookup_hash_entry(sig_kmer_t sig_kmers[], unsigned long long encodedK);
    kmer_memory_image_t *load_raw_kmers(char *file,unsigned long long num_entries, unsigned long long *alloc_sz);
    void advance_past_ambig(unsigned char **p,unsigned char *bound);

    void process_set_of_hits(std::shared_ptr<std::vector<KmerCall>> calls,
				       std::shared_ptr<KmerOtuStats> otu_stats);
    void gather_hits(const char *pseq,
		     unsigned char *pIseq,
		     std::shared_ptr<std::vector<KmerCall>> calls,
		     std::function<void(const hit_in_sequence_t &)> hit_cb,
		     std::shared_ptr<KmerOtuStats> otu_stats);
	
    void process_seq(const char *id,const char *data,
		     std::shared_ptr<std::vector<KmerCall>> calls,
		     std::function<void(const hit_in_sequence_t &)> hit_cb,
		     std::shared_ptr<KmerOtuStats> otu_stats);
	
    void process_aa_seq(const std::string &id, const std::string &seq,
			std::shared_ptr<std::vector<KmerCall>> calls,
			std::function<void(const hit_in_sequence_t&)> hit_cb,
			std::shared_ptr<KmerOtuStats> otu_stats);

    void process_aa_seq_hits(const std::string &id, const std::string &seq,
			std::shared_ptr<std::vector<KmerCall>> calls,
			std::shared_ptr<std::vector<hit_in_sequence_t>> hits,
			std::shared_ptr<KmerOtuStats> otu_stats);

    std::shared_ptr<kmer_handle> init_kmers(const char *dataD);

    std::shared_ptr<KmerImage> image_;

    const char *function_at_index(int i) {
	if (i < 0 || i >= kmersH->function_count)
	    return "INVALID_OFFSET";
	else
	    return kmersH->function_array[i];
    }

    std::string format_call(const KmerCall &c);
    std::string format_hit(const hit_in_sequence_t &h);
    std::string format_otu_stats(const std::string &id, size_t size, KmerOtuStats &otu_stats);

    void find_best_call(std::vector<KmerCall> &calls, FunctionIndex &function_index, std::string &function, float &score, float &weighted_score, float &score_offset);
};

/* redefinition
inline std::ostream &operator<<(std::ostream &os, const KmerCall &c)
{
    os << "KmerCall(" << c.start << "-" << c.end << ": " << c.count << ", " << c.function_index << ", " << c.weighted_hits << ")";
    return os;
}

*/

#endif /* _kmer_guts_h */
