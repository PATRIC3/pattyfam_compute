#include "kguts.h"

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
#include <string.h>

#include <boost/program_options.hpp>
#include "global.h"

static const  char genetic_code[64] = {
    'K','N','K','N','T','T','T','T','R','S','R','S','I','I','M','I',
    'Q','H','Q','H','P','P','P','P','R','R','R','R','L','L','L','L',
    'E','D','E','D','A','A','A','A','G','G','G','G','V','V','V','V',
    '*','Y','*','Y','S','S','S','S','*','C','W','C','L','F','L','F'
};

static const  char prot_alpha[20] = {
    'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y' 
};

KmerGuts::KmerGuts(const std::string &data_dir, std::shared_ptr<KmerImage> image) :
    param_map_(
	{
	    { "order_constraint", order_constraint },
	    { "min_hits", min_hits },
	    { "min_weighted_hits", min_weighted_hits },
	    { "max_gap", max_gap }
	}),
    image_(image)
{
    kmersH = init_kmers(data_dir.c_str());
    
    tot_lookups = 0;
    retry  = 0;

    debug = 0;
    aa = 0;
    hits_only = 0;
    size_hash = image->image()->num_sigs;

    set_default_parameters();

    do_init();
}

void KmerGuts::do_init()
{
    pIseq = (unsigned char *) malloc(MAX_SEQ_LEN / 3);
    cdata = (char *) malloc(MAX_SEQ_LEN);
    data = (char *) malloc(MAX_SEQ_LEN);
    pseq = (char *) malloc(MAX_SEQ_LEN / 3);

    num_hits = 0;
    num_oI = 0;

    current_id[0] = 0;
}

/*
 * Alternative constructor used for loading hash table from
 * kmers and saving to disk.
 */
KmerGuts::KmerGuts(const std::string &data_dir, unsigned long long num_buckets) :
    size_hash(num_buckets),
    kmers_loaded_(0),
    image_(0)
{

    /*
     * Initialize our image pointer using the num_buckets passed in.
     *
     * This is code hoisted from the original init_kmer call. We make
     * the operation of writing hash tables explicit.
     */

    unsigned long long image_size = sizeof(kmer_memory_image_t) + (sizeof(sig_kmer_t) * num_buckets);
    kmer_memory_image_t *image = static_cast<kmer_memory_image_t *>(malloc(image_size));
    if (image == 0)
    {
	std::cerr << "Memory allocation of size " << image_size << " failed\n";
	exit(1);
    }
    std::cerr << "Allocated image of size " << image_size << "\n";
    memset(image, 0, image_size);
    image->num_sigs = num_buckets;
    image->entry_size = sizeof(sig_kmer_t);
    image->version = (long long) KMER_VERSION;

    std::shared_ptr<kmer_handle_t> handle = std::make_shared<kmer_handle_t>();
    handle->kmer_table = (sig_kmer_t *) (image + 1);
    handle->num_sigs   = image->num_sigs;

    kmersH = handle;

    sig_kmer_t *sig_kmers = (sig_kmer_t *) (image + 1);

    for (unsigned long long i =0; (i < size_hash); i++)
	sig_kmers[i].which_kmer = MAX_ENCODED + 1;

    kmer_image_for_loading_ = image;
    sig_kmers_for_loading_ = sig_kmers;
    image_size_for_loading_ = image_size;

    set_default_parameters();
    do_init();
}

KmerGuts::KmerGuts(KmerGuts &k)
{
    std::cerr << "copy construct kmer_guts to create " << this << " from " << &k << "\n";

    debug = k.debug;
    aa = k.aa;
    hits_only = k.hits_only;
    size_hash = k.size_hash;
    num_hits = k.num_hits;
    for (int i = 0; i < num_hits; i++)
	hits[i] = k.hits[i];
    num_oI = k.num_oI;
    for (int i = 0; i < num_oI; i++)
	oI_counts[i] = k.oI_counts[i];

    kmers_loaded_ = k.kmers_loaded_;
    kmer_image_for_loading_ = k.kmer_image_for_loading_;
    image_size_for_loading_= k.image_size_for_loading_;
    sig_kmers_for_loading_ = k.sig_kmers_for_loading_;

    current_fI = k.current_fI;
    strcpy(current_id, k.current_id);
    current_length_contig = k.current_length_contig;
    current_strand = k.current_strand;
    current_prot_off = k.current_prot_off;
    order_constraint = k.order_constraint;
    min_hits = k.min_hits;
    min_weighted_hits = k.min_weighted_hits;
    max_gap = k.max_gap;

    param_map_ = k.param_map_;

    pIseq = 0;
    data = 0;
    cdata = 0;
    pseq = 0;

    tot_lookups = k.tot_lookups;
    retry = k.retry;

    kmersH = k.kmersH;

    image_ = k.image_;

    // long long alloc_sz = sizeof(kmer_memory_image_t) + (sizeof(sig_kmer_t) * k.size_hash);

    do_init();
}

inline unsigned long long KmerGuts::find_empty_hash_entry(sig_kmer_t sig_kmers[], unsigned long long encodedK) {
    unsigned long long hash_entry = encodedK % size_hash;
    while (sig_kmers[hash_entry].which_kmer <= MAX_ENCODED)
      hash_entry = (hash_entry+1)%size_hash;
    return hash_entry;
}
KmerGuts::~KmerGuts() {
    std::cerr << "delete kguts " << this << "\n";
    for (auto x: to_be_freed)
	free(x);
    if (pIseq)
	free(pIseq);
    if (cdata)
	free(cdata);
    if (data)
	free(data);
    if (pseq)
	free(pseq);
}

void KmerGuts::insert_kmer(const std::string &kmer,
			   FunctionIndex function_index, OTUIndex otu_index, unsigned short avg_offset,
			   float function_weight)
{
    unsigned long long encodedK = encoder_.encoded_aa_kmer(kmer.c_str());
    insert_kmer(encodedK, function_index, otu_index, avg_offset, function_weight);
}

void KmerGuts::insert_kmer(unsigned long long encodedK,
			   FunctionIndex function_index, OTUIndex otu_index, unsigned short avg_offset,
			   float function_weight)
{
    if (encodedK > MAX_ENCODED)
    {
	// std::cerr << "Warn: " << kmer << " contains invalid characters\n";
	return;
    }
    unsigned long long hash_entry = find_empty_hash_entry(sig_kmers_for_loading_, encodedK);
    kmers_loaded_++;
    if (kmers_loaded_ >= (size_hash / 2)) {
	fprintf(stderr,"Your Kmer hash is half-full; use -s (and -w) to bump it\n");
	exit(1);
    }
    sig_kmers_for_loading_[hash_entry].which_kmer     = encodedK;
    sig_kmers_for_loading_[hash_entry].avg_from_end   = avg_offset;
    sig_kmers_for_loading_[hash_entry].function_index = function_index;
    sig_kmers_for_loading_[hash_entry].otu_index      = otu_index;
    sig_kmers_for_loading_[hash_entry].function_wt    = function_weight;
}

void KmerGuts::save_kmer_hash_table(const std::string &file)
{
    FILE *fp = fopen(file.c_str(),"w");
    if (fp == NULL) { 
	fprintf(stderr,"could not open %s for writing: %s ",file.c_str(), strerror(errno));
	exit(1);
    }
    std::cerr << "Writing " << image_size_for_loading_ <<" bytes\n";
    fwrite(kmer_image_for_loading_, image_size_for_loading_, 1, fp);
    fclose(fp);
}

void KmerGuts::set_default_parameters()
{
    order_constraint = 0;
    min_hits = 5;
    min_weighted_hits = 0;
    max_gap  = 200;
}

void KmerGuts::set_parameters(const std::map<std::string, std::string> &params)
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
    /*
    std::cerr << "Parameters:\n";
    for (auto p: param_map_)
    {
	std::cerr << p.first << "=" << p.second << "\n";
    }
    */
}


/* =========================== end of reduction global variables ================= */

inline unsigned char KmerGuts::to_amino_acid_off(char c) {
  switch (c)
    {
      case 'A':
	return 0;

      case 'C':
	return 1;

      case 'D':
	return 2;

      case 'E':
	return 3;

      case 'F':
	return 4;

      case 'G':
	return 5;

      case 'H':
	return 6;

      case 'I':
	return 7;

      case 'K':
	return 8;

      case 'L':
	return 9;

      case 'M':
	return 10;

      case 'N':
	return 11;

      case 'P':
	return 12;

      case 'Q':
	return 13;

      case 'R':
	return 14;

      case 'S':
	return 15;

      case 'T':
	return 16;

      case 'V':
	return 17;

      case 'W':
	return 18;

      case 'Y':
	return 19;

      default:
	return 20;
    }
}

char KmerGuts::comp(char c)
{
    switch (c)
    {
      case 'a':
	return 't';
      case 'A':
	return 'T';

      case 'c':
	return 'g';
      case 'C':
	return 'G';

      case 'g':
	return 'c';
      case 'G':
	return 'C';

      case 't':
      case 'u':
	return 'a';
      case 'T':
      case 'U':
	return 'A';

      case 'm':
	return 'k';
      case 'M':
	return 'K';

      case 'r':
	return 'y';
      case 'R':
	return 'Y';

      case 'w':
	return 'w';
      case 'W':
	return 'W';

      case 's':
	return 'S';
      case 'S':
	return 'S';

      case 'y':
	return 'r';
      case 'Y':
	return 'R';

      case 'k':
	return 'm';
      case 'K':
	return 'M';

      case 'b':
	return 'v';
      case 'B':
	return 'V';

      case 'd':
	return 'h';
      case 'D':
	return 'H';

      case 'h':
	return 'd';
      case 'H':
	return 'D';

      case 'v':
	return 'b';
      case 'V':
	return 'B';

      case 'n':
	return 'n';
      case 'N':
	return 'N';

      default:
	return c;
    }
}

void KmerGuts::rev_comp(const char *data,char *cdata) {

    size_t n = strlen(data);
    const char *p  = data + (n-1);
    char *pc = cdata;
    while (n--) {
	*(pc++) = (char) compl(*(p--));
    }
    *pc = 0;
}

unsigned long long KmerGuts::encoded_kmer(unsigned char *p) {
//    std::cout << "encode '" << p << "'\n";
  unsigned long long encodedK = *p;
  int i;
  for (i=1; (i <= KMER_SIZE-1); i++) {
    encodedK = (encodedK * 20) + *(p+i);
  }

  if (encodedK > MAX_ENCODED) {
    fprintf(stderr,"bad encoding - input must have included invalid characters\n");
    for (i=0; (i < KMER_SIZE); i++) {
      fprintf(stderr,"%d ",*(p+i));
    }
    fprintf(stderr,"\n");
    exit(2);
  }
  return encodedK;
}

unsigned long long KmerGuts::encoded_aa_kmer(const char *p)
{
    unsigned char aa_off[KMER_SIZE];
    int j;
    for (j=0; (j < KMER_SIZE); j++) {
	char prot_c = *(p+j);
	aa_off[j] = to_amino_acid_off(prot_c);
	/*
	 * If we have an invalid char, return an invalid encoding.
	 */
	if (aa_off[j] >= 20)
	    return MAX_ENCODED + 1;
    }
    return encoded_kmer(aa_off);
}

void KmerGuts::decoded_kmer(unsigned long long encodedK,char *decoded) {
  
  int i;
  *(decoded+KMER_SIZE) = '\0';
  unsigned long long x = encodedK;

  for (i=KMER_SIZE-1; (i >= 0); i--) {
    *(decoded+i) = prot_alpha[x % 20];
    x = x / 20;
  }
}


int KmerGuts::dna_char(char c)
{
    switch (c)
    {
      case 'a':
      case 'A':
	return 0;

      case 'c':
      case 'C':
	return 1;

      case 'g':
      case 'G':
	return 2;

      case 't':
      case 'u':
      case 'T':
      case 'U':
	return 3;

      default:
	return 4;;
    }
}

void KmerGuts::translate(const char *seq,int off,char *pseq, unsigned char *pIseq) {

  size_t max = strlen(seq) - 3;
  char *p = pseq;
  unsigned char *pI = pIseq;
  for (size_t i = (size_t) off; (i <= max); ) {
    int c1 = dna_char(seq[i++]);
    int c2 = dna_char(seq[i++]);
    int c3 = dna_char(seq[i++]);
    if ((c1 < 4) && (c2 < 4) && (c3 < 4)) {
      int I = (c1 * 16) + (c2 * 4) + c3;
      char prot_c = genetic_code[I];
      *(p++) = prot_c;
      *(pI++) = to_amino_acid_off(prot_c);
    }
    else {
      *(p++)  = 'x';
      *(pI++) = 20;
    }
  }
  *p = 0;
  *pI = 21;
  if (debug >= 3) {
    fprintf(stderr,"len-seq=%d max=%lu p=%ld\n",(int) strlen(seq),max,p-pseq);
  }
}

#define MAX_FUNC_OI_INDEX 1000000
#define MAX_FUNC_OI_VALS  100000000

char **KmerGuts::load_indexed_ar(const char *filename,int *sz) {
    char **index_ar = (char **) malloc(MAX_FUNC_OI_INDEX * sizeof(char *));
    if (index_ar == 0) {
	std::cerr << "Malloc failed at " __FILE__ " " << __LINE__ << "\n";
	exit(1);
    }
    char *vals      = (char *)malloc(MAX_FUNC_OI_VALS);
    if (vals == 0) {
	std::cerr << "Malloc failed at " __FILE__ " " << __LINE__ << "\n";
	exit(1);
    }
    to_be_freed.push_back(vals);
  char *p         = vals;
  FILE *ifp      = fopen(filename,"r");
  if (ifp == NULL) { 
    fprintf(stderr,"could not open %s\n",filename);
    exit(1);
  }

  *sz = 0;
  int j;
  while ((fscanf(ifp,"%d\t",&j) == 1) && fgets(p,1000,ifp)) {
    if (*sz != j) {
      fprintf(stderr,"Your index must be dense and in order (see line %ld, should be %d)\n",p-vals,*sz);
      exit(1);
    }
    /* fprintf(stderr,"%d is %s\n",*sz,p); */
    index_ar[*sz] = p;               /* the fgets leaves the \n at the end of each line */
    p += strlen(index_ar[*sz]) -1;
    *(p++) = '\0';
    if ((*sz >= MAX_FUNC_OI_INDEX) || ((p-vals) > (MAX_FUNC_OI_VALS - 1000))) {
      fprintf(stderr,"Your function or oI index arrays are too small; bump MAX_FUNC_OI_INDEX and MAX_FUNC_OI_VALS\n");
      exit(1);
    }

    *sz += 1;
  }
  fclose(ifp);
  return index_ar;
}

char **KmerGuts::load_functions(const char *file, int *sz) {
  return load_indexed_ar(file,sz);
}

char **KmerGuts::load_otus(const char *file, int *sz) {
  return load_indexed_ar(file,sz);
}

unsigned long long KmerGuts::lookup_hash_entry(sig_kmer_t sig_kmers[],unsigned long long encodedK) {
    unsigned long long  hash_entry = encodedK % size_hash;
    // std::cerr << "lookup " << hash_entry << " " << size_hash << "\n";
    while ((sig_kmers[hash_entry].which_kmer != encodedK) && (sig_kmers[hash_entry].which_kmer <= MAX_ENCODED)) {
      hash_entry = (hash_entry+1)%size_hash;
/*
      hash_entry++;
      if (hash_entry == size_hash)
	hash_entry = 0;
*/
    }
    if (sig_kmers[hash_entry].which_kmer > MAX_ENCODED) {
      return ULLONG_MAX;
    }
    else {
      return hash_entry;
    }
}

kmer_memory_image_t *KmerGuts::load_raw_kmers(char *file,unsigned long long num_entries, unsigned long long *alloc_sz) {
  /*
   * Allocate enough memory to hold the kmer_memory_image_t header plus the hash table itself.
   */
  *alloc_sz = sizeof(kmer_memory_image_t) + (sizeof(sig_kmer_t) * num_entries);

  kmer_memory_image_t *image = (kmer_memory_image_t *) malloc(*alloc_sz);
  if (image == 0)
  {
      std::cerr << "Malloc failed at " __FILE__ << ":" << __LINE__ << "\n";
      exit(1);
  }

  /*
   * Initialize our table pointer to the first byte after the header.
   */
  image->num_sigs = num_entries;
  image->entry_size = sizeof(sig_kmer_t);
  image->version = (long long) KMER_VERSION;

  sig_kmer_t *sig_kmers = (sig_kmer_t *) (image + 1);

  FILE *ifp      = fopen(file,"r");
  if (ifp == NULL) { 
    fprintf(stderr,"could not open %s",file);
    exit(1);
  }

  for (unsigned long long i = 0; (i < size_hash); i++)
    sig_kmers[i].which_kmer = MAX_ENCODED + 1;

  char kmer_string[32] { 0 };
  unsigned short end_off;
  FunctionIndex fI;
  float f_wt;
  OTUIndex oI;
  unsigned long long loaded = 0;
  while (fscanf(ifp,"%31s\t%hu\t%hu\t%f\t%hu",
		kmer_string,&end_off,&fI,&f_wt,&oI) >= 4) {
      if (strlen(kmer_string) != KMER_SIZE)
      {
	  fprintf(stderr, "Invalid kmer read (size=%ld != %d)\n", strlen(kmer_string), KMER_SIZE);
	  exit(1);
      }

      unsigned long long encodedK = encoded_aa_kmer(kmer_string);
      unsigned long long hash_entry = find_empty_hash_entry(sig_kmers,encodedK);
      loaded++;
      if (loaded >= (size_hash / 2)) {
	  fprintf(stderr,"Your Kmer hash is half-full; use -s (and -w) to bump it\n");
	  exit(1);
      }
      sig_kmers[hash_entry].which_kmer     = encodedK;
      sig_kmers[hash_entry].avg_from_end   = end_off;
      sig_kmers[hash_entry].function_index = fI;
      sig_kmers[hash_entry].otu_index      = oI;
      sig_kmers[hash_entry].function_wt    = f_wt;
  }
  if (debug >= 2)
      fprintf(stderr,"loaded %lld kmers\n",loaded);
  
  fclose(ifp);
  return image;
}

std::shared_ptr<KmerGuts::kmer_handle> KmerGuts::init_kmers(const char *dataD)
{
    std::shared_ptr<KmerGuts::kmer_handle> handle = std::make_shared<KmerGuts::kmer_handle>();
    
    kmer_memory_image_t *image = image_->image();
    
    char file[300];
    strcpy(file,dataD);
    strcat(file,"/function.index");
    handle->function_array = load_functions(file, &handle->function_count);

     strcpy(file,dataD);
     strcat(file,"/otu.index");
     handle->otu_array      = load_otus(file, &handle->otu_count);

     size_hash = image->num_sigs;
     handle->num_sigs = size_hash;
     handle->kmer_table = (sig_kmer_t *) (image + 1);
     // std::cerr << "loaded; size_hash=" << size_hash << "\n";
     return handle;
}


void KmerGuts::advance_past_ambig(unsigned char **p,unsigned char *bound) {

   if (KMER_SIZE == 5) {
     while (((*p) < bound) &&
	    ((*(*p) == 20)     || 
	     (*((*p)+1) == 20) || 
	     (*((*p)+2) == 20) || 
	     (*((*p)+3) == 20) || 
	     (*((*p)+4) == 20) )) {
       (*p)++;
     }
   }
   else {   /*  ##### ASSUMING KMER_SIZE == 8 #### */
     int bad = 1;
     while ((*p < bound) && (bad == 1)) {
       bad = 0;
       if      (*((*p)+7) == 20) {
	 bad = 1;
	 (*p) += 8;
       }
       else if (*((*p)+6) == 20) {
	 bad = 1;
	 (*p) += 7;
       }
       else if (*((*p)+5) == 20) {
	 bad = 1;
	 (*p) += 6;
       }
       else if (*((*p)+4) == 20) {
	 bad = 1;
	 (*p) += 5;
       }
       else if (*((*p)+3) == 20) {
	 bad = 1;
	 (*p) += 4;
       }
       else if (*((*p)+2) == 20) {
	 bad = 1;
	 (*p) += 3;
       }
       else if (*((*p)+1) == 20) {
	 bad = 1;
	 (*p) += 2;
       }
       else if (*((*p)+0) == 20) {
	 bad = 1;
	 (*p) += 1;
       }
     } 
   }
 }

 void KmerGuts::process_set_of_hits(std::shared_ptr<std::vector<KmerCall>> calls,
				    std::shared_ptr<KmerOtuStats> otu_stats)
 {
     if (!otu_stats && !calls)
	 return;

     int fI_count = 0;
     float weighted_hits = 0;
     int last_hit=0;
     for (int i = 0; i < num_hits; i++)
     {
	 if (hits[i].fI == current_fI) {
	     last_hit = i;
	     fI_count++;
	     weighted_hits += hits[i].function_wt;
	 }
     }
     if ((fI_count >= min_hits) && (weighted_hits >= min_weighted_hits))
     {
	 if (calls)
	     calls->push_back({ hits[0].from0_in_prot, hits[last_hit].from0_in_prot+(KMER_SIZE-1), fI_count, current_fI, weighted_hits });

	 /* once we have decided to call a region, we take the kmers for fI and
	    add them to the counts maintained to assign an OTU to the sequence */

	 if (otu_stats)
	 {
	     for (int i=0; (i <= last_hit); i++)
	     {
		 if (hits[i].fI == current_fI)
		 {
		     otu_stats->otu_map[hits[i].oI]++;
		 }
	     }
	 }
     }

     if ((hits[num_hits-2].fI != current_fI) && (hits[num_hits-2].fI == hits[num_hits-1].fI)) {
	 current_fI = hits[num_hits-1].fI;
	 hits[0] = hits[num_hits - 2];
	 hits[1] = hits[num_hits - 1];
	 num_hits                 = 2;
     }
     else {
	 num_hits = 0;
     }
 }

void KmerGuts::gather_hits(const char *pseq,
			    unsigned char *pIseq,
			    std::shared_ptr<std::vector<KmerCall>> calls,
			    std::function<void(const hit_in_sequence_t &)> hit_cb,
			    std::shared_ptr<KmerOtuStats> otu_stats)
 {
     unsigned char *p = pIseq;
     /* pseq and pIseq are the same length */

     unsigned char *bound = pIseq + strlen(pseq) - KMER_SIZE;
     advance_past_ambig(&p,bound);
     unsigned long long encodedK=0;
     if (p < bound) {
	 encodedK = encoded_kmer(p);
     }
     while (p < bound) {
	 unsigned long long  where = lookup_hash_entry(kmersH->kmer_table,encodedK);

	 // pLoc is the offset of this kmer in the protein
	 
	 unsigned int pLoc = (unsigned int) (p - pIseq);
	 // std::cerr << p << " " << encodedK << " " << where << "\n";
	 //for (int i = 0; i < KMER_SIZE; i++)
	 //    std::cerr << prot_alpha[p[i]];
	 //std::cerr << " " << encodedK << " " << where << "\n";
	 if (where != ULLONG_MAX) {
	     sig_kmer_t *kmers_hash_entry = &(kmersH->kmer_table[where]);
	     unsigned short avg_off_end = kmers_hash_entry->avg_from_end;
	     FunctionIndex fI        = kmers_hash_entry->function_index;
	     OTUIndex oI          = kmers_hash_entry->otu_index;
	     float f_wt      = kmers_hash_entry->function_wt;
	     if (hit_cb != nullptr)
		 hit_cb(hit_in_sequence_t(*kmers_hash_entry, pLoc));

	     /*
	      * If we are in a run of hits, and the gap between this hit and the start
	      * of the is greater than our allowed gap, process the hits.
	      */
	     if ((num_hits > 0) && (hits[num_hits-1].from0_in_prot + (unsigned int) max_gap) < pLoc)
	     {
		 if (num_hits >= min_hits)
		 {
		     process_set_of_hits(calls, otu_stats);
		 }
		 else
		 {
		     num_hits = 0;
		 }
	     }

	     if (num_hits == 0)
	     {
		 current_fI = fI;   /* if this is the first, set the current_fI */
	     }

	     if ((! order_constraint) || (num_hits == 0) ||
		 ((fI == hits[num_hits-1].fI) &&
		  (labs((pLoc - hits[num_hits-1].from0_in_prot) -
		       (hits[num_hits-1].avg_off_from_end - avg_off_end)
		      ) <= 20)))
	     {
		 /* we have a new hit, so we add it to the global set of hits */
		 hits[num_hits].oI = oI;
		 hits[num_hits].fI = fI;
		 hits[num_hits].from0_in_prot = pLoc;
		 hits[num_hits].avg_off_from_end = avg_off_end;
		 hits[num_hits].function_wt = f_wt;
		 if (num_hits < MAX_HITS_PER_SEQ - 2)
		     num_hits++;
		 if ((num_hits > 1) && (current_fI != fI) &&           /* if we have a pair of new fIs, it is time to */
		     (hits[num_hits-2].fI == hits[num_hits-1].fI))    /* process one set and initialize the next */
		 {
		     process_set_of_hits(calls, otu_stats);
		 }
	     }
	 }
	 p++;
	 if (p < bound) {
	     if (*(p+KMER_SIZE-1) < 20) {
		 encodedK = ((encodedK % CORE) * 20L) + *(p+KMER_SIZE-1);
	     }
	     else {
		 p += KMER_SIZE;
		 advance_past_ambig(&p,bound);
		 if (p < bound) {
		     encodedK = encoded_kmer(p);
		 }
	     }
	 }
     }
     if (num_hits >= min_hits) {
	 process_set_of_hits(calls, otu_stats);
     }
     num_hits = 0;
 }

 void KmerGuts::process_aa_seq_hits(const std::string &id, const std::string &seq,
			       std::shared_ptr<std::vector<KmerCall>> calls,
			       std::shared_ptr<std::vector<hit_in_sequence_t>> hits,
			       std::shared_ptr<KmerOtuStats> otu_stats)
 {
     auto cb = [this, hits](const hit_in_sequence_t &k) { hits->push_back(k); };
     process_aa_seq(id, seq, calls, cb, otu_stats);
 }

 void KmerGuts::process_aa_seq(const std::string &idstr, const std::string &seqstr,
			       std::shared_ptr<std::vector<KmerCall>> calls,
			       std::function<void(const hit_in_sequence_t &)> hit_cb,
			       std::shared_ptr<KmerOtuStats> otu_stats)
 {
     const char *id = idstr.c_str();
     const char *pseq = seqstr.c_str();
     size_t ln = seqstr.size();

     strcpy(current_id,id);
     current_length_contig = ln;
     current_strand        = '+';
     current_prot_off      = 0;
     for (size_t i=0; (i < ln); i++)
	 pIseq[i] = to_amino_acid_off(*(pseq+i));

     // std::cerr << "'" << id << "' '" << pseq << "' " << ln << "\n";
     gather_hits(pseq,pIseq, calls, hit_cb, otu_stats);
     if (otu_stats)
	 otu_stats->finalize();
 }

 void KmerGuts::process_seq(const char *id,const char *data,
			    std::shared_ptr<std::vector<KmerCall>> calls,
			    std::function<void(const hit_in_sequence_t &)> hit_cb,
			    std::shared_ptr<KmerOtuStats> otu_stats)

 {
     strcpy(current_id,id);
     size_t ln = strlen(data);
     current_length_contig = ln;
     unsigned short i;
     for (i=0; (i < 3); i++)
     {
	 translate(data,i,pseq,pIseq);
	 current_strand   = '+';
	 current_prot_off = i;
	 gather_hits(pseq,pIseq, calls, hit_cb, otu_stats);
     }
     rev_comp(data,cdata);
     for (i=0; (i < 3); i++) {
	 translate(cdata,i,pseq,pIseq);

	 current_strand   = '-';
	 current_prot_off = i;
	 gather_hits(pseq,pIseq, calls, hit_cb, otu_stats);
     }
     if (otu_stats)
	 otu_stats->finalize();
 }

 std::string KmerGuts::format_call(const KmerCall &c)
 {
     std::ostringstream oss;
     oss << "CALL\t" << c.start << "\t" << c.end << "\t" << c.count;
     oss << "\t" << c.function_index << "\t" << function_at_index(c.function_index);
     oss << "\t" << c.weighted_hits << "\n";

     return oss.str();
 }

 std::string KmerGuts::format_hit(const hit_in_sequence_t &h)
 {
     std::ostringstream oss;

     char dc[KMER_SIZE + 1];
     decoded_kmer(h.hit.which_kmer, dc);

     oss << "HIT\t" << h.offset << "\t" << dc << "\t" << h.hit.avg_from_end << "\t" << function_at_index(h.hit.function_index) << "\t" << h.hit.function_wt << "\t" << h.hit.otu_index << "\n";
    
    return oss.str();
}

std::string KmerGuts::format_otu_stats(const std::string &id, size_t size, KmerOtuStats &otu_stats)
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

template<class T>
struct less_second : std::binary_function<T,T,bool>
{
    inline bool operator()( const T& lhs, const T& rhs )
	{
	    return rhs.second < lhs.second;
	}
};  

struct FuncScore
{
    int count;
    float weighted;
    FuncScore() : count(0), weighted(0.0) {}
    FuncScore(int c, float w) : count(c), weighted(w) {}
    FuncScore(const FuncScore &f) : count(f.count), weighted(f.weighted) {}
    FuncScore& operator=(const FuncScore& other)
        {
	    // check for self-assignment
	    if(&other == this)
		return *this;
	    count = other.count;
	    weighted = other.weighted;
	    return *this;
	}
};

/*
 * Find the best call from this set of calls.
 *
 * This code replicates the amino acid version of the SEED pipeline
 * km_process_hits_to_regions | km_pick_best_hit_in_peg
 */
void KmerGuts::find_best_call(std::vector<KmerCall> &calls, FunctionIndex &function_index, std::string &function, float &score, float &weighted_score, float &score_offset)
{
    function_index = UndefinedFunction;
    function = "";
    score = 0.0;
    weighted_score = 0.0;

    if (calls.size() == 0)
    {
	return;
    }
    
    /*
     * First merge adjacent hits that have the same function.
     */
    std::vector<KmerCall> collapsed;

    auto comp = calls.begin();

    while (comp != calls.end())
    {
	collapsed.push_back(*comp);
	comp++;
	KmerCall &cur = collapsed.back();

	while (comp != calls.end() && cur.function_index == comp->function_index)
	{
	    cur.end = comp->end;
	    cur.count += comp->count;
	    cur.weighted_hits += comp->weighted_hits;
	    comp++;
	}
    }
    #if 0
    std::cout << "after collapse:\n";
    for (auto iter = collapsed.begin(); iter != collapsed.end(); iter++)
    {
	std::cout << format_call(*iter);
    }
    #endif

    /*
     *
     * Merge hits when we have a case with
     * +------+--+-------+--+-----------+
     * |  F1  |  |   F2  |  |   F1      |
     * +------+--+-------+--+-----------+
     *
     * where the score for F2 is below 5 and the combined scores for
     * the two F1 hits is 10 or more.
     *
     * If that is the case we discard the F2 hit and combine
     * the F1 hits.
     */

    std::vector<KmerCall> merged;

    int merge_interior_thresh = 5;
    int merge_exterior_thresh = 10;

    comp = collapsed.begin();
    while (comp != collapsed.end())
    {
	merged.push_back(*comp);
	comp++;
	auto comp2 = comp + 1;
	KmerCall &cur = merged.back();
	while (comp != collapsed.end() && comp2 != collapsed.end() &&
	       cur.function_index == comp2->function_index &&
	       comp->count < merge_interior_thresh &&
	       (cur.count + comp2->count) >= merge_exterior_thresh)
	{
	    cur.end = comp2->end;
	    cur.count += comp2->count;
	    cur.weighted_hits += comp2->weighted_hits;
	    comp += 2;
	    comp2 = comp + 1;
	}
    }
   
    #if 0
    std::cerr << "after merge:\n";
    for (auto iter = merged.begin(); iter != merged.end(); iter++)
    {
	std::cerr << format_call(*iter);
    }
    #endif

    /*
     * Determine best call.
     *
     * In the current perl kmer_search (km_pick_best_hit_in_peg) code we just take the best
     * function in terms of weighted score. However, that allows tied scores to be called
     * arbitrarily, and close-to-tied scores to be settled based on insignificant differences.
     *
     * In the original kmerv1 code, we required a score threshold difference (typically 5 hits)
     * between the best function and next best function. We resurrect that here.
     *
     */

    typedef std::map<int, FuncScore> map_t;

    map_t by_func;

    for (auto c: merged)
    {
	auto it = by_func.find(c.function_index);
	if (it == by_func.end())
	{
	    by_func.insert(std::make_pair(c.function_index, FuncScore(c.count, c.weighted_hits)));
	}
	else
	{
	    it->second.count += c.count;
	    it->second.weighted += c.weighted_hits;
	}
    }

    //typedef map_t::value_type ent_t;
    typedef std::pair<FunctionIndex, FuncScore> ent_t;

    std::vector<ent_t> vec;
    for (auto it = by_func.begin(); it != by_func.end(); it++)
	vec.push_back(*it);

//    std::cerr << "vec len " << vec.size() << "\n";
    if (vec.size() > 1)
    {
	std::partial_sort(vec.begin(), vec.begin() +  2, vec.end(), 
			  [](const ent_t& s1, const ent_t& s2) {
			      return s1.second.weighted > s2.second.weighted; });
    }
    
    #if 0
    for (auto x: vec)
    {
	std::cerr << x.first << " " << x.second.count << " " << x.second.weighted << " ";
	std::cerr << function_at_index(x.first) << "\n";
    }
    #endif
    
    if (vec.size() == 1)
	score_offset = (float) vec[0].second.count;
    else
	score_offset = (float) (vec[0].second.count - vec[1].second.count);
    
    // std::cerr << "Offset=" << score_offset << "\n";

    if (score_offset >= 5.0)
    {
	auto best = vec[0];
	function_index = best.first;
	function = function_at_index(function_index);
	score = (float) best.second.count;
	weighted_score = best.second.weighted;
    }
    else
    {
	function_index = UndefinedFunction;
	function = "";
	score = 0.0;

	/*
	 * Try to compute a fallback function naming the two best hits if there are two hits within the
	 * threshold but greater than the next hit.
	 */
	if (vec.size() >= 2)
	{
	    std::string f1 = function_at_index(vec[0].first);
	    std::string f2 = function_at_index(vec[1].first);
	    if (f2 > f1)
		std::swap(f1, f2);

	    if (vec.size() == 2)
	    {
		function = f1 + " ?? " + f2;
		score = (float) vec[0].second.count;
	    }
	    else if (vec.size() > 2)
	    {
		float pair_offset = (float) (vec[1].second.count - vec[2].second.count);
		if (pair_offset > 5.0)
		{
		    function = f1 + " ?? " + f2;
		    score = (float) vec[0].second.count;
		    score_offset = pair_offset;
		    weighted_score = vec[0].second.weighted;
		}
	    }
	}
    }
}
