#ifndef _kmer_params_h
#define _kmer_params_h


#define KMER_SIZE 8
#define MAX_SEQ_LEN 500000000

#if KMER_SIZE == 5 
const CORE = 20L*20L*20L*20L;
#endif
#if KMER_SIZE == 8
const unsigned long long CORE = 20L*20L*20L*20L*20L*20L*20L;
#endif
#if KMER_SIZE == 12
const unsigned long long CORE = 20L*20L*20L*20L*20L*20L*20L*20L*20L*20L*20L;
#endif

#define MAX_ENCODED CORE*20L 

#define MAX_HITS_PER_SEQ 40000

#define OI_BUFSZ 5

#endif // _kmer_params_h
