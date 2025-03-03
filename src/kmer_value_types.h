#ifndef _kmer_value_types_h
#define _kmer_value_types_h

/**
 * Function indexes are use to reference the entries in function.index
 * which represent function strings assigned to proteins.
 */
typedef unsigned short FunctionIndex;

/**
 * Value representing a missing or undefined function.
 */
const FunctionIndex UndefinedFunction = USHRT_MAX;

/**
 * OTU indexes are use to reference the entries in otu.index
 * which represent OTUs associated with kmers.
 */
typedef unsigned short OTUIndex;

/**
 * Value representing a missing or undefined function.
 */
const FunctionIndex UndefinedOTU = USHRT_MAX;



#endif
