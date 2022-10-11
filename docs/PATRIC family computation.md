# PATRIC family computation

The family computation comprises a number of phases. We maintain a single directory with all the related data; we call this the `data-dir` below. 

A primary component of the data directory is the `genus.data` directory which contains a subdirectory for each genus included in the family build.

## Create base data directory

Create base data directory with `pf-construct-genus-data-for-genera`
This uses pf-build-nr to compute nonredundant databases, and pf-compute-nr-seqs
to create sets of sequences corresponding to the reference sequence in the nonredundant
sets.

We typically create a family release directory (e.g /vol/patric3/fams/patric-fams-2021-0427) and place the genomic data in a subdirectory genus.data; that subdirectory is what is passed to `pf-construct-genus-data-for-genera`.

```
$ pf-construct-genus-data-for-genera
pf-construct-genus-data-for-genera.pl [-bdghMmp] [long options...] data-dir
        --rank STR                       Use the given taxon rank for
                                         grouping (defaults to genus)
        -g STR... --genus STR...         Limit to the given genus. May be
                                         repeated for multiple genera.
        --genomes STR                    Limit to the genomes in this file
        -b STR... --bad-genome STR...    A bad genome. Don't use it.
        -m INT --min-genomes INT         Minimum number of genomes required
                                         in a genus to build families
        -M INT --max-genomes INT         Maximum number of genomes to be
                                         placed in a genus
        -d STR --solr-url STR            Solr API url
        -p INT --parallel INT            Run with this many procs
        --genome-dir STR                 Directory holding PATRIC genome data
        -h --help                        Show this help message
```



## Compute Local Families

Local families are computed on a per-genus basis using the `pf-compute-local-families`script. It uses the following parameters to control the family creation behavior:

Use pf-submit-all-local families to start the local fam run.



### Notes on scripts from original source under FIGdisk

#### p3x-clean-small-families

Prune family files (global or local) to only families above a minimum size.

#### compute_merge_for_family.pl


Given a family.dist file with the kmer distances, a peg.map file,
and an inflation parameter, compute the clustering and generate
the merge definition for the original families.

#### get_merged_families_p

Given a set of family directories, and the corresponding kmer data directory,
compute a merge of the families.

We merge based on function. We assign function numbers based on the function.index
file in the kmer data directory.

For each function, we construct a fasta file containing the representative pegs
from the families for that function, from each family directory.
The ids given to the proteins in the fasta files are of the form
genus-id:fam-id:peg so that join decisions later on can be made strictly on
the results from the clustering on these fasta files.

The number of pegs to use as representatives is a parameter to this program.
The method for selecting pegs is also a parameter.

Given the set of fasta files, we use the kmer distance clustering algorithm to
compute suggested family merges.

Creates the peg.map file which contains columns

  rep_id, genus, family, func-index, function, genome name

Detailed program flow

 * Create list @dirs of pairs [genus-dir, total-sequence-size]
 * Annotate the work list with per-worker output directories for chunked fasta data
 * LPT-schedule in parallel this:
   * For each genus, collect family fasta
 * Merge the chunked data into larger files
 * Create a work list of an item for each function
 * LPT-schedule based on size of fasta:
    * For each function:
      * Start a kser
      * Load fasta data using the /add route
      * Dump the matrix using /matrix

##### Collect family 

  * Load sequences from genus-dir/nr-seqs into %seq
  * Load family definitions into %fam_info
  * For each of the families we are processing (we may keep only familes for a list of desired roles), emit fasta:
    * Pick a set of representatives from the local family
    * For each rep, write an entry into the peg map, per-family fasta dir, and all-peg fasta dir

#### get_merged_families_2

Clip overly large groups to a minimum kmer score

Run mcl on the groups

#### get_merged_families_3

Load the local family data from the genus directories, including lists of peg ids that are in the families.

Load the peg.map (map from rep id to local family info)


## NOTES

Add a p3x-compute-kmer-distance program that encapsulates the start-kser / load / pull matrix functionality

## Alternative decomposition of the problem

The early versions of this pipeline were targeted to run in parallel on a single node. As the size of the collection has grown this has become insufficient, especially as we wish to run the entire process more often.

There are two distinct phases of the problem. In the first, we are operating on a per-genus basis to compute local families within each genus. After that is complete, we switch to computing on a per-function basis to compute the candidates for the merged global families.

We arrange the end of the per-group computation to initialize the per-family computation with the data that is already local. In our shared directory structure we have the following layout:
 
 families/genus.data/GroupName-Taxid              All data for the group
                                    /peg.map      Table mapping rep sequence to metadata
          function.data/Function-id               All data for the given function ID          function.data/Function-id/aa-fasta      Fasta data of family reps for the function


As a result we can define a refactored pipeline via scripts that perform the following operations:

 * Create an initial directory containing a subdirectory for each genome group to be computed on. These correspond to genera in the bacteria and archaea, and probably families in viruses.
 * For each genome group, we have scripts that perform the following operations:
   * Download sequence data (AA & DNA), look up gene names, scan for contig sizes and mark potentially truncated genes. (currently in pf-construct-genus-data-for-genera)
   * Compute NRs for group (AA & DNA)
   * Compute local families and alignments. (pf-compute-local-families - this is already in the proper form)
   * Write the peg.map for the reps of the group
   * write the function.data/Function-id/aa-fasta/Groupname-Taxid file with the sequence data of the reps
 * 


== The pipeline

```
pf-initialize-data -g Buchnera -g 'Candidatus Carsonella' -g Candidatus\ Riesia --check-quality /vol/bvbrc-families/fams/test-02/genus.data

pf-load-group-data /vol/bvbrc-families/fams/test-02/genus.data/Candidatus\ Carsonella-114185

pf-load-group-data /vol/bvbrc-families/fams/test-02/genus.data/Buchnera-32199


```