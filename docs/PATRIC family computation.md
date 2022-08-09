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

