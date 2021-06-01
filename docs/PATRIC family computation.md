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



