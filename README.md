# PATRIC family computation

 1. Create base data directory with pf-construct-genus-data-for-genera
    This uses pf-build-nr to compute nonredundant databases, and pf-compute-nr-seqs
    to create sets of sequences corresponding to the reference sequence in the nonredundant
    sets.

 2. Use pf-submit-all-local families to start the local fam run.


# Some perf notes

olson@bio-compute-02:/vol/bvbrc-families/fams/test.2024-1206/genus.data$ rm -r /dev/shm/musc; mkdir -p /dev/shm/musc; time pf-compute-local-family-alignments --parallel 30 Buchnera-32199 ../kmers /dev/shm/musc Buchnera-32199  --overload 1 --aligner mafft

real	 7m10.889s
user	 30m36.985s
sys	 161m26.406s

muscle

real	11m35.496s
user	110m13.439s
sys	150m33.937s

muscle 60 proc:
real   8m13.195s
user   130m38.191s
sys    173m23.799s