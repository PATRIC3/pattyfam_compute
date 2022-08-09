#
# Pull the hypotheticals of a given length range from
# the annotations/0 directory of the given sims build.
#

use File::Slurp;
use strict;
use gjoseqlib;

@ARGV == 4 or die "Usage: $0 kmer-dir deleted-ids min max\n";

my $dir = shift;
my $deleted_ids = shift;
my $min_size = shift;
my $max_size = shift;

my $adir = "$dir/Annotations/0";
my $sdir = "$dir/Seqs";

my %del;
if (my @del = read_file($deleted_ids))
{
    chomp @del;
    $del{$_} = 1 foreach @del;
}
opendir(D, $adir) or die "cannot open $adir: $!";

for my $g (sort { $a <=> $b } grep { /^\d+\.\d+$/ && -f "$adir/$_" } readdir(D))
{
    open(ANNO, "<", "$adir/$g") or die "Cannot open $adir/$g: $!";
    my %want;
    while (<ANNO>)
    {
	if (/^(\S+)\thypothetical protein$/ && !$del{$1})
	{
	    $want{$1} = 1;
	}
    }
    close(ANNO);

    open(SEQ, "<", "$sdir/$g") or die "Cannot open $sdir/$g: $!";
    while (my($id, $def, $seq) = read_next_fasta_seq(\*SEQ))
    {
	next unless $want{$id};
	my $l = length($seq);
	if ($l >= $min_size && $l <= $max_size)
	{
	    print_alignment_as_fasta([$id, $def, $seq]);
	}
    }
    close(SEQ);
}
