#
# test the aligner wrapper
#

use P3AlignmentComputeToDisk;
use strict;
use gjoseqlib;
use Data::Dumper;

my $inp = shift;
my $n = shift;

my @seqs = read_fasta($inp);
@seqs or die "No seqs provided in '$inp'\n";

my $p = P3AlignmentComputeToDisk->new("/tmp", aligner => 'muscle');
for (my $rep = 0; $rep < $n; $rep++)
{
    $p->process_data("genus", "fam", \@seqs, 'aa', 11);
}
