#
# Given two id-clust files as output from svr_corresponding_sequences runs, find
# clusters that overlap.
#
# Pull the protein data from the given fasta file and rerun svr_corresponding_sequences
# on the cluster.
#

use strict;
use Data::Dumper;
use File::Path qw(make_path);
use File::Slurp;

@ARGV == 6 or die "Usage: $0 deleted.fids source-proteins1.fa source-proteins2.fa clust-1 clust-2 output-dir\n";

my $deleted_ids = shift;
my $src1_fa = shift;
my $src2_fa = shift;
my $clust1 = shift;
my $clust2 = shift;
my $out = shift;

my %del;
if (my @del = read_file($deleted_ids))
{
    chomp @del;
    $del{$_} = 1 foreach @del;
}

open(SRC1, "<", $src1_fa) or die "Cannot open $src1_fa: $!";
open(SRC2, "<", $src2_fa) or die "Cannot open $src2_fa: $!";
open(C1, "<", $clust1) or die "Cannot open $clust1: $!";
open(C2, "<", $clust2) or die "Cannot open $clust2: $!";


#-d $out and die "Output dir $out already exists\n";
make_path($out);

#
# Read clust1 into a hash whose keys are peg ids and values are the cluster number
#

my $fam1 = 1;
my %c1;
while (<C1>)
{
    chomp;
    my(@pegs) = grep { !$del{$_} } split(/\t/);
    $c1{$_} = $fam1 foreach @pegs;
    $fam1++;
}

die Dumper(\%c1);
