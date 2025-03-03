#
# Given the output of running kmers-call-functions on the uniprot data, find good
# calls and show mappings for the core_cid hypotheticals.
#

use strict;
use Data::Dumper;

@ARGV == 3 or die "Usage: $0 uniprot-sprot.fasta calls kmer-dir\n";

my $uniprot = shift;
my $calls = shift;
my $kmer_dir = shift;

my %uniprot;
open(U, "<", $uniprot) or die;
while (<U>)
{
    if (/^>(\S+)\s+(.*?)[A-Z]+=/)
    {
	$uniprot{$1} = $2};
}

close(U);

my %hypo;
for my $f (<$kmer_dir/recall.report.d/*>)
{
    open(F, "<", $f) or die $f;
    while (<F>)
    {
	chomp;
	my($fid, $old1, $old2, $new, $fnid, $score) = split(/\t/);
	if ($new =~ /core_cid/)
	{
	    push(@{$hypo{$new}}, [$fid, $score]);
	}
    }
    close(F);
}

open(C, "<", $calls) or die;
while (<C>)
{
    chomp;
    my($id, $fn, $fnid, $score) = split(/\t/);

    if ($fn =~ /core_cid/ && $score > 100)
    {
	my $ufn = $uniprot{$id};
	print join("\t", $id, $score, $fn, $ufn), "\n";
	for my $m (@{$hypo{$fn}})
	{
	    print join("\t", '', @$m), "\n";
	}
    }
}
