#
# Look thru singletons in mmseqs data and find clusters in blast data
#

use strict;

@ARGV == 2 or die "Usage: mmseq.clusters blast.clusters\n";

my $mmseq = shift;
my $blast = shift;

open(M, "<", $mmseq) or die;

my %fam;
my %singleton;
while (<M>)
{
    chomp;
    my @c = grep { /fig/ } split(/\t/);
    $fam{$_} = \@c foreach @c;
    if (@c == 1)
    {
	$singleton{$c[0]} = 1;
    }
}
close(M);

open(B, "<", $blast) or die;

while (<B>)
{
    chomp;
    my @c = grep { /fig/ } split(/\t/);
    next unless @c > 2;

    for my $id (@c)
    {
	if ($singleton{$id})
	{
	    print "---\n";
	    print join("\t", $id, @c), "\n";
	    print join("\t", undef, @{$fam{$_}}) . "\n" foreach @c;
	}
    }
}

    
