#
# Look up the given pegs in the local families in the current dir
#

use Data::Dumper;
use strict;

my %want;
my @pegs;
if (@ARGV)
{
    @pegs = @ARGV;
    $want{$_} = 1 foreach @ARGV;
}
else
{
    while (<STDIN>)
    {
	s/^\s*//;
	s/\s*$//;
	push(@pegs, $_);
	$want{$_} = 1;
    }
}
open(F, "<", "local.family.members.expanded") or die;
my %dat;
while (<F>)
{
    chomp;
    my($fam, $fid, $len, $zs) = split(/\t/);
    if (delete($want{$fid}))
    {
	$dat{$fid} = { fam => $fam, len => $len, zs => $zs};
    }
}

open(F, "<", "local.family.defs") or die;
my %fam;
while (<F>)
{
    chomp;
    my @a = split(/\t/);
    my($fam, $fn, $sf) = @a;
    $fam{$fam} = [@a];
}

for my $peg (@pegs)
{
    my $h = $dat{$peg};
    print join("\t", $peg, $h->{fam}, $fam{$h->{fam}}->[1]), "\n";
}
