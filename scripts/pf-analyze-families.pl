#
# Given local families, pull the annotations from BV-BRC and compute the most common
# annotation there and emit a table
#
#  fam-id size fam-function brc-function
#

use strict;
use Data::Dumper;
use P3DataAPI;

my $api = P3DataAPI->new();

my $fam_defs = "local.family.defs";
my $fam_members = "local.family.members";

my %fams;
my @fams;

open(D, "<", $fam_defs) or die "cannot read $fam_defs: $!";
while (<D>)
{
    chomp;
    my($id, $function, $sub, $mean, $dev, $gene, $size) = split(/\t/);
    next if $size < 4;
    push(@fams, $id);
    $fams{$id} = { id => $id, function => $function, size => int($size) };
}
close(D);

my %all_pegs;

open(M, "<", $fam_members) or die "Cannot read $fam_members: $!";
while (<M>)
{
    chomp;
    my($id, $peg, $len, $zsc, $gene) = split(/\t/);
    my $f = $fams{$id};
    next unless $f;
    $all_pegs{$peg}++;
    push(@{$f->{members}}, { peg => $peg });
}

my $fns = $api->function_of([keys %all_pegs]);

for my $id (@fams)
{
    my $fam = $fams{$id};
    my @pegs = map { $_->{peg} } @{$fam->{members}};
    my @fns = map { $fns->{$_} } @pegs;
    my %fn_count;
    $fn_count{$_}++ foreach @fns;
    my @by_count = sort { $b->{count} <=> $a->{count} } map { { function => $_, count => $fn_count{$_} } } keys %fn_count;
    print join("\t",
	       $id, $fam->{size}, $fam->{function},
	       $by_count[0]->{function}, $by_count[0]->{count},
	       $by_count[1]->{function}, $by_count[1]->{count},
	       $fam->{members}->[0]->{peg},
	      ), "\n";
}
