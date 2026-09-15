#
# Based on peg id, find the members of the given family in the new fams file and print the
# matching entries from the old fams file
#

use strict;
use Data::Dumper;

@ARGV == 3 or die "Usage: $0 family-id new-fams old-fams\n";

my $target_fam = shift;
my $new_file = shift;
my $old_file = shift;

my %fids;

open(F, "<", $new_file) or die "Cannot open $new_file: $!";
open(O, "<", $old_file) or die "Cannot open $old_file: $!";


while (<F>)
{
    my($fam, $fid) = /^(\S+)\t(\S+)/;
    last if $fam ne $target_fam && %fids;
    $fids{$fid} = 1 if $fam eq $target_fam;
}

printf STDERR "%d members\n", scalar keys %fids;

my %fams;
my %fids2;

while (<O>)
{
    my(@a) = split(/\t/);
    if ($fids{$a[3]})
    {
	print;
	$fams{$a[0]} = 1;
    }
}
print "---------\n";
seek(O, 0, 0) or die "seek failed: $!";
print STDERR "Search fams matching fids\n";
while (<O>)
{
    my(@a) = split(/\t/);
    if ($fams{$a[0]})
    {
	print;
	$fids2{$a[3]} = 1;
    }
}
print "---------\n";

seek(F, 0, 0) or die "seek failed: $!";
print STDERR "mapping back to new\n";

while (<F>)
{
    my($fid) = /^\S+\t(\S+)/;
    print if ($fids2{$fid});
}
