#
# Given a set of local families, pull the proteins for them from the NR file
# in merge-dir/prop-dir/nr/genus.fa
#

use strict;
use gjoseqlib;
use Data::Dumper;
use File::Basename;

@ARGV >= 3 or die "Usage: $0 nr-file out-dir family [family...]\n";

my $nr_file = shift;
my $out_dir = shift;
my @fams = @ARGV;

my %fams = map { $_ => 1 } @fams;

my $genus = basename($nr_file, ".fa");

open(FA, "<", $nr_file) or die "Cannot read $nr_file: $!";

my %out_fh;
for my $fam (@fams)
{
    my $file = sprintf("$out_dir/%08d.fa", $fam);
    open($out_fh{$fam}, ">", $file) or die "Cannot write $file: $!";
}

while (my($id, $def, $seq) = read_next_fasta_seq(\*FA))
{
    my($pgf, $plf, $g) = split(/\s+/, $def);
    if ($fams{$plf})
    {
	print "found $id\n";
	write_fasta($out_fh{$plf}, [$id, $def, $seq]);
    }
}
