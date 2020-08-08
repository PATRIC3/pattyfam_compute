#
# This is a SAS Component
#


=head1 NAME

    pf-inflate-families

=head1 SYNOPSIS

    pf-inflate-families fam-dir fams-file id-column gene-name-column

=head1 DESCRIPTION

This is a pattyfam computation support routine.

Given a families file that was computed using the nonredundant protein set, 
inflate it to include all the identical proteins from the original data. This
is done using the nr-refs file computed during the early part of the local
family computation process.

=cut
    
use strict;
use Data::Dumper;
use Getopt::Long;

use Getopt::Long::Descriptive;

my($opt, $usage) = describe_options("%c %o fam-dir fam-file id-column gene-name-column",
				    ["output|o=s" => "Write output here"],
				    ["help|h" => "Show this help message."]);

print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 4;

my $fam_dir = shift;
my $fam_file = shift;
my $id_col = shift;
my $gn_col = shift;
$id_col--;
$gn_col--;

my $out_fh;
if ($opt->output)
{
    open($out_fh, ">", $opt->output) or die "Cannot write " . $opt->output . ": $!\n";
}
else
{
    $out_fh = \*STDOUT;
}

#
# Read nr-refs file to find the proteins to be mapped.
#
my %expand;
open(NRR, "<", "$fam_dir/nr-refs") or die "Cannot read $fam_dir/nr-refs: $!";
while (<NRR>)
{
    chomp;
    my($key, @vals) = split(/\t/);
    $expand{$key} = [@vals];
}
close(NRR);

#
# Read families and inflate.
#
# We first read the gene.names file so that when we inflate, we can add gene names to the inflated data.
#

my %gene_name;
open(GN, "<", "$fam_dir/gene.names") or die "Cannot read $fam_dir/gene.names: $!";
while (<GN>)
{
    chomp;
    my($id, $name) = split(/\t/);
    $gene_name{$id} = $name;
}
close(GN);

open(FAMS, "<", $fam_file) or die "Cannot read $fam_file: $!";

while (<FAMS>)
{
    print $out_fh $_;
    chomp;
    my @cols = split(/\t/);
    my $id = $cols[$id_col];
    my $rest = $expand{$id};
    my $gn = $gene_name{$id};
    
    for my $syn (@$rest)
    {
	$cols[$id_col] = $syn;
	$cols[$gn_col] = $gn;
	print $out_fh join("\t", @cols), "\n";
    }
}

close(FAMS);
