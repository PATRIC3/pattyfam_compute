
=head1 NAME

    pf-compute-blast-clusters

=head1 SYNOPSIS

    pf-compute-blast-clusters work-dir gene-names unclassified-proteins-file blast-clusters

=head1 DESCRIPTION

This is a pattyfam computation support routine.

Compute the local family clusters using blast.


=cut

use strict;
use POSIX;
use Time::HiRes 'gettimeofday';
use File::Slurp;
use Getopt::Long::Descriptive;
use KmerGutsNet;
use gjoseqlib;
use IPC::Run 'run';
use IO::Handle;
use Data::Dumper;
use Proc::ParallelLoop;
use File::Path qw(make_path);
use File::Slurp;
use Carp::Always;
use LPTScheduler;
use SeedUtils;

my($opt, $usage) = describe_options("%c %o work-dir gene-names seqs-dir unclassified-proteins-file output-families",
				    ["identity|i=s" => "Identity", { default => 0.5 }],
				    ["parallel=i" => "Parallel threads", { default => 1 }],
				    ["help|h" => "Show this help message."]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 5;

my $work_dir = shift;
my $gene_names = shift;
my $seqs_dir = shift;
my $unclassified_file = shift;
my $output_families = shift;

-d $work_dir or die "Work directory $work_dir is missing\n";

my $map_file = "$work_dir/map";
my $fasta_dir = "$work_dir/fasta";

my %genomes;
my %needed_pegs;
open(IDS, "<", $unclassified_file) or die "Cannot open $unclassified_file: $!";
while (defined(my $peg = <IDS>))
{
    chomp $peg;
    my $g = &SeedUtils::genome_of($peg);
    $genomes{$g} = 1;
    $needed_pegs{$peg} = 1;
}
close(IDS);

my $seq_data = "$work_dir/unclassified.fa";
open(SEQ, ">", $seq_data) or die "cannot write $seq_data: $!";

for my $g (keys %genomes)
{
    open(S, "<", "$seqs_dir/$g") or die "Cannot open $seqs_dir/$g: $!";
    while (my $ent = read_next_fasta_seq(\*S))
    {
	if (delete $needed_pegs{$ent->[0]})
	{
	    write_fasta(\*SEQ, $ent);
	}
    }
    close(S);
}
close(SEQ);

my $fams_file = "$work_dir/blast.fams";
my @cmd = ("svr_representative_sequences", "-s", $opt->identity, "-b", "-f", $fams_file);
my $ok = run(\@cmd, "<", $seq_data, ">", "/dev/null");
$ok or die "Failed running @cmd: $?\n";

#
# Read the gene names file and try to deduce gene name for local fams.
#
my %gene_names;
if (open(GN, "<", $gene_names))
{
    while (<GN>)
    {
	chomp;
	my($id, $name) = split(/\t/);
	$gene_names{$id} = $name;
    }
    close(GN);
}

open(OUT, ">", $output_families) or die "Cannot write $output_families: $!";
open(B, "<", $fams_file) or die "Cannot open $work_dir/blast.fams: $!";
my $idx = 1;
while (<B>)
{
    chomp;
    my %names;
    my @pegs = split(/\t/);
    for my $peg (@pegs)
    {
	my $gn = $gene_names{$peg};
	push(@{$names{$gn}}, $peg) if $gn;
    }
    my $name;
    my @seen = keys %names;
    if (@seen == 1)
    {
	$name = $seen[0];
    }
    for my $peg (@pegs)
    {
	print OUT "hypothetical protein\t$idx\t$peg\t$gene_names{$peg}\t$name\n";
    }
    $idx++;
}
