
=head1 NAME

    pf-compute-kmer-clusters

=head1 SYNOPSIS

    pf-compute-kmer-clusters work-dir gene-names nr-refs output-families unclassified-proteins-file

=head1 DESCRIPTION

This is a pattyfam computation support routine.

Compute the local family clusters using MCL.

We split the computation into two groups. Families smaller in size than our 
cutoff threshold are grouped together and run in parallel batches with sequential MCL.

Families larger than that are run sequentially with multithreaded MCL. My experiments
have shown we get speedups of around 4 with 8 processors which makes it worth it for the 
large families (which tend to otherwise create a long tail of computation when run
with the others in the batches).

When the MCL pass is complete, we write out the families in the standard format.

We also write a file of proteins that were not accounted for in this process;  these
will be sent to a BLAST-based clustering tool.

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

my($opt, $usage) = describe_options("%c %o work-dir gene-names nr-refs output-families unclassified-proteins-file",
				    ["inflation|I=s" => "MCL inflation", { default => 3.0 }],
				    ["parallel=i" => "Parallel threads", { default => 1 }],
				    ["large-family-cutoff=i" => "Cutoff for large families (to switch to parallel MCL)",
				     { default => => 1_000_000_000 }],
				    ["mcl=s" => "MCL executable", { default => "mcl" }],
				    ["help|h" => "Show this help message."]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 5;


my $work_dir = shift;
my $gene_names = shift;
my $nr_refs = shift;
my $output_families = shift;
my $unclassified_file = shift;

-d $work_dir or die "Work directory $work_dir is missing\n";

my $map_file = "$work_dir/map";
my $fasta_dir = "$work_dir/fasta";
my $dist_dir = "$work_dir/dist";
my $mcl_dir = "$work_dir/mcl";

my $n_buckets = $opt->parallel * 10;
my $scheduler_small = LPTScheduler->new($n_buckets);
my $scheduler_large = LPTScheduler->new($n_buckets / 4);

#
# Partition the families that have non-zero counts in the distances directory.
#

my @large_families;
opendir(D, $dist_dir) or die "Cannot opendir $dist_dir: $!";
while (my $f = readdir(D))
{
    next if $f =~ /^\./;
    my $path = "$dist_dir/$f";
    my $mcl_path = "$mcl_dir/$f";

    # next if (-s "$mcl_path.mcl.out");
    next unless -f $path;
    my $sz = -s _;
    if ($sz == 0)
    {
	open(F, ">", $mcl_path) or die "Cannot write $mcl_path: $!";
	close(F);
	next;
    }
    else
    {
	my $work = [$f, $path, $sz, "$mcl_dir/$f"];
	if ($sz < $opt->large_family_cutoff)
	{
	    $scheduler_small->add_work($work, $sz);
	}
	else
	{
	    $scheduler_large->add_work($work, $sz);
	}
    }
}

#print Dumper($scheduler_small, $scheduler_large);

#
# Run the small families.
#
$scheduler_small->run(sub {}, sub {
    my($global, $ent) = @_;
    process_family(@$ent, 1);
}, $opt->parallel);

#
# Run the large families
#
$scheduler_large->run(sub {}, sub {
    my($global, $ent) = @_;
    process_family(@$ent, 4);
}, int($opt->parallel / 4));


#
# Write families in standard format
#

my $next_idx = 1;

#
# Read the family map to get functions.
#

open(F, "<", $map_file) or die $!;
my %fun;
my @fams;
while (<F>)
{
    chomp;
    my($id, $fun) = split(/\t/);
    $fun{$id} = $fun;
    push(@fams, $id);
}
close(F);

#
# Read the hypo.fams file to get the list of all hypothetical non-called IDs, so we can
# clear out IDs that we find in the kmer clusters and generated the
# unclassified file for doing non-kmer-based processing.
#

my %hypo;
open(HP, "<", "$work_dir/hypo.prots") or die "Cannot read $work_dir/hypo.prots: $!";
while (<HP>)
{
    if (/^(\S+)/)
    {
	$hypo{$1} = 1;
    }
}
close(HP);

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

my $out_fh;
open($out_fh, ">", $output_families) or die "Cannot open $output_families: $!";

for my $fam (@fams)
{
    if (!open(F, "<", "$mcl_dir/$fam"))
    {
	warn "cannot open $mcl_dir/$fam: $!";
	next;
    }

    if (-s F > 0)
    {
	while (<F>)
	{
	    chomp;
	    my @pegs = split(/\t/);

	    my $idx = $next_idx++;

	    my %names;
	    for my $peg (@pegs)
	    {
		my $gn = $gene_names{$peg};
		push(@{$names{$gn}}, $peg) if $gn;

		delete $hypo{$peg};
	    }
	    my $name;
	    my @seen = keys %names;
	    if (@seen == 1)
	    {
		$name = $seen[0];
	    }
	    for my $peg (@pegs)
	    {
		print $out_fh "$fun{$fam}\t$idx\t$peg\t$gene_names{$peg}\t$name\n";
	    }

	}
	close(F);
    }
    else
    {
	warn "Empty family $fam\n";
    }
}

#
# Unfound hypos go to really-missed file
#
open(R, ">", $unclassified_file) or die "Cannot open $unclassified_file: $!";
print R "$_\n" foreach sort keys %hypo;
close(R);

close($out_fh);

sub process_family
{
    my($fam, $dist_path, $dist_size, $mcl_path, $threads) = @_;

    my $in_filter;

    if ($dist_size > 10_000_000_000)
    {
	$in_filter = ['awk', '-v', "OFS=\t", "-F", "\t",
		      ' $3 > 4 { print $1, $2, $4 }', $dist_path]
    }
    else
    {
	$in_filter = ['cut', '-f', '1,2,4', $dist_path];
    }

    my @cmd = ($in_filter,
	       '|',
	       [$opt->mcl, '-', '-t', $threads, '-o', $mcl_path, '--abc', '-I', $opt->inflation],
	       '2>', "$mcl_path.mcl.out");

    my $ok = run(@cmd);
    $ok or die "cmd failed with $?: " . Dumper(\@cmd);
}
