
=head1 NAME

    pf-submit-all-local-families

=head1 SYNOPSIS

    pf-submit-all-local-families dir

=head1 DESCRIPTION

This is a pattyfam computation support routine.

Given a top level PATRIC families direcotry (which contains a genus.data directory
with all genus data to be computed), submit, via sbatch, jobs to compute
them all.

We will submit parallel jobs. Genera under the "large genus" cutoff are run with 8 processors, 
while genomes over that cutoff are run with 16.

The measure for large genus is the count of genomes in that genus.


=cut

use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Cwd 'abs_path';
use IPC::Run 'run';

my($opt, $usage) = describe_options("%c %o kmer-dir family-dir",
				    ["large-genus-cutoff=i" => "Cutoff for large genus size, for scheduling", { default => 5000 }],
				    ["small-genus-process-count=i" => "Number of processors for small genera", { default => 8 }],
				    ["large-genus-process-count=i" => "Number of processors for large genera", { default => 16 }],
				    ["small-memory=s" => "Memory required for small genera", { default => "100G" }],
				    ["large-memory=s" => "Memory required for large genera", { default => "500G" }],
				    ["genus-index=s" => "Do not recompute genus data; resubmit the given indices"],
				    ["account=s" => "SLURM account", { default => $ENV{USER} }],
				    ["identity=s" => "Identity for BLAST fallback", { default => 0.5 }],
				    ["inflation=s" => "MCL inflation", { default => 3.0 }],
				    ["good-cutoff=s" => "Fraction of members with unique genomes required to be 'good'", { default => 0.9 }],
				    ["genome-dir=s", "Directory holding PATRIC genome data", { default => "/vol/patric3/downloads/genomes" }],
				    ["kser=s" => "Path to kser executable", { default => "kser" }],
				    ["container=s" => "Container to use", { default => '/vol/patric3/production/containers/bvbrc-dev-258.sif' }],
				    ["dev-container=s" => "Dev container path", { default => "/home/olson/P3/dev-families/dev_container" }],
				    ["just-one" => "Only submit one"],
				    ["help|h" => "Show this help message."]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 2;

my $kmer_dir = abs_path(shift);
my $fam_dir = abs_path(shift);

#
# Count genomes in the genera. Returns a list sorted in reverse by count
# of pairs [genus-name, count].
#

my @genus_counts = count_genera($fam_dir);

#
# If we are running with --genus-index set, read the existing index file.
#

my $index_map;
my $to_compute;
if ($opt->genus_index)
{
    $index_map = read_index_file("$fam_dir/genus.index");
    open(TC, "<", $opt->genus_index) or die "Cannot open " . $opt->genus_index . ": $!";
    $to_compute = [];
    while (<TC>)
    {
	if (/(\d+)/)
	{
	    push(@$to_compute, $1);
	    last if $opt->just_one;
	}
    }
    close(TC);
}
else
{
    #
    # index file maps from job index to genus to process.
    #
    $index_map = write_index_file(\@genus_counts, "$fam_dir/genus.index");
}

#
# Submit jobs. First find the index which denotes teh
# breakpoint between small and large genera. Remember list is sorted
# largest to smallest.
#

my $break = 0;
while ($break < @genus_counts && $genus_counts[$break]->[1] >= $opt->large_genus_cutoff)
{
    $break++;
}

submit_jobs($kmer_dir, $fam_dir, $index_map, \@genus_counts, 0, $break - 1, $opt->large_genus_process_count, $opt->large_memory, $opt, $to_compute);
submit_jobs($kmer_dir, $fam_dir, $index_map, \@genus_counts, $break, $#genus_counts, $opt->small_genus_process_count, $opt->small_memory, $opt, $to_compute);

sub count_genera
{
    my($fam_dir) = @_;

    opendir(D, "$fam_dir/genus.data") or die "Cannot opendir $fam_dir/genus.data: $!";
    my @out;
    while (my $genus = readdir(D))
    {
	if (open(G, "<", "$fam_dir/genus.data/$genus/genomes"))
	{
	    my $n = 0;
	    while (<G>)
	    {
		$n++;
	    }
	    close(G);
	    push @out, [$genus, $n];
	    last if $opt->just_one;
	}
    }
    return sort { $b->[1] <=> $a->[1] } @out;
}

sub write_index_file
{
    my($counts, $index_file) = @_;
    open(I, ">", $index_file) or die "Cannot write index file $index_file: $!";

    my $idx = 1;
    my $map = {};
    for my $ent (@$counts)
    {
	print I "$idx\t$ent->[0]\n";
	$map->{$ent->[0]} = $idx;
	$idx++;
    }
    close(I);
    return $map;
}

sub read_index_file
{
    my($index_file) = @_;
    open(I, "<", $index_file) or die "Cannot read index file $index_file: $!";

    my $map = {};
    while (<I>)
    {
	chomp;
	my($idx, $val) = split(/\t/);
	$map->{$val} = $idx;
    }
    close(I);
    return $map;
}

#
# Submit jobs. We'll run pf-compute-local families in a job script.
# We assume that we have access to the genus data directory via NFS.
# We also arrange to propagate our environment across so PATH, PERL5LIB, etc.
# are correct.
#
# We submit an array job which will use the index_map to find the right position.
# 
sub submit_jobs
{
    my($kmer_dir, $fam_dir, $index_map, $genus_list, $start, $end, $process_count, $mem, $opt, $to_compute) = @_;

    $start++;
    $end++;

    my @ranges;
    if ($to_compute)
    {
	print "Submitting subset: @$to_compute\n";
	for my $i (@$to_compute)
	{
	    if ($i >= $start && $i <= $end)
	    {
		push(@ranges, [$i, $i]);
	    }
	}
    }
    else
    {
	print "Submit $start-$end $process_count\n";
	push(@ranges, [$start, $end]);
    }

#    $end = $start;

    my $account = $opt->account;
    my $inflation = $opt->inflation;
    my $identity = $opt->identity;
    my $good_cutoff = $opt->good_cutoff;
    my $kser = $opt->kser;
    my $genome_dir = $opt->genome_dir;
    my $container = $opt->container;
    my $dev_container = $opt->dev_container;

    for my $range (@ranges)
    {
	my($start, $end) = @$range;
	next if $end < $start;
	my $batch = <<END;
#!/bin/bash
#SBATCH --account $account
#SBATCH --job-name local-fams
#SBATCH --time 4-0
#SBATCH --mem $mem
#SBATCH --nodes 1-1 --ntasks $process_count
#SBATCH --array $start-$end
#SBATCH --export ALL
#SBATCH --constraint sim

export TMPDIR=/disks/tmp/slurm-\$SLURM_ARRAY_JOB_ID-\$SLURM_ARRAY_TASK_ID
mkdir -p \$TMPDIR
export TEMPDIR=\$TMPDIR

genus=`awk -F \$'\\t' "\\\\\$1 == \$SLURM_ARRAY_TASK_ID { print \\\\\$2}" $fam_dir/genus.index`
hostname
echo genus=\$genus
singularity exec -B /home,/vol,/disks $container sh -c ". $dev_container/user-env.sh; 
pf-compute-local-families --identity $identity --inflation $inflation \\
    --good-cutoff $good_cutoff --parallel \$SLURM_JOB_CPUS_PER_NODE \\
    --genome-dir $genome_dir --kser $kser \\
    $kmer_dir '$fam_dir/genus.data/\$genus'"
END
    print "$batch\n\n";

	my $ok = run(["sbatch", "--parsable"], '<', \$batch);
	print "Sub: ok=$ok\n";
    }
}
