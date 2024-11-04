
=head1 NAME

    pf-compute-kmer-distances

=head1 SYNOPSIS

    pf-compute-kmer-distances matrix-url seqs-dir calls-file uncalled-ids-file work-dir

=head1 DESCRIPTION

This is a pattyfam computation support routine.

Use a loaded kser via the given matrix URL to compute all to all kmer distances
for the given sequences.     

Generate a local ID so we do not create filenames with functions.

We create a file to hold the mapping, and create fasta files for
each family in the families.fasta directory.

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

my($opt, $usage) = describe_options("%c %o matrix-url seqs-dir calls-file uncalled-ids-file work-dir",
				    ["parallel=i" => "Parallel threads", { default => 1 }],
				    ["num-fams|n=i" => "Only process N families (for debugging)"],
				    ["truncated-pegs=s" => "File containing possibly truncated pegs"],
				    ["remove-truncated-threshold=i" =>
				     "Family size threshold at which truncated (at the end of contigs) features are removed", { default =>  30_000 }],
				    ["help|h" => "Show this help message."]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 5;

my $matrix_url = shift;
my $seqs_dir = shift;
my $calls_file = shift;
my $uncalled_ids_file = shift;
my $work_dir = shift;

-d $work_dir or die "Work directory $work_dir is missing\n";

my $map_file = "$work_dir/map";
my $fasta_dir = "$work_dir/fasta";
my $dist_dir = "$work_dir/dist";
my $mcl_dir = "$work_dir/mcl";

make_path($fasta_dir, $dist_dir, $mcl_dir);

-d $seqs_dir or die "Seqs directory $seqs_dir does not exist\n";
-f $calls_file or die "Called file $calls_file does not exist\n";

open(MAP, ">", $map_file) or die "Cannot open $map_file for writing: $!";

#
# Read sequence data, populating %seqs.
#

my %seqs;
opendir(DH, $seqs_dir) or die "Cannot opendir $seqs_dir: $!";

while (my $p = readdir(DH))
{
    next if $p =~ /^\./;
    my $f = "$seqs_dir/$p";
    open(F, "<", $f) or die "Cannot open $f: $!";
    while (my($id, $def, $seq) = read_next_fasta_seq(\*F))
    {
	$seqs{$id} = $seq;
    }
    close(F);
}
closedir(DH);

#
# Load the truncated-pegs data if present
#
my %truncated_pegs;

if ($opt->truncated_pegs)
{
    if (open(TRUNC, "<", $opt->truncated_pegs))
    {
	while (<TRUNC>)
	{
	    chomp;
	    my($id) = split(/\t/);
	    $truncated_pegs{$id} = 1;
	}
	close(TRUNC);
    }
    else
    {
	warn "Cannot open truncated pegs file " . $opt->truncated_pegs . ": $!\n";
    }
}

#
# We sort the calls data on function so we can compute distance on a per-function basis.
#
# For each function, we write a file named LFIGXXXXX to store the sequence
# data for the proteins with the given function.
#
# For very large families, we will remove any truncated pegs.
# We collect the set where this has been done; if the families are still too large
# we stratify on size (0-100 AA, 101-200 AA, 201 AA and above) creating new groups.
#

$ENV{LANG} = "C";
open(C, "-|", "sort", "-t", "\t", "-k", 2, $calls_file) or die "Cannot open sort on $calls_file: $!";

my @to_stratify;

my $last;
my $next_id = "LFIG00000";
my $cur_fam;
my %fam_func;
my %fam_sequence_size;
my %fam_sequence_count;
while (<C>)
{
    chomp;
    my($fid, $func) = split(/\t/);
    if ($func ne $last)
    {
	if ($last)
	{
	    close(FH);
	    # print STDERR "Finish $cur_fam, size=$fam_sequence_count{$cur_fam}\n";
	    if (0 && $fam_sequence_count{$cur_fam} > $opt->remove_truncated_threshold)
	    {
		my $n = remove_truncated("$fasta_dir/$cur_fam");
		if ($n > $opt->remove_truncated_threshold)
		{
		    push(@to_stratify, $cur_fam);
		}
	    }
	    my($n) = $next_id =~ /(\d+)/;
	    last if $opt->num_fams && $n > $opt->num_fams;
	}
	open(FH, ">", "$fasta_dir/$next_id") or die "Cannot open $fasta_dir/$next_id for writing: $!";
	$cur_fam = $next_id;
	$fam_func{$cur_fam} = $func;
	print MAP "$next_id\t$func\n";
	$next_id++;
	$last = $func;
    }
    if (!$truncated_pegs{$fid})
    {
	print_alignment_as_fasta(\*FH, [$fid, $func, $seqs{$fid}]);
	$fam_sequence_size{$cur_fam} += length($seqs{$fid});
	$fam_sequence_count{$cur_fam}++;
    }
    else
    {
	print STDERR "Skip truncated $fid\n";
    }
    
}
# print STDERR "Finish $cur_fam, size=$fam_sequence_count{$cur_fam}\n";
if (0 && $fam_sequence_count{$cur_fam} > $opt->remove_truncated_threshold)
{
    my $n = remove_truncated("$fasta_dir/$cur_fam");
    if ($n > $opt->remove_truncated_threshold)
    {
	push(@to_stratify, $cur_fam);
    }
}
close(C);
close(FH);

my %hypo_fam;

if (-s $uncalled_ids_file)
{
    #
    # Create one more for the sequences for which kmers were not able to
    # determine a function. Note that these sequences may still have signature
    # kmers; it may not have been possible for the function assignment
    # algorithm to compute an unambiguous function.
    #
    
    my $cur_fam = $next_id++;
    $next_id++;
    
    open(FH, ">", "$fasta_dir/$cur_fam");
    $fam_func{$cur_fam} = "hypothetical protein";

    $hypo_fam{$cur_fam} = 1;
    
    #
    # Need to initialize here so the key exists so it gets processed.
    #
    $fam_sequence_size{$cur_fam} = $fam_sequence_count{$cur_fam} = 0;

    print MAP "$cur_fam\thypothetical protein\n";
    open(M, "<", $uncalled_ids_file) or die "cannot open missed file $uncalled_ids_file: $!";
    while (<M>)
    {
	chomp;
	if (/(fig\|\d+\.\d+\.peg\.\d+)/ &&
	    !$truncated_pegs{$1})
	{
	    print_alignment_as_fasta(\*FH, [$1, "hypothetical protein", $seqs{$1}]);
	    $fam_sequence_size{$cur_fam} += length($seqs{$1});
	    $fam_sequence_count{$cur_fam}++;
	}
    }
    close(M);
}
close(FH);

# print STDERR "Finish $cur_fam, size=$fam_sequence_count{$cur_fam}\n";
if ($fam_sequence_count{$cur_fam} > $opt->remove_truncated_threshold)
{
    my $n = remove_truncated("$fasta_dir/$cur_fam");
    if ($n > $opt->remove_truncated_threshold)
    {
	push(@to_stratify, $cur_fam);
    }
}

#
# Perform further stratification on families that need it.
#
for my $fam (@to_stratify)
{
    #
    # Rename the current fasta/fam file to fam.pre_stratify
    # Create three new families - fasta/fam plus two new ones.
    #

    my $old = "$fasta_dir/$fam.pre_stratify";
    my $file1 = "$fasta_dir/$fam";
    my $fam2 = $next_id++;
    my $fam3 = $next_id++;
    my $file2 = "$fasta_dir/$fam2";
    my $file3 = "$fasta_dir/$fam3";

    if ($hypo_fam{$fam})
    {
	$hypo_fam{$fam2} = $hypo_fam{$fam3} = 1;
    }

    print STDERR "stratify $fam into $fam2 $fam3\n";

    $fam_func{$fam2} = $fam_func{$fam};
    $fam_func{$fam3} = $fam_func{$fam};
    
    rename("$fasta_dir/$fam", $old);

    open(F1, ">", $file1) or die "cannot open $file1: $!";
    open(F2, ">", $file2) or die "cannot open $file2: $!";
    open(F3, ">", $file3) or die "cannot open $file3: $!";

    print MAP "$fam2\t$fam_func{$fam}\n";
    print MAP "$fam3\t$fam_func{$fam}\n";

    for my $f ($fam, $fam2, $fam3)
    {
	$fam_sequence_count{$f} = $fam_sequence_size{$f} = 0;
    }

    open(IN, "<", $old) or die "cannot open $old: $!";
    while (my $ent = read_next_fasta_seq(\*IN))
    {
	my $sz = length($ent->[2]);
	if ($sz <= 100)
	{
	    $fam_sequence_count{$fam}++;
	    $fam_sequence_size{$fam} += $sz;
	    write_fasta(\*F1, $ent);
	}
	elsif ($sz <= 200)
	{
	    $fam_sequence_count{$fam2}++;
	    $fam_sequence_size{$fam2} += $sz;
	    write_fasta(\*F2, $ent);
	}
	else
	{
	    $fam_sequence_count{$fam3}++;
	    $fam_sequence_size{$fam3} += $sz;
	    write_fasta(\*F3, $ent);
	}
    }
    close(F1);
    close(F2);
    close(F3);
    close(IN);
}

undef %seqs;

close(MAP);

#
# Write the hypothetical family proteins for later use to determine the
# proteins that we need to use BLAST to cluster.
#

open(HP, ">", "$work_dir/hypo.prots") or die "Cannot write $work_dir/hypo.prots: $!";
for my $fam (keys %hypo_fam)
{
    open(F, "<", "$fasta_dir/$fam") or die "cannot open hypo fam $fasta_dir/$fam: $!";
    while (<F>)
    {
	if (/^>(\S+)/)
	{
	    print HP "$1\n";
	}
    }
    close(F);
}
close(HP);


#
# Now that we have a directory of per-function sequence data, we will
# invoke the kser engine via kmer_guts_net to compute the
# kmer distances between each member of each family. We will do this in
# parallel.

my $kgn = KmerGutsNet->new($matrix_url);

#
# Construct a pool of work 10x our parallelism, and schedule
# families into it using LPT algorithm.
#

my @fams_in_order = sort { $fam_sequence_size{$b} <=> $fam_sequence_size{$a} } keys %fam_sequence_size;

my $n_buckets = $opt->parallel * 10;
my $scheduler = LPTScheduler->new($n_buckets);

for my $f (@fams_in_order)
{
    $scheduler->add_work([$f, $fam_func{$f}], $fam_sequence_size{$f});
}

print STDERR "$0 pid $$ starting parallel\n";

$scheduler->run(sub {}, \&work, $opt->parallel);

sub work
{
    my($global, $ent) = @_;

    my($fam, $fun) = @$ent;
    print STDERR "Proc $fam $fun\n";
    
    #
    # If family has a single sequence, the MCL file just has the sequence ID. Singleton.
    #
    
    if ($fam_sequence_count{$fam} == 0)
    {
	open(M, ">", "$mcl_dir/$fam") or die "Cannot write $mcl_dir/$fam: $!";
	open(MO, ">", "$mcl_dir/$fam.mcl.out") or die "Cannot write $mcl_dir/$fam.mcl.out: $!";
	print MO "Not running MCL for empty family $fam $fun\n";
	close(M);
	close(MO);
    }
    elsif ($fam_sequence_count{$fam} == 1)
    {
	open(M, ">", "$mcl_dir/$fam") or die "Cannot write $mcl_dir/$fam: $!";
	open(MO, ">", "$mcl_dir/$fam.mcl.out") or die "Cannot write $mcl_dir/$fam.mcl.out: $!";
	open(S, "<", "$fasta_dir/$fam");
	my $l = <S>;
	if ($l =~ /^>(\S+)/)
	{
	    print M "$1\n";
	    print MO "Not running MCL for singleton $1\n";
	}
	else
	{
	    die "Singleton fam file $fasta_dir/fam not of expected format: $l";
	}
	close(M);
	close(S);
	close(MO);
    }
    else
    {
	
	my $dist_size = -s "$dist_dir/$fam";
	my $dist_exists = -f _;
	if (! $dist_exists)
	{
	    print STDERR "Processing $dist_dir/$fam as it does not exist\n";
	    
	    my $data;
	    my $res;
	    eval {
		$data = read_file("$fasta_dir/$fam");
		# print "Got data $data\n";
		$res = $kgn->request($data, { d => 1, a => 1 });
	    };
	    if ($@)
	    {
		die "Error requesting data for $fasta_dir/$fam: $@";
	    }
	    # print "writing size=" . length($res) . "\n";
	    write_file("$dist_dir/$fam", $res // "");

	    $dist_size = -s "$dist_dir/$fam";;
	    print STDERR "$dist_dir/$fam now $dist_size bytes\n";
	}
	
	#
	# Verify at least the first line if the file isn't empty.
	#
	if ($dist_size > 0)
	{
	    open(my $dfh, "<", "$dist_dir/$fam") or die "Cannot open $dist_dir/$fam: $!";
	    my $l = <$dfh>;
	    
	    if ($l !~ /^fig/)
	    {
		my $s = substr($l, 0, 80);
		my $len = length($l);
		die "Dist file $dist_dir/$fam does not start with 'fig': '$s' (line length $len)";
	    }
	}
    }
}

#
# Remove truncated pegs from the given fasta file.
#
sub remove_truncated
{
    my($file) = @_;
    if (!%truncated_pegs)
    {
	warn "No truncated pegs file so not removing from $file\n";
	return;
    }
    my $new = "$file.new";
    open(IN, "<", $file) or die "Cannot read $file: $!";
    open(OUT, ">", $new) or die "Cannot write $file: $!";
    my $n_in = 0;
    my $n_out = 0;
    while (my $ent = read_next_fasta_seq(\*IN))
    {
	$n_in++;
	if (!$truncated_pegs{$ent->[0]})
	{
	    $n_out++;
	    write_fasta(\*OUT, $ent);
	}
    }
    close(IN);
    close(OUT);
    rename($file, "$file.orig") or die "Error renaming $file to $file.org: $!";
    rename($new, $file) or die "Error renaming $new to $file: $!";
    print STDERR "Truncate of $file complete. $n_in input sequences, $n_out output sequences\n";
    return $n_out;
}
