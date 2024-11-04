
=head1 NAME

    pf-compute-local-families

=head1 SYNOPSIS

    pf-compute-local-families genus-data-dir

=head1 DESCRIPTION

This is a pattyfam computation support routine.

Given a genus data directory, compute the local families.

=cut

use strict;
use POSIX;
use File::Path qw(make_path);
use File::Basename;
use IPC::Run;
use Time::HiRes 'gettimeofday';
use File::Spec;
use File::Temp;
use File::Slurp;
use Getopt::Long::Descriptive;
use Data::Dumper;

my($opt, $usage) = describe_options("%c %o kmer-dir genus-data-dir",
				    ["logfile|l=s" => "Log file"],
				    ["identity=s" => "Identity for BLAST fallback", { default => 0.5 }],
				    ["inflation=s" => "MCL inflation", { default => 3.0 }],
				    ["good-cutoff=s" => "Fraction of members with unique genomes required to be 'good'", { default => 0.9 }],
				    ["parallel=i" => "Parallel threads", { default => 1 }],
				    ["tmpdir=s" => "Temp dir"],
				    ["genome-dir=s", "Directory holding PATRIC genome data", { default => "/vol/patric3/downloads/genomes" }],
				    ["kser=s" => "Path to kser executable", { default => "kser" }],
				    ["help|h" => "Show this help message."]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 2;

my $kmer_dir = shift;
my $fam_dir = shift;

my $tmpdir = $opt->tmpdir // "$fam_dir/tmp";

if (! -d $tmpdir)
{
    mkdir($tmpdir);
}


#
# Use a system-local temp space for processing the large kmer data sets.
#
my $work_dir = File::Spec->tmpdir() . "/work." . basename($fam_dir);

-d $work_dir and die "Work directory $work_dir already exists\n";
make_path($work_dir);


if ($opt->logfile)
{
    open(LOG, ">>", $opt->logfile) or die "Cannot open " . $opt->logfile . " for append: $!";
}
else
{
    open(LOG, ">&STDERR");
}
LOG->autoflush(1);

my $overall_tstart = gettimeofday;
print LOG "start overall_run $overall_tstart\n";

#
# Create NR if necessary.
#
if (! -s "$fam_dir/nr/peg.synonyms")
{
    -d "$fam_dir/nr" or mkdir("$fam_dir/nr") or die "cannot mkdir $fam_dir/nr: $!";
    my $rc = system("pf-build-nr", "$fam_dir/sources", "$fam_dir/nr/nr",
		    "$fam_dir/nr/peg.synonyms", "$fam_dir/nr/nr-len.btree", "$fam_dir/nr/figids");
    if ($rc == 0)
    {
	$rc = system("pf-compute-nr-seqs", "--genome-dir", $opt->genome_dir,  $fam_dir);
	if ($rc != 0)
	{
	    warn "pf-compute-nr-seqs $fam_dir failed: $rc\n";
	}
    }
    else
    {
	warn "build_nr $fam_dir failed: $rc\n";
    }
}

#
# Start the stateful kmer server
#

my $kser_port_file = "$tmpdir/kser_port";
my @kser_cmd = ($opt->kser, 
		"--reserve-mapping", 10000000,
		"--n-kmer-threads", $opt->parallel,
		"--listen-port-file", $kser_port_file, 0, $kmer_dir);

unlink($kser_port_file);

my $kser_pid = fork;
if (!defined($kser_pid))
{
    die "fork failed: $!";
}
 if ($kser_pid == 0)
{
    print "Starting kser_cmd in $$: @kser_cmd\n";
    exec(@kser_cmd);
    die "Kser did not start $?: @kser_cmd\n";
}

my $kser_port;
while (1)
{
    my $kid = waitpid($kser_pid, WNOHANG);
    if ($kid)
    {
	die "kser server did not start\n";
    }
    if (-s $kser_port_file)
    {
	$kser_port = read_file($kser_port_file);
	print "Read '$kser_port' from $kser_port_file\n";
	chomp $kser_port;
	if ($kser_port !~ /^\d+$/)
	{
	    kill(1, $kser_pid);
	    kill(9, $kser_pid);
	    die "Invalid kser port '$kser_port'\n";
	}
	last;
    }
    sleep(0.1);
}

print "Kser running pid $kser_pid on port $kser_port\n";

my $add_url = "http://localhost:$kser_port/add";
my $matrix_url = "http://localhost:$kser_port/matrix";

$SIG{__DIE__} = sub {
    kill 1, $kser_pid;
};

#
# If we have an nr-seqs directory, use that for sequences. Initialize pegs_to_inflate and
# pegs_for_singleton_fams from the nr-refs file.

my %pegs_to_inflate;
my %pegs_for_singleton_fams;

my $seqs_dir;
if (-d "$fam_dir/nr-seqs")
{
    open(REFS, "<", "$fam_dir/nr-refs") or die "Cannot open $fam_dir/nr-refs: $!";
    while (<REFS>)
    {
	chomp;
	my($ref, @rest) = split(/\t/);
	$pegs_to_inflate{$ref} = [@rest];
	$pegs_for_singleton_fams{$ref} = 1;
    }
    close(REFS);
    $seqs_dir = "$fam_dir/nr-seqs";
}
else
{
    $seqs_dir = "$fam_dir/Seqs";
}

#
# Use pf-annotate-seqs to annotate the chosen sequences.
#
# We pass the /add URL from the kser since that will call the
# proteins with a side effect of adding the signature kmers to its
# internal table that maps pegs and kmers.

my $tstart = gettimeofday;
print LOG "start pf-annotate-seqs $tstart\n";

my $calls_file = "$tmpdir/calls";
my $uncalled_ids_file = "$tmpdir/missed";

my @anno_cmd = ('pf-annotate-seqs',
		"--parallel", $opt->parallel,
		"--truncated-pegs", "$fam_dir/truncated.genes",
		$add_url, $fam_dir, $seqs_dir, $calls_file, $uncalled_ids_file);
print "@anno_cmd\n";
my $rc = system(@anno_cmd);

my $tend = gettimeofday;
my $elap = $tend - $tstart;
print LOG "finish pf-annotate-seqs $tend $elap\n";

$rc == 0 or die "annotation failed with $rc: @anno_cmd\n";

#
# Compute kmer distances.
#
my $unclassified = "$tmpdir/unclassified";

my $tstart = gettimeofday;
print LOG "start pf-compute-kmer-distances $tstart\n";
my $cmd = ["pf-compute-kmer-distances",
	   "--parallel", $opt->parallel,
#	   "--remove-truncated-threshold", 1000,
	   "--truncated-pegs", "$fam_dir/truncated.genes",
	   $matrix_url, $seqs_dir,
	   $calls_file, $uncalled_ids_file,
	   $work_dir];
print STDERR "Run: @$cmd\n";
run($cmd) or die "pf-compute-kmer-distances failed: $?";

my $tend = gettimeofday;
my $elap = $tend - $tstart;
print LOG "finish pf-compute-kmer-distances $tend $elap\n";
kill 1, $kser_pid;

#
# Compute kmer clusters with MCL
#

my $tstart = gettimeofday;
print LOG "start pf-compute-kmer-clusters $tstart\n";

my $cmd = ["pf-compute-kmer-clusters",
	   "--parallel", $opt->parallel,
	   $work_dir, "$fam_dir/gene.names", "$fam_dir/nr-refs", "$fam_dir/kmer.clusters", "$fam_dir/unclassified"];
print STDERR "Run: @$cmd\n";
run($cmd) or die "pf-compute-kmer-clusters failed: $?";

my $tend = gettimeofday;
my $elap = $tend - $tstart;
print LOG "finish pf-compute-kmer-clusters $tend $elap\n";

#
# Process unclassified sequences using BLAST.
#

my $tstart = gettimeofday;
print LOG "start pf-compute-blast-clusters $tstart\n";

if (-s "$fam_dir/unclassified")
{
    my $cmd = ["pf-compute-blast-clusters",
	       "--parallel", $opt->parallel,
	       $work_dir, "$fam_dir/gene.names", "$fam_dir/nr-seqs",
	       "$fam_dir/unclassified", "$fam_dir/blast.clusters"];
    print STDERR "Run: @$cmd\n";
    run($cmd) or die "pf-compute-blast-clusters failed: $?";
}
else
{
    print STDERR "No unclassified proteins to process\n";
    open(F, ">", "$fam_dir/blast.clusters") or die "Cannot write $fam_dir/blast.clusters: $!";
    close(F);
}

my $tend = gettimeofday;
my $elap = $tend - $tstart;
print LOG "finish pf-compute-blast-clusters $tend $elap\n";

my $tstart = gettimeofday;
print LOG "start compute singletons $tstart\n";

#
# Determine if we have to mark any singleton families.
#
# Read our cluster files and remove the pegs therein from our
# list of pegs. Anything left goes to singleton.fams.
#

for my $fam_file ("$fam_dir/blast.clusters", "$fam_dir/kmer.clusters")
{
    open(I, "<", $fam_file) or die "Cannot open $fam_file: $!";

    while (<I>)
    {
	my($fam, $subfam, $peg) = split(/\t/);
	delete $pegs_for_singleton_fams{$peg};
    }
    close(I);
}

#
# First process the called proteins for singletons.
#

my %subfam;
open(CALLS, "<", "$tmpdir/calls") or die "Cannot open $tmpdir/calls: $!";
open(B, ">", "$fam_dir/singleton.fams") or die "Cannot write $fam_dir/singleton.fams: $!";
while (<CALLS>)
{
    chomp;
    my($peg, $fun) = split(/\t/);
    if ($pegs_for_singleton_fams{$peg})
    {
	my $sub;
	if (exists($subfam{$fun}))
	{
	    $subfam{$fun}++;
	    $sub = $subfam{$fun};
	}
	else
	{
	    $sub = 1;
	    $subfam{$fun} = 1;
	}
	
	print B "$fun\t$sub\t$peg\n";
	delete $pegs_for_singleton_fams{$peg};
	if (!%pegs_for_singleton_fams)
	{
	    last;
	}
    }
}
close(CALLS);

#
# Anything else is hypothetical.
#
my $fun = "hypothetical protein";
my $sub;
if (exists($subfam{$fun}))
{
    $subfam{$fun}++;
    $sub = $subfam{$fun};
}
else
{
    $sub = 1;
    $subfam{$fun} = 1;
}
for my $peg (sort keys %pegs_for_singleton_fams)
{
    print B "$fun\t$sub\t$peg\n";
    $sub++;
}
close(B);

my $tend = gettimeofday;
my $elap = $tend - $tstart;
print LOG "finish singletons $tend $elap\n";

my $tstart = gettimeofday;
print LOG "start pf-finalize-local-families $tstart\n";

run(["pf-finalize-local-families",
     "--truncated-pegs", "$fam_dir/truncated.genes",
     $fam_dir, "$fam_dir/blast.clusters", "$fam_dir/kmer.clusters", "$fam_dir/singleton.fams"]);
run(["pf-inflate-families", 
     "-o", "$fam_dir/local.family.members.expanded",
     $fam_dir, "$fam_dir/local.family.members", 2, 5]);

my $tend = gettimeofday;
my $elap = $tend - $tstart;
print LOG "finish pf-finalize-local-families $tend $elap\n";

#
# With the local families finalized, we can compute alignments.
#
# We use /dev/shm as a work directory due to the heavy small-file I/O
# incurred by this step.
#

my $align_work_dir = File::Temp->newdir(DIR => "/dev/shm", TEMPLATE => "pfalign-XXXXXX");

my $tstart = gettimeofday;
print LOG "start pf-compute-local-family-alignments.pl $tstart\n";

run(["pf-compute-local-family-alignments",
     "--genome-dir", $opt->genome_dir,
     "--parallel", $opt->parallel,
     $fam_dir, $align_work_dir, $fam_dir]);

my $tend = gettimeofday;
my $elap = $tend - $tstart;
print LOG "finish pf-compute-local-family-alignments.pl $tend $elap\n";

my $overall_elap = $tend - $overall_tstart;

print LOG "finish overall_run $tend $overall_elap\n";

system("rm", "-rf", $work_dir);

#
# Wrapper around IPC::Run::run that does failure checking. Die on failure.
sub run
{
    my(@cmd) = @_;
    my $ok = IPC::Run::run(@cmd);
    if (!$ok)
    {
	die "Failure $? running " . Dumper(\@cmd);
    }
    return $ok;
			   
}
