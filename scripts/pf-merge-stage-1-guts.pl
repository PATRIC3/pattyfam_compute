
#
#
# This is the worker script that performs one chunk work as determined by
# pf-merge-stage-1 --save-work work-file
#
# We are given the work file and the work index to be computed on.
#
# Given the set of fasta files, we use the kmer distance clustering algorithm to
# compute suggested family merges.
#
# We also run mcl to compute the clusters based on that distance.
#

use DBI;
use strict;
use Data::Dumper;
use Getopt::Long::Descriptive;
use IO::Handle;
use JSON::XS;
use File::Path 'make_path';
use File::Copy;
use File::Slurp;
use Statistics::Descriptive;
use File::Basename;
use gjoseqlib;
use Algorithm::Numerical::Shuffle qw /shuffle/;
use POSIX ':sys_wait_h';
use IPC::Run;
use List::Util qw(max min reduce);
use Time::HiRes qw(gettimeofday);

my($opt, $usage) = describe_options("%c %o work-file work-index",
				    ["parallel=i", "Number of threads", { default => 1 }],
				    ["min-kmers-in-common=i", "Minimum kmers in common required", { default => 5 }],
				    ["min-protein-len=i", "Minimum protein length to consider", { default => 0 }],
				    ["log|l=s", "Trace logfile"],
				    ["help|h", "Show this help message"]);

print($usage->text), exit if $opt->help;
die($usage->text) if @ARGV != 2;

if ($opt->log)
{
    open(LOG, ">>", $opt->log) or die "Cannot open logfile " . $opt->log . ": $!";
    LOG->autoflush(1);
}
else
{
    open(LOG, ">", "/dev/null");
}

my $work_file = shift;
my $work_index = shift;
my($kmer_dir, $merge_dir, $dist_dir, $inflation);
my $work_item;

if (open(F, "<", $work_file))
{
    my $txt = read_file(\*F);
    my $dat = decode_json($txt);
    $dat or die "Cannot decode work data from $work_file\n";
    $merge_dir = $dat->{merge_dir};
    $kmer_dir = $dat->{kmer_dir};
    $dist_dir = $dat->{dist_dir};
    for my $k (qw(merge_dir kmer_dir dist_dir))
    {
	if ($dat->{$k})
	{
	    if (! -d $dat->{$k})
	    {
		die "Data dir $dat->{$k} does not exist\n";
	    }
	}
	else
	{
	    die "Work file $work_file missing key '$k'\n";
	}
    }
    $inflation = $dat->{inflation} or die "Work file $work_file missing key 'inflation'\n";
    
    (my $size, $work_item) = @{$dat->{bins}->[$work_index - 1]};
}
else
{
    die "Cannot open $work_file: $!";
}

my $tmpdir = "$merge_dir/tmp";
make_path($tmpdir);

my $mcl_dir = "$merge_dir/mcl/$inflation";
make_path($mcl_dir);

my $stat_dir = "$merge_dir/kmer_stats";
make_path($stat_dir);

# $work_item is a list: [ [ $fn_index, [ list of fa files ] ], size ]

for my $item (@$work_item)
{
    do_work(undef, $item->[0]);
}

sub do_work
{
    my($glob, $item) = @_;
    
    my($fn, $fa_list) = @$item;
    my $n_fa = @$fa_list;
    print "$$ $fn: $n_fa files\n";

    my $t1 = gettimeofday;

    my $dist_file = "$dist_dir/$fn";
    my $mcl_file = "$mcl_dir/$fn";

    my $tmp = File::Temp->new(TEMPLATE => "kmer_inputs_${fn}_XXXXXX", TMPDIR => 1);
    print $tmp "$_\n" foreach @$fa_list;
    close($tmp);

    # my $tmp_dist = File::Temp->new(TEMPLATE => "kmer_dist_${fn}_XXXXXX", TMPDIR => 1);
    # $dist_file = "$tmp_dist";

    $dist_file = "/dev/fd/1";

    my @dist_cmd = ("kmers-matrix-distance-files",
		    "--n-threads", $opt->parallel,
		    "--min-kmers-in-common", $opt->min_kmers_in_common,
		    "--min-protein-len", $opt->min_protein_len,
		    "--kmer-stats", "$stat_dir/$fn",
		    "--input-file-list", "$tmp",
		    $kmer_dir,
		    $dist_file,
	      );

#    my $mcl_tmp = File::Temp->new(TEMPLATE => "mcl_${fn}_XXXXXX", DIR => "/dev/shm");
#    close($mcl_tmp);
    my $ok = IPC::Run::run(\@dist_cmd,
			   '|',
			   ["mcl", "-",
			    "-o", $mcl_file,
			    "--abc",
			    "-I", $inflation,
			    "-te", $opt->parallel],
			   "2>", "$mcl_file.mcl.out");
    $ok or die "mcl failed with $?\n";
    

    # print STDERR "Run @cmd\n";
    # my $rc = system(@cmd);

    # my $t2 = gettimeofday;
    # my $elap = $t2 - $t1;
    # my $sz = 0;
    # $sz += -s $_ foreach @$fa_list;
    # print join("\t", "ITEM", $fn, $sz, $elap), "\n";

    # if ($rc != 0)
    # {
    # 	die "Failure $rc running @cmd\n";
    # }

    # if (-s $dist_file == 0)
    # {
    # 	open(my $fh, ">", "$mcl_file");
    # 	close($fh);
    # 	open(my $fh, ">", "$mcl_file.mcl.out");
    # 	print $fh "Wrote zero length file due to zero-length input $dist_file\n";
    # 	close($fh);
    # }
    # else
    # {
    # 	my $ok = IPC::Run::run(["cut", "-f", "1,2,4", $dist_file],
    # 			       '|',
    # 			       ["mcl", "-",
    # 				"-o", $mcl_file,
    # 				"--abc",
    # 				"-I", $inflation,
    # 				"-te", $opt->parallel], "2>", "$mcl_file.mcl.out");
    # 	$ok or die "mcl failed with $?\n";
    # }

}

