
#
# Analyze the max memory usage per max job item for merge guts runs.
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

my($opt, $usage) = describe_options("%c %o work-file sacct-data",
				    ["help|h", "Show this help message"]);

print($usage->text), exit if $opt->help;
die($usage->text) if @ARGV != 2;


my $work_file = shift;
my $sacct_data = shift;

#
# $sacct_data 
# 15159034_19.batch  COMPLETED   00:11:46     10.18G
#

my($kmer_dir, $merge_dir, $dist_dir, $inflation);
my $work_item;

open(S, "<", $sacct_data) or die "Cannot open $sacct_data: $!";
open(F, "<", $work_file) or die "Cannot open $work_file: $!";

my $txt = read_file(\*F);
my $dat = decode_json($txt);
$dat or die "Cannot decode work data from $work_file\n";

my $bins = $dat->{bins};

while (<S>)
{
    chomp;
    my($job, $idx, $hr, $min, $sec, $mem) =
	/^\s*(\d+)_(\d+)\S+\s+\S+\s+(\d+):(\d+):(\d+)(?:\s+([.\d]+)G)?/;
    # print Dumper($job, $idx, $hr, $min, $sec, $mem);
    
    $job or die "Cannot parse '$_'\n";

    next unless $mem;

    my $elap = (3600 * $hr + 60 * $min + $sec) / 60;
    my $bdata = $bins->[$idx - 1];
    my($total, $elts) = @$bdata;

    my $max_size = 0;
    my $max_files = 0;
    for my $ent (@$elts)
    {
	my($work, $size) = @$ent;
	my($fn, $files) = @$work;
	# print "fn=$fn size=$size\n";
	$max_size = $size if $size > $max_size;
	$max_files = @$files if @$files > $max_files;
    }
    print join("\t", $idx, $total, $max_size, $max_files, scalar @$elts, $elap, $mem), "\n";
}
