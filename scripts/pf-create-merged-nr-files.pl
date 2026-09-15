#
# Given the output of pf-create-nrs-for-annotation, create a smaller set of large
# merged NR files for loading into the annotation server.
#

use strict;
use Data::Dumper;
use Getopt::Long::Descriptive;
use File::Path qw(make_path);
use File::Copy qw(copy);
use Proc::ParallelLoop;
use LPTScheduler;
use gjoseqlib;

my($opt, $usage) = describe_options("%c %o num-nrs nr-base-path input-fasta-files",
				    ["threads|j=i" => "Number of threads for writing output", { default => 1 }],
				    ["help|h" => "Show this help message."]);

print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV < 3;

my $num_nrs = shift;
my $nr_base_path = shift;
my @input_fasta = @ARGV;

#
# Scan input files, and use LPTScheduler to compute an even distribution.
#

my $sched = LPTScheduler->new($num_nrs);

for my $fasta (@input_fasta)
{
    $sched->add_work($fasta, -s $fasta);
}

sub bootstrap
{
    my($num) = @_;
    my $fh;
    open($fh, ">", "$nr_base_path.$num") or die "Cannot open $nr_base_path.$num: $!";
    return [$fh, $num];
}

sub compute
{
    my($glob, $item) = @_;
    my($fh, $num) = @$glob;
    print "Copy $item for $num\n";
    copy($item, $fh);
}

$sched->run(\&bootstrap, \&compute);
