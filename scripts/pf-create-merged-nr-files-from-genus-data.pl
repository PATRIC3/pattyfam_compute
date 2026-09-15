#
# Create a set of NR files for use in loading the annotation server.
# We first read the families file to map peg id to family information.
# Then we read the AA NR data from each genus in genus.data and annotate the
# generated fasta files to have
#
# >peg-id PGF_XXXXXXXX lfam-id GenusName-id
#
# We use the data in genus.data/genus.sizes as created by
# du -sm */nr-seqs > genus.sizes
# to equalize the size of the generated NR files
#

use strict;
use Data::Dumper;
use Getopt::Long::Descriptive;
use File::Path qw(make_path);
use File::Copy qw(copy);
use Proc::ParallelLoop;
use LPTScheduler;
use gjoseqlib;
use POSIX;

my($opt, $usage) = describe_options("%c %o num-nrs nr-base-path fam-file genus-dir work-dir",
				    ["threads|j=i" => "Number of threads for writing output", { default => 1 }],
				    ["help|h" => "Show this help message."]);

print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 5;
    
my $num_nrs = shift;
my $nr_base_path = shift;
my $fam_file = shift;
my $genus_dir = shift;
my $work_dir = shift;

#goto x;
if (-d $work_dir)
{
    die "Work directory $work_dir must not exist\n";
}
make_path($work_dir) or die "Cannot make $work_dir: $!";

my @sort_fields = ("-k2,2", "-k4,4", "-k3,3n");

#
# Load genus.sizes, and use LPTScheduler to compute an even distribution.
#

my $sched = LPTScheduler->new($num_nrs);
open(GS, "<", "$genus_dir/genus.sizes") or die "Cannot read $genus_dir/genus.sizes: $!";
while (<GS>)
{
    chomp;
    my($size, $genus) = m,^(\d+)\s+(.*-\d+)/nr-seqs$,; 
    if (!$genus)
    {
	warn "Error parsing $genus_dir/genus.sizes at line $.: $_\n";
    }
    $sched->add_work($genus, $size);
}

#
# Read family file
#
my %info;
open(F, "<", $fam_file) or die "Cannot read family file $fam_file: $!";
print STDERR "Reading $fam_file\n";
while (<F>)
{
    chomp;
    my($pgf, undef, undef, $peg, undef, undef, $lfam, $genus) = split(/\t/);
    $info{$peg} = [$pgf, $lfam, $genus];
}
print STDERR "Done\n";
close(F);

sub bootstrap
{
    my($num) = @_;
#    my $fh;
#    open($fh, ">", "$nr_base_path.$num") or die "Cannot open $nr_base_path.$num: $!";
    return $num;
}

sub compute
{
    my($glob, $item) = @_;
    # my($fh, $num) = @$glob;
    my $num = $glob;
    my $base = "$genus_dir/$item/nr-seqs";
    opendir(D, $base) or die "Cannot opendir $base: $!";
    my $skip;
    #
    # We run this through a sort on the family IDs so we can easily
    # compute modes later.
    #
    # We save the output in the work directory since we need to merge all of them later
    #
    my $out = "$work_dir/int.$num.$item";
    open(SORT, "|-", "sort", @sort_fields, "-t", "\t", "-S", "50G",
	 "--parallel", 10, "-T", "/disks/tmp", "-o", $out)
	or die "Cannot open sort: $!";
    my $printed;
    for my $d (sort grep { -f "$base/$_" } readdir(D))
    {
	my $path = "$base/$d";
	# print STDERR "Copy $path for $num\n";
	open(F, "<", $path) or die "Cannot read $path:$ !";
	while (<F>)
	{
	    if (my($peg) = /^>(\S+)/)
	    {
		if ($printed)
		{
		    print SORT "\n";
		    $printed = 0;
		}
		my $info = $info{$peg};
		if (!$info)
		{
		    # print STDERR "No info for $peg\n";
		    $skip = 1;
		}
		else
		{
		    my ($pgf, $lfam, $genus) = @$info;
		    # print $fh ">$peg $pgf $lfam $genus\n";
		    print SORT "$peg\t$pgf\t$lfam\t$genus\t";
		    $skip = 0;
		    $printed = 1;
		}
	    }
	    elsif (!$skip)
	    {
		s/\s//g;
		print SORT $_;
		# print $fh $_;
	    }
	}
	if ($printed)
	{
	    print SORT "\n";
	    $printed = 0;
	}
    }
    close(SORT);
}

$sched->run(\&bootstrap, \&compute);
x:
    
#
# we now have $work_dir filled with a bunch of sequence-per-row files. These must
# be merged then rewritten as fasta files. We split into batches at that point.
#

my @files = <$work_dir/int*>;
#
# compute size estimate
#
my $total;
for my $f (@files)
{
    $total += -s $f;
}

my $per_file = ceil($total / $num_nrs);
print STDERR "Per file size $per_file\n";
#
# Merge. Read from the merge and write our output
#

open(MERGE, "-|", "sort", "--merge", @sort_fields, "-t", "\t", "-T", "/disks/tmp",
     "-S", "200G", "--parallel", 30, @files) or die "Cannot open merge: $!";

my $cur_id = 0;
my $cur_size = 0;
my $cur_fh;
my $cur_fam;
open($cur_fh, ">", "$nr_base_path.$cur_id") or die "Cannot open $nr_base_path.$cur_id: $!";

while (<MERGE>)
{
    chomp;
    my($peg, $pgf, $lfam, $genus, $seq) = split("\t");

    my $this_fam = "$pgf $lfam $genus";
    #
    # Only close file on a family boundary
    #
    if ($cur_fam && $this_fam ne $cur_fam)
    {
	if ($cur_size > $per_file)
	{
	    close($cur_fh);
	    undef $cur_fh;
	    $cur_id++;
	    open($cur_fh, ">", "$nr_base_path.$cur_id") or die "Cannot open $nr_base_path.$cur_id: $!";
	    $cur_size = 0;
	}
    }
    
    $cur_fam = $this_fam;    
    write_fasta($cur_fh, [$peg, $this_fam, $seq]);
    $cur_size += length($seq);

}
close(MERGE);
close($cur_fh);
   
    
