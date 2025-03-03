#
# Submit the set of merges to the slumr scheduler
#

use strict;
use Data::Dumper;
use Getopt::Long::Descriptive;
use JSON::XS;
use File::Path 'make_path';
use Cwd qw(abs_path);
use File::Copy;
use File::Slurp;
use Statistics::Descriptive;
use File::Basename;
use gjoseqlib;
use IPC::Run qw(run);
use Template;

my($opt, $usage) = describe_options("%c %o work-file slurm-account",
				    ["parallel=i", "Number of threads per job", { default => 1 }],
				    ["range=s", "Submit these array elements"],
				    ["max-tasks=i", "Only submit this many jobs"],
				    ["memory=s", "Memory request", { default => "32G" }],
				    ["reservation=s", "Slurm reservaton"],
				    ["nodelist|w=s", "Slurm nodelist"],
				    ["help|h", "Show this help message"]);

print($usage->text), exit if $opt->help;
die($usage->text) if @ARGV != 2;

my $work_file = shift;
my $slurm_account = shift;

my $n_tasks;
my $merge_dir;

if (open(F, "<", $work_file))
{
    my $txt = read_file(\*F);
    my $dat = decode_json($txt);
    $dat or die "Cannot decode work data from $work_file\n";
    $merge_dir = $dat->{merge_dir};
    -d $merge_dir or die "Merge_dir $merge_dir does not exist\n";
    $n_tasks = @{$dat->{bins}};
}
else
{
    die "Cannot open $work_file: $!";
}

if ($opt->max_tasks && $n_tasks > $opt->max_tasks)
{
    $n_tasks = $opt->max_tasks;
}

my $output_dir = "$merge_dir/slurm.out";
make_path($output_dir);

my %vars = (account => $slurm_account,
	    output_dir => $output_dir,
	    mem => $opt->memory,
	    reservation => $opt->reservation,
	    nodelist => $opt->nodelist,
	    work_file => abs_path($work_file),
	    ntasks => $opt->parallel,
	    );

my $templ = new Template;
my $batch = "batch.merge.$$";
$templ->process(\*DATA, \%vars, $batch) or die $templ->error();
print "Wrote batch to $batch\n";
my $range = $opt->range // "1-$n_tasks";

run(["sbatch", "-a", $range, $batch]);

__DATA__
#!/bin/bash -x
#SBATCH --account [% account %]
#SBATCH --output "[% output_dir %]/slurm-%A_%a.out"
#SBATCH --job-name "GFM"
#SBATCH --time 4-0
#SBATCH --mem [% IF mem %][% mem %][% ELSE %] 64G [% END %]
#SBATCH --nodes 1-1
#SBATCH --ntasks [% ntasks %]
#SBATCH --export NONE
[% IF reservation -%]
#SBATCH --reservation [% reservation %]
[% END -%]
[% IF nodelist -%]
#SBATCH --nodelist [% nodelist %]
[% END -%]
   
. /home/olson/P3/dev-families/dev_container/user-env.sh

pf-merge-stage-1-guts \
    --parallel $SLURM_JOB_CPUS_PER_NODE \
    [% work_file %] $SLURM_ARRAY_TASK_ID

    
