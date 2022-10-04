#
# Given the output of
#
# fastANI --ql PATHS-TO-GENUS-GENOMES --rl PATHS-TO-REF-GENOMES -o compare-to-ref.out -t 16
#
# find the best match for each query, and emit a list of queries that don't have good references.
#

use strict;
use Data::Dumper;
use File::Slurp;
use Getopt::Long::Descriptive;
use File::Basename;
use Cwd qw(abs_path);

my($opt, $usage) = describe_options("%c %o genus-dir",
				    ["write" => "Write the updated active.genomes file based on the ANI data"],
				    ["parsable" => "Write stats as a tab-separated line"],
				    ["hit-threshold=i" => "minimum ANI distance to a reference required for a genome to not be added tot he active list",
				     { default => 97 }],
				    ["help|h" => "Show this help message"]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 1;

my $hit_thresh = $opt->hit_threshold;
my $genus_dir = shift;
my $ani_file = "$genus_dir/ani-lookup.out";

my $which = basename(abs_path($genus_dir));

open(A, "<", $ani_file) or die "Cannot open fastani-output $ani_file: $!\n";

my %by_query;

while (<A>)
{
    chomp;
    my($qry_str, $ref_str, $ani, $frag_mappings, $fragment_count) = split(/\t/);

    my($qry) = $qry_str =~ /(\d+\.\d+)/;
    my($ref) = $ref_str =~ /(\d+\.\d+)/;

    next if $qry eq $ref;
    push(@{$by_query{$qry}}, [$qry, $ref, $ani, $frag_mappings, $fragment_count]);
}
close(A);

my @missing;
my @good;

while (my($qry, $dat) = each(%by_query))
{
    my @hits = sort { $b->[2] <=> $a->[2] } @$dat;
    my $best_hit = $hits[0];
    my (undef, $ref, $ani, $frag_mappings, $fragment_count) = @$best_hit;
    if ($ani > $hit_thresh)
    {
	push(@good, $best_hit);
    }
    else
    {
	push(@missing, $best_hit);
    }
}

my @refs = read_file("$genus_dir/ref-genomes.ids");
@refs or die "No references found in $genus_dir/ref-genomes.ids";
chomp @refs;

my @all = read_file("$genus_dir/genomes");
@all or die "No genomes found in $genus_dir/genomes";
chomp @all;

#
# Determine the genomes from the all-list that aren't in either good or missing.
#
my %found = map { $_->[0] => 1 } (@good, @missing);
$found{$_} = 1 foreach @refs;

my @not_found = grep { !$found{$_} } @all;
my $n_not_found = @not_found;

my $n_all_genomes = @all;
my $n_refs = @refs;
my $n_good = @good;
my $n_missing = @missing;
my $total = $n_good + $n_missing + $n_refs + $n_not_found;
if ($total != $n_all_genomes)
{
    warn "Total genomes $n_all_genomes != refs+good+missing+not+found $total\n";
}

my $frac = sprintf "%.1f", 100 * ($n_good / $total);

if ($opt->parsable)
{
    print join("\t", $which, $hit_thresh, $frac, $n_refs, $n_good, $n_missing, $n_not_found, $n_all_genomes), "\n";
}
else
{
    printf "hit_threshold=$hit_thresh refs=$n_refs good=$n_good ($frac%) missing=$n_missing\n";
}

#
# Write our active genomes list. We include the ref/rep set plus any missing.
#

if ($opt->write)
{
    open(ACT, ">", "$genus_dir/active.genomes") or die "Cannot write $genus_dir/active.genomes: $!";
    print ACT "$_\n" foreach @refs, @not_found;
    print ACT "$_->[0]\n" foreach @missing;
    close(ACT);
}

#print "missing: " . join(" ", map { $_->[0] } @missing), "\n";
