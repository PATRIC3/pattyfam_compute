#
# Given a genus data directory, retrieve (or link to download files) contigs
# for all genomes in the genus, creating a list to past to fastANI.
#
# Also construct a list of the reference and representative genomes in the genus,
# and create a list for fastANI.
#
# Run fastANI with those lists, saving the results to the genus data directory.
#

use strict;
use Getopt::Long::Descriptive;
use P3DataAPI;
use Data::Dumper;
use File::Slurp;
use File::Path qw(make_path);
use IPC::Run qw(run);
use List::Util qw(min);
use Proc::ParallelLoop;

my($opt, $usage) = describe_options("%c %o genus-dir",
				    ["parallel|j=i" => "Parallelism for fastANI", { default => 1 }],
				    ["help|h" => "Show this help message"],
				   );

print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 1;

my $dir = shift;

my $contig_cache = "$dir/contigs";
make_path($contig_cache);

my @genomes = read_file("$dir/genomes");
chomp @genomes;

my $api = P3DataAPI->new();

my @refs = find_references(\@genomes);
my %refs = map { $_ => 1 } @refs;

#
# Save list of references.
#
my $ref_genomes_file = "$dir/ref-genomes.ids";
open(F, ">", $ref_genomes_file) or die "Cannot write $ref_genomes_file: $!";
print F "$_\n" foreach @refs;
close(F);

#
# Now locate genomes either in the downloads directory or if missing there, create contigs in cache directory.
#

my $all_contigs = "$dir/all-contigs.lst";
my $ref_contigs = "$dir/ref-contigs.lst";

my $genome_dir = "/vol/patric3/downloads/genomes";

open(ALL, ">", $all_contigs) or die "Cannot write $all_contigs: $!";
open(REF, ">", $ref_contigs) or die "Cannot write $ref_contigs: $!";

#
# In parallel, load from API if needed.
#
my $n_loaders = min(4, $opt->parallel);

my @work;

for my $gid (@genomes)
{
    my $gpath = "$genome_dir/$gid/$gid.fna";
    if (! -s $gpath)
    {
	$gpath = "$contig_cache/$gid.fna";

	if (! -s $gpath)
	{
	    #
	    # Make our work list, and skip printing to REF/ALL until we know what has
	    # succeeded.
	    #
	    push(@work, [$gid, $gpath]);
	    next;
	}

    }

    if ($refs{$gid})
    {
	print REF "$gpath\n";
    }
    else
    {
	print ALL "$gpath\n";
    }
}
print Dumper(\@work);
print "Loading from API for " . join(" ", map { $_->[0] } @work) . "\n";
pareach(\@work, sub {
    my($w) = @_;
    my($gid, $gpath) = @$w;
    my $tmp = "$gpath.tmp";
    print STDERR "Loading $gid from API\n";
    my $ok = run(["p3-genome-fasta", "--contig", $gid], ">", $tmp);
    if (!$ok || ! -s $tmp)
    {
	warn "Failure loading $gid from API\n";
    }
    rename($tmp, $gpath);

}, { Max_Workers => $n_loaders });

for my $w (@work)
{
    my($gid, $gpath) = @$w;
    if (! -s $gpath)
    {
	print STDERR "Failed to load '$gid' '$gpath'\n" . Dumper($w);
	next;
    }
    
    if ($refs{$gid})
    {
	print REF "$gpath\n";
    }
    else
    {
	print ALL "$gpath\n";
    }
}

close(REF);
close(ALL);

my $fastani_output = "$dir/ani-lookup.out";
my @cmd = ('fastANI', '--ql', $all_contigs, '--rl', $ref_contigs,  '-o', $fastani_output, "-t", $opt->parallel);

print "@cmd\n";
my $ok = run(\@cmd,
	     ">", "$dir/ani-lookup.stdout",
	     "2>", "$dir/ani-lookup.stderr");
$ok or die "Error $? running @cmd\n";

sub find_references
{
    my($genomes) = @_;

    my @work = @$genomes;

    my @refs;
    while (@work)
    {
	my @chunk = splice(@work, 0, 1000);

        my $q = join(" OR ", map { "genome_id:\"$_\"" } @chunk);
        my $res = $api->solr_query_list("genome",
					[q => $q, fq => "reference_genome:(Reference OR Representative)", fl => "genome_id,reference_genome"]);
        for my $ent (@$res)
        {
	    push(@refs, $ent->{genome_id});
        }
    }
    return @refs;
}
