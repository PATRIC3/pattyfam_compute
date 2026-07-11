#
# Given a fams file as the output of pf-merge-stage-3, split it into a defs and members file
#
# Optionally (with --write-merged-fams-file) also rewrite the original file. Used
# when we have translated GFXXXXXXXX to PGF_XXXXXXXX
#

use strict;
use Getopt::Long::Descriptive;
use gjoseqlib;

my($opt, $usage) = describe_options("%c %o fams-file defs-file members-file singleton-defs singleton-members",
				    ["write-merged-fams-file=s", "Write the merged fams file as well, with singletons removed"],
				    ["cache-dir=s", "Genome cache dir. Use to patch in missing sequence lengths"],
				    ["translate-ids", "Translate GFXXXXXXXX to PGF_XXXXXXXX"],
				    ["sort-parallel=i", "Sort threads", { default => 20 } ],
				    ["genome-dir=s", "Directory holding PATRIC genome data", { default => "/vol/patric3/downloads/genomes" }],
				    ["help|h" => "Show this help message."]);

print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 5;

my $fams_file = shift;
my $defs_file = shift;
my $members_file = shift;
my $singleton_defs_file = shift;
my $singleton_members_file = shift;

my $merged_fams_file = $opt->write_merged_fams_file;

my $mem = "10G";

my @defs_sort = ("sort", "-t", "\t", "--parallel", $opt->sort_parallel, "-S", $mem, "-k", "1,1", "-o");
my @members_sort = ("sort", "-t", "\t", "--parallel", $opt->sort_parallel, "-S", $mem, "-k", "1,1", "-k", "6,6", "-k", "5,5n",
		    "-k", "3,3n", "-k", "2,2", "-o");
my @merged_sort = ("sort", "-t", "\t", "--parallel", $opt->sort_parallel, "-S", $mem, "-k", "1,1", "-k", "8,8", "-k", "7,7n",
		   "-k", "5,5n", "-k", "4,4", "-o");

open(F, "<", $fams_file) or die "Cannot read $fams_file: $!";
open(DEFS, "|-", @defs_sort, $defs_file) or die "Cannot write $defs_file: $!";
open(MEMBERS, "|-", @members_sort, $members_file) or die "Cannot write $members_file: $!";
open(SDEFS, "|-", @defs_sort, $singleton_defs_file) or die "Cannot write $singleton_defs_file: $!";
open(SMEMBERS, "|-", @members_sort, $singleton_members_file) or die "Cannot write $singleton_members_file: $!";
if ($merged_fams_file)
{
    open(MERGED, "|-", @merged_sort, $merged_fams_file) or die "Cannot write $merged_fams_file: $!";
}

my $cur_fam;
my %cur_genera;
my $cur_fun;
my $cur_peg_count = 0;
my $cur_txt = '';
my $cur_merged = '';
while (<F>)
{
    chomp;

    s/^GF(\d)/PGF_$1/ if $opt->translate_ids;
    
    my($fam, $localfam_count, $genus_count, $peg, $len, $fun, $localfam, $genus) = split(/\t/);

    if (!$len && $opt->cache_dir)
    {
	$len = compute_length($peg, $opt->cache_dir);
    }

    if ($fam ne $cur_fam)
    {
	if ($cur_fam)
	{
	    my @genera = sort keys %cur_genera;
	    my $defline = join("\t", $cur_fam, $cur_fun, $cur_peg_count, scalar @genera, join(",", @genera)) .  "\n";
	    
	    if ($cur_peg_count == 1)
	    {
		print SDEFS $defline;
		print SMEMBERS $cur_txt;
	    }
	    else
	    {
		print DEFS $defline;
		print MEMBERS $cur_txt;
		print MERGED $cur_merged if $merged_fams_file;
	    }
	}
	$cur_fam = $fam;
	$cur_fun = $fun;
	%cur_genera = ();
	$cur_peg_count = 0;
	$cur_txt = '';
	$cur_merged = '';
    }
    $cur_genera{$genus}++;
    $cur_peg_count++;
    # print MEMBERS join("\t", $fam, $peg, $len, $fun, $localfam, $genus), "\n";
    $cur_txt .= join("\t", $fam, $peg, $len, $fun, $localfam, $genus) . "\n";
    if ($merged_fams_file) {
	$cur_merged .= join("\t", $fam, $localfam_count, $genus_count, $peg, $len, $fun, $localfam, $genus) . "\n";
    }
}

if ($cur_fam)
{
    my @genera = sort keys %cur_genera;
    my $defline = join("\t", $cur_fam, $cur_fun, $cur_peg_count, scalar @genera, join(",", @genera)) .  "\n";
    
    if ($cur_peg_count == 1)
    {
	print SDEFS $defline;
	print SMEMBERS $cur_txt;
    }
    else
    {
	print DEFS $defline;
	print MEMBERS $cur_txt;
	print MERGED $cur_merged if $merged_fams_file;
    }
}
close(DEFS) or die "Error closing $defs_file: $!";
close(MEMBERS) or die "Error closing $members_file: $!";
close(SDEFS) or die "Error closing $singleton_defs_file $!";
close(SMEMBERS) or die "Error closing $singleton_members_file: $!";
close(MERGED);

my %length_cache;
my %length_genome_loaded;

sub compute_length
{
    my($peg, $cache_dir) = @_;

    my($genome, $tax) = $peg =~ /fig\|((\d+)\.\d+)/;

    if (!$length_genome_loaded{$genome})
    {
	my $base = sprintf("%03d", $tax % 1000);
	my $dir = "$cache_dir/$base/$genome";

	#
	# If there's a cache dir, load from there. Otherwise we
	# need to load from the genome FAA file
	#
	my $where;
	if (-d $dir)
	{
	    if (open(my $lens, "<", "$dir/seq_lens.txt"))
	    {
		$where = "$dir/seq_lens.txt";
		while (<$lens>)
		{
		    my ($fid, $len) = /^(\S+)\t(\d+)/;
		    if ($len) {
			$length_cache{$fid} = $len;
		    }
		}
	    }
	    else
	    {
		die "Cannot load $dir/seq_lens.txt: $!";
	    }
	}
	else
	{
	    my $prots = $opt->genome_dir . "/$genome/$genome.PATRIC.faa";
	    $where = $prots;
	    if (open(P, "<", $prots))
	    {
		while (my($fid, $def, $seq) = read_next_fasta_seq(\*P))
		{
		    $length_cache{$fid} = length($seq);
		}
	    }
	    else
	    {
		die "Cannot load $prots: $!";
	    }
	      
	    close(P);
	}

	print STDERR "Loaded $genome length data\n";
	$length_genome_loaded{$genome} = 1;

	if (!$length_cache{$peg})
	{
	    print STDERR "Did not find length of $peg after loading $where\n";
	}
    }
    return $length_cache{$peg};
}
