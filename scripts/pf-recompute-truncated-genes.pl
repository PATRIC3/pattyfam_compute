=head1 NAME

    pf-recompute-truncated-genes
    
=head1 SYNOPSIS

    pf-recompute-truncated-genes group-dir

=head1 DESCRIPTION

Hack of load-group-data that only recomputes truncated genes with a different threshold.
    
=cut

use lib '/home/olson/perl5/lib/perl5';

use strict;
use File::Path 'make_path';
use File::Slurp;
use LWP::UserAgent::Determined;
use Getopt::Long::Descriptive;
use List::MoreUtils 'first_index';
use Data::Dumper;
use URI;
use URI::Escape;
use JSON::XS;
use DBI;
use IPC::Run 'run';
use gjoseqlib;
use DB_File;
use P3AuthToken;
use PFCache qw(compute_cache_path);
use File::Copy qw(copy);
use IO::File;
use LPTScheduler;

my $token = P3AuthToken->new;
my @auth = (Authorization => $token->token);

my($opt, $usage) = describe_options("%c %o data-dir cache-dir",
				    ["data-dir is the genus directory we are loading with data"],
				    ["cache-dir is the cache of API-generated data used as a fallback, created by pf-load-genomes-to-cache"],
				    [],
				    ["genome-dir=s", "Directory holding PATRIC genome data", { default => "/vol/patric3/downloads/genomes" }],
				    ["solr-url|d=s", "Solr API url", { default => 'https://www.patricbrc.org/api' }],
				    ["parallel|j=i", "Run tasks in parallel", { default => 1 }],
				    ["help|h", "Show this help message"],
				    );

print($usage->text), exit if $opt->help;
$usage->die() if @ARGV != 2;

my $group_dir = shift;
my $cache_base_dir = shift;

-d $group_dir or die "Data directory $group_dir does not exist\n";
-d $cache_base_dir or die "Cache directory $cache_base_dir does not exist\n";

my $json = JSON::XS->new->pretty(1);

my @genome_ids;
if (-s "$group_dir/active.genomes")
{
    @genome_ids = read_file("$group_dir/active.genomes");
}
else
{
    @genome_ids = read_file("$group_dir/genomes");
}

@genome_ids or die "Could not read genome IDs for $group_dir\n";
chomp @genome_ids;

print STDERR "Loading " . scalar(@genome_ids) . " genomes for $group_dir\n";

my $seqs_dir = "$group_dir/Seqs";
my $na_seqs_dir = "$group_dir/NASeqs";

open(TRUNC, ">", "$group_dir/truncated.genes.new") or die "Cannot write $group_dir/truncated.genes.new: $!";


#
# Set up for parallel run. We need to accumulate the four files' data in temps, and then
# write the full file when done.
# The global list will contain four filehandles, one for each of the files above.
#

my $tmpdir = "$group_dir/tmp.$$";
make_path($tmpdir);

my $bootstrap = sub {
    my $trunc_fh = IO::File->new("$tmpdir/trunc.$$", "w");
    return [$trunc_fh];
};

my $sched = LPTScheduler->new($opt->parallel);

$sched->add_work($_, 1) foreach @genome_ids;
$sched->run($bootstrap, sub {
    my($fh_list, $gid) = @_;

    process_genome($gid, $fh_list);
});

#
# Now collect our output and write.
#

opendir(D, $tmpdir) or die "Cannot opendir $tmpdir: $!";
while (my $f = readdir(D))
{
    if ($f =~ /^trunc/)
    {
	copy("$tmpdir/$f", \*TRUNC);
    }
    unlink("$tmpdir/$f") if -f "$tmpdir/$f";
}
rmdir($tmpdir) or die "Error removing $tmpdir: $!";
	
close(TRUNC) or die "Error closing: $!";


sub process_genome
{
    my($gid, $fh_list) = @_;
    my $precomp_dir = $opt->genome_dir . "/$gid";
    my $prots = "$precomp_dir/$gid.PATRIC.faa";

    my($trunc_fh) = @$fh_list;
	
    if (-s $prots)
    {
	load_genome_from_precomputed_data($gid, $precomp_dir, $trunc_fh);
    }
    else
    {
	load_genome_from_cache($gid, $cache_base_dir, $trunc_fh);
    }
}

sub load_genome_from_precomputed_data
{
    my($gid, $precomp_dir, $trunc_fh) = @_;

    my $prots = "$precomp_dir/$gid.PATRIC.faa";
    open(P, "<", $prots) or die "Cannot read $precomp_dir/$gid.PATRIC.faa: $!";
    
    #
    # Read the features.tab file to look up the gene names for possible
    # hypothetical family naming.
    #
    # We also scan for contig sizes so we can tag truncated genes. 
    #
    my $tab = "$precomp_dir/$gid.PATRIC.features.tab";
    my %contig_len;
    my %truncated;

    if (open(TAB, "<", $tab))
    {
	my $key = <TAB>;
	chomp $key;
	my @hdrs = split(/\t/, $key);
	my $gene_idx = first_index { $_ eq 'gene' } @hdrs;
	my $id_idx = first_index { $_ eq 'patric_id' } @hdrs;
	my $ftype_idx = first_index { $_ eq 'feature_type' } @hdrs;
	my $acc_idx = first_index { $_ eq 'accession' } @hdrs;
	my $start_idx = first_index { $_ eq 'start' } @hdrs;
	my $end_idx = first_index { $_ eq 'end' } @hdrs;
	my $len_idx = first_index { $_ eq 'na_length' } @hdrs;
	
	if ($gene_idx < 0)
	{
	    warn "Could not find gene column in $tab. Headers are '$key'\n";
	}
	#
	# Scan once for contig lengths and gene nam
	while (<TAB>)
	{
	    chomp;
	    my @vals = split(/\t/);
	    my $id = $vals[$id_idx];
	    if (defined($gene_idx))
	    {
		my $gene = $vals[$gene_idx];
	    }
	    if ($vals[$ftype_idx] eq 'source')
	    {
		$contig_len{$vals[$acc_idx]} = $vals[$len_idx];
	    }
	}
	seek(TAB, 0, 0);
	my $tried;
	while (<TAB>)
	{
	    chomp;
	    my @vals = split(/\t/);
	    my($id, $type, $start, $end, $acc) = @vals[$id_idx, $ftype_idx, $start_idx, $end_idx, $acc_idx];
	    
	    next unless $type eq 'CDS';
	    
	    my $len = $contig_len{$acc};
	    
	    if (!$len && !$tried)
	    {
		$tried = 1;
		#
		# Length missing means we need to read the contigs.
		#
		my $contigs = $opt->genome_dir . "/$gid/$gid.fna";
		
		#
		# Try the contigs cache dir created during ANI run.
		#
		if (! -s $contigs)
		{
		    $contigs = "$group_dir/contigs/$gid.fna";
		}
		if (open(C, "<", $contigs) )
		{
		    while (my($id, $def, $seq) = read_next_fasta_seq(\*C))
		    {
			$contig_len{$id} = length($seq);
		    }
		    close(C);
		}
		else
		{
		    warn "Cannot open $contigs: $!";
		}
		$len = $contig_len{$acc};
		next unless $len;
	    }
	    
	    if ($start < 100 || $end < 100 || $start > $len - 100 || $end > $len - 100)
	    {
		print $trunc_fh join("\t", $id, $start, $end, $len), "\n";
		$truncated{$id} = 1;
	    }
	}
	
	close(TAB);
    }
    else
    {
	warn "Could not open $tab: $!";
    }
}    

sub load_genome_from_cache
{
    my($gid, $cache_base_dir, $trunc_fh) = @_;

    my $cache_dir = compute_cache_path($cache_base_dir, $gid);

    -d $cache_dir or die "Cache $cache_dir not available for $gid";

    my $dna_seq = "$cache_dir/dna.fa";
    my $aa_seq = "$cache_dir/aa.fa";
    my $bad_seqs = "$cache_dir/bad_seqs.txt";
    my $trunc = "$cache_dir/trunc_seqs.txt";
    my $gene_names = "$cache_dir/gene_names.txt";
    my $seq_lens = "$cache_dir/seq_lens.txt";
    
    #
    # Copy gene name and truncated data.
    #
    copy($trunc, $trunc_fh);
}
