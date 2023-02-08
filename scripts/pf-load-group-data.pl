=head1 NAME

    pf-load-group-data
    
=head1 SYNOPSIS

    pf-load-group-data group-dir

=head1 DESCRIPTION

Load the per-group data for the given group directory. We use the genomes file
to define the genomes to load.
    
=cut

use strict;
use File::Path 'make_path';
use File::Slurp;
use LWP::UserAgent;
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

open(SOURCES, ">", "$group_dir/sources") or die "Cannot write $group_dir/sources: $!";
open(GENE_NAMES, ">", "$group_dir/gene.names") or die "Cannot write $group_dir/gene_names: $!";
open(TRUNC, ">", "$group_dir/truncated.genes") or die "Cannot write $group_dir/truncated.genes: $!";
open(BAD_SEQS, ">", "$group_dir/bad.seqs") or die "Cannot write $group_dir/bad.seqs: $!";

#
# Set up for parallel run. We need to accumulate the four files' data in temps, and then
# write the full file when done.
# The global list will contain four filehandles, one for each of the files above.
#

my $tmpdir = "$group_dir/tmp.$$";
make_path($tmpdir);

my $bootstrap = sub {
    my $src_fh = IO::File->new("$tmpdir/src.$$", "w");
    my $name_fh = IO::File->new("$tmpdir/name.$$", "w");
    my $trunc_fh = IO::File->new("$tmpdir/trunc.$$", "w");
    my $bad_fh = IO::File->new("$tmpdir/bad.$$", "w");
    my $seq_len_fh = IO::File->new("$tmpdir/seq_len.$$", "w");
    return [$src_fh, $name_fh, $trunc_fh, $bad_fh, $seq_len_fh];
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

my %seq_len;
my $seq_len_db = "$group_dir/seq_len.db";
unlink($seq_len_db);
my $tied = tie %seq_len, 'DB_File', $seq_len_db, O_RDWR | O_CREAT, 0664, $DB_BTREE;
$tied or die "Cannot create $seq_len_db: $!";

opendir(D, $tmpdir) or die "Cannot opendir $tmpdir: $!";
while (my $f = readdir(D))
{
    if ($f =~ /^src/)
    {
	copy("$tmpdir/$f", \*SOURCES);
    }
    elsif ($f =~ /^name/)
    {
	copy("$tmpdir/$f", \*GENE_NAMES);
    }
    elsif ($f =~ /^trunc/)
    {
	copy("$tmpdir/$f", \*TRUNC);
    }
    elsif ($f =~ /^bad/)
    {
	copy("$tmpdir/$f", \*BAD_SEQS);
    }
    elsif ($f =~ /^seq_len/)
    {
	open(SL, "<", "$tmpdir/$f") or die "Cannot open $tmpdir/$f: $!";
	while (<SL>)
	{
	    chomp;
	    my($id, $len) = split(/\t/);
	    $seq_len{$id} = $len;
	}
	close(SL);
	$tied->sync();
    }
    unlink("$tmpdir/$f") if -f "$tmpdir/$f";
}
rmdir($tmpdir) or die "Error removing $tmpdir: $!";
untie %seq_len;
undef $tied;
	
close(BAD_SEQS) or die "Error closing: $!";
close(GENE_NAMES) or die "Error closing: $!";
close(SOURCES) or die "Error closing: $!";
close(TRUNC) or die "Error closing: $!";

mkdir("$group_dir/nr") or die "cannot mkdir $group_dir/nr: $!";
my $rc = system("pf-build-nr", "$group_dir/sources", "$group_dir/nr/nr",
		"$group_dir/nr/peg.synonyms", "$group_dir/nr/nr-len.btree", "$group_dir/nr/figids");
if ($rc == 0)
{
    $rc = system("pf-compute-nr-seqs", "--genome-dir", $opt->genome_dir,  $group_dir);
    if ($rc != 0)
    {
	warn "pf-compute-nr-seqs $group_dir failed: $rc\n";
    }
}
else
{
    warn "build_nr $group_dir failed: $rc\n";
}

sub process_genome
{
    my($gid, $fh_list) = @_;
    my $precomp_dir = $opt->genome_dir . "/$gid";
    my $prots = "$precomp_dir/$gid.PATRIC.faa";

    my($src_fh, $name_fh, $trunc_fh, $bad_fh, $seq_len_fh) = @$fh_list;
	
    if (-s $prots)
    {
	load_genome_from_precomputed_data($gid, $precomp_dir, $name_fh, $trunc_fh, $bad_fh, $seq_len_fh);
    }
    else
    {
	load_genome_from_cache($gid, $cache_base_dir, $name_fh, $trunc_fh, $bad_fh, $seq_len_fh);
    }
    print $src_fh "$seqs_dir/$gid\n";
}

sub load_genome_from_precomputed_data
{
    my($gid, $precomp_dir, $name_fh, $trunc_fh, $bad_fh, $seq_len_fh) = @_;

    my $prots = "$precomp_dir/$gid.PATRIC.faa";
    open(P, "<", $prots) or die "Cannot read $precomp_dir/$gid.PATRIC.faa: $!";
    open(SEQS, ">", "$seqs_dir/$gid") or die "Cannot write $seqs_dir/$gid: $!";

    
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
		print $name_fh "$id\t$gene\n" if $gene;
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
	    
	    if ($start < 10 || $end < 10 || $start > $len - 10 || $end > $len - 10)
	    {
		print $trunc_fh join("\t", $id, $start, $end, $len), "\n";
		$truncated{$id} = 1;
	    }
	}
	
	close(TAB);

	while (my($id, $def, $seq) = read_next_fasta_seq(\*P))
	{
	    my $fid;
	    if (($fid) = $id =~ /^(fig\|[^|]+)/)
	    {
		if ($seq =~ /X{10}/)
		{
		    print $bad_fh "Skipping bad sequence $id from $prots at $.\n";
		}
		elsif ($truncated{$fid})
		{
		    print $bad_fh "Skipping truncated $id froM $prots at $.\n";
		}
		else
		{
		    print_alignment_as_fasta(\*SEQS, [$fid, '', $seq]);
		    print $seq_len_fh join("\t", $fid, length($seq)) . "\n";
		}
	    }
	    else
	    {
		warn "Cannot parse $id from $prots at $.\n";
	    }
	}
	close(SEQS);
	close(P);
    }
    else
    {
	warn "Could not open $tab: $!";
    }
}    

sub load_genome_from_cache
{
    my($gid, $cache_base_dir, $name_fh, $trunc_fh, $bad_fh, $seq_len_fh) = @_;

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
    copy($gene_names, $name_fh);
    copy($trunc, $trunc_fh);
    copy($bad_seqs, $bad_fh);
    copy($seq_lens, $seq_len_fh);

    #
    # Copy sequence data.
    #

    copy($aa_seq, "$seqs_dir/$gid");
    copy($dna_seq, "$na_seqs_dir/$gid");
}
