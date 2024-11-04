#
# Load a set of genomes from the API into a cache directory.
#

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
use Proc::ParallelLoop;
use LPTScheduler;
use PFCache qw(compute_cache_path);

my $token = P3AuthToken->new;
my @auth = (Authorization => $token->token);

my($opt, $usage) = describe_options("%c %o genome-file cache-dir",
				    ["genome-dir=s", "Directory holding PATRIC genome data", { default => "/vol/patric3/downloads/genomes" }],
				    ["solr-url|d=s", "Solr API url", { default => 'https://www.patricbrc.org/api' }],
				    ["parallel|j=i", "Run tasks in parallel", { default => 1 }],
				    ["serial", "Run without any parallelism support"],
				    ["sequence-batch-size=i", "Size of sequence retrieval batches", { default => 2000 }],
				    ["raw-solr", "Server is raw Solr not data API"],
				    ["help|h", "Show this help message"],
				    );

print($usage->text), exit if $opt->help;
die($usage->text) if @ARGV != 2;

my $genome_file = shift;
my $cache_base_dir = shift;

my @genomes = read_file($genome_file);
@genomes or die "Could not read genomes from $genome_file\n";
chomp @genomes;

-d $cache_base_dir or die "Cache dir $cache_base_dir does not exist\n";

if ($opt->serial)
{
    for my $gid (@genomes)
    {
	load_genome_from_api($gid);
    }
    exit;
}

my $sched = LPTScheduler->new($opt->parallel);

$sched->add_work($_, 1) foreach @genomes;

$sched->run(sub {}, sub {
    my($glob, $item) = @_;
    load_genome_from_api($item);
});

sub load_genome_from_api
{
    my($gid) = @_;

    my $prots = $opt->genome_dir . "/$gid/$gid.PATRIC.faa";
    if (-s $prots)
    {
	print STDERR "$gid not needed\n";
	return;
    }

    my $cache_dir = compute_cache_path($cache_base_dir, $gid);

    my $dna_seq = "$cache_dir/dna.fa";
    my $aa_seq = "$cache_dir/aa.fa";
    my $bad_seqs = "$cache_dir/bad_seqs.txt";
    my $trunc = "$cache_dir/trunc_seqs.txt";
    my $gene_names = "$cache_dir/gene_names.txt";
    my $seq_lens = "$cache_dir/seq_lens.txt";
    
    if (-s $dna_seq && -s $aa_seq)
    {
	print STDERR "$gid already cached\n";
	return;
    }

    #
    # Query for contig lengths.
    #
    my %seq_len;
    open(SEQ_LEN, ">", $seq_lens) or die "Cannot write $seq_lens: $!";
    
    my $ua = LWP::UserAgent->new;
    my $block = 25000;
    my $start_idx = 0;
    while (1)
    {
	my @q = (q => "genome_id:$gid",
		 fl => "sequence_type,length,accession",
		 rows => $block,
		 start => $start_idx,
		);
	
	my $url = core_url("genome_sequence");

	my $res = core_post($ua, "genome_sequence", \@q);
	# my $res = $ua->post($url,
	# 		   'Content-type' => 'application/solrquery+x-www-form-urlencoded',
	# 		   'Accept' => 'application/solr+json',
	# 		    @auth,
	# 		    Content => \@q
	# 		  );
	if (!$res->is_success)
	{
	    die "Failed: $url " . $res->status_line;
	}
	
	my $range = $res->header('content-range');
	# print "Range: $range\n";
	my($tstart, $tstop, $tlast) = $range =~ m,(\d+)-(\d+)/(\d+),;
	
	my $r = $res->content;
	my $data = decode_json($r);

	my $items = $data;
	# print Dumper($items);
	for my $item (@$items)
	{
	    $seq_len{$item->{accession}} = $item->{length};
	}
	if ($tstop < $tlast)
	{
	    $start_idx = $tstop;
	}
	else
	{
	    last;
	}
    }
    #
    # Now we can query the features.
    #

    open(GENE_NAMES, ">", $gene_names) or die "Cannot write $gene_names: $!";
    open(TRUNC, ">", $trunc) or die "Cannot write $trunc: $!";
    open(BAD_SEQS, ">", $bad_seqs) or die "Cannot write $bad_seqs: $!";

    $start_idx = 0;
    my %aa_seq_id;
    my %na_seq_id;
    print STDERR "Load features from API for $gid\n";
    while (1)
    {
	my @q = (q => "genome_id:$gid feature_type:CDS annotation:PATRIC",
		 fl => "patric_id,start,end,accession,gene,na_length,aa_sequence_md5,na_sequence_md5",
		 rows => $block,
		 start => $start_idx,
		);
	
	my $url = core_url("genome_feature");

	my $res = core_post($ua, "genome_feature", \@q);
	# my $res = $ua->post($url,
	# 		   'Content-type' => 'application/solrquery+x-www-form-urlencoded',
	# 		   'Accept' => 'application/solr+json',
	# 		    @auth,
	# 		    Content => \@q
	# 		  );
	if (!$res->is_success)
	{
	    die "Failed: " . $res->status_line;
	}
	
	my $range = $res->header('content-range');
	# print "Range: $range\n";
	my($tstart, $tstop, $tlast) = $range =~ m,(\d+)-(\d+)/(\d+),;
	
	my $r = $res->content;
	my $data = decode_json($r);
	
	my $items = $data;
	
	for my $item (@$items)
	{
	    my($id, $start, $end, $acc, $gene, $aa_md5, $na_md5) = @$item{qw(patric_id start end accession gene aa_sequence_md5 na_sequence_md5)};
	    # $id =~ s/\.CDS\./.peg./;
	    # print "id=$id start=$start end=$end ac=$acc gene=$gene $aa_md5 $na_md5\n";
	    $aa_seq_id{$aa_md5} = $id;
	    $na_seq_id{$na_md5} = $id;
	    
	    my $len = $seq_len{$acc};
	    
	    print GENE_NAMES "$id\t$gene\n" if $gene;
	    
	    die Dumper($item) unless defined($len);
	    if ($start < 10 || $end < 10 || $start > $len - 10 || $end > $len - 10)
	    {
		print TRUNC join("\t", $id, $start, $end, $len), "\n";
	    }
	}
	
	if ($tstop < $tlast)
	{
	    $start_idx = $tstop;
	}
	else
	{
	    last;
	}
    }
    close(GENE_NAMES);
    close(TRUNC);
    #
    # And finally the sequences, in batches.
    #
    print STDERR "Load sequence data for $gid from API\n";
    my $batch_size = $opt->sequence_batch_size;
    my @seqs = (keys %na_seq_id, keys %aa_seq_id);
    
    open(AA_SEQS, ">", $aa_seq) or die "Cannot write $aa_seq: $!";
    open(NA_SEQS, ">", $dna_seq) or die "Cannot write $dna_seq: $!";
    while (@seqs)
    {
	my @batch = splice(@seqs, 0, $batch_size);
	
	$start_idx = 0;
	while (1)
	{
	    my @q = (q => "*:*",
		     fq => "md5:(@batch)",
		     fl => "md5,sequence",
		     rows => $block,
		     start => $start_idx,
		    );
	    
	    my $url = core_url("feature_sequence");
	    
	    my $start = time;
	    my $last_print;
	    my $res;
	    while (1)
	    {
		# print STDERR "$url\n";
		$res = core_post($ua, "feature_sequence", \@q);
		# $res = $ua->post($url,
		# 		 'Content-type' => 'application/solrquery+x-www-form-urlencoded',
		# 		 'Accept' => 'application/solr+json',
		# 		 @auth,
		# 		 Content => \@q,
		# 		);
		last if $res->is_success;
		sleep(10);
		if (!$last_print || time - $last_print > 120)
		{
		    print STDERR "Still retrying request $url\n";
		    $last_print = time;
		}
	    }
	    
	    my $range = $res->header('content-range');
	    # print "Range: $range\n";
	    my($tstart, $tstop, $tlast) = $range =~ m,(\d+)-(\d+)/(\d+),;
	    
	    my $r = $res->content;
	    my $data = decode_json($r);
	    
	    my $items = $data;
	    
	    for my $item (@$items)
	    {
		my($md5, $seq) = @$item{qw(md5 sequence)};
		if (my $id = $aa_seq_id{$md5})
		{
		    if ($seq =~ /X{10}/)
		    {
			print BAD_SEQS "Skipping bad sequence $id\n";
		    }
		    else
		    {
			print_alignment_as_fasta(\*AA_SEQS, [$id, '', $seq]);
			print SEQ_LEN join("\t", $id, length($seq)),"\n";

		    }
		}
		else
		{
		    print_alignment_as_fasta(\*NA_SEQS, [$na_seq_id{$md5}, '', $seq]);
		}
	    }
	    
	    if ($tstop < $tlast)
	    {
		$start_idx = $tstop;
	    }
	    else
	    {
		last;
	    }
	}
    }
    close(SEQ_LEN);
    close(DNA_SEQS);
    close(AA_SEQS);
}


sub core_url
{
    my($core) = @_;

    if ($opt->raw_solr)
    {
	return $opt->solr_url . "/$core/select";
    }
    else
    {
	return $opt->solr_url . "/$core";
    }
}

sub core_post
{
    my($ua, $core, $qry) = @_;

    my $url = core_url($core);

    my $type;
    if ($opt->raw_solr)
    {
	$type = "application/x-www-form-urlencoded";
    }
    else
    {
	$type = 'application/solrquery+x-www-form-urlencoded';
    }

    return $ua->post($url,
		     'Content-type' => $type,
		     'Accept' => 'application/json', #'application/solr+json',
		     @auth,
		     Content => $qry,
		    );
}
     
    
