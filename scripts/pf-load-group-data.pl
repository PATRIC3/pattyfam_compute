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

my $token = P3AuthToken->new;
my @auth = (Authorization => $token->token);

my($opt, $usage) = describe_options("%c %o data-dir",
				    ["genome-dir=s", "Directory holding PATRIC genome data", { default => "/vol/patric3/downloads/genomes" }],
				    ["solr-url|d=s", "Solr API url", { default => 'https://www.patricbrc.org/api' }],
				    ["cache-dir=s", "Directory holding API-loaded sequence data cache"],
				    ["help|h", "Show this help message"],
				    );

print($usage->text), exit if $opt->help;
die($usage->text) if @ARGV != 1;

my $group_dir = shift;

-d $group_dir or die "Data directory $group_dir does not exist\n";

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
my $na_seqs_dir = "$group_dir/nr-seqs-dna";

open(SOURCES, ">", "$group_dir/sources") or die "Cannot write $group_dir/sources: $!";
open(GENE_NAMES, ">", "$group_dir/gene.names") or die "Cannot write $group_dir/gene_names: $!";
open(TRUNC, ">", "$group_dir/truncated.genes") or die "Cannot write $group_dir/truncated.genes: $!";
open(BAD_SEQS, ">", "$group_dir/bad.seqs") or die "Cannot write $group_dir/bad.seqs: $!";

my %seq_len;
my $tied = tie %seq_len, 'DB_File', "$group_dir/seq_len.db", O_RDWR | O_CREAT, 0664, $DB_BTREE;
$tied or die "Cannot create $group_dir/seq_len.db: $!";

my %contig_len;
tie %contig_len, 'DB_File', "$group_dir/contig_len.btree", O_RDWR | O_CREAT, 0664, $DB_BTREE
    or die "Cannot tie $group_dir/contig_len.btree:$ !";

my @missed;
for my $gid (@genome_ids)
{
    my $prots = $opt->genome_dir . "/$gid/$gid.PATRIC.faa";
    if (-s $prots && open(P, "<", $prots))
    {
	open(SEQS, ">", "$seqs_dir/$gid") or die "Cannot write $seqs_dir/$gid: $!";
	my $skip;
	while (my($id, $def, $seq) = read_next_fasta_seq(\*P))
	{
	    if ($id =~ /^(fig\|[^|]+)/)
	    {
		if ($seq =~ /X{10}/)
		{
		    print BAD_SEQS "Skipping bad sequence $id from $prots at $.\n";
		}
		    else
		    {
			print_alignment_as_fasta(\*SEQS, [$1, '', $seq]);
			$seq_len{$1} = length($seq);
		    }
	    }
	    else
	    {
		warn "Cannot parse $id from $prots at $.\n";
	    }
	}
	close(SEQS);
	close(P);
	
	#
	# Read the features.tab file to look up the gene names for possible
	# hypothetical family naming.
	#
	# We also scan for contig sizes so we can tag truncated genes. 
	#
	my $tab = $opt->genome_dir . "/$gid/$gid.PATRIC.features.tab";
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
		    print GENE_NAMES "$id\t$gene\n" if $gene;
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
		    print TRUNC join("\t", $id, $start, $end, $len), "\n";
		}
	    }
	    
	    close(TAB);
	}
	else
	{
	    warn "Could not open $tab: $!";
	}
    }
    else
    {
	load_genome_from_api($gid);
    }
    $tied->sync();
    print SOURCES "$seqs_dir/$gid\n";
}

untie %seq_len;
untie %contig_len;
undef $tied;
close(GFILE);
close(GNAME);
close(SOURCES);
close(TRUNC);

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

if (@missed)
{
    print STDERR "Missed the following genomes:\n";
    print STDERR "  $_\n" foreach @missed;
}

sub make_query
{
    my(@list) = @_;
    
    my @q;
    while (@list)
    {
	my($k, $v) = splice(@list, 0, 2);
	push(@q, "$k=$v");
    }
    return join("&", @q);
}

sub load_genome_from_api
{
    my($gid) = @_;

    my($dna_seq, $aa_seq, $bad_seqs, $trunc, $gene_names, $cache_dir);

    if ($opt->cache_dir)
    {
	$cache_dir = $opt->cache_dir . "/$gid";
	make_path($cache_dir);
	
	$dna_seq = "$cache_dir/dna.fa";
	$aa_seq = "$cache_dir/aa.fa";
	$bad_seqs = "$cache_dir/bad_seqs.txt";
	$trunc = "$cache_dir/trunc_seqs.txt";
	$gene_names = "$cache_dir/gene_names.txt";

	if (-s $dna_seq && -s $aa_seq)
	{
	    copy($aa_seq, "$seqs_dir/$gid");
	    copy($dna_seq, "$na_seqs_dir/$gid");
	    copy($bad_seqs, \*BAD_SEQS);
	    copy($gene_names, \*GENE_NAMES);
	    copy($trunc, \*TRUNC);
	    return;
	}
    }
    #
    # Query for contig lengths.
    #
    my %seq_len;
    
    my $ua = LWP::UserAgent->new;
    my $block = 25000;
    my $start_idx = 0;
    while (1)
    {
	my $q = make_query(q => "genome_id:$gid",
			   fl => "sequence_type,length,accession",
			   rows => $block,
			   start => $start_idx,
			  );
	
	my $url = $opt->solr_url . "/genome_sequence/?$q";
	
	# print STDERR "$url\n";
	my $res = $ua->get($url,
			   'Content-type' => 'application/solrquery+x-www-form-urlencoded',
			   'Accept' => 'application/solr+json',
			   @auth
			  );
	if (!$res->is_success)
	{
	    die "Failed: $url " . $res->status_line;
	}
	
	my $range = $res->header('content-range');
	# print "Range: $range\n";
	my($tstart, $tstop, $tlast) = $range =~ m,(\d+)-(\d+)/(\d+),;
	
	my $r = $res->content;
	my $data = decode_json($r);
	
	my $items = $data->{response}->{docs};
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
    $start_idx = 0;
    my %aa_seq_id;
    my %na_seq_id;
    print STDERR "Load features from API for $gid\n";
    while (1)
    {
	my $q = make_query(q => "genome_id:$gid feature_type:CDS annotation:PATRIC",
			   fl => "patric_id,start,end,accession,gene,na_length,aa_sequence_md5,na_sequence_md5",
			   rows => $block,
			   start => $start_idx,
			  );
	
	my $url = $opt->solr_url . "/genome_feature/?$q";
	
	# print STDERR "$url\n";
	my $res = $ua->get($url,
			   'Content-type' => 'application/solrquery+x-www-form-urlencoded',
			   'Accept' => 'application/solr+json',
			   @auth
			  );
	if (!$res->is_success)
	{
	    die "Failed: " . $res->status_line;
	}
	
	my $range = $res->header('content-range');
	# print "Range: $range\n";
	my($tstart, $tstop, $tlast) = $range =~ m,(\d+)-(\d+)/(\d+),;
	
	my $r = $res->content;
	my $data = decode_json($r);
	
	my $items = $data->{response}->{docs};
	
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
    #
    # And finally the sequences, in batches.
    #
    print STDERR "Load sequence data for $gid from API\n";
    my $batch_size = 2000;
    my @seqs = (keys %na_seq_id, keys %aa_seq_id);
    
    open(AA_SEQS, ">", "$seqs_dir/$gid") or die "Cannot write $seqs_dir/$gid: $!";
    open(NA_SEQS, ">", "$na_seqs_dir/$gid") or die "Cannot write $na_seqs_dir/$gid: $!";
    while (@seqs)
    {
	my @batch = splice(@seqs, 0, $batch_size);
	
	$start_idx = 0;
	while (1)
	{
	    my @q = (q => "md5:(@batch)",
		     fl => "md5,sequence",
		     rows => $block,
		     start => $start_idx,
		    );
	    
	    my $url = $opt->solr_url . "/feature_sequence";
	    
	    my $start = time;
	    my $last_print;
	    my $res;
	    while (1)
	    {
		# print STDERR "$url\n";
		$res = $ua->post($url,
				 'Content-type' => 'application/solrquery+x-www-form-urlencoded',
				 'Accept' => 'application/solr+json',
				 @auth,
				 Content => \@q,
				);
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
	    
	    my $items = $data->{response}->{docs};
	    
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
    close(SEQS);
}
