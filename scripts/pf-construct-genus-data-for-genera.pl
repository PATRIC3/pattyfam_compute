=head1 NAME

    pf-construct-genus-data-for-genera

=head1 SYNOPSIS

    pf-construct-genus-data-for-genera [opts] data-dir

=head1 DESCRIPTION


Build the data directory for a family computation based on the list
of genera pulled using the PATRIC data api.

This is the first step in a pattyfam computation. It builds the data directory
hierarchy used for the computation of the local families.

We build a directory for each bacterial/archael genus that has a sufficient
number of genomes.

We also build a directory for all phage genomes, and one for each genus of 
non-phage viruses.

=cut



#
# NEED to pull genus mapping data at this time as well.
#
#" http://localhost:9006/solr/taxonomy/select?q=taxon_rank%3Agenus&fl=taxon_name%2Ctaxon_id&wt=csv&csv.separator=%09&rows=1000000&sort=taxon_name%20asc
    

use strict;
use File::Path 'make_path';
use File::Slurp;
use LWP::UserAgent;
use Getopt::Long::Descriptive;
use List::MoreUtils 'first_index';
use Data::Dumper;
use URI;
use JSON::XS;
use DBI;
#
# We use the pl module here; it is a copy of Proc::ParallelLoop with a memory/fd leak fixed.
#
use pl;
#use Proc::ParallelLoop;
use IPC::Run 'run';
use gjoseqlib;
use DB_File;

my($opt, $usage) = describe_options("%c %o data-dir",
				    ['rank=s', "Use the given taxon rank for grouping (defaults to genus)",
				 		{ default => 'genus' }],
				    ['genus|g=s@', 'Limit to the given genus. May be repeated for multiple genera.'],
				    ['genomes=s', 'Limit to the genomes in this file'],
				    ['bad-genome|b=s@', "A bad genome. Don't use it.", { default => ['340189.4'] }],
				    ["min-genomes|m=i", "Minimum number of genomes required in a genus to build families", { default => 4 }],
				    ["solr-url|d=s", "Solr API url", { default => 'https://www.patricbrc.org/api' }],
				    ["parallel|p=i", "Run with this many procs", { default => 1 }],
				    ["genome-dir=s", "Directory holding PATRIC genome data", { default => "/vol/patric3/downloads/genomes" }],
				    ["help|h", "Show this help message"],
				    );

print($usage->text), exit if $opt->help;
die($usage->text) if @ARGV != 1;

my $base_data_dir = shift;


my %gnames;
my @genome_data = get_genomes($opt, \%gnames);

pareach \@genome_data, sub {
    my($ent) = @_;
    my($genus, $genome_ids) = @$ent;

    my $data_dir = "$base_data_dir/$genus";
    die "Data directory $data_dir already exists\n" if -d $data_dir;

    my $seqs_dir = "$data_dir/Seqs";

    make_path($data_dir, $seqs_dir);

    open(GFILE, ">", "$data_dir/genomes") or die "Cannot write $data_dir/genomes: $!";
    open(GNAME, ">", "$data_dir/genome.names") or die "Cannot write $data_dir/genome.names: $!";
    open(SOURCES, ">", "$data_dir/sources") or die "Cannot write $data_dir/sources: $!";
    open(GENE_NAMES, ">", "$data_dir/gene.names") or die "Cannot write $data_dir/gene_names: $!";
    open(TRUNC, ">", "$data_dir/truncated.genes") or die "Cannot write $data_dir/truncated.genes: $!";

    my %seq_len;
    my $tied = tie %seq_len, 'DB_File', "$data_dir/seq_len.db", O_RDWR | O_CREAT, 0664, $DB_BTREE;
    $tied or die "Cannot create $data_dir/seq_len.db: $!";

    my %contig_len;
    tie %contig_len, 'DB_File', "$data_dir/contig_len.btree", O_RDWR | O_CREAT, 0664, $DB_BTREE
	or die "Cannot tie $data_dir/contig_len.btree:$ !";
    
    my @missed;
    for my $gid (@$genome_ids)
    {
	my $prots = $opt->genome_dir . "/$gid/$gid.PATRIC.faa";
	if (!open(P, "<", $prots))
	{
	    warn "Cannot open prots file $prots: $!";
	    next;
	}
	open(SEQS, ">", "$seqs_dir/$gid") or die "Cannot write $seqs_dir/$gid: $!";
	my $skip;
	while (my($id, $def, $seq) = read_next_fasta_seq(\*P))
	{
	    if ($id =~ /^(fig\|[^|]+)/)
	    {
		if ($seq =~ /X{10}/)
		{
		    warn "Skipping bad sequence $id from $prots at $.\n";
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

	$tied->sync();
	print GFILE "$gid\n";
	print GNAME "$gid\t$gnames{$gid}\n";
	print SOURCES "$seqs_dir/$gid\n";
    }

    untie %seq_len;
    untie %contig_len;
    undef $tied;
    close(GFILE);
    close(GNAME);
    close(SOURCES);
    close(TRUNC);

    mkdir("$data_dir/nr") or die "cannot mkdir $data_dir/nr: $!";
    my $rc = system("build_nr_md5", "$data_dir/sources", "$data_dir/nr/nr",
		    "$data_dir/nr/peg.synonyms", "$data_dir/nr/nr-len.btree", "$data_dir/nr/figids");
    if ($rc == 0)
    {
	$rc = system("pf-compute-nr-seqs", "--genome-dir", $opt->genome_dir,  $data_dir);
	if ($rc != 0)
	{
	    warn "pf-compute-nr-seqs $data_dir failed: $rc\n";
	}
    }
    else
    {
	warn "build_nr $data_dir failed: $rc\n";
    }
    
    
    if (@missed)
    {
	print STDERR "Missed the following genomes:\n";
	print STDERR "  $_\n" foreach @missed;
    }
}, { Max_Workers => $opt->parallel };

sub get_genomes
{
    my($opt, $gnames) = @_;

    # curl 'https://www.patricbrc.org/api/genome/?q=taxon_lineage_ids:33882&rows=250&start=1&http_content-type=application/solrquery+x-www-form-urlencoded&http_accept=application/solr+json'

    # curl 'http://macleod.vbi.vt.edu:8080/solr/genome/select?q=!genus:%22%22&fl=genus,genome_id,genome_name&wt=csv&csv.separator=%09&rows=100000&sort=genus+asc'

    my $ua = LWP::UserAgent->new;

    my $start = 0;
    my $block = 25000;
    my @out;
    my $lg;
    my $cur;

    my $rank = $opt->rank;
    my $have_bad = $opt->bad_genome;
    my %bad_genome;
    if ($have_bad)
    {
	$bad_genome{$_} = 1 foreach @{$opt->bad_genome};
    }
    # print Dumper(\%bad_genome);

    my $limit_genomes;
    if ($opt->genomes)
    {
	$limit_genomes = {};
	open(G, "<", $opt->genomes) or die "Cannot read " . $opt->genomes . ": $!";
	while (<G>)
	{
	    if (/(\d+\.\d+)/)
	    {
		$limit_genomes->{$1} = 1;
	    }
	    else
	    {
		die "Cannot parse line $. in " . $opt->genomes;
	    }
	}
	close(G);
    }

    my %by_genus;
    #
    # We hardcode a policy change here where anything in viruses gets
    # stuffed in the Viruses directory.
    #
    # To make this workable, collect up all the data and postprocess.
    #

    while (1)
    {
	my $q = make_query(q => "$rank:* public:1",
			   fl => "$rank,genome_id,genome_name,domain,kingdom,genome_quality,genome_status",
			   rows => $block,
			   start => $start,
			   sort => "$rank asc",
			  );

	my $url = $opt->solr_url . "/genome/?$q";
	
	print STDERR "$url\n";
	my $res = $ua->get($url,
			   'Content-type' => 'application/solrquery+x-www-form-urlencoded',
			   'Accept' => 'application/solr+json',
			  );
	if (!$res->is_success)
	{
	    die "Failed: " . $res->status_line;
	}

	my $range = $res->header('content-range');
	print "Range: $range\n";
	my($tstart, $tstop, $tlast) = $range =~ m,(\d+)-(\d+)/(\d+),;

	my $r = $res->content;
	my $data = decode_json($r);

	# print Dumper($data);

	my $items = $data->{response}->{docs};

	my $limit_genera = defined($opt->genus);

	for my $item (@$items)
	{
	    my($genus, $gid, $name, $kingdom, $quality, $status) = @$item{$rank, qw(genome_id genome_name kingdom genome_quality genome_status)};

	    next if $limit_genomes && !$limit_genomes->{$gid};

	    if ($quality ne 'Good' && $kingdom ne 'Viruses')
	    {
		warn "Skipping $gid due to quality=$quality and kingdom=$kingdom\n";
		next;
	    }

	    $gnames->{$gid} = $name;
	    if ($have_bad && $bad_genome{$gid})
	    {
		print STDERR "Skipping bad genome $gid\n";
		print STDERR Dumper($item);
		next;
	    }
	    if ($limit_genera && $kingdom ne 'Viruses')
	    {
		my $ok;
		for my $x (@{$opt->genus})
		{
		    if ($genus =~ /^$x$/)
		    {
			$ok = 1;
			last;
		    }
		}
		next unless $ok;
	    }
	    
	    next if $genus eq '$genus';
	    next if $genus eq '""';
	    next if $genus eq '';

	    my $ref_genus = ($kingdom eq 'Viruses') ? 'Viruses' : $genus;

	    push(@{$by_genus{$ref_genus}}, $gid);
	}

	if ($tstop < $tlast)
	{
	    $start = $tstop;
	}
	else
	{
	    last;
	}
    }

    for my $genus (sort keys %by_genus)
    {
	my $list = $by_genus{$genus};

	if (@$list >= $opt->min_genomes)
	{
	    push(@out, [$genus, $list]);
	}
    }
	
    return @out;
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
