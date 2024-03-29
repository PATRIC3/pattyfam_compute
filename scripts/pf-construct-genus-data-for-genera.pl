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
number of genomes. The directories are names Genus-GenusTaxId because there the
genus names are not unique across the taxonomic tree (e.g. Ponticoccus is both 983507
519422).

We also build a directory for all phage genomes, and one for each genus of 
non-phage viruses.

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
#
# We use the pl module here; it is a copy of Proc::ParallelLoop with a memory/fd leak fixed.
#
use pl;
#use Proc::ParallelLoop;
use IPC::Run 'run';
use gjoseqlib;
use DB_File;


my @auth;
if (open(P, "<", "$ENV{HOME}/.patric_token"))
{
    my $t = <P>;
    chomp $t;
    close(P);
    @auth = ("Authorization", $t);
}

my($opt, $usage) = describe_options("%c %o data-dir",
				    ['rank=s', "Use the given taxon rank for grouping (defaults to genus)",
				 		{ default => 'genus' }],
				    ['genus|g=s@', 'Limit to the given genus. May be repeated for multiple genera.'],
				    ['genomes=s', 'Limit to the genomes in this file'],
				    ['phages!', 'Collect phages', { default => 1 }],
				    ['check-quality!', 'Only use genomes with good quality', { default => 1 }],
				    ['bad-genome|b=s@', "A bad genome. Don't use it.", { default => ['340189.4'] }],
				    ["min-cds=i", "Minimum number of CDS features required to build families"],
				    ["min-genomes|m=i", "Minimum number of genomes required in a genus to build families", { default => 4 }],
				    ["max-genomes|M=i", "Maximum number of genomes to be placed in a genus", { default => 8000 }],
				    ["solr-url|d=s", "Solr API url", { default => 'https://www.patricbrc.org/api' }],
				    ["parallel|p=i", "Run with this many procs", { default => 1 }],
				    ["genome-dir=s", "Directory holding PATRIC genome data", { default => "/vol/patric3/downloads/genomes" }],
				    ["help|h", "Show this help message"],
				    );

print($usage->text), exit if $opt->help;
die($usage->text) if @ARGV != 1;

my $base_data_dir = shift;

-d $base_data_dir or die "Data directory $base_data_dir does not exist\n";

my $json = JSON::XS->new->pretty(1);

my %gnames;
my ($genome_data, $all_genomes) = get_genomes($opt, \%gnames);

pareach $genome_data, sub {
    my($ent) = @_;
    my($genus, $genome_ids, $genus_taxon) = @$ent;

    my $data_dir = "$base_data_dir/$genus";
    die "Data directory $data_dir already exists\n" if -d $data_dir;

    my $seqs_dir = "$data_dir/Seqs";
    my $na_seqs_dir = "$data_dir/nr-seqs-dna";

    make_path($data_dir, $seqs_dir, $na_seqs_dir);

    if ($genus_taxon)
    {
	open(P, ">", "$data_dir/GENUS_TAXON_ID") or die "cannot write $data_dir/GENUS_TAXON_ID: $!";
	print P "$genus_taxon\n";
	close(P);
    }

    open(GFILE, ">", "$data_dir/genomes") or die "Cannot write $data_dir/genomes: $!";
    open(GNAME, ">", "$data_dir/genome.names") or die "Cannot write $data_dir/genome.names: $!";
    open(SOURCES, ">", "$data_dir/sources") or die "Cannot write $data_dir/sources: $!";
    open(GENE_NAMES, ">", "$data_dir/gene.names") or die "Cannot write $data_dir/gene_names: $!";
    open(TRUNC, ">", "$data_dir/truncated.genes") or die "Cannot write $data_dir/truncated.genes: $!";
    open(BAD_SEQS, ">", "$data_dir/bad.seqs") or die "Cannot write $data_dir/bad.seqs: $!";

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
	if (open(P, "<", $prots))
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
	    my $batch_size = 200;
	    my @seqs = (keys %na_seq_id, keys %aa_seq_id);
	    
	    open(AA_SEQS, ">", "$seqs_dir/$gid") or die "Cannot write $seqs_dir/$gid: $!";
	    open(NA_SEQS, ">", "$na_seqs_dir/$gid") or die "Cannot write $na_seqs_dir/$gid: $!";
	    while (@seqs)
	    {
		my @batch = splice(@seqs, 0, $batch_size);
		
		$start_idx = 0;
		while (1)
		{
		    my $q = make_query(q => "md5:(@batch)",
				       fl => "md5,sequence",
				       rows => $block,
				       start => $start_idx,
				      );
		
		    my $url = $opt->solr_url . "/feature_sequence/?$q";
		    
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

    #
    # Move this to the local family compute script so it can be run in parallel.
    #
    # mkdir("$data_dir/nr") or die "cannot mkdir $data_dir/nr: $!";
    # my $rc = system("pf-build-nr", "$data_dir/sources", "$data_dir/nr/nr",
    # 		    "$data_dir/nr/peg.synonyms", "$data_dir/nr/nr-len.btree", "$data_dir/nr/figids");
    # if ($rc == 0)
    # {
    # 	$rc = system("pf-compute-nr-seqs", "--genome-dir", $opt->genome_dir,  $data_dir);
    # 	if ($rc != 0)
    # 	{
    # 	    warn "pf-compute-nr-seqs $data_dir failed: $rc\n";
    # 	}
    # }
    # else
    # {
    # 	warn "build_nr $data_dir failed: $rc\n";
    # }
    
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
    my $gquery = "";
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
	$gquery = "genome_id:(" . join(" ", keys %$limit_genomes) . ")";
    }

    my %by_genus;

    #
    # We support viruses and phages by creating a single data directory Phages
    # that incorporates all genomes in the NCBI division 'Phages'.
    #
    # Other viruses go into their own genus directories.
    #

    #
    # We begin by pulling all of the taxonomic data for genomes at the genus level.
    #

    my %genus_data;

    my @genera_query;
    if ($opt->genus)
    {
	@genera_query = (fq => "$rank:(" . join(" OR ", map { "\"$_\"" }  @{$opt->genus}) . ")" );
    }

    while (1)
    {
	my $q = make_query(q => "taxon_rank:$rank",
			   fl => "taxon_id,taxon_name,genetic_code,division",
			   rows => $block,
			   start => $start,
			  );

	my $url = $opt->solr_url . "/taxonomy/?$q";

	print STDERR "$url\n";
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
	print "Range: $range\n";
	my($tstart, $tstop, $tlast) = $range =~ m,(\d+)-(\d+)/(\d+),;

	my $r = $res->content;
	my $data = decode_json($r);

	my $items = $data->{response}->{docs};

	for my $item (@$items)
	{
	    $genus_data{$item->{taxon_id}} = $item;
	}

	#print Dumper($items);

	if ($tstop < $tlast)
	{
	    $start = $tstop;
	}
	else
	{
	    last;
	}
    }
    open(G, ">", "$base_data_dir/genus_data.json") or die "Cannot write $base_data_dir/genus_data.json: $!";
    print G $json->encode(\%genus_data);
    
    #
    # Now scan the genomes. With the genus data collected above, we can process and classify
    # the genomes in a single pass
    #

    my %gmap;
    my %all_genomes;

    my $start_idx = 0;
    while (1)
    {
	my $q = make_query(defined($gquery) ? (q => "public:* $gquery") : (q => "$rank:* public:1"),
			   fl => "$rank,genome_id,genome_name,domain,kingdom,genome_quality,genome_status,genetic_code,taxon_lineage_ids,owner,public,contigs,species,genome_quality_flags,cds",
			   @genera_query,
			   # ($opt->min_cds ? ("fq" => "cds:[" . $opt->min_cds . " TO *]") : ()),
			   rows => $block,
			   start => $start_idx,
			   sort => "$rank asc",
			  );

	print "Q=$q\n";
	$q = uri_escape($q, '()" :');
	print "Q=$q\n";
	my $url = $opt->solr_url . "/genome/?$q";
	
	print STDERR "$url\n";
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
	print "Range: $range\n";
	my($tstart, $tstop, $tlast) = $range =~ m,(\d+)-(\d+)/(\d+),;

	my $r = $res->content;
	my $data = decode_json($r);

	# print Dumper($data);

	my $items = $data->{response}->{docs};

	my $limit_genera = defined($opt->genus);
	# print STDERR Dumper($items);

	for my $item (@$items)
	{
	    my($genus, $gid, $name, $kingdom, $quality, $status, $genetic_code, $lineage, $public, $owner, $quality_flags, $cds) = @$item{$rank, qw(genome_id genome_name kingdom genome_quality genome_status genetic_code taxon_lineage_ids public owner genome_quality_flags cds)};

	    if ($opt->min_cds && $cds < $opt->min_cds)
	    {
		warn "Skip $gid $name due to $cds < " . $opt->min_cds . "\n";
		next;
	    }

	    next unless $public || $owner eq 'sars2_bvbrc@patricbrc.org';
#	    if ($owner eq 'sars2_bvbrc@patricbrc.org')
#	    {
#		print Dumper($item);
	    #	    }

	    # print "$genus\t$gid\t$name\n";

	    $all_genomes{$gid} = $item;

	    next if $limit_genomes && !$limit_genomes->{$gid};

	    if (($opt->check_quality && $quality ne 'Good') && $kingdom ne 'Viruses')
	    {
		#
		# If the only reason it's poor quality is that it is short, let it through.
		#
		if (ref($quality_flags) && @$quality_flags == 1 && $quality_flags->[0] eq 'Genome too short')
		{
		    # warn "allowing short genome\n";
		}
		else
		{
		    # warn "Skipping $gid $name due to quality=$quality and kingdom=$kingdom\n";
		    next;
		}
	    }

	    $gnames->{$gid} = $name;
	    if ($have_bad && $bad_genome{$gid})
	    {
		print STDERR "Skipping bad genome $gid\n";
		print STDERR Dumper($item);
		next;
	    }
	    #
	    # We used to have kingdom ne Viruses; that ended up collecting
	    # all the viruses even if we limited genera. Seems wrong.
	    #
	    if ($limit_genera)
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

	    my $gdata;
	    for my $tid (@{$lineage})
	    {
		$gdata = $genus_data{$tid};
		last if $gdata;
	    }
	    if (!$gdata)
	    {
		warn "Skipping $gid $name: no genus data \n";
		next;
	    }

	    my $genus_taxid = $gdata->{taxon_id};
	    my $ref_genus = "$genus-$genus_taxid";
	    if ($gdata->{division} eq 'Phages' && $opt->phages)
	    {
		$ref_genus = 'Phages';
	    }
	    else
	    {
		if (exists($gmap{$ref_genus}))
		{
		    if ($gmap{$ref_genus} ne $gdata->{taxon_id})
		    {
			die "$genus is not unique $gmap{$ref_genus} $gdata->{taxon_id}";
		    }
		}
		$gmap{$ref_genus} = $gdata->{taxon_id};
	    }
	    push(@{$by_genus{$ref_genus}}, $gid);
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

    for my $genus (sort keys %by_genus)
    {
	my $list = $by_genus{$genus};

	next if (@$list < $opt->min_genomes);

	#
	# If the genus is too large, we need to shrink by choosing
	# genomes with the smallest number of contigs.
	#
	
	if (@$list > $opt->max_genomes)
	{
	    print "Limiting $genus to " . $opt->max_genomes . "\n";
	    my %by_species;
	    my @sorted = map { $all_genomes{$_} } sort { $all_genomes{$a}->{contigs} <=> $all_genomes{$b}->{contigs} } @$list;
	    for my $s (@sorted)
	    {
		push(@{$by_species{$s->{species}}}, $s);
		# print join("\t", $s, $all_genomes{$s}->{contigs}, $all_genomes{$s}->{species}, $all_genomes{$s}->{genome_name}), "\n";
	    }
	    my @new_list;
	OUTER:
	    while (keys %by_species)
	    {
		for my $species (keys %by_species)
		{
		    my $list = $by_species{$species};
		    my $ent = shift @$list;
		    push(@new_list, $ent);
		    if (@$list == 0)
		    {
			delete $by_species{$species};
		    }
		    last OUTER if @new_list >= $opt->max_genomes;
		}
	    }
	    # for my $ent (@new_list)
	    # {
	    # 	print join("\t", $ent->{genome_id}, $ent->{contigs}, $ent->{species},
	    # 		   $ent->{genome_name}), "\n";
	    # }
	    @$list = map { $_->{genome_id} } @new_list;
	}
	
	push(@out, [$genus, $list, $gmap{$genus}]);
    }

    return \@out, \%all_genomes;
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
