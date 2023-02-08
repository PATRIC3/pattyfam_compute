=head1 NAME

    pf-initialize-data

=head1 SYNOPSIS

    pf-initialize-data [opts] data-dir

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

We do not retrieve the actual genomic data to be computed on. That is left
for pf-load-group-data.    
    

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
use P3DataAPI;


#
# We use the pl module here; it is a copy of Proc::ParallelLoop with a memory/fd leak fixed.
#
use pl;
#use Proc::ParallelLoop;
use IPC::Run 'run';
use gjoseqlib;
use DB_File;

my $api = P3DataAPI->new;

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
				    ['genus-file|G=s', 'Limit to the genera listed in this file.'],
				    ['genomes=s', 'Limit to the genomes in this file.'],
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

my %bad_genome;
my $have_bad;
if ($opt->bad_genome)
{
    $have_bad = 1;
    $bad_genome{$_} = 1 foreach @{$opt->bad_genome};
}

my $rank = $opt->rank;

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

#
# Pull taxonomic data.
#
my %genus_data;

$api->solr_query_list('taxonomy',
		      [q => "taxon_rank:$rank",
		       fl => "taxon_id,taxon_name,genetic_code,division",
		       ],
		      undef,
		      sub {
			  my($doc) = @_;
			  for my $item (@{$doc->{response}->{docs}})
			  {
			      $genus_data{$item->{taxon_id}} = $item;
			  }
		      });
open(GDATA, ">", "$base_data_dir/genus_data.json") or die "Cannot write $base_data_dir/genus_data.json: $!";
print GDATA $json->encode(\%genus_data);
close(GDATA);


my %gnames;
my ($genome_data, $all_genomes) = get_genomes($opt, \%gnames);

for my $ent (@$genome_data)
{
    my($genus, $genome_ids, $genus_taxon) = @$ent;

    my $data_dir = "$base_data_dir/$genus";
    die "Data directory $data_dir already exists\n" if -d $data_dir;

    my $seqs_dir = "$data_dir/Seqs";
    my $na_seqs_dir = "$data_dir/NASeqs";

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
    
    my @missed;
    for my $gid (@$genome_ids)
    {
	print GFILE "$gid\n";
	print GNAME "$gid\t$gnames{$gid}\n";
	print SOURCES "$seqs_dir/$gid\n";
    }

    close(GFILE);
    close(GNAME);
    close(SOURCES);

    if (@missed)
    {
	print STDERR "Missed the following genomes:\n";
	print STDERR "  $_\n" foreach @missed;
    }
}

sub get_genomes
{
    my($opt, $gnames) = @_;

    my @out;
    my %by_genus;

    #
    # We support viruses and phages by creating a single data directory Phages
    # that incorporates all genomes in the NCBI division 'Phages'.
    #
    # Other viruses go into their own genus directories.
    #

    my @genus_list;
    if ($opt->genus)
    {
	@genus_list = @{$opt->genus};
    }
    if ($opt->genus_file)
    {
	my @x = read_file($opt->genus_file);
	chomp @x;
	s/^\s+// foreach @x;
	s/\s+$// foreach @x;
	push(@genus_list, @x);
    }
    my @genera_query = (fq => "$rank:(" . join(" OR ", map { "\"$_\"" }  @genus_list) . ")" );

    #
    # Scan the genomes. With the genus data collected above, we can process and classify
    # the genomes in a single pass
    #

    my %gmap;
    my %all_genomes;

    my @query;
    if ($gquery ne '')
    {
	push(@query, q => "public:* $gquery");
    }
    else
    {
	push(@query, q => "$rank:* public:1");
    }
    
    push(@query,
	 fl => "$rank,genome_id,genome_name,domain,kingdom,genome_quality,genome_status,genetic_code,taxon_lineage_ids,owner,public,contigs,species,genome_quality_flags,cds,taxon_lineage_names",
	 @genera_query,
	 sort => "$rank asc",
	 fq => "superkingdom:(Bacteria OR Archaea)",
	 fq => "$rank:[* TO *]",
	 fq => "NOT $rank:null",
	 fq => "NOT $rank:\"\"",
	);
    print STDERR "Qry: " . Dumper(@query);

    $api->solr_query_list('genome', \@query, undef, sub {
	my($doc) = @_;

	my $items = $doc->{response}->{docs};
	print "$items->[0]->{genome_name}\n";
	for my $item (@$items)
	{
	    my($genus, $gid, $name, $kingdom, $quality, $status, $genetic_code, $lineage, $public, $owner, $quality_flags, $cds, $lineage_names) =
		@$item{$rank, qw(genome_id genome_name kingdom genome_quality genome_status genetic_code taxon_lineage_ids public owner genome_quality_flags cds taxon_lineage_names)};

	    if (!basic_validity_checks($item))
	    {
		warn "Failed basic validity: " . Dumper($item);
		next;
	    }

	    $all_genomes{$gid} = $item;

	    if (($opt->check_quality && $quality ne 'Good') && $kingdom ne 'Viruses')
	    {
		#
		# If the only reason it's poor quality is that it is short, let it through.
		#
		if (ref($quality_flags) && @$quality_flags == 1 && $quality_flags->[0] eq 'Genome too short')
		{
		    warn "allowing short genome\n";
		}
		else
		{
		    warn "Skipping $gid $name due to quality=$quality and kingdom=$kingdom\n";
		    next;
		}
	    }

	    $gnames->{$gid} = $name;

	    my $gdata;
	    for my $tid (@{$lineage})
	    {
		$gdata = $genus_data{$tid};
		last if $gdata;
	    }
	    if (!$gdata)
	    {
		warn "Skipping $gid $name: no genus data (@$lineage_names)\n" . Dumper($item);
		next;
	    }

	    my $genus_taxid = $gdata->{taxon_id};
	    my $ref_genus = "$genus-$genus_taxid";
	    if (exists($gmap{$ref_genus}))
	    {
		if ($gmap{$ref_genus} ne $gdata->{taxon_id})
		{
		    die "$genus is not unique $gmap{$ref_genus} $gdata->{taxon_id}";
		}
	    }
	    $gmap{$ref_genus} = $gdata->{taxon_id};
	    push(@{$by_genus{$ref_genus}}, $gid);
	}
    });

    for my $genus (sort keys %by_genus)
    {
	my $list = $by_genus{$genus};

	next if (@$list < $opt->min_genomes);

	#
	# If the genus is too large, we need to shrink by choosing
	# genomes with the smallest number of contigs.
	#
	
	if ($opt->max_genomes && @$list > $opt->max_genomes)
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

=head3

Retrieve phage genomes. For our purposes we are using all genomes in the Viruses
superkingdom that have "phage" in the name.

=cut

sub get_phage_genomes
{
    my($opt, $gnames) = @_;

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

    #
    # Now scan the genomes. With the genus data collected above, we can process and classify
    # the genomes in a single pass
    #

    my %gmap;
    my %all_genomes;

    my @query;
    if (defined($gquery))
    {
	push(@query, q => "public:* $gquery");
    }
    else
    {
	push(@query, q => "$rank:* public:1");
    }
    push(@query, fq => "superkingdom:Viruses", fq => "genome_name:phage");
    push(@query, fl => "$rank,genome_id,genome_name,domain,kingdom,genome_quality,genome_status,genetic_code,taxon_lineage_ids,owner,public,contigs,species,genome_quality_flags,cds");
    push(@query, sort => "$rank asc");
    print "Query: @query\n";

    $api->solr_query_list('genome', \@query, undef, sub {
	my($doc) = @_;
	for my $item (@{$doc->{response}->{docs}})
	{
	    my($genus, $gid, $name, $kingdom, $quality, $status, $genetic_code, $lineage, $public, $owner, $quality_flags, $cds) =
		@$item{$rank, qw(genome_id genome_name kingdom genome_quality genome_status genetic_code taxon_lineage_ids public owner genome_quality_flags cds)};

	    if (!basic_phage_validity_checks($item))
	    {
		# warn "Failed basic validity: " . Dumper($item);
		next;
	    }

	    $all_genomes{$gid} = $item;
	    
	    $gnames->{$gid} = $name;
	    
	    my $ref_genus = 'Phages';
	    push(@{$by_genus{$ref_genus}}, $gid);
	}
	return undef;
    });

    my @out;
    for my $genus (sort keys %by_genus)
    {
	my $list = $by_genus{$genus};

	next if (@$list < $opt->min_genomes);

	#
	# We take all phages.
	#
	
	push(@out, [$genus, $list, $gmap{$genus}]);
    }

    return \@out, \%all_genomes;

}

sub basic_validity_checks
{
    my($item) = @_;
    my($genus, $gid, $name, $kingdom, $quality, $status, $genetic_code, $lineage, $public, $owner, $quality_flags, $cds) =
	@$item{$rank, qw(genome_id genome_name kingdom genome_quality genome_status genetic_code taxon_lineage_ids public owner genome_quality_flags cds)};

    if ($opt->min_cds && $cds < $opt->min_cds)
    {
	warn "Skip $item->{genome_id} $item->{genome_name} due to $item->{cds} < " . $opt->min_cds . "\n";
	return 0
	}

    return 0 unless $public || $owner eq 'sars2_bvbrc@patricbrc.org';
	    
    return 0 if $limit_genomes && !$limit_genomes->{$gid};
    
    return 0 if $have_bad && $bad_genome{$gid};
    
    return 0 if $genus eq '$genus' || $genus eq '""' || $genus eq '';

    return 1;
}

#
# Phage checks do not include genome limit or genus name checks.
#
sub basic_phage_validity_checks
{
    my($item) = @_;
    my($genus, $gid, $name, $kingdom, $quality, $status, $genetic_code, $lineage, $public, $owner, $quality_flags, $cds) =
	@$item{$rank, qw(genome_id genome_name kingdom genome_quality genome_status genetic_code taxon_lineage_ids public owner genome_quality_flags cds)};

    if ($opt->min_cds && $cds < $opt->min_cds)
    {
	warn "Skip $item->{genome_id} $item->{genome_name} due to $item->{cds} < " . $opt->min_cds . "\n";
	return 0
    }

    return 0 unless $public || $owner eq 'sars2_bvbrc@patricbrc.org';
	    
    return 0 if $have_bad && $bad_genome{$gid};

    return 1;
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
