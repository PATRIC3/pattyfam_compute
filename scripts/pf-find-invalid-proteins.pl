#
# Find proteins in the given file that do not have valid start codons.
#

use strict;
use Data::Dumper;
use Getopt::Long::Descriptive;
use NCBI_genetic_code;
use gjoseqlib;
use List::MoreUtils 'first_index';
use JSON::XS;
use P3DataAPI;
use File::Slurp;
use Proc::ParallelLoop;
use File::Path qw(make_path);
use IPC::Run qw(run start finish);
use IO::Handle;

my($opt, $usage) = describe_options("%c %o fam-dir genus-dir genome-id [genome-id...]",
				    ["parallel|j=i" => "Number of processes", { default => 1 }],
				    ["genome-dir=s", "Directory holding PATRIC genome data", { default => "/vol/patric3/downloads/genomes" }],
				    ["help|h" => "Show this help message."]);

print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV < 3;

my $fam_dir = shift;
my $genus_dir = shift;
my @genomes = @ARGV;
my $genetic_code;

my $taxon_id = read_file("$genus_dir/GENUS_TAXON_ID");
chomp $taxon_id;

if (open(F, "$genus_dir/../genus_data.json"))
{
    my $genus_data = decode_json(scalar read_file(\*F));
    close(F);
    my $ginfo = $genus_data->{$taxon_id};
    $genetic_code = 11;
    if ($ginfo)
    {
	$genetic_code = $ginfo->{genetic_code};
    }
    else
    {
	die "Cannot find code for $genus_dir\n";
    }
}

my $contig_cache = "$genus_dir/contigs";
my $call_cache = "$genus_dir/prodigal_calls";
my $pseudo_cache = "$genus_dir/pseudogenes";
make_path($contig_cache, $call_cache, $pseudo_cache);

my $api = P3DataAPI->new();

#
# Look up pseudogenes in the genomes.
#

my $gid_query = "(" . join(",", @genomes) . ")";

my %pseudo;
my %pseudo_by_genome;

my $pseudo_cb = sub {
    my($data) = @_;
    for my $ent (@$data)
    {
	my($id, $acc, $start, $end, $strand, $genome, $product) = @$ent{qw(patric_id accession start end strand genome_id product)};
	push(@{$pseudo{$genome, $acc}}, [$id, $acc, $start, $end, $strand, $genome, $product]);
	push(@{$pseudo_by_genome{$genome}}, [$id, $acc, $start, $end, $strand, $genome, $product]);
    }
};

my $genomes_query = $api->query_cb('genome_feature',
				   $pseudo_cb,
				   ['in', 'genome_id', $gid_query],
				   ['eq', 'feature_type', 'pseudogene'],
				   ['select', 'start,end,strand,accession,genome_id,patric_id,product']
				  );

print Dumper(\%pseudo);

if ($opt->parallel > 1)
{
    pareach(\@genomes, sub {
	my($genome) = @_;
	process_genome($genome);
    }, { Max_Workers => $opt->parallel });
}
else
{
    for my $genome (@genomes)
    {
	process_genome($genome);
    }
}

#
# Process a genome:
#  Find contigs
#  Find features
#  Run prodigal
#  Match features against prodigal, looking for mismatches and prodigal-called features that are
#  marked partial.
# 

sub process_genome
{
    my($gid) = @_;

    my $contigs = $opt->genome_dir . "/$gid/$gid.fna";
    if (! -s $contigs)
    {
	$contigs = "$contig_cache/$gid.fna";
	if (! -s $contigs)
	{
	    my $ok = run(["p3-genome-fasta", "--contig", $gid], ">", $contigs);
	    if (! $ok)
	    {
		die "Cannot load contigs for $gid\n";
	    }
	}
    }
    # print "Found contigs $contigs\n";

    #
    # Find call data.
    #
    my @features;
    my %features;
    my %features2;		# features by accession
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
	my $strand_idx = first_index { $_ eq 'strand' } @hdrs;
	my $product_idx = first_index { $_ eq 'product' } @hdrs;
	
	while (<TAB>)
	{
	    chomp;
	    my @vals = split(/\t/);
	    my $id = $vals[$id_idx];
	    my $start = $vals[$start_idx];
	    my $end = $vals[$end_idx];
	    my $acc = $vals[$acc_idx];
	    my $product = $vals[$product_idx];

	    next if $vals[$ftype_idx] ne 'CDS';
	    my $strand = $vals[$strand_idx];
	    push(@features, [$id, $acc, $start, $end, $product]);
	    push(@{$features{$acc, $strand, $end}}, [$id, $start,$product]);
	    push(@{$features2{$acc}}, [$id, $strand, $start, $end,$product]);
	}
	close(TAB);
    }
    else
    {
	my $inp = "it\n$gid\n";

	my $out;
	my $h = run(["p3-get-genome-features",
		     "--eq", "feature_type,CDS",
		     "--attr", "accession,patric_id,start,end,strand,product"],
		    "<", \$inp,
		    ">", \$out);
	open(S, "<", \$out);
	$_ = <S>;
	while (<S>)
	{
	    chomp;
	    my($gid, $acc, $id, $start, $end, $strand, $product) = split(/\t/);
	    push(@features, [$id, $acc, $start, $end,$product]);
	    push(@{$features{$acc, $strand, $end}}, [$id, $start,$product]);
	    push(@{$features2{$acc}}, [$id, $strand, $start, $end,$product]);
	}
    }

    #
    # Scan for features overlapped by pseudogene
    #

    while (my($acc, $flist) = each(%features2))
    {
	my $plist = $pseudo{$gid, $acc};
	
	for my $pg (@$plist)
	{
	    my($id, undef, $start, $end, $strand, undef, $product) = @$pg;

	    my @hits;
	    for my $f (@$flist)
	    {
		my($fid, $fstrand, $fstart, $fend) = @$f;

		# If feature is entirely contained in pseudo
		if ($fstart >= $start && $fend <= $end)
		{
		    push(@hits, $f);
		    # print join("\t", $fid, $fstart, $fend, $id, $start, $end), "\n";
			
		}
	    }
	    if (@hits)
	    {
		print Dumper($pg, \@hits);
	    }
	}
    }
    next;

    my $prodigal_out = "$call_cache/$gid.gff";
    if (! -s $prodigal_out)
    {
	my $ok = run(["prodigal", "-g", $genetic_code, "-f", "gff", "-o", $prodigal_out, "-i", $contigs]);
	$ok or die "Error running prodigal\n";
    }

    if (!open(PROD, "<", $prodigal_out))
    {
	die "Cannot open prodigal output $prodigal_out: $!";
    }

    while (<PROD>)
    {
	next if /^\#/;
	chomp;
	my($acc, $what, $type, $start, $end, undef, $strand, undef, $attrs) = split(/\t/);
	my %attrs = map { split(/=/, $_) } split(/;/, $attrs) ;

	if ($strand eq '-')
	{
	    ($end, $start) = ($start, $end);
	}
	my $match = delete $features{$acc, $strand, $end};

	next unless $type eq 'CDS';
	next if $attrs{partial} eq '00';
	if ($match)
	{
	    for my $m (@$match)
	    {
		my($id, $start) = @$m;
		print join("\t", $id, @attrs{qw(partial score rscore tscore uscore start_type)}), "\n";
	    }
	}
    }

    for my $k (sort keys %features)
    {
	my($acc, $strand, $end) = split(/$;/, $k);

	
	my $mlist = $features{$k};
	for my $m (@$mlist)
	{
	    my($id, $start) = @$m;
	    print STDERR join("\t", $id, $start, $end, $strand), "\n";
	}
    }
    
    close(PROD);
}
