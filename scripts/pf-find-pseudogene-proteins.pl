#
# Find proteins in the set of given genomes which are in range of a pseudogene.
# Mark sets that appear to be frameshifted clusters with the same annotation.
#

use strict;
use Data::Dumper;
use Getopt::Long::Descriptive;
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

my $json = JSON::XS->new->pretty->canonical;

my $taxon_id = read_file("$genus_dir/GENUS_TAXON_ID");
chomp $taxon_id;

my $pseudo_cache = "$genus_dir/pseudogenes";
make_path($pseudo_cache);

my $api = P3DataAPI->new();

#
# Look up pseudogenes in the genomes.
#

my $gid_query = "(" . join(",", @genomes) . ")";

my %pseudo;
my %pseudo_by_genome;
my @req;

my $pseudo_cb = sub {
    my($data) = @_;
    for my $ent (@$data)
    {
	my($id, $acc, $start, $end, $strand, $genome, $product) = @$ent{qw(patric_id accession start end strand genome_id product)};
	push(@{$pseudo{$genome, $acc}}, [$id, $acc, $start, $end, $strand, $genome, $product]);
	push(@{$pseudo_by_genome{$genome}}, [$id, $acc, $start, $end, $strand, $genome, $product]);
	push(@req, [$genome, $acc, $start, $end]);
    }
};

my $genomes_query = $api->query_cb('genome_feature',
				   $pseudo_cb,
				   ['in', 'genome_id', $gid_query],
				   ['eq', 'feature_type', 'pseudogene'],
				   ['select', 'start,end,strand,accession,genome_id,patric_id,product']
				  );

my @genomes_with_pseudo = grep { exists $pseudo_by_genome{$_} } @genomes;

printf "%d genomes, %d with pseudogenes\n", scalar @genomes, scalar @genomes_with_pseudo;
pareach(\@genomes_with_pseudo, sub {
    my($genome) = @_;
    my @req;
    for my $ent (@{$pseudo_by_genome{$genome}})
    {
	my($id, $acc, $start, $end, $strand, $genome, $product) = @$ent;
	push(@req, [$genome, $acc, $start, $end]);
    }
    
    my @gir = $api->genes_in_region_bulk(\@req);

    my %pseudo_by_function;
    for my $ent (@gir)
    {
	my($list, $beg, $end) = @$ent;
	push(@{$pseudo_by_function{$_->{product}}}, $_) foreach grep { $_->{feature_type} eq 'CDS'} @$list;
    }

    write_file("$pseudo_cache/$genome.json", $json->encode(\@gir));
    write_file("$pseudo_cache/$genome.by-function.json", $json->encode(\%pseudo_by_function));

}, { Max_Workers => $opt->parallel });

__END__


for my $set (@gir)
{
    my($features, $start, $end) = @$set;
    my @prots = grep { $_->{patric_id} =~ /\.peg\./ } @$features;
    print Dumper(\@prots);
}

