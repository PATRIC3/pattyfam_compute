=head1 NAME

    pf-construct-genus-data-for-coreseed

=head1 SYNOPSIS

    pf-construct-genus-data-for-coreseed kmer-dir

=head1 DESCRIPTION


Build a special-purpose data directory for the purpose of calling tight families on coreseed genomes.

=cut

use strict;
use File::Path 'make_path';
use File::Slurp;
use LWP::UserAgent;
use Getopt::Long::Descriptive;
use List::MoreUtils 'first_index';
use Data::Dumper;
use URI;
use JSON::XS;
use P3DataAPI;

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

my($opt, $usage) = describe_options("%c %o kmer-dir data-dir",
				    ["parallel|p=i", "Run with this many procs", { default => 1 }],
				    ["help|h", "Show this help message"],
				    );

print($usage->text), exit if $opt->help;
die($usage->text) if @ARGV != 2;

my $kmer_dir = shift;
my $base_data_dir = shift;

-d $base_data_dir or die "Data directory $base_data_dir does not exist\n";

my $json = JSON::XS->new->pretty(1);

#
# Fill in required data from kmer build dir.
#
my $genus = "core";
opendir(D, "$kmer_dir/Seqs") or die "cannot open $kmer_dir/Seqs: $!";
my $genome_ids = [sort grep { -f "$kmer_dir/Seqs/$_" && /^\d+\.\d+$/ } readdir(D)];
closedir(D);
@$genome_ids or die "Did not properly read genome ids";

my $genus_taxon;

my $api = P3DataAPI->new;
my $n1 = $api->genome_name($genome_ids);

my $gnames = {};
for my $gid (@$genome_ids)
{
    if ($n1->{$gid})
    {
	$gnames->{$gid} = $n1->{$gid}->[0];
    }
    elsif (-f "/vol/core-seed/FIGdisk/FIG/Data/Organisms/$gid/GENOME")
    {
	$gnames->{$gid} = read_file("/vol/core-seed/FIGdisk/FIG/Data/Organisms/$gid/GENOME");
    }
    else
    {
	$gnames->{$gid} = read_file("$kmer_dir/gnames/$gid");
    }
    chomp $gnames->{$gid};
}

#
# This code originated form the genus code.
#
if (1)
{
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
	my $prots = "$kmer_dir/Seqs/$gid";
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
	}
	else
	{
	    die "Could not read $prots: $!";
	}
	$tied->sync();
	print GFILE "$gid\n";
	print GNAME "$gid\t$gnames->{$gid}\n";
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
    my $rc = system("pf-build-nr", "$data_dir/sources", "$data_dir/nr/nr",
		    "$data_dir/nr/peg.synonyms", "$data_dir/nr/nr-len.btree", "$data_dir/nr/figids");
    if ($rc == 0)
    {
	$rc = system("pf-compute-nr-seqs", "--skip-dna",  $data_dir);
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
}

