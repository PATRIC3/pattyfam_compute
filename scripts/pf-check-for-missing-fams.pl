#
# Given a genus directory with an existing-fam-membership subdir that has files
# genome-id.plfam and genome-id.pgfam, determine the families that are found in
# genomes not in the active.genomes list.
#

use strict;
use File::Slurp;
use P3DataAPI;
use Data::Dumper;

@ARGV == 1 or die "Usage: $0 genus-dir\n";

my $api = P3DataAPI->new;

my $genus_dir = shift;

my @all_genomes = read_file("$genus_dir/genomes");
@all_genomes or die "Could not read $genus_dir/genomes";
chomp @all_genomes;

my @active_genomes = read_file("$genus_dir/active.genomes");
@active_genomes or die "Could not read $genus_dir/active.genomes";
chomp @active_genomes;
my %active_genomes = map { $_ => 1 } @active_genomes;

my %all_missing;
for my $fam (qw(plfam pgfam))
{
    my %active_fams;
    for my $g (@active_genomes)
    {
	my $fam_file = "$genus_dir/existing-fam-membership/$g.$fam";
	my @fams = read_file($fam_file, err_mode => 'quiet');
	if (!@fams)
	{
	    warn "Could not read active $fam_file\n";
	    next;
	}
	chomp @fams;
	$active_fams{$_} = 1 foreach @fams;
    }

    for my $g (grep { ! $active_genomes{$_} } @all_genomes)
    {
	my $fam_file = "$genus_dir/existing-fam-membership/$g.$fam";
	my @fams = read_file($fam_file, err_mode => 'quiet');
	if (!@fams)
	{
	    warn "Could not read $fam_file\n";
	    next;
	}
	chomp @fams;
	my @missing = grep { ! $active_fams{$_} } @fams;
	$all_missing{$_} = 1 foreach @missing;
	print STDERR "$g: missing " . scalar(@missing) . "\n";
    }

    my @work = keys %all_missing;
    my %def;
    while (@work)
    {
	my @chunk = splice(@work, 0, 500);
	my $q = join(",", @chunk);
	my @res = $api->query("protein_family_ref", ['in', 'family_id', "($q)"], ["select", "family_id,family_product"]);
	$def{$_->{family_id}} = $_->{family_product} foreach @res;
    }
    for my $fam (sort keys %all_missing)
    {
	print "$fam\t$def{$fam}\n";
    }
    exit;
}

