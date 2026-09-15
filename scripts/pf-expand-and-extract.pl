#
# Given a set of genome IDs and a genus, plus a families file,
# extract the families, expanding the NR sequences for that geenus.
#

use strict;
use Text::CSV_XS qw(csv);
use Data::Dumper;

@ARGV >= 4 or die "Usage: genus-dir families-file genus genome-id ...\n";

my $genus_dir = shift;
my $fams = shift;
my $genus = shift;
my @genomes = @ARGV;

my $pegsyn = "$genus_dir/$genus/nr/peg.synonyms";
open(PS, "<", $pegsyn) or die "Cannot open $pegsyn\n";

open(FAMS, "<", $fams) or die "annot open families $fams\n";

my %genomes = map { $_ => 1 } @genomes;

my %expand;
while (<PS>)
{
    chomp;
    my($md5, $list) = split(/\t/);
    my($ref, @pegs) = map { s/,\d+$//; $_ } split(/;/, $list);

    @pegs = grep { my($g) = /fig\|(\d+\.\d+)/; $genomes{$g} } ($ref, @pegs);
    next unless @pegs;
    # print "$ref => @pegs\n";
    $expand{$ref} = \@pegs;
}
close(PS);

while (<FAMS>)
{
    chomp;
    # PGF_00016046	1628	 831	 fig|568815.3.peg.3185		305		L-asparaginase (EC 3.5.1.1)			2481	Brucella-234

    my($pgf, undef, undef, $peg, $len, $func, $lfid, $fgenus) = split(/\t/);

    next unless $genus eq $fgenus;

    my $expand = $expand{$peg};

    for my $epeg (@$expand)
    {
	my($genus_name, $genus_id) = $genus =~ /^(.*)-(\d+)$/;
	my $plf = sprintf("PLF_${genus_id}_%08d", $lfid);
	print join("\t", $epeg, $pgf, $plf, $func), "\n";
    }
}
close(FAMS);
