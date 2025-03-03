#
# Read the coreseed local fams, and write an overrides file based on all local fams there
# larger than the given size.
#

use strict;
use Data::Dumper;
use Getopt::Long::Descriptive;

my($opt, $usage) = describe_options("%c %o local-fam-dir output-file",
				    ["min-fam-size=i" => "Minimum size of family to include", { default => 5 }],
				    ["help|h" => "Show this help message"]);
print($usage->text), exit 1 if $opt->help;
die($usage->text) if @ARGV != 2;

my $fam_dir = shift;
my $output_file = shift;

open(OUT, ">", $output_file) or die "Cannot write $output_file: $!";

my %fams;

my $next_id = 1;
open(DEF, "<", "$fam_dir/local.family.defs") or die "Cannot read $fam_dir/local.family.defs: $!";
while (<DEF>)
{
    chomp;
    my @ent = split(/\t/);
    my ($fam_id, $fam_fn, $cluster, $mean, $dev, $gene, $fam_size) = @ent;
    if ($fam_size < $opt->min_fam_size)
    {
	last;
    }

    next if $fam_fn ne 'hypothetical protein';

    my $tag = substr($cluster, 0, 1);
    my $new_name = sprintf "hypothetical protein core_cid:$tag-%06d", $next_id++;
	
    my $fament =  { fam_id => $fam_id, fam_fn => $fam_fn, cluster => $cluster, mean => $mean, dev => $dev, gene => $gene, fam_size => $fam_size, members => [], name => $new_name };
    $fams{$fam_id} = $fament;
}
close(DEF);

open(FAMS, "<", "$fam_dir/local.family.members") or die "Cannot read $fam_dir/local.family.members: $!";
while (<FAMS>)
{
    chomp;
    my($fam_id, $peg, $len, $zscore) = split(/\t/);
    my $fam = $fams{$fam_id};
    next unless $fam;
    
    print OUT join("\t", $peg, $fam->{name}), "\n";
}
close(FAMS);
close(OUT);
