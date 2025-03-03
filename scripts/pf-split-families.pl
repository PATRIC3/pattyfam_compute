#
# Given a families file, split into per-genus family files.
#

use strict;
use Getopt::Long::Descriptive;
use IO::Handle;

my($opt, $usage) = describe_options("%c %o fam-file output-dir max-id-file",
				    ["map-genera-from=s" => "Map bare genus names to genus-genusid based on the given directory of genus-genusid names"],
				    ["help|h" => "Show this help message"],
				    );
print($usage->text), exit 0 if $opt->help;
die($usage->text), if @ARGV != 3;

my $fam_file = shift;
my $output_dir = shift;
my $max_id_file = shift;

my %genus_to_fh;
my %genus_max;

-d $output_dir or die "Output directory $output_dir does not exist\n";

my %genus_map;
if ($opt->map_genera_from)
{
    opendir(D, $opt->map_genera_from) or die "Cannot opendir " . $opt->map_genera_from . ": $!";
    while (my $p = readdir(D))
    {
	if ($p =~ /^(\S+)-(\d+)$/)
	{
	    $genus_map{$1} = $p;
	}
    }
    closedir(D);
}
open(F, "<", $fam_file) or die "Cannot open families file $fam_file: $!\n";
open(MAX, ">", $max_id_file) or die "Cannot write $max_id_file: $!\n";

while (<F>)
{
    chomp;
    my($gid, $genus, $num) = /^(?:PGF_|GF)(\d+).*?\t([^\t]+)\t(\d+)$/;
    my @cols = split(/\t/);
    defined($gid) or die "cannot parse $_";

    if (my $m = $genus_map{$genus})
    {
	$genus = $m;
	$cols[7] = $m;
    }
    my $fh = $genus_to_fh{$genus};
    $genus_max{GLOBAL} = $gid if $gid > $genus_max{GLOBAL};
    $genus_max{$genus} = $num if $num > $genus_max{$genus};
    if (!$fh)
    {
	$fh = new IO::Handle;
	open($fh, ">", "$output_dir/$genus") or die "Cannot write $output_dir/$genus: $!";
	$genus_to_fh{$genus} = $fh;
    }
    print $fh join("\t", @cols), "\n";
}

for my $genus (sort keys %genus_to_fh)
{
    my $fh = $genus_to_fh{$genus};
    close($fh) or warn "Error closing genus file for $genus: $!\n";
}

for my $genus (sort keys %genus_max)
{
    print MAX "$genus\t$genus_max{$genus}\n";
}
close(MAX);

