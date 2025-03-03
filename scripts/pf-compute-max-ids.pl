#
# Given a families file, compute the max ids seen for global and local families.
#

use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use IO::Handle;

my($opt, $usage) = describe_options("%c %o fam-file max-id-file",
				    ["help|h" => "Show this help message"],
				    );
print($usage->text), exit 0 if $opt->help;
die($usage->text), if @ARGV != 2;

my $fam_file = shift;
my $max_id_file = shift;

my %genus_to_fh;
my %genus_max;

open(F, "<", $fam_file) or die "Cannot open families file $fam_file: $!\n";
open(MAX, ">", $max_id_file) or die "Cannot write $max_id_file: $!\n";

while (<F>)
{
    my($gid, $genus, $num) = /^(?:PGF_|GF)(\d+).*?\t([^\t]+)\t(\d+)$/;
    $genus_max{GLOBAL} = $gid if $gid > $genus_max{GLOBAL};
    $genus_max{$genus} = $num if $num > $genus_max{$genus};
}

for my $genus (sort keys %genus_max)
{
    print MAX "$genus\t$genus_max{$genus}\n";
}
close(MAX);

