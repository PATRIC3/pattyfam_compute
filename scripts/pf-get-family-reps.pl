#
# Retrieve the reps from the given family and genus dir.
#

use strict;
use DB_File;

@ARGV == 2 or die "Usage: $0 fam-dir fam-id\n";

my $dir = shift;
my $fam = shift;

my %reps;

my $rep_file = "$dir/reps.aa.btree";

tie %reps, DB_File => $rep_file, O_RDONLY, 0, $DB_BTREE or die "Can't tie $rep_file: $!";

my $reps = $reps{$fam};

print $reps;
