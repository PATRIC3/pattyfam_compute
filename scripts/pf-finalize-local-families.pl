#
# This is a SAS Component
#

=head1 NAME

    pf-finalize-local-families

=head1 SYNOPSIS

    pf-finalize-local-families

=head1 DESCRIPTION

This is a pattyfam computation support routine.

After the core computational stages of the local families build are complete, combine the 
kmer clusters, blast clusters, and singleton families into a combined families file
and compute statistics.

=cut
    
=head2 Output Format

Output is written to STDOUT and constitutes the derived protein families (which
include singletons).  An 8-column, tab-separated table is written:

    FamilyID - an integer
    Function - function assigned to family
    SubFunction - the Function and an integer (SubFunction) together uniquely
                  determine the FamilyID.  Another way to look at it is

                    a) each family is assigned a unique ID and a function
                    b) multiple families can have the same function (consider
                       "hypothetical protein")
                    c) the Function+SubFunction uniquely determine the FamilyID
    PEG
    LengthProt - the length of the translated PEG
    Mean       - the mean length of PEGs in the family
    StdDev     - standard deviation of lengths for family
    Z-sc       - the Z-score associated with the length of this PEG

=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use SeedEnv;
use gjoseqlib;
use gjostat;
use DB_File;

use Getopt::Long::Descriptive;

my($opt, $usage) = describe_options("%c %o fam-dir fam-file [fam-file...]",
				    ["parallel=i" => "Parallel threads", { default => 1 }],
				    ["help|h" => "Show this help message."]);

print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV < 2;

my $fam_dir = shift;
my @fam_files = @ARGV;

#
# we write to fam-dir/local.family.members and local.family.defs
#

open(FMEM, ">", "$fam_dir/local.family.members") or die "cannot write $fam_dir/local.family.members: $!";
open(FDEF, ">", "$fam_dir/local.family.defs") or die "cannot write $fam_dir/local.family.members: $!";
    
#
# Read family files and break into family sets.
#

my %seq_len;
tie %seq_len, DB_File => "$fam_dir/seq_len.db", O_RDONLY, 0, $DB_BTREE or die "cannot tie $fam_dir/seq_len.db: $!";

my @sets;
my $last_key;
my $cur_set;
my $src_key = 'A';
for my $fam (@fam_files)
{
    open(F, "<", $fam) or die "Cannot open $fam: $!";
    while (<F>)
    {
	chomp;
	my($func, $fam, $peg, $gene_name, $consensus_gene_name) = split(/\t/);
	$fam = "$src_key$fam";
	my $key = "$func$;$fam";
	if ($key ne $last_key)
	{
	    $cur_set = [];
	    push(@sets, $cur_set);
	    $last_key = $key;
	}
	push(@$cur_set, [$func, $fam, $peg, $gene_name, $consensus_gene_name, $seq_len{$peg}]);
    }
    close(F);
    $src_key++;
}
untie %seq_len;

my $famN = 1;
foreach my $set (sort { (@$b <=> @$a) or 
			($a->[0] cmp $b->[0]) or
		        ($a->[1] <=> $b->[1])
		      } @sets)
{
    my($mean,$std_dev) = &mean_stddev(map { $_->[5] } @$set);

    my($fam, $subfam, undef, undef, $consensus_gene_name) = @{$set->[0]};
    
    print FDEF join("\t", $famN, $fam, $subfam, 
		    sprintf("%0.3f",$mean), sprintf("%0.3f",$std_dev),
		    $consensus_gene_name,
		    scalar @$set), "\n";

    foreach my $tuple (@$set)
    {
	my($fam, $subfam, $peg, $gene_name, $consensus_gene_name, $ln) = @$tuple;

	my $z_sc = ($ln - $mean) / ($std_dev + 0.000000001);
	print FMEM join("\t",
			$famN, $peg, $ln,
			sprintf("%0.3f",$z_sc),
			$gene_name), "\n";
    }
    $famN++;
}
close(FDEF);
close(FMEM);

sub mean_stddev
{
    my $n = @_;
    return $n ? ( shift, undef) : ( undef, undef ) if $n < 2;
    my ( $sum, $sum2 ) = ( 0, 0 );
    foreach ( @_ ) { $sum += $_; $sum2 += $_ * $_ }
    my $x = $sum / $n;
    my $to_check = ( $sum2 - ($n * $x * $x )) / ( $n - 1 );
    ( $x, ($to_check > 0) ? sqrt($to_check) : 0);
}
