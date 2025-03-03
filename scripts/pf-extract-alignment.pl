#
# Extract an alignment from a local family btree
#

use strict;
use Getopt::Long::Descriptive;
use Bio::SeqIO;
use DB_File;

my($opt, $usage) = describe_options("%c %o btree-file family-id",
				    ["blast-ids" => "Convert ids to BLAST-compatible form"],
				    ["asn1" => "Write NCBI ASN.1 format"],
				    ["output|o=s" => "Write output to this file"],
				    ["help|h" => "Show this help message"]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) unless @ARGV == 2;

my $btree = shift;
my $fam_id = shift;

my %aligns;
tie %aligns, DB_File => $btree, 0, 0, $DB_BTREE or
    die "Cannot tie $btree: $!";

my $align = $aligns{$fam_id};
if ($opt->blast_ids)
{
    $align =~ s/^>fig\|/>fig-/mg;
}

my $out_fh;
if ($opt->output)
{
    open($out_fh, ">", $opt->output) or die "Cannot write " . $opt->output . ": $!";
}
else
{
    $out_fh = \*STDOUT;
}

if ($opt->asn1)
{
    my $seqio_in = Bio::SeqIO->new(-string => $align, -format => 'fasta');

    my $seqio_out = Bio::SeqIO->new(-fh => $out_fh, -format => 'ASN1');

    while (my $seq = $seqio_in->next_seq) {
	print STDERR "Processing sequence: ", $seq->id, "\n";

	# Write the sequence to the ASN.1 file
	$seqio_out->write_seq($seq);
    }

}
else
{
    print $out_fh $align;
}
    
