#
# Given a genome ID and the location of final families, create a comparison
# between the genome in /vol/patric3/downloads/ and the new annotations.
#
# /vol/bvbrc-families/fams/big.2024-1221/merge.full/prop.4/final
#

use strict;
use Text::CSV_XS qw(csv);
use Data::Dumper;

@ARGV == 2 or @ARGV == 3 or die "Usage: genome-id final-fams [propagated-file]\n";

my $gid = shift;
my $famdir = shift;
my $prop_file = shift;

my $bvbrc = "/vol/patric3/downloads/genomes/$gid/$gid.PATRIC.features.tab";

#
# the anno readers return a hash from fid => anno object
#

my $brc_annos = read_brc_annos($gid);
my $denovo_annos = read_denovo_annos("$famdir/by.genome/$gid");

my $prop_annos;
if ($prop_file)
{
    $prop_annos = read_propagated_annos($prop_file);
}

#
# Loop over the brc annos since those have all pegs.
# For each, find the corresponding other anno, compare, and save.
#

my $null_anno = new Anno();

my @table;
for my $brc_anno (sort { $a->order <=> $b->order } values %$brc_annos)
{
    my $fid = $brc_anno->fid;
    my $denovo = $denovo_annos->{$fid} // $null_anno;;
    my $prop = $prop_annos->{$fid} // $null_anno;;
    my $code_bd = $brc_anno->compare($denovo);
    my $code_bp = $brc_anno->compare($prop);
    my $code_dp = $denovo->compare($prop);
    push(@table, [$fid, $code_bd, $code_bp, $code_dp, $brc_anno->order, $brc_anno, $denovo, $prop]);
}

@table = sort { $b->[1] cmp $a->[1] or  $a->[2] <=> $b->[2] } @table;

print join("\t", "ID", "BRC vs Denovo", "BRC vs Prop", "Denovo vs Prop", "Order",
	   "BRC Func", "Denovo Func", "Prop Func",
	   "BRC PGF", "Denovo PGF", "Prop PGF",
	   "BRC PLF", "Denovo PLF", "Prop PLF"), "\n";

for my $ent (@table)
{
    my($fid, $code_bd, $code_bp, $code_dp, $order, $brc, $denovo, $prop) = @$ent;
    print join("\t", $fid, $code_bd, $code_bp, $code_dp, $order,
	       $brc->func, $denovo->func, $prop->func,
	       $brc->pgf, $denovo->pgf, $prop->pgf,
	       $brc->plf, $denovo->plf, $prop->plf), "\n";
}

sub read_denovo_annos
{
    my($famfile) = @_;
    
    open(F, "<", $famfile) or die "cannot open $famfile: $!";
    
    my %dat;
    while (<F>)
    {
	chomp;
	
	# fig|568815.3.peg.3098	PGF_00001235	PLF_234_00001935	DUF1178 protein RSP_2387
	
	
	my($peg, $pgf, $plf, $func) = split(/\t/);
	
	$dat{$peg} = new Anno($peg, $func, $pgf, $plf);
    }
    close(F);
    return \%dat;
}

sub read_brc_annos
{
    my($gid) = @_;
    
# genome_id	genome_name	accession	annotation	feature_type	patric_id	refseq_locus_tag	start	end	strand	na_length	gene	productplfam_id	pgfam_id
# 568815.3	Brucella microti CCM 4915	NC_013118	PATRIC		source						1	1220319	+		1220319

    my $bfeats = csv(in => $bvbrc, headers => 'auto', sep_char => "\t");
    return { map { my $a = Anno->new_from_brc($_); ( $a->fid => $a ) } grep { $_->{feature_type} eq 'CDS' } @$bfeats };
}    
      
sub read_propagated_annos
{
    my($prop_file) = @_;

    open(P, "<", $prop_file) or die "Cannot read $prop_file: $!";
    my %ret;
    while (<P>)
    {
	chomp;
	my($fid, $pgf, undef, $plf, undef, $func) = split(/\t/);
	($fid) = $fid =~ /^(fig\|\d+\.\d+\.peg\.\d+)/;
	next unless $fid;
	$ret{$fid} = new Anno($fid, $func, $pgf, $plf);
    }

    return \%ret;
}    
      
package Anno;

use strict;
use Data::Dumper;

BEGIN {
use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(pgf plf func fid order));
}

sub new_from_brc
{
    my($class, $brc_feature) = @_;
    return $class->new(@$brc_feature{qw(patric_id product pgfam_id plfam_id)});
}

sub new
{
    my($class, $fid, $func, $pgf, $plf) = @_;
    my $self = {
	fid => $fid,
	func => $func,
	pgf => $pgf,
	plf => $plf,
    };
    ($self->{order}) = $fid =~ /\.peg\.(\d+)/;
    bless $self, $class;
    return $self;
}

sub compare
{
    my($self, $other) = @_;

    my $func_code = ($self->func ne $other->func) ? 1 : 0;
    my $pgf_code = ($self->pgf ne $other->pgf) ? 1 : 0;
    my $plf_code = ($self->plf ne $other->plf) ? 1 : 0;

    my $code = join('', "C", $func_code, $pgf_code, $plf_code);

    return wantarray ?  ($code, $func_code, $plf_code, $plf_code) : $code;
}
