
=head1 NAME

    pf-compute-nr-seqs

=head1 SYNOPSIS

    pf-compute-nr-seqs genus-data-dir

=head1 DESCRIPTION

This is a pattyfam computation support routine.

Given a genus data directory, compute a directory of sequence data files named
with the genome ID. These are derived from the NR and are used to accelerate
the family computation. As such we also write a tab-delimited file nr-refs containing
the pegs in the NR; the first peg on each row is the reference and the rest are
the duplicate pegs. This is derived from the NR peg.synonyms file and is used
later on to inflate the generated families file to include all pegs, including
the duplicates.

=cut

use strict;
use gjoseqlib;
use Getopt::Long::Descriptive;
use Data::Dumper;
use IO::File;

my($opt, $usage) = describe_options("%c %o genus-data-dir",
				    ["skip-dna", "Skip computing DNA NR"],
				    ["genome-dir=s", "Directory holding PATRIC genome data", { default => "/vol/patric3/downloads/genomes" }],
				    ["help|h" => "Show this help message."]);
print($usage->text), exit 1 if $opt->help;
die($usage->text) if @ARGV != 1;

my $fam_dir = shift;

my $refs_file = "$fam_dir/nr-refs";

if (! -f "$fam_dir/nr/peg.synonyms")
{
    die "Missing $fam_dir/nr/peg.synonyms\n";
}

if (! -s "$fam_dir/nr/peg.synonyms")
{
    die "Zero length $fam_dir/nr/peg.synonyms\n";
}

open(PEGSYN, "<", "$fam_dir/nr/peg.synonyms") or die "Cannot open $fam_dir/nr/peg.synonyms: $!";
open(NR, "<", "$fam_dir/nr/nr") or die "Cannot open $fam_dir/nr/nr: $!";

open(REFS, ">", $refs_file) or die "Cannot open refs file $refs_file: $!";

#
# In the event files already exist in nr-seqs, remove them to eliminate duplicates.
#

my $sdir = "$fam_dir/nr-seqs";
-d $sdir or mkdir($sdir);
my $dna_dir = "$fam_dir/nr-seqs-dna";
-d $dna_dir or mkdir($dna_dir);

opendir(DH, $sdir) or die "Cannot opendir $sdir: $!";
while (my $f = readdir(DH))
{
    next if $f =~ /^\./;
    my $p = "$sdir/$f";
    if (-s $p)
    {
	print STDERR "Unlinking existing NR file $p\n";
	unlink($p) or die "Cannot unlink $p: $!";
    }
}
closedir(DH);

my %f;
my %g;

#
# %genomes is a hash mapping from the genome ID to hash mapping feature ID to 1.
# It is used to pull the DNA sequences from the genomes repository for later use in the computation.
#
my %genomes;

# gnl|md5|ffdb15eded2cbc4a89b6e9ece42b30e9,108  fig|1005475.3.peg.3229,108;fig|1182692.3.peg.3088,108;fig|1182722.3.peg.1758,108
while (<PEGSYN>)
{
    if (/^([^,]+),\d+\t(fig\|(\d+\.\d+)\.(?:CDS|peg)\.\d+),\d+(.*)/)
    {
	my $md5 = $1;
	my $ref = $2;
	my $genome = $3;
	my $rest = $4;

	$genomes{$genome}->{$ref} = 1;
	
	$f{$md5} = $ref;
	$g{$md5} = $genome;

	my @pegs;
	while ($rest =~ /(fig\|\d+\.\d+\.(?:CDS|peg)\.\d+)/mg)
	{
	    push(@pegs, $1);
	}
	print REFS join("\t", $ref, @pegs), "\n";
    }
}
close(PEGSYN);

my $open_genome;
my $open_fh;

while (<NR>)
{
    if (/^>(\S+)(.*)/)
    {
	my $fid = $f{$1};
	my $genome = $g{$1};
	
	if (!$genome)
	{
	    warn "Cannot map $1 to genome\n";
	    next;
	}
	
	if ($genome ne $open_genome)
	{
	    close($open_fh) if $open_fh;
	    open($open_fh, ">>", "$sdir/$genome") or die "Cannot append $sdir/$genome: $!";
	    $open_genome = $genome;
	}		    
	
	print $open_fh ">$f{$1}$2\n";
    }
    else
    {
	print $open_fh $_;
    }
}
close(NR);
close($open_fh);

# 
# Ensure we have a file for all organisms, even if empty.
# Empty would mean all sequences were covered in NR from another
# organism.
#
opendir(DH, "$fam_dir/Seqs") or die "Cannot opendir $fam_dir/Seqs: $!";
while (my $ent = readdir(DH))
{
    if ($ent =~ /^\d+\.\d+$/)
    {
	my $f = "$fam_dir/nr-seqs/$ent";
	if (! -f $f)
	{
	    print STDERR "Writing empty nr-seq file $f\n";
	    open(TF, ">", $f) or die "Cannot open $f: $!";
	    close(TF);
	}
    }
}
closedir(DH);

#
# Read and copy the NR dna sequences from the genome directory.
#

if (!$opt->skip_dna)
{
    while (my($genome, $ids) = each %genomes)
    {
	next if -s "$dna_dir/$genome";
	
	my $dna1 = $opt->genome_dir . "/$genome/$genome.PATRIC.ffn";
	
	my $fh = IO::File->new($dna1, "r");
	if (!$fh)
	{
	    my $dna2 = "$fam_dir/NASeqs/$genome";
	    $fh = IO::File->new($dna2, "r");
	    if (!$fh)
	    {
		warn "could not read DNA from either $dna1 or $dna2\n";
	    }
	}
	if ($fh)
	{
	    open(O, ">", "$dna_dir/$genome") or die "cannot write $dna_dir/$genome: $!";
	    
	    while (my($idx, $def, $seq) = read_next_fasta_seq($fh))
	    {
		my($id) = $idx =~ (/(fig\|\d+\.\d+\.[^.]+\.\d+)/);
		
		if (delete $ids->{$id})
		{
		    write_fasta(\*O, [$id, $def, $seq]);
		}
	    }
	}
	
	close(O);
    }
}
