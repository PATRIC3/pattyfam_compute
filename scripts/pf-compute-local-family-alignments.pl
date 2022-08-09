

=head1 NAME

    pf-compute-local-family-alignments

=head1 SYNOPSIS

    pf-compute-local-family-alignments 

=head1 DESCRIPTION

This is a pattyfam computation support routine.

Given a genus data directory,  compute the set of multiple sequence alignments for the families there.

=cut

use P3DataAPI;
use Proc::ParallelLoop;
use File::Copy 'copy';
use IPC::Run 'run';
use IO::File;
use File::Basename;
use strict;
use Data::Dumper;
use Getopt::Long::Descriptive;
use gjoseqlib;
use Digest::MD5 'md5_hex';
use LPTScheduler;
use P3AlignmentComputeToDisk;
use DB_File;
use File::Slurp;

my($opt, $usage) = describe_options("%c %o genus-dir work-dir dest-dir",
				    ["dna-md5-map=s" => "Save DNA md5 mapping here"],
				    ["max-sim=s" => "Maximium identity to retained sequences", { default => 0.8 }],
				    ["genome-dir=s", "Directory holding PATRIC genome data", { default => "/vol/patric3/downloads/genomes" }],
				    ["parallel|p=i" => "run this many parallel alignments", { default => 1 }],
				    ["help|h" => "show this help message"]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) unless @ARGV == 3;

my $api = P3DataAPI->new();

my $dir = shift;
my $work_dir = shift;
my $dest_dir = shift;

-d $dest_dir or die "Destination dir $dest_dir does not exist\n";

my $genus = basename($dir);

my $genus_tax_id;
if ($genus =~ /^(.*)-(\d+)$/)
{
    $genus = $1;
    $genus_tax_id = $2;
}
    

my $bootstrap = sub {
    my $p = P3AlignmentComputeToDisk->new($work_dir);
    return $p;
};

open(FAM, "<", "$dir/local.family.members") or die "$dir/local.family.members does not exist\n";

#
# %fidfam maps from genome id to a hash mapping from feature id to the containing family
#
my %fidfam;
my $last;
my %fam_aa_data;
my %fam_na_data;
my %fam_size;
my %fam_aa_size;
my %fam_na_md5;

while (<FAM>)
{
    chomp;
    my($lid, $fid) = split(/\t/);
    my($g) = $fid =~ /^fig\|(\d+\.\d+)\./;
    $fidfam{$g}->{$fid} = $lid;
    $fam_size{$lid}++;
    if ($lid ne $last)
    {
	if (defined($last) && $fam_size{$last} > 1)
	{
	    my $fa = '';
	    my $na = '';
	    $fam_aa_data{$last} = \$fa;
	    $fam_na_data{$last} = \$na;
	}
	$last = $lid;
    }
}
if (defined($last) && $fam_size{$last} > 1)
{
    my $fa = '';
    my $na = '';
    $fam_aa_data{$last} = \$fa;
    $fam_na_data{$last} = \$na;
}

my $map;
if ($opt->dna_md5_map)
{
    $map = IO::File->new($opt->dna_md5_map, ">");
}

my @genomes = sort keys %fidfam;

my $genetic_code = 11;
my @qry;
if ($genus_tax_id)
{
    @qry = (["eq", "taxon_id", $genus_tax_id]);

}
else
{
    @qry = (["eq", "taxon_rank", "genus"], ["eq", "taxon_name", $genus]);
}
my @res = $api->query("taxonomy", @qry, ["select", "taxon_name,taxon_id,genetic_code"]);
print STDERR  "Genus tax data for '$genus': " . Dumper(\@res);
if (@res)
{
    $genetic_code = $res[0]->{genetic_code};
}

#
# Read fasta data into memory on a per-family basis
#
for my $g (@genomes)
{
    my $fidh = $fidfam{$g};
    if (open(S, "<", "$dir/nr-seqs/$g"))
    {
	my $inref;
	while (<S>)
	{
	    if (/^>(\S+)/)
	    {
		if (defined(my $fam = $fidh->{$1}))
		{
		    $fam_aa_size{$fam}++;
		    $inref = $fam_aa_data{$fam};
		    $$inref .= $_ if $inref;
		}
		else
		{
		    undef $inref;
		}
	    }
	    elsif ($inref)
	    {
		$$inref .= $_;
	    }
	}
	close(S);
    }
    #
    # First try to load from nr-seq-dna data.
    #
    my $nr_seq_dna = "$dir/nr-seqs-dna/$g";
    my $gfna = $opt->genome_dir . "/$g/$g.PATRIC.ffn";
    my $opened;

    if (open(S, "<", $nr_seq_dna))
    {
	$opened++;
    }
    else
    {
	warn "$nr_seq_dna missing, trying $gfna\n";
	if (open(S, "<", $gfna))
	{
	    $opened++;
	}
	else
	{
	    warn "No DNA found\n";
	}
    }
    if ($opened)
    {
	my $inref;
	my $digest;
	my $dna;
	my $id;
	my $fam;
	while (<S>)
	{
	    if (/^>(fig\|\d+\.\d+\.[^.]+\.\d+)/)
	    {
		my $new_id = $1;
		
		#
		# process last fam if we were in one.
		#
		if ($inref)
		{
		    my $md5 = $digest->hexdigest();
		    print $map "$id\t$md5\n" if $map;
		    if (!$fam_na_md5{$fam}->{$md5})
		    {
			$fam_na_md5{$fam}->{$md5} = 1;
			$$inref .= ">$id\n$dna";
			undef $dna;
		    }
		}

		$id = $new_id;
		$fam = $fidh->{$id};

		#
		# Now start new seq processing.
		#
		$inref = $fam_na_data{$fam};		
		if ($inref)
		{
		    $digest = Digest::MD5->new();
		    $dna = '';
		}
	    }
	    else
	    {
		if ($inref)
		{
		    s/\s//g;
		    next if $_ eq '';
		    $digest->add(uc($_));
		    $dna .= $_ . "\n";
		}
	    }
	}

	if ($inref)
	{
	    my $md5 = $digest->hexdigest();
	    print $map "$id\t$md5\n" if $map;
	    if (!$fam_na_md5{$fam}->{$md5})
	    {
		$fam_na_md5{$fam}->{$md5} = 1;
		$$inref .= ">$id\n$dna\n";
		undef $dna;
	    }
	}
	
	close(S);
    }
}
close($map) if $map;

for my $k (keys %fam_aa_size)
{
    if ($fam_aa_size{$k} <= 1)
    {
	# print "Remove fam $k size=$fam_aa_size{$k}\n";
	delete $fam_aa_data{$k};
    }
}
while (my($k, $h) = each %fam_na_md5)
{
    if (keys %$h <= 1)
    {
	# print "Remove fam $k from dna\n";
	delete $fam_na_data{$k};
    }
}

my $sched = LPTScheduler->new($opt->parallel * 10);

my %fams_to_compute;


while (my($fam, $seqp) = each %fam_aa_data)
{
    $fams_to_compute{$fam} = 1;
    my $cost = length($$seqp);
    $cost = 5000 if $cost < 5000;
    $sched->add_work([$fam, $seqp, 'aa', $genus, $genetic_code], $cost);
}

while (my($fam, $seqp) = each %fam_na_data)
{
    $fams_to_compute{$fam} = 1;
    my $cost = length($$seqp);
    $cost = 5000 if $cost < 5000;
    $sched->add_work([$fam, $seqp, 'dna', $genus, $genetic_code], $cost);
}

$sched->run($bootstrap, \&align_one, $opt->parallel);

#
# When alignments are complete, roll the stats, reps and alignments up into btree databases.
#


for my $type ('aa', 'dna')
{
    my(%stats, %reps, %alignments_clean, %alignments_raw);

    tie %stats, DB_File => "$work_dir/stats.$type.btree", O_RDWR | O_CREAT, 0644, $DB_BTREE
	or die "Cannot tie $work_dir/stats.$type.btree: $!";
    
    tie %reps, DB_File => "$work_dir/reps.$type.btree", O_RDWR | O_CREAT, 0644, $DB_BTREE
	or die "Cannot tie $work_dir/reps.$type.btree: $!";

    tie %alignments_clean, DB_File => "$work_dir/alignments.clean.$type.btree", O_RDWR | O_CREAT, 0644, $DB_BTREE
	or die "Cannot tie $work_dir/alignments.clean.$type.btree: $!";
    tie %alignments_raw, DB_File => "$work_dir/alignments.raw.$type.btree", O_RDWR | O_CREAT, 0644, $DB_BTREE
	or die "Cannot tie $work_dir/alignments.raw.$type.btree: $!";

    for my $fam (keys %fams_to_compute)
    {
	load(\%stats, $fam, "stats.$type.json");
	load(\%reps, $fam, "reps.$type");
	load(\%alignments_clean, $fam, "clean_$type.fa");
	load(\%alignments_raw, $fam, "raw_$type.fa");
    }

    #
    # Untie and copy into final place.
    #
    untie %stats;
    untie %reps;
    untie %alignments_clean;
    untie %alignments_raw;
    for my $file ("stats.$type.btree", "reps.$type.btree",
		  "alignments.clean.$type.btree", "alignments.raw.$type.btree")
    {
	print STDERR "Copy $work_dir/$file to $dest_dir/$file\n";
	copy("$work_dir/$file", "$dest_dir/$file");
    }
}

#
# Load data file into btree.
#
sub load
{
    my($hash, $fam, $file) = @_;
    if (open(X, "<", "$work_dir/$fam/$file"))
    {
	$hash->{$fam} = read_file(\*X);
	close(X);
    }
    else
    {
	warn "Cannot open $work_dir/$fam/$file: $!";
    }
}

sub align_one
{
    my($p3align, $work) = @_;
    my($key, $seqp, $tag, $genus, $genetic_code) = @$work;

    open(STR, "<", $seqp);
    my @data = read_fasta(\*STR);
    close(STR);

    $p3align->process_data($genus, $key, \@data, $tag, $genetic_code);
}

