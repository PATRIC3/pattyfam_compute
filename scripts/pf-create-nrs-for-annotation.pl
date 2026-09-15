use strict;
use Data::Dumper;
use Getopt::Long::Descriptive;
use BerkeleyDB;
use File::Path qw(make_path);
use File::Copy qw(copy);
use Proc::ParallelLoop;

use gjoseqlib;

my($opt, $usage) = describe_options("%c %o genus-data-dir prop-out-dir fam-members-file",
				    ["threads|j=i" => "Number of threads for writing output", { default => 1 }],
				    ["help|h" => "Show this help message."]);

print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 3;

my $genus_data_dir = shift;
my $prop_out_dir = shift;
my $fams_file = shift;

open(FAMS, "<", $fams_file) or die "Cannot open $fams_file: $!";
-d $genus_data_dir or die "Genus data dir $genus_data_dir does not exist\n";
-d $prop_out_dir or die "Propagation output data dir $prop_out_dir does not exist\n";

my $ali_cache_dir = "/dev/shm/ali.$$";
my $ali_out_dir = "$prop_out_dir/alignments/raw.aa";
my $nr_out_dir = "$prop_out_dir/nr";
make_path($ali_out_dir, $nr_out_dir, $ali_cache_dir);

open(LMAP, "<", "$prop_out_dir/map.local_fams") or die "Cannot open $prop_out_dir/map.local_fams: $!";

#
# Lmap is structured as $lmap{$genus}->{$lfam} = $mapped
#
# This is because if we have genera that were not in the last build, we use
# the local family IDs as is.
#
print STDERR "Load map $prop_out_dir/map.local_fams\n";
my %lmap;
while (<LMAP>)
{
    my($g1, $lf1, $g2, $lf2) = /^([^.]+)\.(\d+)\t([^\.]+)\.(\d+)$/;
    defined($g1) or die "Cannot parse $_";

    if ($g1 ne $g2)
    {
	die "Unexpected $g1 != $g2 at $.\n";
    }
    $lmap{$g1}->{$lf1} = $lf2;
}
close(LMAP);

#
# Read family data and populate peg to fam map as well as the set of
# genera seen (so we can effectively do testing on a subset of the fams file).
#

my %peg_to_fams;
my %genera;

print STDERR "Load fams $fams_file\n";
while (<FAMS>)
{
    chomp;

    print STDERR "$.\n" if $. % 1000000 == 0;

    my($fam, $peg, undef, undef, $lfam, $genus) = split(/\t/);

    $peg_to_fams{$peg} = [$fam, $lfam, $genus];
    $genera{$genus} = 1;

}

close(FAMS);

#
# Process our genera to create the NR files
#

if ($opt->threads > 1)
{
    pareach([sort keys %genera], \&process_genus, { Max_Workers => $opt->threads });
}
else
{
    for my $genus (sort keys %genera)
    {
	process_genus($genus);
    }
}

sub process_genus
{
    my($genus) = @_;
    my $ecount = 0;

    print STDERR "Process $genus\n";

    my $genus_lmap = $lmap{$genus};

    my $nr_out = "$nr_out_dir/$genus.fa";
    open(my $nr_fh, ">", $nr_out) or die "Cannot write $nr_out: $!";
    
    my $ali_in = "$genus_data_dir/$genus/alignments.raw.aa.btree";
    my $ali_cache = "$ali_cache_dir/$genus.btree";
    my $ali_out = "$ali_out_dir/$genus.btree";
    my $db = new BerkeleyDB::Btree(-Filename => $ali_in);
    $db or die "Cannot open $ali_in: $! $BerkeleyDB::Error";

    my $new_db;
    #
    # Only write new db if we are mapping. Otherwise just copy
    if ($genus_lmap)
    {
	$new_db = new BerkeleyDB::Btree(-Filename => $ali_cache, -Flags => DB_CREATE);
	$new_db or die "Cannot open $ali_cache: $! $BerkeleyDB::Error";
	print STDERR "Mapping lfams for $genus\n";
    }
    else
    {
	copy($ali_in, $ali_out);
	print STDERR "Not mapping lfams for $genus\n";
    }

    my $cursor = $db->db_cursor();
    my $lfam;
    my $peg_list;

    while ($cursor->c_get($lfam, $peg_list, DB_NEXT) == 0)
    {
	my $nlfam = $genus_lmap ? $genus_lmap->{$lfam} : $lfam;
	
	if (!defined($nlfam))
	{
	    die "Could not map $genus\t$lfam\n";
	}

	open(my $pfh, "<", \$peg_list) or die "Cannot open filehandle on peglist: $!";

	while (my($id, $def, $seq) = read_next_fasta_seq($pfh))
	{
	    my $fam_info = $peg_to_fams{$id};
	    if (!$fam_info)
	    {
		warn "No fam info for $id genus $genus $lfam $nlfam\n";
		last;
	    }
	    my($peg_gfam, $peg_lfam, $peg_genus) = @$fam_info;
	    if ($genus ne $peg_genus || $peg_lfam ne $nlfam)
	    {
		die "Family mismatch on $id: $genus $peg_genus $nlfam $peg_lfam";
	    }
	    $seq =~ s/-//g;
	    write_fasta($nr_fh, [$id, "$peg_gfam $peg_lfam $peg_genus", $seq]);
	}
	$new_db->db_put($nlfam, $peg_list) if $new_db;
    }
    $cursor->c_close();
    $db->db_close();
    if ($new_db)
    {
	$new_db->db_close();
	print STDERR "Copy $ali_cache to $ali_out\n";
	copy($ali_cache, $ali_out) or die "error copying $ali_cache ali_out: $!";
	unlink($ali_cache) or warn "Error unlinking $ali_cache: $!";
    }
    close($nr_fh);
}
