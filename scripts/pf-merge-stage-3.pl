
use strict;
use Data::Dumper;
use Getopt::Long::Descriptive;
use File::Path 'make_path';
use File::Slurp;
use Statistics::Descriptive;
use File::Basename;
use gjoseqlib;
use IPC::Run;
use Proc::ParallelLoop;
use Graph;
use DB_File;
#
# Use the output of the MCL runs to compute the sets of local families to be merged, and generate
# the global families.
#


my($opt, $usage) = describe_options("%c %o kmer-dir genus-data-dir merge-dir inflation > families",
				    ["genera-from=s", "Merge genera from the given file"],
				    ["genus|g=s@", "Merge the given genus (may be repeated)"],
				    ["help|h", "Show this help message"]);

print($usage->text), exit if $opt->help;
die($usage->text) if @ARGV != 4;

my $kmer_dir = shift;
my $genus_data_dir = shift;
my $merge_dir = shift;
my $inflation = shift;

my($fn_to_id, $id_to_fn) = read_function_index($kmer_dir);

my @genera;
if ($opt->genera_from)
{
    open(F, "<", $opt->genera_from) or die "Cannot open " . $opt->genera_from . ": $!\n";
    while (<F>)
    {
	chomp;
	s/^\s*//;
	s/\s*$//;
	if (! -d "$genus_data_dir/$_")
	{
	    die "Genus directory $merge_dir/$_ does not exist\n";
	}
	push(@genera, $_);
    }
    close(F);
}
if ($opt->genus)
{
    push(@genera, @{$opt->genus});
}

if (@genera == 0)
{
    my @fams = <$genus_data_dir/*/local.family.defs>;
    @genera = map { basename(dirname($_)) } @fams;
}
my %genera = map { $_ => 1 } @genera;
# die Dumper(\@genera)
#
# Load genus family data.
#

print STDERR "Read genus families\n";
my %local_fams;
for my $genus (@genera)
{
    open(FAM_DEFS, "<", "$genus_data_dir/$genus/local.family.defs") or die "Cannot open $genus_data_dir/$genus/local.family.defs: $!";
    open(FAM_MEM, "<", "$genus_data_dir/$genus/local.family.members") or die "Cannot open $genus_data_dir/$genus/local.family.members: $!";

    #
    # First read the defs to map fam ID to role.
    #
    my %fam_info;
    while (<FAM_DEFS>)
    {
	chomp;
	my($id, $func, $subfam, $mean, $stdev, $gene, $size) = split(/\t/);
	$fam_info{$id} = { func => $func, size => $size, subfam => $subfam,
			   mean => $mean, stdev => $stdev, gene => $gene};
    }
    close(FAM_DEFS);

    my $last;
    my $lastfun;
    my $pegs = [];
    while (<FAM_MEM>)
    {
	chomp;
	my($fam, $peg, $len) = split(/\t/);
	my $fun = $fam_info{$fam}->{func};
	if ($fam ne $last)
	{
	    if ($lastfun)
	    {
		# my $fidx = $fn_to_id->{$lastfun};
		# defined($fidx) or die "Cannot map $fun to an index\n";
		$local_fams{$genus, $last} = $pegs;
	    }
	    $pegs = [];
	    $last = $fam;
	    $lastfun = $fun;
	}
	push(@$pegs, [$peg, $len, $fun, $fam]);
    }
    if ($lastfun)
    {
	# my $fidx = $fn_to_id->{$lastfun};
	# defined($fidx) or die "Cannot map $fun to an index\n";
	$local_fams{$genus, $last} = $pegs;
    }
    close(F);
}

print STDERR "Genus families loaded\n";

#
# Read peg mapping file
#

my %peginfo;

if (tie %peginfo, 'DB_File', "$merge_dir/peg.map.btree", O_RDONLY, 0, $DB_BTREE)
{
    print STDERR "Tied to btree\n";
}
else
{
    print STDERR "Reading peg.map\n";

    open(PM, "<", "$merge_dir/peg.map") or die "Cannot open $merge_dir/peg.map: $!";
    while (<PM>)
    {
	chomp;
	my @a = split(/\t/);
	$peginfo{$a[0]} = \@a;
    }
    close(PM);
}

print STDERR "Peg mapping loaded\n";
# print STDERR Dumper(\%local_fams);

#
# For each family in the merge_dir/mcl directory, use the peg map to map to
# determine the set of local families to merge.
#

my %merged;

my %merged_in_fam;

my $gfam = 'GF00000000';
for my $mcl (<$merge_dir/mcl/$inflation/*>)
{
    my $fam = basename($mcl);
    next unless $fam =~ /^\d+$/;

    # HACK
    # next unless $fam == 2196;
    
    open(F, "<", $mcl) or die "Cannot open $mcl: $!";
    print STDERR "$mcl $gfam\n";

    #
    # Read MCL families and map to the local family identifiers they join.
    # Construct a graph where vertices are local family ids and the initial paths
    # are the connections formed by MCL.
    # After reading the MCL connections for a function, the connected subgraphs
    # form the merged families. In other words, we're both merging across local
    # families based on MCL, and merging MCL clusters based on the contents of
    # the local families.
    #

    my $graph = Graph::Undirected->new();
    while (<F>)
    {
	chomp;
	my @pegs = split(/\t/);
	next unless @pegs > 1;
	my $fam = [];
	my %lf;
	for my $peg (@pegs)
	{
	    my $pi = $peginfo{$peg};
	    if (!defined($pi))
	    {
		warn "No peginfo for $pi\n";
		next;
	    }
	    my($rep, $genus, $fam, $fidx, $fun);
	    if (ref($pi))
	    {
		($rep, $genus, $fam, $fidx, $fun) = @$pi;
	    }
	    else
	    {
		($rep, $genus, $fam, $fidx, $fun) = split(/\t/, $pi);
	    }
		
	    my $lfamid = join($;, $genus, $fam);
	    $lf{$lfamid} = 1;
	}
	$graph->add_path(keys %lf);
    }
    close(F);

    print STDERR Dumper($graph->connected_components);
    for my $group ($graph->connected_components)
    {
	print STDERR "merge\t@$group\n";
	my %gmerge;
	my $found;
	for my $lf (@$group)
	{
	    my($genus, $kgfam) = split(/$;/, $lf);
	    $gmerge{$genus}++;
	}

	my $ng = keys %gmerge;
	my $nf = @$group;

	for my $lf (@$group)
	{
	    my($genus, $kgfam) = split(/$;/, $lf);
	    my $pl = delete $local_fams{$lf};
	    
	    if (ref($pl))
	    {
		for my $pi (@$pl)
		{
		    my($peg, $len, $fun, $fam) = @$pi;
		    print join("\t", $gfam, $nf, $ng, $peg, $len, $fun, $fam, $genus, $kgfam), "\n";
		}
		$merged_in_fam{$lf} = $gfam;
		$found++;
	    }
	    else
	    {
		warn "No local fam for $lf ($merged_in_fam{$lf})\n";
	    }
	}
	$gfam++ if $found;
    }
    next;

    while (<F>)
    {
	chomp;
	my @pegs = split(/\t/);
	next unless @pegs > 1;
	my %merge_these;
	my %gmerge;
	for my $peg (@pegs)
	{
	    my $pi = $peginfo{$peg};
	    if (!ref($pi))
	    {
		warn "No peginfo for $pi\n";
		next;
	    }
	    # print $map_fh "$rep\t$genus\t$fam\t$fun_idx\t$fun\t$gname->{$gid}\n";
	    my($rep, $genus, $fam, $fidx, $fun) = @$pi;
	    next unless $genera{$genus};
	    $gmerge{$genus}++;
	    $merge_these{$genus, $fam} = 1;
	}
	my $n = keys %merge_these;
	if ($n > 1)
	{
	    print STDERR "Merging:\n";
	    for my $kk (keys %merge_these)
	    {
		my($genus, $fam) = split(/$;/, $kk);
		
		my $lfams = ref($local_fams{$kk}) ? join(" ", map { $_->[0] } @{$local_fams{$kk}}) : 'NA';

		print STDERR "  $genus\t$fam\t$merged_in_fam{$kk}\t$lfams\n";
	    }
	    my $found;
	    for my $k (keys %merge_these)
	    {
		my($genus, $kgfam) = split(/$;/, $k);
		my $pl = delete $local_fams{$k};
		my $ng = keys %gmerge;
		if (ref($pl))
		{
		    for my $pi (@$pl)
		    {
			my($peg, $len, $fun, $fam) = @$pi;
			print join("\t", $gfam, $n, $ng, $peg, $len, $fun, $fam, $genus, $kgfam), "\n";
		    }
		    $merged_in_fam{$k} = $gfam;
		    $found++;
		}
		else
		{
		    warn "No local fam for $k ($merged_in_fam{$k})\n";
		}
	    }
	    $gfam++ if $found;
	}
    }
}

for my $lk (sort keys %local_fams)
{
    my $pl = $local_fams{$lk};
    my($genus, $kgfam) = split(/$;/, $lk);

    for my $pi (@$pl)
    {
	my($peg, $len, $fun, $fam) = @$pi;
	print join("\t", $gfam, 1, 1, $peg, $len, $fun, $fam, $genus, $kgfam), "\n";
    }
    $gfam++;
}

sub read_function_index
{
    my($dir) = @_;

    my $idx = {};
    my $from_id = [];

    open(F, "<", "$dir/function.index") or die "Cannot open $dir/function.index:$ !";
    while (<F>)
    {
	chomp;
	my($i, $fun) = split(/\t/);
	$idx->{$fun} = $i;
	$from_id->[$i] = $fun;
    }
    close(F);
    
    return($idx, $from_id);
}
