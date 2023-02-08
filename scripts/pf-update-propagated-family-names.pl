#
# Read the output logs from p3x-propagate-family-names runs and apply the renames
# to the input families file.
#
# We will read all logs (one for global propagation, one for each genus-specific
# local propagation) and cache the data, and apply the changes to the input families
# file(s). 
#
# Emit a record of the family splits and merges.
#

use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use IO::Handle;
use File::Basename;

my($opt, $usage) = describe_options("%c %o [log-files] ::: [family-files]",
				    ["output-directory|o=s" => "If specified write output files to this directory named using the name of the input family file. Otherwise write to standard output"],
				    ["max-id-file=s" => "File containing family max used ID values (required)"],
				    ["local-fam-map-file=s" => "Write this file with local family id mappings"],
				    ["global-fam-map-file=s" => "Write this file with global family id mappings"],
				    ["help|h" => "Show this help message"],
				    );

print($usage->text), exit 1 if $opt->help;
die($usage->text) if @ARGV < 3 || !$opt->max_id_file;

#
# Find logfiles before the ::: indicator and family files after the ::: indicator.
#
my @log_files;
my @fam_files;
my $list = \@log_files;
for my $ent (@ARGV)
{
    if ($ent eq ':::')
    {
	if ($list == \@fam_files)
	{
	    die "Only one instance of ::: allowed in the argument list\n";
	}
	else
	{
	    $list = \@fam_files;
	}
	next;
    }
    if (! -f $ent)
    {
	die "File $ent not readable\n";
    }
    push(@$list, $ent);
}
@log_files or die "At least one logfile must be provided\n";
@fam_files or die "At least one family file must be provided\n";

main();

sub main
{
    my %lf_id_map;
    my %gf_id_map;
    my @splits;
    my @joins;

    #
    # Read the max ID file; increment each by one to create the next id to be assigned
    # for our global and local families.
    #

    my %next_id;
    my %max_new_id;
    open(MID, "<", $opt->max_id_file) or die "Max id file " . $opt->max_id_file . " cannot be opened: $!";
    while (<MID>)
    {
	chomp;
	my($fam, $val) = split(/\t/);
	$next_id{$fam} = 1 + $val;
    }
    close(MID);
    
    for my $lf (@log_files)
    {
	read_log_file($lf, \%lf_id_map, \%gf_id_map, \@splits, \@joins, \%next_id, \%max_new_id);
    }

#    print Dumper(\@splits, \@joins);
#    print Dumper(\%lf_id_map, \%gf_id_map);

    if ($opt->local_fam_map_file)
    {
	open(F, ">", $opt->local_fam_map_file);
	while (my($k, $v) = each %lf_id_map)
	{
	    print F "$k\t$v\n";
	}
	close(F);
    }

    if ($opt->global_fam_map_file)
    {
	open(F, ">", $opt->global_fam_map_file);
	while (my($k, $v) = each %gf_id_map)
	{
	    print F "$k\t$v\n";
	}
	close(F);
    }

    my %fam_fun;
    for my $ff (@fam_files)
    {
	process_fam_file($ff, $opt->output_directory, \%lf_id_map, \%gf_id_map, \%fam_fun, \%next_id);
    }
}

sub process_fam_file
{
    my($file, $out_dir, $lf_id_map, $gf_id_map, $fam_fun, $next_id) = @_;

    open(F, "<", $file) or die "Cannot open $file: $!";
    my $out_fh = \*STDOUT;
    if ($out_dir)
    {
	my $out_file = $out_dir . "/" . basename($file);
	$out_fh = IO::Handle->new();
	open($out_fh, ">", $out_file) or die "Cannot write $out_file: $!";
    }

    while (<F>)
    {
	chomp;
	my($gf, $n, $ng, $peg, $len, $fun, $fam, $genus) = split(/\t/);
	my $lf = "$genus.$fam";

	my $ngf = $gf_id_map->{$gf};
	my $nlf = $lf_id_map->{$lf};
	$ngf or die "Cannot map $gf\n";

	#
	# If there is no local family mapping, see if this is a new genus.
	# A new genus will not have an entry in the $next_id hash that maintains the
	# next id to be assigned for a local family. In this case we just carry across
	# the newly created ID from the incoming family file.
	#

	if (!$nlf)
	{
	    if (exists($next_id->{$genus}))
	    {
		die "Cannot map $lf\n";
	    }
	    else
	    {
		$nlf = $lf;
	    }
	}
	$fam_fun->{$ngf} = $fam_fun->{$nlf} = $fun;
	my($ngenus, $nfam) = $nlf =~ /^(.*?)\.(\d+)$/;
	print $out_fh join("\t", $ngf, $n, $ng, $peg, $len, $fun, $nfam, $ngenus, $nfam), "\n";
    }

    if ($out_dir)
    {
	close($out_fh);
    }
}

sub read_log_file
{
    my($file, $lf_id_map, $gf_id_map, $splits, $joins, $next_id, $max_new_id) = @_;

    open(F, "<", $file) or die "Cannot open $file: $!";

    my $find_unmapped;
    while (<F>)
    {
	chomp;
	if (/^Unmapped\s+new/)
	{
	    $find_unmapped++;
	    last;
	}
	elsif (/^(.*?)\s+NOW\s+(.*?)(\s+weight=.*)?$/)
	{
	    my($new, $old) = ($1, $2);
	    my $map;
	    my $genus;
	    if ($new =~ /^(.*?)\.(\d+)/)
	    {
		$genus = $1;
		$map = $lf_id_map;
	    }
	    else
	    {
		$genus = "GLOBAL";
		$map = $gf_id_map;
	    }
	    # print "'$_' '$new' '$old' '$genus'\n";
	    my $old_upd = $old;
	    if ($old =~ /NEW_(\d+)/)
	    {
		my $val = $next_id->{$genus} + $1;
		$max_new_id->{$genus} = $1 if $1 > $max_new_id->{$genus};
		if ($genus eq 'GLOBAL')
		{
		    $old_upd = sprintf("GF%08d", $val);
		}
		else
		{
		    $old_upd = "$genus.$val";
		}
	    }

	    defined($map->{$new}) && die "map entry for $new already exists with $old";
	    # print "MAP '$new' '$old_upd'\t'$_'\n";
	    $map->{$new} = $old_upd;
	}
	elsif (/^SPLIT\s+O\s+(.*?)\s+=>\s+N\s+(.*)$/)
	{
	    my @newids;
	    my $oid = $1;
	    my $newids = $2;
	    if ($oid =~ /^(.*?)\.(\d+)$/)
	    {
		my $g = $1;
		@newids = $newids =~ /($g\.\d+)/mg;
	    }
	    else
	    {
		@newids = split(/\s+/, $2);
	    }
	    push(@$splits, [$oid, [@newids]]);
	}
	elsif (/^JOIN\s+(.*?)\s+=>\s+(.*)$/)
	{
	    my @oids = split(/\s+/, $1);
	    my $nid = $2;
	    push(@$joins, [[@oids], $nid]);
	}
	else
	{
	    warn "Unmatched: $_\n";
	}
    }

    if ($find_unmapped)
    {
	my %next_new;
	for my $genus (keys %$max_new_id)
	{
	    $next_new{$genus} = 1 + $max_new_id->{$genus};
	}
	while (<F>)
	{
	    if (/^\t([^\t]+)\t/)
	    {
		my $new = $1;
		#
		# This is a new id that was not mapped to an old ID; we need to assign a new
		# family id.
		#
		my $map;
		my $old;
		my $genus;
		if ($new =~ /^(.*)\.(\d+)$/)
		{
		    $genus = $1;
		    my $id = $next_new{$genus}++ + $next_id->{$genus};
		    $map = $lf_id_map;
		    $old = "$genus.$id";
		}
		else
		{
		    $genus = "GLOBAL";
		    my $id = $next_new{$genus}++ + $next_id->{$genus};
		    $map = $gf_id_map;
		    $old = sprintf("GF%08d", $id);
		}
		defined($map->{$new}) && die "map entry for $new already exists with $old in $file at $.";
		$map->{$new} = $old;
	    }
	}
    }
}



sub is_localfam_id
{
    my($id) = @_;
    return $id =~ /^\S+\.\d+$/;
}
