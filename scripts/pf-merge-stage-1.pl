
#
# Given a set of family directories, and the corresponding kmer data directory,
# compute a merge of the families.
#
# We merge based on function. We assign function numbers based on the function.index
# file in the kmer data directory.
#
# For each function, we construct a fasta file containing the representative pegs
# from the families for that function, from each family directory.
# The ids given to the proteins in the fasta files are of the form
# genus-id:fam-id:peg so that join decisions later on can be made strictly on
# the results from the clustering on these fasta files.
#
# The number of pegs to use as representatives is a parameter to this program.
# The method for selecting pegs is also a parameter. 
#
# Given the set of fasta files, we use the kmer distance clustering algorithm to
# compute suggested family merges.
#

use DBI;
use strict;
use Data::Dumper;
use Getopt::Long::Descriptive;
use IO::Handle;
use File::Path 'make_path';
use File::Slurp;
use Statistics::Descriptive;
use File::Basename;
use gjoseqlib;
use Algorithm::Numerical::Shuffle qw /shuffle/;
use POSIX ':sys_wait_h';
use IPC::Run;
use Proc::ParallelLoop;
use List::Util qw(max min);

my %valid_rep_mechs = (random => \&peg_rep_random, 'lfam-alignment' => \&peg_rep_alignment);

my($opt, $usage) = describe_options("%c %o kmer-dir merge-dir famdir famdir [famdir ...]",
				    ['subsystem=s@', "Merge roles that occur in this subsystem only", { default => [] }],
				    ["n-reps|n=i", "Number of pegs per family to use as representative", { default => 10 }],
				    ['rep-type|r=s@', "Mechanism to use to select reps", { default => ['random'] }],
				    ["parallel|p=i", "Number of jobs to run in parallel", { default => 1 }],
				    ["serial", "Force serial execution"],
				    ["log|l=s", "Trace logfile"],
				    ["lfam-db-name=s", "Name of local family database"],
				    ["lfam-db-host=s", "Database host for local family database"],
				    ["lfam-db-user=s", "Username for local family database"],
				    ["lfam-db-password=s", "Username for local family database"],
				    ["dist-dir=s", "Directory into which distances are written"],
				    ["help|h", "Show this help message"]);

print($usage->text), exit if $opt->help;
die($usage->text) if @ARGV < 4;

if ($opt->log)
{
    open(LOG, ">>", $opt->log) or die "Cannot open logfile " . $opt->log . ": $!";
    LOG->autoflush(1);
}
else
{
    open(LOG, ">", "/dev/null");
}

my @rep_mechanisms;
my @lfam_connect;
for my $rep (@{$opt->rep_type})
{
    my $rep_mechanism = $valid_rep_mechs{$rep};
    if (!defined($rep_mechanism))
    {
	my $vals = join(" ", sort keys %valid_rep_mechs);
	die "Invalid rep selection mechanism $rep (valid values are $vals)\n";
    }
    push @rep_mechanisms, $rep_mechanism;

    if ($rep eq 'lfam-alignment')
    {
	#
	# Check connection info.
	#
	@lfam_connect = ("dbi:mysql:host=" . $opt->lfam_db_host . ";dbname=" . $opt->lfam_db_name,
			 $opt->lfam_db_user, $opt->lfam_db_password);
	my $dbh = DBI->connect(@lfam_connect);
	if (!$dbh)
	{
	    die "Cannot connect to lfam database at @lfam_connect: $DBI::errstr";
	}
	$dbh->disconnect();
    }
}

my $kmer_dir = shift;
my $merge_dir = shift;
my @fam_dirs = @ARGV;

my %desired_roles;
my($function_to_index, $index_to_function) = read_function_index($kmer_dir);

for my $subsystem (@{$opt->subsystem})
{
    my $ss_dir = "/vol/core-seed/FIGdisk/FIG/Data/Subsystems/$subsystem";
    open(S, "<",  "$ss_dir/spreadsheet") or die "Cannot open subsystem spreadsheet $ss_dir/spreadsheet: $!";
    
    while (<S>)
    {
	chomp;
	last if m,^//,;
	my($abbr, $role) = split(/\t/);
	$desired_roles{$role} = $abbr;
	if (!defined($function_to_index->{$role}))
	{
	    warn "Role $role not present in function index\n";
	}
	
    }
    close(S);
}

my $tmpdir = "$merge_dir/tmp";
make_path($tmpdir);

my $fa_dir = "$merge_dir/family.fasta";

my $dist_dir = $opt->dist_dir;
if (!$dist_dir)
{
    $dist_dir = "$merge_dir/family.dist";
}

make_path($fa_dir, $dist_dir);

my $all_fasta = "$merge_dir/all.fasta";
my $map_file = "$merge_dir/peg.map";


if (! -s $all_fasta)
{
    my @dirs;
    for my $fam_dir (@fam_dirs)
    {
	my $sz = 0;
	my $seq_subdir = "$fam_dir/nr-seqs";
	opendir(DH, $seq_subdir) or die "Cannot opendir $seq_subdir: $!";
	while (my $f = readdir(DH))
	{
	    next if $f =~ /^\./;
	    my $seqfile = "$seq_subdir/$f";
	    $sz += -s $seqfile;
	}
	closedir(DH);
	push(@dirs, [$fam_dir, $sz]);
    }
    
    my @dirs = sort { $b->[1] <=> $a->[1] } @dirs;
    
    my @work_assigned;

    my $n_buckets = $opt->parallel * 4;
    $work_assigned[$_] = [0, []] for 0..$n_buckets - 1;
    
    for my $ent (@dirs)
    {
	# find slot with smallest load.
	
	my $idxMin = 0;
	$work_assigned[$idxMin]->[0] < $work_assigned[$_]->[0] or $idxMin = $_ for 1 .. $#work_assigned;
	
	my $e = $work_assigned[$idxMin];
	$e->[0] += $ent->[1];
	push(@{$e->[1]}, $ent->[0]);
    }

    make_path("$all_fasta.d", "$map_file.d", "$fa_dir.d");
    my @work;
    for my $idx (0..$#work_assigned)
    {
	my($sz, $list) = @{$work_assigned[$idx]};

	my $all = "$all_fasta.d/$idx";
	my $map = "$map_file.d/$idx";
	my $fa_dirp = "$fa_dir.d/$idx";
	mkdir $fa_dirp;


	push(@work, [$list, $all, $map, $fa_dirp, \@rep_mechanisms, $opt]);
    }

    print Dumper(\@work);

    pareach \@work, sub {
	my($ent) = @_;
	my($list, $all, $map_file, $fa_dirp, $rep_mechanisms, $opt) = @$ent;
	
	#
	# If we are going to use the rep database, connect here in the new subprocess.
	#
	my $lfam_dbh;
	if (@lfam_connect)
	{
	    $lfam_dbh = DBI->connect(@lfam_connect);
	}
	
	my $fh_set = {};

	my $afa_fh;
	open($afa_fh, ">", $all) or die "Cannot open $all: $!";

	for my $fam_dir (@$list)
	{
	    collect_family_fasta($fam_dir, $afa_fh, $fa_dirp, $map_file, $rep_mechanisms, $opt, $fh_set, $lfam_dbh);
	}
	close(AFA);
	
	for my $fh (values(%$fh_set))
	{
	    close($fh) or die "Error closing $fh: $!";
	}
    }, { Max_Workers => $opt->parallel };
    
    #
    # Now, in parallel, merge the pegmap, all.fasta, and individual family files.
    #

    my @work;
    for my $f ($all_fasta, $map_file)
    {
	my @in;
	for my $idx (0..$n_buckets - 1)
	{
	    push(@in, "$f.d/$idx") if -f "$f.d/$idx";
	}
	push(@work, [[[@in], $f]]);
    }
    my @unit;
    for my $fam (0..$#$index_to_function)
    {
	my @in;
	for my $idx (0..$n_buckets - 1)
	{
	    if (-f "$fa_dir.d/$idx/$fam")
	    {
		push(@in, "$fa_dir.d/$idx/$fam");
	    }
	}
	push(@unit, [[@in], "$fa_dir/$fam"]);
	if (@unit >= 10)
	{	    
	    push(@work, [@unit]);
	    @unit = ();
	}
    }
    push(@work, [@unit]);
    
    pareach \@work, sub {
	my($work_list) = @_;

	for my $ent (@$work_list)
	{
	    my($in_list, $out) = @$ent;
	    open(my $out_fh,">", $out) or die "Cannot write $out: $!";
	    for my $in (@$in_list)
	    {
		open(my $in_fh, "<", $in) or die "Cannot read $in: $!";
		my $buf;
		my $n;
		while ($n = read($in_fh, $buf, 32768))
		{
		    print $out_fh $buf;
		}
		
		if (!defined($n))
		{
		    die "Error reading $in: $!";
		}
		close($in_fh);
	    }
	    close($out_fh) or die "Error closing $out for write: $!";;
	}
	for my $ent (@$work_list)
	{
	    my($in_list, $out) = @$ent;
	    unlink(@$in_list);
	}
    
    }, { Max_Workers => $opt->parallel };

}

#
# Given the set of fasta data we've collected, divide into roughly equal groups
# for the amount of work we have.
#
# (LPT algorithm)
#

my @files = grep { $_->[1] > 0 } map { my $s = -s $_; [$_, (($s > 0 && $s < 10000) ? 10000 : $s )] } <$fa_dir/*>;

my @files = grep { my $n = basename($_->[0]); (! -f "$dist_dir/$n") || -s _ == 0 } @files;

my @files = sort { $b->[1] <=> $a->[1] } @files;

my @work_assigned;

my $n_buckets = $opt->parallel * 10;
$work_assigned[$_] = [0, []] for 0..$n_buckets - 1;

for my $ent (@files)
{
    # find slot with smallest load.

    my $idxMin = 0;
    $work_assigned[$idxMin]->[0] < $work_assigned[$_]->[0] or $idxMin = $_ for 1 .. $#work_assigned;

    my $e = $work_assigned[$idxMin];
    $e->[0] += $ent->[1];
    push(@{$e->[1]}, $ent);
}

my @work;
for my $idx (0..$#work_assigned)
{
    my($sz, $list) = @{$work_assigned[$idx]};
    
    if ($sz > 0)
    {
	my $tmp = "$tmpdir/$idx";
	make_path($tmp);
	push(@work, [$list, $sz, $idx, $tmp]);
    }
}
#die Dumper(\@work);
my $real_par = $opt->parallel;

print STDERR "Work assigned:\n";
for (my $i = 0; $i < @work; $i++)
{
    my($list, $sz, $idx, $tmp) = @{$work[$i]};
    print STDERR "$i\t$sz\t$idx\t$tmp\n";
    my $c = 4;
    print STDERR "    ";
    for my $ent (@$list)
    {
	my($file, $fsz) = @$ent;
	my $str = "$file - $fsz ";
	print STDERR $str;
	$c += length($str) + 1;
	if ($c > 100)
	{
	    print STDERR "\n    ";
	    $c = 4;
	}
    }
    print STDERR "\n";
}


#die "Real_par=$real_par\n";

if ($real_par > 1 && !$opt->serial)
{
    pareach \@work, \&do_work, { Max_Workers => $real_par };
}
else
{
    for my $ent (@work)
    {
	do_work($ent);
    }
}

sub do_work
{
    my($item) = @_;
    my($group, $group_size, $index, $tmpdir) = @$item;

    #
    # Load protein data
    #
    
    for my $gent (@$group)
    {
	my($fa_file, $fa_size) = @$gent;
	my $fam = basename($fa_file);

	if (! -s $fa_file)
	{
	    print "Not adding zero size $fa_file\n";
	    next;
	}

	my($kser_pid, $kser_port) = start_kser($tmpdir, $kmer_dir);
	
	my $add_url = "http://localhost:$kser_port/add";
	my $matrix_url = "http://localhost:$kser_port/matrix";
	print "curl to $add_url\n";
	
	$SIG{__DIE__} = sub {
	    kill 1, $kser_pid;
	};
	
	print "Add $fa_file\n";
	my $calls_out = "$tmpdir/calls.$fam";
	$calls_out = "/dev/null"; # don't need this data
	my $rc = system("curl", "-s", "-S", "-o", $calls_out, "-H", "kmer-options: -a -d 1", "--data-binary", "\@$fa_file", $add_url);
	if ($rc != 0)
	{
	    die "Error loading fasta data\n";
	}

	my $fam = basename($fa_file);
	my @cmd = ("kmer_guts_net", "-u", $matrix_url, "-d", 1, "-a");
	my $ok = IPC::Run::run(\@cmd, '<', $fa_file, '>', "$dist_dir/$fam");
	if (!$ok)
	{
	    die "Failed @cmd < $fa_file > $dist_dir/$fam\n";
	}

	print "Killing $kser_pid\n";
	kill 1, $kser_pid;
	print "Waiting on $kser_pid\n";
	waitpid $kser_pid, 0;
	print "Done waiting\n";
    }
}

sub collect_family_fasta
{
    my($fam_dir, $all_fasta_fh, $fa_dir, $map_file, $rep_mechanisms, $opt, $fh_set, $lfam_dbh) = @_;

    open(FAM_DEFS, "<", "$fam_dir/local.family.defs") or die "Cannot open $fam_dir/local.family.defs: $!";
    open(FAM_MEM, "<", "$fam_dir/local.family.members") or die "Cannot open $fam_dir/local.family.members: $!";
    
    my $stat = Statistics::Descriptive::Full->new();

    my $genus = basename($fam_dir);

    my %gname;
    open(G, "<", "$fam_dir/genome.names") or die "Cannot open $fam_dir/genome.names: $!";
    while (<G>)
    {
	chomp;
	my($id, $n) = split(/\t/);
	$gname{$id} = $n;
    }
    close(G);

    open(MAP, ">>", $map_file) or die "Cannot append to $map_file: $!";
    
    my $last;
    my @pegs;

    my %seq;
    my $seq_subdir = "$fam_dir/nr-seqs";
    opendir(DH, $seq_subdir) or die "Cannot opendir $seq_subdir: $!";
    while (my $f = readdir(DH))
    {
	next if $f =~ /^\./;
	my $seqfile = "$seq_subdir/$f";

	open(S, "<", $seqfile) or die "Cannot open $seqfile $!";
	print STDERR "Loading $seqfile\n";
	while (my($id, $def, $seq) = read_next_fasta_seq(\*S))
	{
	    $seq{$id} = $seq;
	}
	close(S);
    }

    #
    # Because we are pulling sequences from the NR, we need to deduplicate the
    # pegs that appear in the families file. The way we're doing this is to only
    # consider pegs that appear as the first synonym in the peg.synonyms file.
    #
    # NB - in newer version, our local.family.members file is constructed from sequences
    # only in the NR.
    #
    # We don't need to actually read peg.synonyms to do this - we need to just
    # only consider proteins that we have sequences for by reading the nr sequences
    # above.
    #

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
    
    while (<FAM_MEM>)
    {
	chomp;
        my @item = split(/\t/);
        my($famId, $peg, $len, $zscore, $gene) = @item;

	next unless $seq{$peg};

	my $finfo = $fam_info{$famId};
	if (%desired_roles)
	{
	    next unless $desired_roles{$finfo->{func}};
	}

	if ($famId != $last->[0])
	{
	    if (@pegs)
	    {
		my($lfamId, $lpeg, $llen, $lzscore) = @$last;
		emit_family_fasta($genus, $all_fasta_fh, $fa_dir, \*MAP, \%seq, $lfamId, $fam_info{$lfamId}->{func},
				  $stat, \@pegs, $rep_mechanisms, $opt, \%gname, $fh_set, $lfam_dbh);
		@pegs = ();
		$stat->clear();
	    }
	    $last = \@item;
	}
	push(@pegs, [$peg, $len]);
	$stat->add_data($len);
    }

    if ($last)
    {
	my($lfamId, $lpeg, $llen, $lzscore) = @$last;
	emit_family_fasta($genus, $all_fasta_fh, $fa_dir, \*MAP, \%seq, $lfamId, $fam_info{$lfamId}->{func},
			  $stat, \@pegs, $rep_mechanisms, $opt, \%gname, $fh_set, $lfam_dbh);
    }
}

sub emit_family_fasta
{
    my($genus, $all_fasta_fh, $fa_dir, $map_fh, $seqs, $fam, $fun, $stat, $pegs, $rep_mechanisms, $opt, $gname, $fh_set, $lfam_dbh) = @_;

    my $fun_idx = $function_to_index->{$fun};
    defined($fun_idx) or die "Cannot map '$fun' to index\n";

    my @reps;
    for my $mech (@$rep_mechanisms)
    {
	@reps = $mech->($stat, $pegs, $fun, $opt, $genus, $fam, $lfam_dbh);
	last if @reps;
    }

    my $fh = $fh_set->{$fun_idx};
    if (!$fh)
    {
	open($fh, ">>", "$fa_dir/$fun_idx") or die "Cannot open $fa_dir/$fun_idx: $!";
	$fh_set->{$fun_idx} = $fh;
    }

    # open(my $fh, ">>", "$fa_dir/$fun_idx") or die "Cannot open $fa_dir/$fun_idx: $!";

    for my $rep (@reps)
    {
	my($gid) = $rep =~ /^fig\|(\d+\.\d+)/;
	print $map_fh "$rep\t$genus\t$fam\t$fun_idx\t$fun\t$gname->{$gid}\n";
	print_alignment_as_fasta($fh, [$rep, $fun, $seqs->{$rep}]);
	print_alignment_as_fasta($all_fasta_fh, [$rep, $fun, $seqs->{$rep}]);
    }
    # close($fh);
}

#
# Find the reps for the given family based on the database of alignment-based reps
#
sub peg_rep_alignment
{
    my($stat, $pegs, $fun, $opt, $genus, $fam, $lfam_dbh) = @_;

    my $res = $lfam_dbh->selectall_arrayref(qq(SELECT af.feature_id
					       FROM alignment a JOIN alignment_feature af ON a.akey = af.alignment
					       WHERE a.genus = ? AND a.family = ? AND a.clean_flag = TRUE AND a.sequence_type = 'D'),
					    undef, $genus, $fam);
    if (!$res || @$res == 0)
    {
	print LOG "No alignment-based reps found for $genus $fam\n" if $opt->log;
	return ();
    }
    else
    {
	my @out = map { $_->[0] } @$res;
	print LOG "Ali found for $genus $fam: @out\n" if $opt->log;

	return @out;
    }
}

sub peg_rep_random
{
    my($stat, $pegs, $fun, $opt, $genus, $fam) = @_;

    my @shuf = shuffle(map { $_->[0] } @$pegs);

    my $n = $opt->n_reps;

    if ($fun eq 'hypothetical protein')
    {
	$n = min($n, 3);
    }

    my @out = @shuf > $n ? @shuf[0..$n-1] : @shuf;
    print LOG "Rand found for $genus $fam: @out\n" if $opt->log;
    return @out;
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

sub start_servers
{
    my($tmpdir, $dataD) = @_;

    my($kser_pid, $kser_port) = start_kser($tmpdir, $dataD);
    
    return ($kser_pid, $kser_port);
}

sub start_kser
{
    my($tmpdir, $dataD) = @_;
    
    #
    # Start the stateful kmer server
    #
    
    my $kser_port_file = "$tmpdir/kser_port";
    my @kser_cmd = ("/scratch/olson/close_kmers/kser",
		    "--listen-port-file", $kser_port_file, 0, $dataD);
    unlink($kser_port_file);
    
    my $kser_pid = fork;
    if (!defined($kser_pid))
    {
	die "fork failed: $!";
    }
    if ($kser_pid == 0)
    {
	print "Starting kser_cmd in $$: @kser_cmd\n";
	exec(@kser_cmd);
	die "Kser did not start $?: @kser_cmd\n";
    }
    
    my $kser_port;
    while (1)
    {
	my $kid = waitpid($kser_pid, WNOHANG);
	if ($kid)
	{
	    die "kser server did not start\n";
	}
	if (-s $kser_port_file)
	{
	    $kser_port = read_file($kser_port_file);
	    chomp $kser_port;
	    if ($kser_port !~ /^\d+$/)
	    {
		kill(1, $kser_pid);
		kill(9, $kser_pid);
		die "Invalid kser port '$kser_port'\n";
	    }
	    last;
	}
	sleep(0.1);
    }
    
    print "Kser running pid $kser_pid on port $kser_port\n";
    
    return ($kser_pid, $kser_port);
}





