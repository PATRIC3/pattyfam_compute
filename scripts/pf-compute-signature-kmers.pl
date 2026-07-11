#!/usr/bin/env perl

use Data::Dumper;
use strict;
use JSON::XS;
use File::Slurp;
use Cwd qw(abs_path);
use Getopt::Long::Descriptive;
use File::Temp;
use File::Copy;
use File::Basename;
use File::Path qw(make_path);
use IPC::Run qw(run);
use gjoseqlib;

my($opt, $usage) = describe_options("%c %o kmer-dir",
				    ["debug_mode", hidden => { one_of => [
								["gdb", "Run with gdb"],
								["valgrind", "Run with valgrind"]
									  ] }],
				    ['threads=i', "Threads to use", { default => 8 }],
				    ['nudb-file=s', "NuDB database file"],
				    ['temp-data=s', "Use this as the temp data dir"],
				    ['recall-output=s', 'Recall output', { default => 'recall.out' }],
				    ["protect-subsystem-roles", "Keep kmers for functions that contain a subsystem role"],
				    ["function-override-file=s", "Use this file for function overrides"],
				    ['organism=s@', "Limit to organism(s)", { default => [] }],
				    ['virus-dir=s', "Add a hierarchy of virus genomes"],
				    ['keep-viral-roles', 'Keep all virus functions in generated kmers'],
				    ['out-dir=s', "Write output here", { default => "out" }],
				    ['phages', 'Use phages'],
				    ['dry-run', "Don't run, just build data and show command"],
				    ['dry-run-file=s', "Save the build command here"],
				    ["help|h", "Show help message"],
				   );


print($usage->text), exit 1 if $opt->help;
die($usage->text) if @ARGV != 1;

#my $bin = "/home/olson/P3/close_kmers";

my $kmer_dir = abs_path(shift);
-d $kmer_dir or die "Kmer dir $kmer_dir does not exist\n";

my $seed_base = "/coreseed";
my $seed_bin = "$seed_base/seed-build/dev_container/bin";
my $seed = "$seed_base/FIGdisk";

my @orgs = @{$opt->organism};

#my @orgs = qw(83333.1 198215.1 83334.1 405955.9 316407.3);

my @betacoronavirus = qw(11128.886 31631.1194 290028.403 502102.6 694006.94 694007.76
    694008.40 694009.690 1160968.8 1263720.4 1335626.1148 1385427.6 1541205.6 1590370.8 1892416.31 2697049.8112);

#@orgs = @betacoronavirus[0..5];

my @good_roles;
if ($opt->protect_subsystem_roles)
{
    if (! -s "$kmer_dir/subsystem.roles")
    {
	my $ok = run(["$seed_bin/write-subsystem-roles", "$kmer_dir/subsystem.roles"]);
	$ok or die "Error writing subsystem roles\n";
    }
    @good_roles = ("--good-roles", "$kmer_dir/subsystem.roles");
}

if ($opt->phages)
{
    open(F, "-|", "$seed_bin/fig", "genomes") or die;
    @orgs = ();
    while (<F>)
    {
	if (/^(\d+\.\d+).*phage/i)
	{
	    push(@orgs, $1);
	}
    }
    close(F);
}
elsif (@orgs == 0)
{
    open(F, "-|", "$seed_bin/fig", "genomes") or die "Cannot run $seed_bin/fig";
    while (<F>)
    {
	if (/^(\d+\.\d+)/)
	{
	    push(@orgs, $1);
	}
    }
    close(F);
}

open(DEL, ">", "$kmer_dir/deleted_fids") or die "Cannot write $kmer_dir/deleted_fids: $!\n";

if (! -f "$kmer_dir/pegs.to.ignore")
{
    my $ok = run(["$seed_bin/find-overlapping-phage-pegs"], ">", "$kmer_dir/pegs.to.ignore");
    $ok or die "Error $? running $seed_bin/find-overlapping-phage-pegs";
}

if (open(PH, "<", "$kmer_dir/pegs.to.ignore"))
{
    # print STDERR "Reading $kmer_dir/pegs.to.ignore\n";
    while (<PH>)
    {
	if (my($id) = /^([^\t]+)\t/)
	{
	    print DEL "$id\n";
	}
    }
    close(PH);
}

my $min_reps_required = 3;
$min_reps_required = @orgs if @orgs < $min_reps_required;

my $orgdir = "$seed/FIG/Data/Organisms";

my @fasta_files;
my @assignment_files;

my $dir = $opt->temp_data // File::Temp->newdir;
my $def_dir = "$dir/defs";
my $fasta_dir = "$dir/fasta";
my $keep_fasta_dir = "$dir/keep_fasta";

make_path($def_dir, $fasta_dir, $keep_fasta_dir);


if ($opt->virus_dir)
{
    my %viral_roles;
    for my $fasta_file (glob($opt->virus_dir . "/fasta/*"))
    {
	my $org = basename($fasta_file);
	my $anno_file = $opt->virus_dir . "/anno/$org";
	if (! -f $anno_file)
	{
	    die "Missing anno file $anno_file\n";
	}
	copy($fasta_file, "$fasta_dir/$org") or die "cannot symlink $fasta_file $fasta_dir/$org: $!";
	copy($anno_file, "$def_dir/$org") or die "cannot symlink $anno_file $def_dir/$org: $!";

	if ($opt->keep_viral_roles)
	{
	    open(VR, "<", $anno_file);
	    while (<VR>)
	    {
		chomp;
		my($id, $fn) = split(/\t/);
		$viral_roles{$fn} = 1;
	    }
	    close(VR);
	}
    }
    if ($opt->keep_viral_roles)
    {
	open(VRF, ">", "$kmer_dir/viral.roles");
	print VRF "$_\n" foreach sort keys %viral_roles;
	close(VRF);
	push(@good_roles, "--good-roles", "$kmer_dir/viral.roles");
    }
}

for my $org (@orgs)
{
    my $fasta = "$orgdir/$org/Features/peg/fasta";
    my $assign = "$orgdir/$org/assigned_functions";
    if (-f $fasta && -f $assign)
    {
	# copy($fasta, "$fasta_dir/$org") or die "cannot symlink $fasta $fasta_dir/$org: $!";
	# copy($assign, "$def_dir/$org") or die "cannot symlink $assign $def_dir/$org: $!";

	open(NFA, ">", "$fasta_dir/$org") or die "cannot write $fasta_dir/$org: $!";
	open(NAS, ">", "$def_dir/$org") or die "cannot write $def_dir/$org: $!";

	my %del;
	if (open(DPEGS, "<", "$orgdir/$org/Features/peg/deleted.features"))
	{
	    # print STDERR "Reading $orgdir/$org/Features/peg/deleted.features\n";
	    while (<DPEGS>)
	    {
		if (/^(\S+)/)
		{
		    print DEL "$1\n";
		    $del{$1} = 1;
		}
	    }
	    close(DPEGS);
	}
	open(FA, "<", $fasta) or die "Cannot open $fasta: $!";
	while (my($id, $def, $seq) = read_next_fasta_seq(\*FA))
	{
	    if (!$del{$id})
	    {
		write_fasta(\*NFA, [$id, $def, $seq]);
	    }
	}
	close(FA);
	close(NFA);

	#
	# We need to read the assigned functions file and create the mapping
	# from id to function. This is because we can likely have multiple
	# entries for the same ID as the ID has been given new assignments.
	# The C++ code currently does not handle that case and as such
	# may compute some statistics incorrectly.
	#
	my %id_func;
	open(AS, "<", $assign) or die "Cannot open $assign: $!";
	while (<AS>)
	{
	    chomp;
	    my($id, $func) = split(/\t/);
	    next if $del{$id};
	    next unless $id =~ /\.peg\./;
	    $id_func{$id} = $func;
	    # print NAS $_ unless $del{$id};
	}
	close(AS);
	for my $id (sort keys %id_func)
	{
	    print NAS "$id\t$id_func{$id}\n";
	}
	close(NAS);
	    
    }
    else
    {
	die "Missing one of $fasta $assign\n";
    }
}
close(DEL);

#
# Check out the repo with additional fasta for definition and the ignored-functions file.
#

my $rc = system("git", "clone", 'https://github.com/TheSEED/kmer-build-data', "$kmer_dir/kmer-build-data");
$rc == 0 or die "git checkout failed with rc=$rc\n";
my @other = <$kmer_dir/kmer-build-data/*.txt>;
my $ignored_functions = "$kmer_dir/kmer-build-data/functions.to.ignore";
for my $f (@other)
{
    my $b = basename($f);
    symlink($f, "$keep_fasta_dir/$b") or die "cannot symlink $f $keep_fasta_dir/$b";
}

#
# Save the functions associated with dlits.
#

my $ok = run(["$seed_bin/write-dlit-functions"], ">", "$kmer_dir/dlit.functions");
$ok or die "error $? running $seed_bin/write-dlit-functions";
push(@good_roles, "--good-roles", "$kmer_dir/dlit.functions");


make_path($opt->recall_output);

my $exe = "kmers-build-signatures";
my @nudb;
if ($opt->nudb_file)
{
    @nudb = ("--nudb-file", $opt->nudb_file);
}

my @override;
if ($opt->function_override_file)
{
    @override = ('--function-override-file', $opt->function_override_file);
}

my @cmd = ($exe,
	   "-D", $def_dir,
	   "-F", $fasta_dir,
	   "-K", $keep_fasta_dir,
	   "--kmer-data-dir", $opt->out_dir,
	   "--deleted-features-file", "$kmer_dir/deleted_fids",
	   "--ignored-functions-file", $ignored_functions,
	   "--min-reps-required", $min_reps_required,
	   "--n-threads", $opt->threads,
	   "--perfect-hash", "$kmer_dir/kmer_data.mph",
	   "--perfect-hash-data", "$kmer_dir/kmer_data.dat",
	   @override,
	   @nudb,
	   @good_roles,
	   );

if ($opt->debug_mode eq 'valgrind')
{
    unshift(@cmd, "/home/olson/P3/valgrind/install/bin/valgrind");
}
elsif ($opt->debug_mode eq 'gdb')
{
    print "command: @cmd\n";
    @cmd = ("/disks/patric-common/runtime/gcc-9.3.0/bin/gdb", $exe);
}

if ($opt->dry_run)
{
    print "Would run: @cmd\n";
    if ($opt->dry_run_file)
    {
	write_file($opt->dry_run_file, JSON::XS->new->pretty->canonical->encode(\@cmd));
    }
}
else
{
    print "Run: @cmd\n";

    my $rc = system(@cmd);
    
    if ($rc != 0)
    {
	die "build failed with $rc: @cmd\n";
    }
}
