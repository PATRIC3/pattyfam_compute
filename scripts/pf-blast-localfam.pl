# Given a fam id and local fam file in a families genus directory, pull
# the sequences and run an all to all blast.

use File::Temp;
use IPC::Run;
use File::Slurp;
use strict;
use Getopt::Long::Descriptive;
use SeedUtils;
use gjoseqlib;
use Data::Dumper;
use File::Basename;
use DBI;
use Cwd 'abs_path';

my($opt, $usage) = describe_options("%c %o fam-file fam-id",
				    ["align-out=s" => "Write alignment here"],
				    ["ids" => "show ids only"],
				    ["min-size=i" => "Only include sequences at least this long"],
				    ["max-size=i" => "Only include sequences less than or equal to this length"],
				    ["save=s" => "save sequence data here"],
				    ["call-families" => "call patric fams"],
				    ["blast" => "show blast results"],
				    ["mafft" => "show mafft results"],
				    ["compare" => "Show compare regions URL"],
				    ["help|h" => "Show this help message."]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 2;

my $fam_file = shift;
my $fam_id = shift;

my @fam_ids;
if ($fam_id =~ /^\d+$/)
{
    @fam_ids = ($fam_id);
}
else
{
    if (open(IDFILE, "<", $fam_id))
    {
	@fam_ids = <IDFILE>;
	chomp @fam_ids;
    }
}
my %fam_ids = map { $_ => 1 } @fam_ids;

my $fam_dir = dirname(abs_path($fam_file));

my %del;
if (my @del = read_file("deleted.ids"))
{
    chomp @del;
    $del{$_} = 1 foreach @del;
}


open(F, "<", $fam_file) or die "Cannot open $fam_file: $!";

my $from_clusters = $fam_file =~ /\.clusters/;

my @ids;
my $fam_function;

if ($from_clusters)
{
    while (<F>)
    {
	chomp;
	my($fn, $fam, $peg) = split(/\t/);
	if ($fam_ids{$fam})
	{
	    next if $del{$peg};
	    $fam_function = $fn;
	    push(@ids, $peg);
	}
	elsif (@fam_ids == 1 && $fam > $fam_ids[0])
	{
	    last;
	}
    }
}
else
{
    while (<F>)
    {
	chomp;
	my($fam, $peg) = split(/\t/);
	if ($fam_ids{$fam})
	{
	    next if $del{$peg};
	    push(@ids, $peg);
	}
	elsif (@fam_ids == 1 && $fam > $fam_ids[0])
	{
	    last;
	}
    }
}

if ($opt->ids)
{
    print "$_\n" foreach @ids;
    exit;
}
elsif ($opt->compare)
{
    my $s = join("&", map { "feature=$_" } @ids);
    print "http://core.theseed.org/FIG/seedviewer.cgi?page=Regions&$s\n";
    exit;
}

my $tmp = File::Temp->new();
my $assigns = {};

if (1)
{
    my $dbh = DBI->connect("DBI:mysql:database=fig_core_seed_local;host=aspen.mcs.anl.gov", "seed");
    $dbh or die;
    my @to_run = @ids;
    my %all = map { $_ => 1 } @ids;
    while (@to_run)
    {
	my @batch = splice(@to_run, 0, 1000);

	my $prots = join(", ", map { $dbh->quote($_) } @batch);
	my $r = $dbh->selectall_arrayref(qq(SELECT prot, assigned_function
					    FROM assigned_functions
					    WHERE prot IN ($prots)));
	for my $ent (@$r)
	{
	    my($id, $func) = @$ent;
	    $assigns->{$id} = $func;
	    delete $all{$id};
	}
    }
    for my $id (keys %all)
    {
	my $g = genome_of($id);
	if (open(G, "/home/olson/P3/close_kmers/virus.out/anno.orig/$g"))
	{
	    while (<G>)
	    {
		chomp;
		my($x, $func) = split(/\t/);
		if ($x eq $id)
		{
		    $assigns->{$id} = $func;
		    last;
		}
	    }
	    close(G);
	}
    }
}

my %len;

if ($opt->save)
{
    open(S, ">", $opt->save) or die "Cannot write " . $opt->save . ": $!\n";
}

my %by_genome;
for my $id (@ids)
{
    my $genome = genome_of($id);
    $by_genome{$genome}->{$id} = 1;
}

for my $genome (sort keys %by_genome)
{
    my $ids = $by_genome{$genome};

    print STDERR "Load from $genome\n";
    open(SEQS, "<", "$fam_dir/Seqs/$genome") or die "Cannot open $fam_dir/Seqs/$genome: $!";
    while (my($this_id, $def, $seq) = read_next_fasta_seq(\*SEQS))
    {
	if ($ids->{$this_id})
	{
	    my $l = length($seq);
	    if (!(defined($opt->min_size) && $l < $opt->min_size) &&
		!(defined($opt->max_size) && $l > $opt->max_size))
	    {
		print_alignment_as_fasta($tmp, [$this_id, $def, $seq]);
		print_alignment_as_fasta(\*S, [$this_id, $def, $seq]) if $opt->save;
	    }
	    $len{$this_id} = $l;
	    delete $ids->{$this_id};
	    last if %$ids == 0;
	}
    }
    close(SEQS);
	   
}
close($tmp);
close(S) if $opt->save;

if ($opt->blast)
{
    system("makeblastdb", "-dbtype", "prot", "-in", "$tmp");
    
    if (open(BL, "-|", "blastp", "-query", $tmp, "-db", $tmp, "-outfmt", "6 std qlen slen qcovs qcovhsp"))
    {
	while (<BL>)
	{
	    my(@a) = split(/\t/);
	    next if $a[0] eq $a[1];
	    print ;
	}
	close(BL);
    }
}

if ($opt->mafft)
{
    IPC::Run::run(["mafft", "--quiet", "--reorder", "$tmp"],
		  ($opt->align_out ? (">", $opt->align_out) : ()));
}

if ($opt->call_families)
{
    IPC::Run::run(["curl", "--data-binary", "\@$tmp", "http://pear:6100/lookup?find_best_match=1"],
		  ">pipe", \*CURL);
    while (<CURL>)
    {
	chomp;
	my($id, $pgf, $pgf_score, $plf, $plf_score, $fun, $score) = split(/\t/);
	print join("\t", $id, $len{$id}, $assigns->{$id}, ($pgf ||  "          "), $fun), "\n";
    }
}

system("rm $tmp.* 2>/dev/null");
