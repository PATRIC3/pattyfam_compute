
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

my($opt, $usage) = describe_options("%c %o kmer-dir merge-dir",
				    ['html-dir=s', 'If specified, write HTML summaries here'],
				    ['inflation|I=s@', 'Run MCL with this inflation parameter', { default => [2.0]} ],
				    ["parallel|p=i", "Run with this many procs", { default => 1 }],
				    ["help|h", "Show this help message"]);

print($usage->text), exit if $opt->help;
die($usage->text) if @ARGV != 2;

my $kmer_dir = shift;
my $merge_dir = shift;
my $mcl = '/scratch/olson/FuncGroups/bin/mcl';

my $html_dir = $opt->html_dir;

make_path($html_dir) if $html_dir;

my($fn_to_id, $id_to_fn) = read_function_index($kmer_dir);

my @functions;
for my $dist (<$merge_dir/family.dist/*>)
{
    my $f = basename($dist);
    push(@functions, $f);

    #
    # Clip the hypothetical protein group to have a minimum kmer score.
    #
    if (-s $dist > 2_000_000_000)
    {
	my $bak = "$merge_dir/dist.$f.orig";
	if (-s $bak)
	{
	    print STDERR "hypo already fixed ($dist $bak)\n";
	}
	else
	{
	    print STDERR "fix up hypothetical dist file $dist => $bak\n";
	    rename($dist, $bak) or die "cannot rename $dist to $bak: $!";
	    open(my $old, "<", $bak) or die "Cannot open $bak: $!";
	    open(my $new, ">", $dist) or die "Cannot open $dist: $!";
	    while (<$old>)
	    {
		my($id1, $id2, $count, $score) = split(/\t/);
		print $new $_ if $count >= 5;
	    }
	    close($new);
	    close($old);
	}
    }
}

my @items;
for my $inflation (@{$opt->inflation})
{
    my $mcl_dir = "$merge_dir/mcl/$inflation";
    make_path($mcl_dir);
    
    for my $dist (<$merge_dir/family.dist/*>)
    {
	my $f = basename($dist);
	my $fn = $id_to_fn->[$f];
# 	if ($fn eq 'hypothetical protein')
# 	{
# 	    # we will deal with these separately
# 	    next;
# 	}
	if (-s $dist == 0)
	{
	    open(my $fh, ">", "$mcl_dir/$f");
	    close($fh);
	    open(my $fh, ">", "$mcl_dir/$f.mcl.out");
	    print $fh "Wrote zero length file due to zero-length input $dist\n";
	    close($fh);
	}
	else
	{
	    if (! -s "$mcl_dir/$f")
	    {
		my $item = [$f, $fn, $inflation, $dist, "$mcl_dir/$f", ($html_dir ? "$html_dir/$inflation.$f.html" : ())];
		push(@items, [$item, -s $dist]);
	    }
	}
    }
}

@items = sort { $b->[1] <=> $a->[1] } @items;
my @work_assigned;
my $n_buckets = $opt->parallel * 10;
$work_assigned[$_] = [0, []] for 0..$n_buckets - 1;

for my $ent (@items)
{
    # find slot with smallest load.
    
    my $idxMin = 0;
    $work_assigned[$idxMin]->[0] < $work_assigned[$_]->[0] or $idxMin = $_ for 1 .. $#work_assigned;
    
    my $e = $work_assigned[$idxMin];
    $e->[0] += $ent->[1];
    push(@{$e->[1]}, $ent->[0]);
}

my %peginfo;
open(PM, "<", "$merge_dir/peg.map") or die "Cannot open $merge_dir/peg.map: $!";
while (<PM>)
{
    chomp;
    my @a = split(/\t/);
    $peginfo{$a[0]} = \@a;
}
close(PM);

my @work = map { $_->[1] } @work_assigned;
print Dumper(\@work, \@work_assigned);

my $n = @work;
print "Starting $n work items\n";

pareach \@work, sub {
    my $work_list = shift;
    for my $work (@$work_list)
    {
	process(@$work);
    }
}, { Max_Workers => $opt->parallel };

if ($html_dir)
{
    open(I, ">", "$html_dir/index.html") or die "cannot write $html_dir/index.html: $!";
    print I "<table>\n";
}

for my $inflation (@{$opt->inflation})
{
    my $mcl_dir = "$merge_dir/mcl/$inflation";
    
    for my $f (@functions)
    {
	my $mcl = "$mcl_dir/$f";
	my $fun = $id_to_fn->[$f];
	open(F, "<", $mcl) or next;

	my $link = "$inflation.$f.html";
	print I "<tr><th colspan=6><a href='$link' target='_blank'>$fun inflation=$inflation</a></th></tr>\n" if $html_dir;

	my $fidx = 1;
	print  I "<tr><th>family</th><th># genera</th><th># genomes</th><th># lfams</th><th># pegs</th></tr>\n" if $html_dir;
	while (<F>)
	{
	    chomp;
	    my(%genus, %genome, %lfam );
	    my @pegs = split(/\t/);
	    for my $peg (@pegs)
	    {
		my $pi = $peginfo{$peg};
		ref($pi) eq 'ARRAY' or die "No peginfo for '$peg' in file $mcl line $.";
		my(undef, $genus, $lfam, $fidx, $fn, $genome) = @$pi;
		$genus{$genus}++;
		$genome{$genome}++;
		$lfam{$lfam}++;
	    }
	    my $ngenus = scalar keys %genus;
	    my $ngenome = scalar keys %genome;
	    my $nlfam = scalar keys %lfam;
	    my $npegs = scalar @pegs;
	    print I "<tr><td>$fidx</td><td>$ngenus</td><td>$ngenome</td><td>$nlfam</td><td>$npegs</td></tr>\n" if $html_dir;

	    $fidx++;
	}
    }
}
if ($html_dir)
{
    print I "</table>\n";
    close(I);    
}

sub process
{
    my($fun, $function_name, $inflation, $dist_file, $mcl_file, $html_file) = @_;

    print "Process @_\n";
    if (0)
    {
	my $pipe = IO::Handle->new();
	open(my $fh, "<", $dist_file) or die "Cannot open $dist_file: $!";
	my $h = IPC::Run::start([$mcl, "-", "-o", $mcl_file, "--abc", "-I", $inflation],
				"<pipe", $pipe,
				"2>", "$mcl_file.mcl.out");
	$h or die "mcl failed with $?\n";
	
	while (<$fh>)
	{
	    chomp;
	    my($id1, $id2, $c, $s) = split(/\t/);
	    print $pipe "$id1\t$id2\t$s\n";
	}
	close($fh);
	close($pipe);
	my $ok = $h->finish();
	$ok or die "mcl failed with $?\n";

    }
    else
    {
	my $ok = IPC::Run::run(["cut", "-f", "1,2,4", $dist_file],
			       '|',
			       [$mcl, "-", "-o", $mcl_file, "--abc", "-I", $inflation], "2>", "$mcl_file.mcl.out");
	$ok or die "mcl failed with $?\n";
    }
    
    #
    # Postprocess to generate webpage
    #

    return unless $html_file;

    my $fh;
    
    open($fh, ">", $html_file) or die "Cannot open $html_file: $!";
    print $fh "<title>$inflation $function_name</title>\n";
    print $fh "<h1>$function_name</h1>\n";
    print $fh "<h1>Inflation=$inflation</h1>\n";
    
    open(my $m, "<", $mcl_file) or die "Cannot open $mcl_file: $!";

    my @all_pegs;
    my $fidx = 1;
    print $fh "<table>\n" if $fh;
    
    while (<$m>)
    {
	chomp;

	print $fh "<tr><th colspan=3>Family $fidx</td></tr>\n" ;
	print $fh "<tr><th>peg</th><th>lfam</th><th>genome</th></tr>\n";
	my @pegs = split(/\t/);
	for my $peg (@pegs)
	{
	    my $pi = $peginfo{$peg};
	    ref($pi) eq 'ARRAY' or die "No peginfo for '$peg'";
	    my(undef, $genus, $lfam, $fidx, $fn, $genome) = @$pi;
	    print $fh "<tr><td>$peg</td><td>$lfam</td><td>$genome</td></tr>\n";
	}

	push(@all_pegs, @pegs);
	$fidx++;
    }
    print $fh "</table>";

    my $s = join("&", map { "feature=$_" } @all_pegs);
    my $n = @all_pegs;
    print $fh "<a target='_blank' href='http://pseed.theseed.org/seedviewer.cgi?page=Regions&$s'>compared region: $n pegs</a>\n";

    close($fh);
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
