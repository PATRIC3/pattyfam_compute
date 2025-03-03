#
# Given output from mmseqs on the coreseed unclassified fasta data (from a pf-compute-local-families run),
# find the clusters and annotate with coreseed functions.
#

use strict;
use Data::Dumper;
use Getopt::Long::Descriptive;
use DBI;
use Template;

my($opt, $usage) = describe_options("%c %o kmer-dir mmseqs-clusters kmer.clusters",
				    ["minimum-length=i" => "minimum protein length to include", { default => 150 }],
				    ["tsv=s" => "write tsv file here"],
				    ["override=s" => "write override file here"],
				    ["password|p=s" => "seed db password"],
				    ["help|h" => "Show this help message"]);
print($usage->text), exit 1 if $opt->help;
die($usage->text) if @ARGV != 3;

my $kmer_dir = shift;
my $fasta = shift;
my $kmer_clusters = shift;

open(FA, "<", $fasta) or die "Cannot open $fasta: $!";

open(DEL, "<", "$kmer_dir/deleted_fids") or die "Cannot open $kmer_dir/deleted_fids: $!";
my @del = <DEL>;
chomp @del;
my %del = map { $_ => 1 } @del;

my %links;

while (<FA>)
{
    chomp;
    my ($id1, $id2) = split /\t/;
    next if $del{$id1} || $del{$id2};
    push @{$links{$id1}}, $id2;
}

my @all_sets;
{
    my @sets = values %links;
    
    @sets = grep { @$_ > 2 } @sets;
    @all_sets = map { { tag => 'B', set => $_} } @sets;
}

print STDERR "Found " . scalar(@all_sets) . " blast sets\n";

#
# Add hypothetical kmer clusters to the sets
#

my @cur;
my $cur_idx;
open(KC, "<", $kmer_clusters) or die "Cannot open $kmer_clusters: $!";

while (<KC>)
{
    chomp;
    my($func, $idx, $fid) = split(/\t/);
    next unless $func eq 'hypothetical protein' || $func =~ /core_cid/;
    if ($idx ne $cur_idx)
    {
	if (@cur > 2)
	{
	    push(@all_sets, { tag => 'K', set => [@cur] });
	}
	@cur =();
	$cur_idx = $idx;
    }
    push(@cur, $fid);
}
if (@cur > 2)
{
    push(@all_sets, { tag => 'K', set => [@cur] });
}
close(KC);
print STDERR "Found " . scalar @all_sets . " sets\n";
my $dbh = DBI->connect("DBI:mysql:database=fig_core_seed;host=arborvitae;mysql_ssl=1",
		       "coreseed", $opt->password);

my $override_fh;
my $tsv_fh;
if ($opt->override)
{
    open($override_fh, ">", $opt->override) or die "Cannot write " . $opt->override . ": $!";
}
if ($opt->tsv)
{
    open($tsv_fh, ">", $opt->tsv) or die "Cannot write " . $opt->tsv . ": $!";
}

my @show_sets;
my $n = 1;
OUTER:
for my $ent (sort { @{$b->{set}} <=> @{$a->{set}} } @all_sets)
{
    my $set = $ent->{set};
    my $tag = $ent->{tag};
    my $in = join(",", map { "?" } @$set);
    my $sizes = $dbh->selectall_hashref(qq(SELECT id, slen FROM protein_sequence_seeks where id IN ($in)),
					'id', undef, @$set);

    my $funcs = $dbh->selectall_hashref(qq(SELECT prot, assigned_function FROM assigned_functions where prot IN ($in)),
					'prot', undef, @$set);

    my $aliases = $dbh->selectall_hashref(qq(SELECT id, aliases, gname FROM features f JOIN genome g ON f.genome = g.genome where id IN ($in)),
					'id', undef, @$set);


    my @rows;
    for my $ent (@$set)
    {
	my $len = $sizes->{$ent}->{slen};
	my $func = $funcs->{$ent}->{assigned_function};
	my $alist = $aliases->{$ent}->{aliases};
	my @alist = split(/,/, $alist);
	my $alias;
	my $no_prefix;
	my $gene_name;
	for my $a (@alist)
	{
            if ($a =~ /^([a-zA-Z]{4})$/)
	    {
		$gene_name = $1;
	    }
	    if ($a =~ /^locus\|(.*)/)
	    {
		$alias = $1;
		last;
	    }
	    elsif ($a =~ /^LocusTag:(.*)/)
	    {
		$alias = $1;
		last;
	    }
	    elsif ($a !~ /[:|]/)
	    {
		$no_prefix = $a;
	    }
	}
	# print STDERR Dumper($ent, \@alist, $no_prefix, $alias);
	$alias //= $no_prefix;
	my $html_alias = $alias // "<i>none</i>";
#	if ($tag eq "B") {
#	    print STDERR "B\t$ent\t$len\n";
#	}
	next OUTER if $len < $opt->minimum_length;
	push(@rows, { id => $ent, len => $len, func => $func, alias => $html_alias, gene_name => $gene_name,
			  raw_alias => $alias,
			  genome => $aliases->{$ent}->{gname} });

    }

    my $new_name = sprintf "hypothetical protein core_cid:$tag-%06d", $n;

    for my $row (@rows)
    {
	if ($tsv_fh)
	{
	    print $tsv_fh join("\t", $n, $row->{id}, $row->{genome}, $row->{len},
			       $row->{gene_name}, $row->{raw_alias}, $row->{func}, $new_name), "\n";
	}
	if ($override_fh)
	{
	    print $override_fh "$row->{id}\t$new_name\n";
	}
    }

    my $s = join("&", map { "feature=$_" } @$set);
    my $compare = "http://core.theseed.org/FIG/seedviewer.cgi?page=Regions&$s";

    push(@show_sets, { compare_url => $compare,
			  tag => $tag,
			  size => scalar @rows,
			  idx => $n, prev_idx => ($n - 1), next_idx => ($n + 1),
			  entries => \@rows });
    $n++;
}

close($tsv_fh) if $tsv_fh;
close($override_fh) if $override_fh;

my %vars = ( sets => \@show_sets );

my $tt = Template->new(ABSOLUTE => 1);
$tt->process("/home/olson/P3/dev-families/dev_container/modules/pattyfam_compute/scripts/core.tmpl", \%vars, \*STDOUT)
or die "Error processing template: " . $tt->error();
