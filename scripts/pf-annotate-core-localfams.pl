#
# Given output from local family compute on coreseed proteins, emit a HTML report.
#

use strict;
use Data::Dumper;
use Getopt::Long::Descriptive;
use DBI;
use File::Slurp;
use Template;
use JSON::XS;
use File::Basename;
use Cwd qw(abs_path);

my $json = JSON::XS->new->pretty->canonical;

my($opt, $usage) = describe_options("%c %o kmer-dir family-dir output-dir file-base",
				    ["hypothetical" => "show hypo families only"],
				    ["uniprot=s" => "Use this uniprot mapping file"],
				    ["minimum-length=i" => "minimum protein length to include", { default => 150 }],
				    ["tsv=s" => "write tsv file here"],
				    ["override=s" => "write override file here"],
				    ["max-proteins-per-page=i" => "maximum number of protein entries per web page", { default => 10000 }],
				    ["min-fam-size=i" => "minimum family size to include", { default => 5 }],
				    ["password|p=s" => "seed db password"],
				    ["help|h" => "Show this help message"]);
print($usage->text), exit 1 if $opt->help;
die($usage->text) if @ARGV != 4;

my $dbh = DBI->connect("DBI:mysql:database=fig_core_seed;host=chestnut;mysql_ssl=1",
		       "coreseed", $opt->password);

my $kmer_dir = shift;
my $fam_dir = shift;
my $output_dir = shift;
my $file_base = shift;

$fam_dir = abs_path($fam_dir);

open(DEL, "<", "$kmer_dir/deleted_fids") or die "Cannot open $kmer_dir/deleted_fids: $!";
my @del = <DEL>;
chomp @del;
my %del = map { $_ => 1 } @del;

#
# Read uniprot mapping if provided
#
my %uniprot;

if ($opt->uniprot)
{
    open(M, "<", $opt->uniprot) or die "Cannot open uniprot mapping " . $opt->uniprot . ": $!";
    while (<M>)
    {
	chomp;
	my($id, $score, $fn, $ufn) = split(/\t/);
	next unless $id;
	my($sid) = $id =~ /^sp\|([^|]+)/;
	push(@{$uniprot{$fn}}, [$ufn, $sid]);
    }
}

#
# Read family defs to paginate
#

my $cur_page = [];
my $cur_size = 0;
my @pages = ($cur_page);

my $last_fam_to_display;
my %id_to_fam;

open(DEF, "<", "$fam_dir/local.family.defs") or die "Cannot read $fam_dir/local.family.defs: $!";
while (<DEF>)
{
    chomp;
    my @ent = split(/\t/);
    my ($fam_id, $fam_fn, $cluster, $mean, $dev, $gene, $fam_size) = @ent;
    if ($fam_size < $opt->min_fam_size)
    {
	$last_fam_to_display = $fam_id;
	last;
    }

    next if ($opt->hypothetical && !($fam_fn eq 'hypothetical protein' || $fam_fn =~ /core_cid/));
	
    if ($cur_size + $fam_size > $opt->max_proteins_per_page)
    {
	$cur_page = [];
	$cur_size = 0;
	push(@pages, $cur_page );
    }
    $cur_size += $fam_size;
    my $fament =  { fam_id => $fam_id, fam_fn => $fam_fn, cluster => $cluster, mean => $mean, dev => $dev, gene => $gene, fam_size => $fam_size, members => [] };
    $id_to_fam{$fam_id} = $fament;
    push(@$cur_page, $fament);
}
close(DEF);

my @to_lookup;

open(FAMS, "<", "$fam_dir/local.family.members") or die "Cannot read $fam_dir/local.family.members: $!";
while (<FAMS>)
{
    chomp;
    my($fam_id, $peg, $len, $zscore) = split(/\t/);
    last if $fam_id eq $last_fam_to_display;

    if (my $fam = $id_to_fam{$fam_id})
    {
	my $ent = { id => $peg, len => $len, zscore => $zscore};
	push(@{$fam->{members}}, $ent);
	push(@to_lookup, $ent);
	if (@to_lookup > 100)
	{
	    do_lookup(\@to_lookup);
	    @to_lookup = ();
	}
    }
    
}
do_lookup(\@to_lookup);

#
# Generate a HTML document per page, with links to previous & next pages.
#

my $tt = Template->new(ABSOLUTE => 1);

my @index;

for my $pagenum (0..$#pages)
{
    my $dispnum = $pagenum + 1;
    my $page_file = sprintf "%s-%04d.html", $file_base, $dispnum;
    my $prev_page_file = sprintf "%s-%04d.html", $file_base, $dispnum - 1;
    my $next_page_file = sprintf "%s-%04d.html", $file_base, $dispnum + 1;

    my $fams = $pages[$pagenum];

    push(@index, { page_file => $page_file, page => $dispnum, first => $fams->[0], last => $fams->[-1] });

    for my $famidx (0..$#$fams)
    {
	my $fam = $fams->[$famidx];
	$fam->{prev_fam} = $fams->[$famidx - 1]->{fam_id};
	$fam->{next_fam} = $fams->[$famidx + 1]->{fam_id};

	my @ids = @{$fam->{members}};
	@ids = @ids[0..99] if @ids > 100;
	my $s = join("&", map { "feature=$_->{id}" } @ids);
	$fam->{compare_regions_all} = "http://core.theseed.org/FIG/seedviewer.cgi?page=Regions&$s";
    }

    my $genus_dir = basename($fam_dir);
    my $tag_dir = dirname(dirname($fam_dir));
    my $tag_name = basename($tag_dir);
    if ($tag_name eq "core.clusters")
    {
	$tag_name = basename(dirname($tag_dir));
    }
	
    my %vars = (page => $dispnum,
		prev => ($dispnum - 1),
		next => ($dispnum + 1),
		prev_link => $prev_page_file,
		next_link => $next_page_file,
		num_pages => scalar @pages,
		uniprot => \%uniprot,
		align_tag => $tag_name,
		align_genus => $genus_dir,
		feature_base_url => "https://core.theseed.org/FIG/seedviewer.cgi?page=Annotation&feature=",
		fams => $fams);

    open(OUT, ">", "$output_dir/$page_file") or die "Cannot write $output_dir/$page_file: $!";

    $tt->process("/home/olson/P3/dev-families/dev_container/modules/pattyfam_compute/scripts/core-fams.tmpl", \%vars, \*OUT)
	or die "Error processing template: " . $tt->error();
    close(OUT);

    my $raw = sprintf "%s-%04d.json", $file_base, $dispnum;
    write_file("$output_dir/$raw", $json->encode(\%vars));
}

open(OUT, ">", "$output_dir/index.html") or die "Cannot write $output_dir/index.html: $!";
my %idx_vars = (pages =>  \@index );

$tt->process("/home/olson/P3/dev-families/dev_container/modules/pattyfam_compute/scripts/core-fams-index.tmpl", \%idx_vars, \*OUT)
	or die "Error processing template: " . $tt->error();
close(OUT);


#
# Look up coreseed data
#
sub do_lookup
{
    my($ents) = @_;
    my @ids = map { $_->{id} } @$ents;
    my $qs = join( ", ", map { "?" } @$ents);;

    
    my $funcs = $dbh->selectall_hashref(qq(SELECT prot, assigned_function FROM assigned_functions where prot IN ($qs)),
					'prot', undef, @ids);
    $_->{core_function} = $funcs->{$_->{id}}->{assigned_function} foreach @$ents;

    my $aliases = $dbh->selectall_hashref(qq(SELECT id, aliases, gname FROM features f JOIN genome g ON f.genome = g.genome where id IN ($qs)),
					'id', undef, @ids);

    $_->{genome} = $aliases->{$_->{id}}->{gname} foreach @$ents;

}

