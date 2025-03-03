#
# Given output from local family compute on a given genus, write a HTML report.
#

use strict;
use Data::Dumper;
use Getopt::Long::Descriptive;
use DBI;
use Fcntl;
use File::Slurp;
use Template;
use JSON::XS;
use File::Basename;
use Cwd qw(abs_path);
use DB_File;
use P3DataAPI;
use File::Path qw(make_path);
use Cache::FastMmap;



my $json = JSON::XS->new->pretty->canonical;

my($opt, $usage) = describe_options("%c %o kmer-dir family-dir output-dir file-base",
				    ["hypothetical" => "show hypo families only"],
				    ["families-for-pegs=s" => "Show only families including the pegs in this file"],
				    ["write-mappings=s" => "Write mappings of the pegs in --families-for-pegs to this file"],
				    ["full" => "use expanded members list"],
				    ["uniprot=s" => "Use this uniprot mapping file"],
				    ["minimum-length=i" => "minimum protein length to include", { default => 150 }],
				    ["tsv=s" => "write tsv file here"],
				    ["cache=s" => "Function lookup cache file"],
				    ["max-proteins-per-page=i" => "maximum number of protein entries per web page", { default => 10000 }],
				    ["min-fam-size=i" => "minimum family size to include", { default => 5 }],
				    ["help|h" => "Show this help message"]);
print($usage->text), exit 1 if $opt->help;
die($usage->text) if @ARGV != 4;

my $func_cache = Cache::FastMmap->new(share_file => $opt->cache,
				      serializer => '',
				      unlink_on_exit => 0,
				     );

my $api = P3DataAPI->new;

my $kmer_dir = shift;
my $fam_dir = shift;
my $output_dir = shift;
my $file_base = shift;

$fam_dir = abs_path($fam_dir);

make_path($output_dir);

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

my $members = $opt->full ? "$fam_dir/local.family.members.expanded" : "$fam_dir/local.family.members";

open(FAMS, "<", $members) or die "Cannot read $members: $!";
#
# If we are looking for families for given pegs, we need to scan the file once to register the family IDs.
#

if ($opt->write_mappings)
{
    open(FAMS_FOR_PEGS, ">", $opt->write_mappings) or die "cannot write " . $opt->write_mappings . ":$ !";
}
my %chosen_fams;
my %selected_pegs;
if ($opt->families_for_pegs)
{
    print STDERR "Read $members for chosen pegs\n";
    my @pegs = read_file($opt->families_for_pegs);
  
    %selected_pegs = map { chomp $_; $_ => 1 } @pegs;
    while (<FAMS>)
    {
	chomp;
	my($fam_id, $peg, $len, $zscore) = split(/\t/);

	if ($selected_pegs{$peg})
	{
	    $chosen_fams{$fam_id} = 1;
	}
    }
    seek(FAMS, 0, 0);
    print STDERR "Read $members for chosen pegs...done\n";
}

open(DEF, "<", "$fam_dir/local.family.defs") or die "Cannot read $fam_dir/local.family.defs: $!";


my(%fam_fn);
while (<DEF>)
{
    chomp;
    my @ent = split(/\t/);
    my ($fam_id, $fam_fn, $cluster, $mean, $dev, $gene, $fam_size) = @ent;

    next if ($opt->families_for_pegs && !$chosen_fams{$fam_id});

    $fam_fn{$fam_id} = $fam_fn;
    
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

my %reps;
tie %reps, DB_File => "$fam_dir/reps.aa.btree", 0, 0, $DB_BTREE or warn "Can't tie $fam_dir/reps.aa.btree: $!";
my %is_rep;

while (<FAMS>)
{
    chomp;
    my($fam_id, $peg, $len, $zscore) = split(/\t/);
    last if $fam_id eq $last_fam_to_display;

    next if ($opt->families_for_pegs && !$chosen_fams{$fam_id});
    if ($opt->write_mappings && $selected_pegs{$peg})
    {
	print FAMS_FOR_PEGS join("\t", $peg, $fam_id, $fam_fn{$fam_id}) . "\n";
    }

    #
    # Parse fam reps if not done yet.
    #
    if (!exists($is_rep{$fam_id}))
    {
	$is_rep{$_} = 1 foreach split(/\n/, $reps{$fam_id});
    }

    if (my $fam = $id_to_fam{$fam_id})
    {
	my $ent = { id => $peg, len => $len, zscore => $zscore, is_rep => $is_rep{$peg}};
	push(@{$fam->{members}}, $ent);
	push(@to_lookup, $ent);
	if (@to_lookup > 1000)
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
my %fam_index;

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

	push(@{$fam_index{$fam->{fam_fn}}}, [$page_file, $fam->{fam_id}]);

	# coreseed
	if (0) {
	    @ids = @ids[0..99] if @ids > 100;
	    my $s = join("&", map { "feature=$_->{id}" } @ids);
	    $fam->{compare_regions_all} = "http://core.theseed.org/FIG/seedviewer.cgi?page=Regions&$s";
	}
	else
	{
	    my @feats = map { my $id = $_->{id}; $id =~ s/\|/%7C/; $id } grep { $_->{is_rep} } @ids;
	    my $joined = join(",", @feats);
	    my $url = "https://www.bv-brc.org/view/FeatureList/?in(patric_id,($joined))#view_tab=features";
	    $fam->{compare_regions_reps} = $url;

	    my @feats = map { my $id = $_->{id}; $id =~ s/\|/%7C/; $id } @ids;
	    my $joined = join(",", @feats);
	    my $url = "https://www.bv-brc.org/view/FeatureList/?in(patric_id,($joined))#view_tab=features";
	    $fam->{compare_regions_all} = $url;
	}
    }

    my $genus_dir = basename($fam_dir);
    my $tag_name = basename(dirname(dirname($fam_dir)));
	
    my %vars = (page => $dispnum,
		prev => ($dispnum - 1),
		next => ($dispnum + 1),
		prev_link => $prev_page_file,
		next_link => $next_page_file,
		num_pages => scalar @pages,
		uniprot => \%uniprot,
		align_tag => $tag_name,
		align_genus => $genus_dir,
		fams => $fams,
		feature_base_url => "https://www.bv-brc.org/view/Feature/",
	       );

    open(OUT, ">", "$output_dir/$page_file") or die "Cannot write $output_dir/$page_file: $!";

    $tt->process("/home/olson/P3/dev-families/dev_container/modules/pattyfam_compute/scripts/core-fams.tmpl", \%vars, \*OUT)
	or die "Error processing template: " . $tt->error();
    close(OUT);

    my $raw = sprintf "%s-%04d.json", $file_base, $dispnum;
    write_file("$output_dir/$raw", $json->encode(\%vars));
}

open(OUT, ">", "$output_dir/index.html") or die "Cannot write $output_dir/index.html: $!";
my %idx_vars = (pages =>  \@index, fam_index => \%fam_index );

$tt->process("/home/olson/P3/dev-families/dev_container/modules/pattyfam_compute/scripts/core-fams-index.tmpl", \%idx_vars, \*OUT)
	or die "Error processing template: " . $tt->error();
close(OUT);


#
# Look up current BV-BRC annotations
#
sub do_lookup
{
    my($ents) = @_;

    my @to_lookup;
    for my $ent (@$ents)
    {
	my $id = $ent->{id};
	my $func = $func_cache->get($id);
	if (defined($func))
	{
	    $ent->{core_function} = $func;
	    # print STDERR "Found $id = $func\n";
	}
	else
	{
	    push(@to_lookup, $ent);
	}
    }

    # print STDERR "Look up " . Dumper(\@to_lookup);
    my $funcs = $api->function_of([map { $_->{id} } @to_lookup]);

    for my $ent (@to_lookup)
    {
	my $id = $ent->{id};
	my $func = $funcs->{$id};
	$ent->{core_function} = $func;
	$func_cache->set($id, $func);
    }
}

