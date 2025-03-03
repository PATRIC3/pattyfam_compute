
#
# Given a recall output directory, write an excel report.
#
# We create three tabs:
#  One for changed, non-hypothetical annotations
#  One for non-hypothetical to hypothetical
#  One for hypothetical to non-hypothetical
#
# In the recall data, we will see no-calls. We treat those as hypothetical
# for placing data into the sets above.
#

use strict;
use File::Basename;
use Excel::Writer::XLSX;
use Getopt::Long::Descriptive;
use Data::Dumper;
use SeedUtils;

my($opt, $usage) = describe_options("%c %o recall-dir kmer-dir output-xlsx",
				    ["help|h" => "Show this help message"]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 3;

my $dir = shift;
my $kmer_dir = shift;
my $xlsx = shift;

my $fidx = "$kmer_dir/function.index";
my $deleted_fids_file = "$kmer_dir/deleted_fids";

my %kmer_function;
if (open(F, "<", $fidx))
{
    while (<F>)
    {
	chomp;
	my($idx, $func) = split(/\t/);
	$kmer_function{$func} = $idx;
    }
    close(F);
}
else
{
    die "Cannot open $fidx: $!";
}

my %deleted;
open(F, "<", $deleted_fids_file) or die "Cannot read $deleted_fids_file: $!";
while (<F>)
{
    chomp;
    $deleted{$_} = 1;
}
close(F);

my %tax;
my %name;
if (open(F, "<", "$ENV{KB_TOP}/modules/pattyfam_compute/taxonomy.txt"))
{
    while (<F>)
    {
	chomp;
	my($gid, $tax) = split(/\t/);
	$tax{$gid} = $tax;
    }
    close(F);
}
if (open(F, "<", "$ENV{KB_TOP}/modules/pattyfam_compute/genomes.txt"))
{
    while (<F>)
    {
	chomp;
	my($gid, $name) = split(/\t/);
	$name{$gid} = $name;
    }
    close(F);
}

my $workbook = Excel::Writer::XLSX->new($xlsx);
my $tab_changed = $workbook->add_worksheet("Changed");
my $tab_added = $workbook->add_worksheet("Added");
my $tab_lost = $workbook->add_worksheet("Lost");
my @tabs = ($tab_changed, $tab_added, $tab_lost);

my %tabs = (changed => $tab_changed,
	    added => $tab_added,
	    lost => $tab_lost);

my %next;
$next{changed} = $next{added} = $next{lost} = 2;

my $hypo = 'hypothetical protein';

for my $tab (@tabs)
{
    $tab->set_column(0, 0, 20);
    $tab->set_column(1, 2, 75);
}

my %count;

my $hdr_big_format = $workbook->add_format(bold => 1,
					   size => 14,
					   align => 'left');

my $hdr_format = $workbook->add_format(bold => 1);
my $format_italic = $workbook->add_format(italic => 1);

my %url_count;
my @row_hdrs = qw(ID Old New Score);

for my $rfile (sort <$dir/*>)
{
    open(RF, "<", $rfile) or die "Cannot open $rfile: $!";

    my $gid = basename($rfile);
    while (my($name, $tab) = each %tabs)
    {
	my $row = $next{$name}++;
	$tab->write($row, 0, $gid, $hdr_big_format);
	$tab->write($row, 1, $name{$gid}, $hdr_big_format);
	$row = $next{$name}++;

	for my $i (0 .. $#row_hdrs)
	{
	    $tab->write($row, $i, $row_hdrs[$i], $hdr_format);
	}
    }

    my @dat;
    while (<RF>)
    {
	chomp;
	my @f = split(/\t/);
	next if $deleted{$f[0]};
	push(@dat, \@f);
    }
    close (RF);
	
    for my $ent (sort { SeedUtils::by_fig_id($a->[0], $b->[0]) } @dat)
    {
	my($fid, $old, $old_strip, $new, $fidx, $score) = @$ent;

	$new = $hypo if $new eq '';
	$count{total}++;
	if ($old eq $new)
	{
	    $count{unchanged}++;
	    next;
	}
	my $tab_name = ($old eq $hypo) ? 'added' : (($new eq $hypo) ? 'lost' : 'changed');

	$count{$tab_name}++;
	my $tab = $tabs{$tab_name};
	my $row = $next{$tab_name}++;
	# print "$tab_name $row\n";
	if ($url_count{$tab_name}++ >= 65530)
	{
	    $tab->write($row, 0, $fid);
	} else {
	    $tab->write_url($row, 0,  "https://core.theseed.org/FIG/seedviewer.cgi?page=Annotation&feature=$fid", undef, $fid);
	}
	$tab->write($row, 1, $old, $kmer_function{$old_strip} ? undef : $format_italic);
	$tab->write($row, 2, $new);
	$tab->write($row, 3, $score);
    }
}

print "Counts: total=$count{total} unchanged=$count{unchanged} changed=$count{changed} added=$count{added} lost=$count{lost}\n";
print "% Counts: ";
for my $k (qw(total unchanged changed added lost))
{
    my $p = 100 * $count{$k} / $count{total};
    printf "$k=%.1f ", $p;
}
print "\n";
$workbook->close();    
