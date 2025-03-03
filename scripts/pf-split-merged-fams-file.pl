#
# Given a fams file as the output of pfpf-merge-stage-3, split it into a defs and members file
#

use strict;
use Getopt::Long::Descriptive;

my($opt, $usage) = describe_options("%c %o fams-file defs-file members-file singleton-defs singleton-members",
				    ["help|h" => "Show this help message."]);

print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 5;

my $fams_file = shift;
my $defs_file = shift;
my $members_file = shift;
my $singleton_defs_file = shift;
my $singleton_members_file = shift;

open(F, "<", $fams_file) or die "Cannot read $fams_file: $!";
open(DEFS, ">", $defs_file) or die "Cannot read $defs_file: $!";
open(MEMBERS, ">", $members_file) or die "Cannot read $members_file: $!";
open(SDEFS, ">", $singleton_defs_file) or die "Cannot read $singleton_defs_file: $!";
open(SMEMBERS, ">", $singleton_members_file) or die "Cannot read $singleton_members_file: $!";

my $cur_fam;
my %cur_genera;
my $cur_fun;
my $cur_peg_count = 0;
my $cur_txt = '';
while (<F>)
{
    chomp;
    my($fam, $localfam_count, $genus_count, $peg, $len, $fun, $localfam, $genus) = split(/\t/);

    if ($fam ne $cur_fam)
    {
	if ($cur_fam)
	{
	    my @genera = sort keys %cur_genera;
	    my $defline = join("\t", $cur_fam, $cur_fun, $cur_peg_count, scalar @genera, join(",", @genera)) .  "\n";
	    
	    if ($cur_peg_count == 1)
	    {
		print SDEFS $defline;
		print SMEMBERS $cur_txt;
	    }
	    else
	    {
		print DEFS $defline;
		print MEMBERS $cur_txt;
	    }
	}
	$cur_fam = $fam;
	$cur_fun = $fun;
	%cur_genera = ();
	$cur_peg_count = 0;
	$cur_txt = '';
    }
    $cur_genera{$genus}++;
    $cur_peg_count++;
    # print MEMBERS join("\t", $fam, $peg, $len, $fun, $localfam, $genus), "\n";
    $cur_txt .= join("\t", $fam, $peg, $len, $fun, $localfam, $genus) . "\n";
}

if ($cur_fam)
{
    my @genera = sort keys %cur_genera;
    my $defline = join("\t", $cur_fam, $cur_fun, $cur_peg_count, scalar @genera, join(",", @genera)) .  "\n";
    
    if ($cur_peg_count == 1)
    {
	print SDEFS $defline;
	print SMEMBERS $cur_txt;
    }
    else
    {
	print DEFS $defline;
	print MEMBERS $cur_txt;
    }
}
close(DEFS) or die "Error closing $defs_file: $!";
close(MEMBERS) or die "Error closing $members_file: $!";
close(SDEFS) or die "Error closing $singleton_defs_file $!";
close(SMEMBERS) or die "Error closing $singleton_members_file: $!";
