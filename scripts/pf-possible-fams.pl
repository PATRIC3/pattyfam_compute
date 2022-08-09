#
# Pick out possible interesting families from the
# representative sequences output.
#
# Parameter is minimum number of proteins. We don't count mulitple proteins from the same genome.
#

use strict;
use Data::Dumper;
use IPC::Run qw(run);
use SeedUtils;
use File::Slurp;
use DBI;

my $dbh = DBI->connect("DBI:mysql:database=fig_core_seed_local;host=aspen.mcs.anl.gov", "seed");
$dbh or die;

my $sth = $dbh->prepare(qq(SELECT assigned_function FROM assigned_functions WHERE prot IN ?));

@ARGV == 2 or die "Usage: $0 deleted-pegs min\n";
my $del_file = shift;
my $n = shift;

my %del;
if (my @del = read_file($del_file))
{
    chomp @del;
    $del{$_} = 1 foreach @del;
}
print "<table border=1>\n";

my @fams;
while (<>)
{
    chomp;
    my @pegs = grep { !$del{$_} } split(/\t/);
    my %g = map { genome_of($_) => 1 } @pegs;
    if (%g >= $n)
    {
	push(@fams, [scalar @pegs, \@pegs]);
    }
}

for my $fam (sort { $b->[0] <=> $a->[0] } @fams)
{
    my($n, $pegs) = @$fam;
    my @pegs = @$pegs;
    
	my $s = join("&", map { "feature=$_" } @pegs);
	my $n = @pegs;
	my $fun;
	#run(["/vol/core-seed/FIGdisk/FIG/bin/function_of"], '<', \$pegs[0], '>', \$fun);
	#$fun =~ s/\S+\s+//;
	
 	# $sth->execute("($pegs[0])");
	# my($fun) = $sth->fetchrow_array();
	my $q = join(", ", map { $dbh->quote($_) } @pegs);
	my $res = $dbh->selectall_hashref(qq(SELECT prot, assigned_function FROM assigned_functions WHERE prot IN ($q)), 'prot');

	next if $res->{$pegs[0]}->{assigned_function} ne 'hypothetical protein';

	my $rlen = $dbh->selectall_hashref(qq(SELECT id, slen FROM protein_sequence_seeks WHERE id IN ($q)), 'id');

	for my $peg (@pegs)
	{
	    if ($peg eq $pegs[0])
	    {
		print "<tr><td><a target='compare' href='http://core.theseed.org/FIG/seedviewer.cgi?page=Regions&$s'>$n pegs</a></td>";
	    }
	    else
	    {
		print "<tr><td></td>";
	    }
	    print "<td>$res->{$peg}->{assigned_function}</td><td>$rlen->{$peg}->{slen}</td><td><a href='https://core.theseed.org/FIG/seedviewer.cgi?page=Annotation&feature=$peg'>$peg</a></td></tr>\n";
	}
}
