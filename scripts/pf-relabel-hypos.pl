#
# Read the family defs to find hypothetical families.
# Read the members file to find non-deleted proteins that are within
# a standard deviation of the mean.
# If there enough, write the pegs with the next hypothetical number.
#

use strict;
use File::Slurp;
use Data::Dumper;

my $min_family_size = 20;

my @defs = read_file("local.family.defs");
my %fams = map { $_->[0] => 1 } grep { $_->[1] eq 'hypothetical protein' && $_->[6] >= $min_family_size } map { chomp; [split(/\t/)] } @defs;
#print STDERR "$_\n" foreach sort keys %fams;

my %del;
if (my @del = read_file("deleted.ids"))
{
    chomp @del;
    $del{$_} = 1 foreach @del;
}

open(M, "<", "local.family.members") or die "cannot read local.family.defs: $!";

my $last_fam;
my @group;
my $next_hypo = 1;

while (<M>)
{
    chomp;
    my($fam, $peg, $len, $dev) = split(/\t/);
    next if $del{$peg};
    next unless $fams{$fam};

    if ($fam ne $last_fam)
    {
	handle_group(\@group);
	@group = ();
	$last_fam = $fam;
    }
    push(@group, $peg) if abs($dev) <= 1;
}
handle_group(\@group);

sub handle_group
{
    my($group) = @_;
    if (@$group > $min_family_size)
    {
	my $id = sprintf("hyp-%06d", $next_hypo);
	$next_hypo++;
	print "$_\t$id\n" foreach @$group;
    }
}
    


