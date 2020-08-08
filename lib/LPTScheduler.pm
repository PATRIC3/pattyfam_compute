#
# Schedule jobs into M machines using the LPT algorithm.
#
# input : n jobs with processing time {p1, p2, . . . , pn} and m machines.
# output: The assignment of jobs to machines.
# 
#  1. Order the jobs in descending order according to their processing times.
#  2. In this order, assign each job to the machine that currently has the
#     least work assigned to it.
#

package LPTScheduler;

use Proc::ParallelLoop;
use strict;
use base 'Class::Accessor';
use List::Util qw(shuffle);
use Data::Dumper;

__PACKAGE__->mk_accessors(qw(m work));

sub new
{
    my($class, $m) = @_;

    $m > 0 or die "LPTScheduler: m must be > 0";
    my $self = {
	m => $m,
	work => [],
    };
    return bless $self, $class;
}

sub add_work
{
    my($self, $item, $time) = @_;
    push(@{$self->work}, [$item, $time]);
}

sub compute_order
{
    my($self) = @_;
    my @sorted = sort { $b->[1] <=> $a->[1] } @{$self->work};
    my @bins;

    return () if @sorted == 0;

    $bins[$_] = [0, []] for 0 .. $self->m - 1;
    for my $w (@sorted)
    {
	my $idxMin = 0;
	$bins[$idxMin]->[0] < $bins[$_]->[0] or $idxMin = $_ for 1..$#bins;
	my $e = $bins[$idxMin];
	$e->[0] += $w->[1];
	push(@{$e->[1]}, $w);
    }

    if ($self->{print_bins})
    {
	for (my $b = 0; $b < @bins; $b++)
	{
	    print STDERR "$b $bins[$b]->[0]\n";
	}
    }
    
    return @bins;
}

sub run
{
    my($self, $bootstrap_cb, $cb, $nprocs) = @_;

    my @bins = $self->compute_order();

    $nprocs //= scalar @bins;

    return if @bins == 0;

    pareach \@bins, sub {
	my($bin) = @_;
	my($size, $list) = @$bin;

	#
	# Shuffle the entries in list so we are not running all
	# of the tiny jobs at the same time in all processes.
	#
	# my @shuffled = shuffle(@$list);

	my $global = $bootstrap_cb->();
	
	for my $ent (@$list)
	{
	    my($item, $time) = @$ent;
	    # print "RUN time=$time\n";
	    $cb->($global, $item);
	}
    }, { Max_Workers => $nprocs };
    
}

1;
