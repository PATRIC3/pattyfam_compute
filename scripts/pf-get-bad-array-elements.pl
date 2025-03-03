#
# Given a job number, return the list of array elements that failed
#

use strict;
@ARGV >= 1 or die "Usage: $0 jobid\n";

my @jobs = map { ("-j", "$_.batch") } @ARGV;


open(P, "-|", "sacct", @jobs, "-o", "jobid%30,state,maxrss", "--units", "g") or die;

my $all;
while (<P>)
{
    next unless /FAIL|CANCEL/;
    print;
    my($jobid) = /^\s+\d+_(\d+)/;
    $all .= "$jobid\n";
}

$all =~ s/(\d+)\n(?=(\d+))/ $1+1==$2 ? "$1-" : $& /ge;
$all =~s/-.*-/-/g;
$all =~ s/\n/,/g;
print "$all\n";
