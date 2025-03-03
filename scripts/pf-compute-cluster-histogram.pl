use strict;

my %hist;
while (<>)
{
    chomp;
    my $c = grep { /fig/ } split(/\t/);
    $hist{$c}++;
}

for my $k (sort { $a <=> $b} keys %hist)
{
    print "$k\t$hist{$k}\n";
}
