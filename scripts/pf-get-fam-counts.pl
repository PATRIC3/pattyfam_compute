#
# Given a file with family id, look up the count of features with that id
#

use strict;
use Data::Dumper;
use P3DataAPI;

my $api = P3DataAPI->new;
my(%feature_count, %genome_count, %product);

my @work;
while (<>)
{
    chomp;
    my($id) = split(/\t/);
    push(@work, $id);
    if (@work > 100)
    {
	query(\@work);
	@work =();
    }
}
query(\@work);

for my $plf (sort keys %genome_count)
{
    print join("\t", $plf, scalar keys %{ $genome_count{$plf} }, $product{$plf}), "\n";
}

sub query
{
    my($list) = @_;
    my $q = join(",", @$list);

    my @res = $api->query("genome_feature", ['eq', 'public', 1], ['in', 'plfam_id', "($q)"], ['select', 'patric_id,product,plfam_id,genome_id']);
    for my $ent (@res)
    {
	$product{$ent->{plfam_id}} = $ent->{product};
	$genome_count{$ent->{plfam_id}}->{$ent->{genome_id}}++;
	$feature_count{$ent->{plfam_id}}++;
    }
}
