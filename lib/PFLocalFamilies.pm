
package PFLocalFamilies;

use strict;

sub load
{
    my($class, $fam_dir) = @_;

    my $self = {
	by_id => {},
    };

    open(DEF, "<", "$fam_dir/local.family.defs") or die "Cannot open $fam_dir/local.family.defs: $!";
    while (<DEF>)
    {
	chomp;

	my $f = [ split(/\t/) ];
	bless $f, 'PFLocalFamily';
	$self->{by_id}->{$f->id} = $f;
    }
    bless $self, $class;
    return $self;

}

sub family
{
    my($self, $id) = @_;
    return $self->{by_id}->{$id};
}

package PFLocalFamily;

sub id { return $_[0]->[0] };
sub function { return $_[0]->[1] };
sub subfam { return $_[0]->[2] };
sub mean { return $_[0]->[3] };
sub std_gev { return $_[0]->[4] };
sub gene_name { return $_[0]->[5] };
sub size { return $_[0]->[6] };

1;
