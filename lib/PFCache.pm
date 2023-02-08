
package PFCache;

use base 'Exporter';

use File::Path qw(make_path);

our @EXPORT_OK = qw(compute_cache_path read_pegmap);



#
# Some common cache utility code.
#

sub compute_cache_path
{
    my($cache_base_dir, $gid) = @_;
    my $xgid = "000" . $gid;
    my($subdir) = $xgid =~ /(\d{3})\./;
    my $path = "$cache_base_dir/$subdir/$gid";
    make_path($path);
    return $path;
}

sub read_pegmap
{
    my($merge_dir, $genus_dir, $peginfo) = @_;

    if (open(PM, "<", "$merge_dir/peg.map"))
    {
	read_one_pegmap(\*PM, \$peginfo);
	close(PM);
    }
    for my $pm (glob("$genus_dir/*/peg.map"))
    {
	print STDERR "Read $pm\n";
	open(PM, "<", $pm) or die "Could not open $pm: $!";
	read_one_pegmap(\*PM, $peginfo);
	close(PM);
    }
}

sub read_one_pegmap
{
    my($fh, $map) = @_;

    while (<$fh>)
    {
	chomp;
	my @a = split(/\t/);
	$map->{$a[0]} = \@a;
    }
}



1;
