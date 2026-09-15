#
# We have truncated seq_lens.txt files in the genome cache.
# Reprocess the aa.fa files to correct them
#

use strict;
use gjoseqlib;

@ARGV == 1 or @ARGV == 2 or die "Usage: $0 cache-dir\n";

my $cache_dir = shift;
my $glist = shift;

if ($glist)
{
    open(G, "<", $glist) or die "cannot open $glist: $!";

    while (my $g = <G>)
    {
	chomp $g;

	my($suf) = $g =~ /(\d\d\d)\./;

	my $dir = "$cache_dir/$suf/$g";
	if (! -d $dir)
	{
	    warn "No dir $dir\n";
	    next;
	}
	fixup($_);
    }
}
else
{
    opendir(D, $cache_dir) or die "Cannot opendir $cache_dir: $!";
    
    for my $top (sort readdir(D))
    {
	next unless $top =~ /^\d+$/;
	
	opendir(G, "$cache_dir/$top") or die "Cannot opendir $cache_dir/$top: $!";
	while (my $g = readdir(G))
	{
	    next unless $g =~ /\d/;
	    my $dir = "$cache_dir/$top/$g";
	    next unless -d $dir;
	    # print "$dir\n";
	    
	    fixup($dir);
	}
    }
}

sub fixup
{
    my($dir) = @_;
    if (!open(LEN, "<", "$dir/seq_lens.txt"))
    {
	if (-f "$dir/aa.fa")
	{
	    die "Cannot open $dir/aa.fa: $!";
	}
	else
	{
	    return;
	}
    }
    local $/;
    undef $/;
    seek(LEN, 2, -10);
    my $txt = <LEN>;
    next if (substr($txt, -1, 1) eq "\n");
    print "Need to fix $dir\n";
    if (!open(FA, "<", "$dir/aa.fa"))
    {
	warn "Cannot open $dir/aa.fa: $!";
	next;
    }
    if (! -s FA)
    {
	system("ls", "-l", "$dir");
	next;
    }
    open(NLEN, ">", "$dir/seq_lens.fixed.txt") or die "Cannot write $dir/seq_lens.fixed.txt: $!";
    while (my($id, $def, $seq) = read_next_fasta_seq(\*AA))
    {
	print NLEN "$id\t" . length($seq) . "\n";
    }
    close(LNEN);
    print "Wrote $dir\n";
    exit;
}
