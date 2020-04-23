
=head1 NAME

    pf-find-truncated-genes

=head1 SYNOPSIS

    pf-find-truncated-genes genus-dir work-dir

=head1 DESCRIPTION

This is a pattyfam computation support routine.

Given a genus directory and the associated work directory, compute the potentially
truncated genes (those with stops or starts within 10 base pairs of the end of the contig).

=cut

use strict;
use POSIX;
use File::Temp;
use File::Path qw(make_path);
use File::Copy 'copy';
use File::Basename;
use IPC::Run qw(run);
use Time::HiRes 'gettimeofday';
use File::Spec;
use File::Slurp;
use Getopt::Long::Descriptive;
use List::MoreUtils 'first_index';
use gjoseqlib;
use LPTScheduler;

my($opt, $usage) = describe_options("%c %o genus-dir work-dir",
				    ["parallel|p=i" => "Parallelism", { default => 1 }],
				    ["help|h" => "Show this help message."]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 2;

my $genus_dir = shift;
my $work_dir = shift;

open(G, "<", "$genus_dir/genomes") or die "Cannot open $genus_dir/genomes: $!";

open(TRUNC, ">", "$work_dir/truncated") or die "Cannot write $work_dir/truncated: $!";

my $scheduler = LPTScheduler->new($opt->parallel * 10);

my $out_dir = File::Temp->newdir(CLEANUP => 1);

while (my $genome = <G>)
{
    chomp $genome;

    my $gdir = "/vol/patric3/downloads/genomes/$genome";
    if (! -d $gdir)
    {
	warn "Missing genome dir $gdir\n";
	next;
    }
    my $tab = "$gdir/$genome.PATRIC.features.tab";
    if (! -f $tab)
    {
	warn "Tag file $tab missing\n";
	next;
    }
    
    $scheduler->add_work([$genome, $gdir, $tab], -s $tab);
}

$scheduler->run(sub {
    my $fh;
    if (open($fh, ">", "$out_dir/out.$$"))
    {
	return $fh;
    }
    else
    {
	die "Cannot open $out_dir/out.$$: $!";
    }
}, \&work, $opt->parallel);

opendir(D, $out_dir) or die "Cannot opendir $out_dir: $!";
while (my $f = readdir(D))
{
    next unless $f =~ /^out/;
    copy("$out_dir/$f", \*TRUNC);
}
close(TRUNC);

sub work
{
    my($out_fh, $ent) = @_;
    my($genome, $gdir, $tab) = @$ent;

    if (!open(TBL, "<", $tab))
    {
	die "Cannot open $gdir/$genome.PATRIC.features.tab: $!";
    }
    # find the source records
    my $key = <TBL>;
    chomp $key;
    my @hdrs = split(/\t/, $key);
    my $ftype_idx = first_index { $_ eq 'feature_type' } @hdrs;
    my $id_idx = first_index { $_ eq 'patric_id' } @hdrs;
    my $acc_idx = first_index { $_ eq 'accession' } @hdrs;
    my $start_idx = first_index { $_ eq 'start' } @hdrs;
    my $end_idx = first_index { $_ eq 'end' } @hdrs;
    my $len_idx = first_index { $_ eq 'na_length' } @hdrs;
    my %contig_len;
    while (<TBL>)
    {
	chomp;
	my @vals = split(/\t/);
	if ($vals[$ftype_idx] eq 'source')
	{
	    $contig_len{$vals[$acc_idx]} = $vals[$len_idx];
	}
    }
    seek(TBL, 0, 0);
    my $tried;
    while (<TBL>)
    {
	chomp;
	my @vals = split(/\t/);
	my($id, $type, $start, $end, $acc) = @vals[$id_idx, $ftype_idx, $start_idx, $end_idx, $acc_idx];

	next unless $type eq 'CDS';
	
	my $len = $contig_len{$acc};

	if (!$len && !$tried)
	{
	    $tried = 1;
	    #
	    # Length missing means we need to read the contigs.
	    #
	    if (open(C, "<", "$gdir/$genome.fna") )
	    {
		while (my($id, $def, $seq) = read_next_fasta_seq(\*C))
		{
		    $contig_len{$id} = length($seq);
		}
		close(C);
	    }
	    else
	    {
		warn "Cannot open $gdir/$genome.fna: $!";
	    }
	    $len = $contig_len{$acc};
	    next unless $len;
	}

	if ($start < 10 || $end < 10 || $start > $len - 10 || $end > $len - 10)
	{
	    print $out_fh join("\t", $id, $start, $end, $len), "\n";
	}
    }
}
