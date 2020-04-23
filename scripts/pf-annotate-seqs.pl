
=head1 NAME

    pf-annotate-seqs

=head1 SYNOPSIS

    pf-annotate-seqs kmer-server-url genus-data-dir sequences-dir calls-file uncalled-ids-file

=head1 DESCRIPTION

This is a pattyfam computation support routine.

Annotate the sequences in the given sequences directory using the given kmers. Write
Write the calls to calls-file and the list of IDs for which calls could not be
made to uncalled-ids-file.
    
Note that the URL we'll be passed to do the processing may be the one that accumulates
kmer hits for doing later kmer similarity processing. The kser code that handles the
add-hits processing appears to be multithread safe so we can run in parallel here.

=cut

use strict;
use Getopt::Long::Descriptive;
use File::Temp;
use File::Copy;
use Data::Dumper;
use LPTScheduler;

my($opt, $usage) = describe_options("%c %o kmer-server-url genus-data-dir sequences-dir calls-file uncalled-ids-file",
				    ["parallel=s" => "Parallel threads", { default => 1 }],
				    ["help|h" => "Show this help message."]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 5;

my $url = shift;
my $fam_dir = shift;
my $seqs_dir = shift;
my $calls_file = shift;
my $uncalled_file = shift;

#
# Construct the LPT scheduler where work elements are genomes to be annotated. Work size is
# the size of the file. Since we are splitting the work up, we will
# need a temp dir into which each chunk will write its output. We'll merge them
# when all are done.
#
my $anno_temp = File::Temp->newdir(CLEANUP => 1);
my $scheduler = LPTScheduler->new(10 * $opt->parallel);

open(GENOMES, "<", "$fam_dir/genomes") or die "Cannot open $fam_dir/genomes: $!";

while (my $genome = <GENOMES>)
{
    chomp $genome;

    my $seqs = "$seqs_dir/$genome";

    if (! -f $seqs)
    {
	warn "Missing sequences file $seqs";
    }
    else
    {
	my $size = -s $seqs;
	if ($size)
	{
	    # print "Work: $genome $size\n";
	    $scheduler->add_work([$genome, "$seqs_dir/$genome", $anno_temp], $size);
	}
    }
}

open(CALLS, ">", $calls_file) or die "Cannot write $calls_file: $!";
open(UNCALLED, ">", $uncalled_file) or die "Cannot write $uncalled_file: $!";

#
# Run. We use the global callback to open up a per-worker file.
#

$scheduler->run(sub {
    my($called, $uncalled);
    open($called, ">", "$anno_temp/called.$$") or die "cannot write $anno_temp/called.$$: $!";
    open($uncalled, ">", "$anno_temp/uncalled.$$") or die "cannot write $anno_temp/called.$$: $!";
    return [$called, $uncalled];
},
		\&work, $opt->parallel);

#
# Collect output.
#

opendir(D, $anno_temp) or die "Cannot open $anno_temp: $!";
while (my $p = readdir(D))
{
    if ($p =~ /^called/)
    {
	copy("$anno_temp/$p", \*CALLS);
    }
    elsif ($p =~ /^uncalled/)
    {
	copy("$anno_temp/$p", \*UNCALLED);
    }
}
closedir(D);
close(CALLS);
close(UNCALLED);
    
sub work
{
    my($global, $ent) = @_;
    my($called, $uncalled) = @$global;

    my($genome, $seqs, $anno_temp) = @$ent;

    #
    # Annotate the sequences using a curl to the kmer server URL.
    # We scan the file first
    # so we can tell what sequences were not annotated.
    #
    my $tmp = File::Temp->new(UNLINK => 1);
    close($tmp);
    my @cmd = ("curl", "-s", "--data-binary", '@-', "-o", "$tmp", $url);
    # my @cmd = ('kmer_search', @url, "-d", $kmer_dir, "-a", "-o", "$tmp");
    
    open(SEQS, "<", $seqs) or die "Cannot open $seqs: $!";
    open(KS, "|-", @cmd) or die "Cannot open @cmd: $!";
    my %seen;
    while (<SEQS>)
    {
	if (/^>(\S+)/)
	{
	    $seen{$1} = 1;
	}
	print KS $_;
    }
    close(SEQS);
    close(KS) or die "Error on close of @cmd: $! $?";
    
    #
    # Read generated calls to enable logging of uncalled IDs
    #
    # The raw kmer server output will include all the low level
    # CALL / OTU-COUNTS / etc output, but also a BEST-CALL
    # line that incorporates the logic originally from kmer_search.pl.
    #
    # We use this here to remove the process-call overhead that
    # we have with a lot of small files.
    #
    
    open(TMP, "<", "$tmp") or die "Cannot read $tmp: $!";
    while (<TMP>)
    {
	# print $called $_;
	if (my($id, $fn, $score) = /^BEST-CALL\t([^\t]+)\t([^\t]+)\t([^\t]+)\t/)
	{
	    #
	    # Score will be zeroish if no kmer call. Worry about FP precision and comparisons to 0.0 so just
	    # mark uncalled if <1.
	    if ($score < 1)
	    {
		print $uncalled "$id\t$fn\n";
	    }
	    else
	    {
		print $called "$id\t$fn\n";
	    }
	    
	    my $x = delete $seen{$id};
	    if (!$x)
	    {
		warn "Delete '$id' from seen fails\n";
	    }
	}
    }
    close(TMP);
    
    for my $id (keys %seen)
    {
	print $uncalled "$id\n";
    }
}


