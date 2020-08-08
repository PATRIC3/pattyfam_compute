package P3AlignmentComputeToDisk;

use strict;
use gjoseqlib;
use gjoalignment;
use IPC::Run;
use IO::File;
use IO::Handle;
use Getopt::Long::Descriptive;
use Data::Dumper;
use NCBI_genetic_code;
use Statistics::Descriptive;
use Time::HiRes 'gettimeofday';
use JSON::XS;

use base 'Class::Accessor';

__PACKAGE__->mk_accessors(qw(length_cutoff max_sim dest_dir json));

sub new
{
    my($class, $dest_dir, %params) = @_;

    my $self = {
	length_cutoff => 2,
	max_sim => 0.9,
	dest_dir => $dest_dir,
	json => JSON::XS->new->pretty(1),
	%params,
    };
    return bless $self, $class;
}

#
# Process one block of data, passed as a list of sequence triples.
#
sub process_data
{
    my($self, $genus, $family, $data, $seq_type, $genetic_code) = @_;

    my $is_start_codon = NCBI_genetic_code::is_start_codon($genetic_code);
    my $is_stop_codon = NCBI_genetic_code::is_stop_codon($genetic_code);

    my $stype = uc(substr($seq_type,0,1));
    my $raw_key = join("-", $genus, $family, $stype, 'R');
    my $clean_key = join("-", $genus, $family, $stype, 'C');

    my $fam_dir = $self->dest_dir. "/$family";
    mkdir($fam_dir);

    #
    # Discard entries without valid stop/start.
    # Compute statistics on remaining seqs.
    #
    my $stats = Statistics::Descriptive::Sparse->new;
    my $raw_stats = Statistics::Descriptive::Sparse->new;
    my @ali2;
    for my $ent (@$data)
    {
	my($id, $def, $seq) = @$ent;
	
	my $l = length($seq);
	$raw_stats->add_data($l);
	
	if ($seq_type eq 'dna')
	{
	    my($start) = $seq =~ /^([a-zA-Z]{3})/;
	    my($stop) = $seq =~ /([a-zA-Z]{3})$/;
	    # print "$id $start $stop\n";
	    
	    if (!$is_start_codon->{$start})
	    {
		# print STDERR "Skip $id: invalid start $start\n";
		next;
	    }
	    if (!$is_stop_codon->{$stop})
	    {
		# print STDERR "Skip $id: invalid stop $stop\n";
		next;
	    }
	}
	
	push(@ali2, [@$ent, $l]);
	
	$stats->add_data($l);
    }

    if (@ali2 == 0)
    {
	# print STDERR "$genus $family: No sequences left after bad start/stop trimming. Emitting empty file\n";
	$self->write_stats({ error => 'No sequences left after bad start/stop trimming' }, $seq_type, $fam_dir);
	return;
    }
    
    # print STDERR "mean=" . $stats->mean() . "\n";
    # print STDERR "dev=" . $stats->standard_deviation() . "\n";
    
    #
    # Remove sequences more than 2 devs away from mean.
    #
    my $min = $stats->mean() - $self->length_cutoff * $stats->standard_deviation();
    my $max = $stats->mean() + $self->length_cutoff * $stats->standard_deviation();

    my $trimmed_stats = Statistics::Descriptive::Sparse->new;

    my @ali3;
    for my $ent (@ali2)
    {
	if ($ent->[3] >= $min && $ent->[3] <= $max)
	{
	    push(@ali3, $ent);
	    $trimmed_stats->add_data($ent->[3]);
	}
    }

    if (@ali3 == 0)
    {
	# print STDERR "$genus $family: No sequences left after length trim. Emitting empty file\n";
	$self->write_stats({ error => 'No sequences left after length trim.' }, $seq_type, $fam_dir);
	return;
    }
    
    # printf STDERR "Before trim: %d after frag trim: %d after len trim: %d\n", scalar @$data, scalar @ali2, scalar @ali3;

    #
    # Align.
    #
    # If we just have a single sequence, copy in to out.
    #
    my $mafft_elap;
    my @ali4;
    my @ali5;
    
    if (@ali3 == 1)
    {
	@ali4 = @ali3;
	@ali5 = @ali3;
	$mafft_elap = 0;
    }
    else
    {
	my $to_mafft = IO::Handle->new;
	my($mafft_out, $mafft_err);
	my $t1 = gettimeofday;
	my $h = IPC::Run::start(["mafft", "--reorder", "--anysymbol", "-"],
				"<pipe", $to_mafft,
				">", \$mafft_out,
				"2>", \$mafft_err,
				IPC::Run::timeout(1800),
			       );
	$to_mafft->blocking(1);
	for my $ent (@ali3)
	{
	    my($id, $def, $seq) = @$ent;
	    my $rc = print $to_mafft ">$id $def\n" . pack_seq($seq)  . "\n";
	    if (!$rc)
	    {
		$rc or die "failed $!";
	    }
	}
	close($to_mafft) or die "close failed: $!";

	my $ok;
	eval {
	    $ok = $h->finish();
	};

	if ($@ || !$ok)
	{
	    my $msg = $@ ? $@ : "mafft failed: $mafft_err";
	    print STDERR "Mafft failed: $msg\n";
	    $self->write_stats({ error => "mafft failed: $msg" }, $seq_type, $fam_dir);
	    return;
	}
	
	my $t2 = gettimeofday;
	$mafft_elap = $t2 - $t1;
	
	open(S, "<", \$mafft_out);
	
	@ali4 = read_fasta(\*S);
    
	#
	# Compute representative alignment.
	#
	
	@ali5 = gjoalignment::representative_alignment(\@ali4, { max_sim => $self->max_sim });
    }
    my $raw_ali;
    open(A, ">", \$raw_ali);
    eval {
	write_fasta(\*A, @ali4);
    };
    if ($@)
    {
	die "Error writing ali4: $@\n" . Dumper(\@ali4);
    }
    close(A);

    my $clean_stats = Statistics::Descriptive::Sparse->new;

    open(A, ">", "$fam_dir/clean_$seq_type.fa") or die "Cannot write $fam_dir/clean_$seq_type.fa: $!";
    write_fasta(\*A, \@ali5);
    close(A);

    open(A, ">", "$fam_dir/raw_$seq_type.fa") or die "Cannot write $fam_dir/raw_$seq_type.fa: $!";
    write_fasta(\*A, \@ali4);
    close(A);
    
    for my $ent (@ali5)
    {
	my $len = $ent->[2] =~ tr/-//c;
	$clean_stats->add_data($len);
    }

    open(R, ">", "$fam_dir/reps.$seq_type") or die "Cannot write $fam_dir/reps.$seq_type: $!";
    print R "$_->[0]\n" foreach @ali5;
    close(R);

    my $stat_summary = {
	raw_size => $raw_stats->count(),
	raw_mean => $raw_stats->mean(),
	raw_sd => $raw_stats->standard_deviation(),
	frags_trimmed_size => $stats->count(),
	frags_trimmed_mean => $stats->mean(),
	frags_trimmed_sd => $stats->standard_deviation(),
	clean_size => $clean_stats->count(),
	clean_mean => $clean_stats->mean(),
	clean_sd => $clean_stats->standard_deviation(),
	compute_time => $mafft_elap,
    };
    $self->write_stats($stat_summary, $seq_type, $fam_dir);
    
}

sub write_stats
{
    my($self, $stats, $seq_type, $fam_dir) = @_;
    open(my $fh, ">", "$fam_dir/stats.$seq_type.json") or die "Cannot write $fam_dir/stats.$seq_type.json: $!";
    print $fh $self->json->encode($stats);
    close($fh);
}
1;
