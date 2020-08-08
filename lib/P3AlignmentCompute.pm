package P3AlignmentCompute;

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
use DBI;

use base 'Class::Accessor';

__PACKAGE__->mk_accessors(qw(length_cutoff max_sim db_name db_host db_user db_pass dbh));

sub new
{
    my($class, %params) = @_;

    my $dbh = DBI->connect("DBI:mysql:host=$params{db_host};dbname=$params{db_name}",
			   $params{db_user}, $params{db_pass},
		       { RaiseError => 1, AutoCommit => 1 });
    $dbh->do(qq(SET sql_mode='traditional'));
    my $self = {
	length_cutoff => 2,
	max_sim => 0.9,
	dbh => $dbh,
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

    my $dbh = $self->dbh;
    
    my $is_start_codon = NCBI_genetic_code::is_start_codon($genetic_code);
    my $is_stop_codon = NCBI_genetic_code::is_stop_codon($genetic_code);

    my $stype = uc(substr($seq_type,0,1));
    my $raw_key = join("-", $genus, $family, $stype, 'R');
    my $clean_key = join("-", $genus, $family, $stype, 'C');

    $dbh->do(qq(INSERT INTO family_status (genus, family, sequence_type, processing_start, failure_message, complete, failed)
		VALUES (?, ?, ?, CURRENT_TIMESTAMP, NULL, FALSE, FALSE)
		ON DUPLICATE KEY UPDATE processing_start = CURRENT_TIMESTAMP, complete = false, failed = false, processing_end = '0000-00-00 00:00:00', failure_message = NULL),
	     undef, $genus, $family, $stype);

    eval {
	my $str;
	my $fh;
	open($fh, ">", \$str);
	write_fasta($fh, @$data);
	close($fh);
	$dbh->do(qq(INSERT INTO family_sequence VALUES (?,?,?,?)), undef,
		 $genus, $family, $stype, $str);
    };
    if ($@)
    {
	warn "Error writing family seq for $genus $family $seq_type: $@";
    }

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
	print STDERR "$genus $family: No sequences left after bad start/stop trimming. Emitting empty file\n";
	$dbh->do('INSERT INTO alignment (akey, genus, family, clean_flag, sequence_type) VALUES (?,?,?,1,?)',
		 undef, $clean_key, $genus, $family, $stype);
	$dbh->do('INSERT INTO alignment (akey, genus, family, clean_flag, sequence_type) VALUES (?,?,?,0,?)',
		 undef, $raw_key, $genus, $family, $stype);

	$dbh->do(qq(UPDATE family_status
		    SET processing_end = CURRENT_TIMESTAMP,
		        complete = TRUE
		    WHERE genus = ? AND family = ? AND sequence_type = ?), undef,
		 $genus, $family, $stype);
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
	print STDERR "$genus $family: No sequences left after length trim. Emitting empty file\n";
	$dbh->do('INSERT INTO alignment (akey, genus, family, clean_flag, sequence_type) VALUES (?,?,?,1,?)',
		 undef, $clean_key, $genus, $family, $stype);
	$dbh->do('INSERT INTO alignment (akey, genus, family, clean_flag, sequence_type) VALUES (?,?,?,0,?)',
		 undef, $raw_key, $genus, $family, $stype);
	$dbh->do(qq(UPDATE family_status
		    SET processing_end = CURRENT_TIMESTAMP,
		        complete = TRUE
		    WHERE genus = ? AND family = ? AND sequence_type = ?), undef,
		 $genus, $family, $stype);
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
	    $dbh->do(qq(UPDATE family_status
			SET processing_end = CURRENT_TIMESTAMP,
		        complete = TRUE,
			failed = TRUE,
			failure_message = ?
			WHERE genus = ? AND family = ? AND sequence_type = ?), undef,
		     $msg, $genus, $family, $stype);
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

    my $clean_ali;
    open(A, ">", \$clean_ali);
    eval {
	write_fasta(\*A, @ali5);
    };
    if ($@)
    {
	die "Error writing ali5: $@\n" . Dumper(\@ali5);
    }
    close(A);

    my $clean_stats = Statistics::Descriptive::Sparse->new;

    my $sth = $dbh->prepare(qq(INSERT INTO alignment_feature (genus, family, alignment, feature_id) VALUES (?, ?, ?, ?)));
    for my $ent (@ali5)
    {
	my $len = $ent->[2] =~ tr/-//c;
	$clean_stats->add_data($len);
	$sth->execute($genus, $family, $clean_key, $ent->[0]);
    }

    for my $ent (@ali4)
    {
	$sth->execute($genus, $family, $raw_key, $ent->[0]);
    }

    $dbh->do(qq(INSERT INTO family_stats (genus, family, sequence_type,
					  raw_size, raw_mean, raw_sd,
					  frags_trimmed_size, frags_trimmed_mean, frags_trimmed_sd,
					  trimmed_size, trimmed_mean, trimmed_sd,
					  reps_size, reps_mean, reps_sd, compute_time) VALUES (?, ?, ?,
										 ?, ?, ?,
										 ?, ?, ?,
										 ?, ?, ?,
										 ?, ?, ?, ?)),
	     undef,
	     $genus, $family, $stype,
	     $raw_stats->count(), $raw_stats->mean(), $raw_stats->standard_deviation(),
	     $stats->count(), $stats->mean(), $stats->standard_deviation(),
	     $trimmed_stats->count(), $trimmed_stats->mean(), $trimmed_stats->standard_deviation(),
	     $clean_stats->count(), $clean_stats->mean(), $clean_stats->standard_deviation(),
	     $mafft_elap);
    $dbh->do(qq(INSERT INTO alignment (akey, genus, family, clean_flag, sequence_type, alignment_text)
		VALUES (?, ?, ?, 0, ?, ?)), undef,
	     $raw_key, $genus, $family, $stype, $raw_ali);
    $dbh->do(qq(INSERT INTO alignment (akey, genus, family, clean_flag, sequence_type, alignment_text)
		VALUES (?, ?, ?, 1, ?, ?)), undef,
	     $clean_key, $genus, $family, $stype, $clean_ali);

    $dbh->do(qq(UPDATE family_status
		SET processing_end = CURRENT_TIMESTAMP,
		complete = TRUE
		WHERE genus = ? AND family = ? AND sequence_type = ?), undef,
	     $genus, $family, $stype);
    
}

sub clear_genus_data
{
    my($self, $genus) = @_;

    my $dbh = $self->dbh;

    for my $tbl (qw(family_stats alignment alignment_feature family_sequence
		    local_family family_status feature_family family_map))
    {
	print STDERR "Clear $tbl for $genus\n";
	$dbh->do(qq(DELETE FROM $tbl WHERE genus = ?), undef, $genus);
    }
}
   
    
1;

__END__

CREATE TABLE family_stats
(
	genus varchar(255),
	family integer,
	sequence_type char,
	raw_size integer,
	raw_mean float,
	raw_sd float,
	frags_trimmed_size integer,
	frags_trimmed_mean float,
 	frags_trimmed_sd float,
	trimmed_size integer,
	trimmed_mean float,
	trimmed_sd float,
	reps_size integer,
	reps_mean float,
	reps_sd float,
	compute_time float,
	index(genus, family)
) ENGINE = InnoDB;

CREATE TABLE alignment
(
	akey varchar(255) primary key,
	genus varchar(255),
	family integer,
	clean_flag boolean,
	sequence_type char,
	alignment_text longtext,
	index (genus, family)
) ENGINE = InnoDB;

CREATE TABLE alignment_feature
(
	genus varchar(255),
	family integer,
	alignment varchar(255),
	feature_id varchar(255),
	index (genus, family),
	index (alignment),
	index (feature_id)
) ENGINE = InnoDB;

CREATE TABLE family_sequence
(
	genus varchar(255),
	family integer,
	sequence_type char,
	sequence_text longtext,
	PRIMARY KEY (genus, family, sequence_type)  
) ENGINE = InnoDB;

CREATE TABLE local_family
(
	genus varchar(255),
	family integer,
	family_function text,
	family_size integer,
	PRIMARY KEY (genus, family),
	INDEX (family_function(500))
) ENGINE = InnoDB;

CREATE TABLE family_status
(
	genus varchar(255),
	family integer,
	sequence_type char,
	processing_start timestamp DEFAULT 0,
	processing_end timestamp DEFAULT 0,
	failure_message text,
	complete boolean,
	failed boolean,
	PRIMARY KEY (genus, family, sequence_type)
) ENGINE = InnoDB;

CREATE TABLE feature_family
(
	feature_id varchar(255) PRIMARY KEY,
	pgfam_id varchar(12),
	plfam_id integer,
	family_function text,
	genus varchar(255),
	INDEX (family_function(100)),
	INDEX (pgfam_id),
	INDEX (plfam_id),
	INDEX (genus)
) ENGINE = InnoDB;	

CREATE TABLE family_map
(
	pub_pgfam_id varchar(12),
	denovo_pgfam_id varchar(12),
	genus varchar(255),
	pub_plfam_id integer,
	denovo_plfam_id integer,
	INDEX (pub_pgfam_id)
) ENGINE = InnoDB;

CREATE TABLE published_fam_to_alignment
(
	pub_pgfam_id varchar(12),
	alignment_key varchar(255),
	INDEX (pub_pgfam_id)
) ENGINE = InnoDB;

CREATE TABLE published_fam_to_aa_alignment
(
	pub_pgfam_id varchar(12),
	alignment_key varchar(255),
	INDEX (pub_pgfam_id)
) ENGINE = InnoDB;

