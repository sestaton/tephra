package Tephra::Genome::SoloLTRSearch;

use 5.014;
use Moose;
use File::Spec;
use File::Find;
use File::Basename;
use File::Copy;
use File::Temp          qw(tempfile);
use File::Path          qw(make_path remove_tree);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use Time::HiRes         qw(gettimeofday);
use List::UtilsBy       qw(nsort_by);
use Cwd                 qw(abs_path);
use Sort::Naturally;
use Path::Class::File;
use Bio::DB::HTS::Kseq;
use Bio::AlignIO;
use Bio::SearchIO;
use Sort::Naturally;
use Set::IntervalTree;
use Carp 'croak';
use Tephra::Config::Exe;
use namespace::autoclean;
#use Data::Dump::Color;

with 'Tephra::Role::File',
     'Tephra::Role::Run::Any',
     'Tephra::LTR::Role::Utils';

=head1 NAME

Tephra::Genome::SoloLTRSearch - Find solo-LTRs in a refence genome

=head1 VERSION

Version 0.13.1

=cut

our $VERSION = '0.13.1';
$VERSION = eval $VERSION;

has dir => (
      is       => 'ro',
      isa      => 'Path::Class::Dir',
      required => 1,
      coerce   => 1,
);

has genome => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has report => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 0,
    coerce   => 1,
);

has outfile => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has seqfile => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 0,
    coerce   => 1,
);

has percentid => (
    is        => 'ro',
    isa       => 'Num',
    predicate => 'has_percentid',
    lazy      => 1,
    default   => 39,
);

has percentcov => (
    is        => 'ro',
    isa       => 'Num',
    predicate => 'has_percentcov',
    lazy      => 1,
    default   => 80,
);

has matchlen => (
    is        => 'ro',
    isa       => 'Num',
    predicate => 'has_matchlen',
    lazy      => 1,
    default   => 80,
);

has clean => (
    is       => 'ro',
    isa      => 'Bool',
    required => 0,
    default  => 1,
);

has fullanalysis => (
    is       => 'ro',
    isa      => 'Bool',
    required => 0,
    default  => 0,
);

has numfamilies => (
    is        => 'ro',
    isa       => 'Int',
    predicate => 'has_numfamilies',
    lazy      => 1,
    default   => 20,
);

has threads => (
    is        => 'ro',
    isa       => 'Int',
    predicate => 'has_threads',
    lazy      => 1,
    default   => 1,
);

has debug => (
    is         => 'ro',
    isa        => 'Bool',
    predicate  => 'has_debug',
    lazy       => 1,
    default    => 0,
);

#
# Methods
#
sub find_soloLTRs {
    my $self = shift;
    my $anno_dir = $self->dir->absolute->resolve;
    my $genome   = $self->genome->absolute->resolve;
    my $report   = $self->report; 
    my $seqfile  = $self->seqfile;

    my $write_seqs = $seqfile ? 1 : 0;

    my @sfs;
    find( sub { push @sfs, $File::Find::name if -d && /_copia\z|_gypsy\z|_unclassified\z/ }, $anno_dir);
    croak "\n[ERROR]: Could not find the expected sub-directories ending in 'copia', 'gypsy' and 'unclassified'. ".
	"Please check input. Exiting.\n" unless @sfs; #== 2;

    my $forks = @sfs;
    my $thr = $forks > $self->threads ? $self->threads : $forks;

    my $pm = Parallel::ForkManager->new($thr);
    local $SIG{INT} = sub {
        warn "Caught SIGINT; Waiting for child processes to finish.";
        $pm->wait_all_children;
        exit 1;
    };

    my (@soloLTR_reps, @soloLTR_seqs, @model_dirs);
    $pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
			      if (defined $data_ref) { 
				  my ($model_dir, $reports, $seqfiles) = @{$data_ref}{qw(model_dir reports seqfiles)};
				  if (defined $model_dir && defined $reports) { 
				      push @model_dirs, $model_dir;
				      push @soloLTR_reps, @$reports;
				  }
				  push @soloLTR_seqs, @$seqfiles 
				      if defined $seqfiles && $write_seqs;
			      }
			} );

    my $dirct = 0;
    for my $dir (nsort @sfs) {
	$dirct++;
	my $sf = (split /_/, $dir)[-1];
	$pm->start($sf) and next;
	$SIG{INT} = sub { $pm->finish };
	my $data_ref = $self->run_parallel_model_search($genome, $dir, $write_seqs, $thr);
	my ($model_dir, $reports, $seqfiles) = @{$data_ref}{qw(model_dir reports seqfiles)};
	$pm->finish(0, { model_dir => $model_dir, reports => $reports, seqfiles => $seqfiles });
    }
    $pm->wait_all_children;

    my $soloct = $self->collate_sololtr_reports(\@soloLTR_reps, \@soloLTR_seqs, $report, $write_seqs);

    if ($soloct > 0) {
	$self->write_sololtr_gff($report);
    }
    else {
	say STDERR "\n[WARNING]: No solo-LTRs were found so none will be reported.\n";
	unlink $report;
	unlink $seqfile if -e $seqfile;
    }

    if ($self->clean) {
	for my $dir (@model_dirs) { 
	    remove_tree( $dir, { safe => 1} );                                                                
	}
    }

    return;
}

sub run_parallel_model_search {
    my $self = shift;
    my ($genome, $dir, $write_seqs, $thread_level) = @_;
    my $numfams = $self->numfamilies;
    my $allfams = $self->fullanalysis;
    my $threads = $self->threads;
    my $outfile = $self->report;

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($hmm3bin) = @{$config}{qw(hmmer3bin)};
    my $nhmmer = File::Spec->catfile($hmm3bin, 'nhmmer');

    my $model_dir = File::Spec->catdir($dir, 'Tephra_LTR_exemplar_models');
    unless ( -e $model_dir ) {
        make_path( $model_dir, {verbose => 0, mode => 0771,} );
    }

    my $ltrfams = $self->get_exemplar_ltrs_for_sololtrs({ input_dir => $dir, full_analysis => $allfams });
    #dd $ltrfams and exit;
    return unless defined $ltrfams && %$ltrfams;
    
    my (@parsed_seqfiles, @parsed_reports);
    my %reports;
    my $t0 = gettimeofday();

    my ($gname, $gpath, $gsuffix) = fileparse($genome, qr/\.[^.]*/);
    my $logfile = File::Spec->catfile($model_dir, 'all_solo-ltr_searches.log');
    open my $log, '>>', $logfile or die "\n[ERROR]: Could not open file: $logfile\n";

    my $thr;
    if ($threads > 1) { 
	if ($threads % 3 == 0) {
	    $thr = sprintf("%.0f",$threads/3);
	}
	elsif (+($threads-1) % 3 == 0) {
	    $thr = sprintf("%.0f",($threads-1)/3);
	}
	elsif (+($threads-2) % 3 == 0) {
	    $thr = sprintf("%.0f",($threads-2)/3);
	}
	else {
	    $thr = 1;
	}
    }
    else {
	$thr = 1;
    }

    my $pm = Parallel::ForkManager->new($thr, $gpath);
    local $SIG{INT} = sub {
        warn "Caught SIGINT; Waiting for child processes to finish.";
        $pm->wait_all_children;
        exit 1;
    };

    $pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
			      my ($query, $nhmmer_out, $aln_stats) = @{$data_ref}{qw(query report aln_stats)};
                              if (defined $nhmmer_out && -e $nhmmer_out && -s $nhmmer_out) {
                                  my ($parsed_report, $parsed_seqs) = 
				      $self->write_nhmmer_report($query, $aln_stats, $nhmmer_out, $write_seqs);
				  push @parsed_reports, $parsed_report
				      if defined $parsed_report;
				  push @parsed_seqfiles, $parsed_seqs 
				      if defined $parsed_report && $write_seqs;
			      }
                              #unlink $_ for keys %$data_ref;
                              my $t1 = gettimeofday();
                              my $elapsed = $t1 - $t0;
                              my $time = sprintf("%.2f",$elapsed/60);
                              say $log basename($ident),
			      " just finished with PID $pid and exit code: $exit_code in $time minutes";
                        } );


    my $fam_ct = 0;
    for my $element (nsort keys %$ltrfams) {
	$fam_ct++;
        $pm->start($element) and next;
        $SIG{INT} = sub { $pm->finish };
	my ($query, $nhmmer_out, $aln_stats) = 
	    $self->run_serial_model_search($nhmmer, $gname, $genome, $element, $ltrfams, $model_dir, $fam_ct);
	%reports = ( query => $query, report => $nhmmer_out, aln_stats => $aln_stats );

        $pm->finish(0, \%reports);
    }

    $pm->wait_all_children;

    my $t2 = gettimeofday();
    my $total_elapsed = $t2 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);

    say $log "\n========> Finished HMMER search on $fam_ct solo-LTR models in $final_time minutes";
    close $log;

    # Each process stores the child's data structure and I've found this can cause issues with the
    # default /tmp directory in some cases, so we store this data in the working directory and clean up the
    # files manually. 
    my @process_logs;
    find( sub { push @process_logs, $File::Find::name if -f and /Parallel-ForkManager*.txt/ }, $gpath );
    unlink @process_logs;

    return { reports => \@parsed_reports, seqfiles => \@parsed_seqfiles, model_dir => $model_dir };
}

sub run_serial_model_search {
    my $self = shift;
    my ($nhmmer, $gname, $genome, $element, $ltrfams, $model_dir, $fam_ct) = @_;
    my $numfams = $self->numfamilies;
    my $allfams = $self->fullanalysis;  

    if ($allfams || (!$allfams && $numfams >= $fam_ct)) {
	my $ltrfile = File::Spec->catfile($model_dir, $element.'_ltrseqs.fasta');
	open my $ltrfh, '>', $ltrfile or die "\n[ERROR]: Could not open file: $ltrfile\n";
	
	for my $pair (@{$ltrfams->{$element}}) {
	    say $ltrfh join "\n", ">".$pair->{id}, $pair->{seq};
	}
	close $ltrfh;
	
	my ($name, $path, $suffix) = fileparse($ltrfile, qr/\.[^.]*/);
	my $tre = File::Spec->catfile( abs_path($path), $name.'.dnd' );
	my $aln = File::Spec->catfile( abs_path($path), $name.'_muscle-out.aln' );
	my $log = File::Spec->catfile( abs_path($path), $name.'_muscle-out.log' );

	my $config  = Tephra::Config::Exe->new->get_config_paths;
	my ($muscle) = @{$config}{qw(muscle)};

	my $muscmd = "$muscle -quiet -clwstrict -in $ltrfile -out $aln -log $log";
	say STDERR "DEBUG: $muscmd" if $self->debug;
	my $status;
	if ($allfams) {
	    $status = $self->capture_cmd($muscmd);
	}
	else {
	    if ($numfams >= $fam_ct) {
		$status = $self->capture_cmd($muscmd);
	    }
	}
	#unlink $log, $ltrfile;

	if (defined $status && $status =~ /failed/i) {
	    say STDERR "\n[ERROR]: $muscmd failed. Removing $ltrfile";
	    #unlink $ltrfile;
	    return undef;
	}
	
	if (defined $aln && -e $aln && -s $aln) {
	    my $query = $aln =~ s/\.aln$//r;
	    $query = basename($query);
	    #my $seq_stats = $self->get_seq_len($ltrfile);
	    #my $aln_stats;# = $self->get_aln_len($aln);
	    my ($stdout, $stderr, $aln_stats) = capture { 
		my $ltr_retro = basename($aln);
		$ltr_retro =~ s/\.aln//;
		my $aln_in = Bio::AlignIO->new(-file => $aln, -format => 'clustalw');
		$aln_in->alphabet('dna');

		my $aln_obj = $aln_in->next_aln;
		#$aln_stats{$ltr_retro} = $aln_obj->length;
		return { $ltr_retro => $aln_obj->length };
                #return $self->get_aln_len($aln); }; 
	    };
	    my ($len) = values %$aln_stats;
	    return undef unless $len =~ /\d+/;
	    my $nhmmer_out = $self->search_with_models($gname, $aln, $nhmmer, $genome);

	    return ($query, $nhmmer_out, $aln_stats);
	}
	else {
	    return undef;
	}
    }
    else {
	return undef;
    }
}

sub search_with_models {
    my $self = shift;
    my ($gname, $aln, $nhmmer, $genome) = @_;

    my $nhmmer_out = $aln;
    $nhmmer_out =~ s/\.aln.*$//;
    $nhmmer_out .= '_'.$gname.'.nhmmer';

    my $nhmmer_cmd = "$nhmmer --cpu 1 -o $nhmmer_out $aln $genome";
    say STDERR "DEBUG: $nhmmer_cmd" if $self->debug;
    $self->run_cmd($nhmmer_cmd);

    return $nhmmer_out;
}

sub write_nhmmer_report {
    my $self = shift;
    my ($query, $aln_stats, $search_report, $write_seqs) = @_;
    my $match_pid  = $self->percentid;
    my $match_len  = $self->matchlen;
    
    my $parsed_report = $search_report =~ s/\.nhmmer/_parsed.txt/r;
    my $parsed_seqs = $search_report =~ s/\.nhmmer/_seqs.fasta/r;
    
    my $num_hits = 0;
    my (%results, %ids);
    open my $fh, '<', $search_report or die "\n[ERROR]: Could not open file: $search_report\n";;

    while (my $line = <$fh>) {
	chomp $line;
	$line =~ s/^\s+|\s+$//;
	next if $line =~ /^#|^[-]+|E-value|Scores/;
	last if $line =~ /Annotation for each hit/i;
	last if $line =~ /No hits detected/i;
	my $q;
	#Query:       RLC_family5_LTR_retrotransposon595_ltrseqs_muscle-out  [M=158]
	#Scores for complete hits:
	#  E-value  score  bias  Sequence    start      end  Description
	#  ------- ------ -----  --------    -----    -----  -----------
	#  0.00081   25.6   6.6  Chr5     20120128 20120038  
	if ($line =~ /^Query:\s+(\S+)\s+\S+/) { 
	    $q = $1;
	    $ids{query}{$q} = 1;
	}
	else { 
	    my @f = split /\s+/, $line;
	    if (@f == 6) {
		$num_hits++;
		$ids{subject}{$f[3]} = 1;
	    }
	}
    }
    close $fh;

    if ($num_hits == 0) {
	say STDERR "\n[WARNING]: No hits found for: $search_report.\n" if $self->debug;
	unlink $search_report;
	return undef;
    }

    open my $outfh, '>', $parsed_report or die "\n[ERROR]: Could not open file: $parsed_report\n";
    open my $seqfh, '>', $parsed_seqs or die "\n[ERROR]: Could not open file: $parsed_seqs\n"
        if $write_seqs;

    {
	open my $in, '<', $search_report or die "\n[ERROR]: Could not open file: $search_report\n";
	local $/ = "\n>>";
	
	while (my $line = <$in>) {
	    chomp $line;
	    next if $line =~ /^#/;
	    my @line = split /\n/, $line;
	    my ($qid, $qstring, $sstring, $sid, $len, $previd, $score, $bias, $evale, $hmmfrom, $hmmto, 
		$alifrom, $alito, $envfrom, $envto, $seqlen, $acc, $seq);
	    for my $l (@line) {
		$l =~ s/^\s+//;
		if ($l =~ /(\S+)\s+?.*/) {
		    my $id = $1;
		    $sid = $id if exists $ids{subject}{$id};
		}
		next if $l =~ /^score|^------\s+|^alignment/i;
		if ($l =~ /^(?:\!|\?)\s+\d+/) {
		    ##score  bias    Evalue   hmmfrom    hmm to     alifrom    ali to      envfrom    env to       sq len      acc 
		    my @info = split /\s+/, $l;
		    ($score, $bias, $evale, $hmmfrom, $hmmto, $alifrom, $alito, $envfrom, $envto, $seqlen, $acc) = 
			@info[1..5,7..8,10..11,13..14];
		}
		if ($l =~ /^(\S+)\s+\d+\s+([atcgnATCGNRYSWKMBDHVU.*-]+)\s+\d+/) {
		    my $id = $1;
		    my $seqstr = $2;
		    if ($id eq $sid) {
			$sstring .= $seqstr;
		    }
		    if ($id eq $query && exists $ids{query}{$id}) {
			$qid = $id;
			$qstring .= $seqstr;
		    }
		    
		}
		if ($l =~ /[atcgnATCGN.*-]+/) {
		    next if $l =~ /\d+/;
		    $seq .= $l;
		}
	    }
	    ## for DEBUG 
	    #unless (defined $seq) {
	        #say STDERR "ERROR: seq not defined";
		#dd \@line and exit;
	    #}
	    #unless (defined $sstring) {
		#say STDERR "ERROR: sstring not defined";
		#dd \@line and exit;
	    #}
	    #unless (defined $qid) {
		#say STDERR "query: $query";
		#dd \@line and exit;
	    #}
	    ## 
	    my $identical = ($seq =~ tr/a-zA-Z//);
	    my $qlen = $aln_stats->{$qid};
	    my $pid = sprintf("%.2f", $identical/$qlen * 100);	
	    my $hsplen = $alito-$alifrom+1;
	    my $percent_q_coverage = sprintf("%.2f", $hsplen/$qlen * 100);
	    
	    if ($hsplen >= $match_len && $pid >= $match_pid) {    
		my $key = join "||", $qid, $alifrom;
		push @{$results{$key}}, 
		    join "||", $qid, $qlen, $sid, $percent_q_coverage, $hsplen, $pid, $hmmfrom, $hmmto, $alifrom, $alito, $qstring, $sstring;
	    }
	    undef $seq;
	    undef $qstring;
	    undef $sstring;
	}
	close $in;
    }

    my $solo_ct = 0;
    for my $qid (nsort_by { m/\S+\|\|(\d+)/ and $1 } keys %results) {
	for my $res (@{$results{$qid}}) {
	    $solo_ct++;
	    my ($qid, $qlen, $sid, $percent_q_coverage, $hsplen, $pid, $hmmfrom, $hmmto, $alifrom, $alito, $qstring, $sstring) = 
		split /\|\|/, $res;
	    $qid =~ s/_ltrseqs_muscle-out//;
	    say $outfh join "\t", $qid, $qlen, $num_hits, $sid, $percent_q_coverage, 
	        $hsplen, $pid, $hmmfrom, $hmmto, $alifrom, $alito;
	    
	    if ($write_seqs) {
		my $seqid = join '_', '>'.$qid, $sid, $alifrom, $alito;
		## We show the location of the solo-LTR in the header, as in the GFF3.
		$qstring = uc($qstring);
		$qstring =~ s/\./N/g;
		say $seqfh join "\n", $seqid, $qstring;
	    }

	}
    }
    close $outfh;
    close $seqfh if $write_seqs;
    #unlink $search_report;

    return ($parsed_report, $parsed_seqs);
}

sub write_sololtr_gff {
    my $self = shift;
    my ($report) = @_;
    my $genome  = $self->genome->absolute->resolve;
    my $outfile = $self->outfile;;

    my $seqlen = $self->get_seq_len($genome);
    open my $in, '<', $report or die "\n[ERROR]: Could not open file: $report\n";
    open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";

    if (! -s $outfile) {
	say $out '##gff-version 3';
	for my $id (nsort keys %$seqlen) {
	    say $out join q{ }, '##sequence-region', $id, '1', $seqlen->{$id};
	}
    }

    my $ct = 0;
    my (%features, %intervals, %ltr_lengths);
    while (my $line = <$in>) {
	chomp $line;
	next if $line =~ /^#/;
	next unless $line =~ /\S/;
	my ($model_num, $query, $query_length, $number_of_hits, $hit_name, $perc_coverage, $hsp_length, 
	    $hsp_perc_ident, $hsp_query_start, $hsp_query_end, $hsp_hit_start, $hsp_hit_end) = split /\t/, $line;

	if (exists $seqlen->{$hit_name}) {
	    $ct++;
	    say $out join "\t", $hit_name, 'Tephra', 'solo_LTR', $hsp_hit_start, $hsp_hit_end, 
	        '.', '?', '.', "ID=$model_num;match_id=$query;Name=solo_LTR;Ontology_term=SO:0001003";
	}
	else {
	    croak "\n[ERROR]: $hit_name not found in $genome. This is a bug, please report it. Exiting.\n";
	}
    }
    close $in;
    close $out; 

    return;
}

sub collate_sololtr_reports {
    my $self = shift;
    my ($soloLTR_reps, $soloLTR_seqs, $report, $write_seqs) = @_;

    my (%seen, %parsed_alns);
    for my $file (grep { -s $_ } @$soloLTR_reps) {
	open my $fh_in, '<', $file or die "\n[ERROR]: Could not open file: $file\n";
	my $first = <$fh_in>;
	my ($qid, $qlen, $num_hits, $sid, $percent_q_coverage,
	    $hsplen, $pid, $hmmfrom, $hmmto, $alifrom, $alito) = split /\t/, $first;
	my ($last_chr, $last_start, $last_end) = ($sid, $alifrom, $alito);

	while (my $line = <$fh_in>) {
	    chomp $line;
	    next if $line =~ /^#/;
	    next unless $line =~ /\S/;
	    ($qid, $qlen, $num_hits, $sid, $percent_q_coverage,
	     $hsplen, $pid, $hmmfrom, $hmmto, $alifrom, $alito) = split /\t/, $line;
	    if (exists $seen{$sid}{$alifrom}) {
		my ($prev_len, $prev_pid) = split /\|\|/, $seen{$sid}{$alifrom};
		if ($alito-$alifrom+1 > $prev_len && $pid > $prev_pid) {
		    $parsed_alns{$sid}{$alifrom} = join "||", ($qid, $qlen, $num_hits, $sid, $percent_q_coverage,
							       $hsplen, $pid, $hmmfrom, $hmmto, $alifrom, $alito);
		}
	    }
	    if ($last_chr eq $sid) { 
		if ($alifrom > $last_start && $alifrom > $last_end && $alito > $last_end) {
		    $parsed_alns{$sid}{$alifrom} = join "||", ($qid, $qlen, $num_hits, $sid, $percent_q_coverage,
							       $hsplen, $pid, $hmmfrom, $hmmto, $alifrom, $alito);
		    $seen{$sid}{$alifrom} = join "||", $alito-$alifrom+1, $pid;
		}
	    }
	    else {
		# solo-LTRs on a new chromosome
		$parsed_alns{$sid}{$alifrom} = join "||", ($qid, $qlen, $num_hits, $sid, $percent_q_coverage,
							   $hsplen, $pid, $hmmfrom, $hmmto, $alifrom, $alito);
		$seen{$sid}{$alifrom} = join "||", $alito-$alifrom+1, $pid;
	    }
	    ($last_chr, $last_start, $last_end) = ($sid, $alifrom, $alito);
	}
	close $fh_in;

	# if there was only one line in the report, store that now
	unless (%parsed_alns) {
	    $parsed_alns{$sid}{$alifrom} = join "||", ($qid, $qlen, $num_hits, $sid, $percent_q_coverage,
						       $hsplen, $pid, $hmmfrom, $hmmto, $alifrom, $alito);
	}
    }

    open my $repfh, '>>', $report or die "\n[ERROR]: Could not open file: $report\n";
    say $repfh join "\t", "#model_num", "query", "query_length", "number_of_hits", "hit_name",
        "perc_coverage","hsp_length", "hsp_perc_ident","hsp_query_start", "hsp_query_end", "hsp_hit_start", "hsp_hit_end";

    my $solo_ct = 0;
    my %id_map;
    for my $chr (nsort keys %parsed_alns) {
	for my $start (sort { $a <=> $b } keys %{$parsed_alns{$chr}}) {
	    $solo_ct++;
	    my @f = split /\|\|/, $parsed_alns{$chr}{$start};
	    chomp $f[10];
	    say $repfh join "\t", "solo_LTR$solo_ct", @f;
	    $id_map{ $f[0].'_'.$chr.'_'.$start.'_'.$f[10] } = "solo_LTR$solo_ct".'_'.$f[0].'_'.$chr.'_'.$start.'_'.$f[10];
	}
    }
    #dd \%id_map;
    #dd $soloLTR_seqs;

    if ($write_seqs && @$soloLTR_seqs) { 
        my $seqfile = $self->seqfile;
        open my $seqfh, '>', $seqfile or die "\n[ERROR]: Could not open file: $seqfile";

	for my $file (@$soloLTR_seqs) {
	    my $kseq = Bio::DB::HTS::Kseq->new($file);
	    my $iter = $kseq->iterator;
	    while (my $seqo = $iter->next_seq) {
		my $id = $seqo->name;
		my $seq = $seqo->seq;
		$seq =~ s/.{60}\K/\n/g;
		if (exists $id_map{$id}) { 
		    # Do not warn if the IDs were not found; that means they were overlapping and removed.
		    say $seqfh join "\n", '>'.$id_map{$id}, $seq;
		}
	    }
	}
    }

    return $solo_ct;
}

sub check_report_summary {
    my $self = shift;
    my ($hmmsearch_summary) = @_;

    open my $in, '<', $hmmsearch_summary or die "\n[ERROR]: Could not open file: $hmmsearch_summary\n";
    my $ct = 0;

    while (my $line = <$in>) {
	chomp $line;
	next if $line =~ /^#/;
	$ct++;
    }

    return $ct;
}
    
sub get_seq_len {
    my $self = shift;
    my ($genome) = @_;
    
    my %len;

    my $kseq = Bio::DB::HTS::Kseq->new($genome);
    my $iter = $kseq->iterator();

    while ( my $seq = $iter->next_seq() ) {
	my $id  = $seq->name;
	my $seq = $seq->seq;
	$len{$id} = length($seq);
    }       

    return \%len;
}

sub get_aln_len {
    my $self = shift;
    my ($aln) = @_;
    
    my %aln_stats;

    my $ltr_retro = basename($aln);
    $ltr_retro =~ s/\.aln//;
    my $aln_in = Bio::AlignIO->new(-file => $aln, -format => 'clustalw');
    $aln_in->alphabet('dna');

    while ( my $aln_obj = $aln_in->next_aln() ) {
	$aln_stats{$ltr_retro} = $aln_obj->length;
    }	

    return \%aln_stats;
}

=head1 AUTHOR

S. Evan Staton, C<< <evan at evanstaton.com> >>

=head1 BUGS

Please report any bugs or feature requests through the project site at 
L<https://github.com/sestaton/tephra/issues>. I will be notified,
and there will be a record of the issue. Alternatively, I can also be 
reached at the email address listed above to resolve any questions.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Tephra::Genome::SoloLTRSearch


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
