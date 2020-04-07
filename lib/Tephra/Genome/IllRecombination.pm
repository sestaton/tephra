package Tephra::Genome::IllRecombination;

use 5.014;
use Moose;
use File::Spec;
use File::Find;
use File::Basename;
use File::Path   qw(make_path remove_tree);
use File::Temp   qw(tempfile);
use Time::HiRes  qw(gettimeofday);
use Log::Any     qw($log);
use Scalar::Util qw(openhandle);
use Cwd          qw(getcwd abs_path);
use Path::Class::File;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::HTS::Kseq;
use Bio::AlignIO;
use Bio::SearchIO;
use Sort::Naturally;
use Parallel::ForkManager;
use Statistics::Descriptive;
use Tephra::Config::Exe;
#use Data::Dump::Color;
use namespace::autoclean;

with 'Tephra::Role::File',
     'Tephra::Role::Run::Any';

=head1 NAME

Tephra::Genome::IllRecombination - Calculate illegitimate recombination in a genome

=head1 VERSION

Version 0.12.6

=cut

our $VERSION = '0.12.6';
$VERSION = eval $VERSION;

has infile => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 1,
      coerce   => 1,
);

has outfile => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has allstatsfile => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has illrecstatsfile => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has threads => (
    is        => 'ro',
    isa       => 'Int',
    predicate => 'has_threads',
    lazy      => 1,
    default   => 1,
);

has family_size => (
    is        => 'ro',
    isa       => 'Int',
    predicate => 'has_family_size',
    lazy      => 1,
    default   => 100,
);

has clean => (
    is       => 'ro',
    isa      => 'Bool',
    required => 0,
    default  => 1,
);

sub find_illegitimate_recombination {
    my $self = shift;
    my $statsfile = $self->allstatsfile;

    my ($alignments, $model_dir) = $self->align_features;
 
    my $all_gap_stats = {};
    for my $aln_file (@$alignments) {
	unless (defined $aln_file && -s $aln_file) {
	    unlink $aln_file if -e $aln_file;
	    next;
	}
	$self->find_align_gaps($all_gap_stats, $aln_file);
    }

    $self->collate_gap_stats($all_gap_stats, $statsfile);
    remove_tree( $model_dir, { safe => 1 } );

    return;
}


sub align_features {
    my $self = shift;
    my $threads = $self->threads;
    my $infile = $self->infile->absolute->resolve;

    my ($name, $path, $suffix) = fileparse($infile, qr/\.[^.]*/);
    my $model_dir = File::Spec->catdir($path, 'Tephra_LTR_illrecomb_models');
    unless ( -e $model_dir ) {
	make_path( $model_dir, {verbose => 0, mode => 0771,} );
    }

    my $args = $self->collect_align_args($model_dir);
    my $t0 = gettimeofday();
    my $doms = 0;

    my $logfile = File::Spec->catfile( $model_dir, 'all_illrecomb_muscle_reports.log' );
    open my $log, '>>', $logfile or die "\n[ERROR]: Could not open file: $logfile\n";

    my $pm = Parallel::ForkManager->new($threads);
    local $SIG{INT} = sub {
        $log->warn("Caught SIGINT; Waiting for child processes to finish.");
        $pm->wait_all_children;
        exit 1;
    };

    $pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
			      unlink $data_ref->{data}{seqs}, unlink $data_ref->{data}{log};
			      if ($data_ref->{status} =~ /failed/i) {
				  say $log "\n[WARNING]: ",basename($ident), " failed with exit code: $exit_code\n";
			      }
			      else {
				  my $t1 = gettimeofday();
				  my $elapsed = $t1 - $t0;
				  my $time = sprintf("%.2f",$elapsed/60);
				  say $log basename($ident),
			              " just finished with PID $pid and exit code: $exit_code in $time minutes";
			      }
			} );

    for my $name (keys %$args) {
	$doms++;
	$pm->start($name) and next;
	$SIG{INT} = sub { $pm->finish };
	my $status = $self->capture_cmd($args->{$name}{args});
	$pm->finish(0, { data => $args->{$name}, status => $status });
    }

    $pm->wait_all_children;
    
    my $t2 = gettimeofday();
    my $total_elapsed = $t2 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);

    say $log "\n========> Finished running MUSCLE on $doms families in $final_time minutes";
    close $log;
    unlink $logfile if $self->clean;

    my @aligns;
    for my $k (keys %$args) {
	push @aligns, $args->{$k}{aln};
    }

    return (\@aligns, $model_dir);
}

sub collect_align_args {
    my $self = shift;
    my ($model_dir) = @_;
    my (@full, @aln, %aln_args);

    my $config  = Tephra::Config::Exe->new->get_config_paths;
    my ($muscle) = @{$config}{qw(muscle)};

    my $seqstore = $self->_filter_families_by_size($model_dir);

    for my $fam (keys %$seqstore) {
	my ($name, $path, $suffix) = fileparse($seqstore->{$fam}, qr/\.[^.]*/);
	my $aln = File::Spec->catfile( $model_dir, $name.'_muscle-out.fas' );
	my $log = File::Spec->catfile( $model_dir, $name.'_muscle-out.log' );
	
	my $muscmd  = "$muscle -quiet -in $seqstore->{$fam} -out $aln -log $log";
	$aln_args{$name} = { seqs => $seqstore->{$fam}, args => $muscmd, aln => $aln, log => $log };
    }

    if ($self->clean) {
	unlink $_ for keys %$seqstore;
    }

    return \%aln_args;
}

sub find_align_gaps {
    my $self = shift;
    my ($gap_stats, $aln_file) = @_;

    my $outfile = $self->outfile;;
    my $illrecstatsfile = $self->illrecstatsfile;;

    my ($pos, $gap, $del) = (0, 0, 0);
    my (@indels, @flanking_seqs);

    my $cwd = getcwd();
    open my $illrecstat_fh, '>>', $illrecstatsfile or die "\n[ERROR]: Could not open file: $!";
    open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";

    my ($seqs_in_aln, $count) = $self->split_aln($aln_file);

    for my $fas (@$seqs_in_aln) {
	my $aln_in = Bio::AlignIO->new(-file => $fas, -format => 'fasta');

	my ($fname, $fpath, $fsuffix) = fileparse($fas, qr/\.[^.]*/);
	my $seq_out = File::Spec->catfile( abs_path($fpath), $fname.'_gap_flanking_sequences.fasta' );
	open my $each_out, '>>', $seq_out or die "\n[ERROR]: Could not open file: $seq_out\n";

	while ( my $aln = $aln_in->next_aln() ) {
	    my $aln_len = $aln->length;
	    for my $hash (@{$aln->gap_col_matrix}) {
		for my $key (keys %$hash) {
		    $pos++;
		    if ($hash->{$key} == 1) { # gap columns are coded as '1'
			$gap++;
			push @indels, $pos;
		    }
		}
	    }

	    my ($indel_lengths, $indel_ranges) = $self->get_indel_range(@indels);
	    next unless @$indel_lengths;
	    @indels = ();
	    $self->get_stats($gap_stats, $fas, $pos, $gap, $indel_lengths, $fname);

	    for my $line (@$indel_ranges) {
		$del++;
		my ($indel_len, $indel_spos, $indel_epos) = split /\|\|/, $line;

		# Only analyze repeats flanking deletions > 10 bp; Ma et al. 2004,Genome Research
		next unless $indel_len >= 10; 

		# What we want to do is search 20 bp on either side of the gap.
		# This is based on the 1-15 direct repeats found in Arabidopsis... (Devos et al. 2002)
		# Other plants may be different, so we search 20 bp.
		my $upstream_spos   = $indel_spos - 20;
		my $upstream_epos   = $indel_spos - 1;
		my $downstream_spos = $indel_epos + 1;
		my $downstream_epos = $indel_epos + 20;

		if ($upstream_spos < 1 || $upstream_epos > $aln_len) {
		    #say STDERR "[WARNING]: Deletion $del has a flanking repeat out of bounds.";
		    next;
		}
		if ($downstream_epos > $aln_len) {
		    #say STDERR "[WARNING]: Deletion $del has the downstream flanking repeat out of bounds.";
		    last;
		}

		for my $seq ($aln->each_seq) { # $seq "isa" a Bio::LocatableSeq object 
		    my $upstr_seq   = $seq->subseq($upstream_spos,   $upstream_epos);
		    my $downstr_seq = $seq->subseq($downstream_spos, $downstream_epos);

		    my $upstr_id   = $seq->id.'_upstr-del-'.$del.'_'.$upstream_spos.'-'.$upstream_epos;
		    my $downstr_id = $seq->id.'_downstr-del-'.$del.'_'.$downstream_spos.'-'.$downstream_epos;

		    ## Devos et al. 2002 (I think) defined a threshold of 2bp of non-matching bases
		    ## following a gap
		    my @ugaps = ($upstr_seq =~ /(\-+$)/g);
		    my @dgaps = ($downstr_seq =~ /(^\-+)/g);
		    if (@ugaps < 2 && @dgaps < 2) {
			$self->bl2seq_compare({ indel_length   => $indel_len, 
						upstream_seq   => $upstr_seq, 
						downstream_seq => $downstr_seq, 
						upstream_id    => $upstr_id, 
						downstream_id  => $downstr_id, 
						fas_name       => $fname, 
						fas_path       => $fpath, 
						seqs_fh        => $each_out,
					        stats_fh       => $illrecstat_fh });
		    }
		}
	    }
	}
	$pos = 0;
	$del = 0;
	close $each_out;
	$self->collate($seq_out, $out);
	unlink $fas;
	unlink $seq_out;
    }
    close $out;
    close $illrecstat_fh;
    unlink $aln_file; # if $self->clean;

    return;
}

sub split_aln {
    my $self = shift;
    my ($input) = @_;

    my ($iname, $ipath, $isuffix) = fileparse($input, qr/\.[^.]*/);

    my $kseq = Bio::DB::HTS::Kseq->new($input);
    my $iter = $kseq->iterator();

    my $count  = 0;
    my $fcount = 1;
    my @split_files;
    $iname =~ s/\.fa.*//;
    my $cwd = getcwd();
    
    my $tmpiname = $iname.'_'.$fcount.'_XXXX';
    my ($tmp_fh, $tmp_filename) = tempfile( TEMPLATE => $tmpiname, DIR => $ipath, SUFFIX => '.fasta', UNLINK => 0 );
    push @split_files, $tmp_filename;

    while (my $seq = $iter->next_seq) {
	if ($count % 1 == 0 && $count > 0) {
	    $fcount++;
	    $tmpiname = $iname.'_'.$fcount.'_XXXX';
	    ($tmp_fh, $tmp_filename) = tempfile( TEMPLATE => $tmpiname, DIR => $ipath, SUFFIX => '.fasta', UNLINK => 0 );
	    push @split_files, $tmp_filename;
	}
	say $tmp_fh join "\n", ">".$seq->name, $seq->seq;
	$count++;
    }

    return (\@split_files, $count);
}

sub get_indel_range {
    my $self = shift;
    my @indels = @_;

    my (@indel_ranges, @indel_lengths);
    my $gap_range = [ $indels[0] ];

    for my $indel (@indels[1..$#indels]) {
	if ($indel == $gap_range->[-1] + 1) {
	    push @$gap_range, $indel;
	}
	else {
	    my $gap_length = ($gap_range->[-1] - $gap_range->[0]) + 1;
	    push @indel_lengths, $gap_length;
	    push @indel_ranges, join "||", $gap_length, $gap_range->[0], $gap_range->[-1];
	    $gap_range = [ $indel ];
	}
    }

    return (\@indel_lengths, \@indel_ranges);
}

sub bl2seq_compare {
    my $self = shift;
    my ($args) = @_;
    my ($indel_len, $upstr_seq, $downstr_seq, $upstr_id, $downstr_id, $fname, $fpath, $each_out, $illrecstat_fh) =
	@{$args}{qw(indel_length upstream_seq downstream_seq upstream_id downstream_id fas_name fas_path seqs_fh stats_fh)};

    my $tmpiname = $fname.'_XXXX';
    my ($qtmp_fh, $qtmp_filename) = tempfile( TEMPLATE => $tmpiname, DIR => $fpath, SUFFIX => '.fasta', UNLINK => 0 );
    my ($rtmp_fh, $rtmp_filename) = tempfile( TEMPLATE => $tmpiname, DIR => $fpath, SUFFIX => '.fasta', UNLINK => 0 );
    my ($tmp_outfh, $tmp_outfile) = tempfile( TEMPLATE => $tmpiname, DIR => $fpath, SUFFIX => '.blo',   UNLINK => 0 );

    my $config  = Tephra::Config::Exe->new->get_config_paths;
    my ($blastbin) = @{$config}{qw(blastpath)};
    my $blastn  = File::Spec->catfile($blastbin, 'blastn');

    say $qtmp_fh join "\n", ">".$upstr_id, $upstr_seq;
    say $rtmp_fh join "\n", ">".$downstr_id, $downstr_seq;
    close $qtmp_fh;
    close $rtmp_fh;

    my $blcmd = "$blastn -query $qtmp_filename -subject $rtmp_filename -word_size 4 -outfmt 5 -out $tmp_outfile 2>&1 | ";
    $blcmd .= "egrep -v \"CFastaReader|Error\" | ";
    $blcmd .= "sed \'/^ *\$/d\'"; ## This is a hack for these warnings. I'm throwing away the errors which is why it's a hack but
                                  ## they are not informative. Need to find a workaround for this.
    $self->run_cmd($blcmd);
    unlink $qtmp_filename, $rtmp_filename;

    my $bl2seq_report = Bio::SearchIO->new(-file => $tmp_outfile, -format => 'blastxml');
    
    while ( my $result = $bl2seq_report->next_result ) {
	if ($result->num_hits == 0) {
	    #warn "[WARNING]: No hits found.";
	    unlink $tmp_outfile;
	    return;
	}
	else {
	    #say "HIT FOUND in $upstr_seq and $downstr_seq";# if $debug
	    my $query = $result->query_name();
	    my $qlen  = $result->query_length();
	    while ( my $hit = $result->next_hit ) {
		my $hitid = $hit->name();
		my $hlen  = $hit->length;
		while ( my $hsp = $hit->next_hsp ) {
		    my $hsplen    = $hsp->length('total');
		    my $hstart    = $hsp->start('hit');
		    my $hstop     = $hsp->end('hit');
		    my $qstart    = $hsp->start('query');
		    my $qstop     = $hsp->end('query');
		    my $qstring   = $hsp->query_string;
		    my $hstring   = $hsp->hit_string;
		    my $homstring = $hsp->homology_string;
		    my $hpid      = $hsp->percent_identity;

		    my $qend = $qlen - $qstop;

		    # As defined in Ma et al. 2004 and Devos et al. 2002, 
		    # we want a match of at least 2 bp ($hsplen),
		    # no more than 2 bp from the upstream gap terminus ($qend)
		    # and no more than 2bp from the start of the downstream gap terminus ($hstart).
		    #
		    # Finally, if we have the minimum 2bp match in a 20bp string, 
		    # we should get >= 10% identity ($hpid). This needs some attention bc the match
		    # string is shorter.
		    if ($hsplen > 2 && $hpid >= 10 && $qend <= 2 && $hstart <= 2) {    
			my $qustr  = $upstr_seq;
			my $hitstr = $downstr_seq;
			
			## just testing the indexing for the hit, which is off for some matches
			#my $que_start = index($qustr, $qstring);
			#my $hit_start = index($hitstr, $hstring);
			
			say $each_out join "\n", ">$upstr_id", $qstring;
			say $each_out join "\n", ">$downstr_id", $hstring;
			
			say $illrecstat_fh join "\t", "Query_ID","Hit_ID","HSP_len","Hit_start",
			    "Hit_stop","Query_start","Query_stop","HSP_PID";
			say $illrecstat_fh join "\t", $upstr_id, $downstr_id, $hsplen, $hstart, 
			    $hstop, $qstart, $qstop, $hpid,"\n";

			say $illrecstat_fh "Gap length        : $indel_len";
			say $illrecstat_fh "Query len         : $qlen";
			say $illrecstat_fh "Hit len           : $hlen";
			say $illrecstat_fh "Query string      : $qustr";
			say $illrecstat_fh "Hit string        : $hitstr\n";
			say $illrecstat_fh "Query match string: ",uc($qstring);
			say $illrecstat_fh "Homology string   : $homstring";
			say $illrecstat_fh "Hit match string  : ",uc($hstring),"\n";
		    }
		}
	    }
	}
    }
    unlink $tmp_outfile;

    return;
}

sub get_stats {
    my $self = shift;
    my ($gap_stats, $fas, $pos, $gap, $indel_lengths, $fname) = @_;

    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@$indel_lengths);

    my ($fasname, $faspath, $fassuffix) = fileparse($fas, qr/\.[^.]*/);
    $fasname =~ s/_muscle-out.*//;
    my $count = $stat->count;
    if ($count > 0) {
	my $gap_percent = sprintf("%.2f",$gap/$pos);
	my $ave = $stat->mean; # // 0;
	my $mean = sprintf("%.2f",$ave);
	my $min = $stat->min; # // 0;
	my $max = $stat->max; # // 0;
	my $sum = $stat->sum;

	push @{$gap_stats->{$fasname}}, join "||", $gap, $count, $gap_percent, $mean, $min, $max, $sum;
    }

    return;
}

sub collate_gap_stats {
    my $self = shift;
    my ($gap_stats, $statsfile) = @_;

    open my $gap_stats_fh_out, '>', $statsfile or die "\n[ERROR]: Could not open file: $statsfile\n";
    
    say $gap_stats_fh_out join "\t", "Family_name", "Total_fam_gap_count", "Total_length_of_gaps", 
        "Mean_fam_gap_count (stddev)", "Mean_fam_gap_size (stddev)", "Mean_fam_gap_perc (stdev)", 
        "Mean_fam_gap_size (stddev)", "Mean_gap_min_size (stddev)", "Mean_gap_max_size (stddev)";

    my (@repeat_names, @total_gap_char, @diff_gap_sizes, @gap_char_perc, 
	@mean_gap_size, @min_gap_size, @max_gap_size, %fam_gap_stats);

    for my $fasname (keys %$gap_stats) {
	for my $gap (@{$gap_stats->{$fasname}}) {
	    my @all_gap_stats = split /\|\|/, $gap;

	    push @{$fam_gap_stats{ $fasname }{ total_gap_char }},  $all_gap_stats[0];
	    push @{$fam_gap_stats{ $fasname }{ diff_gap_sizes }},  $all_gap_stats[1];
	    push @{$fam_gap_stats{ $fasname }{ gap_char_perc }},   $all_gap_stats[2];
	    push @{$fam_gap_stats{ $fasname }{ mean_gap_size }},   $all_gap_stats[3];
	    push @{$fam_gap_stats{ $fasname }{ min_gap_size }},    $all_gap_stats[4];
	    push @{$fam_gap_stats{ $fasname }{ max_gap_size }},    $all_gap_stats[5];
	    push @{$fam_gap_stats{ $fasname }{ sum_gap_lengths }}, $all_gap_stats[6];
	}
    }

    for my $family (nsort keys %fam_gap_stats) {
	my $fam_name = $family;
	$fam_name =~ s/\_\d+\_.*//;
	$fam_name =~ s/^all\_//;
	
	my $total_gap_char_stats = Statistics::Descriptive::Full->new;
	my $gap_size_stats       = Statistics::Descriptive::Full->new;
	my $gap_char_perc_stats  = Statistics::Descriptive::Full->new;
	my $mean_gap_size_stats  = Statistics::Descriptive::Full->new;
	my $min_gap_size_stats   = Statistics::Descriptive::Full->new;
	my $max_gap_size_stats   = Statistics::Descriptive::Full->new;
	my $gap_length_sum_stats = Statistics::Descriptive::Full->new;

	$total_gap_char_stats->add_data(@{$fam_gap_stats{$family}{total_gap_char}});
	$gap_size_stats->add_data(@{$fam_gap_stats{$family}{diff_gap_sizes}});
	$gap_char_perc_stats->add_data(@{$fam_gap_stats{$family}{gap_char_perc}});
	$mean_gap_size_stats->add_data(@{$fam_gap_stats{$family}{mean_gap_size}});
	$min_gap_size_stats->add_data(@{$fam_gap_stats{$family}{min_gap_size}});
	$max_gap_size_stats->add_data(@{$fam_gap_stats{$family}{max_gap_size}});
	$gap_length_sum_stats->add_data(@{$fam_gap_stats{$family}{sum_gap_lengths}});

	my $grand_mean_fam_count = $total_gap_char_stats->count;
	my $grand_gap_char_mean  = sprintf("%.2f",$total_gap_char_stats->mean);
	my $grand_gap_size_mean  = sprintf("%.2f",$gap_size_stats->mean);
	my $grand_gap_char_perc  = sprintf("%.2f",$gap_char_perc_stats->mean);
	my $grand_mean_gap_size  = sprintf("%.2f",$mean_gap_size_stats->mean);
	my $grand_gap_size_min   = sprintf("%.2f",$min_gap_size_stats->mean);
	my $grand_gap_size_max   = sprintf("%.2f",$max_gap_size_stats->mean);
	my $gap_lengths_sum      = $gap_length_sum_stats->sum;
	
	my $grand_gap_char_sd      = sprintf("%.2f",$total_gap_char_stats->standard_deviation);
	my $grand_gap_size_sd      = sprintf("%.2f",$gap_size_stats->standard_deviation);
	my $grand_gap_char_perc_sd = sprintf("%.2f",$gap_char_perc_stats->standard_deviation);
	my $grand_gap_size_mean_sd = sprintf("%.2f",$mean_gap_size_stats->standard_deviation);
	my $grand_gap_size_min_sd  = sprintf("%.2f",$min_gap_size_stats->standard_deviation);
	my $grand_gap_size_max_sd  = sprintf("%.2f",$max_gap_size_stats->standard_deviation);
	
	say $gap_stats_fh_out join "\t", $fam_name, $grand_mean_fam_count, $gap_lengths_sum, 
	    $grand_gap_char_mean.' ('.$grand_gap_char_sd.')',
	    $grand_gap_size_mean.' ('.$grand_gap_size_sd.')',
	    $grand_gap_char_perc.' ('.$grand_gap_char_perc_sd.')',
	    $grand_mean_gap_size.' ('.$grand_gap_size_mean_sd.')',
	    $grand_gap_size_min.' ('.$grand_gap_size_min_sd.')',
	    $grand_gap_size_max.' ('.$grand_gap_size_max_sd.')';
    }
    close $gap_stats_fh_out;

    return;
}

sub _filter_families_by_size {
    my $self = shift;
    my ($model_dir) = @_;
    my $infile  = $self->infile;
    my $sthresh = $self->family_size;
    my $lthresh = 1.2e4; # elements over 12kb cannot be aligned reliably

    my ($name, $path, $suffix) = fileparse($infile, qr/\.[^.]*/);
    my $kseq = Bio::DB::HTS::Kseq->new($infile);
    my $iter = $kseq->iterator();

    my (@seqs, $foutfile, $fout);
    while (my $seqobj = $iter->next_seq) {
        my $id = $seqobj->name;
        my $seq = $seqobj->seq;
        my ($family) = ($id =~ /^(RL[CGX]_family\d+)_LTR/);
        if (defined $family) {
            $foutfile = File::Spec->catfile($model_dir, $family.'.fas');
            if (-e $foutfile && openhandle($fout)) {
                say $fout join "\n",">$id", $seq;
            }
            else {
                open $fout, '>', $foutfile or die "\n[ERROR]: Could not open file: $foutfile\n";
                say $fout join "\n", ">$id", $seq;
                push @seqs, $foutfile;
            }
        }
    }

    my ($largest, $seqct) = (0, 0);
    my (%seqstore, @over_thresh);

    for my $seqfile (@seqs) {
	my ($family) = ($seqfile =~ /\/(\w+_family\d+)/);
	my $kseq = Bio::DB::HTS::Kseq->new($seqfile);
	my $iter = $kseq->iterator();
	while (my $seqobj = $iter->next_seq) {
	    $seqct++;
	    my $seq = $seqobj->seq;
	    my $len = length($seq);
	    $largest = $len if $len > $largest;
	}
	if ($seqct <= $sthresh && $largest <= $lthresh) {
	    $seqstore{ $family } = $seqfile;
	}
	else {
	    push @over_thresh, $seqfile;
	}
	($seqct, $largest) = (0, 0);
    }
    unlink @over_thresh; # remove temp files that cannot be aligned

    return \%seqstore;
}

sub _remove_singletons {
    my $self = shift;
    my ($args) = @_;

    my @singles;
    my $seqct = 0;
    for my $name (keys %$args) {
	my $db = $args->{$name}{seqs};
	next unless defined $db && -e $db;
	my $kseq = Bio::DB::HTS::Kseq->new($db);
	my $iter = $kseq->iterator();
	while (my $seqobj = $iter->next_seq) { $seqct++ if defined $seqobj->seq; }
	if ($seqct < 2) {
	    push @singles, $name;
	    unlink $db;
	}
	$seqct = 0;
    }

    delete $args->{$_} for @singles;

    return;
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

    perldoc Tephra::Genome::IllRecombination


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
