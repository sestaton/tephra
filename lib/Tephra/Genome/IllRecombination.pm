package Tephra::Genome::IllRecombination;

use 5.010;
use Moose;
use Cwd;
use File::Spec;
use File::Find;
use File::Basename;
use File::Temp;
use Path::Class::File;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SearchIO;
use Parallel::ForkManager;
use Statistics::Descriptive;
use Time::HiRes qw(gettimeofday);
use Log::Any qw($log);
use Tephra::Config::Exe;
use namespace::autoclean;

with 'Tephra::Role::Util';

=head1 NAME

Tephra::Genome::IllRecombination - Calculate illigetimate recombination in a genome

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';
$VERSION = eval $VERSION;

has dir => (
      is       => 'ro',
      isa      => 'Path::Class::Dir',
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

has clean => (
    is       => 'ro',
    isa      => 'Bool',
    required => 0,
    default  => 1,
);

sub align_features {
    my $self = shift;
    my $threads = $self->threads;
    my $dir = $self->dir;

    my $args = $self->collect_align_args($dir);
    $self->_remove_singletons($args);

    my $t0 = gettimeofday();
    my $doms = 0;
    my $logfile = File::Spec->catfile($dir, 'all_muscle_reports.log');
    open my $log, '>>', $logfile or die "\nERROR: Could not open file: $logfile\n";

    my $pm = Parallel::ForkManager->new($threads);
    $pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
			      my $t1 = gettimeofday();
			      my $elapsed = $t1 - $t0;
			      my $time = sprintf("%.2f",$elapsed/60);
			      say $log basename($ident),
			          " just finished with PID $pid and exit code: $exit_code in $time minutes";
			} );

    for my $name (keys %$args) {
	$doms++;
	$pm->start($name) and next;
	$self->capture_cmd($args->{$name}{args});
	$pm->finish(0);
    }

    $pm->wait_all_children;
    
    my $t2 = gettimeofday();
    my $total_elapsed = $t2 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);

    say $log "\n========> Finished running MUSCLE on $doms families in $final_time minutes";
    close $log;

    my @aligns;
    for my $k (keys %$args) {
	push @aligns, $args->{$k}{aln};
    }
    return \@aligns;
}

sub collect_align_args {
    my $self = shift;
    my ($dir) = @_;
    my (@full, @aln, %aln_args);
    find( sub { push @full, $File::Find::name if -f and /(?:family\d+|singleton_families).fasta$/ }, $dir);

    for my $fam (@full) {
	my ($name, $path, $suffix) = fileparse($fam, qr/\.[^.]*/);
	my $aln = File::Spec->catfile($path, $name."_muscle-out.fas");
	my $log = File::Spec->catfile($path, $name."_muscle-out.log");
	
	my $clwcmd  = "muscle -in $fam -out $aln 2>$log";
	$aln_args{$name} = { seqs => $fam, args => $clwcmd, aln => $aln };
    }

    return \%aln_args;
}

sub find_illigetimate_recombination {
    my $self = shift;

    my $alignments = $self->align_features;

    for my $aln_file (@$alignments) {
	$self->find_align_gaps($aln_file);
    }
}

sub find_align_gaps {
    my $self = shift;
    my ($aln_file) = @_;

    my $outfile   = $self->outfile;
    my $statsfile = $self->allstatsfile;
    my $illrecstatsfile = $self->illrecstatsfile;

    my $dr_pid;
    my $gap_stats;
    my $help;

    my $pos = 0;
    my $gap = 0;
    my $del = 0;
    my @indels;
    my @flanking_seqs;
    
    $dr_pid = defined($dr_pid) ? $dr_pid : '10'; ## make class attribute
    my $cwd = getcwd();
    open my $illrecstat_fh, '>>', $illrecstatsfile or die "\nERROR: Could not open file: $!";
    
    my $statstmp = $statsfile.".tmp";
    open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n";
    open my $stats_out_tmp, '>>', $statstmp or die "\nERROR: Could not open file: $statstmp\n";
    say $stats_out_tmp join "\t", "#LTR_retro_name", "Total_num_gap_sites", "Num_diff_gap_sizes", 
        "Percent_gap_sites", "Mean_gap_size", "Min_gap_size", "Max_gap_size";

    my ($seqs_in_aln, $count) = $self->split_aln($aln_file);
    
    for my $fas (@$seqs_in_aln) {
	my $aln_in = Bio::AlignIO->new(-fh => \*$fas, -format => 'fasta');

	my ($fname, $fpath, $fsuffix) = fileparse($fas, qr/\.[^.]*/);
	my $seq_out = File::Spec->catfile($fpath, $fname."_gap_flanking_sequences.fasta");
	open my $each_out, '>>', $seq_out or die "\nERROR: Could not open file: $seq_out\n";

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
	    @indels = ();

	    $gap_stats = $self->get_stats($fas,$pos,$gap,$indel_lengths,$fname) if $statsfile;

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
		    warn "Deletion $del has a flanking repeat out of bounds.";
		    next;
		}
		if ($downstream_epos > $aln_len) {
		    warn "Deletion $del has the downstream flanking repeat out of bounds.";
		    last;
		}

		for my $seq ($aln->each_seq) { # $seq "isa" a Bio::LocatableSeq object 
		    my $upstr_seq   = $seq->subseq($upstream_spos,   $upstream_epos);
		    my $downstr_seq = $seq->subseq($downstream_spos, $downstream_epos);

		    my $upstr_id   = $seq->id."_upstr-del-".$del."_".$upstream_spos."-".$upstream_epos;
		    my $downstr_id = $seq->id."_downstr-del-".$del."_".$downstream_spos."-".$downstream_epos;

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
	$self->collate($gap_stats, $stats_out_tmp);
	unlink $fas;
	unlink $seq_out;
	unlink $gap_stats;
    }
    close $out;
    close $stats_out_tmp;
    close $illrecstat_fh;
    
    $self->collate_gap_stats($statstmp, $statsfile);
    unlink $statstmp;
}

sub split_aln {
    my $self = shift;
    my ($input) = @_;

    my ($iname, $ipath, $isuffix) = fileparse($input, qr/\.[^.]*/);

    my $seq_in  = Bio::SeqIO->new(-file  => $input, -format => 'fasta');
    my $count = 0;
    my $fcount = 1;
    my @split_files;
    $iname =~ s/\.fa.*//;
    my $cwd = getcwd();
    
    my $tmpiname = $iname."_".$fcount."_XXXX";
    my $fname = File::Temp->new( TEMPLATE => $tmpiname,
				 DIR      => $ipath,
				 UNLINK   => 0,
				 SUFFIX   => ".fasta");

    my $seq_out = Bio::SeqIO->new(-file   => ">$fname", -format => 'fasta');

    push @split_files, $fname;
    while (my $seq = $seq_in->next_seq) {
	if ($count % 1 == 0 && $count > 0) {
	    $fcount++;
	    $tmpiname = $iname."_".$fcount."_XXXX";
	    $fname = File::Temp->new( TEMPLATE => $tmpiname,
				      DIR      => $ipath,
				      UNLINK   => 0,
				      SUFFIX   => ".fasta");

	    $seq_out = Bio::SeqIO->new(-file   => ">$fname",
				       -format => 'fasta');

	    push @split_files, $fname;
	}
	$seq_out->write_seq($seq);
	$count++;
    }

    return (\@split_files, $count);
}

sub get_indel_range {
    my $self = shift;
    my @indels = @_;
    my @indel_ranges;
    my @indel_lengths;

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
    my ($indel_len, $upstr_seq, $downstr_seq, $upstr_id, $downstr_id, $fname, $fpath, $each_out, $illrecstat_fh)
	= @{$args}{qw(indel_length upstream_seq downstream_seq upstream_id downstream_id fas_name fas_path seqs_fh stats_fh)};

    my $qname = File::Temp->new( TEMPLATE => $fname."_XXXX",
				 DIR      => $fpath,
				 UNLINK   => 0,
				 SUFFIX   => ".fasta");

    my $rname = File::Temp->new( TEMPLATE => $fname."_XXXX",
				 DIR      => $fpath,
				 UNLINK   => 0,
				 SUFFIX   => ".fasta");

    my $out = File::Temp->new( TEMPLATE => $fname."_XXXX",
			       DIR      => $fpath,
			       UNLINK   => 0,
			       SUFFIX   => ".blo");

    my $outfile = $out->filename;
    my $qfile   = $qname->filename;
    my $rfile   = $rname->filename;
    my $config  = Tephra::Config::Exe->new->get_config_paths;
    my ($blastbin) = @{$config}{qw(blastpath)};
    my $blastn = File::Spec->catfile($blastbin, 'blastn');

    say $qname join "\n", ">".$upstr_id, $upstr_seq;
    say $rname join "\n", ">".$downstr_id, $downstr_seq;
    my $blcmd = "$blastn -query $qname -subject $rname -word_size 4 -outfmt 5 -out $outfile 2>&1 | ";
    $blcmd .= "grep -v \"CFastaReader: Hyphens are invalid and will be ignored\""; ## what is worse, this hack or the stupid warnings?
    #say STDERR $blcmd;
    $self->run_cmd($blcmd);
    unlink $qname, $rname;

    my $bl2seq_report = Bio::SearchIO->new(-file => $outfile, -format => 'blastxml');
    
    while ( my $result = $bl2seq_report->next_result ) {
	if ($result->num_hits == 0) {
	    #warn "No hits found.";
	    unlink $outfile;
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
			my $que_start = index($qustr, $qstring);
			my $hit_start = index($hitstr, $hstring);
			
			say $each_out join "\n", ">$upstr_id", $qstring;
			say $each_out join "\n", ">$downstr_id", $hstring;
			
			say $illrecstat_fh join "\t", "Query_ID","Hit_ID","HSP_len","Hit_start",
			    "Hit_stop","Query_start","Query_stop","HSP_PID";
			say $illrecstat_fh join "\t", $query,$hitid,$hsplen,$hstart,$hstop,$qstart,$qstop,$hpid,"\n";

			say $illrecstat_fh "Gap length        : $indel_len";
			say $illrecstat_fh "Query len         : $qlen";
			say $illrecstat_fh "Hit len           : $hlen";
			say $illrecstat_fh "Query string      : $qustr\n";
			say $illrecstat_fh "Query match string: ",uc($qstring);
			say $illrecstat_fh "Homology string   : $homstring";
			say $illrecstat_fh "Hit match string  : ",uc($hstring),"\n";
			say $illrecstat_fh "Hit string        : $hitstr\n";
		    }
		}
	    }
	}
    }
    unlink $outfile;
}

sub get_stats {
    my $self = shift;
    my ($fas, $pos, $gap, $indel_lengths, $fname) = @_;

    my $gap_stats = $fname."_gap_stats.txt";
    open my $statsout, '>', $gap_stats or die "\nERROR: Could not open file: $gap_stats\n";
    my $stat = Statistics::Descriptive::Full->new;

    $stat->add_data(@$indel_lengths);

    my ($fasname, $faspath, $fassuffix) = fileparse($fas, qr/\.[^.]*/);
    my $count = $stat->count;
    my $gap_percent = sprintf("%.2f",$gap/$pos);
    my $mean = sprintf("%.2f",$stat->mean);
    my $min = $stat->min;
    my $max = $stat->max;

    say $statsout join "\t", $fasname, $gap, $count, $gap_percent, $mean, $min, $max;

    close $statsout;
    return $gap_stats;
}

sub collate {
    my $self = shift;
    my ($file_in, $fh_out) = @_;
    my $lines = do { 
	local $/ = undef; 
	open my $fh_in, '<', $file_in or die "\nERROR: Could not open file: $file_in\n";
	<$fh_in>;
    };
    print $fh_out $lines;
}

sub collate_gap_stats {
    my $self = shift;
    my ($statstmp, $statsfile) = @_;

    open my $gap_stats_fh_in, '<', $statstmp or die "\nERROR: Could not open file: $statstmp\n";
    open my $gap_stats_fh_out, '>', $statsfile or die "\nERROR: Could not open file: $statsfile\n";

    my (@repeat_names, @total_gap_char, @diff_gap_sizes, @gap_char_perc, 
	@mean_gap_size, @min_gap_size, @max_gap_size);

    while (<$gap_stats_fh_in>) {
	if (/^\#/) {
	    print $gap_stats_fh_out $_;
	}
	else {
	    my @all_gap_stats = split /\t/, $_;
	    push @repeat_names,   $all_gap_stats[0];
	    push @total_gap_char, $all_gap_stats[1];
	    push @diff_gap_sizes, $all_gap_stats[2];
	    push @gap_char_perc,  $all_gap_stats[3];
	    push @mean_gap_size,  $all_gap_stats[4];
	    push @min_gap_size,   $all_gap_stats[5];
	    push @max_gap_size,   $all_gap_stats[6];

	    print $gap_stats_fh_out $_;
	}
    }

    my $fam_name = pop @repeat_names;
    $fam_name =~ s/\_\d+\_.*//;
    $fam_name =~ s/^all\_//;

    my $total_gap_char_stats = Statistics::Descriptive::Full->new;
    my $gap_size_stats       = Statistics::Descriptive::Full->new;
    my $gap_char_perc_stats  = Statistics::Descriptive::Full->new;
    my $mean_gap_size_stats  = Statistics::Descriptive::Full->new;
    my $min_gap_size_stats   = Statistics::Descriptive::Full->new;
    my $max_gap_size_stats   = Statistics::Descriptive::Full->new;

    $total_gap_char_stats->add_data(@total_gap_char);
    $gap_size_stats->add_data(@diff_gap_sizes);
    $gap_char_perc_stats->add_data(@gap_char_perc);
    $mean_gap_size_stats->add_data(@mean_gap_size);
    $min_gap_size_stats->add_data(@min_gap_size);
    $max_gap_size_stats->add_data(@max_gap_size);

    my $grand_mean_fam_count = $total_gap_char_stats->count;
    my $grand_gap_char_mean  = sprintf("%.2f",$total_gap_char_stats->mean);
    my $grand_gap_size_mean  = sprintf("%.2f",$gap_size_stats->mean);
    my $grand_gap_char_perc  = sprintf("%.2f",$gap_char_perc_stats->mean);
    my $grand_mean_gap_size  = sprintf("%.2f",$mean_gap_size_stats->mean);
    my $grand_gap_size_min   = sprintf("%.2f",$min_gap_size_stats->mean);
    my $grand_gap_size_max   = sprintf("%.2f",$max_gap_size_stats->mean);

    my $grand_gap_char_sd      = sprintf("%.2f",$total_gap_char_stats->standard_deviation);
    my $grand_gap_size_sd      = sprintf("%.2f",$gap_size_stats->standard_deviation);
    my $grand_gap_char_perc_sd = sprintf("%.2f",$gap_char_perc_stats->standard_deviation);
    my $grand_gap_size_mean_sd = sprintf("%.2f",$mean_gap_size_stats->standard_deviation);
    my $grand_gap_size_min_sd  = sprintf("%.2f",$min_gap_size_stats->standard_deviation);
    my $grand_gap_size_max_sd  = sprintf("%.2f",$max_gap_size_stats->standard_deviation);

    say $gap_stats_fh_out "\nFamily_name\tTotal_fam_gap_count\tMean_fam_gap_count ".
	"(stddev)\tMean_fam_gap_size (stddev)\tMean_fam_gap_perc ".
	"(stdev)\tMean_fam_gap_size (stddev)\tMean_gap_min_size ".
	"(stddev)\tMean_gap_max_size (stddev)";

    say $gap_stats_fh_out $fam_name."\t".$grand_mean_fam_count."\t".
	$grand_gap_char_mean." (".$grand_gap_char_sd.")"."\t".
	$grand_gap_size_mean." (".$grand_gap_size_sd.")"."\t".
	$grand_gap_char_perc." (".$grand_gap_char_perc_sd.")"."\t".
	$grand_mean_gap_size." (".$grand_gap_size_mean_sd.")"."\t".
	$grand_gap_size_min." (".$grand_gap_size_min_sd.")"."\t".
	$grand_gap_size_max." (".$grand_gap_size_max_sd.")"."\n";


    close $gap_stats_fh_in;
    close $gap_stats_fh_out;
}

sub _cnv_align_fmt {
    my $self = shift;
    my ($aln) = @_;

    my ($name, $path, $suffix) = fileparse($aln, qr/\.[^.]*/);
    my $fas = File::Spec->catfile($path, $name.".fas");

    my $seqin  = Bio::AlignIO->new(-file => $aln, -format => 'clustalw');
    my $seqout = Bio::AlignIO->new(-file => ">$fas", -format => 'fasta');
    
    my $seqct = 0;
    my $index = 0; # this is a horrible hack...need to get original ids

    while (my $seqobj = $seqin->next_aln) {
	for my $seq ($seqobj->each_seq) {
	    my $id = $seq->id;
	    $id .= "_$index";
	    $seq->id($id);
	    $index++;
	}
	$seqout->write_aln($seqobj);
    }
    
    return $fas;
}

sub _remove_singletons {
    my $self = shift;
    my ($args) = @_;

    my @singles;
    my $seqct = 0;
    for my $name (keys %$args) {
	my $db = $args->{$name}{seqs};
	my $seqio = Bio::SeqIO->new( -file => $db, -format => 'fasta' );
	while (my $seqobj = $seqio->next_seq) { $seqct++ if defined $seqobj->seq; }
	if ($seqct < 2) {
	    push @singles, $name;
	    unlink $db;
	}
	$seqct = 0;
    }

    delete $args->{$_} for @singles;
}

=head1 AUTHOR

S. Evan Staton, C<< <statonse at gmail.com> >>

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
