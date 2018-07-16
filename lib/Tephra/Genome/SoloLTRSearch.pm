package Tephra::Genome::SoloLTRSearch;

use 5.014;
use Moose;
use File::Spec;
use File::Find;
use File::Basename;
use File::Copy;
use File::Path          qw(make_path remove_tree);
use IPC::System::Simple qw(capture system);
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

with 'Tephra::Role::Run::Any',
     'Tephra::LTR::Role::Utils';

=head1 NAME

Tephra::Genome::SoloLTRSearch - Find solo-LTRs in a refence genome

=head1 VERSION

Version 0.11.1

=cut

our $VERSION = '0.11.1';
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
    my $seqfile  = $self->seqfile;
    my $report   = $self->report; 

    my @sfs;
    find( sub { push @sfs, $File::Find::name if -d && /_copia\z|_gypsy\z|_unclassified\z/ }, $anno_dir);
    croak "\n[ERROR]: Could not find the expected sub-directories ending in 'copia', 'gypsy' and 'unclassified'. ".
	"Please check input. Exiting.\n" unless @sfs; #== 2;

    my $forks = @sfs;
    my $pm = Parallel::ForkManager->new($forks);
    local $SIG{INT} = sub {
        warn "Caught SIGINT; Waiting for child processes to finish.";
        $pm->wait_all_children;
        exit 1;
    };

    $pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
			      if (defined $data_ref) { 
				  my ($model_dir, $reports) = @{$data_ref}{qw(model_dir reports)};
				  $self->_collate_sololtr_reports($reports, $report);

				  if ($self->clean) {
				      remove_tree( $model_dir, { safe => 1} );
				  }
			      }
			} );

    my $dirct = 0;
    for my $dir (nsort @sfs) {
	my $sf = (split /_/, $dir)[-1];
	$pm->start($sf) and next;
	$SIG{INT} = sub { $pm->finish };
	my $data_ref = $self->run_sf_search($sf, $dir, $dirct, $report);
	$dirct++;
	$pm->finish(0, $data_ref);
    }
    $pm->wait_all_children;

    my $soloct = $self->_check_report_summary($report);

    if ($soloct > 0) {
	$self->write_sololtr_gff($report);
    }
    else {
	say STDERR "\n[WARNING]: No solo-LTRs were found so none will be reported.\n";
	unlink $report, $seqfile;
    }

    return;
}

sub run_sf_search {
    my $self = shift;
    my ($sf, $dir, $dirct, $hmmsearch_summary) = @_;
    my $genome = $self->genome->absolute->resolve;
    my $debug  = $self->debug;
    
    my ($hmmbuild, $hmmsearch) = $self->_find_hmmer;

    print STDERR "Getting LTR alignments for ",ucfirst($sf),"...." if $debug;
    my $ltr_aln_files = $self->_get_ltr_alns($dir);    
    say STDERR "done with alignments." if $debug;
    #dd $ltr_aln_files and exit;

    ## need masked genome here
    if (!defined $ltr_aln_files || @$ltr_aln_files < 1) {
	say STDERR "\n[ERROR]: The expected alignments files were not found for ",ucfirst($sf).
	    ". Skipping solo-LTR search for this superfamily.";
	return undef;
    }
        
    print STDERR "Getting alignment statistics for ",ucfirst($sf),"..." if $debug;
    my $aln_stats = $self->_get_aln_len($ltr_aln_files); # return a hash-ref
    say STDERR "done with alignment statistics." if $debug;
    #dd $aln_stats;

    # make one directory
    my $model_dir = File::Spec->catdir($dir, 'Tephra_LTR_exemplar_models');
    unless ( -e $model_dir ) {
	make_path( $model_dir, {verbose => 0, mode => 0771,} );
    }
        
    print STDERR "Building LTR exemplar models for ",ucfirst($sf),"..." if $debug;
    for my $ltr_aln (nsort_by { m/family(\d+)/ and $1 } @$ltr_aln_files) {
	$self->build_model($ltr_aln, $model_dir, $hmmbuild);        
	unlink $ltr_aln;    
    }
    say STDERR "done with models." if $debug;
    
    my @ltr_hmm_files;
    find( sub { push @ltr_hmm_files, $File::Find::name if -f and /\.hmm$/ }, $model_dir);
        
    print STDERR "Search genome with models for ",ucfirst($sf),"..." if $debug;
    $self->do_parallel_search($hmmsearch, $genome, \@ltr_hmm_files, $model_dir, $aln_stats);
    say STDERR "done searching with LTR models." if $debug;

    print STDERR "Collating results and writing GFF..." if $debug;
    my @reports;
    find( sub { push @reports, $File::Find::name if -f and /\.txt$/ }, $model_dir );

    return { reports => \@reports, model_dir => $model_dir };
}
    

sub do_parallel_search {
    my $self = shift;
    my ($hmmsearch, $genome, $ltr_hmm_files, $model_dir, $aln_stats) = @_;
    my $threads = $self->threads;

    my %reports;
    my $t0 = gettimeofday();
    my $model_num = @$ltr_hmm_files;

    my ($gname, $gpath, $gsuffix) = fileparse($genome, qr/\.[^.]*/);
    my $logfile = File::Spec->catfile($model_dir, 'all_solo-ltr_searches.log');
    open my $log, '>>', $logfile or die "\n[ERROR]: Could not open file: $logfile\n";

    my $thr;
    if ($threads % 2 == 0) {
        $thr = sprintf("%.0f",$threads/2);
    }
    elsif ($threads-1 % 2 == 0) {
        $thr = sprintf("%.0f",$threads-1/2);
    }
    else {
        $thr = 1;
    }

    my $pm = Parallel::ForkManager->new($threads);
    local $SIG{INT} = sub {
        warn "Caught SIGINT; Waiting for child processes to finish.";
        $pm->wait_all_children;
        exit 1;
    };

    $pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
			      for my $hmmrep (sort keys %$data_ref) {
				  unlink $hmmrep;
			      }
			      my $t1 = gettimeofday();
			      my $elapsed = $t1 - $t0;
			      my $time = sprintf("%.2f",$elapsed/60);
			      say $log basename($ident),
			        " just finished with PID $pid and exit code: $exit_code in $time minutes";
			} );

    for my $hmm (nsort @$ltr_hmm_files) {
	$pm->start($hmm) and next;
	$SIG{INT} = sub { $pm->finish };
	
	my $hmmsearch_out = $self->search_with_models($gname, $hmm, $hmmsearch, $genome, $aln_stats);
	$reports{$hmmsearch_out} = 1;

	$pm->finish(0, \%reports);
    }

    $pm->wait_all_children;

    my $t2 = gettimeofday();
    my $total_elapsed = $t2 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);

    say $log "\n========> Finished hmmsearch on $model_num solo-LTR models in $final_time minutes";
    close $log;
}

sub write_sololtr_gff {
    my $self = shift;
    my ($hmmsearch_summary) = @_;
    my $genome  = $self->genome->absolute->resolve;
    my $outfile = $self->outfile;;

    my $seqlen = $self->_get_seq_len($genome);
    open my $in, '<', $hmmsearch_summary or die "\n[ERROR]: Could not open file: $hmmsearch_summary\n";
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
	my ($query, $query_length, $number_of_hits, $hit_name, $perc_coverage, $hsp_length, 
	    $hsp_perc_ident, $hsp_query_start, $hsp_query_end, $hsp_hit_start, $hsp_hit_end,
	    $search_type) = split /\t/, $line;
	$query =~ s/_ltrs_muscle-out//;

	if (exists $seqlen->{$hit_name}) {
	    $ct++;
	    say $out join "\t", $hit_name, 'Tephra', 'solo_LTR', $hsp_hit_start, $hsp_hit_end, 
	        '.', '?', '.', "ID=solo_LTR$ct;match_id=$query;Name=solo_LTR;Ontology_term=SO:0001003";
	}
	else {
	    croak "\n[ERROR]: $hit_name not found in $genome. This is a bug, please report it. Exiting.\n";
	}
    }
    close $in;
    close $out; 

    return;
}

sub build_model {
    my $self = shift;
    my ($aln, $model_dir, $hmmbuild) = @_;

    my $hmmname = basename($aln);
    $hmmname =~ s/\.aln$/\.hmm/;
    my $model_path = File::Spec->catfile($model_dir, $hmmname);
    # hmmer3 has --dna and --cpu options and not -g (global)
    my $hmm_cmd = "$hmmbuild -f --nucleic $model_path $aln";
    say STDERR "DEBUG: $hmm_cmd" if $self->debug;

    my $status = $self->capture_cmd($hmm_cmd);
    unlink $aln && return if $status =~ /failed/i;

    return;
}

sub search_with_models {
    my $self = shift;
    my ($gname, $hmm, $hmmsearch, $fasta, $aln_stats) = @_;

    my $hmmsearch_out = $hmm;
    $hmmsearch_out =~ s/\.hmm.*$//;
    $hmmsearch_out .= '_'.$gname.'.hmmer';

    ## no -o option in hmmer2
    my $hmmsearch_cmd = "$hmmsearch --cpu 1 $hmm $fasta > $hmmsearch_out";
    say STDERR "DEBUG: $hmmsearch_cmd" if $self->debug;
    $self->run_cmd($hmmsearch_cmd);

    if (-e $hmmsearch_out) {   
	$self->write_hmmsearch_report($aln_stats, $hmmsearch_out);
    }
    unlink $hmm;

    return $hmmsearch_out;
}

sub write_hmmsearch_report {
    my $self = shift;
    my ($aln_stats, $search_report) = @_;
    my $genome     = $self->genome->absolute->resolve;
    my $match_pid  = $self->percentid;
    my $match_len  = $self->matchlen;
    #my $match_pcov = $self->percentcov;
    #dd $aln_stats;

    my $parsed = $search_report;
    $parsed =~ s/\.hmmer$/\_hmmer_parsed.txt/;
    my $element = $search_report;

    my $model_type = 'local';

    my ($gname, $gpath, $gsuffix) = fileparse($genome, qr/\.[^.]*/);
    my ($name, $path, $suffix)    = fileparse($search_report, qr/\.[^.]*/);
    $element =~ s/_$gname*//;
    open my $out, '>>', $parsed or die "\n[ERROR]: Could not open file: $parsed\n";

    my ($seq, $seqfile);
    if ($self->seqfile) {
	$seqfile = $self->seqfile;
	open $seq, '>>', $seqfile or die "\n[ERROR]: Could not open file: $seqfile";
    }

    my $hmmerin = Bio::SearchIO->new(-file => $search_report, -format => 'hmmer');

    my ($positions, $matches) = (0, 0);
    while ( my $result = $hmmerin->next_result() ) {
	my $query    = $result->query_name();
	my $num_hits = $result->num_hits();
	my $qlen     = $aln_stats->{$query};
	die "\n[ERROR]: Could not determine query length for: $query" 
	    unless defined $qlen;
	while ( my $hit = $result->next_hit() ) {
	    my $hitid = $hit->name();
	    while ( my $hsp = $hit->next_hsp() ) {
		my $percent_q_coverage = sprintf("%.2f", $hsp->length('query')/$qlen * 100);
		my $hspgaps = $hsp->gaps;
		my $hsplen  = $hsp->length('total');
		my $hstart  = $hsp->start('hit');
		my $hstop   = $hsp->end('hit');
		my $qstart  = $hsp->start('query');
		my $qstop   = $hsp->end('query');
		my $qstring = $hsp->query_string;
		my $qhsplen = $hsp->length('query');
		my $pid     = $hsp->percent_identity;

		if (exists $aln_stats->{$query}) {
		    if ($hsplen >= $match_len && $pid >= $match_pid) { #$hsplen >= $aln_stats->{$query} * ($match_pcov/100) ) {
			my $qid = $query =~ s/_ltrseqs_muscle-out//r; # non-destructive substitution in v5.14+
			$matches++;
			say $out join "\t", $qid, $aln_stats->{$query}, $num_hits, $hitid, 
			    $percent_q_coverage, $hsplen, $pid, $qstart, $qstop, $hstart, $hstop, $model_type;
			
			if ($seqfile) {
			    my $seqid = join '_', '>'.$qid, $hitid, $hstart, $hstop; 
			    ## It makes more sense to show the location of the hit
			    ## Also, this would pave the way for creating a gff of solo-LTRs
			    ## my $seqid = ">".$query."|".$hitid."_".$hstart."-".$hstop
			    say $seq join "\n", $seqid, $qstring;
			}
		    }
		}
	    }
	}
    }
    close $out;
    close $seq;
    #unlink $seqfile, $parsed if $matches == 0;
    return;
}

sub _get_ltr_alns {
    my $self = shift;
    my ($dir) = @_;
    my $numfams = $self->numfamilies;
    my $allfams = $self->fullanalysis;

    my (@ltrseqs, @aligns);

    my $ltrseqs = $self->get_exemplar_ltrs_for_sololtrs({ input_dir => $dir, full_analysis => $allfams });
    return unless defined $ltrseqs && @$ltrseqs;

    # This is where families are filtered by size. Since largest families come first,
    # a simple sort will filter the list.
    my $aln_ct = 0;
    for my $ltrseq (nsort_by { m/family(\d+)/ and $1 } @$ltrseqs) {
	$aln_ct++;
	my ($name, $path, $suffix) = fileparse($ltrseq, qr/\.[^.]*/);
	my $tre = File::Spec->catfile( abs_path($path), $name.'.dnd' );
	my $aln = File::Spec->catfile( abs_path($path), $name.'_muscle-out.aln' );
	my $log = File::Spec->catfile( abs_path($path), $name.'_muscle-out.log' );
     
	my $muscmd = "muscle -quiet -clwstrict -in $ltrseq -out $aln -log $log";
        if ($allfams) {
	    my $status = $self->capture_cmd($muscmd);
	    unlink $ltrseq && next if $status =~ /failed/i;
	    unlink $log;
	    push @aligns, $aln;
	}
	else {
            if ($numfams >= $aln_ct) {
		my $status = $self->capture_cmd($muscmd);
		unlink $ltrseq && next if $status =~ /failed/i;
		unlink $log;
		push @aligns, $aln;
            }
        }
    }
    unlink $_ for @$ltrseqs;

    return \@aligns;
}

sub _collate_sololtr_reports {
    my $self = shift;
    my ($files, $outfile) = @_;

    ##TODO: use lexical vars instead for array indexes
    my (%seen, %parsed_alns);
    for my $file (@$files) {
	open my $fh_in, '<', $file or die "\n[ERROR]: Could not open file: $file\n";
	while (my $line = <$fh_in>) {
	    chomp $line;
	    next if $line =~ /^#/;
	    my @f = split /\t/, $line;
	    if (exists $seen{$f[3]}{$f[9]}) {
		my ($len, $pid) = split /\|\|/, $seen{$f[3]}{$f[9]};
		if ($f[10]-$f[9] > $len && $f[6] > $pid) {
		    $parsed_alns{$f[3]}{$f[9]} =  join "||", @f;
		}
	    }
	    $parsed_alns{$f[3]}{$f[9]} =  join "||", @f;
	    $seen{$f[3]}{$f[9]} = join "||", $f[10]-$f[9], $f[6];
	}
	close $fh_in;
    }
   
    open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
    if (! -s $outfile) { 
	say $out join "\t", "#query", "query_length", "number_of_hits", "hit_name",
            "perc_coverage","hsp_length", "hsp_perc_ident","hsp_query_start", "hsp_query_end",
            "hsp_hit_start", "hsp_hit_end","search_type";

    }

    for my $chr (nsort keys %parsed_alns) {
	for my $start (sort { $a <=> $b } keys %{$parsed_alns{$chr}}) {
	    my @f = split /\|\|/, $parsed_alns{$chr}{$start};
	    say $out join "\t", @f;
	}
    }
    close $out;

    return;
}

sub _check_report_summary {
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
    
sub _find_hmmer {
    my $self = shift; 

    my $config      = Tephra::Config::Exe->new->get_config_paths;
    my ($hmmer2bin) = @{$config}{qw(hmmer2bin)};
    my $hmmbuild    = File::Spec->catfile($hmmer2bin, 'hmmbuild');
    my $hmmsearch   = File::Spec->catfile($hmmer2bin, 'hmmsearch');

    if (-e $hmmbuild && -x $hmmbuild &&
	-e $hmmsearch && -x $hmmsearch) {
	return ($hmmbuild, $hmmsearch);
    }
    else {
	croak "\n[ERROR]: Could not get HMMERv2 PATH. This indicates that Tephra was not configured correctly. Exiting.\n";
    }
}

sub _get_seq_len {
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

sub _get_aln_len {
    my $self = shift;
    my ($aln_files) = @_;
    
    my %aln_stats;

    for my $aligned (@$aln_files) {
	my $ltr_retro = basename($aligned);
	$ltr_retro =~ s/\.aln//;
	my $aln_in = Bio::AlignIO->new(-file => $aligned, -format => 'clustalw');

	while ( my $aln = $aln_in->next_aln() ) {
	    $aln_stats{$ltr_retro} = $aln->length;
	}	
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
