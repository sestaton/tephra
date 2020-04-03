package Tephra::NonLTR::NonLTRSearch;

use 5.014;
use Moose;
use File::Spec;
use File::Find;
use File::Basename;
use File::Path qw(make_path remove_tree);
use Cwd        qw(abs_path);
use Bio::DB::HTS::Kseq;
use Parallel::ForkManager;
use Tephra::NonLTR::RunHMM;
use Tephra::NonLTR::Postprocess;
use Tephra::NonLTR::QValidation;
use Tephra::NonLTR::SeqUtils;
use Tephra::Config::Exe;
use namespace::autoclean;
#use Data::Dump::Color;

=head1 NAME

Tephra::NonLTR::NonLTRSearch - Search a genome for non-LTR retrotransposons

=head1 VERSION

Version 0.12.5

=cut

our $VERSION = '0.12.5';
$VERSION = eval $VERSION;

has genome  => ( is => 'ro', isa => 'Maybe[Str]', required => 1 );
has outdir  => ( is => 'ro', isa => 'Maybe[Str]', required => 0 );
has pdir    => ( is => 'ro', isa => 'Maybe[Str]', required => 0 );
has threads => ( is => 'ro', isa => 'Int',  predicate => 'has_threads', lazy => 1, default => 1 );
has verbose => ( is => 'ro', isa => 'Bool', predicate => 'has_debug',   lazy => 1, default => 0 );

sub find_nonltrs {
    my $self = shift;
    my $main_data_dir = $self->outdir;
    #my $program_dir   = $self->pdir;
    my $genome = $self->genome;
    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($phmm_dir) = @{$config}{qw(modeldir)};

    my ($gname, $gpath, $gsuffix) = fileparse($genome, qr/\.[^.]*/);
    my $genome_dir = File::Spec->catdir( abs_path($gpath), $gname.'_genome' );
    $main_data_dir //= File::Spec->catdir( abs_path($gpath), $gname.'_nonLTRs' );

    unless ( -d $genome_dir ) {
	make_path( $genome_dir, {verbose => 0, mode => 0771,} );
    }

    unless ( -d $main_data_dir ) {
	make_path( $main_data_dir, {verbose => 0, mode => 0771,} );
    }

    my $plus_out_dir  = File::Spec->catdir($main_data_dir, 'f');
    my $minus_out_dir = File::Spec->catdir($main_data_dir, 'b');
    my $minus_dna_dir = $genome_dir.'_b';

    unless (-e $minus_dna_dir) {
	make_path( $minus_dna_dir, {verbose => 0, mode => 0771,} );
    }

    # set up parallel processing of both strands
    my $pm = Parallel::ForkManager->new(2);
    local $SIG{INT} = sub {
        warn("Caught SIGINT; Waiting for child processes to finish.");
        $pm->wait_all_children;
        exit 1;
    };

    # create data directories for forward strand
    my $fasfiles = $self->_split_genome($genome, $genome_dir);
    die "\n[ERROR]: No FASTA files found in genome directory. Sequences must be over 50kb and less than 50% gaps. Exiting.\n" 
	if @$fasfiles == 0;

     # create data directories for reverse strand
    my $sequtils = Tephra::NonLTR::SeqUtils->new;
    my $revfasfiles = $sequtils->invert_seq($genome_dir, $minus_dna_dir);
    die "\n[ERROR]: No FASTA files found in genome directory. Sequences must be over 50kb and less than 50% gaps. Exiting.\n"
        if @$revfasfiles == 0;

    my %data_dirs = ( plus  => { files => $fasfiles,    outdir => $plus_out_dir  }, 
		      minus => { files => $revfasfiles, outdir => $minus_out_dir } );

    #dd \%data_dirs and exit;
    my ($is_pos_success, $is_rev_success) = (0, 0);
    $pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
			      
			      for my $key (%$data_ref) {
				  $is_pos_success = $data_ref->{$key} if $key eq 'plus';
				  $is_rev_success = $data_ref->{$key} if $key eq 'minus';
			      }
			      #say STDERR "Finished running HMM search on $$data_ref..." if $self->verbose; 
			      #my $t1 = gettimeofday();
			      #my $elapsed = $t1 - $t0;
			      #my $time = sprintf("%.2f",$elapsed/60);
			      #say $fmlog "$ident just finished with PID $pid and exit code: $exit_code in $time minutes";
			} );

    my %results;
    for my $strand (keys %data_dirs) {
	$pm->start($strand) and next;
	$SIG{INT} = sub { $pm->finish };
	
	my $pp_result = $self->run_model_search_and_postprocess( $strand, 
								 $data_dirs{$strand}{files}, 
								 $genome_dir, 
								 $data_dirs{$strand}{outdir}, 
								 $phmm_dir );
	$results{$strand} = $pp_result;

	$pm->finish(0, \%results); #\$strand);
    }

    $pm->wait_all_children;

    ##say STDERR "Running forward..." if $self->verbose;
    if ($is_pos_success && $is_rev_success) {
	# Validation for Q value
	my $pp2 = Tephra::NonLTR::QValidation->new( 
	    outdir  => $main_data_dir, 
	    phmmdir => $phmm_dir, 
	    fasta   => $genome_dir,
	    threads => $self->threads );
	$pp2->validate_q_score;
	
	return ($genome_dir, $main_data_dir);
    }
    else {
	remove_tree( $main_data_dir, { safe => 1} );                                                         
        remove_tree( $genome_dir, { safe => 1} );                                                                                            
        remove_tree( $minus_out_dir, { safe => 1} );                                                                                         
	remove_tree( $minus_dna_dir, { safe => 1} );

        return (undef, undef);                                                                                                                
    }
}

sub run_model_search_and_postprocess {
    my $self = shift;
    my ($strand, $fasfiles, $genome_dir, $out_dir, $phmm_dir) = @_;
    my $threads = $self->threads;
    my $program_dir = $self->pdir;

    # adjust requested threads since we are running analysis on both strands
    my $thr;
    if ($threads % 2 == 0) {
	$thr = sprintf("%.0f",$threads/2);
    }
    elsif (+($threads-1) % 2 == 0) {
	$thr = sprintf("%.0f",$threads-1/2);
    }
    else {
	$thr = 1;
    }

    my $reverse = defined $strand && $strand eq 'plus' ? 0 : 1;

    #say "\nDEBUG: fasfiles\n";
    #dd $fasfiles;
    # run HMM model search
    for my $file (sort @$fasfiles) {
        my $run_hmm = Tephra::NonLTR::RunHMM->new(
            fasta    => $file,
            fastadir => $genome_dir,
            outdir   => $out_dir,
            phmmdir  => $phmm_dir,
            pdir     => $program_dir,
            strand   => $strand,
            threads  => $threads,
            verbose  => $self->verbose );
        $run_hmm->run_mgescan;
    }

    # run postprocessing
    my $pp = Tephra::NonLTR::Postprocess->new(
        fastadir => $genome_dir,
        outdir   => $out_dir,
        reverse  => $reverse );
    my $pp_result = $pp->postprocess;

    return $pp_result;
}

sub _split_genome {
    my $self = shift;
    my ($genome, $genome_dir) = @_;

    my $length_thresh = 1e4; # 10kb
    my $nperc_thresh  = 50;  # 50%

    my @fasfiles;
    my $kseq = Bio::DB::HTS::Kseq->new($genome);
    my $iter = $kseq->iterator;

    while (my $seqobj = $iter->next_seq) {
	my $id = $seqobj->name;
	my $seq = $seqobj->seq;
	# filter by length (over 10kb) and N-percent (reject over 50%) to speed up search
	# and reduce spurious matches
	my $seqlength = length($seq);
	my $n_count = ($seq =~ tr/Nn//);
	my $n_perc = sprintf("%.2f",($n_count/$seqlength)*100);
	if ($seqlength >= $length_thresh && $n_perc <= $nperc_thresh) {
	    my $outfile = File::Spec->catfile($genome_dir, $id.'.fasta');
	    open my $out, '>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
	    say $out join "\n", ">".$id, $seq;
	    close $out;
	    push @fasfiles, $outfile;
	}
    }

    return \@fasfiles;
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

    perldoc Tephra::NonLTR::NonLTRSearch


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
