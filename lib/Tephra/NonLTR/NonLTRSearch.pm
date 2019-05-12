package Tephra::NonLTR::NonLTRSearch;

use 5.014;
use Moose;
use File::Spec;
use File::Find;
use File::Basename;
use File::Path qw(make_path remove_tree);
use Cwd        qw(abs_path);
use Bio::DB::HTS::Kseq;
use Tephra::NonLTR::RunHMM;
use Tephra::NonLTR::Postprocess;
use Tephra::NonLTR::QValidation;
use Tephra::NonLTR::SeqUtils;
use Tephra::Config::Exe;
use namespace::autoclean;

=head1 NAME

Tephra::NonLTR::NonLTRSearch - Search a genome for non-LTR retrotransposons

=head1 VERSION

Version 0.12.3

=cut

our $VERSION = '0.12.3';
$VERSION = eval $VERSION;

has genome  => ( is => 'ro', isa => 'Maybe[Str]', required => 1 );
has outdir  => ( is => 'ro', isa => 'Maybe[Str]', required => 0 );
has pdir    => ( is => 'ro', isa => 'Maybe[Str]', required => 0 );
has verbose => ( is => 'ro', isa => 'Bool', predicate  => 'has_debug', lazy => 1, default => 0 );

sub find_nonltrs {
    my $self = shift;
    my $main_data_dir = $self->outdir;
    my $program_dir   = $self->pdir;
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

    # Forward strand
    my $fasfiles = $self->_split_genome($genome, $genome_dir);
    die "\n[ERROR]: No FASTA files found in genome directory. Sequences must be over 50kb and less than 50% gaps. Exiting.\n" 
	if @$fasfiles == 0;

    say STDERR "Running forward..." if $self->verbose;
    for my $file (sort @$fasfiles) {    
	my $run_hmm = Tephra::NonLTR::RunHMM->new( 
	    fasta    => $file, 
	    fastadir => $genome_dir,
	    outdir   => $plus_out_dir, 
	    phmmdir  => $phmm_dir, 
	    pdir     => $program_dir,
	    strand   => 'plus',
	    verbose  => $self->verbose );
	$run_hmm->run_mgescan;
    }

    my $pp = Tephra::NonLTR::Postprocess->new( 
	fastadir => $genome_dir, 
	outdir   => $plus_out_dir, 
	reverse  => 0 );
    my $fpp_result = $pp->postprocess;

    unless ($fpp_result) {
	say STDERR "\n[WARNING]: No non-LTR elements were found on the forward strand. Will search reverse strand.\n";
    }

    # Backward strand
    say STDERR "Running backward..." if $self->verbose;

    my $sequtils = Tephra::NonLTR::SeqUtils->new;
    my $revfasfiles = $sequtils->invert_seq($genome_dir, $minus_dna_dir);
    die "\n[ERROR]: No FASTA files found in genome directory. Sequences must be over 50kb and less than 50% gaps. Exiting.\n"
        if @$revfasfiles == 0;
    
    for my $file (sort @$revfasfiles) {
	my $run_rev_hmm = Tephra::NonLTR::RunHMM->new( 
	    fasta    => $file, 
	    fastadir => $minus_dna_dir,
	    outdir   => $minus_out_dir, 
	    phmmdir  => $phmm_dir, 
	    pdir     => $program_dir,
	    strand   => 'minus',
	    verbose  => $self->verbose );
	$run_rev_hmm->run_mgescan;
    }

    my $pp_rev = Tephra::NonLTR::Postprocess->new( 
	fastadir => $minus_dna_dir, 
	outdir   => $minus_out_dir, 
	reverse  => 1 );
    my $bpp_result = $pp_rev->postprocess;

    unless ($bpp_result) {
        say STDERR "\n[WARNING]: No non-LTR elements were found on the reverse strand.\n";
    }
 
    if (!$fpp_result && !$bpp_result) { 
	remove_tree( $main_data_dir, { safe => 1} );
	remove_tree( $genome_dir, { safe => 1} );
	remove_tree( $minus_dna_dir, { safe => 1} );

	return (undef, undef);
    }
    else {
	# Validation for Q value
	my $pp2 = Tephra::NonLTR::QValidation->new( 
	    outdir  => $main_data_dir, 
	    phmmdir => $phmm_dir, 
	    fasta   => $genome_dir );
	$pp2->validate_q_score;
	
	return ($genome_dir, $main_data_dir);
    }
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
