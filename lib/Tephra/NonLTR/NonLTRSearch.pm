package Tephra::NonLTR::NonLTRSearch;

use 5.010;
use Moose;
use Cwd;
use File::Spec;
use File::Find;
use File::Path qw(make_path);
use Tephra::NonLTR::RunHMM;
use Tephra::NonLTR::Postprocess;
use Tephra::NonLTR::QValidation;
use Tephra::NonLTR::SeqUtils;
use namespace::autoclean;

has fastadir => ( is => 'ro', isa => 'Path::Class::File', required => 1, coerce => 1 );
has outdir   => ( is => 'ro', isa => 'Path::Class::Dir',  required => 1, coerce => 1 );
has pdir     => ( is => 'ro', isa => 'Path::Class::Dir',  required => 0, coerce => 1 );

sub find_nonltrs {
    my $self = shift;
    my $main_data_dir = $self->outdir;
    my $genome_dir    = $self->fastadir;
    my $program_dir   = $self->pdir;
    $program_dir //= File::Spec->catdir($ENV{HOME}, '.tephra');

    my $phmm_dir      = File::Spec->catdir($program_dir, 'pHMM');
    my $plus_out_dir  = File::Spec->catdir($main_data_dir, 'f');
    my $minus_out_dir = File::Spec->catdir($main_data_dir, 'b');
    my $minus_dna_dir = $genome_dir."_b";

    unless (-e $minus_dna_dir) {
	make_path( $minus_dna_dir, {verbose => 0, mode => 0771,} );
    }

    # Forward strand
    my @fasfiles;
    find( sub { push @fasfiles, $File::Find::name if -f and /\.fa.*?$/ }, $genome_dir );
    die "\nERROR: No FASTA files found in genome directory. Exiting.\n" if @fasfiles == 0;

    printf STDERR "Running forward...\n";
    for my $file (sort @fasfiles) {    
	my $run_hmm = Tephra::NonLTR::RunHMM->new( 
	    fasta   => $file, 
	    outdir  => $plus_out_dir, 
	    phmmdir => $phmm_dir, 
	    pdir    => $program_dir );
	$run_hmm->run_mgescan;
    }

    my $pp = Tephra::NonLTR::Postprocess->new( 
	fastadir => $genome_dir, 
	outdir   => $plus_out_dir, 
	reverse  => 0 );
    $pp->postprocess;

    # Backward strand
    printf "Running backward...\n";

    my $sequtils = Tephra::NonLTR::SeqUtils->new;
    $sequtils->invert_seq($genome_dir, $minus_dna_dir);

    my @revfasfiles;
    find( sub { push @revfasfiles, $File::Find::name if -f and /\.fa.*$/ }, $minus_dna_dir );
    
    for my $file (sort @revfasfiles) {
	my $run_rev_hmm = Tephra::NonLTR::RunHMM->new( 
	    fasta   => $file, 
	    outdir  => $minus_out_dir, 
	    phmmdir => $phmm_dir, 
	    pdir    => $program_dir );

	$run_rev_hmm->run_mgescan;
    }

    my $pp_rev = Tephra::NonLTR::Postprocess->new( 
	fastadir => $minus_dna_dir, 
	outdir   => $minus_out_dir, 
	reverse  => 1 );

    $pp_rev->postprocess;
    
    # Validation for Q value
    my $pp2 = Tephra::NonLTR::QValidation->new( 
	outdir  => $main_data_dir, 
	phmmdir => $phmm_dir, 
	fasta   => $genome_dir );

    $pp2->validate_q_score;
}

__PACKAGE__->meta->make_immutable;

1;
