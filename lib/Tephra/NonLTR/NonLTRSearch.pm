package Tephra::NonLTR::NonLTRSearch;

use 5.010;
use Moose;
use Cwd;
use File::Spec;
use File::Find;
use File::Basename;
use File::Path qw(make_path);
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

Version 0.03.0

=cut

our $VERSION = '0.03.0';
$VERSION = eval $VERSION;

has genome => ( is => 'ro', isa => 'Maybe[Str]', required => 1 );
has outdir => ( is => 'ro', isa => 'Maybe[Str]',  required => 0 );
has pdir   => ( is => 'ro', isa => 'Maybe[Str]',  required => 0 );

sub find_nonltrs {
    my $self = shift;
    my $main_data_dir = $self->outdir;
    my $program_dir   = $self->pdir;
    my $genome = $self->genome;
    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($phmm_dir) = @{$config}{qw(modeldir)};

    my ($gname, $gpath, $gsuffix) = fileparse($genome, qr/\.[^.]*/);
    my $genome_dir = File::Spec->catdir($gpath, $gname.'_genome');
    $main_data_dir //= File::Spec->catdir($gpath, $gname.'_nonLTRs');

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
    $self->_split_genome($genome, $genome_dir);
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

    return ($genome_dir, $main_data_dir);
}

sub _split_genome {
    my $self = shift;
    my ($genome, $genome_dir) = @_;

    my $kseq = Bio::DB::HTS::Kseq->new($genome);
    my $iter = $kseq->iterator;
    while (my $seqobj = $iter->next_seq) {
	my $id = $seqobj->name;
	my $outfile = File::Spec->catfile($genome_dir, $id.'.fasta');
	open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";
	say $out join "\n", ">".$id, $seqobj->seq;
	close $out;
    }
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

    perldoc Tephra::NonLTR::NonLTRSearch


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
