package Tephra::NonLTR::RunHMM;

use 5.014;
use Moose;
use MooseX::Types::Path::Class;
use autodie;
use File::Path          qw(make_path);
use IPC::System::Simple qw(capture EXIT_ANY);
use Cwd                 qw(abs_path);
use File::Basename;
use File::Spec;
use File::Find;
use Bio::SearchIO;
use Try::Tiny;
use Tephra::NonLTR::SeqUtils;
use Tephra::Config::Exe;
use namespace::autoclean;

with 'Tephra::NonLTR::Role::PathUtils';

=head1 NAME

Tephra::NonLTR::RunHMM - Search for non-LTR coding domains (adapted from MGEScan-nonLTR)

=head1 VERSION

Version 0.12.3

=cut

our $VERSION = '0.12.3';
$VERSION = eval $VERSION;

has fasta    => ( is => 'ro', isa => 'Path::Class::File', required => 1, coerce => 1 );
has fastadir => ( is => 'ro', isa => 'Path::Class::File', required => 1, coerce  => 1 );
has outdir   => ( is => 'ro', isa => 'Path::Class::Dir',  required => 1, coerce => 1 );
has phmmdir  => ( is => 'ro', isa => 'Path::Class::Dir',  required => 1, coerce => 1 );
has pdir     => ( is => 'ro', isa => 'Path::Class::Dir',  required => 1, coerce => 1 );

has verbose  => ( is => 'ro', isa => 'Bool', predicate  => 'has_verbose', lazy => 1, default => 0 );
has strand   => ( is => 'ro', isa => 'Str',  required => 1, default  => 'plus' );

sub run_mgescan {
    my $self = shift;
    my $dna_file = $self->fasta->absolute->resolve;
    my $fas_dir  = $self->fastadir;
    my $out_dir  = $self->outdir;
    my $phmm_dir = $self->phmmdir->absolute->resolve;
    my $pdir     = $self->pdir->absolute->resolve;
    my $strand   = $self->strand;

    my $sequtils = Tephra::NonLTR::SeqUtils->new( verbose => $self->verbose );
    my ($dna_name, $dna_path, $dna_suffix) = fileparse($dna_file, qr/\.[^.]*/);
    my $outf_dir = File::Spec->catdir($out_dir, 'out1');
    my $pos_dir  = File::Spec->catdir($out_dir, 'pos');
    unless ( -d $outf_dir ) {
	make_path( $outf_dir, {verbose => 0, mode => 0771,} );
    }
    
    unless ( -d $pos_dir ) {
	make_path( $pos_dir, {verbose => 0, mode => 0771,} );
    }
    
    # get signal for some state of ORF1, RT, and APE
    say STDERR "Getting signal..." if $self->verbose;
    say STDERR "    Protein sequence..." if $self->verbose;
    my $pep_file = File::Spec->catfile($out_dir, $dna_name.$dna_suffix.'.pep');
    my $result = $sequtils->translate($fas_dir, $dna_file, $pep_file, $strand);
    unless ($result) {
	say STDERR "\n[ERROR]: 'transeq' failed. Please check input.\n";
	return;
    }

    say STDERR "    RT signal..." if $self->verbose;
    my $phmm_file = File::Spec->catfile($phmm_dir, 'ebi_ds36752_seq.hmm');
    my $domain_rt_pos_file = File::Spec->catfile($pos_dir, $dna_name.$dna_suffix.'.rt.pos');
    $self->get_signal_domain($pep_file, $phmm_file, $domain_rt_pos_file);
    
    say STDERR "    APE signal..." if $self->verbose;
    $phmm_file = File::Spec->catfile($phmm_dir, 'ebi_ds36736_seq.hmm');
    my $domain_ape_pos_file = File::Spec->catfile($pos_dir, $dna_name.$dna_suffix.'.ape.pos');
    $self->get_signal_domain($pep_file, $phmm_file, $domain_ape_pos_file);
    
    # generate corresponsing empty domains files if either of them does not exist 
    if (-e $domain_rt_pos_file || -e $domain_ape_pos_file ) {
	#print $dna_name."\n" if $self->verbose;	
	if (! -e $domain_rt_pos_file){
	    open my $out, '>', $domain_rt_pos_file or die "\n[ERROR]: Could not open file: $domain_rt_pos_file\n";
	    print $out "";
	    close $out;
	}
	elsif (! -e $domain_ape_pos_file){
	    open my $out, '>', $domain_ape_pos_file or die "\n[ERROR]: Could not open file: $domain_ape_pos_file\n";
	    print $out "";
	    close $out;
	}

	# run hmm
	say STDERR "Running HMM..." if $self->verbose;

	my $mgescan  = File::Spec->catfile($pdir, 'hmm', 'tephra-MGEScan');
	my $out_file = File::Spec->catfile($outf_dir, $dna_name.$dna_suffix);
	my $chrhmm   = File::Spec->catfile($pdir, 'hmm', 'chr.hmm');
	my $ldir = $pdir.'/';
	$outf_dir .= '/';
	my $tephra_dir = $ENV{TEPHRA_DIR} // File::Spec->catfile($ENV{HOME}, '.tephra');
	$ENV{PATH} = join ':', $ENV{PATH}, File::Spec->catfile($tephra_dir, 'EMBOSS-6.5.7', 'bin');
	my $cmd = "$mgescan -m $chrhmm -s $dna_file -r $domain_rt_pos_file -a $domain_ape_pos_file -o $out_file -p $ldir -d $outf_dir";
	say STDERR "CMD: $cmd" if $self->verbose;
	system($cmd);
	#exit;
    }

    # not implemented
    #unless ($self->debug) {
    unlink $pep_file if -e $pep_file;
    #}

    return;
}

sub get_signal_domain {
    my $self = shift;
    my ($pep_file, $phmm_file, $domain_rt_pos_file) = @_;

    # debugging
    my ($pname, $ppath, $psuffix) = fileparse($phmm_file, qr/\.[^.]*/);
    my $signal_out = $pep_file.'_'.$pname.'_signal_searchout.txt';

    my %domain_start;
    my %domain_end;
    my %domain_pos;
    my $evalue;
    my $temp_file   = $domain_rt_pos_file.'temp';
    my $stemp_file  = $domain_rt_pos_file.'temp_sorted';
    my $output_file = $domain_rt_pos_file;

    # run hmmsearch to find the domain and save it in the temprary file
    my $hmmsearch   = $self->find_hmmsearch;
    my @hmm_results = capture([0..5], $hmmsearch, '-E', '0.00001', $phmm_file, $pep_file);
    $self->_parse_hmmsearch(\@hmm_results, $signal_out, $temp_file);

    if (-s $temp_file) {	
	$self->_sort_matches($temp_file, $stemp_file);
        my ($start, $end) = (-1, -1);
        my @pre = (-1000, -1000, -1000, -1000, -1000, -1000);
        open my $in, '<', $stemp_file or die "\n[ERROR]: Could not open file: $stemp_file\n";
        open my $out, '>', $output_file or die "\n[ERROR]: Could not open file: $output_file\n";

	## NB hard thresholds: length = 300, evalue = 1e-5 
        while (my $each_line = <$in>) {
	    chomp $each_line;
            my @temp = split /\t/, $each_line;
	    # seq_start seq_end domain_start domain_end score evalue
            if ($temp[0] - $pre[1] < 300) {
                $end = $temp[1];
                $evalue = $evalue * $temp[5];
            }
	    else {
                if ($start >= 0 && $evalue <= 0.00001) {
                    say $out join "\t", $start, $end, @pre[4..5];
                }
                ($start, $end, $evalue) = @temp[0,1,5];
            }
            @pre = @temp;
        }

        if ($start >= 0 && $evalue <= 0.00001) { 
            say $out join "\t", $start, $end, @pre[4..5];
        }
        close $in;
        close $out;
	unlink $stemp_file;
    }
    unlink $temp_file;

    return;
}

sub _parse_hmmsearch {
    my $self = shift;
    my ($hmm_results, $signal_out, $outfile) = @_;

    open my $o, '>', $signal_out or die "\n[ERROR]: Could not open file: $signal_out\n";
    print $o @$hmm_results;
    close $o;

    if (-s $signal_out) {
	open my $out, '>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
	my $hmmer_in = Bio::SearchIO->new(-file => $signal_out, -format => 'hmmer');
	
	my @evalues;
	while ( my $result = $hmmer_in->next_result ) {    
	    while ( my $hit = $result->next_hit ) {
		while ( my $hsp = $hit->next_hsp ) {
		    my $hstart = $hsp->start('hit');
		    my $hstop  = $hsp->end('hit');
		    my $qstart = $hsp->start('query');
		    my $qstop  = $hsp->end('query');
		    my $score  = $hsp->score;
		    my $e_val  = $hsp->evalue;
		    say $out join "\t", eval($hstart*3), eval($hstop*3), $qstart, $qstop, $score, $e_val;
		}
	    }
	}
	close $out;
    }
    unlink $signal_out;

    return;
}

sub _sort_matches {
    my $self = shift;
    my ($unsorted, $sorted) = @_;
    my %hash;

    ##NB: This method returns the parsed HMM matches in order by coordinate. 
    # Input:
    #46927035 46928451 -87.4 39e-08
    #84073929 84075108 -115.7 6.1e-07
    #217494303 217495593 -75.2 1.2e-08
    #
    # Output:
    #217494303 217495593 -75.2 1.2e-08
    #46927035 46928451 -87.4 3.9e-08
    #84073929 84075108 -115.7 6.1e-07

    open my $in, '<', $unsorted or die "\n[ERROR]: Could not open file: $unsorted\n";
    open my $out, '>', $sorted or die "\n[ERROR]: Could not open file: $sorted\n";

    while (my $l = <$in>) {
	chomp $l;
	my @f = split /\t/, $l;
	my $start = shift @f;
	$hash{$start} = join "||", @f;
    }
    close $in;

    for my $sk (sort { $a <=> $b } keys %hash) {
	my @doms = split /\|\|/, $hash{$sk};
	say $out join "\t", $sk, @doms;
    }
    close $out;

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

    perldoc Tephra::NonLTR::RunHMM


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
