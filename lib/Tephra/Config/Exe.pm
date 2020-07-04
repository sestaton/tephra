package Tephra::Config::Exe;

use 5.014;
use Moose;
use MooseX::Types::Path::Class;
use File::Spec;
use File::Path qw(make_path);
use Log::Any   qw($log);
use File::Basename;
use namespace::autoclean;

=head1 NAME

Tephra::Config::Exe - Class for setting up PATHs for Tephra dependencies

=head1 VERSION

Version 0.13.1

=cut

our $VERSION = '0.13.1';

has basedir => (
    is       => 'ro',
    isa      => 'Path::Class::Dir',
    required => 0,
    coerce   => 1,
    default  => sub {
	return $ENV{TEPHRA_DIR} // Path::Class::Dir->new($ENV{HOME}, '.tephra')
    },
);

sub get_config_paths {
    my $self = shift;
    my $root = $self->basedir; #->absolute->resolve;

    my $bindir = File::Spec->catdir($root, 'bin');
    unless (-e $bindir) {
        make_path($bindir, {verbose => 0, mode => 0711,});
    }

    ## Binaries that are part of, or used by, Tephra
    my $baseml   = File::Spec->catfile($bindir, 'baseml');
    my $vmatch   = File::Spec->catfile($bindir, 'vmatch');
    my $mkvtree  = File::Spec->catfile($bindir, 'mkvtree');
    my $cleanpp  = File::Spec->catfile($bindir, 'cleanpp.sh');
    my $gt       = File::Spec->catfile($bindir, 'gt');
    my $musbin   = File::Spec->catfile($bindir, 'muscle');
    my $transeq  = File::Spec->catfile($bindir, 'transeq');
    my $blastn   = File::Spec->catfile($bindir, 'blastn');
    my $mblastdb = File::Spec->catfile($bindir, 'makeblastdb');
    my $mgescan  = File::Spec->catfile($bindir, 'tephra-MGEScan');
    my $transla  = File::Spec->catfile($bindir, 'tephra-translate');

    ## Required libraries and deps that cannot be easily dropped into the same bin directory
    my $gtdata  = File::Spec->catdir($root,  'gtdata');
    my $htsdir  = File::Spec->catdir($root,  'htslib-1.3.1', 'htslib');
    my $hmm2bin = File::Spec->catdir($root,  'hmmer-2.3.2', 'bin');
    my $hmm3bin = File::Spec->catdir($root,  'hmmer-3.1b2', 'binaries');
    my $hscan   = File::Spec->catfile($root, 'helitronscanner', 'HelitronScanner', 'HelitronScanner.jar');

    ## Databases and models used by Tephra 
    my $tephradb = File::Spec->catfile($root,     'TephraDB');
    my $trnadb   = File::Spec->catfile($tephradb, 'eukaryotic-tRNAs.fa');
    my $hmmdb    = File::Spec->catfile($tephradb, 'transposable+element.hmm');
    my $chrhmm   = File::Spec->catfile($tephradb, 'mgescan-chr.hmm');
    my $moddir   = File::Spec->catdir($root,      'pHMM');

    return ({
	tephrabin   => $bindir,
	tephradb    => $tephradb,
        gt          => $gt,
	gtdata      => $gtdata,
	vmatch      => $vmatch,
	mkvtree     => $mkvtree,
	cleanpp     => $cleanpp,
        hscanjar    => $hscan,
        hmmer2bin   => $hmm2bin,
	hmmer3bin   => $hmm3bin,
        modeldir    => $moddir,
	trnadb      => $trnadb,
	hmmdb       => $hmmdb,
        hmmdir      => $tephradb,
	chrhmm      => $chrhmm,
        mgescan     => $mgescan,
        transcmd    => $transla,
        baseml      => $baseml,
        transeq     => $transeq,
        blastn      => $blastn,
	makeblastdb => $mblastdb,
	htslibdir   => $htsdir,
	muscle      => $musbin });
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

    perldoc Tephra::Config::Exe


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut 

__PACKAGE__->meta->make_immutable;

1;
