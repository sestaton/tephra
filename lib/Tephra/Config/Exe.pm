package Tephra::Config::Exe;

use 5.010;
use Moose;
use MooseX::Types::Path::Class;
use File::Spec;
use File::Path qw(make_path remove_tree);
use File::Basename;
use Log::Any qw($log);
use namespace::autoclean;

=head1 NAME

Tephra::Config::Exe - Class for setting up PATHs for Tephra dependencies

=head1 VERSION

Version 0.03.4

=cut

our $VERSION = '0.03.4';

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
    my $root = $self->basedir;

    unless (-e $root) {
        make_path($root, {verbose => 0, mode => 0711,});
    }

    # we don't want to reconfigure every time the tests run 
    my $gt       = File::Spec->catfile($root,   'gt', 'bin', 'gt');
    my $hscan    = File::Spec->catfile($root,   'helitronscanner', 'HelitronScanner', 'HelitronScanner.jar');
    my $hmm2bin  = File::Spec->catdir($root,    'hmmer-2.3.2', 'bin');
    my $hmm3bin  = File::Spec->catdir($root,    'hmmer-3.1b2-linux-intel-x86_64', 'binaries');
    my $moddir   = File::Spec->catdir($root,    'pHMM');
    my $chrdir   = File::Spec->catdir($root,    'hmm');
    my $mgescan  = File::Spec->catfile($chrdir, 'tephra-MGEScan');
    my $transla  = File::Spec->catfile($chrdir, 'tephra-translate');
    my $pamlbin  = File::Spec->catdir($root,    'paml4.8', 'bin');
    my $transeq  = File::Spec->catdir($root,    'EMBOSS-6.5.7', 'bin', 'transeq');
    my $blastph  = File::Spec->catdir($root,    'ncbi-blast-2.3.0+', 'bin');
    my $htsdir   = File::Spec->catdir($root,    'htslib-1.3.1', 'htslib');

    # this is to avoid building each time
    my @path = split /:|;/, $ENV{PATH};    
    for my $p (@path) {
	my $texe  = File::Spec->catfile($p, 'transeq');
	my $bexe  = File::Spec->catfile($p, 'blastn');
	if (-e $texe && -x $texe) {
	    $transeq = $texe;
	}
	if (-e $bexe && -x $bexe) {
	    $blastph = $p;
	}
    }

    return ({
        gt         => $gt,
        hscanjar   => $hscan,
        hmmer2bin  => $hmm2bin,
	hmmer3bin  => $hmm3bin,
        modeldir   => $moddir,
        hmmdir     => $chrdir,
        mgescan    => $mgescan,
        transcmd   => $transla,
        pamlbin    => $pamlbin,
        transeq    => $transeq,
        blastpath  => $blastph,
	htslibdir  => $htsdir });
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

    perldoc Tephra::Config::Exe


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut 

__PACKAGE__->meta->make_immutable;

1;
