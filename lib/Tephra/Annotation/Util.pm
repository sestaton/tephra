package Tephra::Annotation::Util;

use 5.014;
use Moose;
use namespace::autoclean;

=head1 NAME

Tephra::Annotation::Util - Utility methods for working with transposon annotations

=head1 VERSION

Version 0.08.1

=cut

our $VERSION = '0.08.1';
$VERSION = eval $VERSION;

has 'debug' => (
    is         => 'ro',
    isa        => 'Bool',
    predicate  => 'has_debug',
    lazy       => 1,
    default    => 0,
);

=head2 map_superfamily_name

 Title   : map_superfamily_name

 Usage   : my $superfamily_name = $self->map_superfamily_name($match);
           
 Function: Get the superfamily common name based on the 3-letter code
                                                                            Return_type
 Returns : The TE superfamily name                                          Scalar
            
                                                                            Arg_type
 Args    : A BLAST hit or any sequence identifier                           Scalar

=cut

## NB: This method is pulled straight from Transposome (github.com/sestaton/Transposome).
sub map_superfamily_name {
    my $self = shift;
    my ($id) = @_;

    my ($sfamily_code) = ($id =~ /(^[A-Z]{3})_?/);
    unless (defined $sfamily_code) {
	say STDERR "\n[WARNING]: Could not get 3-letter code from: $id. Skipping.\n"
	    if $self->debug;
	return 0;
    }

    my %sfcode_table = (
        'RYD' => 'DIRS',
        'RYN' => 'Ngaro',
        'RYX' => 'Unknown_DIRS',
        'RYV' => 'VIPER',
        'RII' => 'I',
        'RIJ' => 'Jockey',
        'RIL' => 'L1',
        'RIR' => 'R2',
        'RIT' => 'RTE',
        'RIX' => 'Unknown_LINE',
        'RLB' => 'Bel/Pao',
        'RLG' => 'Gypsy',
        'RLC' => 'Copia',
        'RLE' => 'ERV',
        'RLR' => 'Retrovirus',
        'RLX' => 'Unknown_LTR',
        'PPP' => 'Penelope',
        'RPX' => 'Unknown_PLE',
        'RSS' => '5S',
        'RSL' => '7SL',
        'RST' => 'tRNA',
        'RSX' => 'Unknown_SINE',
        'RXX' => 'Unknown_retrotransposon',
        'DYC' => 'Crypton',
        'DYX' => 'Unknown_Crypton',
        'DTC' => 'CACTA',
        'DTA' => 'hAT',
        'DTE' => 'Merlin',
        'DTM' => 'Mutator',
        'DTP' => 'P',
        'DTH' => 'PIF/Harbinger',
        'DTB' => 'PiggyBac',
        'DTT' => 'Tc1/Mariner',
        'DTR' => 'Transib',
        'DTX' => 'Unknown_TIR',
        'DXX' => 'Unknown_DNA_transposon',
        'DHH' => 'Helitron',
        'DHX' => 'Unknown_Helitron',
        'DMM' => 'Maverick',
        'DMX' => 'Unknown_Maverick',
	'RST' => 'SINE2/tRNA' );
    
    if (exists $sfcode_table{$sfamily_code}) {
	return $sfcode_table{$sfamily_code};
    }
    else {
	say STDERR "\n[WARNING]: No 3-letter code could be found for: $sfamily_code\n"
	    if $self->debug;
	return 0;
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

    perldoc Tephra::Annotation::Util


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
