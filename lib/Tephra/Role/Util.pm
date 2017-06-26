package Tephra::Role::Util;

use 5.014;
use Moose::Role;
use Bio::DB::HTS::Faidx;
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use Try::Tiny;
use Log::Log4perl;
use Log::Any::Adapter;
use namespace::autoclean;

=head1 NAME

Tephra::Role::Util - Helper methods for running programs

=head1 VERSION

Version 0.08.2

=cut

our $VERSION = '0.08.2';
$VERSION = eval $VERSION;

sub run_cmd {
    my $self = shift;
    my ($cmd) = @_;

    try {
	system([0..5], $cmd);
    }
    catch {
	die "\nERROR: $cmd failed. Here is the exception: $_\n";
    };
}

sub capture_cmd {
    my $self = shift;
    my ($cmd) = @_;

    try {
	my @out = capture { system([0..5], $cmd) };
    }
    catch {
	my $err = $_;
	if ($err =~ /SEGV/) {
	    return 'failed';
	}
	die "\nERROR: $cmd failed. Here is the exception: $_\n";
    };
}

sub index_ref {
    my $self = shift;
    my ($fasta) = @_;

    my $index = Bio::DB::HTS::Faidx->new($fasta);
    return $index;
}

sub get_tephra_logger {
    my $self = shift;
    my ($logfile) = @_;

    my $conf = qq{
    log4perl.category.Tephra      = INFO, Logfile, Screen
    log4perl.appender.Logfile          = Log::Log4perl::Appender::File
    log4perl.appender.Logfile.filename = $logfile
    log4perl.appender.Logfile.layout   = Log::Log4perl::Layout::PatternLayout
    log4perl.appender.Logfile.layout.ConversionPattern = %m%n
    log4perl.appender.Screen         = Log::Log4perl::Appender::Screen
    log4perl.appender.Screen.stderr  = 1
    log4perl.appender.Screen.layout  = Log::Log4perl::Layout::SimpleLayout
    };

    Log::Log4perl::init( \$conf );
    Log::Any::Adapter->set('Log4perl');
    my $log = Log::Any->get_logger( category => "Tephra" );

    return $log;
}
sub get_SO_terms {
    my $self = shift;

    my %table = (
        'LTR_retrotransposon'     => 'SO:0000186',
        'non_LTR_retrotransposon' => 'SO:0000189',
        
        'U_box'                => 'SO:0001788',
        'RR_tract'             => 'SO:0000435',
        'long_terminal_repeat' => 'SO:0000286',
        'inverted_repeat'      => 'SO:0000294',
        'primer_binding_site'  => 'SO:0005850',
        'protein_match'        => 'SO:0000349',
        
        'terminal_inverted_repeat_element' => 'SO:0000208',
        'terminal_inverted_repeat'         => 'SO:0000481',
        'helitron'                         => 'SO:0000544',
        'MITE'                             => 'SO:0000338',
        'DNA_transposon'                   => 'SO:0000182' );

    return \%table;
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

    perldoc Tephra::Role::Util


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1;
