package Tephra::Role::Logger;

use 5.014;
use Moose::Role;
use Log::Any;
use Log::Any::Adapter;
use Log::Log4perl;
use DateTime;
use POSIX       qw(strftime);
use Time::HiRes qw(gettimeofday);
use Lingua::EN::Inflect;
use namespace::autoclean;

=head1 NAME

Tephra::Role::Logger - Helper methods for logging Tephra commands

=head1 VERSION

Version 0.11.0

=cut

our $VERSION = '0.11.0';
#$VERSION = eval $VERSION;

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

sub init_tephra {
    my $self = shift;
    my ($config) = @_;

    #my $log_file = File::Spec->catfile($config->{output_directory}, $config->{run_log_file});
    my $conf = qq{
    log4perl.category.Tephra      = INFO, Logfile, Screen
    log4perl.appender.Logfile          = Log::Log4perl::Appender::File
    log4perl.appender.Logfile.filename = $config->{logfile}
    log4perl.appender.Logfile.layout   = Log::Log4perl::Layout::PatternLayout
    log4perl.appender.Logfile.layout.ConversionPattern = %m%n
    log4perl.appender.Screen         = Log::Log4perl::Appender::Screen
    log4perl.appender.Screen.stderr  = 1
    log4perl.appender.Screen.layout  = Log::Log4perl::Layout::SimpleLayout
    };

    Log::Log4perl::init( \$conf );
    Log::Any::Adapter->set('Log4perl');
    my $log = Log::Any->get_logger( category => "Tephra" );
    
    my $t0 = [Time::HiRes::gettimeofday()];
    my $ts = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("======== Tephra version: $VERSION (started at: $ts) ========");
    $log->info("Configuration - Log file for monitoring progress and errors: $config->{logfile}");
    $log->info("Configuration - Genome file:                                 $config->{genome}");
    $log->info("Configuration - Repeat database:                             $config->{repeatdb}");
    $log->info("Configuration - Number of threads:                           $config->{threads}");

    return ($t0, $log);
}

sub log_interval {
    my $self = shift;
    my ($t0, $log) = @_;
    
    #load_classes('DateTime', 'Time::HiRes', 'Lingua::EN::Inflect', 'POSIX');

    my $t1    = [Time::HiRes::gettimeofday()];
    my $t0_t1 = Time::HiRes::tv_interval($t0, $t1);
    my $dt    = DateTime->from_epoch( epoch => 0 );

    $dt = $dt->add( seconds => $t0_t1 );
    $dt = $dt - DateTime->from_epoch( epoch => 0 );
    
    my @time;
    push @time, $dt->days . Lingua::EN::Inflect::PL_N( ' day', $dt->days ) if $dt->days;
    push @time, $dt->hours . Lingua::EN::Inflect::PL_N( ' hour', $dt->hours ) if $dt->hours;
    push @time, $dt->minutes . Lingua::EN::Inflect::PL_N( ' minute', $dt->minutes ) if $dt->minutes;
    push @time, $dt->seconds . Lingua::EN::Inflect::PL_N( ' second', $dt->seconds ) if $dt->seconds;
    my $timestr = join ', ', @time;
    
    my $fs = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("======== Tephra completed at: $fs. Elapsed time: $timestr. ========");
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

    perldoc Tephra::Role::Logger


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1;
