package Tephra::LTR::Role::Config;

use 5.014;
use Moose::Role;

=head1 NAME

Tephra::LTR::Role::Config - Attributes and routines for parsing Tephra LTR configuration. 

=head1 VERSION

Version 0.05.0

=cut

our $VERSION = '0.05.0';
$VERSION = eval $VERSION;

=head1 SYNOPSIS

    use Tephra::LTR::LTRSearch;

    with 'Tephra::LTR::Role::Config'
    ...

=cut

=head1 METHODS

=head2 parse_configuration

 Title    : parse_config

 Usage    : my $config = $trans_obj->parse_configuration;

 Function : The parsed configuration for Tephra LTR options.

                                                           Return_type
 Returns  : A hash containing the user-set configuration   HashRef
            for how transposome is to be executed

 Args     : None. This is a role that can be consumed.

=cut 

sub parse_configuration {
    my $self = shift;
    my ($yaml) = @_;
    my %config;
    my $index = 0;

    # ltrharvest options from config
    $index = 0;
    $config{mintsd}     = $yaml->[0]->{ltrharvest}->[$index]->{mintsd};
    $index++;
    $config{maxtsd}     = $yaml->[0]->{ltrharvest}->[$index]->{maxtsd};
    $index++;
    $config{minlenltr}  = $yaml->[0]->{ltrharvest}->[$index]->{minlenltr};
    $index++;
    $config{maxlenltr}  = $yaml->[0]->{ltrharvest}->[$index]->{maxlenltr};
    $index++;
    $config{mindistltr} = $yaml->[0]->{ltrharvest}->[$index]->{mindistltr};
    $index++;
    $config{maxdistltr} = $yaml->[0]->{ltrharvest}->[$index]->{maxdistltr};
    $index++;
    $config{seedlength} = $yaml->[0]->{ltrharvest}->[$index]->{seedlength};
    $index++;
    $config{tsdradius}  = $yaml->[0]->{ltrharvest}->[$index]->{tsdradius};
    $index++;
    $config{xdrop}      = $yaml->[0]->{ltrharvest}->[$index]->{xdrop};
    $index++;
    $config{swmat}      = $yaml->[0]->{ltrharvest}->[$index]->{swmat};
    $index++;
    $config{swmis}      = $yaml->[0]->{ltrharvest}->[$index]->{swmis};
    $index++;
    $config{swins}      = $yaml->[0]->{ltrharvest}->[$index]->{swins};
    $index++;
    $config{swdel}      = $yaml->[0]->{ltrharvest}->[$index]->{swdel};
    $index++;
    $config{overlaps}   = $yaml->[0]->{ltrharvest}->[$index]->{overlaps};

    # ltrdigest options from config
    $index = 0;
    $config{pptradius}      = $yaml->[0]->{ltrdigest}->[$index]->{pptradius};
    $index++;
    $config{pptlen}         = $yaml->[0]->{ltrdigest}->[$index]->{pptlen};
    $index++;
    $config{pptagpr}        = $yaml->[0]->{ltrdigest}->[$index]->{pptagpr};
    $index++;
    $config{uboxlen}        = $yaml->[0]->{ltrdigest}->[$index]->{uboxlen};
    $index++;
    $config{uboxutpr}       = $yaml->[0]->{ltrdigest}->[$index]->{uboxutpr};
    $index++;
    $config{pbsradius}      = $yaml->[0]->{ltrdigest}->[$index]->{pbsradius};
    $index++;
    $config{pbslen}         = $yaml->[0]->{ltrdigest}->[$index]->{pbslen};
    $index++;
    $config{pbsoffset}      = $yaml->[0]->{ltrdigest}->[$index]->{pbsoffset};
    $index++;
    $config{pbstrnaoffset}  = $yaml->[0]->{ltrdigest}->[$index]->{pbstrnaoffset};
    $index++;
    $config{pbsmaxeditdist} = $yaml->[0]->{ltrdigest}->[$index]->{pbsmaxeditdist};
    $index++;
    $config{pdomevalue}     = $yaml->[0]->{ltrdigest}->[$index]->{pdomevalue};
    $index++;
    $config{pdomcutoff}     = $yaml->[0]->{ltrdigest}->[$index]->{pdomcutoff};
    $index++;
    $config{maxgaplen}      = $yaml->[0]->{ltrdigest}->[$index]->{maxgaplen};

    my $valid_config = $self->_validate_params(\%config);

    return $valid_config;
}

=head2 _validate_params

 Title    : _validate_params

 Usage    : This is a private method, do not use it directly.

 Function : Valiadate whether all of the slots in config file
            have been set.

                                                           Return_type
 Returns  : A hash containing the user-set configuration   HashRef
            for how transposome is to be executed

 Args     : None. This is a role that can be consumed.

=cut 

sub _validate_params {
    my $self = shift;
    my ($config) = @_;

    for my $k (keys %$config) {
	my $v = $config->{$k};
        if (not defined $v) {
            die "[ERROR]: '$k' is not defined after parsing configuration file.\n".
	        "         This indicates there may be a blank line in your configuration file.\n".
	        "         Please check your configuration file and try again. Exiting.\n";
        }
        elsif ($v =~ /^~/) {
            $v =~ s/^~/$ENV{"HOME"}/;
            $config->{$k} = $v;
        }
    }
    
    return $config;
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

    perldoc Tephra::LTR::Role::Config


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015 S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1; 

