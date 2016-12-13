package Tephra::NonLTR::Role::PathUtils;

use 5.010;
use Moose::Role;
use Config;
use File::Spec;
use IPC::System::Simple qw(capture);
use Carp 'croak';
use Tephra::Config::Exe;

=head1 NAME

Tephra::NonLTR::Role::PathUtils - Helper role for setting proper paths to programs

=head1 VERSION

Version 0.04.5

=cut

our $VERSION = '0.04.5';
$VERSION = eval $VERSION;

sub find_hmmsearch {
    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($hmmer2bin) = @{$config}{qw(hmmer2bin)};
    my $hmmsearch   = File::Spec->catfile($hmmer2bin, 'hmmsearch');

    if (-e $hmmsearch && -x $hmmsearch) {
        return $hmmsearch;
    }
    elsif (defined $ENV{HMMER2}) {
        $hmmsearch = join "/", $ENV{HMMER2}, 'src', 'hmmsearch';
        if (-e $hmmsearch && -x $hmmsearch) {
            return $hmmsearch;
        }
        else {
            $hmmsearch = join "/", $ENV{HMMER2}, 'bin', 'hmmsearch';
            if (-e $hmmsearch && -x $hmmsearch) {
                return $hmmsearch;
            }
        }
    }
    else {
	my @path = split /:|;/, $ENV{PATH};
        
        for my $p (@path) {
            $hmmsearch = File::Spec->catfile($p, 'hmmsearch');
            if (-e $hmmsearch && -x $hmmsearch) {
                my @out = capture([0..5], "hmmsearch", "-h");
                my ($version) = grep { /HMMER/ } @out;              
                if ($version =~ /HMMER (\d\.\d\w?\d+?) \(/) {
                    my $release = $1;                  
                    if ($release =~ /^2/) {
                        return $hmmsearch;
                    }
                    else {
                        croak "\nERROR: HMMER version 2 is required but was not found. Exiting.\n";
                    }
                }
            }
        }
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

    perldoc Tephra::NonLTR::Role::PathUtils


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1;
