package Tephra::Role::GFF;

use 5.010;
use Moose::Role;
use Bio::GFF3::LowLevel qw(gff3_parse_feature);
use Path::Class::File
#use Data::Dump::Color;
use namespace::autoclean;

=head1 NAME

Tephra::Role::GFF - Utility methods for working with GFF files

=head1 VERSION

Version 0.02.7

=cut

our $VERSION = '0.02.7';
$VERSION = eval $VERSION;

#
# methods
#
sub collect_gff_features {
    my $self = shift;
    my ($gff) = @_;

    my $header;
    open my $in, '<', $gff or die "\nERROR: Could not open file: $gff\n";
    while (<$in>) {
	chomp;
	next if /^###$/;
	if (/^##?\w+/) {
	    $header .= $_."\n";
	}
	else {
	    last;
	}
    }
    close $in;
    chomp $header;

    open my $gffio, '<', $gff or die "\nERROR: Could not open file: $gff\n";

    my ($start, $end, $region, $key, %features);
    while (my $line = <$gffio>) {
        chomp $line;
        next if $line =~ /^#/;
        my $feature = gff3_parse_feature( $line );
        if ($feature->{type} eq 'repeat_region') {
            $region = @{$feature->{attributes}{ID}}[0];
            ($start, $end) = @{$feature}{qw(start end)};
	    $key = join "||", $region, $start, $end;

        }
	if ($feature->{type} ne 'repeat_region') {
            if ($feature->{start} >= $start && $feature->{end} <= $end) {
		push @{$features{$key}}, $feature;
            }
        }
    }
    close $gffio;

    return ($header, \%features);
}

sub get_parent_coords {
    my $self = shift;
    my ($parent, $coord_map) = @_;

    my ($seq_id, $start, $end) = split /\|\|/, $coord_map->{$parent};
    my $pkey = join "||", $parent, $start, $end;

    return ($seq_id, $pkey);
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

    perldoc Tephra::Role::GFF


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1;
