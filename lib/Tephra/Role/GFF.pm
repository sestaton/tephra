package Tephra::Role::GFF;

use 5.010;
use Moose::Role;
use Bio::Tools::GFF;
use Path::Class::File;
use namespace::autoclean;

=head1 NAME

Tephra::Role::GFF - Utility methods for working with GFF files

=head1 VERSION

Version 0.02

=cut

our $VERSION = '0.02';
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

    my $gffio = Bio::Tools::GFF->new( -file => $gff, -gff_version => 3 );

    my ($start, $end, $region, %features);
    while (my $feature = $gffio->next_feature()) {
	if ($feature->primary_tag eq 'repeat_region') {
	    my @string = split /\t/, $feature->gff_string;
	    ($region) = ($string[8] =~ /ID=?\s+?(repeat_region\d+)/);
	    ($start, $end) = ($feature->start, $feature->end);
	}
	next $feature unless defined $start && defined $end;
	if ($feature->primary_tag ne 'repeat_region') {
	    if ($feature->start >= $start && $feature->end <= $end) {
		push @{$features{$region.".".$start.".".$end}}, join "||", split /\t/, $feature->gff_string;
	    }
	}
    }

    return ($header, \%features);
}

sub get_source {
    my $self = shift;
    my ($ref) = @_;
    for my $feat (@$ref) {
	for my $rfeat (@$feat) {
	    my @feats = split /\|\|/, $rfeat;
	    return ($feats[0], $feats[1]);
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

    perldoc Tephra::Role::GFF


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1;
