package Tephra::Role::GFF;

use 5.010;
use Moose::Role;
use Bio::Tools::GFF;
use Path::Class::File;
use namespace::autoclean;

#
# methods
#
sub collect_gff_features {
    my $self = shift;
    my ($gff) = @_;
    #my $gff = $self->gff;

    my $header;
    open my $in, '<', $gff;
    #my $in = $gff->open('r') or die "\n[ERROR]: Could not open file: $gff\n";
    while (<$in>) {
	chomp;
	if (/^#/) {
	    $header .= $_."\n";
	}
	else {
	    last;
	}
    }
    close $gff;
    #$in->close;
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

1;
