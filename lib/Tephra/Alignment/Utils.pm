package Tephra::Alignment::Utils;

use 5.014;
use Moose;
use File::Basename;
use File::Spec;
use Cwd qw(abs_path);
use Bio::AlignIO;
use Bio::TreeIO;
#use Data::Dump::Color;
use namespace::autoclean;

=head1 NAME

Tephra::Alignment::Utils - Reusable methods for manipulating multiple sequence alignments

=head1 VERSION

Version 0.11.1

=cut

our $VERSION = '0.11.1';
$VERSION = eval $VERSION;

#
# methods
#
sub parse_aln {
    my $self = shift;
    my ($aln, $tre, $dnd) = @_;

    my ($name, $path, $suffix) = fileparse($aln, qr/\.[^.]*/);
    my $phy = File::Spec->catfile( abs_path($path), $name.'.phy' );
    
    my $aln_in  = Bio::AlignIO->new(-file  => $aln,    -format => 'clustalw');
    my $aln_out = Bio::AlignIO->new(-file  => ">$phy", -format => 'phylip', -flag_SI => 1, -idlength => 20);

    while (my $alnobj = $aln_in->next_aln) {
	$aln_out->write_aln($alnobj);
    }

    my $tre_in  = Bio::TreeIO->new(-file => $tre,    -format => 'newick');
    my $tre_out = Bio::TreeIO->new(-file => ">$dnd", -format => 'newick');

    while (my $treobj = $tre_in->next_tree) {
	for my $node ($treobj->get_nodes) {
	    my $id = $node->id;
	    next unless defined $id;
	    my $newid = substr $id, 0, 20;
	    $node->id($newid);
	}
	$tre_out->write_tree($treobj);
    }
    unlink $tre;
    
    return $phy;
}

sub check_divergence {
    my $self = shift;
    my ($phy) = @_;

    my $alnio = Bio::AlignIO->new(-file => $phy, -format => 'phylip', -longid => 1); 
    my $aln = $alnio->next_aln;
    my $pid = $aln->overall_percentage_identity;
    my $div = 100 - $pid;

    return $div;
}

sub revcom {
    my $self = shift;
    my ($seq) = @_;

    if ($seq =~ /[atcg]/i) {
        my $revcom = reverse $seq;
        $revcom =~ tr/ACGTacgt/TGCAtgca/;
        return $revcom;
    }
    else {
        say STDERR "\n[WARNING]: Not going to reverse protein sequence.\n";
        return $seq;
    }
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

    perldoc Tephra::Alignment::Utils


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
