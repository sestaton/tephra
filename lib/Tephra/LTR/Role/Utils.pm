package Tephra::LTR::Role::Utils;

use 5.014;
use Moose::Role;
use File::Spec;
use File::Find;
use File::Basename;
use Bio::DB::HTS::Kseq;
use Carp 'croak';
use namespace::autoclean;
#use Data::Dump::Color;

=head1 NAME

Tephra::LTR::Role::Utils - Common utility methods for working with LTR retrotransposons

=head1 VERSION

Version 0.11.1

=cut

our $VERSION = '0.11.1';
$VERSION = eval $VERSION;

sub get_exemplar_ltrs {
    my $self = shift;
    my ($dir) = @_;

    my ($ltrfile, @ltrseqs, %ltrfams);
    find( sub { $ltrfile = $File::Find::name if -f and /exemplar_repeats.fasta$/ }, $dir);
    unless (defined $ltrfile) {
	say STDERR "\n[WARNING]: Exemplar LTR file not found in $dir.\n";
	return;
    }
	    
    if ($ltrfile =~ /^RL|family\d+/) {
	croak "\n[ERROR]: Expecting a single file of LTR exemplar sequences but it appears this command has ".
	    "been run before. This will cause problems. Please re-run 'classifyltrs' or report this issue. Exiting.\n";
	return;
    }

    my $kseq = Bio::DB::HTS::Kseq->new($ltrfile);
    my $iter = $kseq->iterator();

    while ( my $seq = $iter->next_seq() ) {
	my $id  = $seq->name;
	my $seq = $seq->seq;
	if ($id =~ /^[35]prime_(RL[CGX]_family\d+)_LTR_retrotransposon.*/) {
	    my $family = $1;
	    push @{$ltrfams{$family}}, { id => $id, seq => $seq };
	}
    }

    for my $family (keys %ltrfams) {
	my $outfile = File::Spec->catfile($dir, $family.'_exemplar_ltrseqs.fasta');
	open my $out, '>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
	for my $pair (@{$ltrfams{$family}}) {
	    say $out join "\n", ">".$pair->{id}, $pair->{seq};
	}
	close $out;
	push @ltrseqs, $outfile;
    }

    return \@ltrseqs;
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

    perldoc Tephra::LTR::Role::Utils


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1;
