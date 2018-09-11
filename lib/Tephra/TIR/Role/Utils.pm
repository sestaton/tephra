package Tephra::TIR::Role::Utils;

use 5.014;
use Moose::Role;
use File::Spec;
use File::Find;
use File::Basename;
use File::Path qw(remove_tree);
use Bio::DB::HTS::Kseq;
use Carp 'croak';
use namespace::autoclean;
#use Data::Dump::Color;

=head1 NAME

Tephra::TIR::Role::Utils - Common utility methods for working with TIR transposons

=head1 VERSION

Version 0.12.1

=cut

our $VERSION = '0.12.1';
$VERSION = eval $VERSION;

sub get_exemplar_tirs_for_age {
    my $self = shift;
    my ($dir, $outdir) = @_;

    my (@dirs, @tirseqs);
    find( sub { push @dirs, $File::Find::name if -d && /_hAT\z|_mite\z|_mutator\z|_tc1-mutator\z|_unclassified\z/ }, $dir);
    croak "\n[ERROR]: Could not find the expected sub-directories ending in 'hAT','mite','mutator','tc1-mutator' and 'unclassified'. Please ".
        "check input. Exiting.\n" unless @dirs;

    for my $sfdir (@dirs) {
	my ($tirfile, %tirfams);
	find( sub { $tirfile = $File::Find::name if -f and /exemplar_repeats.fasta$/ }, $sfdir);
	unless (defined $tirfile) {
	    say STDERR "\n[WARNING]: No exemplar TIR file was found in: $sfdir.";
	    say STDERR "This is likely because there were no families identified by the 'classifytirs' command for this superfamily.";
	    say STDERR "You can try the 'age' command again with the --all flag to process all TIRs.\n";
	    remove_tree( $outdir, { safe => 1 } )
		if $self->clean;
	    exit;
	}

	my $kseq = Bio::DB::HTS::Kseq->new($tirfile);
	my $iter = $kseq->iterator();
	
	my $re = qr/terminal_inverted_repeat_element\d+|MITE\d+/;
	while ( my $seq = $iter->next_seq() ) {
	    my $id  = $seq->name;
	    my $seq = $seq->seq;
	    #if ($id =~ /^[35]prime_(RL[CGX]_family\d+)_(?:LTR|TRIM)_retrotransposon.*/) {
	    if ($id =~ /^(?:[35]prime_)?(\w{3}_)?((?:singleton_)?(?:family\d+_))?$re?_\S+_\d+[-_]\d+/) {
		my $code = $1;
		my $family = $2;
		$family =~ s/_//g;
		push @{$tirfams{$code.$family}}, { id => $id, seq => $seq };
	    }
	}

	for my $family (keys %tirfams) {
	    my $outfile = File::Spec->catfile($outdir, $family.'_exemplar_tirseqs.fasta');
	    open my $out, '>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
	    for my $pair (@{$tirfams{$family}}) {
		say $out join "\n", ">".$pair->{id}, $pair->{seq};
	    }
	    close $out;
	    push @tirseqs, $outfile;
	}
    }

    return \@tirseqs;
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

    perldoc Tephra::TIR::Role::Utils


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1;
