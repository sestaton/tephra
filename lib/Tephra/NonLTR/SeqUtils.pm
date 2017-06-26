package Tephra::NonLTR::SeqUtils;

use 5.014;
use Moose;
use File::Find;
use File::Basename;
use File::Path qw(make_path);
use Cwd        qw(abs_path);

=head1 NAME

Tephra::NonLTR::SeqUtils - Minor sequence utilities for non-LTR finding

=head1 VERSION

Version 0.08.2

=cut

our $VERSION = '0.08.2';
$VERSION = eval $VERSION;

sub invert_seq {
    my $self = shift;
    my ($plus_dna_dir, $minus_dna_dir) = @_;

    unless ( -d $minus_dna_dir ) {
        make_path( abs_path($minus_dna_dir), {verbose => 0, mode => 0771,} );
    }

    my @fasfiles;
    find( sub { push @fasfiles, $File::Find::name if -f }, $plus_dna_dir );

    my @revfasfiles;
    for my $file (@fasfiles) {
	my ($name, $path, $suffix) = fileparse($file, qr/\.[^.]*/);

	open my $in, '<', $file or die "\nERROR: Could not open file: $file";
        my @temp = <$in>;
        close $in;

	shift @temp if $temp[0] =~ /^\>/;
        chomp @temp;
        my $seq = join "", @temp;
        my $revseq = reverse $seq;
        $revseq =~ tr/[A,C,G,T,a,c,g,t]/[T,G,C,A,t,g,c,a]/;
        $revseq =~ s/.{60}\K/\n/g;
        my $outfile = File::Spec->catfile($minus_dna_dir, $name.$suffix);
        open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile";;
	say $out join "\n", ">".$file, $revseq;
	close $out;
	push @revfasfiles, $outfile;
    }

    return \@revfasfiles;
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

    perldoc Tephra::NonLTR::SeqUtils


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
