package SeqUtils;

use 5.010;
use Moose;
use File::Find;
use File::Basename;
use File::Path qw(make_path);

sub invert_seq {
    my $self = shift;
    my ($plus_dna_dir, $minus_dna_dir) = @_;

    unless ( -d $minus_dna_dir ) {
        make_path( $minus_dna_dir, {verbose => 0, mode => 0771,} );
    }

    my @fasfiles;
    find( sub { push @fasfiles, $File::Find::name if -f }, $plus_dna_dir );

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
    }
}

__PACKAGE__->meta->make_immutable;

1;
