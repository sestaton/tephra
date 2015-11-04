package SeqUtils;

use 5.010;
use Moose;
use File::Find;
use File::Basename;
use File::Path qw(make_path);

#fasta =>$file, outdir => $plus_out_dir,hmmdir=> $phmm_dir

sub invert_seq {
    my $self = shift;
    my ($plus_dna_dir, $minus_dna_dir) = @_;

    unless ( -d $minus_dna_dir ) {
        make_path( $minus_dna_dir, {verbose => 0, mode => 0771,} );
    }


    my @fasfiles;
    find( sub { push @fasfiles, $File::Find::name if -f }, $plus_dna_dir );

    #opendir(DIRHANDLE, $_[0]) || die ("Cannot open directory ".$_[0]);
    #foreach my $name1 (sort readdir(DIRHANDLE)) {
    for my $file (@fasfiles) {
	my ($name, $path, $suffix) = fileparse($file, qr/\.[^.]*/);
        #my $file = $_[0].$name1;
        #open (IN, $file)||die "Couldn't open ".$file;
	open my $in, '<', $file or die "\nERROR: Could not open file: $file";
        my @temp = <$in>;
        close $in;

        #if ($temp[0] =~ /\>/){
	shift @temp if $temp[0] =~ /^\>/;
        #}
        chomp @temp;
        my $seq1 = join "", @temp;
        my $seq2 = reverse $seq1;
        $seq2 =~ tr/[A,C,G,T,a,c,g,t]/[T,G,C,A,t,g,c,a]/;
        
        my $header = ">".$name;
        my $outfile = File::Spec->catfile($minus_dna_dir, $name);
        open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile";;
        say $out $header;

        my $i = 0;
        my $seq_len = length($seq2);

        while ($i<$seq_len-60){
            say $out substr($seq2, $i, 60);
            $i += 60;
        }
        say $out substr($seq2, $i, $seq_len-$i);
        close $out;
    }
}

__PACKAGE__->meta->make_immutable;

1;
