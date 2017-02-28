package Tephra::NonLTR::Postprocess;

use 5.014;
use Moose;
use MooseX::Types::Path::Class;
use Bio::DB::HTS::Kseq;
use File::Find;
use File::Spec;
use File::Basename;
use File::Path qw(make_path);
use Cwd        qw(abs_path);
use namespace::autoclean;

=head1 NAME

Tephra::NonLTR::Postprocess - Postprocess initial scan for non-LTR coding domains

=head1 VERSION

Version 0.07.0

=cut

our $VERSION = '0.07.0';
$VERSION = eval $VERSION;

has fastadir    => ( is => 'ro', isa => 'Path::Class::File', required => 1, coerce  => 1 );
has outdir      => ( is => 'ro', isa => 'Path::Class::Dir',  required => 1, coerce  => 1 );
has reverse     => ( is => 'ro', isa => 'Bool', required => 0, default => 1 );
has n_threshold => ( is => 'ro', isa => 'Num',  required => 0, default => 0.30 );

sub postprocess {
    my $self = shift;
    my $dna_dir = $self->fastadir->absolute->resolve;
    my $out_dir = $self->outdir->absolute->resolve;
    my $rev     = $self->reverse;

    # identify full and frag
    my $outf_dir = File::Spec->catdir($out_dir, 'out1');
    my $outr_dir = File::Spec->catdir($out_dir, 'out2');
    my $outr_full_file = File::Spec->catfile($out_dir, 'out2', 'full');   # full-length
    my $outr_frag_file = File::Spec->catfile($out_dir, 'out2', 'frag');   # fragmented
    $self->merge_thmm($outf_dir, $outr_dir, $outr_full_file, $outr_frag_file, $dna_dir);
    
    # convert minus coordinates to plus coordinates
    if ($rev == 1) {
	my $full_result = $outr_full_file.'_converted';
	my $frag_result = $outr_frag_file.'_converted';
	$self->convert_minus_to_plus($outr_full_file, $full_result, $dna_dir);
	$self->convert_minus_to_plus($outr_frag_file, $frag_result, $dna_dir);
    }    
 }   

sub convert_minus_to_plus {
    my $self = shift;
    my ($out_file, $result_file, $dna_dir) = @_;

    my %len;
    my @fasfiles;
    find( sub { push @fasfiles, $File::Find::name if -f and /\.fa.*$/ }, $dna_dir );

    for my $file (sort @fasfiles) {
	my ($genome, $head) = $self->get_sequence_id($file);
	my $filename = basename($file);
	$len{$filename} = length($genome);
    }

    open my $out, '>', $result_file or die "\nERROR: Could not open file: $result_file";
    open my $in, '<', $out_file or die "\nERROR: Could not open file: $out_file";

    ## Output: seqname start end length clade
    #Ha1.fasta 172342043 172344615 2572 L1
    #Ha1.fasta 162671813 162674406 2593 L1
    #Ha1.fasta 123785696 123788295 2599 L1

    while (my $line = <$in>) {
	chomp $line;
	my @temp = split /\s+/, $line;
	say $out join "\t", 
	    $temp[0], eval($len{$temp[0]}-$temp[2]), eval($len{$temp[0]}-$temp[1]), @temp[3..4];
    }
    close $in;
    close $out;
}

sub merge_thmm {
    my $self = shift;
    my ($outf_dir, $outr_dir, $outr_full_file, $outr_frag_file, $dna_dir) = @_;
    my $n_thresh = $self->n_threshold;

    my @fasfiles;
    find( sub { push @fasfiles, $File::Find::name if -f and /\.fa.*$/ }, $dna_dir );

    if (-e $outr_dir) {
	my @files;
	find( sub { push @files, $File::Find::name if -f }, $outr_dir );
	
	unlink @files if @files;
    }
    else {
	make_path( $outr_dir, {verbose => 0, mode => 0771,} );
    }

    unless ( -e $outf_dir ) {
	make_path( $outf_dir, {verbose => 0, mode => 0771,} );
    }

    open my $out, '>', $outr_full_file or die "\nERROR: Could not open file: $outr_full_file";
    open my $frag, '>', $outr_frag_file or die "\nERROR: Could not open file: $outr_frag_file";

    my @resfiles;
    find( sub { push @resfiles, $File::Find::name if -f }, $outf_dir );

    for my $file (sort @resfiles) {
	my $filename = basename($file);
	my ($name, $path, $suffix) = fileparse($file, qr/\.[^.]*/);
	my ($chr_file) = grep { /$filename/ } @fasfiles;
	next unless defined $chr_file;
	my $end   = -1000;
	my $start = -1000;
	my $te    = -1;
	my $te_name;
	my $count = 0;
	open my $in, '<', $file or die "\nERROR: Could not open file: $file";
	
	## Input: 
	#303624192 24 26.225023
	#303623033 0 16.655420
	#194759853 24 16.298745

	while (my $each_line = <$in>){
	    chomp $each_line;
	    my @temp = split /\s+/, $each_line;
	    if ($te == $temp[1] + 1 && $temp[1] != 0) {
		$start = $temp[0];
		$te    = $temp[1];
		$count += 1;
	    }
	    elsif ($temp[1] != 0 && $te > 0){
		print $each_line."\n";
	    }
	    elsif ($temp[1] == 0) {
		if ($count == 3 || ($count == 1 && $te > 30)) { # full-length elements
		    if ($te <= 3 ) {
			$te_name = 'Jockey';
		    }
		    elsif ($te <= 6) {
			$te_name = 'I';
		    }
		    elsif ($te <= 9) {
			$te_name = 'CR1';
		    }
		    elsif ($te <= 12) {
			$te_name = 'Rex';
		    }
		    elsif ($te <= 15) {
			$te_name = 'R1';
		    }
		    elsif ($te <= 18) {
			$te_name = 'Tad1';
		    }
		    elsif ($te <= 21) {
			$te_name = 'RTE';
		    }
		    elsif ($te <= 24) {
			$te_name = 'L1';
		    }
		    elsif ($te <= 27) {
			$te_name = 'RandI';
		    }
		    elsif ($te <= 30) {
			$te_name = 'L2';
		    }
		    elsif ($te <= 31) {
			$te_name = 'CRE';
		    }
		    elsif ($te <= 32) {
			$te_name = 'R2';
		    }
		    
		    my $seq_file = File::Spec->catfile($outr_dir, $te_name.'_full'); #$_[1].$te_name."_full";
		    open my $out1, '>>', $seq_file or die "\nERROR: Could not open file: $seq_file";
		    my $header = '>'.$filename.'_'.$temp[0].'_'.$end;

		    my ($genome, $head) = $self->get_sequence_id($chr_file);

		    my $start_pos;
		    my $end_pos;
		    if ($temp[0] < 2000) {
			$start_pos = 0;
		    }
		    else{
			$start_pos = $temp[0]-2000;
		    }
		    if ($end+2000 > length($genome)){
			$end_pos = length($genome)-1;
		    }
		    else{
			$end_pos = $end + 2000;
		    }
		    
		    my $sequence = substr($genome, $start_pos, eval($end_pos-$start_pos+1));
		    $sequence =~ s/.{60}\K/\n/g;
		    say $out1 join "\n", $header, $sequence;
		    say $out join "\t",  $filename, $temp[0], $end, eval($end-$temp[0]), $te_name;

		    close $out1;
		    unlink $seq_file unless -s $seq_file;
		}
		elsif ($count == 1) { # fragmented elements		    
		    if ($te <= 3) {
			$te_name = 'Jockey';
		    }
		    elsif ($te <= 6) {
			$te_name = 'I';
		    }
		    elsif ($te <= 9) {
			$te_name ='CR1';
		    }
		    elsif ($te <= 12) {
			$te_name = 'Rex';
		    }
		    elsif ($te <= 15) {
			$te_name = 'R1';
		    }
		    elsif ($te <= 18) {
			$te_name = 'Tad1';
		    }
		    elsif ($te <= 21) {
			$te_name = 'RTE';
		    }
		    elsif ($te <= 24) {
			$te_name = 'L1';
		    }
		    elsif ($te <= 27) {
			$te_name = 'RandI';
		    }
		    elsif ($te <= 30) {
			$te_name = 'L2';
		    }
		    
		    my $seq_file = File::Spec->catfile($outr_dir, $te_name.'_frag'); #$_[1].$te_name."_frag";
		    open my $out1, '>>', $seq_file or die "\nERROR: Could not open file: $seq_file";
		    my $header = '>'.$filename.'_'.$temp[0].'_'.$end; ;

		    my ($genome, $head) = $self->get_sequence_id($chr_file);

		    my $sequence = substr($genome, $temp[0], eval($end-$temp[0]+1));
		    $sequence =~ s/.{60}\K/\n/g;
		    say $out1 join "\n", $header, $sequence;
		    say $frag join "\t",  $filename, $temp[0], $end, eval($end-$temp[0]), $te_name;

		    close $out1;
		    unlink $seq_file unless -s $seq_file;
		}
		$te = 0;
	    }
	    else {
		$start = $temp[0];
		$end   = $temp[0];
		$te    = $temp[1];
		$count = 1;
	    }
	}
	close $in;
    }
    close $out;
    close $frag;

    unlink $outr_full_file unless -s $outr_full_file;
    unlink $outr_frag_file unless -s $outr_frag_file;
}

sub get_sequence_id {
    my $self = shift;
    my ($chr_file) = @_;

    my ($genome, $head);
    my $kseq = Bio::DB::HTS::Kseq->new($chr_file);
    my $iter = $kseq->iterator;
    while (my $seqobj = $iter->next_seq) {
	$head   = $seqobj->name;
	$genome = $seqobj->seq;
    }

    return ($genome, $head);
}

sub _filterNpercent {
    my $self = shift;
    my ($sequence) = @_;

    my $length  = length($sequence);
    my $n_count = ($sequence =~ tr/Nn//);
    my $n_perc  = sprintf("%.2f",$n_count/$length);

    return $n_perc;
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

    perldoc Tephra::NonLTR::Postprocess


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
