package Postprocess;

use 5.010;
use Moose;
use MooseX::Types::Path::Class;
use Bio::SeqIO;
use File::Find;
use File::Spec;
use File::Path qw(make_path);
use File::Basename;
use Data::Printer;
#with 'SeqUtils';

has fastadir => ( is => 'ro', isa => 'Path::Class::File', required => 1, coerce => 1 );
has outdir   => ( is => 'ro', isa => 'Path::Class::Dir', required => 1, coerce => 1 );
has reverse  => ( is => 'ro', isa => 'Bool', required => 0, default => 1 );

sub postprocess {
    my $self = shift;
   ##########################################################
   # get input parameter of dna file, pep file, output dir
   ##########################################################
    my $dna_dir = $self->fastadir;
    my $out_dir = $self->outdir;
    my $rev     = $self->reverse;

    #my ($dna_dir, $out_dir, $rev);
    #get_parameter(\$dna_dir, \$out_dir, \$rev);
    
    ##########################################################
    # identify full and frag
    ##########################################################
    my $out1_dir  = File::Spec->catdir($out_dir, 'out1');
    my $out2_dir  = File::Spec->catdir($out_dir, 'out2');
    my $out_file1 = File::Spec->catdir($out_dir, 'out2', 'full');   # full-length
    my $out_file2 = File::Spec->catdir($out_dir, 'out2', 'frag');   # fragmented
    $self->merge_thmm($out1_dir, $out2_dir, $out_file1, $out_file2, $dna_dir);
    
    ##########################################################
    # convert minus coordinates to plus coordinates
    ##########################################################
    if ($rev == 1){
	my $result1 = $out_file1."_converted";
	my $result2 = $out_file2."_converted";
	$self->convert_minus_to_plus($out_file1, $result1, $dna_dir);
	$self->convert_minus_to_plus($out_file2, $result2, $dna_dir);
    }
    
 }   

sub convert_minus_to_plus {
    my $self = shift;
    #$_[0]: minus coordinates
    #$_[1]: plus coordinates
    my ($out_file1, $result1, $dna_dir) = @_;

    my %len; # = ();
    my @fasfiles;
    find( sub { push @fasfiles, $File::Find::name if -f and /\.fa.*$/ }, $dna_dir );

    #find the length of seq
    #opendir(DIRHANDLE, $_[2]) || die ("Cannot open directory ".$_[2]);
    #foreach my $name1 (sort readdir(DIRHANDLE)) {
	
    #if ($name1 !~ /^\./ ){    
    for my $file (sort @fasfiles) {
	#my $chr_file = $_[2].$name1;
	my ($genome, $head);
	$self->get_sequence($file, $genome, $head);
	#my $key = substr($, 0, length($name1));  #change when the file name is changed
	$len{$file} = length($genome);
	#}
    }
    #closedir(DIRHANDLE);
    
    #open OUT, ">$_[1]";
    #open (IN, $_[0]);
    open my $out, '>', $result1 or die "\nERROR: Could not open file: $result1";
    open my $in, '<', $out_file1 or die "\nERROR: Could not open file: $out_file1";

    while(my $each_line = <$in>){
	my @temp = split /\s+/, $each_line;
	say $out join "\t", $temp[0], eval($len{$temp[0]}-$temp[2]), eval($len{$temp[0]}-$temp[1]), $temp[3..4];
    }
    close $in;
    close $out;
}

sub merge_thmm {
    my $self = shift;
    my ($out1_dir, $out2_dir, $out_file1, $out_file2, $dna_dir) = @_;

    my @fasfiles;
    find( sub { push @fasfiles, $File::Find::name if -f and /\.fa.*$/ }, $dna_dir );

    if (-e $out2_dir) {
	my @files;
	find( sub { push @files, $File::Find::name if -f }, $out2_dir );
	
	unlink @files if @files;
	    #system("rm ".$_[1]."/*");
	    #unlink $files;
	#}
    }
    else {
	#system("mkdir ".$_[1]);
	make_path( $out2_dir, {verbose => 0, mode => 0771,} );
    }

    unless ( -e $out1_dir ) {
	make_path( $out1_dir, {verbose => 0, mode => 0771,} );
    }
    #open OUT, ">$_[2]";     # full-length
    #open FRAG, ">$_[3]";    # frag
    open my $out, '>', $out_file1 or die "\nERROR: Could not open file: $out_file1";
    open my $frag, '>', $out_file2 or die "\nERROR: Could not open file: $out_file2";

    #opendir(DIRHANDLE, $_[0]) || die ("Cannot open directory ".$_[0]);  # result from HMM
    #foreach my $name (sort readdir(DIRHANDLE)) {
    my @resfiles;
    find( sub { push @resfiles, $File::Find::name if -f }, $out1_dir );
    #p @fasfiles;

    for my $file (sort @resfiles) {
	my ($name, $path, $suffix) = fileparse($file, qr/\.[^.]*/);
	#say STDERR "name is: $name $path suffix at 110 postprocess";
	my ($chr_file) = grep { /$name/ } @fasfiles;
	next unless defined $chr_file;
	#say STDERR "chrfile: $chr_file";
	#for each genome file
	#if ($name !~ /^\./){
	    
	#my $result_file = File::Spec->catfile($out1_dir, $file); #$_[0].$name;
	#my ($end, $start, $te, $te_name,$count);
	my $end   = -1000;
	my $start = -1000;
	my $te    = -1;
	my $te_name; # ="";
	my $count = 0;
	#open (IN, $result_file)|| die "Couldn't open ".$result_file;
	open my $in, '<', $file or die "\nERROR: Could not open file: $file";
	
	while (my $each_line = <$in>){
	    chomp $each_line;
	    my @temp = split /\s+/, $each_line;
	    if ($te == $temp[1] + 1 && $temp[1] != 0){
		$start = $temp[0];
		$te    = $temp[1];
		$count += 1;
	    }
	    elsif ($temp[1] != 0 && $te > 0){
		print $each_line."\n";
	    }
	    elsif ($temp[1] == 0){
		if ($count == 3 || ($count == 1 && $te > 30)){   #full-length elements
		    if ($te <= 3 ){
			$te_name = "Jockey";
		    }elsif($te <= 6){
			$te_name = "I";
		    }elsif($te <= 9){
			$te_name = "CR1";
		    }elsif($te <= 12){
			$te_name = "Rex";
		    }elsif($te <= 15){
			$te_name = "R1";
		    }elsif($te <= 18){
			$te_name = "Tad1";
		    }elsif($te <= 21){
			$te_name = "RTE";
		    }elsif($te <= 24){
			$te_name = "L1";
		    }elsif($te <= 27){
			$te_name = "RandI";
		    }elsif($te <= 30){
			$te_name = "L2";
		    }elsif($te <= 31){
			$te_name = "CRE";
		    }elsif($te <= 32){
			$te_name = "R2";
		    }
		    
		    say $out join "\t",  $name, $temp[0], $end, eval($end-$temp[0]), $te_name;
		    
		    my $seq_file = File::Spec->catfile($out2_dir, $te_name."_full"); #$_[1].$te_name."_full";
		    open my $out1, '>>', $seq_file or die "\nERROR: Could not open file: $seq_file";
		    say $out1 ">".$name."_".$temp[0];
		    
		    #my $chr_file = File::Spec->catfile($dna_dir, $name.$suffix); #$_[4].$name;    
		    #say STDERR "chrfile2: $chr_file";
		    my ($genome, $head) = $self->get_sequence($chr_file);
		    
		    my $start_pos;
		    my $end_pos;
		    if ($temp[0] < 2000){
			$start_pos = 0;
			#$start_pos = $temp[0];
		    }
		    else{
			$start_pos = $temp[0]-2000;
			#$start_pos = $temp[0];
			
		    }
		    if ($end+2000 > length($genome)){
			$end_pos = length($genome)-1;
			#$end_pos = $end;
		    }
		    else{
			$end_pos = $end + 2000;
			#$end_pos = $end;
			
			}
		    say $out1 substr($genome, $start_pos, eval($end_pos-$start_pos+1));
		    close $out1;
		}
		elsif ($count == 1){ #fragmented elements
		    
		    if ($te <= 3 ){
			$te_name = "Jockey";
		    }elsif($te <= 6){
			$te_name = "I";
		    }elsif($te <= 9){
			$te_name = "CR1";
		    }elsif($te <= 12){
			$te_name = "Rex";
		    }elsif($te <= 15){
			$te_name = "R1";
		    }elsif($te <= 18){
			$te_name = "Tad1";
		    }elsif($te <= 21){
			$te_name = "RTE";
		    }elsif($te <= 24){
			$te_name = "L1";
		    }elsif($te <= 27){
			$te_name = "RandI";
		    }elsif($te <= 30){
			$te_name = "L2";
		    }
		    
		    say $frag join "\t", $name, $temp[0], $end, eval($end-$temp[0]), $te_name;
			
		    my $seq_file = File::Spec->catfile($out2_dir, $te_name."_frag"); #$_[1].$te_name."_frag";
		    open my $out1, '>>', $seq_file or die "\nERROR: Could not open file: $seq_file";
		    say $out1 ">".$name."_".$temp[0];
		    
		    #my $chr_file = File::Spec->catfile($dna_dir, $name.$suffix); #$_[4].$name;    
		    #say STDERR "chrfile3: $chr_file";
		    my ($genome, $head) = $self->get_sequence($chr_file);
			
		    say $out1 substr($genome, $temp[0], eval($end-$temp[0]+1));
		    close $out1;
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
}

sub get_sequence {  # file name, variable for seq, variable for head                                   \
    my $self = shift;
    my ($chr_file) = @_;

    my ($genome, $head);
    #say STDERR "chrfile4: $chr_file";
    my $seqio = Bio::SeqIO->new(-file => $chr_file, -format => 'fasta');
    while (my $seqobj = $seqio->next_seq) {
	$head   = $seqobj->id;
	$genome = $seqobj->seq;
    }

    return ($genome, $head);
	#open (GENOME, $_[0])|| die("ERROR: Couldn't open genome_file $_[0]!\n");
    #while( my $each_line=<GENOME>)  {

     #   if ($each_line =~ m/>/){
     #       ${$_[1]} = "";
     #       chomp($each_line);
     #       ${$_[2]} = $each_line;
     #   }else{
     #       chomp($each_line);
     #       ${$_[1]} .= $each_line;
     #   }
    #}
    #close(GENOME);
}

__PACKAGE__->meta->make_immutable;

1;
