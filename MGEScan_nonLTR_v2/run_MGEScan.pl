#!/usr/bin/env perl

use strict;
use Getopt::Long;

############################################
# INPUT
############################################
my $conf_file;
my $main_data_dir;
my $plus_dna_dir;
my $plus_out_dir;
my $minus_dna_dir;
my $minus_out_dir;
my $phmm_dir;
my $program_dir;
my $genome_dir;

print "\n\n";
GetOptions(
           'data=s' => \$main_data_dir,
           'genome=s' => \$genome_dir,
           'program=s' => \$program_dir,
           );

if (length($genome_dir)==0){
    print "ERROR: An input genome directory was not specified.\n";
    print_usage();
    exit;
}elsif (! -e $genome_dir){
    print "ERROR: The input genome directory [$genome_dir] does not exist.\n";
    print_usage();
    exit;
}else{
    if (substr($genome_dir, length($genome_dir)-1, 1) ne "/"){
	$minus_dna_dir = $genome_dir."_b/";
	$plus_dna_dir = $genome_dir."/";
	
    }else {
	$minus_dna_dir = substr($genome_dir, 0, length($genome_dir)-1)."_b/";
	$plus_dna_dir = $genome_dir;
    }
}

if (length($main_data_dir) == 0 ){
    print "ERROR: An output directory was not specified.\n";
    print_usage();
    exit;
}else{
    if (!-e $main_data_dir){
	system("mkdir ".$main_data_dir);
    }
    if (substr($main_data_dir, length($main_data_dir)-1, 1) ne "/"){
	$main_data_dir .= "/";
    }
}

if (length($program_dir)==0){
    print "ERROR: The program directory was not specified.\n";
    print_usage();
    exit;
}elsif (! -e $program_dir){
    print "ERROR: The program directory [$program_dir] does not exist.\n";
    print_usage();
    exit;
}else{
    if (substr($program_dir, length($program_dir)-1, 1) ne "/"){
	$program_dir .= "/";
    }
}

$conf_file = $program_dir."path_conf";
$phmm_dir = $program_dir."pHMM/";
$plus_out_dir=$main_data_dir."f/";
$minus_out_dir=$main_data_dir."b/";

#say STDERR "DEBUG: $conf_file $phmm_dir $plus_out_dir $minus_out_dir";
#exit;

############################################
# Forward strand
############################################ 
printf "Running forward...\n";

opendir(DIRHANDLE, $plus_dna_dir) || die ("Cannot open directory ".$plus_dna_dir);
foreach my $name (sort readdir(DIRHANDLE)) {

    if ($name !~ /^\./){  

	my $plus_dna_file = $plus_dna_dir.$name;
	my $command = $program_dir."run_hmm.pl --dna=".$plus_dna_file."  --out=".$plus_out_dir." --phmm=".$phmm_dir." --pdir=".$program_dir;
	system($command);
    }
}
#system("rm ".$plus_out_dir."out1/aaaaa"); 
#system("rm ".$plus_out_dir."out1/bbbbb");
#system("rm ".$plus_out_dir."out1/ppppp");
#system("rm ".$plus_out_dir."out1/qqqqq");

my $command = $program_dir."post_process.pl --dna=".$plus_dna_dir." --out=".$plus_out_dir." --rev=0";
#print $command."\n";
system($command);


############################################
#Backward strand
############################################
printf "Running backward...\n";
invert_seq($plus_dna_dir, $minus_dna_dir);
opendir(DIRHANDLE, $minus_dna_dir) || die ("Cannot open directory ".$minus_dna_dir);
foreach my $name (sort readdir(DIRHANDLE)) {

    if ($name !~ /^\./){  
	my $minus_dna_file = $minus_dna_dir.$name;
	my $command = $program_dir."run_hmm.pl --dna=".$minus_dna_file." --out=".$minus_out_dir." --phmm=".$phmm_dir." --pdir=".$program_dir;
	system($command);
    }
}

#system("rm ".$minus_out_dir."out1/aaaaa");
#system("rm ".$minus_out_dir."out1/bbbbb");
#system("rm ".$minus_out_dir."out1/ppppp");
#system("rm ".$minus_out_dir."out1/qqqqq");

my $command = $program_dir."post_process.pl --dna=".$minus_dna_dir." --out=".$minus_out_dir." --rev=1";
#print $command."\n";
system($command);

###########################################
#validation for Q value
###########################################

my $command = $program_dir."post_process2.pl --data_dir=".$main_data_dir." --phmm_dir=".$phmm_dir;
system($command);



sub invert_seq{

    if (!-e $_[1]){
	system("mkdir ".$_[1]);
    }    

    opendir(DIRHANDLE, $_[0]) || die ("Cannot open directory ".$_[0]);
    foreach my $name1 (sort readdir(DIRHANDLE)) {

	my $file = $_[0].$name1;
	open (IN, $file)||die "Couldn't open ".$file;
	my @temp=<IN>;
	close(IN);

	if ($temp[0] =~ /\>/){
	    shift(@temp);
	}
	chomp(@temp);
	my $seq1 = join("", @temp);
	my $seq2 = reverse($seq1);
	$seq2 =~ tr/[A,C,G,T,a,c,g,t]/[T,G,C,A,t,g,c,a]/;
	
	my $head = ">".$name1;
	my $file2 = $_[1].$name1;
	open OUT, ">$file2";
	print OUT $head."\n";

	my $i=0;
	my $seq_len = length($seq2);

	while($i<$seq_len-60){
	    print OUT substr($seq2, $i, 60)."\n";
	    $i += 60;
	}
	print OUT substr($seq2, $i, $seq_len-$i)."\n";
	close(OUT);
    }
}
sub print_usage{

    print "USAGE: ./run_MGEScan.pl -genome=[a directory name for genome sequences] -data=[a directory name for output files]  -program=[a directory name for MGEScan-nonLTR]\n\n\n\n";

}
