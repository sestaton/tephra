#!/usr/bin/perl
use strict; 
use Getopt::Long;

my $seq; 
my $phmm_file;
my $out_dir;
my $seq_file;
my $pep_file;
my $command;
my @hmm_results;
my $hmm_result;

GetOptions(    'seq=s' => \$seq,
 	       'hmmfile=s' => \$phmm_file,
	       'odir=s' => \$out_dir,
               );

$seq_file = $out_dir."aaaaa";
$pep_file = $out_dir."bbbbb";

system("echo ".$seq." > ".$seq_file);
system("transeq -frame=f ".$seq_file." -outseq=".$pep_file." 2>/dev/null");
$command = "hmmsearch ".$phmm_file." ".$pep_file;
$hmm_result = `$command`;
 
@hmm_results = split(/\n/, $hmm_result);
for (my $i=0; $i<=$#hmm_results; $i++){
    
    if ($hmm_results[$i] =~ /^\-\-\-\-\-\-\-\-\s\-\-\-\-\-\-\-\s/){
	if ($hmm_results[$i+1] =~ /^\S/){
	    my @temp = split(/\s+/, $hmm_results[$i+1]);
	    print $temp[9];
	}else{
	    print "1";
	}
	last;
    }
}
