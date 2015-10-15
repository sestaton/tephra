#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Data::Dump;
use Data::Printer;

#our $hmmsearch = '/home/statonse/Desktop/hannuus_annotation/non-ltr/MGEScan_nonLTR_v2/hmmer-2.3.2/src/hmmsearch';
##########################################################
# get input parameter of dna file, pep file, output dir
##########################################################
#print "Getting input parameter...\n";
my ($dna_file, $pep_file, $out_dir, $dna_name, $command);
my ($out1_dir, $out_file, $pos_dir, $phmm_dir, $pdir);
get_parameter(\$dna_file, \$out_dir, \$phmm_dir, \$pdir);
get_id(\$dna_file, \$dna_name);

#say STDERR join q{ }, $dna_file, $dna_name;
#exit;

$out1_dir = $out_dir."out1/";
if (-e $out1_dir){
}else{
    system("mkdir ".$out1_dir);
}
$pos_dir = $out_dir."pos/";
if (-e $pos_dir){
}else{
    system("mkdir ".$pos_dir);
}


##########################################################
# get signal for some state of ORF1, RT, and APE
# need HMMSEARCH
##########################################################
print "Getting signal...\n";
my ($phmm_file, $domain_rt_pos_file, $domain_ape_pos_file, $domain_orf1_pos_file);

print "    Protein sequence...\n";
$pep_file = $out_dir.$dna_name.".pep";
$command = $pdir."translate -d ".$dna_file." -h ".$dna_name." -p ".$pep_file;
say STDERR "translate cmd: ";
say STDERR $command;
system($command);
#exit;

print "    RT signal...\n";
$phmm_file = $phmm_dir."ebi_ds36752_seq.hmm3";
$domain_rt_pos_file = $pos_dir.$dna_name.".rt.pos";
get_signal_domain($pep_file, $phmm_file, $domain_rt_pos_file);
#exit;

print "    APE signal...\n";
$phmm_file = $phmm_dir."ebi_ds36736_seq.hmm3";
$domain_ape_pos_file = $pos_dir.$dna_name.".ape.pos";
get_signal_domain($pep_file, $phmm_file, $domain_ape_pos_file);

##############################################################################
# generate corresponsing empty domains files if either of them does not exist 
##############################################################################
if (-e $domain_rt_pos_file  || -e $domain_ape_pos_file ){
    
    print $dna_name."\n";

    if (! -e $domain_rt_pos_file){
	open OUT, ">$domain_rt_pos_file";
	print OUT "";
	close(OUT);
    }elsif (! -e $domain_ape_pos_file){
	open OUT, ">$domain_ape_pos_file";
	print OUT "";
	close(OUT);
    }

    $command = $pdir."match_pos.pl -rt=".$domain_rt_pos_file." -ape=".$domain_ape_pos_file;
    #system($command);
    ###########################################################
    # run hmm
    ###########################################################
    #print "Running HMM...\n";

    $out_file = $out1_dir.$dna_name;
    $command = $pdir."hmm/MGEScan -m ".$pdir."hmm/chr.hmm -s ".$dna_file." -r ".$domain_rt_pos_file." -a ".$domain_ape_pos_file." -o ".$out_file." -p ".$pdir." -d ".$out1_dir;
    #print $command."\n";
    system($command); 
}

 
if (-e $pep_file){
    #system("rm ".$pep_file);
}

###########################################################
#                        SUBROUTINE                       #
###########################################################


sub get_signal_domain{
    my ($pep_file, $phmm_file, $domain_rt_pos_file) = @_;
    #$_[0]: pep seq file
    #$_[1]: domain hmm file
    #$_[2]: output domain dna position file

    my $hmm_command = "hmmsearch  -E 0.00001 $phmm_file $pep_file";;
    my $hmm_result = `$hmm_command`;
    my %domain_start=();
    my %domain_end=();
    my $evalue;
    my $temp_file   = $domain_rt_pos_file."temp";
    my $temp_file2  = $domain_rt_pos_file."temp2";
    my $output_file = $domain_rt_pos_file;

    say STDERR "TEMP: $temp_file";
    # run hmmsearch to find the domain and save it in the temprary file    
    open (OUT, ">$temp_file");
    while ($hmm_result =~ /((\S)+\s+\d+\/\d+\s+\d+\s+\d+\s+(\[|\.)(\]|\.)\s+\d+\s+\d+\s+(\[|\.)(\]|\.)\s+(-)*\d+\.\d+\s+((\d|\-|\.|e)+))\s*/g){

	say STDERR "match: $1";
        my @temp = split(/\s+/, $1);
	p @temp;
	if ($temp[9]<0.001 ){
	    print OUT eval($temp[2]*3)."\t".eval($temp[3]*3)."\t".$temp[5]."\t".$temp[6]."\t".$temp[8]."\t".$temp[9]."\n";
	}
    }
    close(OUT);
    if (-s $temp_file >0){

        system("sort +0 -1n ".$temp_file." > ".$temp_file2);

        my $start = -1;
        my $end = -1;
        my @pre = (-1000, -1000, -1000, -1000, -1000, -1000);
	say STDERR "temp2: $temp_file2";
        open(IN, $temp_file2);
        open OUT, ">$output_file";
        while(my $each_line=<IN>){
            my @temp = split(/\s+/, $each_line);

            if ($temp[0] - $pre[1] < 300 ) {
                $end = $temp[1];
                $evalue = $evalue * $temp[5];
            }else{
                if($start>=0 && $evalue < 0.00001){
                    print OUT $start."\t".$end."\t".$pre[4]."\t".$pre[5]."\n";
                }
                $start = $temp[0];
                $end = $temp[1];
                $evalue = $temp[5];
            }
            @pre = @temp;
        }
        if($start>=0 && $evalue < 0.00001){ 
            print OUT $start."\t".$end."\t".$pre[4]."\t".$pre[5]."\n";
        }
        close(IN);
        close(OUT);
        #system("rm ".$temp_file2);
    }
    #system("rm ".$temp_file);
}

sub get_id{

    my @temp = split(/\//, ${$_[0]});
    ${$_[1]} = $temp[$#temp];
}


sub usage {
    die "Usage: run_hmm.pl --dna=<dna_file_path>  --out=<output_dir> --pdir=<program dir>";
}


sub get_parameter{

    my ($dna, $pep, $out, $phmm, $pdir);

    GetOptions(
               'dna=s' => \$dna,
               'out=s' => \$out,
	       'phmm=s' => \$phmm,
	       'pdir=s' => \$pdir,
               );

    if (! -e $dna){
        print "ERROR: The file $dna does not exist.\n";
        usage();
    }
    if (! -e $phmm){
        print "ERROR: The file $phmm does not exist.\n";
        usage();
    }
    if (! -e $pdir){
        print "ERROR: The file $phmm does not exist.\n";
        usage();
    }
    if (! -d $out){
	system("mkdir ".$out);
    }

    ${$_[0]} = $dna;
    ${$_[1]} = $out;
    ${$_[2]} = $phmm;
    ${$_[3]} = $pdir;
}



sub get_sequence{  # file name, variable for seq, variable for head                                   \
                                                                                                       

    open(GENOME, $_[0])|| die("ERROR: Couldn't open genome_file $_[0]!\n");
    while( my $each_line=<GENOME>)  {

        if ($each_line =~ m/>/){
            ${$_[1]} = "";
            chomp($each_line);
            ${$_[2]} = $each_line;
        }else{
            chomp($each_line);
            ${$_[1]} .= $each_line;
        }
    }
    close(GENOME);
}

