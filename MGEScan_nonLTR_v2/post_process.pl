#!/usr/bin/perl  -w
use strict;
use Getopt::Long;

##########################################################
# get input parameter of dna file, pep file, output dir
##########################################################
my ($dna_dir, $out_dir, $rev);
get_parameter(\$dna_dir, \$out_dir, \$rev);

##########################################################
# identify full and frag
##########################################################
my $out1_dir = $out_dir."out1/";
my $out2_dir = $out_dir."out2/";
my $out_file1 = $out_dir."out2/full";   # full-length
my $out_file2 = $out_dir."out2/frag";   # fragmented
merge_thmm($out1_dir, $out2_dir, $out_file1, $out_file2, $dna_dir);

##########################################################
# convert minus coordinates to plus coordinates
##########################################################
if ($rev==1){
    my $result1 = $out_file1."_converted";
    my $result2 = $out_file2."_converted";
    convert_minus_to_plus($out_file1, $result1, $dna_dir);
    convert_minus_to_plus($out_file2, $result2, $dna_dir);
}


###########################################################
#                        SUBROUTINE                       #
###########################################################

sub convert_minus_to_plus{

    #$_[0]: minus coordinates
    #$_[1]: plus coordinates

    my %len=();
    #find the length of seq
    opendir(DIRHANDLE, $_[2]) || die ("Cannot open directory ".$_[2]);
    foreach my $name1 (sort readdir(DIRHANDLE)) {

        if ($name1 !~ /^\./ ){    
	    my $chr_file = $_[2].$name1;
	    my ($genome, $head);
	    get_sequence($chr_file, \$genome, \$head);
	    my $key = substr($name1, 0,length($name1));  #change when the file name is changed
	    $len{$key} = length($genome);
	}
    }
    closedir(DIRHANDLE);

    open OUT, ">$_[1]";
    open (IN, $_[0]);
    while(my $each_line=<IN>){
	my @temp = split(/\s+/, $each_line);
	print OUT $temp[0]."\t".eval($len{$temp[0]}-$temp[2])."\t".eval($len{$temp[0]}-$temp[1])."\t".$temp[3]."\t".$temp[4]."\n";
    }
    close(IN);
    close(OUT);
}

sub merge_thmm{

    if (-e $_[1]){
	system("rm ".$_[1]."/*");
    }else{
	system("mkdir ".$_[1]);
    }
    
    open OUT, ">$_[2]";     # full-length
    open FRAG, ">$_[3]";    # frag
    
    opendir(DIRHANDLE, $_[0]) || die ("Cannot open directory ".$_[0]);  # result from HMM
    foreach my $name (sort readdir(DIRHANDLE)) {
	
	#for each genome file
	if ($name !~ /^\./){
	    
	    my $result_file = $_[0].$name;
	    my ($end, $start, $te, $te_name,$count);
	    $end = -1000;
	    $start = -1000;
	    $te = -1;
	    $te_name ="";
	    $count=0;
	    open (IN, $result_file)|| die "Couldn't open ".$result_file;
	    while(my $each_line=<IN>){
		chomp($each_line);
		my @temp = split(/\s+/, $each_line);
		if ($te == $temp[1] + 1 && $temp[1] != 0){
		    $start = $temp[0];
		    $te = $temp[1];
		    $count += 1;
		}elsif ($temp[1] != 0 && $te >0){
		    print $each_line."\n";
		}elsif ($temp[1] == 0){
		    if ($count == 3 || ($count ==1 && $te>30)){   #full-length elements
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

			print OUT $name."\t".$temp[0]."\t".$end."\t".eval($end-$temp[0])."\t".$te_name."\n";
			
			my $seq_file = $_[1].$te_name."_full";
			open OUT1, ">>$seq_file";
			print OUT1 ">".$name."_".$temp[0]."\n";
			
			my $chr_file = $_[4].$name;    
 			my ($genome, $head);
			get_sequence($chr_file, \$genome, \$head);
			
			my $start_pos;
			my $end_pos;
			if ($temp[0]<2000){
			    $start_pos = 0;
			    #$start_pos = $temp[0];
			}else{
			    $start_pos = $temp[0]-2000;
			    #$start_pos = $temp[0];

			}
			if ($end+2000> length($genome)){
			    $end_pos = length($genome)-1;
			    #$end_pos = $end;
			}else{
			    $end_pos = $end + 2000;
			    #$end_pos = $end;

			}
			print OUT1 substr($genome, $start_pos, eval($end_pos-$start_pos+1))."\n";
			close(OUT1);
		    }elsif ($count==1){                                      #fragmented elements

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

			print FRAG $name."\t".$temp[0]."\t".$end."\t".eval($end-$temp[0])."\t".$te_name."\n";
			
			my $seq_file = $_[1].$te_name."_frag";
			open OUT1, ">>$seq_file";
			print OUT1 ">".$name."_".$temp[0]."\n";
			
			my $chr_file = $_[4].$name;    
 			my ($genome, $head);
			get_sequence($chr_file, \$genome, \$head);
			
			print OUT1 substr($genome, $temp[0], eval($end-$temp[0]+1))."\n";
			close(OUT1);

		    }

		    $te = 0;
		}else{
		    $start = $temp[0];
		    $end = $temp[0];
		    $te = $temp[1];
		    $count = 1;
		}
	    }
	    close(IN);
	}
    }
    closedir(DIRHANDLE);
    close(OUT);
    close(FRAG);
}


sub usage {
    die "Usage: post_process.pl --dna=<dna_dir_path> --out=<output_dir> --rev=[0|1]";
}


sub get_parameter{

    my ($dna, $pep, $out);

    GetOptions(
               'dna=s' => \$dna,
               'out=s' => \$out,
               'rev=s' => \$rev,
               );

    if (! -d $dna){
        print "ERROR: The directory $dna does not exist.\n";
        usage();
    }

    if (! -d $out){
        print "ERROR: The directory $out does not exist.\n";
        usage();
    }

    if ($rev ne "0" && $rev ne "1"){
        print "ERROR: The \$rev has wrong value.\n";
        usage();
    }

    ${$_[0]} = $dna;
    ${$_[1]} = $out;
    ${$_[2]} = $rev;
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

