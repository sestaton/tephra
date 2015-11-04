#!/usr/bin/perl  -w
use strict;
use Getopt::Long;

our $hmmsearch = '/home/statonse/github/tephra/MGEScan_nonLTR_v2/hmmer-2.3.2/src/hmmsearch';
my @all_clade = ('CR1', 'I', 'Jockey', 'L1', 'L2', 'R1', 'RandI', 'Rex', 'RTE', 'Tad1', 'R2','CRE');
my @en_clade = ('CR1', 'I', 'Jockey', 'L1', 'L2', 'R1', 'RandI', 'Rex', 'RTE', 'Tad1');
my $genome="";   
my $dir;
my $hmm_dir;
my $domain;
my $seq_dir;
my $validation_dir;
my $validation_file;
my $evalue_file;
my $tree_dir;
my $seq;

get_parameter(\$dir, \$hmm_dir);

# copy seq from out2 dir into info dir
system("mkdir ".$dir."info/");
system("mkdir ".$dir."info/full");
get_full_frag($genome, $dir, \@all_clade);

# get domain seq
get_domain_for_full_frag($genome, "en", \@en_clade, $dir, $hmm_dir);
get_domain_for_full_frag($genome, "rt", \@all_clade, $dir, $hmm_dir);

# get Q value after running pHMM for EN in full elements
$validation_dir = $dir."info/validation/";
system("mkdir ".$validation_dir);

$domain = "en";
$seq_dir = $dir."info/full/";
$validation_file = $validation_dir.$domain;
$evalue_file = $validation_dir.$domain."_evalue";

opendir(DIRHANDLE, $seq_dir) || die ("Cannot open directory ".$seq_dir);
foreach my $name (sort readdir(DIRHANDLE)) {
    if ($name !~ /^\./ && $name ne "R2" && $name ne "CRE" ){
	$seq = $seq_dir.$name."/".$name.".".$domain.".pep";
	vote_hmmsearch($seq, $hmm_dir, $domain, $validation_file, $evalue_file, \@en_clade);
    }
}
closedir(DIRHANDLE);

# get Q value after running pHMM for RT in full elements
$domain = "rt";
$seq_dir = $dir."info/full/";
$validation_file = $validation_dir.$domain;
$evalue_file = $validation_dir.$domain."_evalue";

opendir(DIRHANDLE, $seq_dir) || die ("Cannot open directory ".$seq_dir);
foreach my $name (sort readdir(DIRHANDLE)) {
    if ($name !~ /^\./){
	$seq = $seq_dir.$name."/".$name.".".$domain.".pep";
	vote_hmmsearch($seq, $hmm_dir, $domain, $validation_file, $evalue_file, \@all_clade);
    }
}
close(DIRHANDLE);
#system("rm ".$seq_dir."*/*.pep");
#system("rm ".$seq_dir."*/*.rt.*");
#system("rm ".$seq_dir."*/*.en.*");
#system("rm -r ".$dir."b");
#system("rm -r ".$dir."f");


sub get_parameter{

    my ($dir, $hmm_dir);

    GetOptions(
               'data_dir=s' => \$dir,
               'phmm_dir=s' => \$hmm_dir,
               );

    if (! -e $dir){
        print "ERROR: The directory $dir does not exist.\n";
        usage();
    }

    if (! -e $hmm_dir){
        print "ERROR: The directory $hmm_dir does not exist.\n";
        usage();
    }

    ${$_[0]} = $dir;
    ${$_[1]} = $hmm_dir;
}



sub get_domain_for_full_frag{

    my $genome = $_[0];
    my $domain = $_[1];
    my @all_clade = @{$_[2]};
    my $hmm_dir = $_[4];
    my ($dir);
    my ($pep_file, $dna_file);
    my ($phmm_file, $result_pep_file, $result_dna_file);

    for (my $i=0; $i<=$#all_clade; $i++){

	my $clade = $all_clade[$i];

	$dir = $_[3]."info/full/".$clade."/";
	$pep_file = $dir.$clade.".pep";
	$dna_file = $dir.$clade.".dna";
	
	if (-e $pep_file ){

	    $phmm_file = $hmm_dir.$clade.".".$domain.".hmm";
	    $result_pep_file = $dir.$clade.".".$domain.".pe";
	    $result_dna_file = $dir.$clade.".".$domain.".dna";
	    
	    my $flag = 2;  #1: protein-protein, 2: protein-dna 
	    get_domain_pep_seq($pep_file, $phmm_file, $result_pep_file);
	    get_domain_dna_seq($pep_file, $phmm_file, $result_dna_file, $dna_file, $flag);

	    my $command = "sed 's/>/>".$clade."_/' ".$result_pep_file." > ".$result_pep_file."p";
	    system($command);
	    #system("rm ".$result_pep_file);
	    
	}
    }
}




sub get_full_frag{

    my $genome=$_[0];
    my $dir = $_[1];
    my @all_clade = @{$_[2]};
    my $file;
    my $clade_dir;
    my ($dna_file, $pep_file, $file_f, $file_b);


    for (my $i=0; $i<=$#all_clade; $i++){

	# create a clade dir
	$clade_dir = $dir."info/full/".$all_clade[$i]."/";
	$file_f = $dir."f/out2/".$all_clade[$i]."_full";
	$file_b = $dir."b/out2/".$all_clade[$i]."_full";
	if (-e $file_f || -e $file_b){
	    system("mkdir ".$clade_dir);
	}

	# copy full length in + strand
	if (-e $file_f){
	    my $command = "cat ".$file_f." > ".$clade_dir.$all_clade[$i].".dna";
	    system($command);
	}

	# copy full length in - strand
	if (-e $file_b){
	    my $command = "cat ".$file_b." >> ".$clade_dir.$all_clade[$i].".dna";
	    system($command);
	}

	# translate
	$dna_file = $dir."info/full/".$all_clade[$i]."/".$all_clade[$i].".dna";
	$pep_file = $dir."info/full/".$all_clade[$i]."/".$all_clade[$i].".pep";	   
	if (-e $dna_file){
	    my $command = "transeq -frame=f -sequence=".$dna_file." -outseq=".$pep_file." 2>/dev/null";
	    system($command);
	}
    }
}


sub get_domain_pep_seq{
    
    #$_[0]: pep seq file
    #$_[1]: domain hmm file
    #$_[2]: output domain pep seq file 
    
    my $hmm_command = "$hmmsearch  ".$_[1]." ".$_[0];
    my $hmm_result = `$hmm_command`;
    my %domain_start=();
    my %domain_end=();
    my %result_start=();
    my %result_end=();
    my %uniq_head=();


    while ($hmm_result =~ /((\d|\w|\-|\_|\#|\/|\.)+\s+\d+\/\d+\s+\d+\s+\d+\s+(\[|\.)(\]|\.)\s+\d+\s+\d+\s+(\[|\.)(\]|\.)\s+(\-)*\d+\.\d+\s+((\d|\-|\.|e)+))\s*/g){
	    
	my @temp = split(/\s+/, $1);
#	if ($temp[9]<0.000001 ){
	    my $key = substr($temp[0],0,length($temp[0]));
	    my $uniq_key = substr($temp[0],0,length($temp[0])-2);

	    if (exists $uniq_head{$uniq_key}){
	    }else{
		$uniq_head{$uniq_key} = 1;
		$result_start{$key} = $temp[2];
		$result_end{$key} = $temp[3];
	    }
#	}
    }

    my $flag=0;
    my $head="";
    my $seq="";
    open (IN, $_[0]);
    open OUT, ">$_[2]";
    while(my $each_line=<IN>){
	chomp($each_line);
	if ($each_line =~ /\>/){
	    if (length($head)>0 && $flag==1 ){
		print OUT ">".$head."\n";
		print OUT substr($seq, $result_start{$head}, eval($result_end{$head}-$result_start{$head}+1))."\n";
	    }
	    my @temp = split(/\s+/, $each_line);
	    if (exists $result_start{substr($temp[0], 1, length($temp[0])-1)}){
		$flag=1;
		$head = substr($temp[0], 1, length($temp[0])-1);
	    }else{
		$flag=0;
	    }
	    $seq="";
	}else{
	    if($flag==1){
	        $seq .= $each_line;
	    }
	}
    }
    if($flag==1){
	print OUT ">".$head."\n";
	print OUT substr($seq, $result_start{$head}, eval($result_end{$head}-$result_start{$head}+1))."\n";
    }
    close(IN);
    close(OUT);
}


sub get_domain_dna_seq{
    
    #$_[0]: pep seq file
    #$_[1]: domain hmm file
    #$_[2]: output domain dna seq file 
    #$_[3]: dna seq file

    my $hmm_command = "$hmmsearch  ".$_[1]." ".$_[0];
    my $hmm_result = `$hmm_command`;
    my %domain_start=();
    my %domain_end=();
    my %result_start=();
    my %result_end=();
    my %uniq_head=();

    while ($hmm_result =~ /((\d|\w|\-|\_|\#|\/|\.)+\s+\d+\/\d+\s+\d+\s+\d+\s+(\[|\.)(\]|\.)\s+\d+\s+\d+\s+(\[|\.)(\]|\.)\s+(\-)*\d+\.\d+\s+((\d|\-|\.|e)+))\s*/g){
	    
	my @temp = split(/\s+/, $1);

	    my $key = substr($temp[0],0,length($temp[0]));
	    my $uniq_key = substr($temp[0],0,length($temp[0])-2);

	    if (exists $result_start{$uniq_key}){
	    }else{
		$uniq_head{$uniq_key} = 1;
		if ($_[4] ==1){
		    $result_start{$key} = $temp[2];
		    $result_end{$key} = $temp[3];
		}elsif($_[4] ==2){
		    $result_start{$uniq_key} = $temp[2];
		    $result_end{$uniq_key} = $temp[3];
		}
	    }
    }
    my $flag=0;
    my $head="";
    my $seq="";
    open (IN, $_[3]);
    open OUT, ">$_[2]";
    while(my $each_line=<IN>){
	chomp($each_line);
	if ($each_line =~ /\>/){
	    if (length($head)>0 && $flag==1 ){
		print OUT ">".$head."\n";
		print OUT substr($seq, $result_start{$head}*3-3, eval(($result_end{$head}-$result_start{$head}+1)*3+3))."\n";
	    }
	    my @temp = split(/\s+/, $each_line);
	    if (exists $result_start{substr($temp[0], 1, length($temp[0])-1)}){
		$flag=1;
		$head = substr($temp[0], 1, length($temp[0])-1);
	    }else{
		$flag=0;
	    }
	    $seq="";
	}else{
	    if($flag==1){
	        $seq .= $each_line;
	    }
	}
    }
    if($flag==1){
	print OUT ">".$head."\n";
	print OUT substr($seq, $result_start{$head}*3-3, eval(($result_end{$head}-$result_start{$head}+1)*3+3))."\n";
    }
    close(IN);
    close(OUT);
}



sub vote_hmmsearch{

    my @line = @{$_[5]};
    my %evalue;
    my %save_evalue;
    my %clade;
    my %sig;
    my $i;
    my $anno_clade; 

    open (IN, $_[0]);
    while(my $each_line=<IN>){
	if ($each_line=~ /\>/){
	    chomp($each_line);
	    my $uniq_key = substr($each_line,1,length($each_line)-1);
	    $evalue{$uniq_key} = 1000;
	    $save_evalue{$uniq_key} = 1000;
	    $clade{$uniq_key} = "-";
	    $sig{$uniq_key} = 1;
	}
    }
    close(IN);

    for ($i=0; $i<=$#line; $i++){

        my $command = "$hmmsearch ".$_[1].$line[$i].".".$_[2].".hmm ".$_[0];
        my $hmm_result = `$command`;

        while ($hmm_result =~ /((\d|\w|\-|\_|\#|\/|\.)+\s+\d+\/\d+\s+\d+\s+\d+\s+(\[|\.)(\]|\.)\s+\d+\s+\d+\s+(\[|\.)(\]|\.)\s+\-*\d+\.\d+\s+((\d|\-|\.|e)+))\s*/g){
            my @temp = split(/\s+/, $1);
            my $uniq_key = substr($temp[0],0,length($temp[0]));

	    $save_evalue{$uniq_key} = $save_evalue{$uniq_key}."\t".$line[$i]."\t".$temp[9];
#	    print $uniq_key."\t\t".$clade{$uniq_key}."\t".$line[$i]."\t".$temp[9]."\n";
            if ($evalue{$uniq_key} > $temp[9]){
		$sig{$uniq_key} = $temp[9]/$evalue{$uniq_key};
                #print $uniq_key."\t\t".$clade{$uniq_key}."\t".$evalue{$uniq_key}."\t".$line[$i]."\t".$temp[9]."\n";
                $evalue{$uniq_key} = $temp[9];
                $clade{$uniq_key} = $line[$i];
            }
	    elsif ($evalue{$uniq_key}/$temp[9] > $sig{$uniq_key}){
		$sig{$uniq_key} = $evalue{$uniq_key}/$temp[9];
	    }
        }
    }
    open (OUT, ">>$_[3]");
#    open (OUT1, ">>$_[4]");
    if ($_[0] =~ /\/((\w|\d)+)\./){
	$anno_clade = $1;
    }
    print OUT "$anno_clade-------------------------------------\n";
    for my $key (keys %evalue){
        print OUT $key."\t".$clade{$key}."\t".$evalue{$key}."\t";
	printf OUT "%.1e\n", $sig{$key};
#	print OUT1 $key."\t".$save_evalue{$key}."\n";
    }
    close(OUT);
#    close(OUT1);

}
