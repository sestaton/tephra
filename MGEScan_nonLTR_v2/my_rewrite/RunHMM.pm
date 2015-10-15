package RunHMM;

use Moose;
use MooseX::Path::Class::File;
use File::Path qw(make_path);
use File::Basename;

has dna  => ( is => 'ro', isa => 'Path::Class::File', required => 1, coerce => 1 );
has out  => ( is => 'ro', isa => 'Path::Class::Dir', required => 1, coerce => 1 );
has phmm => ( is => 'ro', isa => 'Path::Class::Dir', required => 1, coerce => 1 );
has pdir => ( is => 'ro', isa => 'Path::Class::Dir', required => 1, coerce => 1 );
    
##########################################################
# get input parameter of dna file, pep file, output dir
##########################################################
#print "Getting input parameter...\n";
#my ($dna_file, $pep_file, $out_dir, $dna_name, $command);
#my ($out1_dir, $out_file, $pos_dir, $phmm_dir, $pdir);
#get_parameter(\$dna_file, \$out_dir, \$phmm_dir, \$pdir);
my ($dna_name, $dna_path, $dna_suffix) = fileparse($dna_file, qr/\.[^.]*/);
#get_id(\$dna_file, \$dna_name);

my $out1_dir = File::Spec->catdir($out_dir, "out1");
my $pos_dir  = File::Spec->catdir($out_dir, "pos");
unless ( -d $out1_dir ) {
    make_path( $out1_dir, {verbose => 0, mode => 0771,} );
}

unless ( -d $pos_dir ) {
    make_path( $pos_dir, {verbose => 0, mode => 0771,} );
}}

##########################################################
# get signal for some state of ORF1, RT, and APE
# need HMMSEARCH
##########################################################
print "Getting signal...\n";
my ($phmm_file, $domain_rt_pos_file, $domain_ape_pos_file, $domain_orf1_pos_file);

print "    Protein sequence...\n";
#$pep_file = $out_dir.$dna_name.".pep";
#$command = $pdir."translate -d ".$dna_file." -h ".$dna_name." -p ".$pep_file;
#system($command);
translate_forward($dna_file, $pep_file);

print "    RT signal...\n";
$phmm_file = File::Spec->catfile($self->phmm, "ebi_ds36752_seq.hmm");
$domain_rt_pos_file = $pos_dir.$dna_name.".rt.pos";
get_signal_domain(\$pep_file, \$phmm_file, \$domain_rt_pos_file);

print "    APE signal...\n";
$phmm_file = File::Spec->catfile($self->phmm, "ebi_ds36736_seq.hmm");
$domain_ape_pos_file = $pos_dir.$dna_name.".ape.pos";
get_signal_domain(\$pep_file, \$phmm_file, \$domain_ape_pos_file);

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
sub translate_forward {
    my ($in, $out) = @_;

    my $seqin  = Bio::SeqIO->new(-file => $in, -format => 'fasta');
    my $seqout = Bio::SeqIO->new(-file => ">>$out", -format => 'fasta');

    my $seqobj = $seqin->next_seq;
    for my $frame (0..2) { ## forward 3 frames
	my $prot_obj = $seqobj->translate(-frame => $frame);
	my $id = $prot_obj->id;
	$id .= "_$frame";
	$prot_obj->id($id);
	$seqout->write_seq($prot_obj);
    }
 
   return $out;
}

sub get_signal_domain{

    #$_[0]: pep seq file
    #$_[1]: domain hmm file
    #$_[2]: output domain dna position file

    my $hmm_command = "hmmsearch  -E 0.00001 ".${$_[1]}." ".${$_[0]};
    my $hmm_result = `$hmm_command`;
    my %domain_start=();
    my %domain_end=();
    my $evalue;
    my $temp_file = ${$_[2]}."temp";
    my $temp_file2 = ${$_[2]}."temp2";
    my $output_file = ${$_[2]};

    # run hmmsearch to find the domain and save it in the temprary file    
    open (OUT, ">$temp_file");
    while ($hmm_result =~ /((\S)+\s+\d+\/\d+\s+\d+\s+\d+\s+(\[|\.)(\]|\.)\s+\d+\s+\d+\s+(\[|\.)(\]|\.)\s+(-)*\d+\.\d+\s+((\d|\-|\.|e)+))\s*/g){

        my @temp = split(/\s+/, $1);
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
        system("rm ".$temp_file2);
    }
    system("rm ".$temp_file);
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

