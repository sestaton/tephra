package RunHMM;

use 5.010;
use Moose;
use MooseX::Types::Path::Class;
use File::Path qw(make_path);
use File::Basename;
use File::Spec;
use Bio::SeqIO;
#with 'SeqUtils';

has fasta   => ( is => 'ro', isa => 'Path::Class::File', required => 1, coerce => 1 );
has outdir  => ( is => 'ro', isa => 'Path::Class::Dir', required => 1, coerce => 1 );
has phmmdir => ( is => 'ro', isa => 'Path::Class::Dir', required => 1, coerce => 1 );
has pdir    => ( is => 'ro', isa => 'Path::Class::Dir', required => 1, coerce => 1 );

sub run_mgescan {
    my $self = shift;
    my $dna_file = $self->fasta;
    my $out_dir  = $self->outdir;
    my $phmm_dir = $self->phmmdir;
    my $pdir     = $self->pdir;

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
    }
    
    ##########################################################
    # get signal for some state of ORF1, RT, and APE
    # need HMMSEARCH
    ##########################################################
    print "Getting signal...\n";
    #my ($phmm_file, $domain_rt_pos_file, $domain_ape_pos_file, $domain_orf1_pos_file);

    print "    Protein sequence...\n";
    my $pep_file = File::Spec->catfile($out_dir, $dna_name.".pep");
    #$command = $pdir."translate -d ".$dna_file." -h ".$dna_name." -p ".$pep_file;
    #system($command);
    $self->translate_forward($dna_file, $pep_file);

    print "    RT signal...\n";
    my $phmm_file = File::Spec->catfile($phmm_dir, "ebi_ds36752_seq.hmm");
    my $domain_rt_pos_file = File::Spec->catfile($pos_dir, $dna_name.".rt.pos");
    $self->get_signal_domain($pep_file, $phmm_file, $domain_rt_pos_file);
    
    print "    APE signal...\n";
    $phmm_file = File::Spec->catfile($phmm_dir, "ebi_ds36736_seq.hmm");
    my $domain_ape_pos_file = File::Spec->catfile($pos_dir, $dna_name.".ape.pos");
    $self->get_signal_domain($pep_file, $phmm_file, $domain_ape_pos_file);
    
    ##############################################################################
    # generate corresponsing empty domains files if either of them does not exist 
    ##############################################################################
    if (-e $domain_rt_pos_file  || -e $domain_ape_pos_file ){
	print $dna_name."\n";
	
	if (! -e $domain_rt_pos_file){
	    open my $out, '>', $domain_rt_pos_file;
	    print $out "";
	    close $out;
	}
	elsif (! -e $domain_ape_pos_file){
	    open my $out, '>', $domain_ape_pos_file;
	    print $out "";
	    close $out;
	}
	
	#$command = $pdir."match_pos.pl -rt=".$domain_rt_pos_file." -ape=".$domain_ape_pos_file;
	#system($command);
	###########################################################
	# run hmm
	###########################################################
	#print "Running HMM...\n";

	my $mgescan  = File::Spec->catfile($pdir, 'hmm', 'MGEScan');
	my $out_file = File::Spec->catfile($out1_dir, $dna_name);
	my $chrhmm   = File::Spec->catfile($pdir, 'hmm', 'chr.hmm');
	my $ldir = $pdir."/";
	my $command = "$mgescan -m $chrhmm -s $dna_file -r $domain_rt_pos_file -a $domain_ape_pos_file -o $out_file -p $ldir -d $out1_dir";
	#print $command."\n";
	system($command); 
    }

 
    #if (-e $pep_file){
        #system("rm ".$pep_file);
    #} 
}

sub translate_forward {
    my $self = shift;
    my ($in, $out) = @_;

    say STDERR join q{ }, $in, $out;
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
 
   #return $out;
}

sub get_signal_domain {
    my $self = shift;
    my ($pep_file, $phmm_file, $domain_rt_pos_file) = @_;
    #$_[0]: pep seq file
    #$_[1]: domain hmm file
    #$_[2]: output domain dna position file

    my $hmmsearch = $self->_find_hmmsearch;
    my $hmm_command = "$hmmsearch  -E 0.00001 $phmm_file $pep_file";
    my $hmm_result = `$hmm_command`;
    my %domain_start=();
    my %domain_end=();
    my $evalue;
    my $temp_file   =  $domain_rt_pos_file."temp";
    my $temp_file2  =  $domain_rt_pos_file."temp2";
    my $output_file =  $domain_rt_pos_file;

    # run hmmsearch to find the domain and save it in the temprary file    
    open my $out, '>', $temp_file;
    while ($hmm_result =~ /((\S)+\s+\d+\/\d+\s+\d+\s+\d+\s+(\[|\.)(\]|\.)\s+\d+\s+\d+\s+(\[|\.)(\]|\.)\s+(-)*\d+\.\d+\s+((\d|\-|\.|e)+))\s*/g){

        my @temp = split /\s+/, $1;
	if ($temp[9]< 0.001 ){
	    say $out join "\t", eval($temp[2]*3), eval($temp[3]*3), @temp[5..6], @temp[8..9];
	}
    }
    close $out;
    if (-s $temp_file > 0){

        system("sort +0 -1n ".$temp_file." > ".$temp_file2);

        my $start = -1;
        my $end = -1;
        my @pre = (-1000, -1000, -1000, -1000, -1000, -1000);
        open my $in, '<', $temp_file2;
        open my $out, '>', $output_file;

        while(my $each_line = <$in>){
            my @temp = split /\s+/, $each_line;

            if ($temp[0] - $pre[1] < 300 ) {
                $end = $temp[1];
                $evalue = $evalue * $temp[5];
            }
	    else{
                if ($start >= 0 && $evalue < 0.00001){
                    say $out join "\t", $start, $end, @pre[4..5];
                }
                $start  = $temp[0];
                $end    = $temp[1];
                $evalue = $temp[5];
            }
            @pre = @temp;
        }

        if ($start >= 0 && $evalue < 0.00001){ 
            say $out join "\t", $start, $end, @pre[4..5];
        }
        close $in;
        close $out;
        #system("rm ".$temp_file2);
	unlink $temp_file2;
    }
    #system("rm ".$temp_file);
    unlink $temp_file;
}

#sub get_id {
#    my $self = shift;
#    my @temp = split(/\//, ${$_[0]});
#    ${$_[1]} = $temp[$#temp];
#}


#sub usage {
#    die "Usage: run_hmm.pl --dna=<dna_file_path>  --out=<output_dir> --pdir=<program dir>";
#}


#sub get_parameter{

#    my ($dna, $pep, $out, $phmm, $pdir);

#    GetOptions(
#               'dna=s' => \$dna,
#               'out=s' => \$out,
#	       'phmm=s' => \$phmm,
#	       'pdir=s' => \$pdir,
#               );

#    if (! -e $dna){
#        print "ERROR: The file $dna does not exist.\n";
#        usage();
#    }
#    if (! -e $phmm){
#        print "ERROR: The file $phmm does not exist.\n";
#        usage();
#    }
#    if (! -e $pdir){
#        print "ERROR: The file $phmm does not exist.\n";
#        usage();
#    }
#    if (! -d $out){
#	system("mkdir ".$out);
#    }

#    ${$_[0]} = $dna;
#    ${$_[1]} = $out;
#    ${$_[2]} = $phmm;
#    ${$_[3]} = $pdir;
#}

sub _find_hmmsearch {
    my $self = shift;

    my $hmmsearch = '/home/statonse/github/tephra/MGEScan_nonLTR_v2/hmmer-2.3.2/src/hmmsearch';
    return $hmmsearch;
}

__PACKAGE__->meta->make_immutable;

1;
