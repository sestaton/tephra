package Postprocess2;

use 5.010;
use Moose;
use MooseX::Types::Path::Class;
use Bio::SeqIO;
use File::Find;
use File::Spec;
use File::Path qw(make_path);
use File::Basename;

#with 'SeqUtils';

has fasta   => ( is => 'ro', isa => 'Path::Class::File', required => 1, coerce => 1 );
has outdir  => ( is => 'ro', isa => 'Path::Class::Dir', required => 1, coerce => 1 );
has phmmdir => ( is => 'ro', isa => 'Path::Class::Dir', required => 1, coerce => 1 );

sub postprocess2 {
    my $self = shift;
    my $dir  = $self->outdir;
    my $hmm_dir = $self->phmmdir;

    my $hmmsearch = $self->_find_hmmsearch;

    #our $hmmsearch = '/home/statonse/github/tephra/MGEScan_nonLTR_v2/hmmer-2.3.2/src/hmmsearch';

    my $genome = $self->fasta;
    my @all_clade  = ('CR1', 'I', 'Jockey', 'L1', 'L2', 'R1', 'RandI', 'Rex', 'RTE', 'Tad1', 'R2','CRE');
    my @en_clade   = ('CR1', 'I', 'Jockey', 'L1', 'L2', 'R1', 'RandI', 'Rex', 'RTE', 'Tad1');
    #my $genome="";   
    #my $dir;
    #my $hmm_dir;
    #my $domain;
    #my $seq_dir;
    #my $validation_dir;
    #my $validation_file;
    #my $evalue_file;
    my $tree_dir;
    #my $seq;
    
    #get_parameter(\$dir, \$hmm_dir);
    
    # copy seq from out2 dir into info dir
    #system("mkdir ".$dir."info/");
    #system("mkdir ".$dir."info/full");
    my $infodir = File::Spec->catdir($dir, 'info');
    my $fulldir = File::Spec->catdir($dir, 'info', 'full');
    make_path( $infodir, {verbose => 0, mode => 0771,} );
    make_path( $fulldir, {verbose => 0, mode => 0771,} );

    $self->get_full_frag($genome, $dir, \@all_clade);
    
    # get domain seq
    $self->get_domain_for_full_frag($genome, "en", \@en_clade, $dir, $hmm_dir);
    $self->get_domain_for_full_frag($genome, "rt", \@all_clade, $dir, $hmm_dir);

    # get Q value after running pHMM for EN in full elements
    my $validation_dir = File::Spec->catdir($dir, 'info', 'validation');
    make_path( $validation_dir, {verbose => 0, mode => 0771,} );
    #system("mkdir ".$validation_dir);
    
    my $domain          = 'en';
    my $seq_dir         = File::Spec->catdir($dir, 'info', 'full');
    my $validation_file = File::Spec->catfile($validation_dir, $domain);
    my $evalue_file     = File::Spec->catfile($validation_dir, $domain.'_evalue');

    my @fasfiles;
    find( sub { push @fasfiles, $File::Find::name if -f and /\.fa.*$/ }, $seq_dir );

    #opendir(DIRHANDLE, $seq_dir) || die ("Cannot open directory ".$seq_dir);
    #foreach my $name (sort readdir(DIRHANDLE)) {
	#if ($name !~ /^\./ && $name ne "R2" && $name ne "CRE" ){
    for my $file (@fasfiles) {
	my $name = basename($file);
	my $seq = File::Spec->catfile($seq_dir, $name, $name.".".$domain.'.pep');
	$self->vote_hmmsearch($seq, $hmm_dir, $domain, $validation_file, $evalue_file, \@en_clade);
    }
    #closedir(DIRHANDLE);
    
    # get Q value after running pHMM for RT in full elements
    $domain          = 'rt';
    $seq_dir         = File::Spec->catdir($dir, 'info', 'full');
    $validation_file = File::Spec->catfile($validation_dir, $domain);
    $evalue_file     = File::Spec->catfile($validation_dir, $domain.'_evalue');
    
    @fasfiles = ();
    find( sub { push @fasfiles, $File::Find::name if -f and /\.fa.*$/ }, $seq_dir );

    #opendir(DIRHANDLE, $seq_dir) || die ("Cannot open directory ".$seq_dir);
    #foreach my $name (sort readdir(DIRHANDLE)) {
    #if ($name !~ /^\./){
    for my $file (@fasfiles) {
	my $name = basename($file);
	my $seq = File::Spec->catfile($seq_dir, $name, $name.'.'.$domain.'.pep');
	$self->vote_hmmsearch($seq, $hmm_dir, $domain, $validation_file, $evalue_file, \@all_clade);
	#}
    }
    #close(DIRHANDLE);
    #system("rm ".$seq_dir."*/*.pep");
    #system("rm ".$seq_dir."*/*.rt.*");
    #system("rm ".$seq_dir."*/*.en.*");
    #system("rm -r ".$dir."b");
    #system("rm -r ".$dir."f");
}

sub get_domain_for_full_frag {
    my $self = shift;
    my ($genome, $domain, $en_clade, $dir, $hmm_dir) = @_;

    #my $genome = $_[0];
    #my $domain = $_[1];
    #my @all_clade = @$en_clade;
    #my $hmm_dir = $_[4];
    #my ($dir);

    #my ($pep_file, $dna_file);
    #my ($phmm_file, $result_pep_file, $result_dna_file);

    #for (my $i=0; $i<=$#all_clade; $i++){
    for my $clade (@$en_clade) {
	#my $clade = $all_clade[$i];

	my $resdir   = File::Spec->catdir($dir, 'info', 'full', $clade);
	my $pep_file = File::Spec->catfile($resdir, $clade.'.pep');
	my $dna_file = File::Spec->catfile($resdir, $clade.'.dna');
       
	if (-e $pep_file ){
	    my $phmm_file       = File::Spec->catfile($hmm_dir, $clade.'.'.$domain.'.hmm');
	    my $result_pep_file = File::Spec->catfile($dir, $clade.'.'.$domain.'.pe');
	    my $result_dna_file = File::Spec->catfile($dir, $clade.'.'.$domain.'.dna');
	    
	    my $flag = 2;  #1: protein-protein, 2: protein-dna 
	    $self->get_domain_pep_seq($pep_file, $phmm_file, $result_pep_file);
	    $self->get_domain_dna_seq($pep_file, $phmm_file, $result_dna_file, $dna_file, $flag);

	    my $command = "sed 's/>/>".$clade."_/' ".$result_pep_file." > ".$result_pep_file."p";
	    system($command);
	    #system("rm ".$result_pep_file);
	}
    }
}

sub get_full_frag {
    my $self = shift;
    my ($genome, $dir, $all_clade) = @_;

    #my $genome=$_[0];
    #my $dir = $_[1];
    #my @all_clade = @$clade;
    #my $file;
    #my $clade_dir;
    #my ($dna_file, $pep_file, $file_f, $file_b);


    #for (my $i=0; $i<=$#all_clade; $i++){
    for my $clade (@$all_clade) {
	# create a clade dir
	my $clade_dir = File::Spec->catdir($dir, 'info', 'full', $clade);
	my $file_f    = File::Spec->catdir($dir, 'f', 'out2', $clade.'_full');
	my $file_b    = File::Spec->catdir($dir, 'b', 'out2', $clade.'_full');

	if (-e $file_f || -e $file_b){
	    #system("mkdir ".$clade_dir);
	    make_path( $clade_dir, {verbose => 0, mode => 0771,} );
	}

	# copy full length in + strand
	if (-e $file_f) { # collate here and below
	    my $command = "cat ".$file_f." > ".$clade_dir."/".$clade.".dna";
	    system($command);
	}

	# copy full length in - strand
	if (-e $file_b) {
	    my $command = "cat ".$file_b." >> ".$clade_dir."/".$clade.".dna";
	    system($command);
	}

	# translate
	my $dna_file = File::Spec->catfile($dir, 'info', 'full', $clade, $clade.'.dna');
	my $pep_file = File::Spec->catfile($dir, 'info', 'full', $clade, $clade.'.pep');   
	if (-e $dna_file){
	    my $command = "transeq -frame=f -sequence=$dna_file -outseq=$pep_file 2>/dev/null";
	    system($command);
	}
    }
}


sub get_domain_pep_seq {
    my $self = shift;
    my ($pep_file, $phmm_file, $result_pep_file) = @_;
    #$_[0]: pep seq file
    #$_[1]: domain hmm file
    #$_[2]: output domain pep seq file 

    my $hmmsearch = $self->_find_hmmsearch;

    my $hmm_command = "$hmmsearch $phmm_file $pep_file";
    my $hmm_result = `$hmm_command`;
    my (%domain_start, %domain_end, %result_start, %result_end, %uniq_head);


    while ($hmm_result =~ /((\d|\w|\-|\_|\#|\/|\.)+\s+\d+\/\d+\s+\d+\s+\d+\s+(\[|\.)(\]|\.)\s+\d+\s+\d+\s+(\[|\.)(\]|\.)\s+(\-)*\d+\.\d+\s+((\d|\-|\.|e)+))\s*/g){
	    
	my @temp = split /\s+/, $1;
#	if ($temp[9]<0.000001 ){
	    my $key      = substr($temp[0],0,length($temp[0]));
	    my $uniq_key = substr($temp[0],0,length($temp[0])-2);

	if (not exists $uniq_head{$uniq_key}){
	    $uniq_head{$uniq_key} = 1;
	    $result_start{$key}   = $temp[2];
	    $result_end{$key}     = $temp[3];
	}
    }

    my $flag = 0;
    my $head;
    my $seq;
    open my $in, '<', $pep_file or die "\nERROR: Could not open file: $pep_file";
    #open OUT, ">$_[2]";
    open my $out, '>', $result_pep_file or die "\nERROR: Could not open file: $result_pep_file";
    while (my $each_line = <$in>){
	chomp $each_line;
	if ($each_line =~ /\>/){
	    if (defined $head && length($head) > 0 && $flag == 1){
		say $out ">".$head;
		say $out substr($seq, $result_start{$head}, eval($result_end{$head}-$result_start{$head}+1));
	    }
	    my @temp = split /\s+/, $each_line;
	    if (exists $result_start{substr($temp[0], 1, length($temp[0])-1)}) {
		$flag = 1;
		$head = substr($temp[0], 1, length($temp[0])-1);
	    }
	    else {
		$flag = 0;
	    }
	    $seq = "";
	}
	else{
	    if ($flag == 1) {
	        $seq .= $each_line;
	    }
	}
    }
    if ($flag == 1) {
	say $out ">".$head;
	say $out substr($seq, $result_start{$head}, eval($result_end{$head}-$result_start{$head}+1));
    }
    close $in;
    close $out;
}


sub get_domain_dna_seq {
    my $self = shift;
    my ($pep_file, $phmm_file, $result_dna_file, $dna_file, $flag) = @_;
    #$_[0]: pep seq file
    #$_[1]: domain hmm file
    #$_[2]: output domain dna seq file 
    #$_[3]: dna seq file

    my $hmmsearch = $self->_find_hmmsearch;
    my $hmm_command = "$hmmsearch $phmm_file $pep_file";
    my $hmm_result = `$hmm_command`;
    my (%domain_start, %domain_end, %result_start, %result_end, %uniq_head);

    while ($hmm_result =~ /((\d|\w|\-|\_|\#|\/|\.)+\s+\d+\/\d+\s+\d+\s+\d+\s+(\[|\.)(\]|\.)\s+\d+\s+\d+\s+(\[|\.)(\]|\.)\s+(\-)*\d+\.\d+\s+((\d|\-|\.|e)+))\s*/g){
	    
	my @temp = split /\s+/, $1;

	    my $key      = substr($temp[0],0,length($temp[0]));
	    my $uniq_key = substr($temp[0],0,length($temp[0])-2);

	    if (not exists $result_start{$uniq_key}){
		$uniq_head{$uniq_key} = 1;
		if ($flag == 1){
		    $result_start{$key} = $temp[2];
		    $result_end{$key} = $temp[3];
		}
		elsif($_[4] ==2){
		    $result_start{$uniq_key} = $temp[2];
		    $result_end{$uniq_key} = $temp[3];
		}
	    }
    }
    $flag = 0;
    my $head;
    my $seq;
    #open (IN, $_[3]);
    #open OUT, ">$_[2]";
    open my $in, '<', $dna_file or die "\nERROR: Could not open file: $dna_file";
    open my $out, '>', $result_dna_file or die "\nERROR: Could not open file: $result_dna_file";
    while(my $each_line = <$in>){
	chomp $each_line;
	if ($each_line =~ /\>/){
	    if (defined $head && length($head) > 0 && $flag == 1){
		say $out ">".$head;
		say $out substr($seq, $result_start{$head}*3-3, eval(($result_end{$head}-$result_start{$head}+1)*3+3));
	    }
	    my @temp = split /\s+/, $each_line;
	    if (exists $result_start{substr($temp[0], 1, length($temp[0])-1)}){
		$flag = 1;
		$head = substr($temp[0], 1, length($temp[0])-1);
	    }
	    else{
		$flag = 0;
	    }
	    $seq = "";
	}
	else {
	    if ($flag == 1){
	        $seq .= $each_line;
	    }
	}
    }
    if ($flag == 1){
	say $out ">".$head;
	say $out substr($seq, $result_start{$head}*3-3, eval(($result_end{$head}-$result_start{$head}+1)*3+3));
    }
    close $in;
    close $out;
}



sub vote_hmmsearch {
    my $self = shift;
    my ($seq, $hmm_dir, $domain, $validation_file, $evalue_file, $en_clade) = @_;

    my $hmmsearch = $self->_find_hmmsearch;
    my @line = @$en_clade;
    my %evalue;
    my %save_evalue;
    my %clade;
    my %sig;
    my $i;
    my $anno_clade; 

    #open (IN, $_[0]);
    open my $in, '<', $seq or die "\nERROR: Could not open file: $seq";
    while (my $each_line = <$in>){
	if ($each_line =~ /\>/){
	    chomp $each_line;
	    my $uniq_key = substr($each_line,1,length($each_line)-1);
	    $evalue{$uniq_key}      = 1000;
	    $save_evalue{$uniq_key} = 1000;
	    $clade{$uniq_key}       = "-";
	    $sig{$uniq_key}         = 1;
	}
    }
    close $in;

    for ($i=0; $i<=$#line; $i++){

        my $command = "$hmmsearch ".$hmm_dir."/".$line[$i].".".$domain.".hmm ".$seq;
        my $hmm_result = `$command`;

        while ($hmm_result =~ /((\d|\w|\-|\_|\#|\/|\.)+\s+\d+\/\d+\s+\d+\s+\d+\s+(\[|\.)(\]|\.)\s+\d+\s+\d+\s+(\[|\.)(\]|\.)\s+\-*\d+\.\d+\s+((\d|\-|\.|e)+))\s*/g){
            my @temp = split /\s+/, $1;
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
    #open (OUT, ">>$_[3]");
    open my $out, '>>', $validation_file or die "\nERROR: Could not open file: $validation_file";
#    open (OUT1, ">>$_[4]");
    if ($seq =~ /\/((\w|\d)+)\./){
	$anno_clade = $1;
    }
    say $out "$anno_clade-------------------------------------";
    for my $key (keys %evalue){
        print $out join "\t", $key, $clade{$key}, $evalue{$key}."\t";
	printf $out "%.1e\n", $sig{$key};
#	print OUT1 $key."\t".$save_evalue{$key}."\n";
    }
    close $out;
#    close(OUT1);

}

sub _find_hmmsearch {
    my $self = shift;

    my $hmmsearch = '/home/statonse/github/tephra/MGEScan_nonLTR_v2/hmmer-2.3.2/src/hmmsearch';
    return $hmmsearch;
}

__PACKAGE__->meta->make_immutable;

1;
