package RunHMM;

use 5.010;
use Moose;
use MooseX::Types::Path::Class;
use autodie;
use File::Path          qw(make_path);
use IPC::System::Simple qw(capture EXIT_ANY);
use File::Basename;
use File::Spec;
use File::Find;
use Bio::SeqIO;
use Bio::SearchIO;
use Parallel::ForkManager;
use Try::Tiny;
use Data::Dump;

with 'PathUtils';

has fasta   => ( is => 'ro', isa => 'Path::Class::File', required => 1, coerce => 1 );
has outdir  => ( is => 'ro', isa => 'Path::Class::Dir',  required => 1, coerce => 1 );
has phmmdir => ( is => 'ro', isa => 'Path::Class::Dir',  required => 1, coerce => 1 );
has pdir    => ( is => 'ro', isa => 'Path::Class::Dir',  required => 1, coerce => 1 );

sub run_mgescan {
    my $self = shift;
    my $dna_file = $self->fasta;
    my $out_dir  = $self->outdir;
    my $phmm_dir = $self->phmmdir;
    my $pdir     = $self->pdir;

    my ($dna_name, $dna_path, $dna_suffix) = fileparse($dna_file, qr/\.[^.]*/);
    my $outf_dir = File::Spec->catdir($out_dir, "out1");
    my $pos_dir  = File::Spec->catdir($out_dir, "pos");
    unless ( -d $outf_dir ) {
	make_path( $outf_dir, {verbose => 0, mode => 0771,} );
    }
    
    unless ( -d $pos_dir ) {
	make_path( $pos_dir, {verbose => 0, mode => 0771,} );
    }
    
    # get signal for some state of ORF1, RT, and APE
    print "Getting signal...\n";
    print "    Protein sequence...\n";
    my $pep_file = File::Spec->catfile($out_dir, $dna_name.$dna_suffix.".pep");
    $self->translate_forward($dna_file, $pep_file);

    print "    RT signal...\n";
    my $phmm_file = File::Spec->catfile($phmm_dir, "ebi_ds36752_seq.hmm");
    my $domain_rt_pos_file = File::Spec->catfile($pos_dir, $dna_name.$dna_suffix.".rt.pos");
    $self->get_signal_domain($pep_file, $phmm_file, $domain_rt_pos_file);
    
    print "    APE signal...\n";
    $phmm_file = File::Spec->catfile($phmm_dir, "ebi_ds36736_seq.hmm");
    my $domain_ape_pos_file = File::Spec->catfile($pos_dir, $dna_name.$dna_suffix.".ape.pos");
    $self->get_signal_domain($pep_file, $phmm_file, $domain_ape_pos_file);
    
    # generate corresponsing empty domains files if either of them does not exist 
    if (-e $domain_rt_pos_file || -e $domain_ape_pos_file ){
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

	# run hmm
	print "Running HMM...\n";

	my $mgescan  = File::Spec->catfile($pdir, 'hmm', 'MGEScan');
	my $out_file = File::Spec->catfile($outf_dir, $dna_name.$dna_suffix);
	my $chrhmm   = File::Spec->catfile($pdir, 'hmm', 'chr.hmm');
	my $ldir = $pdir."/";
	$outf_dir .= "/";
	#try {
	    system("$mgescan -m $chrhmm -s $dna_file -r $domain_rt_pos_file -a $domain_ape_pos_file -o $out_file -p $ldir -d $outf_dir");
	    #say STDERR "cmd: $command";
	    #system($command);
	#}
	#catch {
	    #say "mgescan died with: $_";
	#};
    }
 
    #if (-e $pep_file){
        #system("rm ".$pep_file);
    #} 
}

sub translate_forward {
    my $self = shift;
    my ($in, $out) = @_;

    my @parts;
    my ($name, $path, $suffix) = fileparse($out, qr/\.[^.]*/);
    my $seqin = Bio::SeqIO->new(-file => $in, -format => 'fasta');
    
    my $pm = Parallel::ForkManager->new(3);
    my $seqobj = $seqin->next_seq;
    for my $frame (0..2) { ## forward 3 frames
	$pm->start($frame) and next;
	my $outpart    = $out."_frame$frame";
	my $seqoutpart = Bio::SeqIO->new(-file => ">$outpart", -format => 'fasta');
	my $prot_obj   = $seqobj->translate(-frame => $frame);
	#my $id = $prot_obj->id;
	#$id .= "_$frame";
	my $frameid = $frame+1;              # this is a hack to get mgscan working
	my $id = basename($in)."_$frameid";  # //
	$prot_obj->id($id);
	$seqoutpart->write_seq($prot_obj);
	$pm->finish(0);
    }
    $pm->wait_all_children;
    
    find( sub { push @parts, $File::Find::name if -f and /frame[012]$/ }, $path);
    open my $seqout, '>>', $out or die "\nERROR: Could not open file: $out\n";
    for my $part (sort @parts) {
	say "writing $part ...";
	my $lines = do { 
	    local $/ = undef; 
	    open my $fh_in, '<', $part or die "\nERROR: Could not open file: $part\n";
	    <$fh_in>;
	};
	chomp $lines;
	say $seqout $lines;
	unlink $part;
    }
    close $seqout;
}

sub get_signal_domain {
    my $self = shift;
    my ($pep_file, $phmm_file, $domain_rt_pos_file) = @_;

    # debugging
    my ($pname, $ppath, $psuffix) = fileparse($phmm_file, qr/\.[^.]*/);
    my $signal_out = $pep_file."_".$pname."_signal_searchout.txt";

    my %domain_start;
    my %domain_end;
    my %domain_pos;
    my $evalue;
    my $temp_file   =  $domain_rt_pos_file."temp";
    my $stemp_file  =  $domain_rt_pos_file."temp_sorted";
    my $output_file =  $domain_rt_pos_file;

    # run hmmsearch to find the domain and save it in the temprary file
    my $hmmsearch   = $self->find_hmmsearch;
    my @hmm_results = capture([0..5], $hmmsearch, "-E", "0.00001", $phmm_file, $pep_file);
    $self->_parse_hmmsearch(\@hmm_results, $signal_out, $temp_file);

    if (-s $temp_file) {	
	#say "DEBUG temp_file: $temp_file";
        #system("sort +0 -1n ".$temp_file." > ".$temp_file2);
	$self->_sort_matches($temp_file, $stemp_file);
        my ($start, $end) = (-1, -1);
        my @pre = (-1000, -1000, -1000, -1000, -1000, -1000);
        open my $in, '<', $stemp_file or die "\nERROR: Could not open file: $stemp_file\n";
        open my $out, '>', $output_file or die "\nERROR: Could not open file: $output_file\n";

        while (my $each_line = <$in>) {
	    chomp $each_line;
            my @temp = split /\t/, $each_line;
            if ($temp[0] - $pre[1] < 300 ) {
                $end = $temp[1];
                $evalue = $evalue * $temp[5];
            }
	    else {
                if ($start >= 0 && $evalue < 0.00001) {
                    say $out join "\t", $start, $end, @pre[4..5];
                }
                ($start, $end, $evalue) = @temp[0,1,5];
            }
            @pre = @temp;
        }

        if ($start >= 0 && $evalue < 0.00001) { 
            say $out join "\t", $start, $end, @pre[4..5];
        }
        close $in;
        close $out;
	#unlink $stemp_file;
    }
    #unlink $temp_file;
}

sub _parse_hmmsearch {
    my $self = shift;
    my ($hmm_results, $signal_out, $outfile) = @_;

    open my $o, '>', $signal_out or die "\nERROR: Could not open file: $signal_out\n";
    print $o @$hmm_results;
    close $o;

    if (-s $signal_out) {
	say "debug hmmer file: $signal_out";
	open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";
	my $hmmer_in = Bio::SearchIO->new(-file => $signal_out, -format => 'hmmer');
	
	my @evalues;
	while ( my $result = $hmmer_in->next_result ) {    
	    #my $query = $result->query_name;
	    while ( my $hit = $result->next_hit ) {
		#my $hitid  = $hit->name;
		#my $score  = $hit->raw_score;
		#my $signif = $hit->significance;
		#my $bits   = $hit->bits;
		while ( my $hsp = $hit->next_hsp ) {
		    my $hstart = $hsp->start('hit');
		    my $hstop  = $hsp->end('hit');
		    my $qstart = $hsp->start('query');
		    my $qstop  = $hsp->end('query');
		    my $score  = $hsp->score;
		    my $e_val  = $hsp->evalue;
		    #push @evalues, $e_val;
		    #say join q{ }, $hitid, $hstart, $hstop, $qstart, $qstop, $score, $e_val;
		    #say $out join "\t", eval($temp[2]*3), eval($temp[3]*3), @temp[5..6], @temp[8..9]; 
		    say $out join "\t", eval($hstart*3), eval($hstop*3), $qstart, $qstop, $score, $e_val;
		}
	    }
	}
	close $out;
    }
    #unlink $signal_out;
}

sub _sort_matches {
    my $self = shift;
    my ($unsorted, $sorted) = @_;
    my %hash;

    open my $in, '<', $unsorted or die "\nERROR: Could not open file: $unsorted\n";
    open my $out, '>', $sorted or die "\nERROR: Could not open file: $sorted\n";
    while (my $l = <$in>) {
	chomp $l;
	my @f = split /\t/, $l;
	my $start = shift @f;
	$hash{$start} = join "||", @f;
    }
    close $in;

    for my $sk (sort { $a <=> $b } keys %hash) {
	my @doms = split /\|\|/, $hash{$sk};
	say $out join "\t", $sk, @doms;
    }
    close $out;
}


__PACKAGE__->meta->make_immutable;

1;
