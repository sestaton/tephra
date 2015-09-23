package Tephra::LTR::LTRStats;

use 5.010;
use Moose;
use MooseX::Types::Path::Class;
use Statistics::Descriptive;
use Sort::Naturally;
use List::MoreUtils qw(indexes any);
use File::Spec;
use File::Find;
use File::Basename;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Tools::GFF;
use Time::HiRes qw(gettimeofday);
use Parallel::ForkManager;
use Cwd;
use Try::Tiny;
use namespace::autoclean;

use Data::Dump;
use Data::Printer;

with 'Tephra::Role::GFF',
     'Tephra::Role::Util';

has genome => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has outdir => (
    is       => 'ro',
    isa      => 'Path::Class::Dir',
    required => 1,
    coerce   => 1,
);

has threads => (
    is        => 'ro',
    isa       => 'Int',
    predicate => 'has_threads',
    lazy      => 1,
    default   => 1,
);

has gff => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

#
# methods
#
sub collect_ltr_features {
    my $self = shift;
    my $fasta = $self->genome;
    my $dir   = $self->outdir;
    my $gff   = $self->gff;
    
    #my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    #my $comp = File::Spec->catfile($dir, $name."_complete.fasta");
    #my $five_pr_ltrs  = File::Spec->catfile($dir, $name."_5prime-ltrs.fasta");
    #my $three_pr_ltrs = File::Spec->catfile($dir, $name."_3prime-ltrs.fasta");

    #open my $allfh, '>>', $comp;
    #open my $fivefh, '>>', $five_pr_ltrs;
    #open my $threfh, '>>', $three_pr_ltrs;

    my $gffio = Bio::Tools::GFF->new( -file => $gff, -gff_version => 3 );

    my ($start, $end, $elem_id, $key, %feature, %ltrs, %seen);
    while (my $feature = $gffio->next_feature()) {
	if ($feature->primary_tag eq 'LTR_retrotransposon') {
	    my @string = split /\t/, $feature->gff_string;
	    ($elem_id) = ($string[8] =~ /ID=?\s+?(LTR_retrotransposon\d+)/);
	    ($start, $end) = ($feature->start, $feature->end);
	    $key = join ".", $elem_id, $start, $end;
	    #$ltrs{$key}{'full'} = join "-", $string[0], $feature->primary_tag, @string[3..4];
	}
	next unless defined $start && defined $end;
	if ($feature->primary_tag eq 'long_terminal_repeat') {
	    my @string = split /\t/, $feature->gff_string;
	    if ($feature->start >= $start && $feature->end <= $end) {
		my $ltrkey = join "||", $string[0], $feature->primary_tag, @string[3..4];
		push @{$ltrs{$key}{'ltrs'}}, $ltrkey unless exists $seen{$ltrkey};
		$seen{$ltrkey} = 1;
	    }
	}
    }

    #p %ltrs and exit;
    
    my %pdoms;
    my $ltrct = 0;
    for my $ltr (sort keys %ltrs) {
	my ($element, $rstart, $rend) = split /\./, $ltr;
	for my $ltr_repeat (@{$ltrs{$ltr}{'ltrs'}}) {
	    my ($src, $ltrtag, $s, $e) = split /\|\|/, $ltr_repeat;
	    my $ltrs_out = File::Spec->catfile($dir, $ltr."_ltrs.fasta");
	    open my $ltrs_outfh, '>>', $ltrs_out;
	    if ($ltrct) {
		my $fiveprime_tmp = File::Spec->catfile($dir, $ltr."_5prime-ltr.fasta");
		$self->subseq($fasta, $src, $element, $s, $e, $fiveprime_tmp, $ltrs_outfh);
		$ltrct = 0;
	    }
	    else {
		my $threeprime_tmp = File::Spec->catfile($dir, $ltr."_3prime-ltr.fasta");
		$self->subseq($fasta, $src, $element, $s, $e, $threeprime_tmp, $ltrs_outfh);
		$ltrct++;
	    }
	    close $ltrs_outfh;
	    #$ltrct = 0;
	}
	#$ltrct = 0;
    }
    #close $allfh;
    #close $fivefh;
    #close $threfh;

    #for my $file ($comp, $five_pr_ltrs, $three_pr_ltrs) {
	#unlink $file if ! -s $file;
    #}
}

sub collect_feature_args {
    my $self = shift;
    my $dir = $self->outdir;

    my (@ltrs, %aln_args);
    find( sub { push @ltrs, $File::Find::name if -f and /ltrs.fasta$/ }, $dir);

    # ltr
    #my $args = "muscle -in -out ";
    $aln_args{ltrs} = { seqs => \@ltrs };

    return \%aln_args;
}

sub align_features {
    my $self = shift;
    my $dir  = $self->outdir;
    my $threads = $self->threads;

    my $args = $self->collect_feature_args;
    #$self->_remove_singletons($args);

    #use Data::Dump;
    #say STDERR "DEBUG: ";
    #dd $args;

    my $t0 = gettimeofday();
    my $doms = 0;
    my %reports;
    my $outfile = File::Spec->catfile($dir, 'all_aln_reports.txt');
    my $logfile = File::Spec->catfile($dir, 'all_aln_reports.log');
    open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n";
    open my $log, '>>', $logfile or die "\nERROR: Could not open file: $logfile\n";
    
    my $pm = Parallel::ForkManager->new($threads);
    $pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
			      for my $bl (sort keys %$data_ref) {
				  open my $report, '<', $bl or die "\nERROR: Could not open file: $bl\n";
				  print $out $_ while <$report>;
				  close $report;
				  unlink $bl;
			      }
			      my $t1 = gettimeofday();
			      my $elapsed = $t1 - $t0;
			      my $time = sprintf("%.2f",$elapsed/60);
			      say $log basename($ident),
			      " just finished with PID $pid and exit code: $exit_code in $time minutes";
			} );

    for my $type (keys %$args) {
	for my $db (@{$args->{$type}{seqs}}) {
	    $doms++;
	    $pm->start($db) and next;
	    my $mrep = $self->process_aln_args($db);
	    $reports{$mrep} = 1;

	    $pm->finish(0, \%reports);
	}
    }

    $pm->wait_all_children;
    close $out;

    my $t2 = gettimeofday();
    my $total_elapsed = $t2 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);

    say $log "\n========> Finished running PAML on $doms LTRs in $final_time minutes";
    close $log;

    return $outfile;
}

sub run_baseml {
    my $self = shift;
    
    my $name     = shift;
    my $divfile  = $name."-divergence.txt";
    my $outfile  = $name."-paml.out";
    my $phylip   = $name.".phy";
    my $treefile = $name.".dnd";
    
    my $ctl_file = "      seqfile = $phylip 
     treefile = $treefile

      outfile = $outfile       * main result file
        noisy = 0   * 0,1,2,3: how much rubbish on the screen
      verbose = 0   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI 

        model = 1   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
                    * 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu
        
        Mgene = 0   * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff

*        ndata = 5
        clock = 0   * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
    fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below
        kappa = 5  * initial or fixed kappa

    fix_alpha = 0   * 0: estimate alpha; 1: fix alpha at value below
        alpha = 0.5   * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * 1: different alpha's for genes, 0: one alpha
        ncatG = 5   * # of categories in the dG, AdG, or nparK models of rates
        nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK 

        nhomo = 0   * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1   * (0,1,2): rates (alpha>0) or ancestral states

   Small_Diff = 7e-6
    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
*        icode = 0  * (with RateAncestor=1. try GC in data,model=4,Mgene=4)
*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed
       method = 0  * Optimization method 0: simultaneous; 1: one branch a time";

    my $control_file = "baseml.ctl";
    open my $out, '>', $control_file or die "\nERROR: Could not open file: $control_file\n";
    
    print $out $ctl_file;
    close $out;
    my $basemlcmd = "baseml 2>&1 > /dev/null";
    $self->run_cmd($basemlcmd);
	
    #my $pwdivfile = "2base.t";     # using the output now to get kappa and divergence
    open my $divin, '<', $outfile or die "ERROR: Could not open divergence file: $outfile\n";
    open my $divout, '>', $divfile or die "ERROR: Could not open divergence file: $divfile\n";
    
    while (<$divin>) {
	chomp;
	if (/^\w+\s+\d\.\d+\(\s/) {
	    # 3prime_Ung         0.0269( 8.7752)
	    my ($seqid,$divergence_time,$kappa) = split /\s+/;
	    
	    $divergence_time =~ s/\($//;
	    $kappa =~ s/\)$//;
	    my $time = $divergence_time/(1e-8 * 2);   # T=k/2r, k=1.0 10-8
	    
	    # alignID divergence age Ts:Tv
	    say $divout join "\t", $phylip,$divergence_time,$time,$kappa;
	}
	elsif (/^(\d\w+)         (\d\.\d+\()(\d+\.\d+\))/) {
	    # 3prime_RL1         0.0087(999.0000)
	    
	    my $divergence_time = $2;
	    my $kappa = $3;
	    $divergence_time =~ s/\($//;
	    $kappa =~ s/\)$//;
	    my $time = $divergence_time/(1e-8 * 2);
	    
	    # alignID divergence age Ts:Tv
	    say $divout join "\t", $phylip,$divergence_time,$time,$kappa;
	}
	elsif (/^(\d\w+)         (\d\.\d+\()(\-\d+\.\d+\))/) {
	    # 3prime_RL1         0.0017(-0.0025)
	    
	    my $divergence_time = $2;
	    my $kappa = $3;
	    $divergence_time =~ s/\($//;
	    $kappa =~ s/\)$//;
	    my $time = $divergence_time/(1e-8 * 2);
	    
	    # alignID divergence age Ts:Tv
	    say $divout join "\t", $phylip,$divergence_time,$time,$kappa;
	}
    }
    close $divin;
    close $divout;
    my $summary_file_path = "Divergence_time_files";
    copy $divfile, $summary_file_path or die "\nERROR: Copy failed: $!";
    #system("cp $divfile $summary_file_path");
    unlink $control_file;           # instead of overwriting, delete so the last one is not left

    if ($self->clean) {
	# remove the PAML output but keep the summary produced
	# by this script in a separate directory.
	unlink "2base.t", "rub", "rst", "rst1", "lnf", "rates", "in.basemlg", $divfile, $outfile;
    }
}

sub process_align_args {
    my $self = shift;
    my ($db) = @_;

    my ($name, $path, $suffix) = fileparse($db, qr/\.[^.]*/);
    my $aln = File::Spec->catfile($path, $name."_clustal-out.aln");
    my $log = File::Spec->catfile($path, $name."_clustal-out.log");

    my $clwcmd  = "clustalw -infile=$db -outfile=aln 2>$log";
    $self->run_cmd($clwcmd);

    my $phy = $self->parse_aln($aln);
    return $phy;
}

sub parse_aln {
    my $self = shift;
    my ($aln) = @_;

    my ($name, $path, $suffix) = fileparse($aln, qr/\.[^.]*/);
    my $phy = File::Spec->catfile($path, $name.".phy");
    
    my $aln_in = Bio::AlignIO->new(-file   => $aln,
				   -format => 'clustalw');

    my $aln_out = Bio::AlignIO->new(-file  => $phy,
				    -format => 'phylip',
				    -flag_SI => 1);

    while (my $aln = $aln_in->next_aln) {
	$aln_out->write_aln($aln);
    }
    
    return $phy;
}

sub subseq {
    my $self = shift;
    my ($fasta, $loc, $elem, $start, $end, $tmp, $out) = @_;
    my $cmd = "samtools faidx $fasta $loc:$start-$end > $tmp";
    $self->run_cmd($cmd);

    my $id = join "_", $loc, $elem, "$start-$end";
    if (-s $tmp) {
	my $seqio = Bio::SeqIO->new( -file => $tmp, -format => 'fasta' );
	while (my $seqobj = $seqio->next_seq) {
	    my $seq = $seqobj->seq;
	    if ($seq) {
		$seq =~ s/.{60}\K/\n/g;
		say $out join "\n", ">".$id, $seq;
	    }
	}
    }
    unlink $tmp;
}

sub _remove_singletons {
    my $self = shift;
    my ($args) = @_;

    my @singles;
    my ($index, $seqct) = (0, 0);
    for my $type (keys %$args) {
	delete $args->{$type} if ! @{$args->{$type}{seqs}};
	for my $db (@{$args->{$type}{seqs}}) {
	    my $seqio = Bio::SeqIO->new( -file => $db, -format => 'fasta' );
	    while (my $seqobj = $seqio->next_seq) { $seqct++ if defined $seqobj->seq; }
	    if ($seqct < 2) {
		push @singles, $index;
		unlink $db;
	    }
	    $index++;
	    $seqct = 0;
	}

	if (@{$args->{$type}{seqs}}) {
	    if (@singles > 1) {
		for (@singles) {
		    splice @{$args->{$type}{seqs}}, $_, 1;
		    @singles = map { $_ - 1 } @singles; # array length is changing after splice so we need to adjust offsets
		}
	    }
	    else {
		splice @{$args->{$type}{seqs}}, $_, 1 for @singles;
	    }
	}
	else {
	    delete $args->{$type};
	}
	$index = 0;
	@singles = ();
    }
}
	
__PACKAGE__->meta->make_immutable;

1;
