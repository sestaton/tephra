package Tephra::Role::Run::PAML;

use 5.010;
use Moose::Role;
use MooseX::Types::Path::Class;
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use Try::Tiny;
use File::Copy qw(move copy);
use File::Spec;
use File::Find;
use File::Basename;
use Log::Any        qw($log);
use Cwd;
use namespace::autoclean;

has baseml_exec => (
    is        => 'rw',
    isa       => 'Path::Class::File',
    required  => 0,
    coerce    => 1, 
    reader    => 'get_baseml_exec',
    writer    => 'set_baseml_exec',
    predicate => 'has_baseml_exec',
    builder   => '_build_baseml_exec',
);

has genome => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has clean => (
    is       => 'ro',
    isa      => 'Bool',
    required => 0,
    default  => 0,
);

sub run_baseml {
    my $self = shift;
    my ($args) = @_;
    my $wd = $args->{working_dir};
    chdir $wd or die "\nERROR: Could not change directories: $!";
    
    my $baseml = $self->get_baseml_exec;
    my ($stdout, $stderr, $exit);
    try {
	my @out = capture { system([0..5], $baseml) };
    }
    catch {
	$log->error("baseml failed. Here is the exception: $_\nExiting.");
	exit(1);
    };
    return;
}

sub create_baseml_files {
    my $self = shift;
    my ($args) = @_;
    my $phy = $args->{phylip};
    my $dnd = $args->{treefile};
    my $seqfile = basename($phy);
    my $treefile = basename($dnd);
    
    my $cwd = getcwd();
    my ($pname, $ppath, $psuffix) = fileparse($phy, qr/\.[^.]*/);
    my $outfile  = $pname."-paml.out";
    
    my $ctl_file = "      seqfile = $seqfile
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

    my $control_file = File::Spec->catfile($ppath, "baseml.ctl");
    open my $out, '>', $control_file or die "\nERROR: Could not open file: $control_file\n";
    print $out $ctl_file;
    close $out;

    return { outfile         => $outfile,
	     control_file    => $control_file };
}

sub parse_baseml {
    my $self = shift;
    my ($args) = @_;
    my $divergence_file = $args->{divergence_file};
    my $outfile         = $args->{outfile};
    my $phylip          = $args->{phylip};
    my $treefile        = $args->{treefile};
    my $control_file    = $args->{control_file};
    my $results_dir     = $args->{results_dir};

    my $divfile = basename($divergence_file);
    my $out = basename($outfile);
    my $wd = getcwd();
    
    open my $divin, '<', $out or die "ERROR: Could not open outfile: $!\n";
    open my $divout, '>', $divfile or die "ERROR: Could not open divergence file: $!\n";
    
    while (<$divin>) {
	chomp;
	if (/^\d\w+\s+\d\.\d+\(\s/) {
	    # 3prime_Ung         0.0269( 8.7752)
	    my ($seqid,$divergence_time,$kappa) = split /\s+/;

	    $divergence_time =~ s/\($//;
	    $kappa =~ s/\)$//;
	    my $time = $divergence_time/(1e-8 * 2);   # T=k/2r, k=1.0 10-8

	    # alignID divergence age Ts:Tv
	    say $divout join "\t", basename($phylip), $divergence_time, $time, $kappa;
	}
	elsif (/^(\d\w+)         (\d\.\d+\()(\d+\.\d+\))/) {
	    # 3prime_RL1         0.0087(999.0000)

	    my $divergence_time = $2;
	    my $kappa = $3;
	    $divergence_time =~ s/\($//;
	    $kappa =~ s/\)$//;
	    my $time = $divergence_time/(1e-8 * 2);

	    # alignID divergence age Ts:Tv
	    say $divout join "\t", basename($phylip), $divergence_time, $time, $kappa;
	}
	elsif (/^(\d\w+)         (\d\.\d+\()(\-\d+\.\d+\))/) {
	    # 3prime_RL1         0.0017(-0.0025)

	    my $divergence_time = $2;
	    my $kappa = $3;
	    $divergence_time =~ s/\($//;
	    $kappa =~ s/\)$//;
	    my $time = $divergence_time/(1e-8 * 2);

	    # alignID divergence age Ts:Tv
	    say $divout join "\t", basename($phylip), $divergence_time, $time, $kappa;
	}
    }
    close $divin;
    close $divout;

    my $resdir = basename($results_dir);
    my $dest_file = File::Spec->catfile($resdir, $divfile);
    copy($divfile, $dest_file) or die "\nERROR: Move failed: $!";
    unlink basename($control_file);

    if ($self->clean) {
	# remove the PAML output but keep the summary produced
	# by this script in a separate directory.
	unlink "2base.t", "rub", "rst", "rst1", "lnf", "rates", "in.basemlg", $divfile, $out;
    }
    
    return;
}

sub _build_baseml_exec { # this should probably be a separate role 
    my $self  = shift;
    my $blexe = $self->get_baseml_exec; # set during class initialization
    
    # check if the executable path was set first
    if (defined $blexe && -e $blexe && -x $blexe) {
        return $blexe;
    }
    elsif (! defined $blexe) {
	my @path = split /:|;/, $ENV{PATH};
	for my $p (@path) {
	    my $bl = File::Spec->catfile($p, 'baseml');

	    if (-e $bl && -x $bl) {
		$self->set_baseml_exec($bl);
		return $bl;
	    }
	}
    }
    else {
        $log->error("Unable to find 'baseml' executable. Make sure PAML is installed. Exiting.");
        exit(1);
    }
}

1;
