package Tephra::Role::Run::PAML;

use 5.010;
use Moose::Role;
use MooseX::Types::Path::Class;
use File::Spec;
use File::Find;
use File::Basename;
use Try::Tiny;
use Capture::Tiny       qw(capture);
use IPC::System::Simple qw(system);
use File::Path          qw(remove_tree);
use File::Copy          qw(move copy);
use Log::Any            qw($log);
use Cwd;
use Tephra::Config::Exe;
use namespace::autoclean;

=head1 NAME

Tephra::Role::Run::PAML - Helper role for running PAML

=head1 VERSION

Version 0.03.0

=cut

our $VERSION = '0.03.0';
$VERSION = eval $VERSION;

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

has subs_rate => (
    is       => 'ro',
    isa      => 'Num',
    default  => 1e-8,
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
	#system($baseml);
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
    my $seqfile  = basename($phy);
    my $treefile = basename($dnd);
    
    my $cwd = getcwd();
    my ($pname, $ppath, $psuffix) = fileparse($phy, qr/\.[^.]*/);
    my $outfile = $pname.'-paml.out';
    
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

    my $control_file = File::Spec->catfile($ppath, 'baseml.ctl');
    open my $out, '>', $control_file or die "\nERROR: Could not open file: $control_file\n";
    print $out $ctl_file;
    close $out;

    return ({ outfile     => $outfile,
	     control_file => $control_file });
}

sub parse_baseml {
    my $self = shift;
    my ($args) = @_;
    my $subs_rate       = $self->subs_rate;
    my $divergence_file = $args->{divergence_file};
    my $outfile         = $args->{outfile};
    my $phylip          = $args->{phylip};
    my $treefile        = $args->{treefile};
    my $control_file    = $args->{control_file};
    my $results_dir     = $args->{results_dir};

    my $element = basename($phylip);
    $element =~ s/_ltrs_clustal-out.*//;
    my $out = basename($outfile);
    my $wd  = getcwd();
    my $dirobj = Path::Class::Dir->new($wd);
    my $parent = $dirobj->parent;

    open my $divin, '<', $out or die "ERROR: Could not open outfile: $out\n";
    open my $divout, '>', $divergence_file or die "ERROR: Could not open divergence file: $divergence_file\n";

    while (my $line = <$divin>) {
	chomp $line;
	if ($line =~ /^[35]prime?|unk-prime-[fr]?\s+/) {
	    if ($line =~ /(\d+\.\d+)\(\s?(\-?\d+\.\d+)\)/) {
		# 3prime_Ung         0.0269( 8.7752)
		my $divergence_time = $1;
		my $kappa = $2;
		my $time = $divergence_time/($subs_rate * 2);
		
		# alignID divergence age Ts:Tv
		say $divout join "\t", $element, $divergence_time, $time, $kappa;
	    }
	}
    }
    close $divin;
    close $divout;

    my $resdir = basename($results_dir);
    my $dest_file = File::Spec->catfile($parent, $resdir, $divergence_file);
    #say STDERR join q{ }, $divergence_file, $dest_file;
    copy($divergence_file, $dest_file) or die "\nERROR: Move failed: $!";
    chdir $parent or die $!;
    remove_tree($wd, { safe => 1 });
    
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
	my $config = Tephra::Config::Exe->new->get_config_paths;
	my ($pamlbin) = @{$config}{qw(pamlbin)};
	$blexe = File::Spec->catfile($pamlbin, 'baseml');
	if (-e $blexe && -x $blexe) {
	    $self->set_baseml_exec($blexe);
	    return $blexe;
	}

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

=head1 AUTHOR

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests through the project site at 
L<https://github.com/sestaton/tephra/issues>. I will be notified,
and there will be a record of the issue. Alternatively, I can also be 
reached at the email address listed above to resolve any questions.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Tephra::Role::Run::PAML


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1;
