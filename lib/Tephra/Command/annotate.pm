package Chloro::Command::assemble;
# ABSTRACT: Run many chloroplast genome assemblies and pick the best one.

use 5.010;
use strict;
use warnings;
use Chloro -command;
use Cwd                 qw(abs_path);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(:all);
use File::Basename;
use File::Spec;

sub opt_spec {
    return (    
	[ "paired|p=s",   "A file of paired, interleaved chlorplast sequences to be assembled"        ],
	[ "unpaired|u=s", "The file of unpaired, singleton sequences"                                 ], 
	[ "threads|t=i",  "The number of threads (hash steps) to execute simultaneously (Default: 1)" ],
	[ "hashs|s=i",    "The starting hash size (Default: 59)"                                      ],
	[ "hashe|e=i",    "The maximum hash size (Default: 89)"                                       ],
    );
}

sub validate_args {
    my ($self, $opt, $args) = @_;

    my $command = __FILE__;
    if ($self->app->global_options->{man}) {
	system([0..5], "perldoc $command");
    }
    elsif ($self->app->global_options->{help}) {
	$self->help;
    }
    elsif (!$opt->{paired} && !$opt->{unpaired}) {
	say "\nERROR: Required arguments not given.";
	$self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    exit(0) if $self->app->global_options->{man} ||
	$self->app->global_options->{help};

    my $result = _run_assembly($opt);
}

sub _run_assembly {
    my ($opt) = @_;
    my $pairfile   = abs_path($opt->{paired});
    my $singletons = abs_path($opt->{unpaired}); 
    my $hashs      = $opt->{hashs};
    my $hashe      = $opt->{hashe};
    my $thread     = $opt->{threads};

    my $file = __FILE__;
    my $cmd_dir = basename(dirname(abs_path($file)));
    my $hmm_dir = basename(dirname($cmd_dir));
    my $chl_dir = basename(dirname($hmm_dir));
    my $vel_dir = File::Spec->catdir($chl_dir, 'src', 'velvet');
    my $velveth = File::Spec->catfile($chl_dir, 'src', 'velvet', 'velveth');
    my $velvetg = File::Spec->catfile($chl_dir, 'src', 'velvet', 'velvetg');
    my $vo      = File::Spec->catfile(abs_path($chl_dir), 'src', 'VelvetOptimiser', 'VelvetOptimiser.pl');

    local $ENV{PATH} = "$vel_dir:$ENV{PATH}";

    unless (-e $vo) {
        die "\nERROR: 'VelvetOptimiser' not found. Please run the 'install_deps.sh' script before proceeding. Exiting.\n";
    }

    unless (-e $velveth && -e $velvetg) {
        die "\nERROR: Velvet executables not found. Please run the 'install_deps.sh' script before proceeding. Exiting.\n";
    }

    my $exit_value;
    $hashs  //= 59;
    $hashe  //= 89;
    $thread //= 1;
    my ($ifile, $idir, $iext) = fileparse($pairfile, qr/\.[^.]*/);
    my $dirname = "VelvetOpt_k$hashs-k$hashe";

    my @vo_cmd = "$vo ".
	         "-s $hashs ".
		 "-e $hashe ".
		 "-t $thread ".
		 "-p $dirname ".
		 "-d $dirname ".
		 "-f '-fasta -shortPaired $pairfile -fasta -short $singletons'";

    my ($stdout, $stderr, @res) = capture { system([0..5], @vo_cmd); };

    unless (-d $dirname) {
	say "\nERROR: VelvetOptimiser seems to have exited. Here is the message:\n$stderr" if $stderr;
    }
    
    say "\nAssembly results can be found in 'contigs.fa' in the directory: $dirname.\n".
	"See '${dirname}_logfile.txt' for assembly details.\n";

    return $exit_value;
}

sub help {
    print STDERR<<END

USAGE: chloro assemble [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -p|paired     :   A file of paired, interleaved chlorplast sequences to be assembled.
    -u|unpaired   :   The file of unpaired, singleton sequences.

Options:
    -t|threads    :   The number of threads (hash steps) to execute simultaneously (Default: 1).
    -s|hashs      :   The starting hash size (Default: 59).
    -e|hashe      :   The maximum hash size (Default: 89).

END
}


1;
__END__

=pod

=head1 NAME
                                                                       
 chloro assemble - Run many chloroplast genome assemblies and pick the best one

=head1 SYNOPSIS    

 chloro assemble -p reads_paired_interl.fq -u reads_unpaired.fq -s 59 -e 99

=head1 DESCRIPTION

This command will take a file of paired and a file of unpaired reads, for example the files
produced by the 'chloro screenreads' command, and run VelvetOptimiser for a range of hash lengths
and produce the best assembly found.

=head1 AUTHOR 

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -p, --paired

A file of interleaved, paired reads in FASTA format.

=item -u, --unpaired

A file of unpaired reads in FASTA format.

=back

=head1 OPTIONS

=over 2

=item -t, --treads

The number of threads to use with VelvetOptimiser (Default: 1).

=item -s, --hashs

The starting hash length for Velvet (Default: 59).

=item -e, --hashe

The ending hash length for Velvet (Default: 89).

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=back

=cut
