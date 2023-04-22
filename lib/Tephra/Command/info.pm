package Tephra::Command::info;
# ABSTRACT: Show version information for all external programs used by Tephra.

use 5.014;
use strict;
use warnings;
use Pod::Find     qw(pod_where);
use Pod::Usage    qw(pod2usage);
use Cwd           qw(abs_path);
use Capture::Tiny qw(capture_merged capture);
use File::Spec;
use File::Basename;
use Tephra::Config::Exe;
use Tephra -command;
use Data::Dump::Color;

our $VERSION = '0.13.2';

sub opt_spec {
    return (    
	[ "help|h",        "Display the usage menu and exit. "  ],
        [ "man|m",         "Display the full manual. "          ],
    );
}

sub validate_args {
    my ($self, $opt, $args) = @_;

    my $command = __FILE__;
    if ($opt->{man}) {
        system('perldoc', $command) == 0 or die $!;
        exit(0);
    }
    elsif ($opt->{help}) {
        $self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    my $success = _get_command_info($opt);
}

sub _get_command_info {
    my ($opt) = @_;

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($gt, $vmatch, $baseml, $transeq, $blastn, $muscle)
	= @{$config}{qw(gt vmatch baseml transeq blastn muscle)};

    my $emboss_ver = get_emboss_version($transeq);
    my $vmatch_ver = get_vmatch_version($vmatch);
    #my $paml_ver   = get_paml_version($baseml);
    my $gt_ver     = get_gt_version($gt);
    my $blast_ver  = get_blastplus_version($blastn);
    
    my %progs = (
	'HelitronScanner         ' => '                              |  v1.1     |',
	'Vmatch                  ' => "                              |  v$vmatch_ver   |",
	'GenomeTools (LTRharvest, LTRdigest, TIRvish, Tallymer)' => "|  v$gt_ver   |",
	'Muscle                  ' => '                              |  v3.8.31  |',
	'PAML (baseml)           ' => '                              |  v4.10.6  |',
	'BLAST+ (blastn)         ' => "                              |  v$blast_ver   |",
	'HTSlib                  ' => '                              |  v1.3.1   |',
	'EMBOSS (transeq)        ' => "                              |  v$emboss_ver |",
	'HMMERv2 (soloLTR search)' => '                              |  v2.3.2   |',
	'HMMERv3 (LTRdigest)     ' => '                              |  v3.1b    |',
	);

    say STDOUT "\nTephra v$VERSION (using Perl version $^V) has the following the configuration:\n";
    say STDOUT "+------------------------------------------------------+-----------+";
    for my $prog (sort keys %progs) {
	say STDOUT join q{ }, $prog, $progs{$prog};
    }
    say STDOUT "+------------------------------------------------------+-----------+\n";
}

#
# methods
#
sub get_blastplus_version {
    my ($exe) = @_;

    my ($stdout, $stderr, @res) = capture { system($exe, '-version') };
    my ($ver) = ($stdout =~ /blastn: (\d+.\d+.\d+)/g);
    
    return $ver;
}

sub get_emboss_version {
    my ($exe) = @_;

    my ($stdout, $stderr, @res) = capture { system($exe, '-version') };
    my ($ver) = ($stderr =~ /EMBOSS:(\d+.\d+.\d+.\d+)/g);
    
    return $ver;
}
    
sub get_vmatch_version {
    my ($exe) = @_;

    my ($stdout, $stderr, @res) = capture { system($exe, '-version') };
    my ($ver) = ($stdout =~ /\(Vmatch\)\s+(\d+.\d+.\d+)\s+/g);
    
    return $ver;
}

sub get_gt_version {
    my ($exe) = @_;

    my ($stdout, $stderr, @res) = capture { system($exe, '-version') };
    my ($ver) = ($stdout =~ /\(GenomeTools\)\s+(\d+.\d+.\d+)/g);
    
    return $ver;
}

sub get_paml_version {
    my ($exe) = @_;

    my $desc = capture_merged { system($exe) };
    my ($ver) = ($desc =~ /version (\d+.\d+.\d+),/g);

    return $ver;
}

sub help {
    my $desc = capture_merged {
        pod2usage(-verbose => 99, -sections => "NAME|DESCRIPTION", -exitval => "noexit",
		  -input => pod_where({-inc => 1}, __PACKAGE__));
    };
    chomp $desc;
    print STDERR<<END
$desc
USAGE: tephra info [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -c|config     :   The Tephra configuration file.

END
}

1;
__END__

=pod

=head1 NAME
                                                                       
 tephra info - Show version information for all external programs used by Tephra. 

=head1 SYNOPSIS    

 tephra info

=head1 DESCRIPTION

 This command will report the versions of all configured external programs used by Tephra. This subcommand
 is very useful for citation purposes, and for repeatabily of analyses so you can compare between runs
 in an informed manner.

=head1 AUTHOR 

S. Evan Staton, C<< <evan at evanstaton.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=back

=head1 OPTIONS

=over 2

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
