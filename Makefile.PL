use strict;
use warnings FATAL => 'all';
use ExtUtils::MakeMaker;
use Cwd;
use lib qw(lib);
use Tephra::Config::Install;

if (eval ($ExtUtils::MakeMaker::VERSION) <= 6.55) {
    print "\n[ERROR]: ExtUtils::MakeMaker >= 6.5503 is required to build the package correctly.\n";
    print "The installed version is $ExtUtils::MakeMaker::VERSION.\n";
    print "Install the latest ExtUtils::MakeMaker and try again. Exiting.\n";
    exit(1);
}

WriteMakefile(
    NAME             => 'Tephra',
    AUTHOR           => q{S. Evan Staton <statonse@gmail.com>},
    VERSION_FROM     => 'bin/tephra',
    LICENSE          => 'MIT',
    PL_FILES         => {},
    MIN_PERL_VERSION => 5.010,
    CONFIGURE_REQUIRES => {
        'ExtUtils::MakeMaker' => 6.5503,
    },
    BUILD_REQUIRES => {
        'Test::More' => 0,
    },
    PREREQ_PM => {
	'Moose'                      => 0,
	'MooseX::Types::Path::Class' => 0,
	'Parallel::ForkManager'      => 0,
	'IPC::System::Simple'        => 0,
	'Sort::Naturally'            => 0,
	'List::MoreUtils'            => 0,
	'List::UtilsBy'              => 0,
	'Number::Range'              => 0,
	'Set::IntervalTree'          => 0,
	'Statistics::Descriptive'    => 0,
	'App::Cmd'                   => 0,
	'Log::Any'                   => 0,
	'Try::Tiny'                  => 0,
	'Capture::Tiny'              => 0,
	'HTTP::Tiny'                 => 0,
	'HTML::TreeBuilder'          => 0,
	'Net::FTP'                   => 0,
    },
    INST_SCRIPT => 'blib/bin',
    EXE_FILES   => [ 'bin/tephra', 'bin/tephra-getphmm' ],
    dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    test  => { TESTS => 't/*.t' },
    clean => { FILES => [ 't/test_data/*.fai', 'src/*.o', 'src/tephra-*' ] },
);

{
    print STDERR "=====> Configuring Tephra, this may take some time...";
    my $cwd = getcwd();
    my $confobj = Tephra::Config::Install->new( workingdir => $cwd );
    my $config  = $confobj->configure_root;
}

