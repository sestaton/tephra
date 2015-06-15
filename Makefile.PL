use strict;
use warnings FATAL => 'all';
use ExtUtils::MakeMaker;

if (eval ($ExtUtils::MakeMaker::VERSION) <= 6.55) {
    print "\n[ERROR]: ExtUtils::MakeMaker >= 6.5503 is required to build the package correctly.\n";
    print "The installed version is $ExtUtils::MakeMaker::VERSION.\n";
    print "Install the latest ExtUtils::MakeMaker and try again. Exiting.\n";
    exit(1);
}

BEGIN {
    system("chmod +x bin/tephra") == 0 or die $!;
    if (! -d 'src') {
	system("./install_deps.sh") == 0 or die $!;
    }
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
	'Moo'                   => 0,       
        'IPC::System::Simple'   => 0,
	'Parallel::ForkManager' => 0,
	'Try::Tiny'             => 0,
    },
    INST_SCRIPT => 'bin',
    dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    test  => { TESTS => 't/*.t' },
    clean => { FILES => ['bin/tephra', 'src'] },
);

