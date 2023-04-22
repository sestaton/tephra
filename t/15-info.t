#!/usr/bin/env perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use Cwd                 qw(getcwd);
use autodie             qw(open);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use File::Spec;
use Scalar::Util        qw(looks_like_number);

use Test::More tests => 13;

$| = 1;

my $cwd = getcwd();
$ENV{PATH}   = join ":", File::Spec->catdir($cwd, 'blib', 'bin'), $ENV{PATH};
my $cmd      = File::Spec->catfile('blib', 'bin', 'tephra');

{
    my @help_args = ($cmd, 'info', '-h');
    my ($stdout, $stderr, $exit) = capture { system(@help_args) };
    #say STDERR "stderr: $stderr";
    ok($stderr, 'Can execute info subcommand');
}

my ($stdout, $stderr, $exit) = capture { system(@help_args) };

my ($tephra_ver, $perl_ver, $blast_ver, $emboss_ver, $gt_ver, $hmmer2_ver,
    $hmmer3_ver, $htslib_ver, $hscan_ver, $mus_ver, $paml_ver, $vmatch_ver);

for my $line (split /^/, $stdout) {
    if ($line =~ /Tephra v(\d+.\d+.\d+) \(using Perl version v(\d+.\d+.\d+)\)/) {
	($tephra_ver, $perl_ver) = ($1, $2);
    }
    elsif ($line =~ /BLAST+ (blastn)\s+\|  v(\d+.\d+.?\d+?)/) {
	$blast_ver = $1;
    }
    elsif ($line =~ /EMBOSS \(transeq\)\s+\|  v(\d+.\d+.\d+.\d+)/) {
	$emboss_ver = $1;
    }
    elsif ($line =~ /GenomeTools \(LTRharvest\, LTRdigest\, TIRvish\, Tallymer\) \|  v(\d+.\d+.\d+)/) {
	$gt_ver = $1;
    }
    elsif ($line =~ /HMMERv2 \(soloLTR search\)\s+\| v(2.3.2)/) {
	$hmmer2_ver = $1;
    }
    elsif ($line = ~/HMMERv3 \(LTRdigest\)\s+\|  v(3.1b)/) {
	$hmmer3_ver = $1;
    }
    elsif ($line = ~/HTSlib\s+\|  v(\d+.\d+.\d+)/) {
	$htslib_ver = $1;
    }
    elsif ($line = ~/HelitronScanner\s+\|  v(\d+.\d+)/) {
	$hscan_ver = $1;
    }
    elsif ($line = ~/Muscle\s+\|  v(\d+.\d+.\d+)/) {
	$mus_ver = $1;
    }
    elsif ($line = ~/PAML \(baseml\)\s+\|  v\d+.\d+.\d+)/) {
	$paml_ver = $1;
    }
    elsif ($line = ~/Vmatch\s+\|v(\d+.\d+.\d+)/) {
	$vmatch_ver = $1;
    }
}

ok( looks_like_number($gt_ver),     'Can get GenomeTools version from info subcommand' );
ok( looks_like_number($vmatch_ver), 'Can get Vmatch version from info subcommand' );
ok( looks_like_number($mus_ver),    'Can get Muscle version from info subcommand' );
ok( looks_like_number($blast_ver),  'Can get BLAST+ version from info subcommand' );
ok( looks_like_number($emboss_ver), 'Can get EMBOSS version from info subcommand' );
ok( looks_like_number($hmmer2_ver), 'Can get HMMER2 version from info subcommand' );
ok( looks_like_number($hmmer3_ver), 'Can get HMMER3 version from info subcommand' );
ok( looks_like_number($htslib_ver), 'Can get HTSlib version from info subcommand' );
ok( looks_like_number($hscan_ver),  'Can get HelitronScanner version from info subcommand' );
ok( looks_like_number($paml_ver),   'Can get PAML version from info subcommand' );
ok( looks_like_number($perl_ver),   'Can get Perl version from info subcommand' );
ok( looks_like_number($tephra_ver), 'Can get Tephra version from info subcommand' );

done_testing();


