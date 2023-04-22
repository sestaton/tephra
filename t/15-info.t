#!/usr/bin/env perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use Cwd                 qw(getcwd);
use autodie             qw(open);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use File::Spec;
#use Data::Dump::Color;

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

my ($stdout, $stderr, $exit) = capture { system($cmd, 'info') };

my ($tephra_ver, $perl_ver, $blast_ver, $emboss_ver, $gt_ver, $hmmer2_ver,
    $hmmer3_ver, $htslib_ver, $hscan_ver, $mus_ver, $paml_ver, $vmatch_ver);

my %vers;
for my $line (split /^/, $stdout) {
    if ($line =~ /^Tephra/) {
	($tephra_ver, $perl_ver) = ($line =~ /Tephra v(\d+.\d+.\d+)\s+\(using Perl version v(\d+.\d+.\d+)\)/g);
	$vers{tephra} = $tephra_ver;
	$vers{perl} = $perl_ver;
    }
    elsif ($line =~ /^BLAST/) {
	($blast_ver) = ($line =~ /^BLAST\+ \(blastn\)\s+\|\s+v(\d+.\d+.?\d+?)/g);
	$vers{blast} = $blast_ver;
    }
    elsif ($line =~ /^EMBOSS/) {
	($emboss_ver) = ($line =~ /^EMBOSS \(transeq\)\s+\|\s+v(\d+.\d+.\d+.\d+)/g);
	$vers{emboss} = $emboss_ver;
    }
    elsif ($line =~ /^GenomeTools/) {
	($gt_ver) = ($line =~ /^GenomeTools \(LTRharvest\, LTRdigest\, TIRvish\, Tallymer\)\s+\|\s+v(\d+.\d+.\d+)/g);
	$vers{gt} = $gt_ver;
    }
    elsif ($line =~ /^HMMERv2/) {
	($hmmer2_ver) = ($line =~ /^HMMERv2 \(soloLTR search\)\s+\|\s+v(\d+.\d+.\d+)/g);
	$vers{hmmer2} = $hmmer2_ver;
    }
    elsif ($line =~ /^HMMERv3/) {
	($hmmer3_ver) = ($line =~ /^HMMERv3 \(LTRdigest\)\s+\|\s+v(\d+.\d+\w+)/g);
	$vers{hmmer3} = $hmmer3_ver;
    }
    elsif ($line =~ /^HTSlib/) {
	($htslib_ver) = ($line =~ /^HTSlib\s+\|\s+v(\d+.\d+.\d+)/g);
	$vers{htslib} = $htslib_ver;
    }
    elsif ($line =~ /^Helitron/) {
	($hscan_ver) = ($line =~ /^HelitronScanner\s+\|\s+v(\d+.\d+)/g);
	$vers{hscan} = $hscan_ver;
    }
    elsif ($line =~ /^Muscle/) {
	($mus_ver) = ($line =~ /^Muscle\s+\|\s+v(\d+.\d+.\d+)/g);
	$vers{muscle} = $mus_ver;
    }
    elsif ($line =~ /^PAML/) {
	($paml_ver) = ($line =~ /^PAML\s+\(baseml\)\s+\|\s+v(\d+.\d+.\d+)/g);
	$vers{paml} = $paml_ver;
    }
    elsif ($line =~ /^Vmatch/) {
	($vmatch_ver) = ($line =~ /^Vmatch\s+\|\s+v(\d+.\d+.\d+)/g);
	$vers{vmatch} = $vmatch_ver;
    }
}

#dd \%vers and exit;
ok( $vers{gt} =~ /\d+/,     'Can get GenomeTools version from info subcommand' );
ok( $vers{vmatch} =~ /\d+/, 'Can get Vmatch version from info subcommand' );
ok( $vers{muscle} =~ /\d+/, 'Can get Muscle version from info subcommand' );
ok( $vers{blast} =~ /\d+/,  'Can get BLAST+ version from info subcommand' );
ok( $vers{emboss} =~ /\d+/, 'Can get EMBOSS version from info subcommand' );
ok( $vers{hmmer2} =~ /\d+/, 'Can get HMMER2 version from info subcommand' );
ok( $vers{hmmer3} =~ /\d+/, 'Can get HMMER3 version from info subcommand' );
ok( $vers{htslib} =~ /\d+/, 'Can get HTSlib version from info subcommand' );
ok( $vers{hscan} =~ /\d+/,  'Can get HelitronScanner version from info subcommand' );
ok( $vers{paml} =~ /\d+/,   'Can get PAML version from info subcommand' );
ok( $vers{perl} =~ /\d+/,   'Can get Perl version from info subcommand' );
ok( $vers{tephra} =~ /\d+/, 'Can get Tephra version from info subcommand' );

done_testing();


