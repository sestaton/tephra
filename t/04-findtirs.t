#!/usr/bin/env perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use autodie             qw(open);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use File::Find;
use File::Spec;

use Test::More tests => 4;

$| = 1;

my $devtests = 0;
if (defined $ENV{TEPHRA_ENV} && $ENV{TEPHRA_ENV} eq 'development') {
    $devtests = 1;
}

my $cmd      = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir  = File::Spec->catdir('t', 'test_data');
my $genome   = File::Spec->catfile($testdir, 'ref.fas');
my $gff      = File::Spec->catfile($testdir, 'ref_tirs.gff3');
my $fas      = File::Spec->catfile($testdir, 'ref_tirs.fasta');
my $genefile = File::Spec->catfile($testdir, 'devtest_gene_seqs.fas.gz');

{
    my @help_args = ($cmd, 'findtirs', '-h');
    my ($stdout, $stderr, $exit) = capture { system(@help_args) };
    #say STDERR "stderr: $stderr";
    ok($stderr, 'Can execute findtirs subcommand');
}

if ($devtests) {
    $genome = File::Spec->catfile($testdir, 'TAIR10_chr1.fas');
    $gff    = File::Spec->catfile($testdir, 'TAIR10_chr1_tirs.gff3');
    $fas    = File::Spec->catfile($testdir, 'TAIR10_chr1_tirs.fasta');
}

my @find_cmd = "$cmd findtirs -g $genome -o $gff -r $genefile --clean";
#say STDERR join q{ }, @find_cmd;
my @ret = capture { system([0..5], @find_cmd) };

my @files;
find( sub { push @files, $File::Find::name if /tirs_?(?:filtered)?.gff3$/ }, $testdir);
ok( @files == 1, 'Can find some TIRs' ); # only 1 after rename

my ($fasct, $gffct) = (0, 0);
open my $fasin, '<', $fas;
while (<$fasin>) { chomp; $fasct++ if /^>/; }
close $fasin;

open my $gffin, '<', $gff;
while (<$gffin>) { chomp; next if /^#/; my @f = split /\t/; $gffct++ if $f[2] eq 'terminal_inverted_repeat_element'; }
close $gffin;

#say STDERR "TIR FASCT: $fasct";
#say STDERR "TIR GFFCT: $gffct";

if ($devtests) {
    ok( $fasct == 181,    'Found the correct number of TIRs in FASTA' );
    ok( $gffct == $fasct, 'Found the correct number of TIRs in FASTA and GFF3' );
}
else {
    ok( $fasct == 1,      'Found the correct number of TIRs in FASTA' );
    ok( $gffct == $fasct, 'Found the correct number of TIRs in FASTA and GFF3' );
}

## clean up
#unlink $fas;
    
done_testing();
