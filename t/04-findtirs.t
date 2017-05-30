#!/usr/bin/env perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use autodie             qw(open);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use File::Find;
use File::Spec;

use Test::More tests => 2;

$| = 1;

my $cmd     = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir = File::Spec->catdir('t', 'test_data');
my $genome  = File::Spec->catfile($testdir, 'ref.fas');
my $gff     = File::Spec->catfile($testdir, 'ref_tirs.gff3');
my $fas     = File::Spec->catfile($testdir, 'ref_tirs.fasta');
## these are subsets for testing
#my $model   = File::Spec->catfile($testdir, 'te.hmm');

{
    my @help_args = ($cmd, 'findtirs', '-h');
    my ($stdout, $stderr, $exit) = capture { system(@help_args) };
        #say STDERR "stderr: $stderr";
    ok($stderr, 'Can execute findtirs subcommand');
}

my @find_cmd = "$cmd findtirs -g $genome -o $gff --clean";
#say STDERR join q{ }, @find_cmd;
my @ret = capture { system([0..5], @find_cmd) };

my @files;
find( sub { push @files, $File::Find::name if /tirs_?(?:filtered)?.gff3$/ }, $testdir);
ok( @files == 1, 'Can find some tirs' ); # only 1 after rename

## clean up
#my @outfiles;
#find( sub { push @outfiles, $File::Find::name if /^ref_tirs/ && ! /$gff/ }, $testdir);
#unlink @outfiles;
unlink $fas;
    
done_testing();
