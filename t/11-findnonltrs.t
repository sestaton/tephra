#!/usr/bin/env perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use autodie             qw(open);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use File::Path          qw(remove_tree);
use File::Find;
use File::Spec;
use Data::Dump;

use Test::More tests => 2;

#my $bindir = File::Spec->catdir('t', 'gt', 'bin');
#local $ENV{PATH} = "$bindir:$ENV{PATH}";

my $cmd     = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir = File::Spec->catdir('t', 'test_data');
my $genome  = File::Spec->catfile($testdir, 'nonltrgenome'); #, 'Ha1.fa');
my $outdir  = File::Spec->catdir($testdir, 'nonltrdata');
my $devtest = 0;

my @results = capture { system([0..5], "$cmd findnonltrs -h") };
ok( @results, 'Can execute findnonltrs subcommand' );

SKIP: {
    skip 'skip lengthy tests', 1 unless $devtest;
    my $find_cmd = "$cmd findnonltrs -g $genome -o $outdir";
    #say STDERR $find_cmd;

    my ($stdout, $stderr, @ret) = capture { system([0..5], $find_cmd) };
       
    my @files;
    find( sub { push @files, $File::Find::name if /\.gff3$/ }, $outdir);
    ok( @files == 1, 'Can find some non-LTRs' );

    ## clean up
    unlink @files;
    remove_tree( $outdir, { safe => 1 } );
}
    
done_testing();
