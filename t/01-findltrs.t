#!/usr/bin/env perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use autodie             qw(open);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use List::Util          qw(sum);
use File::Find;
use File::Spec;
use Data::Dump;

use Test::More tests => 8;

my $bindir = File::Spec->catdir('t', 'gt', 'bin');
local $ENV{PATH} = "$bindir:$ENV{PATH}";

my $cmd     = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir = File::Spec->catdir('t', 'test_data');
my $genome  = File::Spec->catfile($testdir, 'ref.fas');
my $model   = File::Spec->catfile($testdir, 'te.hmm');
my $trnas   = File::Spec->catfile($testdir, 'trnas.fas');

my @assemb_results = capture { system([0..5], "$cmd findltrs -h") };

ok(@assemb_results, 'Can execute findltrs subcommand');

my $find_cmd = "$cmd findltrs -g $genome -t $trnas -p $model --clean";
say STDERR $find_cmd;

my ($stdout, $stderr, @ret) = capture { system([0..5], $find_cmd) };

my @files;
find( sub { push @files, $File::Find::name if /\.gff3$/ }, $testdir);
ok( @files == 3, 'Can find some ltrs' ); # 2 ltrdigest files + combined file

my $combined;
for my $line (split /^/, $stderr) {
    if ($line =~ /Number of elements filtered/) {
	my @tot = ($line =~ /=(\d)/g);
	ok( @tot == 5, 'Correct number of filters applied to results' );
	my $sum = sum(@tot);
	ok( $sum == 0, 'Correct number of elements filtered' );
    }
    elsif ($line =~ /Number of elements found/) {
	my @tot = ($line =~ /=(\d)/g);
	ok( @tot == 4, 'Correct number of refinement steps' );
	($combined) = ($line =~ /Combined=(\d)/);
	ok( $combined == 6, 'Correct number of combined elements' );
    }
    elsif ($line =~ /Total elements written/) {
	my ($tot) = ($line =~ /\: (\d)/);
	ok( $tot == 6, 'Correct number of total elements' );
	ok( $tot == $combined, 'Correct number of total and combined elements' );
    }
}

## clean up
my @outfiles;
find( sub { push @outfiles, $File::Find::name if /^ref_ltr/ && ! /filtered.gff3$/ }, $testdir);
unlink @outfiles;
    
done_testing();
