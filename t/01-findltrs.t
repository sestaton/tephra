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

use Test::More tests => 9;

my $cmd     = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir = File::Spec->catdir('t', 'test_data');
my $genome  = File::Spec->catfile($testdir, 'ref.fas');
my $model   = File::Spec->catfile($testdir, 'te.hmm');
my $trnas   = File::Spec->catfile($testdir, 'trnas.fas');

my $config  = write_config($testdir);
ok( -e $config, 'Can create config file for testing' );

my @results = capture { system([0..5], "$cmd findltrs -h") };
ok(@results, 'Can execute findltrs subcommand');

my $find_cmd = "$cmd findltrs -c $config -g $genome -t $trnas -d $model --clean";
#say STDERR $find_cmd;

my ($stdout, $stderr, @ret) = capture { system([0..5], $find_cmd) };
      
my @files;
find( sub { push @files, $File::Find::name if /\.gff3$/ }, $testdir);

ok( @files == 3, 'Can find some ltrs' ); # 2 ltrdigest files + combined file

my $combined;
for my $line (split /^/, $stderr) {
    if ($line =~ /Number of elements filtered/) {
	my @tot = ($line =~ /=(\d)/g);
	ok( @tot == 3, 'Correct number of filters applied to results' );
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
find( sub { push @outfiles, $File::Find::name if /^ref_ltr/ && ! /combined_filtered.gff3$/ }, $testdir);
unlink @outfiles;
unlink $config;
    
done_testing();

sub write_config {
    my ($testdir) = @_;
    my $config = File::Spec->catfile($testdir, 'tephra_ltr_config.yml');

my $conf = "ltrharvest:
  - mintsd: 4
  - maxtsd: 6
  - minlenltr: 100
  - maxlenltr: 6000
  - mindistltr: 1500
  - maxdistltr: 25000
  - seedlength: 30
  - tsdradius: 60
  - xdrop: 5
  - swmat: 2 
  - swmis: -2
  - swins: -3
  - swdel: -3
  - overlaps: best
ltrdigest:
  - pptradius: 30
  - pptlen: 8 30
  - pptagpr: 0.25
  - uboxlen: 3 30
  - uboxutpr: 0.91
  - pbsradius: 30
  - pbslen: 11 30
  - pbsoffset: 0 5
  - pbstrnaoffset: 0 5
  - pbsmaxeditdist: 1
  - pdomevalue: 10E-6
  - pdomcutoff: NONE
  - maxgaplen: 50";

    open my $out, '>', $config;
    say $out $conf;
    close $out;

    return $config;
}
