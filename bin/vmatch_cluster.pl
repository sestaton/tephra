#!/usr/bin/env perl

use 5.020;
use warnings;
use File::Find;
use File::Basename;
use IPC::System::Simple qw(system);
use Try::Tiny;
use Getopt::Long;
use Cwd;
use Data::Dump;
use experimental 'signatures';

my %opt;
my %ltrs;

GetOptions(\%opt, 'dir|d=s');

$opt{dir} //= getcwd();
my $vmatch_args = collect_feature_args($opt{dir});
#dd $vmatch_args and exit;
cluster_features($vmatch_args);

sub collect_feature_args ($dir) {
    my (@fiveltrs, @threeltrs, @ppt, @pbs, @pdoms, %vmatch_args); 
    find( sub { push @fiveltrs, $File::Find::name if -f and /5prime-ltrs.fasta$/ }, $dir);
    find( sub { push @threeltrs, $File::Find::name if -f and /3prime-ltrs.fasta$/ }, $dir);
    find( sub { push @ppt, $File::Find::name if -f and /ppts.fasta$/ }, $dir);
    find( sub { push @pbs, $File::Find::name if -f and /pbs.fasta$/ }, $dir);
    find( sub { push @pdoms, $File::Find::name if -f and /pdoms.fasta$/ }, $dir);

    # ltr 
    my $fiveltr = $fiveltrs[0];
    my $fiveargs = "-dbcluster 80 20 dbcluster-5primeseqs -s -p -d -seedlength 10 ";
    $fiveargs .= "-exdrop 7 -l 80 -showdesc 0 -sort ld -best 100 -identity 80";
    $vmatch_args{fiveltr} = { seqs => \@fiveltrs, args => $fiveargs };
    
    my $threeltr = $threeltrs[0];
    my $threeargs = "-dbcluster 80 20 dbcluster-3primeseqs -s -p -d -seedlength 10 ";
    $threeargs .= "-exdrop 7 -l 80 -showdesc 0 -sort ld -best 100 -identity 80";
    $vmatch_args{threeltr} = { seqs => \@threeltrs, args => $threeargs };

    # pbs/ppt
    #$vmatch -dbcluster 90 90 -identity 90 -seedlength 3 -exdrop 2
    ### args=-dbcluster 90 90 dbcluster-pbs -s -p -d -seedlength 3 -exdrop 2 -l 3 -showdesc 0 -sort ld -best 100
    ### args=-dbcluster 90 90 dbcluster-ppt -s -p -d -seedlength 3 -exdrop 2 -l 3 -showdesc 0 -sort ld -best 100

    # pdoms
    #$vmatch -dbcluster 80 -seedlength 10 -exdrop 3
    
    # chromo
    #-dbcluster 80 80 dbcluster-chromo -s -p -d -seedlength 10 -l 40 -exdrop 3 -showdesc 0 -sort ld -best 100
    
    # gag
    # args=-dbcluster 80 80 dbcluster-gag -s -p -d -seedlength 10 -l 40 -exdrop 3 -showdesc 0 -sort ld -best 100 
    
    # rnasheh
    ## args=-dbcluster 80 80 dbcluster-rnaseh -s -p -d -seedlength 10 -l 40 -exdrop 3 -showdesc 0 -sort ld -best 100
    
    # rvt 1/2
    # args=-dbcluster 80 80 dbcluster-rvt2 -s -p -d -seedlength 10 -l 40 -exdrop 3 -showdesc 0 -sort ld -best 100
    # args=-dbcluster 80 80 dbcluster-rvt1 -s -p -d -seedlength 10 -l 40 -exdrop 3 -showdesc 0 -sort ld -best 100
    
    # rve
    # args=-dbcluster 80 80 dbcluster-rve -s -p -d -seedlength 10 -l 40 -exdrop 3 -showdesc 0 -sort ld -best 100

    return \%vmatch_args;
}

sub cluster_features ($args) {
    my $vmatch  = '/usr/local/bioinfo/vmatch/vmatch-2.2.4-Linux_x86_64-64bit/vmatch';
    my $mkvtree = '/usr/local/bioinfo/vmatch/vmatch-2.2.4-Linux_x86_64-64bit/mkvtree';

    for my $type (keys %$args) {
	for my $db (@{$args->{$type}{seqs}}) {
	    my ($name, $path, $suffix) = fileparse($db, qr/\.[^.]*/);
	    my $index = File::Spec->catfile($path, $name.".index");
	    my $vmrep = File::Spec->catfile($path, $name."_vmatch-out.txt");
	    my $log   = File::Spec->catfile($path, $name."_vmatch-out.log");;
	    my $mkvtreecmd = "time $mkvtree -db $db -dna -indexname $index -allout -v -pl 2>&1 > $log";
	    my $vmatchcmd  = "time $vmatch $args->{$type}{args} $index > $vmrep";
	    say STDERR "=====> Running mkvtree on $type";
	    run_cmd($mkvtreecmd);
	    say STDERR "=====> Running vmatch on $type";
	    run_cmd($vmatchcmd);
	    unlink glob "$index*";
	}
    }   
}

sub run_cmd ($cmd) {
    try {
	system([0..5], $cmd);
    }
    catch {
	die "\nERROR: $cmd failed. Here is the exception: $_\n";
    };
}
