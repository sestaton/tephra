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
use Parallel::ForkManager;
use Time::HiRes qw(gettimeofday);
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
    find( sub { push @pdoms, $File::Find::name if -f and /pdom.fasta$/ }, $dir);

    # ltr 
    my $fiveargs = "-dbcluster 80 20 dbcluster-5primeseqs -s -p -d -seedlength 10 ";
    $fiveargs .= "-exdrop 7 -l 80 -showdesc 0 -sort ld -best 100 -identity 80";
    $vmatch_args{fiveltr} = { seqs => \@fiveltrs, args => $fiveargs };
    
    my $threeargs = "-dbcluster 80 20 dbcluster-3primeseqs -s -p -d -seedlength 10 ";
    $threeargs .= "-exdrop 7 -l 80 -showdesc 0 -sort ld -best 100 -identity 80";
    $vmatch_args{threeltr} = { seqs => \@threeltrs, args => $threeargs };

    # pbs/ppt
    my $pbsargs = "-dbcluster 90 90 dbcluster-pbs -s -p -d -seedlength 5 -exdrop 2 ";
    $pbsargs .= "-l 3 -showdesc 0 -sort ld -best 100";
    $vmatch_args{pbs} = { seqs => \@pbs, args => $pbsargs, prefixlen => 5 };

    my $pptargs = "-dbcluster 90 90 dbcluster-ppt -s -p -d -seedlength 5 -exdrop 2 ";
    $pptargs .= "-l 3 -showdesc 0 -sort ld -best 100";
    $vmatch_args{ppt} = { seqs => \@ppt, args => $pptargs, prefixlen => 5 };
    
    # pdoms
    my $pdomargs = "-dbcluster 80 80 dbcluster-pdom -s -p -d -seedlength 10 -exdrop 3 ";
    $pdomargs .= "-l 40 -showdesc 0 -sort ld -best 100";
    $vmatch_args{pdoms} = { seqs => \@pdoms, args => $pdomargs };

    return \%vmatch_args;
}

sub cluster_features ($args) {
    #my $vmatch  = '/usr/local/bioinfo/vmatch/vmatch-2.2.4-Linux_x86_64-64bit/vmatch';
    #my $mkvtree = '/usr/local/bioinfo/vmatch/vmatch-2.2.4-Linux_x86_64-64bit/mkvtree';

    my $t0 = gettimeofday();
    my $doms = 0;
    my %reports;
    my $outfile = 'all_vmatch_reports.txt';
    open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n"; 
    my $pm = Parallel::ForkManager->new(12);
    $pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
			  for my $bl (sort keys %$data_ref) {
			      open my $report, '<', $bl or die "\nERROR: Could not open file: $bl\n";
			      print $out $_ while <$report>;
			      close $report;
			      unlink $bl;
			  }
			  my $t1 = gettimeofday();
			  my $elapsed = $t1 - $t0;
			  my $time = sprintf("%.2f",$elapsed/60);
			  say basename($ident),
			      " just finished with PID $pid and exit code: $exit_code in $time minutes";
			} );
    
    for my $type (keys %$args) {
	for my $db (@{$args->{$type}{seqs}}) {
	    $doms++;
	    $pm->start($db) and next;
	    my $vmrep = process_cluster_args($args, $type, $db);
	    $reports{$vmrep} = 1;
    
	    $pm->finish(0, \%reports);
	}
    }

    $pm->wait_all_children;
    close $out;

    my $t2 = gettimeofday();
    my $total_elapsed = $t2 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);

    say "\n========> Finished running vmatch on $doms domains in $final_time minutes";

}

sub process_cluster_args ($args, $type, $db) {
    my $vmatch  = '/usr/local/bioinfo/vmatch/vmatch-2.2.4-Linux_x86_64-64bit/vmatch';
    my $mkvtree = '/usr/local/bioinfo/vmatch/vmatch-2.2.4-Linux_x86_64-64bit/mkvtree';
    
    my ($name, $path, $suffix) = fileparse($db, qr/\.[^.]*/);
    my $index = File::Spec->catfile($path, $name.".index");
    my $vmrep = File::Spec->catfile($path, $name."_vmatch-out.txt");
    my $log   = File::Spec->catfile($path, $name."_vmatch-out.log");;

    my $mkvtreecmd = "time $mkvtree -db $db -dna -indexname $index -allout -v -pl ";
    if (defined $args->{$type}{prefixlen}) {
	$mkvtreecmd .= "$args->{$type}{prefixlen} ";
    }
    $mkvtreecmd .= "2>&1 > $log";
    my $vmatchcmd  = "time $vmatch $args->{$type}{args} $index > $vmrep";
    #say STDERR "=====> Running mkvtree on $type";
    run_cmd($mkvtreecmd);
    #say STDERR "=====> Running vmatch on $type";
    run_cmd($vmatchcmd);
    unlink glob "$index*";

    return $vmrep;
}

sub run_cmd ($cmd) {
    try {
	system([0..5], $cmd);
    }
    catch {
	die "\nERROR: $cmd failed. Here is the exception: $_\n";
    };
}
