package Tephra::Config::Install;

use 5.014;
use Moose;
use MooseX::Types::Path::Class;
use Path::Class::File;
use File::Spec;
use File::Find;
use File::Basename;
use Cwd           qw(getcwd abs_path);
use Log::Any      qw($log);
use File::Copy    qw(copy move);
use File::Path    qw(make_path remove_tree);
use Capture::Tiny qw(capture);
use HTTP::Tiny;
use HTML::TreeBuilder;
#use Net::FTP;
use Tephra::Config::Exe;
use namespace::autoclean;
#use Data::Dump::Color;

=head1 NAME

Tephra::Config::Install - Class for setting up Tephra dependencies

=head1 VERSION

Version 0.13.1

=cut

our $VERSION = '0.13.1';

has basedir => (
    is       => 'ro',
    isa      => 'Path::Class::Dir',
    required => 0,
    coerce   => 1,
    default  => sub {
	return $ENV{TEPHRA_DIR} // Path::Class::Dir->new($ENV{HOME}, '.tephra')
    },
);

has workingdir => ( 
    is       => 'ro', 
    isa      => 'Path::Class::Dir', 
    required => 0, 
    coerce   => 1 
);

has debug => (
    is         => 'ro',
    isa        => 'Bool',
    predicate  => 'has_debug',
    lazy       => 1,
    default    => 0,
);

sub configure_root {
    my $self = shift;
    my $basedir = $self->basedir; #->absolute->resolve;
    my $debug = $self->debug;

    my $config = Tephra::Config::Exe->new( basedir => $basedir )->get_config_paths;

    unless (-e $config->{transeq}) {
	say STDERR "getting EMBOSS" if $debug;
	$config->{transeq} = $self->fetch_emboss($config->{tephrabin});
	print STDERR ".";
    }

    unless (-e $config->{gt} && -x $config->{gt} &&
	    -e $config->{gtdata} && -d $config->{gtdata}) {
	say STDERR "getting GenomeTools" if $debug;
	$config->{gt} = $self->fetch_gt_exes($config->{tephrabin});
	print STDERR ".";
    }

    unless (-e $config->{vmatch} && -x $config->{vmatch} &&
	    -e $config->{mkvtree} && -x $config->{mkvtree} &&
	    -e $config->{cleanpp} && -x $config->{cleanpp}) {
	say STDERR "getting Vmatch" if $debug;
	($config->{vmatch}, $config->{mkvtree}, $config->{cleanpp}) = $self->fetch_vmatch_exes($config->{tephrabin});
	print STDERR ".";
    }
    
    unless (-e $config->{hscanjar}) {
	say STDERR "getting HSCAN" if $debug;
	$config->{hscanjar} = $self->fetch_hscan;
	print STDERR ".";
    }
    
    unless (-e $config->{hmmer2bin}) {
	say STDERR "getting HMMER2" if $debug;
	$config->{hmmer2bin} = $self->fetch_hmmer2;
	print STDERR ".";
    }
    
    unless (-e $config->{hmmer3bin}) {
	say STDERR "getting HMMER3" if $debug;
        $config->{hmmer3bin} = $self->fetch_hmmer3;
        print STDERR ".";
    }
    
    unless (-e $config->{trnadb}) {
	say STDERR "getting tRNAdb" if $debug;
	$config->{trnadb} = $self->fetch_trnadb;
	print STDERR ".";
    }

    unless (-e $config->{hmmdb}) {
	say STDERR "getting HMMDB" if $debug;
	$config->{hmmdb} = $self->fetch_hmmdb;
	print STDERR ".";
    }

    unless (-e $config->{modeldir}) {
	say STDERR "getting ModelDB" if $debug;
	$config->{modeldir} = $self->fetch_hmm_models;
	print STDERR ".";
    }
    
    unless (-e $config->{chrhmm}) {
	say STDERR "writing Chromosome HMM file" if $debug;
	$config->{chrhmm} = $self->make_chrom_dir($config->{chrhmm});
	print STDERR ".";
    }
    
    unless (-e $config->{mgescan} && -e $config->{transcmd}) {
	say STDERR "getting MGEScan and translate-cmd" if $debug;
	($config->{mgescan}, $config->{transcmd}) = $self->build_mgescan($config->{tephrabin});
	print STDERR ".";
    }
    
    unless (-e $config->{baseml} && -x $config->{baseml}) {
	say STDERR "getting PAML" if $debug;
	$config->{baseml} = $self->fetch_paml($config->{tephrabin});
	print STDERR ".";
    }
    
    unless (-e $config->{blastn} && -x $config->{blastn} &&
	    -e $config->{makeblastdb} && -x $config->{makeblastdb}) {
	say STDERR "getting BLAST" if $debug;
        ($config->{blastn}, $config->{makeblastdb}) = $self->fetch_blast($config->{tephrabin});
	print STDERR ".";
    }

    unless (-e $config->{htslibdir}) {
	say STDERR "getting HTSlib" if $debug;
        $config->{htslibdir} = $self->fetch_htslib;
	print STDERR ".";
    }

    unless (-e $config->{muscle}) {
        say STDERR "getting MUSCLE" if $debug; 
        $config->{muscle} = $self->fetch_muscle($config->{tephrabin});
        print STDERR ".";
    }

    print STDERR "Done.\n";

    return $config;
}

sub fetch_gt_exes {
    my $self   = shift;
    my ($bindir) = @_;

    my $root = $self->basedir->absolute->resolve;
    my $wd   = $self->workingdir->absolute->resolve;
    
    ##TODO: fetch from github  
    my $host = 'http://genometools.org';
    my $dir  = 'pub/binary_distributions';
    my $file = 'gt_distlisting.html';
    $self->fetch_file($file, $host."/".$dir);
    
    my $tree = HTML::TreeBuilder->new;
    $tree->parse_file($file);
    
    my ($dist, $ldist, $ldir);
    for my $tag ($tree->look_down(_tag => 'a')) {
	if ($tag->attr('href')) {
	    if ($tag->as_text =~ /Linux_x86_64-64bit-barebone.tar.gz\z/) {
		$dist = $tag->as_text;
		my $archive = join "/", $host, $dir, $dist;
		$self->fetch_file($dist, $archive);
		
		$ldist = $dist;
		$ldist =~ s/\.tar.gz\z//;
		$ldir = File::Spec->catdir($root, 'gt');
		
		system("tar xzf $dist") == 0 or die $!;
		
		move $ldist, $ldir or die "\n[ERROR]: move failed: $ldist -> $ldir: $!\n";
		unlink $dist;
	    }
	}
    }
    unlink $file;
    
    my $lgt = File::Spec->catfile($ldir, 'bin', 'gt');
    my $tgt = File::Spec->catfile($bindir, 'gt');
    copy $lgt, $tgt or die "\n[ERROR]: copy failed: $lgt -> $tgt: $!\n";
    chmod 0755, $tgt;

    my $gtdata  = File::Spec->catdir($ldir, 'gtdata');
    my $tgtdata = File::Spec->catdir($root, 'gtdata');
    move $gtdata, $tgtdata or die "\n[ERROR]: move failed: $gtdata -> $tgtdata: $!\n";

    remove_tree( $ldir, { safe => 1 } );

    return $tgt;
}

sub fetch_vmatch_exes {
    my $self   = shift;
    my ($bindir) = @_;

    my $root = $self->basedir->absolute->resolve;
    my $wd   = $self->workingdir->absolute->resolve;
    
    ##TODO: fetch from github  
    my $host = 'http://vmatch.de';
    my $dir  = 'distributions';
    my $file = 'vmatch-2.3.0-Linux_x86_64-64bit.tar.gz';
    my $url     = join "/", $host, $dir, $file;
    my $outfile = File::Spec->catfile($root, $file);

    system("wget -q -O $outfile $url 2>&1 > /dev/null") == 0
        or die $!;
    chdir $root;
    my $dist = 'vmatch-2.3.0-Linux_x86_64-64bit';
    my $ldir = File::Spec->catdir($root, $dist);
    system("tar xzf $file") == 0 or die "tar failed: $!";
    
    my $distfile = File::Spec->catfile($root, $file);
    unlink $distfile;
    
    my $vmatch   = File::Spec->catfile($ldir, 'vmatch');
    my $mkvtree  = File::Spec->catfile($ldir, 'mkvtree');
    my $cleanpp  = File::Spec->catfile($ldir, 'cleanpp.sh');
    my $tvmatch  = File::Spec->catfile($bindir, 'vmatch');
    my $tmkvtree = File::Spec->catfile($bindir, 'mkvtree');
    my $tcleanpp = File::Spec->catfile($bindir, 'cleanpp.sh');
    copy $vmatch, $tvmatch or die "\n[ERROR]: copy failed: $!\n";
    copy $mkvtree, $tmkvtree or die "\n[ERROR]: copy failed: $!\n";
    copy $cleanpp, $tcleanpp or die "\n[ERROR]: copy failed: $!\n";
    chmod 0755, $tvmatch, $tmkvtree, $tcleanpp;

    remove_tree( $ldir, { safe => 1 } );

    return ($tvmatch, $tmkvtree, $tcleanpp);
}

sub fetch_hscan {
    my $self = shift;
    my $root = $self->basedir->absolute->resolve;
    my $wd   = $self->workingdir->absolute->resolve;
    
    my $host = 'https://sourceforge.net';
    my $dir  = 'projects/helitronscanner/files/HelitronScanner_V1.0.zip';
    my $ldir = File::Spec->catdir($root, 'helitronscanner');
    make_path( $ldir, {verbose => 0, mode => 0771,} );
    my $file = 'HelitronScanner_V1.0.zip';
    my $path = File::Spec->catfile($ldir, $file);
    #$self->fetch_file($path, $host."/".$dir);
    my $remote = join "/", $host, $dir;
    chdir $ldir or die $!;
    system("wget -q $remote") == 0 or die $!;
    system("unzip $file 2>&1 > /dev/null") == 0 or die $!;
    
    my $cwd   = getcwd();
    my $hscan = File::Spec->catfile($cwd, 'HelitronScanner', 'HelitronScanner.jar');
    chdir $wd;
    
    return $hscan;
}

sub fetch_blast {
    my $self   = shift;
    my ($bindir) = @_;

    my $root = $self->basedir->absolute->resolve;

    chdir $root or die $!;
    my $host = 'https://ftp.ncbi.nlm.nih.gov';
    my $dir  = 'blast/executables/blast+/2.9.0';
    my $file = 'ncbi-blast-2.9.0+-x64-linux.tar.gz';
    my $url  = join "/", $host, $dir, $file; 

    system("wget -q $url") == 0 or die $!;
    system("tar xzf $file 2>&1 > /dev/null") == 0 or die $!;
    my $bdir = 'ncbi-blast+';
    my $ldir = $file;
    $ldir =~ s/\-x64-linux.tar.gz//;

    unlink $file if -e $file;
    my $blastn   = File::Spec->catfile($ldir, 'bin', 'blastn');
    my $mblastdb = File::Spec->catfile($ldir, 'bin', 'makeblastdb');
    my $tblastn   = File::Spec->catfile($bindir, 'blastn');
    my $tmblastdb = File::Spec->catfile($bindir, 'makeblastdb');

    copy $blastn, $tblastn or die "\n[ERROR]: copy failed: l328 $!\n"; # $tblastn means tephra copy, not to be confused with tblastn
    copy $mblastdb, $tmblastdb or die "\n[ERROR]: copy failed: l329 $!\n";
    chmod 0755, $tblastn, $tmblastdb;
    remove_tree( $ldir, { safe => 1 } );

    return ($tblastn, $tmblastdb);
}

sub fetch_hmmer2 {
    my $self = shift;
    my $root = $self->basedir->absolute->resolve;
    my $wd   = $self->workingdir->absolute->resolve;

    ##TODO: fetch from github
    my $urlbase = 'http://eddylab.org'; 
    my $dir     = 'software';
    my $tool    = 'hmmer';
    my $version = '2.3.2';
    my $file    = 'hmmer-2.3.2.tar.gz';
    my $url     = join "/", $urlbase, $dir, $tool, $version, $file;
    my $outfile = File::Spec->catfile($root, $file);
    $self->fetch_file($outfile, $url);

    chdir $root;
    my $dist = 'hmmer-2.3.2';
    system("tar xzf $file") == 0 or die "tar failed: $!";
    chdir $dist;
    my $cwd = getcwd();
    system("./configure --enable-threads --prefix=$cwd 2>&1 > /dev/null") == 0
	or die "configure failed: $!";
    system("make -j4 2>&1 >/dev/null") == 0 
	or die "make failed: $!";
    system("make install 2>&1 >/dev/null") == 0
	 or die "make failed: $!";
    my $hmmer2bin = File::Spec->catdir($cwd, 'bin');
    my $distfile  = File::Spec->catfile($root, $file);
    unlink $distfile;
    chdir $wd;
    
    return $hmmer2bin;
}

sub fetch_hmmer3 {
    my $self = shift;
    my $root = $self->basedir->absolute->resolve;
    my $wd   = $self->workingdir->absolute->resolve;
    
    ##TODO: fetch from github  
    my $urlbase = 'http://eddylab.org'; 
    my $dir     = 'software';
    my $tool    = 'hmmer3';
    my $version = '3.1b2';
    my $file    = 'hmmer-3.1b2-linux-intel-x86_64.tar.gz';
    my $url     = join "/", $urlbase, $dir, $tool, $version, $file;
    my $outfile = File::Spec->catfile($root, $file);
    $self->fetch_file($outfile, $url);

    chdir $root;
    my $dist  = 'hmmer-3.1b2-linux-intel-x86_64';
    my $ldist = 'hmmer-3.1b2';

    system("tar xzf $file") == 0 or die "tar failed: $!";
    move $dist, $ldist or die "\n[ERROR]: move failed: $!\n";
    chdir $ldist;
    my $cwd = getcwd();
    my $hmmer3bin = File::Spec->catdir($cwd, 'binaries');
    my $distfile  = File::Spec->catfile($root, $file);
    unlink $distfile;
    chdir $wd;
    
    return $hmmer3bin;
}

sub fetch_paml {
    my $self   = shift;
    my ($bindir) = @_;

    my $root = $self->basedir->absolute->resolve;
    my $wd   = $self->workingdir->absolute->resolve;

    my $urlbase = 'http://abacus.gene.ucl.ac.uk';
    my $dir     = 'software';
    my $file    = 'pamlX1.3.1+paml4.8a-win32.tgz';
    my $url     = join "/", $urlbase, $dir, $file;
    my $outfile = File::Spec->catfile($root, $file);
    $self->fetch_file($outfile, $url);

    chdir $root;
    my $dist  = 'paml4.8';
    my $xdist = 'pamlX';
    system("tar xzf $file") == 0 or die "tar failed: $!";
    remove_tree( $xdist, { safe => 1 } );
    unlink $file;

    chdir $dist;
    my $cwd = getcwd();
    my $bin = File::Spec->catdir($cwd, 'bin');
    my @exes;
    find( sub { push @exes, $File::Find::name if -f and /\.exe$/ }, $bin );
    unlink @exes;
    chdir 'src';

    my @results = capture { system('make', '-j4') };
    for my $l (split /^/, @results) {
        if ($l =~ /error/i) {
            say STDERR "\n[ERROR]: 'make' failed for PAML. Please report the error below. Exiting.\n";
            say STDERR @results;
        }
    }

    my @exelist = ('yn00', 'baseml', 'basemlg', 'mcmctree', 'pamp', 'evolver', 'infinitesites', 'codeml');

    for my $file (grep { /baseml\z/ } @exelist) {
	my $binfile = File::Spec->catfile($bindir, $file);
	copy $file, $binfile or die "Copy failed: $file -> l446 $!";
	chmod 0755, $binfile;
    }

    chdir $root;
    remove_tree( $dist, { safe => 1 } );
    my $baseml = File::Spec->catfile($bindir, 'baseml');

    return $baseml;
}

sub fetch_emboss {
    my $self   = shift;
    my ($bindir) = @_;

    my $root = $self->basedir->absolute->resolve;
    my $wd   = $self->workingdir->absolute->resolve;

    my $urlbase = 'ftp://emboss.open-bio.org';
    my $dir     = 'pub';
    my $tool    = 'EMBOSS';
    my $release = 'old';
    my $version = '6.5.0';
    my $file    = 'EMBOSS-6.5.7.tar.gz';
    my $url     = join "/", $urlbase, $dir, $tool, $release, $version, $file;
    #my $file = 'emboss-latest.tar.gz';
    #my $url = join "/", $urlbase, $dir, $tool, $file;
    my $outfile = File::Spec->catfile($root, $file);

    system("wget -q -O $outfile $url 2>&1 > /dev/null") == 0
	or die $!;
    chdir $root;
    my $dist = File::Spec->catdir($root, 'EMBOSS-6.5.7');
    system("tar xzf $file") == 0 or die "tar failed: $!";
    chdir $dist;
    my $cwd = getcwd();
    system("./configure --without-x --without-mysql --disable-shared --prefix=$root 2>&1 > /dev/null") == 0
	or die "configure failed: $!";
    system("make -j4 2>&1 > /dev/null") == 0 
	or die "make failed: $!";
    system("make install 2>&1 > /dev/null") == 0
	or die "make failed: $!";

    ## clean up
    chdir $root or die $!;

    my $distfile = File::Spec->catfile($root, $file);
    unlink $distfile;
    ##################################
    my $share = File::Spec->catdir($root, 'share');
    my $inc   = File::Spec->catdir($root, 'include');
    my $lib   = File::Spec->catdir($root, 'lib');
    my $bin   = File::Spec->catdir($root, 'bin');
    my $data  = File::Spec->catdir($root, 'data');

    remove_tree( $dist, { safe => 1 } );
    remove_tree( $lib,  { safe => 1 } );
    remove_tree( $inc,  { safe => 1 } );

    my (@sharedirs, @sharefiles, @binfiles, @datafiles);
    find(sub { push @sharedirs, $File::Find::name if -d and 
		   /doc|index|jemboss|test|AAINDEX|CODONS|JASPAR|OBO|PRINTS|PROSITE|REBASE|TAXONOMY/ }, $share);
    for my $dir (@sharedirs) { 
	remove_tree( $dir, { safe => 1 } );
    }

    #find(sub { push @datafiles, $File::Find::name if -f }, $share);
    #for my $file (@sharefiles) {
	#my ($name, $path, $suffix) = fileparse($file, qr/\.[^.]*/);
	#unlink $file unless $file =~ /EGC\.0$/;
    #}

    find(sub { push @sharefiles, $File::Find::name if -f }, $share);
    #dd \@sharefiles;
    for my $file (@sharefiles) {
	my ($name, $path, $suffix) = fileparse($file, qr/\.[^.]*/);
	next if $file =~ /EGC.0\z/;
	unlink $file unless $name =~ /^transeq|knowntypes|codes/;
    }

    find(sub { push @binfiles, $File::Find::name }, $bin);
    for my $file (@binfiles) {
	my ($name, $path, $suffix) = fileparse($file, qr/\.[^.]*/);
	unlink $file unless $name =~ /^transeq|gt|tephra|baseml|blast|vmatch|mkvtree|muscle|cleanpp/;
    }
    ##################################
    my $transeq  = File::Spec->catdir($root, 'bin', 'transeq');

    return $transeq;
}

sub fetch_htslib {
    my $self = shift;
    my $root = $self->basedir->absolute->resolve;
    my $wd   = $self->workingdir->absolute->resolve;

    my $urlbase = 'https://github.com';
    my $dir     = 'samtools';
    my $tool    = 'htslib';
    my $release = 'releases/download';
    my $version = '1.3.1';
    my $file    = 'htslib-1.3.1.tar.bz2';
    my $url     = join "/", $urlbase, $dir, $tool, $release, $version, $file;
    my $outfile = File::Spec->catfile($root, $file);

    system("wget -q -O $outfile $url 2>&1 > /dev/null") == 0
	or die $!;
    chdir $root;
    my $dist = 'htslib-1.3.1';
    my $libdir = File::Spec->catdir($root, $dist); #, 'htslib');
    system("tar xjf $file") == 0 or die "tar failed: $!";
    chdir $dist;
    my $cwd = getcwd();
    system("./configure --prefix=$cwd 2>&1 > /dev/null") == 0
	or die "configure failed: $!";
    system("make -j4 2>&1 > /dev/null") == 0 
	or die "make failed: $!";
    system("make install 2>&1 > /dev/null") == 0
	or die "make failed: $!";
    
    my $distfile = File::Spec->catfile($root, $file);
    unlink $distfile;
    chdir $wd;

    $ENV{HTSLIB_DIR} = $libdir;
    #system("cpanm -q Bio::DB::HTS") == 0
	#or die "Installing Bio::DB::HTS failed. Here is the HTSLIB_DIR: $libdir. [ERROR]: $!\n";
    #system('cpanm', '-q', '-n', 'Bio::Root::Version') == 0
        #or die "BioPerl install failed: $!";
    #my @results = capture { system('cpanm', '-q', '-n', 'Bio::Root::Version') };
    my @results = capture { system('cpanm', '-q', 'Bio::DB::HTS') };

    return $libdir;
}

sub fetch_muscle {
    my $self   = shift;
    my ($bindir) = @_;

    my $root = $self->basedir->absolute->resolve;
    my $wd   = $self->workingdir->absolute->resolve;
    
    my $host = 'http://www.drive5.com';
    my $dir  = 'muscle';
    my $page = 'downloads.htm';
    my $file = 'muscle_distlisting.html';
    my $endpoint = join "/", $host, $dir, $page;
    $self->fetch_file($file, $endpoint);
    
    my $tree = HTML::TreeBuilder->new;
    $tree->parse_file($file);
    
    my ($dist, $ldist, $ldir, $musbin);
    for my $tag ($tree->look_down(_tag => 'a')) {
        if ($tag->attr('href')) {
            if ($tag->as_text =~ /muscle(\d+\.\d+\.\d+)_i86linux64.tar.gz\z/) {
		my $ver = $1;
		my $vdir = "downloads$ver";
                $dist = $tag->as_text;
                my $archive = join "/", $host, $dir, $vdir, $dist;
                $self->fetch_file($dist, $archive);
                
                $ldist = $dist;
                $ldist =~ s/\.tar.gz\z//;
                $ldir = File::Spec->catdir($root, 'muscle'); #$ldist);
		make_path( $ldir, {verbose => 0, mode => 0771,} );
                $musbin = File::Spec->catfile($ldir, 'muscle');
                system("tar xzf $dist") == 0 or die $!;
                
                move $ldist, $musbin or die "\n[ERROR]: move failed: l590 $!\n";
                unlink $dist;
            }
        }
    }
    unlink $file;
    
    my $muscle = File::Spec->catfile($bindir, 'muscle');
    copy $musbin, $muscle or die "Copy failed: l598 $!";
    chmod 0755, $muscle;

    remove_tree( $ldir, { safe => 1 } );

    return $muscle;
}

sub fetch_trnadb {
    my $self = shift;
    my $root = $self->basedir->absolute->resolve;

    # make db directory
    my $db_dir = File::Spec->catdir($root, 'TephraDB');
    unless ( -e $db_dir ) {
	make_path( $db_dir, {verbose => 0, mode => 0771,} );
    }

    my $file = 'eukaryotic-tRNAs.fa.gz';
    my $trdb = File::Spec->catdir('build', $file);
    copy $trdb, $db_dir or die "Copy failed: $!";

    ## website down as of 10/5/16
    #my $urlbase = 'http://lowelab.ucsc.edu';
    #my $dir     = 'download';
    #my $release = 'tRNAs';
    #my $file    = 'eukaryotic-tRNAs.fa.gz';
    #my $url     = join "/", $urlbase, $dir, $release, $file;
    #my $outfile = File::Spec->catfile($db_dir, $file);

    #system("wget -q -O $outfile $url 2>&1 > /dev/null") == 0
	#or die $!;
    chdir $db_dir;

    system("gunzip $file") == 0 or die "tar failed: $!";
    $file =~ s/\.gz//;
    my $trnadb = File::Spec->catfile($db_dir, $file); 

    return $trnadb;
}

sub fetch_hmmdb {
    my $self = shift;
    my $root = $self->basedir->absolute->resolve;
    my $wd   = $self->workingdir->absolute->resolve;
   
    chdir $wd;
    # make db directory
    my $db_dir = File::Spec->catdir($root, 'TephraDB');
    unless ( -e $db_dir ) {
	make_path( $db_dir, {verbose => 0, mode => 0771,} );
    }

    ## The HMM library was generated with HMMER2GO (https://github.com/sestaton/HMMER2GO) with the following command:
    ##     hmmer2go pfamsearch -s 'transposable element' -d -o transposable+element.hmm
    ##
    ## This file is distributed with Tephra as of v0.03.8
    my $hmms = 'transposable+element.hmm.gz';
    my $hlib = File::Spec->catfile('build', $hmms); # for LTR-RT search with LTRdigest
    copy $hlib, $db_dir or die "Copy failed: $!";
    chdir $db_dir;
    system("gunzip $hmms") == 0 or die $!;

    $hmms =~ s/\.gz//;
    my $hmmdb = File::Spec->catfile($db_dir, $hmms);

    return $hmmdb;
}

sub fetch_hmm_models {
    my $self = shift;
    my $root = $self->basedir->absolute->resolve;
    my $wd   = $self->workingdir->absolute->resolve;
   
    chdir $wd;
    my $file = 'pHMM.tar.gz';
    my $dist = File::Spec->catfile('build', $file); # for non-LTR-RT search

    copy $dist, $root or die "Copy failed: $!";
    chdir $root;
    system("tar xzf $file") == 0 or die $!;
    unlink $file;

    my $dir = File::Spec->catfile($root, 'pHMM');

    return $dir;
}

sub make_chrom_dir {
    my $self = shift;
    my $root = $self->basedir->absolute->resolve;
    my ($chr_file) = @_;

    open my $out, '>', $chr_file or die "\n[ERROR]: Could not open file: $chr_file\n";

    say $out "Symbol= 4";
    say $out "State= 33";
    say $out "Transition= 73";
    say $out join "\t", '0', '1', '0.0455';
    say $out join "\t", '0', '3', '0.0455';
    say $out join "\t", '0', '4', '0.0455';
    say $out join "\t", '0', '6', '0.0455';
    say $out join "\t", '0', '7', '0.0455';
    say $out join "\t", '0', '9', '0.0455';
    say $out join "\t", '0', '10', '0.0455';
    say $out join "\t", '0', '12', '0.0455';
    say $out join "\t", '0', '13', '0.0455';
    say $out join "\t", '0', '15', '0.0455';
    say $out join "\t", '0', '16', '0.0455';
    say $out join "\t", '0', '18', '0.0455';
    say $out join "\t", '0', '19', '0.0455';
    say $out join "\t", '0', '21', '0.0455';
    say $out join "\t", '0', '22', '0.0455';
    say $out join "\t", '0', '24', '0.0455';
    say $out join "\t", '0', '25', '0.0455';
    say $out join "\t", '0', '27', '0.0455';
    say $out join "\t", '0', '28', '0.0455';
    say $out join "\t", '0', '30', '0.0455';
    say $out join "\t", '0', '31', '0.0455';
    say $out join "\t", '0', '32', '0.0455';
    say $out join "\t", '1', '0', '0.5';
    say $out join "\t", '1', '2', '0.5';
    say $out join "\t", '2', '3', '1.0';
    say $out join "\t", '3', '0', '1.0';
    say $out join "\t", '4', '0', '0.5';
    say $out join "\t", '4', '5', '0.5';
    say $out join "\t", '5', '6', '1.0';
    say $out join "\t", '6', '0', '1.0';
    say $out join "\t", '7', '0', '0.5';
    say $out join "\t", '7', '8', '0.5';
    say $out join "\t", '8', '9', '1.0';
    say $out join "\t", '9', '0', '1.0';
    say $out join "\t", '10', '0', '0.5';
    say $out join "\t", '10', '11', '0.5';
    say $out join "\t", '11', '12', '1.0';
    say $out join "\t", '12', '0', '1.0';
    say $out join "\t", '13', '0', '0.5';
    say $out join "\t", '13', '14', '0.5';
    say $out join "\t", '14', '15', '1.0';
    say $out join "\t", '15', '0', '1.0';
    say $out join "\t", '16', '0', '0.5';
    say $out join "\t", '16', '17', '0.5';
    say $out join "\t", '17', '18', '1.0';
    say $out join "\t", '18', '0', '1.0';
    say $out join "\t", '19', '0', '0.5';
    say $out join "\t", '19', '20', '0.5';
    say $out join "\t", '20', '21', '1.0';
    say $out join "\t", '21', '0', '1.0';
    say $out join "\t", '22', '0', '0.5';
    say $out join "\t", '22', '23', '0.5';
    say $out join "\t", '23', '24', '1.0';
    say $out join "\t", '24', '0', '1.0';
    say $out join "\t", '25', '0', '0.5';
    say $out join "\t", '25', '26', '0.5';
    say $out join "\t", '26', '27', '1.0';
    say $out join "\t", '27', '0', '1.0';
    say $out join "\t", '28', '0', '0.5';
    say $out join "\t", '28', '29', '0.5';
    say $out join "\t", '29', '30', '1.0';
    say $out join "\t", '30', '0', '1.0';
    say $out join "\t", '31', '0', '1.0';
    say $out join "\t", '32', '0', '1.0';
    say $out "Pi= 33";
    say $out "0.1";
    print $out "0.05\n" x 31;
    print $out '0.05';
    close $out;

    return $chr_file;
}

sub build_mgescan {
    my $self = shift;
    my ($bindir) = @_;

    my $root = $self->basedir->absolute->resolve;
    my $wd   = $self->workingdir->absolute->resolve;

    my $src_dir = File::Spec->catdir($wd, 'src');
    my $mgexe   = 'tephra-MGEScan';
    my $trexe   = 'tephra-translate';
    my $mgescan = File::Spec->catfile($bindir, $mgexe);
    my $transla = File::Spec->catfile($bindir, $trexe);

    chdir $src_dir;
    system("make clean -f mgescan-makefile 2>&1 >/dev/null") == 0 
	or die "make failed: $!";
    system("make -f mgescan-makefile 2>&1 >/dev/null") == 0 
	or die "make failed: $!";
    system("make clean -f translate-makefile 2>&1 >/dev/null") == 0
        or die "make failed: $!";
    system("make all -f translate-makefile 2>&1 >/dev/null") == 0
        or die "make failed: $!";
    
    copy $mgexe, $bindir or die "Copy failed: $!";
    copy $trexe, $bindir or die "Copy failed: $!";
    chmod 0755, $mgescan, $transla;
    
    return ($mgescan, $transla);
}

sub fetch_file {
    my $self = shift;
    my ($file, $endpoint) = @_;

    unless (-e $file) {
	my $response = HTTP::Tiny->new->get($endpoint);
	unless ($response->{success}) {
	    die "Can't get url $endpoint -- Status: ", $response->{status}, 
	        " -- Reason: ", $response->{reason};
	}
	open my $out, '>', $file;
	print $out $response->{content};
	#sleep 1;
	close $out;
    }
}

=head1 AUTHOR

S. Evan Staton, C<< <evan at evanstaton.com> >>

=head1 BUGS

Please report any bugs or feature requests through the project site at 
L<https://github.com/sestaton/tephra/issues>. I will be notified,
and there will be a record of the issue. Alternatively, I can also be 
reached at the email address listed above to resolve any questions.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Tephra::Config::Install


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut 

__PACKAGE__->meta->make_immutable;

1;
