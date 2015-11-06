 package Tephra::Config;

 use 5.010;
 use Moose;
 use MooseX::Types::Path::Class;
 use Cwd;
 use Config;
 use File::Spec;
 use File::Find;
 use File::Copy qw(copy move);
 use File::Path qw(make_path remove_tree);
 use File::Basename;
 use Path::Class::File;
 use HTML::TreeBuilder;
 use HTTP::Tiny;
 use Log::Any qw($log);
 use namespace::autoclean;

=head1 NAME

Tephra::Config - Class for setting up Tephra dependencies

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

has basedir => (
    is       => 'ro',
    isa      => 'Path::Class::Dir',
    required => 0,
    coerce   => 1,
    default  => sub {
	return Path::Class::Dir->new($ENV{HOME})
    },
);

has workingdir => ( 
    is       => 'ro', 
    isa      => 'Path::Class::Dir', 
    required => 1, 
    coerce   => 1 
);

sub configure_root {
    my $self = shift;
    my $root = $self->basedir;
    
    unless (-e $root) {
	make_path($root, {verbose => 0, mode => 0711,});
    }
    
    # we don't want to reconfigure every time the tests run
    my $gt      = File::Spec->catfile($root, 'gt', 'bin', 'gt');
    my $hscan   = File::Spec->catfile($root, 'helitronscanner', 'HelitronScanner', 'HelitronScanner.jar');
    my $hmmbin  = File::Spec->catdir($root, 'hmmer-2.3.2', 'bin');
    my $moddir  = File::Spec->catdir($root, 'pHMM');
    my $chrdir  = File::Spec->catdir($root, 'hmm');
    my $mgescan = File::Spec->catfile($chrdir, 'tephra-MGEScan');
    my $transla = File::Spec->catfile($chrdir, 'tephra-translate');
    my $clw     = File::Spec->catdir($root, 'clustalw-2.1', 'bin', 'clustalw2');
    my $pamlbin = File::Spec->catdir($root, 'paml4.8', 'bin');
    my $transeq = File::Spec->catdir($root, 'EMBOSS-6.5.7', 'bin');
    
    unless (-e $gt) {
	$gt = $self->get_gt_exes;
    }
    
    unless (-e $hscan) {
	$hscan = $self->get_hscan;
    }
    
    unless (-e $hmmbin) {
	$hmmbin = $self->fetch_hmmer2;
    }
    
    unless (-e $moddir) {
	$moddir = $self->fetch_hmm_models;
    }
    
    unless (-e $chrdir) {
	$chrdir = $self->make_chrom_dir;
    }
    
    unless (-e $mgescan && -e $transla) {
	($mgescan, $transla) = $self->build_mgescan;
    }
    
    unless (-e $clw) {
	$clw = $self->fetch_clustalw2;
    }
    
    unless (-e $pamlbin) {
	$pamlbin = $self->fetch_paml;
    }
    
    unless (-e $transeq) {
	$transeq = $self->fetch_paml;
    }

    return ({
	gt       => $gt, 
	hscandir => $hscan, 
	hmmerbin => $hmmbin, 
	modeldir => $moddir, 
	hmmdir   => $chrdir, 
	mgescan  => $mgescan, 
	transcmd => $transla, 
	clustalw => $clw, 
	pamlbin  => $pamlbin,
	transeq  => $transeq });
}

sub get_gt_exes {
    my $self = shift;
    my $root = $self->basedir;
    my $wd   = $self->workingdir;
    
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
		
		move $ldist, $ldir or die "Move failed: $!";
		unlink $dist;
	    }
	}
    }
    unlink $file;
    
    my $gt = File::Spec->catfile($ldir, 'bin', 'gt');

    return $gt
}

sub get_hscan {
    my $self = shift;
    my $root = $self->basedir;
    my $wd   = $self->workingdir;
    
    my $host = 'http://sourceforge.net';
    my $dir  = 'projects/helitronscanner/files/HelitronScanner_V1.0.zip/download';
    my $ldir = File::Spec->catdir($root, 'helitronscanner');
    make_path( $ldir, {verbose => 0, mode => 0771,} );
    my $file = 'HelitronScanner.zip';
    my $path = File::Spec->catfile($ldir, $file);
    $self->fetch_file($path, $host."/".$dir);
    chdir $ldir or die $!;
    system("unzip $file 2>&1 > /dev/null") == 0 or die $!;
    
    my $cwd   = getcwd();
    my $hscan = File::Spec->catfile($cwd, 'HelitronScanner', 'HelitronScanner.jar');
    chdir $wd;
    
     return $hscan;
}

sub fetch_hmmer2 {
    my $self = shift;
    my $root = $self->basedir;
    my $wd   = $self->workingdir;
    
    my $urlbase = 'http://selab.janelia.org'; 
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
    system("./configure --prefix=$cwd 2>&1 > /dev/null") == 0
	or die "configure failed: $!";
    system("make -j4 2>&1 >/dev/null") == 0 
	or die "make failed: $!";
    system("make install 2>&1 >/dev/null") == 0
	 or die "make failed: $!";
    my $hmmbin = File::Spec->catdir($cwd, 'bin');
    my $distfile = File::Spec->catfile($root, $file);
    unlink $distfile;
    chdir $wd;
    
    return $hmmbin;
}

sub fetch_clustalw2 {
    my $self = shift;
    my $root = $self->basedir;
    my $wd   = $self->workingdir;
    
    my $urlbase = 'http://www.clustal.org';
    my $dir     = 'download';
    my $tool    = 'current';
    my $file    = 'clustalw-2.1.tar.gz';
    my $url     = join "/", $urlbase, $dir, $tool, $file;
    my $outfile = File::Spec->catfile($root, $file);
    $self->fetch_file($outfile, $url);

    chdir $root;
    my $dist = 'clustalw-2.1';
    system("tar xzf $file") == 0 or die "tar failed: $!";
    chdir $dist;
    my $cwd = getcwd();
    system("./configure --prefix=$cwd 2>&1 > /dev/null") == 0
	or die "configure failed: $!";
    system("make -j4 2>&1 > /dev/null") == 0 
	or die "make failed: $!";
    system("make install 2>&1 > /dev/null") == 0
	or die "make failed: $!";
    
    my $clw = File::Spec->catdir($cwd, 'bin', 'clustalw2');
    my $distfile = File::Spec->catfile($root, $file);
    unlink $distfile;
    chdir $wd;

    return $clw;
}

sub fetch_paml {
    my $self = shift;
    my $root = $self->basedir;
    my $wd   = $self->workingdir;

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
    system("make -j4 2>&1 >/dev/null") == 0 
	or die "make failed: $!";

    my @exelist = ('yn00', 'baseml', 'basemlg', 'mcmctree', 'pamp', 'evolver', 'infinitesites', 'codeml');
    my $rootdir = File::Spec->catdir($root, 'paml4.8', 'bin');
    unless ( -d $rootdir ) {
	make_path( $rootdir, {verbose => 0, mode => 0771,} );
    }

    for my $file (@exelist) {
	copy $file, $rootdir or die "Copy failed: $!";
    }

    my @nonexes = map { File::Spec->catfile($rootdir, $_) } @exelist;
    my $cnt = chmod 0755, @nonexes;
    
    if ($cnt == @exelist) {
	return $rootdir;
    }
}

sub fetch_emboss {
    my $self = shift;
    my $root = $self->basedir;
    my $wd   = $self->workingdir;

    my @path = split /:|;/, $ENV{PATH};
    
    for my $p (@path) {
	my $transeq  = File::Spec->catfile($p, 'transeq');
	if (-e $transeq && -x $transeq) {
	    return $transeq;
	}
    }

    my $urlbase = 'ftp://emboss.open-bio.org';
    my $dir     = 'pub';
    my $tool    = 'EMBOSS';
    my $release = 'old';
    my $version = '6.5.0';
    my $file    = 'EMBOSS-6.5.7.tar.gz';
    my $url     = join "/", $urlbase, $dir, $tool, $release, $version, $file;
    my $outfile = File::Spec->catfile($root, $file);
    $self->fetch_file($outfile, $url);

    chdir $root;
    my $dist = 'EMBOSS-6.5.7';
    system("tar xzf $file") == 0 or die "tar failed: $!";
    chdir $dist;
    my $cwd = getcwd();
    system("./configure --prefix=$cwd 2>&1 > /dev/null") == 0
	or die "configure failed: $!";
    system("make -j4 2>&1 > /dev/null") == 0 
	or die "make failed: $!";
    system("make install 2>&1 > /dev/null") == 0
	or die "make failed: $!";
    
    my $transeq = File::Spec->catdir($cwd, 'bin', 'transeq');
    my $distfile = File::Spec->catfile($root, $file);
    unlink $distfile;
    chdir $wd;

    return $transeq;
}

sub fetch_hmm_models {
    my $self = shift;
    my $root = $self->basedir;
    my $wd   = $self->workingdir;
   
    chdir $wd;
    my $file = 'pHMM.tar.gz';
    my $dist = File::Spec->catfile('build', $file);
    copy $dist, $root or die "Copy failed: $!";
    chdir $root;
    system("tar xzf $file") == 0 or die $!;
    unlink $file;

    my $dir = File::Spec->catfile($root, 'pHMM');

    return $dir;
}

sub make_chrom_dir {
    my $self = shift;
    my $root = $self->basedir;

    my $hmm_dir = File::Spec->catdir($root, 'hmm');
    unless ( -d $hmm_dir ) {
        make_path( $hmm_dir, {verbose => 0, mode => 0771,} );
    }

    my $chr_file = File::Spec->catfile($hmm_dir, 'chr.hmm');
    open my $out, '>', $chr_file;
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

    return $hmm_dir;
}

sub build_mgescan {
    my $self = shift;
    my $root = $self->basedir;
    my $wd   = $self->workingdir;

    my $hmm_dir = File::Spec->catdir($root, 'hmm');
    my $src_dir = File::Spec->catdir($wd, 'src');
    my $mgexe   = 'tephra-MGEScan';
    my $trexe   = 'tephra-translate';
    my $mgescan = File::Spec->catfile($hmm_dir, $mgexe);
    my $transla = File::Spec->catfile($hmm_dir, $trexe);

    chdir $src_dir;
    system("make clean -f mgescan-makefile 2>&1 >/dev/null") == 0 
	or die "make failed: $!";
    system("make -f mgescan-makefile 2>&1 >/dev/null") == 0 
	or die "make failed: $!";
    system("make clean -f translate-makefile 2>&1 >/dev/null") == 0
        or die "make failed: $!";
    system("make all -f translate-makefile 2>&1 >/dev/null") == 0
        or die "make failed: $!";
    
    copy $mgexe, $hmm_dir or die "Copy failed: $!";
    copy $trexe, $hmm_dir or die "Copy failed: $!";
    my $cnt = chmod 0755, $mgescan, $transla;
    
    if ($cnt == 2 && -e $mgescan && -e $transla) {
	return ($mgescan, $transla);
    }
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

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests through the project site at 
L<https://github.com/sestaton/tephra/issues>. I will be notified,
and there will be a record of the issue. Alternatively, I can also be 
reached at the email address listed above to resolve any questions.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Tephra::Config


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut 

__PACKAGE__->meta->make_immutable;

1;
