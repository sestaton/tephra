package Tephra::Config;

use 5.010;
use Moose;
use Cwd;
use Config;
use File::Spec;
use File::Find;
use File::Copy qw(copy move);
use File::Path qw(make_path);
use File::Basename;
use Path::Class::File;
use HTML::TreeBuilder;
use HTTP::Tiny;
use Log::Any qw($log);
use namespace::autoclean;

has basedir => (
      is       => 'ro',
      isa      => 'Path::Class::Dir',
      required => 0,
      coerce   => 1,
      default  => sub {
	return Path::Class::Dir->new($ENV{HOME})
    },
);


sub get_gt_exes {
    my $host = 'http://genometools.org';
    my $dir  = 'pub/binary_distributions';
    my $file = 'gt_distlisting.html';
    fetch_file($file, $host."/".$dir);
    
    my $tree = HTML::TreeBuilder->new;
    $tree->parse_file($file);
    
    my ($dist, $ldist, $ldir);
    for my $tag ($tree->look_down(_tag => 'a')) {
	if ($tag->attr('href')) {
	    if ($tag->as_text =~ /Linux_x86_64-64bit-barebone.tar.gz\z/) {
		$dist = $tag->as_text;
		my $archive = join "/", $host, $dir, $dist;
		fetch_file($dist, $archive);
		
		$ldist = $dist;
		$ldist =~ s/\.tar.gz\z//;
		$ldir = 'gt';
		
		system("tar xzf $dist") == 0 or die $!;
		
		move $ldist, $ldir or die "Move failed: $!";
		unlink $dist;
	    }
	}
    }
    unlink $file;

    my $cwd = getcwd();
    my $gt  = File::Spec->catfile($cwd, $ldir, 'bin', 'gt');
    
    return $gt
}

sub get_hscan {
    my $host = 'http://sourceforge.net';
    my $dir  = 'projects/helitronscanner/files/HelitronScanner_V1.0.zip/download';
    my $ldir = 'helitronscanner';
    make_path( $ldir, {verbose => 0, mode => 0771,} );
    my $file = 'HelitronScanner.zip';
    my $path = File::Spec->catfile($ldir, $file);
    fetch_file($path, $host."/".$dir);
    chdir $ldir or die $!;
    system("unzip $file 2>&1 > /dev/null") == 0 or die $!;

    my $cwd   = getcwd();
    my $hscan = File::Spec->catfile($cwd, 'HelitronScanner', 'HelitronScanner.jar');

    return $hscan;
}

sub fetch_file {
    my ($file, $endpoint) = @_;
    unless (-e $file) {
	my $response = HTTP::Tiny->new->get($endpoint);
	unless ($response->{success}) {
	    die "Can't get url $endpoint -- Status: ", $response->{status}, " -- Reason: ", $response->{reason};
	}
	open my $out, '>', $file;
	print $out $response->{content};
	#sleep 1;
	close $out;
    }
}


__PACKAGE__->meta->make_immutable;

1;
