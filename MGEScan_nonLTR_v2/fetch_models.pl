use 5.010;
use strict;
use warnings;
use File::Spec;
use File::Basename;
use WWW::Mechanize;
use HTTP::Tiny;

my $usage = basename($0)." outdir";
my $dir   = shift or die $usage;

fetch_hmm_models($dir);

###
sub fetch_hmm_models {
    my $outdir = shift;

    my $phmm_dir = File::Spec->catdir($outdir, 'pHMM');
    unless ( -d $phmm_dir ) {
	make_path( $phmm_dir, {verbose => 0, mode => 0771,} );
    }

    my $url = 'https://github.com/MGEScan/mgescan/tree/master/mgescan/nonltr/pHMM';
    my $urlbase = 'https://raw.githubusercontent.com/MGEScan/mgescan/master/mgescan/nonltr/pHMM'; 
    
    my %links;
    my $mech = WWW::Mechanize->new();
    $mech->get( $url );
    my @links = $mech->links();
    for my $link ( @links ) {
	next unless defined $link->text;
	if ($link->url =~ /.hmm3$/) {
	    $links{$link->text} = join "/", $urlbase, $link->text;
	}
    }
    
    for my $model (keys %links) {
	my $response = HTTP::Tiny->new->get($links{$model});

	unless ($response->{success}) {
	    die "Can't get url $url -- Status: ", $response->{status}, " -- Reason: ", $response->{reason};
	}

	my $outfile = File::Spec->catfile($outdir, $phmm_dir, $model);
	open my $out, '>', $outfile;
	say $out $response->{content};
	close $out;
    }
}
