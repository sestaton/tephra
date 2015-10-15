package Tephra::NonLTR::NonLTRSearch;

use 5.010;
use Moose;
use Cwd;
use File::Spec;
use File::Find;
use File::Basename;
use File::Path          qw(make_path);
use IPC::System::Simple qw(system EXIT_ANY);
use Log::Any            qw($log);
use Try::Tiny;
use HTTP::Tiny;
use WWW::Mechanize;
use namespace::autoclean;

#with 'Tephra::Role::Run::GT';

sub find_nonltrs {
    my $self = shift;
 
    my $genome = $self->genome->absolute;
    my $outdir = $self->outdir;

    $self->fetch_hmm_models;

    
}

sub fetch_hmm_models {
    my $self = shift;
    my $outdir = $self->outdir;

    my $phmm_dir = File::Spec->catdir($outdir, 'pHMM');
    unless ( -d $phmm_dir ) {
	make_path( $phmm_dir, {verbose => 0, mode => 0771,} );
    }

    my $url = 'https://github.com/MGEScan/mgescan/tree/master/mgescan/nonltr/pHMM';
    my $urlbase = 'https://github.com';

    my %links;
    my $mech = WWW::Mechanize->new();
    $mech->get( $url );
    my @links = $mech->links();
    for my $link ( @links ) {
	next unless defined $link->text;
	if ($link->url =~ /.hmm3$/) {
	    $links{$link->text} = $urlbase.$link->url;
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


__PACKAGE__->meta->make_immutable;

1;
