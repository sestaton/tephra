package Tephra::NonLTR::Utils;

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

sub fetch_hmmer2 {
    my $self = shift;
    my $outdir = $self->outdir;
    
    my $urlbase = 'http://selab.janelia.org'; 
    my $dir     = 'software';
    my $tool    = 'hmmer';
    my $version = '2.3.2';
    my $file    = 'hmmer-2.3.2.tar.gz';
    my $url     = join "/", $urlbase, $dir, $tool, $version, $file;
    my $response = HTTP::Tiny->new->get($url);

    unless ($response->{success}) {
	die "Can't get url $url -- Status: ", $response->{status}, " -- Reason: ", $response->{reason};
    }

    my $outfile = File::Spec->catfile($outdir, $file);
    open my $out, '>', $outfile;
    say $out $response->{content};
    close $out;

    my $ldir = 'hmmer-2.3.2';
    system("tar xzf $file") == 0 or die "tar failed: $!";
    chdir $ldir;
    system("./configure && make -j4 2>&1 >/dev/null") == 0 or die "configure failed: $!";

    return $ldir;
}

sub fetch_hmm_models {
    my $self = shift;
    my $outdir = $self->outdir;

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
        if ($link->url =~ /.hmm$/) {
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


__PACKAGE__->meta->make_immutable;

1;
