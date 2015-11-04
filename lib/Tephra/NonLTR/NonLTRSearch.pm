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


__PACKAGE__->meta->make_immutable;

1;
