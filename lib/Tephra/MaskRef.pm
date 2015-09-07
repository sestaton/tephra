package Tephra::MaskRef;;

use 5.010;
use Moose;
use Cwd;
use File::Spec;
use File::Find;
use File::Basename;
use IPC::System::Simple qw(system EXIT_ANY);
use Sort::Naturally;
use Path::Class::File;
use Log::Any            qw($log);
use Try::Tiny;
use namespace::autoclean;

with 'Tephra::Role::Util';

has genome => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 1,
      coerce   => 1,
);

has repeatdb => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 0,
      coerce   => 1,
);

has clean => (
    is       => 'ro',
    isa      => 'Bool',
    required => 0,
    default  => 1,
    );

sub mask_reference {
    my $self = shift;
    my $genome   = $self->genome->absolute;
    my $repeatdb = $self->repeatdb;
    my (%suf_args, %ltrh_cmd, %ltrd_cmd);

    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
	$name =~ s/$1//;
    }

    my $masked_ref = File::Spec->catfile($path, $name."_masked.fas");
    my $index      = $repeatdb.".index";
    my $vmatch_log = File::Spec->catfile($path, $name."_vmatch.err");
    my $mkvtree = "mkvtree -db $repeatdb -indexname $index -dna -allout -v -pl 2>&1 > /dev/null";
    my $vmatch  = "vmatch -q $genome -l 50 -s 90 -qmaskmatch N $index 1> $masked_ref 2> $vmatch_log";
    $self->run_cmd($mkvtree);
    $self->run_cmd($vmatch);

    $self->clean_index($path) if $self->clean;

    return $masked_ref;
}

sub clean_index {
    my $self = shift;
    my ($dir) = @_;
    
    #my $dir = getcwd();
    my @files;
    find( sub { push @files, $File::Find::name
		    if /\.al1|\.llv|\.ssp|\.bck|\.ois|\.sti1|\.bwt|\.prj|\.suf|\.des|\.sds|\.tis|\.lcp|\.skp/
	  }, $dir);
    unlink @files;
}


__PACKAGE__->meta->make_immutable;

1;
