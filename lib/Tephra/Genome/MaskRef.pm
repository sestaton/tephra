package Tephra::Genome::MaskRef;

use 5.010;
use Moose;
use Cwd;
use File::Spec;
use File::Find;
use File::Basename;
#use Path::Class::File;
use Log::Any qw($log);
use namespace::autoclean;

with 'Tephra::Role::Util';

=head1 NAME

Tephra::Genome::MaskRef - Mask a reference with repeats to reduce false positives

=head1 VERSION

Version 0.03.5

=cut

our $VERSION = '0.03.5';
$VERSION = eval $VERSION;

has genome => (
      is       => 'ro',
      isa      => 'Maybe[Str]',
      required => 1,
      coerce   => 0,
);

has repeatdb => (
      is       => 'ro',
      isa      => 'Maybe[Str]',
      required => 0,
      coerce   => 0,
);

has outfile => (
    is       => 'ro',
    isa      => 'Maybe[Str]',
    required => 1,
    coerce   => 0,
);

has clean => (
    is       => 'ro',
    isa      => 'Bool',
    required => 0,
    default  => 1,
);

sub mask_reference {
    my $self = shift;
    my $genome   = $self->genome;
    my $repeatdb = $self->repeatdb;

    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
	$name =~ s/$1//;
    }

    my $outfile = $self->outfile // File::Spec->catfile($path, $name."_masked.fas");

    my (%suf_args, %ltrh_cmd, %ltrd_cmd);
    my $index      = $repeatdb.".index";
    my $vmatch_log = File::Spec->catfile($path, $name."_vmatch.err");
    my $mkvtree = "mkvtree -db $repeatdb -indexname $index -dna -allout -v -pl 2>&1 > /dev/null";
    my $vmatch  = "vmatch -q $genome -l 50 -s 90 -qmaskmatch N $index 1> $outfile 2> $vmatch_log";
    $self->run_cmd($mkvtree); # need to warn here, not just log errors
    $self->run_cmd($vmatch);

    $self->clean_index($path) if $self->clean;

    return $outfile;
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

=head1 AUTHOR

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests through the project site at 
L<https://github.com/sestaton/tephra/issues>. I will be notified,
and there will be a record of the issue. Alternatively, I can also be 
reached at the email address listed above to resolve any questions.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Tephra::Genome::MaskRef


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
