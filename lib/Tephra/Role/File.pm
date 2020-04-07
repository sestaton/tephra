package Tephra::Role::File;

use 5.010;
use Moose::Role;
use MooseX::Types::Path::Class;
#use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::File;
#use Symbol;
#use Method::Signatures;

=head1 NAME

Tephra::Role::File - File handling methods for Tephra.

=head1 VERSION

Version 0.12.6

=cut

our $VERSION = '0.12.6';
$VERSION = eval $VERSION;

=head1 SYNOPSIS

    use Tephra;

    with 'Tephra::Role::File'
    ...

=cut

#has 'config' => (
#    is            => 'ro',
#    isa           => 'Str',
#    required      => 0,
#    documentation => qq{The Tephra configuration file},
#);

#has 'config' => (
#      is       => 'ro',
#      isa      => 'Path::Class::File',
#      required => 0,
#      coerce   => 1,
#);

has 'file' => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 0,
      coerce   => 1,
);

has 'dir' => (
      is       => 'ro',
      isa      => 'Path::Class::Dir',
      required => 0,
      coerce   => 1,
);

has 'fh' => (
    is         => 'ro',
    #isa        => 'IO::File',
    predicate  => 'has_fh',
    lazy_build => 1,
    builder    => '_build_fh',
);

has 'format' => (
    is        => 'ro',
    isa       => 'Str',
    predicate => 'has_format',
    default   => 'fasta'
);

sub _build_fh {
    my $self = shift;
    my $file = $self->file->absolute;
    my $fh = IO::File->new();

    # make sure zcat/bzcat can be found under different shells
    my @path = split /:|;/, $ENV{PATH};
    unless (@path) {
        $ENV{PATH} = "/usr/bin:/bin";
    }

    if ($file =~ /\.gz$/) {
        open $fh, '-|', 'zcat', $file or die "\nERROR: Could not open file: $file\n";
	#$fh = new IO::Uncompress::Gunzip $file->stringify;
	    #or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
    }
    elsif ($file =~ /\.bz2$/) {
        open $fh, '-|', 'bzcat', $file or die "\nERROR: Could not open file: $file\n";
    }
    elsif ($file =~ /^-$|STDIN/) {
        open $fh, '< -' or die "\nERROR: Could not open STDIN\n";
    }
    else {
	open $fh, '<', $file or die "\nERROR: Could not open file: $file\n";
    }

    return $fh;
}

=head1 METHODS

=head2 get_fh

 Title   : get_fh

 Usage   : my $fh = $trans_obj->file->get_fh;
          
 Function: Gets a filehandle for the associated
           file.

                                                   Return_type
 Returns : An open filehandle for reading          Scalar

 Args    : None. This is a role that can
           be consumed.

=cut

sub get_fh {
    my $self = shift;
    my ($file) = @_;

    # make sure sort can be found under different shells
    my @path = split /:|;/, $ENV{PATH};
    unless (@path) {
        $ENV{PATH} = "/usr/bin:/bin";
    }

    my $fh;
    if ($file =~ /\.gz$/) {
	open $fh, '-|', 'zcat', $file or die "\nERROR: Could not open file: $file: $!\n";
    }
    elsif ($file =~ /\.bz2$/) {
	open $fh, '-|', "bzcat <$file" or die "\nERROR: Could not open file: $file: $!\n";
    }
    elsif ($file =~ /^-$|STDIN/) {
	open $fh, '< -' or die "\nERROR: Could not open STDIN\n";
    }
    else {
	$fh = $self->file->openr;
    }
    return $fh;
}

=head2 collate

 Title   : collate

 Usage   : $tephra_obj->collate($file, $fh_out);
          
 Function: Takes an input file and an output filehandle 
           and appends the file contents to the open file.

                                                   Return_type
 Returns : None. No data is returned.

                                                   Arg_type
 Args    : file - The input file to append to      Scalar
                  the output
           fh_out - The output file handle to      Scalar
                    use for writing the file
                    contents

=cut

sub collate {
    my $self = shift;
    my ($file_in, $fh_out) = @_;

    my $lines = do { 
        local $/ = undef; 
        open my $fh_in, '<', $file_in or die "\n[ERROR]: Could not open file: $file_in\n";
        <$fh_in>;
    };
    print $fh_out $lines;

    return;
}

=head1 AUTHOR

S. Evan Staton, C<< <evan at evanstaton.com> >>

=head1 BUGS

Please report any bugs or feature requests through the project site at 
L<https://github.com/sestaton/Tephra/issues>. I will be notified,
and there will be a record of the issue. Alternatively, I can also be 
reached at the email address listed above to resolve any questions.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Tephra::Role::File


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2013-2019 S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1;
