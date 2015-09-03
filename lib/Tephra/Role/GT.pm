package Tephra::Role::GT;

use 5.010;
use Moose::Role;
use MooseX::Types::Path::Class;
use File::Spec;
#use Types::Standard qw(Str);
use Log::Any        qw($log);

has gt_exec => (
    is        => 'rw',
    isa       => 'Path::Class::File',
    required  => 0,
    coerce    => 1, 
    reader    => 'get_gt_exec',
    writer    => 'set_gt_exec',
    predicate => 'has_gt_exec',
    builder   => '_build_gt_exec',
);

sub _build_gt_exec {
    my $self  = shift;
    my $gtexe = $self->get_gt_exec; # set during class initialization
    #my $gtexe = '/home';
    
    # check if the executable path was set first
    if (defined $gtexe && -e $gtexe && -x $gtexe) {
        return $gtexe;
    }
    elsif (! defined $gtexe) {
	my @path = split /:|;/, $ENV{PATH};
	for my $p (@path) {
	    my $gt = File::Spec->catfile($p, 'gt');

	    if (-e $gt && -x $gt) {
		$self->set_gt_exec($gt);
		return $gt;
	    }
	}
    }
    else {
        #$log->error(
	say STDERR "Unable to find 'gt' executable. Make sure genometools is installed. Exiting."; #);
        exit(1);
    }
}

1;
