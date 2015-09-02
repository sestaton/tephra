package Tephra::Role::GT;

use 5.010;
use Moo::Role;
use File::Spec;
use Types::Standard qw(Str);
use Log::Any        qw($log);

has gt_exec => (
    is        => 'rw',
    isa       => 'Str',
    reader    => 'get_gt_exec',
    writer    => 'set_gt_exec',
    predicate => 'has_gt_exec',
);

sub find_gt {
    my $self  = shift;
    my $gtexe = $self->get_gt_exec; # set during class initialization

    # check if the executable path was set first
    if (defined $gtexe && -e $gtexe && -x $gtexe) {
        return ($gtexe);
    }
    elsif (! -e $gtexe) {
	my @path = split /:|;/, $ENV{PATH};

    for my $p (@path) {
        my $gt = File::Spec->catfile($p, 'gt');

    if (-e $gt && -x $gt) {
        return ($gt);
    }
    else {
        $log->error("Unable to find 'gt' executable. Make sure genometools is installed. Exiting.");
        exit(1);
    }
}
