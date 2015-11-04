package PathUtils;

use 5.010;
use Moose::Role;
use Config;
use File::Spec;
use IPC::System::Simple qw(capture);
use Carp 'croak';

sub find_hmmsearch {
    if (defined $ENV{HMMER2}) {
        my $hmmsearch = join "/", $ENV{HMMER2}, 'src', 'hmmsearch';
        if (-e $hmmsearch && -x $hmmsearch) {
            return $hmmsearch;
        }
        else {
            $hmmsearch = join "/", $ENV{HMMER2}, 'bin', 'hmmsearch';
            if (-e $hmmsearch && -x $hmmsearch) {
                return $hmmsearch;
            }
        }
    }
    else {
	my @path = split /:|;/, $ENV{PATH};
	
	for my $p (@path) {
	    my $hmmsearch = File::Spec->catfile($p, 'hmmsearch');
	    if (-e $hmmsearch && -x $hmmsearch) {
		my @out = capture([0..5], "hmmsearch", "-h");
		my ($version) = grep { /HMMER/ } @out;              
		if ($version =~ /HMMER (\d\.\d\w?\d+?) \(/) {
		    my $release = $1;                  
		    if ($release =~ /^2/) {
			return $hmmsearch;
		    }
		    else {
			croak "\nERROR: HMMER version 2 is required but was not found. Exiting.\n";
		    }
		}
	    }
	}
    }
}

1;
