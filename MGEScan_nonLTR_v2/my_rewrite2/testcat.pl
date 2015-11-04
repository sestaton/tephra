use strict;
use warnings;
use File::Spec;
use Data::Printer;

my $seqdir = File::Spec->catdir('data', 'info', 'full');
my @all_clade = ('CR1', 'I', 'Jockey', 'L1', 'L2', 'R1', 'RandI', 'Rex', 'RTE', 'Tad1', 'R2','CRE');
my @clades = File::Spec->catdir($seqdir, $_), map @all_clade;
p @clades;
