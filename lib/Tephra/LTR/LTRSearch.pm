package Tephra::LTR::LTRSearch;

use 5.014;
use Moose;
use File::Spec;
use File::Find;
use File::Basename;
use File::Path          qw(remove_tree);
use IPC::System::Simple qw(system EXIT_ANY);
use Cwd                 qw(abs_path);
use Log::Any            qw($log);
use Try::Tiny;
use YAML::Tiny;
use namespace::autoclean;
#use Data::Dump;

with 'Tephra::Role::Logger',
     'Tephra::Role::Run::GT';

=head1 NAME

Tephra::LTR::LTRSearch - Find LTR retrotransposons in a reference genome

=head1 VERSION

Version 0.09.3

=cut

our $VERSION = '0.09.3';
$VERSION = eval $VERSION;

has config => (
    is            => 'ro',
    isa           => 'Str',
    required      => 0,
    documentation => qq{The Tephra LTR configuration file},
);

has logfile => (
      is        => 'ro',
      isa       => 'Str',
      predicate => 'has_logfile',
);

sub ltr_search {
    my $self = shift;
    my ($search_obj) = @_;
    my ($config, $index, $mode) = @{$search_obj}{qw(config index mode)};
    
    my $genome  = $self->genome->absolute->resolve;
    my $hmmdb   = $self->hmmdb->absolute->resolve;
    my $trnadb  = $self->trnadb->absolute->resolve;
    my $logfile = $self->logfile;

    ## LTRharvest constraints
    my ($overlaps, $mintsd, $maxtsd, $minlenltr, $maxlenltr, $mindistltr, $maxdistltr, $pdomcutoff, $pdomevalue) = 
	@{$config->{findltrs}}{qw(overlaps mintsd maxtsd minlenltr maxlenltr mindistltr maxdistltr pdomcutoff pdomevalue)};

    my ($seedlength, $tsdradius, $xdrop, $swmat, $swmis, $swins, $swdel) = 
	@{$config->{findltrs}}{qw(seedlength tsdradius xdrop swmat swmis swins swdel)};

    ## LTRdigest constraints
    my ($pptradius, $pptlen, $pptagpr, $uboxlen, $uboxutpr, $pbsradius, $pbslen, $pbsoffset, $pbstrnaoffset, $pbsmaxeditdist, $maxgaplen) = 
	@{$config->{findltrs}}{qw(pptradius pptlen pptagpr uboxlen uboxutpr pbsradius pbslen pbsoffset pbstrnaoffset pbsmaxeditdist maxgaplen)};

    my (%ltrh_cmd, %ltrd_cmd, @ltrh_opts, @ltrh_args, $ltrh_gff, $ltrd_gff);
    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
	$name =~ s/$1//;
    }
    
    if ($mode eq 'strict') {
	$ltrh_gff = File::Spec->catfile( abs_path($path), $name.'_ltrharvest99.gff3' );
	$ltrd_gff = File::Spec->catfile( abs_path($path), $name.'_ltrdigest99.gff3' );

	@ltrh_opts = qw(-seqids -mintsd -maxtsd -minlenltr -maxlenltr -mindistltr 
                       -maxdistltr -motif -similar -vic -index -overlaps -seed -vic 
                       -xdrop -mat -mis -ins -del -gff3);
	@ltrh_args = ("yes",$mintsd,$maxtsd,$minlenltr,$maxlenltr,$mindistltr,$maxdistltr,"tgca","99",
		      "10",$index,$overlaps,$seedlength,$tsdradius,$xdrop,
		      $swmat,$swmis,$swins,$swdel,$ltrh_gff);
    }
    elsif ($mode eq 'relaxed') {
	$ltrh_gff = File::Spec->catfile( abs_path($path), $name.'_ltrharvest85.gff3' );
	$ltrd_gff = File::Spec->catfile( abs_path($path), $name.'_ltrdigest85.gff3' );

	@ltrh_opts = qw(-seqids -mintsd -maxtsd -minlenltr -maxlenltr 
                       -mindistltr -maxdistltr -similar -vic -index -overlaps
                       -seed -vic -xdrop -mat -mis -ins -del -gff3);
	@ltrh_args = ("yes",$mintsd,$maxtsd,$minlenltr,$maxlenltr,$mindistltr,$maxdistltr,"85","10",
		      $index,$overlaps,$seedlength,$tsdradius,$xdrop,$swmat,
		      $swmis,$swins,$swdel,$ltrh_gff);
    }
    else {
	say STDERR "\nERROR: Could not get 'mode' for LTR search. This is a bug, please report it. Exiting.\n";
        exit(1);
    }

    my $log = $self->get_tephra_logger($logfile);

    @ltrh_cmd{@ltrh_opts} = @ltrh_args;
    my $ltrh_succ = $self->run_ltrharvest(\%ltrh_cmd, $log);

    if ($ltrh_succ && -s $ltrh_gff) {
	my $gffh_sort = $self->sort_gff($ltrh_gff, $log);

	my @ltrd_opts = qw(-trnas -hmms -seqfile -matchdescstart -seqnamelen -o 
                           -pdomevalcutoff -pdomcutoff -pptradius -pptlen -pptaprob 
                           -pptgprob -uboxlen -pptuprob -pbsalilen -pbsradius -pbsoffset -pbstrnaoffset
                           -pbsmaxedist -maxgaplen);
	my @ltrd_args = ($trnadb,$hmmdb,$genome,"yes","50",$ltrd_gff,
			 $pdomevalue,$pdomcutoff,$pptradius,$pptlen,$pptagpr,$pptagpr,$uboxlen,
	                 $uboxutpr,$pbslen,$pbsradius,$pbsoffset,$pbstrnaoffset,$pbsmaxeditdist,$maxgaplen);
	@ltrd_cmd{@ltrd_opts} = @ltrd_args;

	my $ltrd_succ = $self->run_ltrdigest(\%ltrd_cmd, $gffh_sort, $log);
	$self->clean_indexes($path) 
	    if $self->clean && $mode eq 'relaxed';
	unlink $ltrh_gff;
	unlink $gffh_sort;
	
	return $ltrd_gff;
    }
    else {
	$self->clean_indexes($path) 
	    if $self->clean && $mode eq 'relaxed';
	unlink $ltrh_gff;

	return 0;
    }
}

=head1 AUTHOR

S. Evan Staton, C<< <evan at evanstaton.com> >>

=head1 BUGS

Please report any bugs or feature requests through the project site at 
L<https://github.com/sestaton/tephra/issues>. I will be notified,
and there will be a record of the issue. Alternatively, I can also be 
reached at the email address listed above to resolve any questions.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Tephra::LTR::LTRSearch


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
