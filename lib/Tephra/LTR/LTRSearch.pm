package Tephra::LTR::LTRSearch;

use 5.010;
use Moose;
use Cwd;
use File::Spec;
use File::Find;
use File::Basename;
use IPC::System::Simple qw(system EXIT_ANY);
use Path::Class::File;
use Log::Any            qw($log);
use Try::Tiny;
use YAML::Tiny;
use namespace::autoclean;
#use Data::Dump;

with 'Tephra::LTR::Role::Config',
     'Tephra::Role::Run::GT';

=head1 NAME

Tephra::LTR::LTRSearch - Find LTR retrotransposons in a reference genome

=head1 VERSION

Version 0.02.3

=cut

our $VERSION = '0.02.3';
$VERSION = eval $VERSION;

has config => (
    is            => 'ro',
    isa           => 'Str',
    required      => 0,
    documentation => qq{The Tephra LTR configuration file},
);

sub ltr_search_strict {
    my $self = shift;
    my ($config, $index) = @_;
    
    my $genome = $self->genome;
    my $hmmdb  = $self->hmmdb;
    my $trnadb = $self->trnadb;

    ## LTR constraints
    my $overlaps   = $config->{overlaps};
    my $mintsd     = $config->{mintsd};
    my $maxtsd     = $config->{maxtsd};
    my $minlenltr  = $config->{minlenltr};
    my $maxlenltr  = $config->{maxlenltr};
    my $mindistltr = $config->{mindistltr};
    my $maxdistltr = $config->{maxdistltr};
    my $pdomcutoff = $config->{pdomcutoff};
    my $pdomevalue = $config->{pdomevalue};

    my $pptradius = $config->{pptradius};
    my $pptlen    = $config->{pptlen};
    my $pptagpr   = $config->{pptagpr};
    my $uboxlen   = $config->{uboxlen};
    my $uboxutpr  = $config->{uboxutpr};
    my $pbsradius = $config->{pbsradius};
    my $pbslen    = $config->{pbslen};
    my $pbsoffset = $config->{pbsoffset};
    my $pbstrnaoffset  = $config->{pbstrnaoffset};
    my $pbsmaxeditdist = $config->{pbsmaxeditdist};
    my $maxgaplen = $config->{maxgaplen};

    my $seedlength = $config->{seedlength};
    my $tsdradius  = $config->{tsdradius};
    my $xdrop      = $config->{xdrop};
    my $swmat      = $config->{swmat};
    my $swmis      = $config->{swmis};
    my $swins      = $config->{swins};
    my $swdel      = $config->{swdel};

    my (%suf_args, %ltrh_cmd, %ltrd_cmd);
    
    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
	$name =~ s/$1//;
    }
    
    my $ltrh_out = File::Spec->catfile($path, $name."_ltrharvest99_pred-all");
    my $ltrh_out_inner = File::Spec->catfile($path, $name."_ltrharvest99_pred-inner");
    my $ltrh_gff = File::Spec->catfile($path, $name."_ltrharvest99.gff3");
    my $ltrg_gff = File::Spec->catfile($path, $name."_ltrdigest99.gff3");
    my $ltrg_out = File::Spec->catfile($path, $name."_ltrdigest99");

    my @ltrh_opts = qw(-longoutput -seqids -tabout -mintsd -maxtsd -minlenltr -maxlenltr -mindistltr 
                       -maxdistltr -motif -similar -vic -index -out -outinner -overlaps -seed -vic 
                       -xdrop -mat -mis -ins -del -gff3);
    my @ltrh_args = ("no","yes","no",$mintsd,$maxtsd,$minlenltr,$maxlenltr,$mindistltr,$maxdistltr,"tgca","99",
		     "10",$index,$ltrh_out,$ltrh_out_inner,$overlaps,$seedlength,$tsdradius,$xdrop,
		     $swmat,$swmis,$swins,$swdel,$ltrh_gff);
    @ltrh_cmd{@ltrh_opts} = @ltrh_args;
    
    my $ltr_succ  = $self->run_ltrharvest(\%ltrh_cmd);
    if (-s $ltrh_gff > 1) {
	my $gffh_sort = $self->sort_gff($ltrh_gff);

	my @ltrd_opts = qw(-trnas -hmms -aliout -aaout -seqfile -matchdescstart -seqnamelen -o 
                           -outfileprefix -pdomevalcutoff -pdomcutoff -pptradius -pptlen -pptaprob 
                           -pptgprob -uboxlen -pptuprob -pbsalilen -pbsradius -pbsoffset -pbstrnaoffset
                           -pbsmaxedist -maxgaplen);
	my @ltrd_args = ($trnadb,$hmmdb,"no","no",$genome,"yes","50",$ltrg_gff,$ltrg_out,
			 $pdomevalue,$pdomcutoff,$pptradius,$pptlen,$pptagpr,$pptagpr,$uboxlen,
	                 $uboxutpr,$pbslen,$pbsradius,$pbsoffset,$pbstrnaoffset,$pbsmaxeditdist,$maxgaplen);
	@ltrd_cmd{@ltrd_opts} = @ltrd_args;
	
	my $ltr_dig = $self->run_ltrdigest(\%ltrd_cmd, $gffh_sort);
	unlink $ltrh_gff;
	unlink $gffh_sort;
    
	return $ltrg_gff;
    }
    else {
	unlink $ltrh_gff;
	return undef;
    }
}

sub ltr_search_relaxed {
    my $self = shift;
    my ($config, $index) = @_;
    
    my $genome = $self->genome;
    my $hmmdb  = $self->hmmdb;
    my $trnadb = $self->trnadb;
    #my $index  = $self->index;

    ## LTR constraints                                                                                         
    my $overlaps   = $config->{overlaps};
    my $mintsd     = $config->{mintsd};
    my $maxtsd     = $config->{maxtsd};
    my $minlenltr  = $config->{minlenltr};
    my $maxlenltr  = $config->{maxlenltr};
    my $mindistltr = $config->{mindistltr};
    my $maxdistltr = $config->{maxdistltr};
    my $pdomcutoff = $config->{pdomcutoff};
    my $pdomevalue = $config->{pdomevalue};
    
    my $seedlength = $config->{seedlength};
    my $tsdradius  = $config->{tsdradius};
    my $xdrop      = $config->{xdrop};
    my $swmat      = $config->{swmat};
    my $swmis      = $config->{swmis};
    my $swins      = $config->{swins};
    my $swdel      = $config->{swdel};

    my $pptradius = $config->{pptradius};
    my $pptlen    = $config->{pptlen};
    my $pptagpr   = $config->{pptagpr};
    my $uboxlen   = $config->{uboxlen};
    my $uboxutpr  = $config->{uboxutpr};
    my $pbsradius = $config->{pbsradius};
    my $pbslen    = $config->{pbslen};
    my $pbsoffset = $config->{pbsoffset};
    my $pbstrnaoffset  = $config->{pbstrnaoffset};
    my $pbsmaxeditdist = $config->{pbsmaxeditdist};
    my $maxgaplen = $config->{maxgaplen};
    #my $genome = $self->genome->absolute;
    #my $hmmdb  = $self->hmmdb;
    #my $trnadb = $self->trnadb;
    my $gt = $self->get_gt_exec;

    my (%suf_args, %ltrh_cmd, %ltrd_cmd);
    
    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
	$name =~ s/$1//;
    }
    
    my $ltrh_out = File::Spec->catfile($path, $name."_ltrharvest85_pred-all");
    my $ltrh_out_inner = File::Spec->catfile($path, $name."_ltrharvest85_pred-inner");
    my $ltrh_gff = File::Spec->catfile($path, $name."_ltrharvest85.gff3");
    my $ltrg_gff = File::Spec->catfile($path, $name."_ltrdigest85.gff3");
    my $ltrg_out = File::Spec->catfile($path, $name."_ltrdigest85");

    my @ltrh_opts = qw(-longoutput -seqids -tabout -mintsd -maxtsd -minlenltr -maxlenltr 
                       -mindistltr -maxdistltr -similar -vic -index -out -outinner -overlaps
                       -seed -vic -xdrop -mat -mis -ins -del -gff3);

    my @ltrh_args = ("no","yes","no",$mintsd,$maxtsd,$minlenltr,$maxlenltr,$mindistltr,$maxdistltr,"85","10",
		     $index,$ltrh_out,$ltrh_out_inner,$overlaps,$seedlength,$tsdradius,$xdrop,$swmat,
                     $swmis,$swins,$swdel,$ltrh_gff);
    @ltrh_cmd{@ltrh_opts} = @ltrh_args;
    
    my $ltr_succ  = $self->run_ltrharvest(\%ltrh_cmd);
    my $gffh_sort = $self->sort_gff($ltrh_gff);

    my @ltrd_opts = qw(-trnas -hmms -aliout -aaout -seqfile -matchdescstart -seqnamelen -o 
                       -outfileprefix -pdomevalcutoff -pdomcutoff -pptradius -pptlen -pptaprob
                       -pptgprob -uboxlen -pptuprob -pbsalilen -pbsradius -pbsoffset 
                       -pbstrnaoffset -pbsmaxedist -maxgaplen);
    my @ltrd_args = ($trnadb,$hmmdb,"no","no",$genome,"yes","50",$ltrg_gff,$ltrg_out,$pdomevalue,
		     $pdomcutoff,$pptradius,$pptlen,$pptagpr,$pptagpr,$uboxlen,$uboxutpr,$pbslen,
		     $pbsradius,$pbsoffset,$pbstrnaoffset,$pbsmaxeditdist,$maxgaplen);
    @ltrd_cmd{@ltrd_opts} = @ltrd_args;
    
    my $ltr_dig = $self->run_ltrdigest(\%ltrd_cmd, $gffh_sort);
    $self->clean_index if $self->clean;
    unlink $ltrh_gff;
    unlink $gffh_sort;

    return $ltrg_gff;
}

sub get_configuration {
    my $self = shift;
    my $configfile   = YAML::Tiny->read( $self->config );
    my $valid_config = $self->parse_configuration( $configfile );
    return $valid_config;
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

    perldoc Tephra::LTR::LTRSearch


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
