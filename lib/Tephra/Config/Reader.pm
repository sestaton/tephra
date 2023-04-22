package Tephra::Config::Reader;

use 5.014;
use Moose;
use YAML::Tiny;
use File::Spec;
use File::Basename;
use File::Temp qw(tempfile);
#use Data::Dump::Color;
use namespace::autoclean;

with 'Tephra::Role::File';

=head1 NAME

Tephra::Config::Reader - Attributes and routines for parsing Tephra configuration. 

=head1 VERSION

Version 0.14.0

=cut

our $VERSION = '0.14.0';
$VERSION = eval $VERSION;

=head1 SYNOPSIS

    use Tephra::Config::Reader;

    ...

=cut

has config => (
    is            => 'ro',
    isa           => 'Str',
    required      => 0,
    documentation => qq{The Tephra configuration file},
);

=head1 METHODS

=head2 get_configuration

=cut

sub get_configuration {
    my $self = shift;
    my $configfile   = YAML::Tiny->read( $self->config );
    my $valid_config = $self->parse_configuration( $configfile );
    return $valid_config;
}

=head2 parse_configuration

 Title    : parse_config

 Usage    : my $config = $trans_obj->parse_configuration;

 Function : The parsed configuration for Tephra options.

                                                           Return_type
 Returns  : A hash containing the user-set configuration   HashRef
            for how transposome is to be executed

 Args     : None. This is a class method.

=cut 

sub parse_configuration {
    my $self = shift;
    my ($yaml) = @_;

    my %config;
    my $index = 0;

    # all options from config
    $config{all}{logfile}   = $yaml->[0]{all}[$index]{logfile};
    $index++;
    $config{all}{genome}    = $yaml->[0]{all}[$index]{genome};
    $index++;
    $config{all}{outfile}   = $yaml->[0]{all}[$index]{outfile};
    $index++;
    $config{all}{repeatdb}  = $yaml->[0]{all}[$index]{repeatdb};
    $index++;
    $config{all}{genefile}  = $yaml->[0]{all}[$index]{genefile};
    $index++;
    $config{all}{trnadb}    = $yaml->[0]{all}[$index]{trnadb};
    $index++;
    $config{all}{hmmdb}     = $yaml->[0]{all}[$index]{hmmdb};
    $index++;
    $config{all}{threads}   = $yaml->[0]{all}[$index]{threads};
    $index++;
    $config{all}{clean}     = $yaml->[0]{all}[$index]{clean};
    $index++;
    $config{all}{debug}     = $yaml->[0]{all}[$index]{debug};
    $index++;
    $config{all}{subs_rate} = $yaml->[0]{all}[$index]{subs_rate};

    # findltrs
    # ltrharvest options from config
    my $ltrh_index = 0;
    $index = 0;
    $config{findltrs}{dedup}      = $yaml->[0]{findltrs}[$index]{dedup};
    $index++;
    $config{findltrs}{tnpfilter}  = $yaml->[0]{findltrs}[$index]{tnpfilter};
    $index++;
    $config{findltrs}{domains_required} = $yaml->[0]{findltrs}[$index]{domains_required};
    $index++;
    $config{findltrs}{mintsd}     = $yaml->[0]{findltrs}[$index]{ltrharvest}[$ltrh_index]{mintsd};
    $ltrh_index++;
    $config{findltrs}{maxtsd}     = $yaml->[0]{findltrs}[$index]{ltrharvest}[$ltrh_index]{maxtsd};
    $ltrh_index++;
    $config{findltrs}{minlenltr}  = $yaml->[0]{findltrs}[$index]{ltrharvest}[$ltrh_index]{minlenltr};
    $ltrh_index++;
    $config{findltrs}{maxlenltr}  = $yaml->[0]{findltrs}[$index]{ltrharvest}[$ltrh_index]{maxlenltr};
    $ltrh_index++;
    $config{findltrs}{mindistltr} = $yaml->[0]{findltrs}[$index]{ltrharvest}[$ltrh_index]{mindistltr};
    $ltrh_index++;
    $config{findltrs}{maxdistltr} = $yaml->[0]{findltrs}[$index]{ltrharvest}[$ltrh_index]{maxdistltr};
    $ltrh_index++;
    $config{findltrs}{seedlength} = $yaml->[0]{findltrs}[$index]{ltrharvest}[$ltrh_index]{seedlength};
    $ltrh_index++;
    $config{findltrs}{tsdradius}  = $yaml->[0]{findltrs}[$index]{ltrharvest}[$ltrh_index]{tsdradius};
    $ltrh_index++;
    $config{findltrs}{xdrop}      = $yaml->[0]{findltrs}[$index]{ltrharvest}[$ltrh_index]{xdrop};
    $ltrh_index++;
    $config{findltrs}{swmat}      = $yaml->[0]{findltrs}[$index]{ltrharvest}[$ltrh_index]{swmat};
    $ltrh_index++;
    $config{findltrs}{swmis}      = $yaml->[0]{findltrs}[$index]{ltrharvest}[$ltrh_index]{swmis};
    $ltrh_index++;
    $config{findltrs}{swins}      = $yaml->[0]{findltrs}[$index]{ltrharvest}[$ltrh_index]{swins};
    $ltrh_index++;
    $config{findltrs}{swdel}      = $yaml->[0]{findltrs}[$index]{ltrharvest}[$ltrh_index]{swdel};
    $ltrh_index++;
    $config{findltrs}{overlaps}   = $yaml->[0]{findltrs}[$index]{ltrharvest}[$ltrh_index]{overlaps};
    $index++;

    # ltrdigest options from config
    my $ltrd_index = 0;
    $config{findltrs}{pptradius}      = $yaml->[0]{findltrs}[$index]{ltrdigest}[$ltrd_index]{pptradius};
    $ltrd_index++;
    $config{findltrs}{pptlen}         = $yaml->[0]{findltrs}[$index]{ltrdigest}[$ltrd_index]{pptlen};
    $ltrd_index++;
    $config{findltrs}{pptagpr}        = $yaml->[0]{findltrs}[$index]{ltrdigest}[$ltrd_index]{pptagpr};
    $ltrd_index++;
    $config{findltrs}{uboxlen}        = $yaml->[0]{findltrs}[$index]{ltrdigest}[$ltrd_index]{uboxlen};
    $ltrd_index++;
    $config{findltrs}{uboxutpr}       = $yaml->[0]{findltrs}[$index]{ltrdigest}[$ltrd_index]{uboxutpr};
    $ltrd_index++;
    $config{findltrs}{pbsradius}      = $yaml->[0]{findltrs}[$index]{ltrdigest}[$ltrd_index]{pbsradius};
    $ltrd_index++;
    $config{findltrs}{pbslen}         = $yaml->[0]{findltrs}[$index]{ltrdigest}[$ltrd_index]{pbslen};
    $ltrd_index++;
    $config{findltrs}{pbsoffset}      = $yaml->[0]{findltrs}[$index]{ltrdigest}[$ltrd_index]{pbsoffset};
    $ltrd_index++;
    $config{findltrs}{pbstrnaoffset}  = $yaml->[0]{findltrs}[$index]{ltrdigest}[$ltrd_index]{pbstrnaoffset};
    $ltrd_index++;
    $config{findltrs}{pbsmaxeditdist} = $yaml->[0]{findltrs}[$index]{ltrdigest}[$ltrd_index]{pbsmaxeditdist};
    $ltrd_index++;
    $config{findltrs}{pdomevalue}     = $yaml->[0]{findltrs}[$index]{ltrdigest}[$ltrd_index]{pdomevalue};
    $ltrd_index++;
    $config{findltrs}{pdomcutoff}     = $yaml->[0]{findltrs}[$index]{ltrdigest}[$ltrd_index]{pdomcutoff};
    $ltrd_index++;
    $config{findltrs}{maxgaplen}      = $yaml->[0]{findltrs}[$index]{ltrdigest}[$ltrd_index]{maxgaplen};

    # classifyltrs
    $index = 0;
    $config{classifyltrs}{percentcov} = $yaml->[0]{classifyltrs}[$index]{percentcov};
    $index++;
    $config{classifyltrs}{percentid}  = $yaml->[0]{classifyltrs}[$index]{percentid};
    $index++;
    $config{classifyltrs}{hitlen}     = $yaml->[0]{classifyltrs}[$index]{hitlen};

    # illrecomb
    $index = 0;
    $config{illrecomb}{repeat_pid} = $yaml->[0]{illrecomb}[$index]{repeat_pid};

    # ltrage
    $config{ltrage}{all} = $yaml->[0]{ltrage}[$index]{all};

    # maskref
    $config{maskref}{percentid} = $yaml->[0]{maskref}[$index]{percentid};
    $index++;
    $config{maskref}{hitlength} = $yaml->[0]{maskref}[$index]{hitlength};
    $index++;
    $config{maskref}{splitsize} = $yaml->[0]{maskref}[$index]{splitsize};
    $index++;
    $config{maskref}{overlap}   = $yaml->[0]{maskref}[$index]{overlap};

    # sololtr
    $index = 0;
    $config{sololtr}{percentid}   = $yaml->[0]{sololtr}[$index]{percentid};
    $index++;
    $config{sololtr}{percentcov}  = $yaml->[0]{sololtr}[$index]{percentcov};
    $index++;
    $config{sololtr}{matchlen}    = $yaml->[0]{sololtr}[$index]{matchlen};
    $index++;
    $config{sololtr}{numfamilies} = $yaml->[0]{sololtr}[$index]{numfamilies};
    $index++;
    $config{sololtr}{allfamilies} = $yaml->[0]{sololtr}[$index]{allfamilies};

    # tirage
    $index = 0;
    $config{tirage}{all} = $yaml->[0]{tirage}[$index]{all};
    
    #dd \%config and exit;
    my $valid_config = $self->_validate_params(\%config);

    return $valid_config;
}

=head2 _validate_params

 Title    : _validate_params

 Usage    : This is a private method, do not use it directly.

 Function : Valiadate whether all of the slots in config file
            have been set.

                                                           Return_type
 Returns  : A hash containing the user-set configuration   HashRef
            for how transposome is to be executed

 Args     : None. This is a class method.

=cut 

sub _validate_params {
    my $self = shift;
    my ($config) = @_;
    
    for my $cmd (keys %$config) {
	for my $opt (keys %{$config->{$cmd}}) {
	    my $v = $config->{$cmd}{$opt};
	    if (defined $v && $v =~ /^~/) {
                $v =~ s/^~/$ENV{"HOME"}/;
                $config->{$cmd}{$opt} = $v;
            }
	    elsif (not defined $v) {
	    #if ($cmd ne 'all' && ! defined $v) {
		#elsif (not defined $v && $cmd ne 'all') {
		die "[ERROR]: '$opt' under '$cmd' is not defined after parsing configuration file.\n".
		    "         This indicates there may be a blank line in your configuration file.\n".
		    "         Please check your configuration file and try again. Exiting.\n";
	    }
	}
    }
    
    return $config;
}

sub get_all_opts {
    my $self = shift;
    my ($config) = @_;

    my ($logfile, $genome, $repeatdb, $genefile, $hmmdb, $trnadb, $outfile, $clean, $debug, 
	$threads, $subs_rate, $genome_compressed, $repeatdb_compressed, $genefile_compressed);

    if (defined $config->{all}{genome} && -e $config->{all}{genome}) {
        ($genome, $genome_compressed) = 
	    $self->check_if_compressed({ infile => $config->{all}{genome}, is_genome => 1, is_repeatdb => 0, is_genefile => 0 });
    }
    else {
        say STDERR "\n[ERROR]: genome file was not defined in configuration or does not exist. Check input. Exiting.\n";
        exit(1);
    }

    if (defined $config->{all}{repeatdb} && -e $config->{all}{repeatdb}) {
        ($repeatdb, $repeatdb_compressed) = 
	    $self->check_if_compressed({ infile => $config->{all}{repeatdb}, is_genome => 0, is_repeatdb => 1, is_genefile => 0 });
    }   
    else {
        say STDERR "\n[ERROR]: repeatdb file was not defined in configuration or does not exist. Check input. Exiting.\n";
        exit(1);
    }

    if (defined $config->{all}{genefile} && -e $config->{all}{genefile}) {
        ($genefile, $genefile_compressed) = 
	    $self->check_if_compressed({ infile => $config->{all}{genefile}, is_genome => 0, is_repeatdb => 0, is_genefile => 1 });
    }
    else {
        say STDERR "\n[ERROR]: gene file was not defined in configuration or does not exist. Check input. Exiting.\n";
        exit(1);
    }

    my ($name, $path, $suffix); # genome file specs
    if (defined $config->{all}{outfile}) {
        $outfile = $config->{all}{outfile};
    }
    else {
        ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
        $outfile = File::Spec->catfile( abs_path($path), $name.'_tephra_transposons.gff3' );
        $logfile = $config->{all}{logfile} // File::Spec->catfile( abs_path($path), $name.'_tephra_full.log' );
    }

    my $execonfig = Tephra::Config::Exe->new->get_config_paths;
    my ($tephra_hmmdb, $tephra_trnadb) = @{$execonfig}{qw(hmmdb trnadb)};
    
    $hmmdb = defined $config->{all}{hmmdb}  && $config->{all}{hmmdb} =~ /tephradb/i ? 
        $tephra_hmmdb : $config->{all}{hmmdb};
    $trnadb = defined $config->{all}{trnadb} && $config->{all}{trnadb} =~ /tephradb/i ? 
        $tephra_trnadb : $config->{all}{trnadb};

    $clean = defined $config->{all}{clean} && $config->{all}{clean} =~ /yes/i ? 1 : 0;
    $debug = defined $config->{all}{debug} && $config->{all}{debug} =~ /yes/i ? 1 : 0;
    $threads = $config->{all}{threads} // 1;
    $subs_rate = $config->{all}{subs_rate} // 1e-8;
    $logfile = $config->{all}{logfile} // File::Spec->catfile( abs_path($path), $name.'_tephra_full.log' );

    return { logfile   => $logfile,
             genome    => $genome, 
             repeatdb  => $repeatdb,
	     genefile  => $genefile,
             hmmdb     => $hmmdb, 
             trnadb    => $trnadb, 
             outfile   => $outfile, 
             clean     => $clean, 
             debug     => $debug, 
             threads   => $threads, 
             subs_rate => $subs_rate, 
             genome_is_compressed   => $genome_compressed,
             repeatdb_is_compressed => $repeatdb_compressed,
             genefile_is_compressed => $genefile_compressed };
}

sub check_if_compressed {
    my $self = shift;
    my ($struc) = @_;
    my ($infile, $is_genome, $is_repeatdb, $is_genefile) =
	@{$struc}{qw(infile is_genome is_repeatdb is_genefile)};

    my $is_compressed = 0;
    if ($infile =~ /\.gz\z|\.bz2\z/) {
	my $fh = $self->get_fh($infile);

	my $tmpfname;
	my ($name, $path, $suffix) = fileparse($infile, qr/\.[^.]*/);
	$name =~ s/\.fa.*$//;
	$tmpfname = $name.'_tephra_tmp_genome_XXXX' if $is_genome;
	$tmpfname = $name.'_tephra_tmp_repeatdb_XXXX' if $is_repeatdb;
	$tmpfname = $name.'_tephra_tmp_genes_XXXX' if $is_genefile;

	my ($outf, $ffilename) = tempfile( TEMPLATE => $tmpfname, DIR => $path, UNLINK => 0, SUFFIX => '.fasta' );
	while (my $line = <$fh>) {
	    chomp $line;
	    say $outf $line;
	}
	close $outf;
	close $fh;

	$infile = $ffilename;
	$is_compressed = 1;
    }

    return ($infile, $is_compressed);
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

    perldoc Tephra::LTR::Role::Config


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015 S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1; 

