package Tephra::Hel::HelSearch;

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
use namespace::autoclean;

with 'Tephra::Role::Run::HelitronScanner';

sub find_helitrons {
    my $self = shift;
    
    my $genome    = $self->genome->absolute;
    my $hscan_dir = $self->helitronscanner_dir;
    my $gff       = $self->outfile;

    my (%scanh_cmd, %scant_cmd, %pair_cmd, %draw_cmd);
    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
	$name =~ s/$1//;
    }
 
    my $g_headlcvs = File::Spec->catfile($path, $name."_hscan_head.lcvs");
    my $g_taillcvs = File::Spec->catfile($path, $name."_hscan_tail.lcvs");
    my $g_paired   = File::Spec->catfile($path, $name."_hscan_paired.txt");
    my $g_helname  = File::Spec->catfile($path, $name."_tephra_hscan_helitrons");
    
    my $jar  = File::Spec->catfile($hscan_dir, "HelitronScanner.jar");
    my $lcvs = File::Spec->catfile($hscan_dir, "TrainingSet", "head.lcvs");
    my $rcvs = File::Spec->catfile($hscan_dir, "TrainingSet", "tail.lcvs");

    #java -jar $jar scanHead -g $db -lf $lcvs -o $g_headlcvs --rc -tl 10 -buffer_size 1000000
    #java -jar $jar scanTail -g $db -lf $rcvs -o $g_taillcvs --rc -tl 10 -buffer_size 1000000
    #java -jar $jar pairends -hs $g_headlcvs -ts $g_taillcvs --rc -o $g_paired
    #java -jar $jar draw -p $g_paired -g $db -o bronze_helitrons_hscan --pure --flanking --ext 
    #-ext5 100 -ext3 100

    my @scanh_opts = qw(-g -lf -o -tl -buffer_size);
    my @scanh_args = ($genome, $lcvs, $g_headlcvs, '10', '1000000');
    @scanh_cmd{@scanh_opts} = @scanh_args;

    my @scant_opts = qw(-g -lf -o -tl -buffer_size);
    my @scant_args = ($genome, $rcvs, $g_taillcvs, '10', '1000000');
    @scant_cmd{@scant_opts} = @scant_args;

    my @pair_opts = qw(-hs -ts -o);
    my @pair_args = ($g_headlcvs, $g_taillcvs, $g_paired);
    @pair_cmd{@pair_opts} = @pair_args;

    my @draw_opts = qw(-p -g -o -ext5 -ext3);
    my @draw_args = ($g_paired, $genome, $g_helname, '100', '100');
    @draw_cmd{@draw_opts} = @draw_args;
    
    $self->run_hscan_headtail(\%scanh_cmd, $jar, 'scanHead');
    $self->run_hscan_headtail(\%scant_cmd, $jar, 'scanTail');
    $self->run_hscan_pair(\%pair_cmd, $jar);
    $self->run_hscan_draw(\%draw_cmd, $jar);

    #$self->make_hscan_gff;
}

__PACKAGE__->meta->make_immutable;

1;
