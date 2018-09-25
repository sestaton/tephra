package Tephra::Classify::Role::LogResults;

use 5.014;
use Moose::Role;
use Statistics::Descriptive;
use namespace::autoclean;
#use Data::Dump::Color;

=head1 NAME

Tephra::Classify::Role::LogResults - Reusable methods for calculating and logging basic classification statistics

=head1 VERSION

Version 0.12.2

=cut

our $VERSION = '0.12.2';
$VERSION = eval $VERSION;

#
# methods
#
sub log_basic_element_stats {
    my $self = shift;
    my ($object) = @_;

    my ($lengths, $type, $log, $pdoms) = @{$object}{qw(lengths type log pdom_ct)};

    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@$lengths);
    my $min   = $stat->min // 0;
    my $max   = $stat->max // 0;
    my $mean  = defined $stat->mean ? sprintf("%.2f", $stat->mean) : 0;
    my $count = $stat->count;

    my $tot_str = sprintf("%-70s %-10s", "Results - Total number of $type elements:", $count);
    my $min_str = sprintf("%-70s %-10s", "Results - Minimum length of $type elements:", $min);
    my $max_str = sprintf("%-70s %-10s", "Results - Maximum length of $type elements:", $max);
    my $ave_str = sprintf("%-70s %-10s", "Results - Mean length of $type elements:", $mean);
    my $pct_str = sprintf("%-70s %-10s", "Results - Number of $type elements with protein matches:", $pdoms);

    $log->info($tot_str);
    $log->info($min_str);
    $log->info($max_str);
    $log->info($ave_str);
    $log->info($pct_str);
    
    return $count;
}

sub write_superfam_pdom_organization {
    my $self = shift;
    my ($dom_obj) = @_;
    
    my ($pdom_index, $domoutfile) = @{$dom_obj}{qw(pdom_index outfile)};
    open my $domf, '>>', $domoutfile or die "\n[ERROR]: Could not open file: $domoutfile\n";

    my %tot_dom_ct;
    say $domf join "\t", "Strand", "Domain_organizaion", "Domain_count";

    for my $strand (keys %$pdom_index) {
	for my $org (reverse sort { $pdom_index->{$strand}{$a} <=> $pdom_index->{$strand}{$b} } keys %{$pdom_index->{$strand}}) {
	    $tot_dom_ct{$org} += $pdom_index->{$strand}{$org};
	    say $domf join "\t", $strand, $org, $pdom_index->{$strand}{$org};
	}
    }
    
    say $domf "#========== Below are the domain architectures summarized for both strands ==========";
    say $domf join "\t", "Domain_organization", "Domain_count";
    for my $domorg (reverse sort { $tot_dom_ct{$a} <=> $tot_dom_ct{$b} } keys %tot_dom_ct) {
	say $domf join "\t", $domorg, $tot_dom_ct{$domorg};
    }
    close $domf;

    return;
}

sub write_fam_pdom_organization {
    my $self = shift;
    my ($dom_obj) = @_;
    
    my ($pdom_famid_map, $outfh, $elemnum, $elemid, $famname) = @{$dom_obj}{qw(famid_map outfh elemnum elemid famname)};
    my $re = qr/helitron\d+|non_LTR_retrotransposon\d+|(?:LTR|TRIM|LARD)_retrotransposon\d+|terminal_inverted_repeat_element\d+|MITE\d+/;
    my ($element) = ($elemnum =~ /($re)/);
 
    say $outfh join "\t", $famname, $elemid, $pdom_famid_map->{$element}{pdoms};

    return;
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

    perldoc Tephra::Classify::Role::LogResults


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut 

1;
