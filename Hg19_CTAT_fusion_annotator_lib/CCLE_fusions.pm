package CCLE_fusions;

use strict;
use warnings;
use Carp;

sub load_data {
    my ($annotations_href, $ccle_fusions) = @_;

    print STDERR "-parsing CCLE_fusions: $ccle_fusions\n";
    
    open (my $fh, $ccle_fusions) or confess "Error, cannot open file $ccle_fusions";
   
    my %fusion_to_tissue_counts;
   
    while (<$fh>) {
        chomp;
        if (/^\#/) { next; }
        unless (/\w/) { next; }
        my @x = split(/\t/);
        unless (scalar @x > 17) { next; }
        my $fusion = join("--", $x[1], $x[2]);
        my $tissue_type = $x[17];
        
        $fusion_to_tissue_counts{$fusion}->{$tissue_type}++;
    }
    close $fh;
    
    foreach my $fusion (keys %fusion_to_tissue_counts) {
        
        my $tissues_href = $fusion_to_tissue_counts{$fusion};

        my @tissues = reverse sort {$tissues_href->{$a} <=> $tissues_href->{$b}} keys %$tissues_href;

        my @annots;
        foreach my $tissue (@tissues) {
            my $count = $tissues_href->{$tissue};
            push (@annots, "$tissue=$count");
        }
        
        my $annot = "{CCLE:" . join(",", @annots) . "}";
        $annot =~ s/\s+/_/g;
        
        $annotations_href->{$fusion}->{$annot} = 1;
        
    }
    
    return;
}


1; #EOM


    
