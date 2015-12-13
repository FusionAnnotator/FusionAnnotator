package Klijn_fusions;

use strict;
use warnings;
use Carp;


# fusions reported in Klijn et al. http://www.ncbi.nlm.nih.gov/pubmed/25485619 
# A comprehensive transcriptional portrait of human cancer cell lines.
# Nat Biotechnol. 2015 Mar;33(3):306-12. doi: 10.1038/nbt.3080. Epub 2014 Dec 8.

sub load_data {
    my ($annotations_href, $klijn_fusions) = @_;

    print STDERR "-parsing Klijn et al. cancer cell line fusions: $klijn_fusions\n";
    
    open (my $fh, $klijn_fusions) or confess "Error, cannot open file $klijn_fusions";
   
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
        
        my $annot = "{Klijn_CCL:" . join(",", @annots) . "}";
        $annot =~ s/\s+/_/g;
        
        $annotations_href->{$fusion}->{$annot} = 1;
        
    }
    
    return;
}


1; #EOM


    
