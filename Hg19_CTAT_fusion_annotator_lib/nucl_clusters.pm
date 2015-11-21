package nucl_clusters;

use strict;
use warnings;
use Carp;

sub load_data {
    my ($annotations_href, $nucl_clusters_dat_file) = @_;

    print STDERR "-parsing nucl clusters\n";
    
    open (my $fh, "$nucl_clusters_dat_file") or confess "Error, cannot read file: $nucl_clusters_dat_file";
    while (<$fh>) {
        chomp;
        my @genes = split(/\s+/);
        for (my $i = 0; $i < $#genes; $i++) {
            my $gene_i = $genes[$i];
            
            for (my $j = $i + 1; $j <= $#genes; $j++) {
                my $gene_j = $genes[$j];
                
                $annotations_href->{"$gene_i--$gene_j"}->{"NUCLCLUST"} = 1;
                $annotations_href->{"$gene_j--$gene_i"}->{"NUCLCLUST"} = 1;
                
            }
        }
    }
    
    
    return;
}


1; #EOM


    
