package Stransky2015_normals;

use strict;
use warnings;
use Carp;

sub load_data {
    my ($annotations_href, $gtex_normals) = @_;

    print STDERR "-parsing $gtex_normals\n";
    
    open (my $fh, $gtex_normals) or confess "Error, cannot open file $gtex_normals";
    
    my $header = <$fh>;
    while (<$fh>) {
        chomp;
        my ($geneA, $geneB, $count, $tissues, $sample_count) = split(/\t/);
                
        my $fusion = join("--", $geneA, $geneB);
        
        my $pct_normals = sprintf("$count/$sample_count=%.2f%%", $count/$sample_count*100);
        
        $annotations_href->{$fusion}->{"GTEx\&tcga_normals:$pct_normals"} = 1;
    }
    close $fh;
    
    return;
}


1; #EOM


    
