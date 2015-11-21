package gene_spans;

use strict;
use warnings;
use Carp;

sub load_data {
    my ($annotations_href, $gene_spans_file) = @_;
    
    print STDERR "- parsing gene chr locations\n";

    open (my $fh, $gene_spans_file) or confess "Error, cannot open file $gene_spans_file";
        
    while (<$fh>) {
        chomp;
            
        my ($gene_id, $chr, $lend, $rend, $orient, $gene_name) = split(/\t/);
        if ($gene_name && $gene_name ne ".") {
            $gene_id = $gene_name;
        }
            
        $annotations_href->{$gene_id}->{"$chr:$lend-$rend:$orient"} = 1;
    }
    close $fh;
    
    
    return;
}


1; #EOM


    
