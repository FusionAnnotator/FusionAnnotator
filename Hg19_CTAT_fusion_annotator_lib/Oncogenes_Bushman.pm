package Oncogenes_Bushman;

use strict;
use warnings;
use Carp;


# oncogenes list from: http://www.bushmanlab.org/links/genelists

sub load_data {
    my ($annotations_href, $oncogenes_file) = @_;

    print STDERR "-parsing Oncogenes listing: $oncogenes_file\n";
    
    open (my $fh, $oncogenes_file) or confess "Error, cannot open file $oncogenes_file";
    my $header = <$fh>;

    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $gene_id = $x[2];
                
        $annotations_href->{$gene_id}->{"Oncogene{Bushman}"} = 1;    
        
    }
    close $fh;
    
        
    return;
}


1; #EOM


    
