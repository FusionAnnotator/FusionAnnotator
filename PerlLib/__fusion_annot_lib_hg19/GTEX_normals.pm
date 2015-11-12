package GTEX_normals;

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
        my ($fusion, $num_normals, $pct_normals) = split(/\t/);
            
        $annotations_href->{$fusion}->{"GTEx:$pct_normals%"} = 1;
    }
    close $fh;
        
    return;
}


1; #EOM


    
