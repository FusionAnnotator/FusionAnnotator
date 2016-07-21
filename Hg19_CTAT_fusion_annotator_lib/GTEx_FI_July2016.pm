package GTEx_FI_July2016;

## bhaas - ran starF / FI on all of GTEx

use strict;
use warnings;
use Carp;

sub load_data {
    my ($annotations_href, $gtex_normals) = @_;

    print STDERR "-parsing GTEx_FI.July2016 $gtex_normals\n";
    
    open (my $fh, $gtex_normals) or confess "Error, cannot open file $gtex_normals";

    while (<$fh>) {
        chomp;
        my ($fusion, $GTEx_FI_annot) = split(/\t/);
            
        $annotations_href->{$fusion}->{"GTEx_FI:{$GTEx_FI_annot}"} = 1;
    }
    close $fh;
    
    return;
}


1; #EOM


    
