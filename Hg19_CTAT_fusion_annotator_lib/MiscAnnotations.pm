package MiscAnnotations;

use strict;
use warnings;
use Carp;

sub load_data {
    my ($annotations_href, $annots_file) = @_;

    print STDERR "-parsing misc annotations: $annots_file\n";
    
    open (my $fh, $annots_file) or confess "Error, cannot open file $annots_file";
    my $header = <$fh>;

    while (<$fh>) {
        chomp;
        if (/^\#/) { next; }
        my ($target, $annot) = split(/\s+/, $_, 2);
        
                
        $annotations_href->{$target}->{"$annot"} = 1;    
        
    }
    close $fh;
    
        
    return;
}


1; #EOM


    
