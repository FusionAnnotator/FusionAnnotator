package cosmic_fusions;

use strict;
use warnings;
use Carp;

sub load_data {
    my ($annotations_href, $cosmic_fusions_file) = @_;

    print STDERR "-parsing cosmic_fusions: $cosmic_fusions_file\n";
    
    open (my $fh, "$cosmic_fusions_file") or confess "Error, cannot read file: $cosmic_fusions_file";
    while (<$fh>) {
        chomp;
        if (/^\#/) { next; }
        my ($fusion, $samples, $mutations, $papers) = split(/\s+/);
        $fusion =~ s|/|--|;
        $annotations_href->{$fusion}->{"Cosmic{samples=$samples,mutations=$mutations,papers=$papers}"} = 1;
        
    }
    close $fh;
    
    return;
}


1; #EOM


    
