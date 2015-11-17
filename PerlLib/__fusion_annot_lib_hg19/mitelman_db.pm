package mitelman_db;

use strict;
use warnings;
use Carp;

sub load_data {
    my ($annotations_href, $mitelman_db_file) = @_;

    print STDERR "-parsing mitelman db: $mitelman_db_file\n";
    
    open (my $fh, "$mitelman_db_file") or confess "Error, cannot read file: $mitelman_db_file";
    while (<$fh>) {
        chomp;
        my ($fusion, $tumor_types) = split(/\t/);
        $annotations_href->{$fusion}->{"Mitelman{$tumor_types}"} = 1;
    }
    close $fh;
    
    return;
}


1; #EOM


    
