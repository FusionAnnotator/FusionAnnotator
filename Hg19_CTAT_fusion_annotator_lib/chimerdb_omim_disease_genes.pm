package chimerdb_omim_disease_genes;

use strict;
use warnings;
use Carp;

sub load_data {
    my ($annotations_href, $chimerdb_omim_txt) = @_;

    print STDERR "-parsing chimerdb omim: $chimerdb_omim_txt\n";
    
    open (my $fh, "$chimerdb_omim_txt") or confess "Error, cannot read file: $chimerdb_omim_txt";
    my $header = <$fh>;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $geneA = $x[2];
        my $geneB = $x[4];
        my $disease = $x[5];
        unless ($disease && $disease =~ /\w/) { next; }
        $disease =~ s/\"//g;
        
        my $fusion = "$geneA--$geneB";
        $annotations_href->{$fusion}->{"chimerdb_omim:{$disease}"} = 1;
    }
    close $fh;
    
    return;
}


1; #EOM


    
