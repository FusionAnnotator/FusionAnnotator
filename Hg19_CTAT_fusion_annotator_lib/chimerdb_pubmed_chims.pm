package chimerdb_pubmed_chims;

use strict;
use warnings;
use Carp;

sub load_data {
    my ($annotations_href, $chimerdb_pubmed_chims) = @_;

    print STDERR "-parsing chimerdb_pubmed_chims: $chimerdb_pubmed_chims\n";

    my %fusion_to_annotations;
    
    open (my $fh, "$chimerdb_pubmed_chims") or confess "Error, cannot read file: $chimerdb_pubmed_chims";
    my $header = <$fh>;
    while (<$fh>) {
        chomp;
        my ($geneA, $geneB, $disease_name, $pubmed) = split(/\t/);
        my $fusion = join("--", $geneA, $geneB);
        if ($disease_name && $disease_name =~ /\w/) {
            $fusion_to_annotations{$fusion}->{$disease_name}++;
        }
    }
    close $fh;
    
    foreach my $fusion (keys %fusion_to_annotations) {
        my @annotations = keys %{$fusion_to_annotations{$fusion}};
        
        my $annot_list = join(",", @annotations);

        $annotations_href->{$fusion}->{"chimerdb_pubmed{$annot_list}"} = 1;
    }
    
    return;
}


1; #EOM


    
