package blast_pairs;

use strict;
use warnings;
use Carp;

sub load_data {
    my ($annotations_href, $blast_pairs_gz_filename) = @_;
    
    print STDERR "-parsing blast pairs: $blast_pairs_gz_filename\n";
    
    open (my $fh, "gunzip -c $blast_pairs_gz_filename | ") or confess "Error, cannot read file via: gunzip -c $blast_pairs_gz_filename";
    while (<$fh>) {
        chomp;
        my ($geneA, $geneB, $per_id, $Evalue) = split(/\s+/);
        if ($geneA eq $geneB) { next; }
        my $token = "{$geneA|$geneB|pID:$per_id|E:$Evalue}";
        $annotations_href->{"$geneA--$geneB"}->{"BLASTMATCH:$token"} = 1;
        $annotations_href->{"$geneB--$geneA"}->{"BLASTMATCH:$token"} = 1;
    }
    close $fh;
    
    return;
}


1; #EOM


    
