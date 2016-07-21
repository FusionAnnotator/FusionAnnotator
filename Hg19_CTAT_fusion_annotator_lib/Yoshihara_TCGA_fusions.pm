package Yoshihara_TCGA_fusions;

use strict;
use warnings;
use Carp;

sub load_data {
    my ($annotations_href, $tcga_fusions_file) = @_;
    
    print STDERR "-parsing TCGA recurrent fusions from Yoshihara Oncogene 2014\n";
    
    open (my $fh, $tcga_fusions_file) or confess "Error, cannot open file $tcga_fusions_file";
   
    my $topline = <$fh>;
    my $header = <$fh>;
    chomp $header;
    my @tumor_type = split(/\t/, $header);
    shift @tumor_type;
    
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $fusion = shift @x;
        $fusion =~ s/__/--/;
        
        my @fusion_annots;
        for (my $i = 0; $i <= 12; $i++) {
            my $num_samples = $x[$i];
            if ($num_samples) {
                my $tumor = $tumor_type[$i];
                push (@fusion_annots, "$tumor:$num_samples");
            }
        }
        
        my $fusion_annot = "YOSHIHARA_TCGA_num_samples[" . join("|", @fusion_annots) . "]";
        
        $annotations_href->{$fusion}->{$fusion_annot} = 1;
        
    }
    close $fh;
    
    
    return;
}


1; #EOM


    
