#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Carp::Assert;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use Data::Dumper;
use JSON::XS;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib", "$FindBin::Bin");
use TiedHash;
use Nuc_translator;
use FusionCodingEffect;


my $DEBUG = 0;

my $usage = <<__EOUSAGE__;

##########################################################################
#
#  --fusions <string>                : fusion predictions
#
#  --genome_lib_dir <string>         : CTAT genome lib
#                                      
#  optional:
#
#  --show_all                        : by default, only the single 'best' fusion is shown
#                                      prioritizizing INFRAME and longer fusion transcripts.
#
##########################################################################


__EOUSAGE__

    ;


my $help_flag;
my $fusions_file;
my $genome_lib_dir;
my $show_all_flag = 0;

&GetOptions ( 'h' => \$help_flag,
              'fusions=s' => \$fusions_file,
              'genome_lib_dir=s' => \$genome_lib_dir,
   
              'show_all' => \$show_all_flag,
    );

unless ($fusions_file && $genome_lib_dir) {
    
    die $usage;
}

my $prot_info_db = "$genome_lib_dir/ref_annot.prot_info.dbm";

unless (-s $prot_info_db) {
    confess "Error, cannot locate $prot_info_db - be sure to use the latest version of the CTAT genome lib";
}


my %priority = ( 'INFRAME' => 1,
                 'FRAMESHIFT' => 2,
                 'NA' => 3,
                 '.' => 4);



main: {

    open (my $fh, $fusions_file) or die "Error, cannot open file $fusions_file";
    my $header = <$fh>;
    
    unless ($header =~ /^\#/) {
        die "Error, fusion_file: $fusions_file has unrecognizable header: $header";
    }
    chomp $header;
   
    my %header_to_index;
    {
        my @fields = split(/\t/, $header);
        for (my $i = 0; $i <= $#fields; $i++) {
            $header_to_index{$fields[$i]} = $i;
        }
        
        # ensure we have the fields we need:
        my @fields_need = ('#FusionName', 'LeftGene', 'LeftBreakpoint', 'RightGene', 'RightBreakpoint');
        my $missing_field = 0;
        foreach my $field (@fields_need) {
            unless (exists $header_to_index{$field}) {
                print STDERR "ERROR, missing column header: $field\n";
                $missing_field = 1;
            }
        }
        if ($missing_field) {
            die "Error, at least one column header was missing";
        }
    }
    

    print $header . "\t" . join("\t", 
                                "CDS_LEFT_ID", 
                                "CDS_LEFT_RANGE",
                                "CDS_RIGHT_ID",
                                "CDS_RIGHT_RANGE",
                                "PROT_FUSION_TYPE",
                                "FUSION_MODEL",
                                "FUSION_CDS",
                                "FUSION_TRANSL",
                                "PFAM_LEFT",
                                "PFAM_RIGHT",
        ) . "\n";
        
    
    my $prot_info_db_tied_hash = new TiedHash( { 'use' => $prot_info_db } );

    my $pfam_domain_dbm = "$genome_lib_dir/pfam_domains.dbm";
    my $pfam_domain_db_tied_hash = undef;
    if (-s $pfam_domain_dbm) {
        $pfam_domain_db_tied_hash = new TiedHash( { 'use' => $pfam_domain_dbm } );
    }
    
    
    while (<$fh>) {
        chomp;        
        my $line = $_;

        my @x = split(/\t/);

        my $fusion_name = $x[ $header_to_index{'#FusionName'} ];
        my $gene_left = $x[ $header_to_index{'LeftGene'}];
        my $break_left = $x[ $header_to_index{'LeftBreakpoint'} ];
                
        my $gene_right = $x[ $header_to_index{'RightGene'} ];
        my $break_right = $x[ $header_to_index{'RightBreakpoint'} ];
        

        # remove gene symbol, just want gene ID
        $gene_left =~ s/^.*\^//;
        $gene_right =~ s/^.*\^//;
        

        my @results = &examine_fusion_coding_effect($gene_left, $break_left, $gene_right, $break_right, $prot_info_db_tied_hash, $pfam_domain_db_tied_hash);
        
        if (@results) {
            ## just take the single 'best' one, arbitrarily chosen as the one with the longest fusion sequence.
            @results = sort { 
                $priority{$a->{prot_fusion_type}} <=> $priority{$b->{prot_fusion_type}}
                ||
                length($b->{prot_fusion_seq}) <=> length($a->{prot_fusion_seq}) } @results;
            
            foreach my $result (@results) {
                print join("\t", 
                           $line,
                           $result->{cds_left_id}, $result->{cds_left_range},
                           $result->{cds_right_id}, $result->{cds_right_range},
                           $result->{prot_fusion_type},
                           $result->{fusion_coding_descr},
                           $result->{cds_fusion_seq},
                           $result->{prot_fusion_seq},
                           
                           $result->{left_domains},
                           $result->{right_domains}
                           
                    ) . "\n";
                
                unless ($show_all_flag) {
                    # only showing first entry.
                    last;
                }
            }
        }
        else {
            print "$line" . ("\t." x 10) . "\n";
        }

    }
    close $fh;
    
    
    exit(0);
}


####
sub examine_fusion_coding_effect {
    my ($gene_left, $break_left, $gene_right, $break_right, $prot_info_db_tied_hash, $pfam_domain_db_tied_hash) = @_;
    
    my $gene_left_aref = &decode_my_json($gene_left, $prot_info_db_tied_hash, 1);
    my $gene_right_aref = &decode_my_json($gene_right, $prot_info_db_tied_hash, 1);
    
    unless ($gene_left_aref && $gene_right_aref) {
        return();
    }

    if (0) {
        ## debugging
        foreach my $cds_obj (@$gene_left_aref, @$gene_right_aref) {
            my @segments = @{$cds_obj->{phased_segments}};
            print &segments_to_string(@segments) . "\n";
        }
        
        return();
    }
    my @results;

    
    foreach my $cds_left_obj (@$gene_left_aref) {
        
        my $cds_left_seq = $cds_left_obj->{cds_seq};
        my $cds_left_id = $cds_left_obj->{cds_id};
        
        foreach my $cds_right_obj (@$gene_right_aref) {
            
            my $cds_right_seq = $cds_right_obj->{cds_seq};
            my $cds_right_id = $cds_right_obj->{cds_id};
            
            my ($left_fuse_segments_aref, $right_fuse_segments_aref) = &try_fuse_cdss($cds_left_obj, $break_left, $cds_right_obj, $break_right);

            if (@$left_fuse_segments_aref && @$right_fuse_segments_aref) {

                my $initial_left_seg = $left_fuse_segments_aref->[0];
                
                ## see if compatible
                my $terminal_left_seg = $left_fuse_segments_aref->[$#$left_fuse_segments_aref];
                my $left_end_phase = $terminal_left_seg->{phase_end};
                my $left_rel_rend = $terminal_left_seg->{rel_rend};

                my $left_cds_part = substr($cds_left_seq, 0, $left_rel_rend);
                
                my $initial_right_seg = $right_fuse_segments_aref->[0];
                my $right_beg_phase = $initial_right_seg->{phase_beg};
                my $right_rel_lend = $initial_right_seg->{rel_lend};
                
                my $right_cds_part = substr($cds_right_seq, $right_rel_lend - 1);

                my $fusion_seq = join("", lc($left_cds_part), uc($right_cds_part));
                my $pep = translate_sequence($fusion_seq, &get_translation_phase($initial_left_seg->{phase_beg}));
                
                my $prot_fusion_type = "NA";
                if ($left_end_phase ne '.' && $right_beg_phase ne '.') {
                    $prot_fusion_type = ( ($left_end_phase + 1) % 3 == $right_beg_phase) ? "INFRAME" : "FRAMESHIFT";
                }
                
                my $left_segs_string = &segments_to_string(@$left_fuse_segments_aref);
                my $right_segs_string = &segments_to_string(@$right_fuse_segments_aref);

                my $left_domains_string = &get_pfam_domains($cds_left_id, $left_rel_rend, "left", $pfam_domain_db_tied_hash) || ".";
                my $right_domains_string = &get_pfam_domains($cds_right_id, $right_rel_lend, "right", $pfam_domain_db_tied_hash) || ".";
                
                
                push (@results, { cds_left_id => $cds_left_id,
                                  cds_right_id => $cds_right_id,
                                  cds_left_range => "1-$left_rel_rend",
                                  cds_right_range => "$right_rel_lend-" . length($cds_right_seq),
                                  prot_fusion_type => $prot_fusion_type,
                                  cds_fusion_seq => $fusion_seq,
                                  prot_fusion_seq => $pep,
                                  fusion_coding_descr => join("<==>", $left_segs_string, $right_segs_string),
                                  left_domains => $left_domains_string,
                                  right_domains => $right_domains_string,
                                  
                      }
                    );
                
            }
        }
    }
    
    return (@results);
}

