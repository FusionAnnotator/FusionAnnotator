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
#  --fusion <string>                 :geneA--geneB  (ex. "ENSG00000143924.14--ENSG00000171094.11")
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
my $fusion_name;
my $genome_lib_dir;
my $show_all_flag = 0;

&GetOptions ( 'h' => \$help_flag,
              'fusion=s' => \$fusion_name,
              'genome_lib_dir=s' => \$genome_lib_dir,
   
              'show_all' => \$show_all_flag,
    );

unless ($genome_lib_dir) {
    
    die $usage;
}

my $prot_info_db = "$genome_lib_dir/ref_annot.prot_info.dbm";

unless (-s $prot_info_db) {
    confess "Error, cannot locate $prot_info_db - be sure to use the latest version of the CTAT genome lib";
}

my $prot_info_db_tied_hash = new TiedHash( { 'use' => $prot_info_db } );

main: {

    if (1) {
        my @genes = $prot_info_db_tied_hash->get_keys();
        foreach my $gene (@genes) {
            &validate_gene_reading_frames($gene);
        }
        
        exit(0);
    }
    
    
    if ($fusion_name) {
        my ($geneA, $geneB) = split(/--/, $fusion_name);
        &examine_gene_pair($geneA, $geneB);
    }
    else {
        my @genes = $prot_info_db_tied_hash->get_keys();
        for (my $i = 0; $i < $#genes; $i++) {

            for (my $j = $i + 1; $j <= $#genes; $j++) {

                my $gene_i = $genes[$i];
                my $gene_j = $genes[$j];
                
                &examine_gene_pair($gene_i, $gene_j);
            }
        }
    }
    

    exit(0);
}


sub examine_gene_pair {
    my ($gene_left, $gene_right) = @_;
    
    my $fusion_name = join("--", $gene_left, $gene_right);
    
    my $gene_left_aref = &decode_my_json($gene_left, $prot_info_db_tied_hash, 1);
    my $gene_right_aref = &decode_my_json($gene_right, $prot_info_db_tied_hash, 1);
    
    unless ($gene_left_aref && $gene_right_aref) {
        die "Error, didn't find both genes in db";
    }

    if (0) {
        ## debugging
        foreach my $cds_obj (@$gene_left_aref, @$gene_right_aref) {
            my @segments = @{$cds_obj->{phased_segments}};
            print &segments_to_string(@segments) . "\n";
        }
    }
    

    my @results;
    
    foreach my $cds_left_obj (@$gene_left_aref) {
        
        my $cds_left_seq = $cds_left_obj->{cds_seq};
        my $cds_left_id = $cds_left_obj->{cds_id};
        
        foreach my $cds_right_obj (@$gene_right_aref) {
            
            my $cds_right_seq = $cds_right_obj->{cds_seq};
            my $cds_right_id = $cds_right_obj->{cds_id};
    

            my @left_breaks = &get_candidate_breaks($cds_left_obj, 'left');
            my @right_breaks = &get_candidate_breaks($cds_right_obj, 'right');
            
            foreach my $break_left (@left_breaks) {
                foreach my $break_right (@right_breaks) {

                    # returned segment lists are in order according of 5' to 3' based on strand.
                    my ($left_fuse_segments_aref, $right_fuse_segments_aref) = &try_fuse_cdss($cds_left_obj, $break_left, $cds_right_obj, $break_right);
            
                    if (@$left_fuse_segments_aref && @$right_fuse_segments_aref) {
                        
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
                        my $pep = translate_sequence($fusion_seq, 1);
                        
                        my $prot_fusion_type = "NA";
                        if ($left_end_phase ne '.' && $right_beg_phase ne '.') {
                            $prot_fusion_type = ( ($left_end_phase + 1) % 3 == $right_beg_phase) ? "INFRAME" : "FRAMESHIFT";
                        }
                        
                        my $left_segs_string = &segments_to_string(@$left_fuse_segments_aref);
                        my $right_segs_string = &segments_to_string(@$right_fuse_segments_aref);
                        
                        push (@results, { cds_left_id => $cds_left_id,
                                          cds_right_id => $cds_right_id,
                                          cds_left_range => "1-$left_rel_rend",
                                          cds_right_range => "$right_rel_lend-" . length($cds_right_seq),
                                          prot_fusion_type => $prot_fusion_type,
                                          cds_fusion_seq => $fusion_seq,
                                          prot_fusion_seq => $pep,
                                          fusion_coding_descr => join("<==>", $left_segs_string, $right_segs_string),
                                          
                              }
                            );
                        
                    }
                }
            }
        }
    }



    if (@results) {
        
        foreach my $result (@results) {

            my $STATUS = "OK";
            if ($result->{prot_fusion_seq} =~ /\*\S/) {
                if ($result->{prot_fusion_type} ne "FRAMESHIFT") {

                    ## see if it has just a single -inframe stop
                    my $prot_fusion_seq = $result->{prot_fusion_seq};
                    chop $prot_fusion_seq;
                    my $num_internal_stops = 0;
                    while ($prot_fusion_seq =~ /\*/g) {
                        $num_internal_stops++;
                    }
                    if ($num_internal_stops == 1) {
                        $STATUS = "INFRAMESTOP";
                    }
                    else {
                        $STATUS = "ERROR-$num_internal_stops";
                    }
                }
            }
            elsif ($result->{prot_fusion_seq} =~ /^M[^\*]+\*/ && $result->{prot_fusion_type} eq "FRAMESHIFT") {
                my $cds_fusion_seq = $result->{cds_fusion_seq};
                $cds_fusion_seq =~ /[a-z]([A-Z]+)$/;
                my $fusion_ptB_len = length($1);
                my $stop_codon = substr($cds_fusion_seq, -3);
                if (length($1) > 100 && $stop_codon =~ /^(TGA|TAA|TAG)$/) {
                    $STATUS = "ERROR-notFrameshift";
                }
            }
            
            elsif ($result->{prot_fusion_type} eq "FRAMESHIFT") {
                #$STATUS = "ERROR";
            }
            
            print join("\t",
                       $STATUS,
                       $fusion_name,
                       $result->{cds_left_id}, $result->{cds_left_range},
                       $result->{cds_right_id}, $result->{cds_right_range},
                       $result->{prot_fusion_type},
                       $result->{fusion_coding_descr},
                       $result->{cds_fusion_seq},
                       $result->{prot_fusion_seq},
                       
                ) . "\n";
        }
    }
    
    
    return;
    
}


####
sub validate_gene_reading_frames {
    my ($gene) = @_;

    my $cds_objs_aref = &decode_my_json($gene, $prot_info_db_tied_hash, 1); 

    foreach my $cds_obj (@$cds_objs_aref) {
        &validate_cds_obj_reading_frames($cds_obj);
    }
}

####
sub validate_cds_obj_reading_frames {
    my ($cds_obj) = @_;

    
