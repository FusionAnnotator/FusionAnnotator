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


my %PROBLEM_GENE;

main: {

    if (1) {

        print STDERR "-Validating each gene separately\n";
        
        # examining each gene separately
        my @genes;
        if ($fusion_name) {
            my ($geneA, $geneB) = split(/--/, $fusion_name);
            @genes = ($geneA, $geneB);
        }
        else {
            @genes = $prot_info_db_tied_hash->get_keys();
        }
        foreach my $gene (@genes) {
            &validate_gene_reading_frames($gene);
        }
        
        
    }


    print STDERR "\n\n-Exploring all fusion pairs\n\n";
    
    # examining fusions
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
    

            push (@left_breaks, &generate_random_breaks($cds_left_obj, 3));
            push (@right_breaks, &generate_random_breaks($cds_right_obj, 3));
            
        
            foreach my $break_left (@left_breaks) {
                foreach my $break_right (@right_breaks) {

                    # returned segment lists are in order according of 5' to 3' based on strand.
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

                        my $pep_left = translate_sequence($left_cds_part, &get_translation_phase($initial_left_seg->{phase_beg})) || "";
                        my $pep_right = translate_sequence($right_cds_part, &get_translation_phase($initial_right_seg->{phase_beg})) || "";
                        
                        
                        my $prot_fusion_type = "NA";
                        if ($left_end_phase ne '.' && $right_beg_phase ne '.') {
                            $prot_fusion_type = ( ($left_end_phase + 1) % 3 == $right_beg_phase) ? "INFRAME" : "FRAMESHIFT";
                        }
                        
                        my $left_segs_string = &segments_to_string(@$left_fuse_segments_aref);
                        my $right_segs_string = &segments_to_string(@$right_fuse_segments_aref);
                        
                        my $result = {cds_left_id => $cds_left_id,
                                          cds_right_id => $cds_right_id,
                                          cds_left_range => "1-$left_rel_rend",
                                          cds_right_range => "$right_rel_lend-" . length($cds_right_seq),
                                          prot_fusion_type => $prot_fusion_type,
                                          cds_fusion_seq => $fusion_seq,
                                          prot_fusion_seq => $pep,
                                          fusion_coding_descr => join("<==>", $left_segs_string, $right_segs_string),
                                          pep_left => $pep_left,
                                          pep_right => $pep_right
                        };

                        
                        my $STATUS = "OK";
                        unless ($PROBLEM_GENE{$gene_left} || $PROBLEM_GENE{$gene_right}) {
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
                                   $result->{pep_left},
                                   $result->{pep_right}
                                   
                            ) . "\n";
                    }
                }
            }
        }
    }
    
    return;
    
}


####
sub validate_gene_reading_frames {
    my ($gene) = @_;

    my $cds_objs_aref = &decode_my_json($gene, $prot_info_db_tied_hash, 1); 

    foreach my $cds_obj (@$cds_objs_aref) {
        &validate_cds_obj_reading_frames($gene, $cds_obj);
    }
}

####
sub validate_cds_obj_reading_frames {
    my ($gene, $cds_obj) = @_;

    my @segments = sort {$a->{lend}<=>$b->{lend}} @{$cds_obj->{phased_segments}};
    my $orient = $segments[0]->{orient};

    if ($orient eq '-') {
        @segments = reverse @segments;
    }
    
    my $cds_seq = $cds_obj->{cds_seq};
    my $phase_seq = "";;
    foreach my $i (0..(length($cds_seq)-1)) {
        $phase_seq .= $i % 3;
    }

    my $prev_seg = undef;
    my $first_pep = "";
    my $initial_phase = $segments[0]->{phase_beg};
    while (@segments)  {
        my $segs_string = &segments_to_string(@segments);
        my $rel_lend = $segments[0]->{rel_lend};
        my $seg_seq = substr($cds_seq, $rel_lend - 1);

        my $seg_phase_seq = substr($phase_seq, $rel_lend -1);
        my $phase_beg = $segments[0]->{phase_beg};
        my $theor_phase = ($rel_lend + $initial_phase - 1) % 3;
        if ($theor_phase != $phase_beg) {
            confess "Error, phase_beg not $theor_phase: " . Dumper($segments[0]);
        }
        my $rel_rend = $segments[0]->{rel_rend};
        my $phase_end = $segments[0]->{phase_end};
        my $theor_phase_end = ($phase_beg + $rel_rend - $rel_lend +1 -1) % 3;
        if ($theor_phase_end != $phase_end) {
            confess "Error, phase_end not $theor_phase_end: " . Dumper($segments[0]);
        }
        
        my $translate_phase = &get_translation_phase($phase_beg);
        
        my $pep = translate_sequence($seg_seq, $translate_phase);
        
        if ($pep) {
            $pep =~ s/\*$//; # trim stop
        }
        else {
            #die "Error, no pep for seq [$seg_seq]";
            $pep = "";
        }
        if (!$first_pep) {
            $first_pep = $pep;
            
            
            if ($first_pep =~ /\*/) {
                $PROBLEM_GENE{$gene} = 1;
            }
        }
        elsif ($first_pep !~ /\*/ && $pep =~ /\*/) {
            die "Error, stop introduced in translation of downstream segments ";
        }
        
        
        print join("\t", $segs_string, "$rel_lend-$rel_rend...", $phase_beg, "?$theor_phase", $seg_phase_seq, $seg_seq, $pep) . "\n";

        if ($prev_seg) {
            my $prev_phase_end = $prev_seg->{phase_end};
            print "Phase transitions: $prev_phase_end <-> $phase_beg\n";

            if ( ($prev_phase_end + 1) % 3 != $phase_beg) {
                die "Error, phases out of order";
            }
        }
        
        $prev_seg = shift @segments;
    }
    print "\n";
}


####
sub generate_random_breaks {
    my ($cds_obj, $num_breaks) = @_;

    
    my @segments = sort {$a->{lend}<=>$b->{lend}} @{$cds_obj->{phased_segments}};

    my $chr = $segments[0]->{chr};
    my $orient = $segments[0]->{orient};
    
    my @breaks;
    
    for(1..$num_breaks) {

        my $rand_segment = $segments[int(rand(scalar(@segments)))];
        my ($lend, $rend) = ($rand_segment->{lend}, $rand_segment->{rend});
        my $seg_len = $rend - $lend + 1;
        my $random_break = $lend + 1 + int(rand($seg_len-1));
        push (@breaks, "$chr:$random_break:$orient");
    }

    return(@breaks);
}

