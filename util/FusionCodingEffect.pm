#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Carp::Assert;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use Data::Dumper;
use JSON::XS;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use TiedHash;
use Nuc_translator;

my $DEBUG = 0;


####
sub try_fuse_cdss {
    my ($cds_left_obj, $break_left, $cds_right_obj, $break_right) = @_;

    
    # get left part
    my @left_fusion_partner_segments = &get_left_fusion_partner_segments($cds_left_obj, $break_left);
        
    # todo: get right part
    my @right_fusion_partner_segments = &get_right_fusion_partner_segments($cds_right_obj, $break_right);
    
        
    ## piece it together.
    
    #print STDERR "Left: " . Dumper(\@left_fusion_partner_segments) . "\nRight: " . Dumper(\@right_fusion_partner_segments);
    
    return (\@left_fusion_partner_segments, \@right_fusion_partner_segments);
    
}


####
sub get_left_fusion_partner_segments {
    my ($cds_obj, $breakpoint_info) = @_;
    
    my ($chr, $breakpoint_coord, $orient) = split(/:/, $breakpoint_info);
    
    # ensure breakpoint overlaps a coding segment
    unless (&breakpoint_overlaps_cds_segment($cds_obj, $breakpoint_coord)) {
        return();
    }
    
    my ($left_segs_aref, $right_segs_aref) = &split_cds_at_breakpoint($cds_obj, $breakpoint_coord);

    if ($orient eq '+') {
        return(@$left_segs_aref);
    }
    else {
        return(reverse @$right_segs_aref);
    }
}

####
sub get_right_fusion_partner_segments {
    my ($cds_obj, $breakpoint_info) = @_;
    
    my ($chr, $breakpoint_coord, $orient) = split(/:/, $breakpoint_info);

    # ensure breakpoint overlaps a coding segment
    unless (&breakpoint_overlaps_cds_segment($cds_obj, $breakpoint_coord)) {
        return();
    }
    
    my ($left_segs_aref, $right_segs_aref) = &split_cds_at_breakpoint($cds_obj, $breakpoint_coord);
    
    if ($orient eq '+') {
        return(@$right_segs_aref);
    }
    else {
        return(reverse @$left_segs_aref);
    }
}



####
sub breakpoint_overlaps_cds_segment {
    my ($cds_obj, $breakpoint_coord) = @_;

    #print STDERR Dumper($cds_obj);
    
    my @segments = sort {$a->{lend}<=>$b->{lend}} @{$cds_obj->{phased_segments}};
    
    foreach my $segment (@segments) {
        if ($segment->{lend} <= $breakpoint_coord && $breakpoint_coord <= $segment->{rend}) {
            return(1);
        }
    }

    return(0); # no overlap
}


####
sub split_cds_at_breakpoint {
    my ($cds_obj, $breakpoint_coord) = @_;
    
    my @segments = sort {$a->{lend}<=>$b->{lend}} @{$cds_obj->{phased_segments}};
    
    my @segs_left;
    my @segs_right;

    
    foreach my $segment (@segments) {
        if ($segment->{rend} <= $breakpoint_coord) {
            push (@segs_left, $segment);
        }
        elsif ($segment->{lend} >= $breakpoint_coord) {
            push (@segs_right, $segment);
        }
        elsif(&overlaps_breakpoint($breakpoint_coord ,$segment->{lend}, $segment->{rend})) {

            ## split the segment at the breakpoint, keep breakpoint coordinate in each piece.
            my $orient = $segment->{orient};
            if ($orient eq '+') {
                my $new_left_segment = { chr => $segment->{chr}, 
                                         lend => $segment->{lend},
                                         rend => $breakpoint_coord,
                                         orient => $orient,
                                         rel_lend => $segment->{rel_lend},
                                         rel_rend => $segment->{rel_lend} + ($breakpoint_coord - $segment->{lend}),
                                         phase_beg => $segment->{phase_beg},
                                         phase_end => ".", # set below
                };
                if ($segment->{phase_beg} ne '.') {
                    $new_left_segment->{phase_end} = ($segment->{phase_beg} + ($breakpoint_coord - $segment->{lend})) % 3;
                }
                
                my $new_right_segment = { chr => $segment->{chr},
                                          lend => $breakpoint_coord,
                                          rend => $segment->{rend},
                                          orient => $orient,
                                          rel_lend => $new_left_segment->{rel_rend},
                                          rel_rend => $segment->{rel_rend},
                                          phase_beg => $new_left_segment->{phase_end},
                                          phase_end => $segment->{phase_end},
                };
                
                   
                push (@segs_left, $new_left_segment);
                push (@segs_right, $new_right_segment);
                
                
                print "Original segment: " . &segments_to_string($segment) . 
                    "\nNew frags: " . &segments_to_string($new_left_segment) . "\t" . &segments_to_string($new_right_segment) . "\n\n" if $DEBUG;
                
                
                if ($segment->{phase_beg} ne '.') {
                    my $supposed_end_phase = ($new_right_segment->{phase_beg} + ($new_right_segment->{rel_rend} - $new_right_segment->{rel_lend})) % 3;
                    print "supposed end phase: $supposed_end_phase\n" if $DEBUG;
                    assert($supposed_end_phase == $new_right_segment->{phase_end});
                }
                
             
            }
            else {
                ## orient eq '-'
                
                my $new_right_segment = { chr => $segment->{chr},
                                          lend => $breakpoint_coord,
                                          rend => $segment->{rend},
                                          orient => $orient,
                                          rel_lend => $segment->{rel_rend} + ($segment->{rend} - $breakpoint_coord),
                                          rel_rend => $segment->{rel_rend},
                                          phase_beg => $segment->{phase_beg},
                                          phase_end => ".", # set below
                };
                
                if ($segment->{phase_beg} ne '.') {
                    $new_right_segment->{phase_end} = ($segment->{phase_beg} + (($segment->{rend} - $breakpoint_coord))) % 3;
                }

                my $new_left_segment = { chr => $segment->{chr},
                                         lend => $segment->{lend},
                                         rend => $breakpoint_coord,
                                         orient => $orient,
                                         rel_lend => $segment->{rel_lend},
                                         rel_rend => $new_right_segment->{rel_lend},
                                         phase_beg => $new_right_segment->{phase_end},
                                         phase_end => $segment->{phase_end},
                };
                
                push (@segs_left, $new_left_segment);
                push (@segs_right, $new_right_segment);
                
                print "Original segment: " . &segments_to_string($segment) . 
                    "\nNew frags: " . &segments_to_string($new_left_segment) . "\t" . &segments_to_string($new_right_segment) . "\n\n" if $DEBUG;
                
                if ($segment->{phase_beg} ne ".") {
                    assert($new_right_segment->{phase_end} == ($new_right_segment->{phase_beg} + $new_right_segment->{rend} - $new_right_segment->{lend}) % 3);
                    
                    assert($new_left_segment->{phase_end} == ($new_left_segment->{phase_beg} + $new_left_segment->{rend} - $new_left_segment->{lend}) % 3);
                }
            }                
            
        }
        else {
            die "Error, shouldn't get here";
        }
    }
    
    return(\@segs_left, \@segs_right);
    
}

####
sub overlaps_breakpoint {
    my ($breakpoint_coord, $lend, $rend) = @_;

    if ($breakpoint_coord >= $lend && $breakpoint_coord <= $rend) {
        return(1);
    }
    else {
        return(0);
    }
}

####
sub segments_to_string {
    my (@segments) = @_;

    @segments = sort {$a->{lend}<=>$b->{lend}} @segments;

    my $chr = $segments[0]->{chr};
    my $orient = $segments[0]->{orient} or die "Error: " . Dumper(\@segments);

    my @coord_text;
    foreach my $segment (@segments) {
        my ($phase_left, $phase_right) = ($orient eq '+') 
            ? ($segment->{phase_beg}, $segment->{phase_end})
            : ($segment->{phase_end}, $segment->{phase_beg});
        
        push (@coord_text, "[$phase_left]" . join("-", $segment->{lend}, $segment->{rend}) . "[$phase_right]");
    }
    
    my $descr_text = join("|", $chr, $orient, @coord_text);

    return($descr_text);
    
}


####
sub get_pfam_domains {
    my ($cds_id, $cds_coord, $left_or_right_side, $pfam_domain_db_tied_hash) = @_;

    unless ($pfam_domain_db_tied_hash) {
        return;
    }

    my $pfam_hits_aref = &decode_my_json($cds_id, $pfam_domain_db_tied_hash);

    unless ($pfam_hits_aref) { return; }
    
    my @pfam_hits = @$pfam_hits_aref;
    
    my @pfam_domains_selected;

    foreach my $pfam_hit (@pfam_hits) {
        my ($start, $end) = ($pfam_hit->{query_start}, $pfam_hit->{query_end});

        if  (($left_or_right_side eq 'left' && $end <= $cds_coord)
             ||
             ($left_or_right_side eq 'right' && $start >= $cds_coord) ) {

            ## domain entirely on the side of the protein included in the fusion.
            
            push (@pfam_domains_selected, $pfam_hit);
        }
        elsif ($start < $cds_coord && $cds_coord < $end) {
            ## overlaps

            ## fragment it and return the fragment.
            my $pfam_hit_copy = &clone($pfam_hit);
            if ($left_or_right_side eq 'left') {
                $pfam_hit_copy->{query_end} = $cds_coord;
                $pfam_hit_copy->{query_end_partial} = 1;
            }
            else {
                # right side
                $pfam_hit_copy->{query_start} = $cds_coord;
                $pfam_hit_copy->{query_start_partial} = 1;
            }
            $pfam_hit_copy->{hmmer_domain} .= "-PARTIAL";
            push (@pfam_domains_selected, $pfam_hit_copy);
            
        }
    }

    ## generate a summary string.
    @pfam_domains_selected = sort {$a->{query_start}<=>$b->{query_start}} @pfam_domains_selected;

    my @pfam_descrs;
    foreach my $pfam_domain (@pfam_domains_selected) {
        #'query_start' => '369',
        #'cds_id' => 'DISP1|ENST00000284476.6',
        #'domain_evalue' => '7.3e-21',
        #'query_end' => '733',
        #'hmmer_domain' => 'Patched'

        if ($pfam_domain->{query_start_partial}) {
            $pfam_domain->{query_start} = "~" . $pfam_domain->{query_start};
        }
        if ($pfam_domain->{query_end_partial}) {
            $pfam_domain->{query_end} = $pfam_domain->{query_end} . "~";
        }
        
        my $descr = join("|", $pfam_domain->{hmmer_domain},
                         $pfam_domain->{query_start} . "-" . $pfam_domain->{query_end},
                         $pfam_domain->{domain_evalue});

        push (@pfam_descrs, $descr);
    }

    my $ret_descr = join("^", @pfam_descrs);

    return($ret_descr);
            
}


####
sub clone {
    my ($hashref) = @_;

    my $clone_ref = {};
    foreach my $key (keys %$hashref) {

        $clone_ref->{$key} = $hashref->{$key};
    }

    return($clone_ref);
}


####
sub decode_my_json {
    my ($key, $prot_info_db_tied_hash, $report_failed_retrievals_flag) = @_;
    
    my $coder = new JSON::XS();

    my $decoded = undef;
    my $json = undef;
    
    eval {
        $json = $prot_info_db_tied_hash->get_value($key);

        if ($json) {
            $decoded = $coder->decode($json);
            #print Dumper($decoded);
        }
        else {
            if ($report_failed_retrievals_flag) {
                print STDERR "WARNING, no entry stored in dbm for [$key]\n";
            }
            return(undef);
        }
    };

    if ($@) {
        print STDERR "WARNING, key: $key returns json: $json and error decoding: $@";
    }

    return($decoded);
}


####
sub get_candidate_breaks {
    my ($gene_obj, $fusion_side) = @_;

    my @phased_segments = sort {$a->{lend} <=> $b->{lend}} @{$gene_obj->{phased_segments}};
    
    my $orient = $phased_segments[0]->{orient};
    my $chr = $phased_segments[0]->{chr};
    my %breaks;

    my $first_segment = $phased_segments[0];
    my $last_segment = $phased_segments[$#phased_segments];
    
    foreach my $segment (@phased_segments) {
        my ($lend, $rend) = ($segment->{lend}, $segment->{rend});

        $lend = "$chr:$lend:$orient";
        $rend = "$chr:$rend:$orient";
        
        if ($orient eq '+') {
            if ($fusion_side eq 'left' && $segment ne $last_segment) {
                $breaks{$rend}++;
            }
            elsif ($segment ne $first_segment) {
                $breaks{$lend}++;
            }
        }
        elsif ($orient eq '-') {
            if ($fusion_side eq 'left' && $segment ne $first_segment) {
                $breaks{$lend}++;
            }
            elsif ($segment ne $last_segment) {
                $breaks{$rend}++;
            }
        }
    }

    return(keys %breaks);
}

        
; #EOM

