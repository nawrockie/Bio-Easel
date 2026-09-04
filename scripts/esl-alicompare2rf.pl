#!/usr/bin/env perl
# 
# esl-alicompare2rf.pl: given a RF annotated Stockholm alignment, output differences between 
#                       each sequence and the RF annotation.
#
# EPN, Fri Jan  8 10:59:21 2021
# 

use strict;
use Getopt::Long;
use Bio::Easel::MSA;

my $version     = "0.18";
my $date        = "Sep 2026";

my $in_alifile  = "";    # name of input MSA file

my $usage;
$usage  = "# esl-alicompare2rf.pl :: output differences between aligned sequences and RF\n";
$usage .= "# Bio-Easel $version ($date)\n";
$usage .= "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
$usage .= "\n";
$usage .= "Usage: esl-alicompare2rf.pl <alignment file with RF anntotation>\n";
#$usage .= "Usage: esl-alicompare2rf.pl [OPTIONS] <alignment file with RF anntotation>\n";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t--alldel    : include terminal deletions, [df: do not]\n";
$usage .= "\t\t--seqrf <s> : set RF to compare to as seq <s>\n";
#$usage .= "\t\t-w       : use sequence weights in the alignment file to weight counts [df: do not]\n";
#$usage .= "\t\t--nc <f> : ignore PPs, remove consensus basepairs for which > <f> fraction are non-canonical\n";
#$usage .= "\t\t--dg <f> : ignore PPs, remove consensus basepairs for which > <f> fraction of seqs are double gaps\n";

my $do_alldel = 0;
my $seqrf     = undef;
&GetOptions( "alldel"  => \$do_alldel,
             "seqrf=s" => \$seqrf);
#             "w"    => \$use_weights, 
#             "nc=s" => \$min_fractnc,
#             "dg=s" => \$min_fractdg);

if(scalar(@ARGV) != 1) { die $usage; }
($in_alifile) = @ARGV;

# validate input args
if(! -e $in_alifile) { die "ERROR $in_alifile does not exist"; }

# open file 
my $msa = Bio::Easel::MSA->new({ fileLocation => $in_alifile,
                                 forceText => 1,
                               });

# check if we have RF, we need it unless --seqrf
my $has_rf = $msa->has_rf;
if((! defined $seqrf) && (! $has_rf)) { die "ERROR, if --seqrf not used, alignment must have RF annotation, it does not"; }

my $alen = $msa->alen;
my $nseq = $msa->nseq; 
my $i    = 0; # counter over sequences 

# get RF or the sequence that you want to use as RF
my $rf_str = undef;
my $gc_rf_str = ($has_rf) ? $msa->get_rf : undef;
if(! defined $seqrf) { 
  $rf_str = $gc_rf_str;
}
  
else { # find the sequence
  for($i = 0; $i < $nseq; $i++) { 
    my $seq_name = $msa->get_sqname($i);
    if($seq_name eq $seqrf) {
      $rf_str = $msa->get_sqstring_aligned($i);
    }
  }
  if(! defined $rf_str) {
    die "ERROR, --seqrf used, but did not find a sequence named $seqrf in alignment.";
  }
}

my @gc_rf_A = ();
my $use_gc_rf = 0;
if($has_rf && (defined $seqrf)) {
  @gc_rf_A = split("", $gc_rf_str);
  $use_gc_rf = 1;
}
my @rf_A = split("", $rf_str);
if(scalar(@rf_A) != $alen) { 
  die "ERROR unexpected alignment length mismatch $alen != %d\n";
}
my $i;

#printf("%-30s  %5s  %5s  %5s  %6s  %6s  description\n", 
#         "#seqname", "rfpos", "sqpos", "apos", "rfchar", "sqchar");
if(defined $seqrf) {
  print("# Reference (RF) set as $seqrf\n");
}
else { 
  print("# Reference (RF) set as #=GC RF from alignment $seqrf\n");
}
my $optional_field = ($use_gc_rf) ? "\tGC_RF_is_gap?" : "";
printf("%s\t%s\t%s\t%s\t%s\t%s\t%s%s\n",
       "#seqname", "rfpos", "sqpos", "apos", "rfchar", "sqchar", "description", $optional_field);

# for each sequence, go through each position and determine differences with RF
# we save these lines to an array so we can go back through the array and combine
# any deletions or insertions with lengths greater than 1 into a single line.
for($i = 0; $i < $nseq; $i++) { 
  my $seq_name = $msa->get_sqname($i);
  my $asqstring = $msa->get_sqstring_aligned($i);
  my @asqstring_A = split("", $asqstring);
  my $rfpos = 0;
  my $sqpos = 0;
  my @out_str_A  = (); # array of output strings, all information except the 'description' field
  my @out_desc_A = (); # array of output descriptions
  my @out_opt_A  = (); # array of optional fields for output, only filled if --seqrf and alignment has GC RF annotation
  
  # determine first and final apos with a nongap residue
  my $apos = 0;
  my $spos = $alen+1;
  my $epos = 0;
  for($apos = 1; $apos <= $alen; $apos++) { 
    if ($asqstring_A[($apos-1)] =~ m/[A-Z]/) {
      $spos = $apos;
      $apos = $alen+1; # breaks loop
    }
  }
  for($apos = $alen; $apos >= 1; $apos--) { 
    if ($asqstring_A[($apos-1)] =~ m/[A-Z]/) {
      $epos = $apos;
      $apos = 0; # breaks loop
    }
  }

  my $gc_rf_is_gap;
  for($apos = 1; $apos <= $alen; $apos++) { 
    my $sqchar    = $asqstring_A[($apos-1)];
    my $rfchar    = $rf_A[($apos-1)];
    my $gc_rfchar = ($use_gc_rf) ? $gc_rf_A[($apos-1)] : undef;
    my $tmp_sqchar = $sqchar;
    my $tmp_rfchar = $rfchar;
    $tmp_sqchar =~ tr/a-z/A-Z/; # uppercase-ize
    $tmp_rfchar =~ tr/a-z/A-Z/; # uppercase-ize
    my $sq_is_gap    = ($tmp_sqchar =~ m/[A-Z]/) ? 0 : 1;
    my $rf_is_gap    = ($tmp_rfchar =~ m/[A-Z]/) ? 0 : 1;
    $gc_rf_is_gap = 0;
    if(defined $gc_rfchar) { 
      $gc_rfchar =~ tr/a-z/A-Z/; # uppercase-ize
      $gc_rf_is_gap = ($gc_rfchar =~ m/[A-Z]/) ? 0 : 1;
    }
    my $desc = undef;
    if(! $rf_is_gap) { $rfpos++; }
    if(! $sq_is_gap) { $sqpos++; }

    if($rf_is_gap) { 
      if(! $sq_is_gap) { 
        $desc = "insert-after-RF-position";
      }
    }
    else { # rf is not a gap
      if($sq_is_gap) { 
        if(($do_alldel) || (($apos >= $spos) && ($apos <= $epos))) { 
          $desc = "deletion";
        }
      }
      elsif($tmp_rfchar ne $tmp_sqchar) { 
        $desc = "substitution";
      }
    }
    if(defined $desc) { 
      # if --seqrf and we have RF annotation, output extra column indicating if GC RF position is a gap or not
      push(@out_str_A, sprintf("%s\t%d\t%d\t%d\t%s\t%s", 
                               $seq_name, $rfpos, $sqpos, $apos, $rfchar, $sqchar));
      push(@out_desc_A, $desc);

      if($use_gc_rf) { 
        push(@out_opt_A, ($gc_rf_is_gap) ? "y" : "n");
      }
    }
  }

  # output for all RF positions, combining indels length > 1 into one line,
  # we can go to alen (rflen must be <= alen)
  my $nlines = scalar(@out_str_A);
  my $prv_desc = undef;
  my ($cur_seq_name, $cur_rfpos,  $cur_sqpos, $cur_apos, $cur_rfchar, $cur_sqchar, $cur_desc);
  my $cur_sqpos_start = undef;
  my $cur_sqpos_end   = undef;
  my $cur_apos_start  = undef;
  my $cur_apos_end    = undef;
  my $cur_rfpos_start = undef;
  my $cur_rfpos_end   = undef;
  my $cur_sqchar_str  = "";
  my $cur_rfchar_str  = "";
  my $cur_opt         = undef; # current optional field
  my $cur_opt_str     = undef; # current optional field
  for(my $l = 0; $l < $nlines; $l++) { 
    my $cur_str = $out_str_A[$l];
    $cur_desc = $out_desc_A[$l];
    $cur_opt  = ($use_gc_rf) ? $out_opt_A[$l] : undef;
    
    # if substitution - output it
    # if deletion, keep going until sqpos changes, then output summary line
    # if insertion, keep going until rfpos changes, then output summary line
    if($cur_desc eq "substitution") {
      # output any insertion or deletion strings we have
      if(defined $cur_sqpos_start) {
        output_summary_insertion_string($cur_seq_name, $cur_rfpos, $cur_sqpos_start, $cur_sqpos_end, $cur_apos_start, $cur_apos_end, $cur_sqchar_str, $cur_opt_str);
        ($cur_sqpos_start, $cur_sqpos_end, $cur_apos_start, $cur_apos_end) = (undef, undef, undef, undef);
        $cur_sqchar_str = "";
      }
      if(defined $cur_rfpos_start) {
        output_summary_deletion_string($cur_seq_name, $cur_rfpos_start, $cur_rfpos_end, $cur_sqpos, $cur_apos_start, $cur_apos_end, $cur_rfchar_str, $cur_opt_str);
        ($cur_rfpos_start, $cur_rfpos_end, $cur_apos_start, $cur_apos_end) = (undef, undef, undef, undef);
        $cur_rfchar_str = "";
      }
      printf($cur_str . "\tsubstitution%s\n", (defined $cur_opt) ? "\t" . $cur_opt : "");
    }
    elsif($cur_desc eq "deletion") {
      # output any insertion strings we have
      if(defined $cur_sqpos_start) {
        output_summary_insertion_string($cur_seq_name, $cur_rfpos, $cur_sqpos_start, $cur_sqpos_end, $cur_apos_start, $cur_apos_end, $cur_sqchar_str, $cur_opt_str);
        ($cur_sqpos_start, $cur_sqpos_end, $cur_apos_start, $cur_apos_end) = (undef, undef, undef, undef);
        $cur_sqchar_str = "";
      }

      my @el_A = split(/\t/, $cur_str);
      if(scalar(@el_A) != 6) { die "ERROR internal coding problem"; }
      ($cur_seq_name, $cur_rfpos,  $cur_sqpos, $cur_apos, $cur_rfchar, $cur_sqchar) = @el_A;

      # build up the deletion line
      if(! defined $cur_rfpos_start) {
        $cur_rfpos_start = $cur_rfpos;
        $cur_apos_start  = $cur_apos;
        $cur_rfchar_str  = $cur_rfchar;
        if(defined $cur_opt) { 
          $cur_opt_str = $cur_opt;
        }
      }
      else {
        $cur_rfchar_str .= $cur_rfchar;
        if(defined $cur_opt) { 
          $cur_opt_str .= $cur_opt;
        }
      }
      $cur_rfpos_end = $cur_rfpos;
      $cur_apos_end  = $cur_apos;
      # sqpos doesn't change
    }
    elsif($cur_desc eq "insert-after-RF-position") {
      # output any deletion strings we have
      if(defined $cur_rfpos_start) {
        output_summary_deletion_string($cur_seq_name, $cur_rfpos_start, $cur_rfpos_end, $cur_sqpos, $cur_apos_start, $cur_apos_end, $cur_rfchar_str, $cur_opt_str);
        ($cur_rfpos_start, $cur_rfpos_end, $cur_apos_start, $cur_apos_end) = (undef, undef, undef, undef);
        $cur_rfchar_str = "";
      }

      my @el_A = split(/\t/, $cur_str);
      if(scalar(@el_A) != 6) { die "ERROR internal coding problem"; }
      ($cur_seq_name, $cur_rfpos,  $cur_sqpos, $cur_apos, $cur_rfchar, $cur_sqchar) = @el_A;

      # build up the insertion line
      if(! defined $cur_sqpos_start) {
        $cur_sqpos_start = $cur_sqpos;
        $cur_apos_start  = $cur_apos;
        $cur_sqchar_str  = $cur_sqchar;
        if(defined $cur_opt) { 
          $cur_opt_str = $cur_opt;
        }
      }
      else {
        $cur_sqchar_str .= $cur_sqchar;
        if(defined $cur_opt) { 
          $cur_opt_str .= $cur_opt;
        }
      }
      $cur_sqpos_end = $cur_sqpos;
      $cur_apos_end  = $cur_apos;
    }
  }
  # output final insertion and/or deletion lines
  if(defined $cur_sqpos_start) {
    output_summary_insertion_string($cur_seq_name, $cur_rfpos, $cur_sqpos_start, $cur_sqpos_end, $cur_apos_start, $cur_apos_end, $cur_sqchar_str, $cur_opt_str);
    ($cur_sqpos_start, $cur_sqpos_end, $cur_apos_start, $cur_apos_end) = (undef, undef, undef, undef);
    $cur_sqchar_str = "";
  }
  if(defined $cur_rfpos_start) {
    output_summary_deletion_string($cur_seq_name, $cur_rfpos_start, $cur_rfpos_end, $cur_sqpos, $cur_apos_start, $cur_apos_end, $cur_rfchar_str, $cur_opt_str);
    ($cur_rfpos_start, $cur_rfpos_end, $cur_apos_start, $cur_apos_end) = (undef, undef, undef, undef);
    $cur_rfchar_str = "";
  }
}

#################################################################
# Subroutine: output_summary_insertion_string()
#
# Arguments:
#   $seq_name:    sequence name
#   $rfpos:       reference position
#   $sqpos_start: start sequence position 
#   $sqpos_end:   end sequence position 
#   $apos_start:  start alignment position
#   $apos_end:    end alignment position
#   $sqchar_str:  sequence string
#   $opt_str:     string for optional output, can be undef
#
#################################################################
sub output_summary_insertion_string { 
  my $sub_name = "output_summary_insertion_string";
  my $nargs_expected = 8;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }

  my ($seq_name, $rfpos, $sqpos_start, $sqpos_end,  $apos_start, $apos_end, $sqchar_str, $opt_str) = @_;

  my $sqpos_str = ($sqpos_start eq $sqpos_end) ? $sqpos_start : $sqpos_start . ".." . $sqpos_end;
  my $apos_str  = ($apos_start  eq $apos_end)  ? $apos_start  : $apos_start  . ".." . $apos_end;
  printf("%s\t%s\t%s\t%s\t%s\t%s\t%s%s\n", 
         $seq_name, $rfpos, $sqpos_str, $apos_str, "-", $sqchar_str, "insert-after-RF-position",
         (defined $opt_str) ? ("\t" . $opt_str) : "");

  return;
}

#################################################################
# Subroutine: output_summary_deletion_string()
#
# Arguments:
#   $seq_name:    sequence name
#   $rfpos_start: start reference position
#   $rfpos_end:   end reference position
#   $sqpos:       sequence position 
#   $apos_start:  start alignment position
#   $apos_end:    end alignment position
#   $rfchar_str:  reference string
#   $opt_str:     optional output string, can be undef
#
#################################################################
sub output_summary_deletion_string { 
  my $sub_name = "output_summary_deletion_string";
  my $nargs_expected = 8;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }

  my ($seq_name, $rfpos_start, $rfpos_end, $sqpos, $apos_start, $apos_end, $rfchar_str, $opt_str) = @_;

  my $rfpos_str = ($rfpos_start eq $rfpos_end) ? $rfpos_start : $rfpos_start . ".." . $rfpos_end;
  my $apos_str  = ($apos_start  eq $apos_end)  ? $apos_start  : $apos_start  . ".." . $apos_end;
  printf("%s\t%s\t%s\t%s\t%s\t%s\t%s%s\n", 
         $seq_name, $rfpos_str, $sqpos, $apos_str, $rfchar_str, "-", "deletion",
         (defined $opt_str) ? ("\t" . $opt_str) : "");

  return;
}
