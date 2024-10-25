#!/usr/bin/env perl
# 
# esl-aliconsensus.pl:  add GC annotation to an alignment summarizing the
#                       conservation in the alignment. Various options allow
#                       the addition of other types of annotation. Use
#                       --describe to see explanation of the annotation.
#                       
# EPN, Wed Sep 11 11:48:00 2024
# 
# This script uses BioEasel's MSA module and adds GC annotation
# to the input alignment.

use strict;
use Getopt::Long;
use Bio::Easel::MSA;

my $version     = "0.16";
my $date        = "Dec 2022";

my $in_alifile    = "";    # name of input MSA file
my $outfile       = "";    # name of output alignment file
my $do_checkonly  = 0;     # set to '1' if --checkonly 
my $do_keepseqs   = 0;     # set to '1' if --keepseqs
my $do_keepsscons = 0;     # set to '1' if --keepsscons
my $do_perposn    = 0;     # set to '1' if --perposn

my $usage;
$usage  = "# esl-aliconsensus.pl :: define per column consensus annotation for an alignment\n";
$usage .= "# Bio-Easel $version ($date)\n";
$usage .= "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
$usage .= "\n";
$usage .= "Usage: esl-aliconsensus.pl [options] <input alignment>\n";
$usage .= "\n";
$usage .= "options related to RF annotation:\n";
$usage .= "\t--rf_no     : do not add RF annotation (if it does not already exist)\n";
$usage .= "\t--rf_ignore : annotate all positions [df: only annotate nongap RF positions (if RF exists)]\n";
$usage .= "\t--rf_cons   : set RF annotation as CONS annotation (even if it already exists)\n";
$usage .= "\t--rf_mis    : set RF annotation as MIS annotation (even if it already exists)\n";
$usage .= "\t--rf_x      : set RF annotation as 'x' and gaps (only if it doesn't already exists)\n";
$usage .= "\t--rf_gapthr : if no RF, set gap threshold for defining RF annotation to add as <x>\n";
$usage .= "\n";
$usage .= "options related to per-column CONS annotation:\n";
$usage .= "\t--cons_thr1 <x> : threshold for fraction of seqs that must be covered by consensus iupac nt [df: 0.5]\n";
$usage .= "\t--cons_thr2 <x> : threshold for making consensus iupac nt uppercase [df: 0.75]\n";
$usage .= "\t--cons_no       : do not add consensus sequence (CONS) annotation [df: do add it]\n";
$usage .= "\t--cons_fract    : add consensus sequence fraction (CONSFRACT) annotation\n";
$usage .= "\n";
$usage .= "options related to other per-column annotation:\n";
$usage .= "\t--info     : add information content (INFO) annotation\n";
$usage .= "\t--relent   : add relative entropy (RELENT) annotation\n";
$usage .= "\t--gapfract : add fraction of seqs that are gaps (GAPFRACT) annotation\n";
$usage .= "\t--mis      : add 'most informative sequence' (MIS) annotation\n";
$usage .= "\n";
$usage .= "options related to not annotating (skipping) positions based on gap frequency:\n";
$usage .= "\t--skip          : do skip positions that are >= <x> fraction gaps\n";
$usage .= "\t--skip_thr <x>  : set gap threshold for skipping positions to <x> [df: 0.05]\n";
$usage .= "\t--skip_char <s> : set character for skipped positions (that are not gap RF) as <s> [df: 'x']\n";
$usage .= "\n";
$usage .= "other options:\n";
$usage .= "\t--weights   : use sequence weights in the alignment\n";
$usage .= "\t--describe  : output descriptions of possible annotation and exit\n";
$usage .= "\t--data <s>  : save tabular output of data underlying annotations to file <s>\n";

# set defaults
my $opt_rf_no      = 0;
my $opt_rf_ignore  = 0;
my $opt_rf_cons    = 0;
my $opt_rf_mis     = 0;
my $opt_rf_x       = 0;
my $opt_rf_gapthr  = undef;
my $df_rf_gapthr   = 0.5;
my $opt_cons_thr1  = undef;
my $df_cons_thr1   = 0.5;
my $opt_cons_thr2  = undef;
my $df_cons_thr2   = 0.75;
my $opt_cons_no    = 0;
my $opt_cons_fract = 0;
my $opt_info       = 0;
my $opt_relent     = 0;
my $opt_gapfract   = 0;
my $opt_mis        = 0;
my $opt_skip       = 0;
my $opt_skip_thr   = undef;
my $df_skip_thr    = 0.05;
my $opt_skip_char  = undef;
my $df_skip_char   = "x";
my $opt_weights    = 0;
my $opt_describe   = 0;
my $opt_data       = undef;

my $opt_nocomment  = 0;  # secret option, added so diffs in testing would be clean (comments include file paths which cause diffs to fail)

my $cmdline = "esl-aliconsensus.pl ". join(" ", @ARGV);

&GetOptions( "rf_no"       => \$opt_rf_no,
             "rf_ignore"   => \$opt_rf_ignore,
             "rf_cons"     => \$opt_rf_cons,
             "rf_mis"      => \$opt_rf_mis,
             "rf_x"        => \$opt_rf_x,
             "rf_gapthr=s" => \$opt_rf_gapthr,
             "cons_thr1=s" => \$opt_cons_thr1,
             "cons_thr2=s" => \$opt_cons_thr2,
             "cons_no"     => \$opt_cons_no,
             "cons_fract"  => \$opt_cons_fract,
             "info"        => \$opt_info,
             "relent"      => \$opt_relent,
             "gapfract"    => \$opt_gapfract,
             "mis"         => \$opt_mis,
             "skip"        => \$opt_skip,
             "skip_thr=s"  => \$opt_skip_thr,
             "skip_char=s" => \$opt_skip_char,
             "weights"     => \$opt_weights,
             "describe"    => \$opt_describe, 
             "data=s"      => \$opt_data, 
             "nocomment"   => \$opt_nocomment);

if(scalar(@ARGV) != 1) { die $usage; }
($in_alifile) = @ARGV;

# validate input args
if(! -e $in_alifile) { die "ERROR $in_alifile does not exist"; }

# check that option values make sense
if($opt_rf_no && $opt_rf_cons) {
  die "ERROR --rf_cons does not make sense in combination with --rf_no";
}
if($opt_cons_no && $opt_rf_cons) {
  die "ERROR --cons_no does not make sense in combination with --rf_cons";
}
if($opt_rf_no && $opt_rf_mis) {
  die "ERROR --rf_mis does not make sense in combination with --rf_no";
}
if($opt_rf_no && $opt_rf_x) {
  die "ERROR --rf_x does not make sense in combination with --rf_no";
}
if($opt_rf_no && (defined $opt_rf_gapthr)) {
  die "ERROR --rf_gapthr does not make sense in combination with --rf_no";
}
if($opt_rf_ignore && $opt_rf_cons) {
  die "ERROR --rf_cons does not make sense in combination with --rf_ignore";
}
if($opt_rf_ignore && $opt_rf_mis) {
  die "ERROR --rf_mis does not make sense in combination with --rf_ignore";
}
if($opt_rf_ignore && $opt_rf_x) {
  die "ERROR --rf_x does not make sense in combination with --rf_ignore";
}
if($opt_rf_ignore && (defined $opt_rf_gapthr)) {
  die "ERROR --rf_gapthr does not make sense in combination with --rf_ignore";
}
if((! $opt_skip) && (defined $opt_skip_thr)) {
  die "ERROR --skip_thr does not make sense without --skip";
}
if((! $opt_skip) && (defined $opt_skip_char)) {
  die "ERROR --skip_char does not make sense without --skip";
}

my $rf_gapthr = (defined $opt_rf_gapthr) ? $opt_rf_gapthr : $df_rf_gapthr;
my $cons_thr1 = (defined $opt_cons_thr1) ? $opt_cons_thr1 : $df_cons_thr1;
my $cons_thr2 = (defined $opt_cons_thr2) ? $opt_cons_thr2 : $df_cons_thr2;
my $skip_thr  = (defined $opt_skip_thr)  ? $opt_skip_thr  : $df_skip_thr;
my $skip_char = (defined $opt_skip_char) ? $opt_skip_char : $df_skip_char;
# check <x> from --cons_thr1 <x> and --cons_thr2 <x> are between 0.0 and 1.0
if($rf_gapthr < 0. || $rf_gapthr > 1.0) { 
  die "ERROR with --rf_gapthr <x>, <x> must be between 0.0 and 1.0";
}
if($cons_thr1 < 0. || $cons_thr1 > 1.0) { 
  die "ERROR with --cons_thr1 <x>, <x> must be between 0.0 and 1.0";
}
if($cons_thr2 < 0. || $cons_thr2 > 1.0) { 
  die "ERROR with --cons_thr2 <x>, <x> must be between 0.0 and 1.0";
}
if($skip_thr < 0. || $skip_thr > 1.0) { 
  die "ERROR with --skip_thr <x>, <x> must be between 0.0 and 1.0";
}
if(length($skip_char) != 1) { 
  die "ERROR with --skip_char <s>, <s> must be a single character";
}
if($skip_char =~ /[\.\-\~]/) { 
  die "ERROR with --skip_char <s>, <s> can't be '.', '-', or '~'";
}
   
if($opt_describe) {
  # describe annotations and exit
  #        "\t--cons_thr1 <x> : threshold for fraction of seqs that must be covered by consensus iupac nt [df: 0.5]\n";
  print("This script adds per-column (GC) annotation to an existing input MSA.\n");
  print("\n");
  print("If the input MSA has existing RF annotation it will be kept (or rewritten if --rf_cons, --rf_mis or\n");
  print("--rf_x, but gap/nongap column definitions will not be changed).\n");
  print("\n");
  print("If the input MSA does not have RF annotation it will be added (unless --rf_no) and RF gap positions\n");
  print("will be defined as positions with > $rf_gapthr fraction gaps (changeable to <x> with --rf_gapthr <x>).\n");
  print("\n");
  print("For all added GC annotation, gap RF positions will be gaps ('.'), unless --rf_ignore used in which case\n");
  print("no gaps will appear in added GC annotation.\n");
  print("\n");
  print("For all added GC annotation, if --skip is used, positions with > $skip_thr gaps (changeable to <x> with\n");
  print("--skip_thr <x>) will be marked as \"$skip_char\" (changeable to <c> with --skip_char <c>).\n");
  print("\n");
  print("Per-column annotation descriptions:\n");
  print("\n");
  print("CONS:      Most specific IUPAC nt that explains <y> > $cons_thr1 (changeable to <x> with --cons_thr1 <x>)\n");
  print("           if lower case, then <y> < $cons_thr2 (changeable to <y> with --cons_thr2 <y>)\n");
  print("           CONS is always added unless --cons_no is used, or added as RF if --rf_cons is used.\n");
  print("\n");
  print("CONSFRACT: <y> value for CONS annotation, encoded as fractional values [0..1] as explained below.\n");
  print("           Only added if --cons_fract used.\n");
  print("\n");
  print("INFO:      Information content of each position, only available for DNA or RNA alignments.\n");
  print("           Values range from 0.0 to 2.0 bits and are encoded like fractional values explained below, only\n");
  print("           with values multiplied by 2 (e.g. '*': [1.90..2.00]).\n");
  print("\n");
  print("RELENT:    Proxy for relative entropy of each position, with background distribution set as background\n");
  print("           frequency in alignment, across all columns. Actual values range from 0.0 to 1.0, and are calc'ed\n");
  print("           as '1 - e^{-1 * relent}' where relent is the relative entropy (or Kullback-Leibler distance) for\n");
  print("           the given column. The values are encoded as fractional values [0..1] as explained below.\n");
  print("\n");
  print("GAPFRACT:  Fraction of seqs that are a gap, encoded as values [0..1] as explained below.\n");
  print("           Only added if --gapfract used.\n");
  print("\n");
  print("MIS:       Most-informative-sequence, IUPAC code that corresponds to all nt above background;\n");
  print("           background calculated as the frequency of nt in alignment, across all columns.\n");
  print("           Only added if --mis used. Or added as RF if --rf_mis used.\n");
  print("\n");
  print("Encoding of fractional values [0.0..1.0] in CONSFRACT and GAPFRACT annotation:\n");
  print("     '*': [0.95..1.00]\n");
  print("     '9': [0.85..0.95)\n");
  print("     '8': [0.75..0.85)\n");
  print("     '7': [0.65..0.75)\n");
  print("     '6': [0.55..0.65)\n");
  print("     '5': [0.45..0.55)\n");
  print("     '4': [0.35..0.45)\n");
  print("     '3': [0.25..0.35)\n");
  print("     '2': [0.15..0.25)\n");
  print("     '1': [0.05..0.15)\n");
  print("     '0': [0.00..0.05)\n");
  print("\n");
  exit 0;
}

# open input file 
my $msa = Bio::Easel::MSA->new({ fileLocation => $in_alifile });
my $alen = $msa->alen;

# open output file, if --data
my $OUT_FH = undef;
if(defined $opt_data) {
  open($OUT_FH, ">", $opt_data) || die "ERROR unable to open $opt_data from --data for writing";
  print $OUT_FH ("#colidx\trfcolidx\ttag\tnumval\tcode\tskipped?\n");
}

# are we adding RF?
my $in_has_rf = $msa->has_rf();
if(($in_has_rf) && ($opt_rf_x)) {
  die "ERROR --rf_x only allowed if msa does not already have RF annotation";
}
my $add_rf = ((! $in_has_rf) && (! $opt_rf_no)) ? 1 : 0;
# are we using RF to define gap positions an all annotation we add?
my $use_rf = (($opt_rf_ignore) || ($opt_rf_no && (! $in_has_rf))) ? 0 : 1;
# are we 'skipping' annotation of gappy columns?
my $do_skip = ($opt_skip) ? 1 : 0;

my $apos;
# determine gap freqs if we need them
my @gap_fract_A      = (); # [0..a..apos-1] frequency of gaps in position a
my @gap_fract_code_A = (); # [0..a..apos-1] code for frequency of gaps in position a
if($add_rf || $do_skip || $opt_gapfract) { # we need to know gap frequencise
  @gap_fract_A = $msa->pos_gap($opt_weights);
  for($apos = 0; $apos < $msa->alen; $apos++) {
    $gap_fract_code_A[$apos] = frequency_to_annotation_code($gap_fract_A[$apos]);
  }
}

my @i_am_rf_A = (); # [0..a..alen-1]: '1' if a is a nongap RF position, '0' if gap

# if we are enforcing RF: update gap_fract_code_A to have gaps at gap RF positions
# we are enforcing RF if ($in_has_rf && (! $opt_rf_ignore)) OR
#                     if ($add_rf) 
if($in_has_rf && $use_rf) { 
  my $gap_str  = ".-~"; 
  my $rfseq = $msa->get_rf();
  my @rfseq_A = split("", $rfseq);
  for($apos = 0; $apos < $alen; $apos++) {
    if($rfseq_A[$apos] =~ /[\Q$gap_str\E]/) { 
      $gap_fract_code_A[$apos] = ".";
      $i_am_rf_A[$apos] = 0;
    }
    else {  
      $gap_fract_code_A[$apos] = frequency_to_annotation_code($gap_fract_A[$apos]);
      $i_am_rf_A[$apos] = 1;
    }
  }
}
elsif($add_rf) {
  my $rfseq = "";
  for($apos = 0; $apos < $alen; $apos++) {
    if($gap_fract_A[$apos] > $rf_gapthr) { 
      $rfseq .= ".";
      $gap_fract_code_A[$apos] = ".";
      $i_am_rf_A[$apos] = 0;
    }
    else {  
      $rfseq .= "x";
      $gap_fract_code_A[$apos] = frequency_to_annotation_code($gap_fract_A[$apos]);
      $i_am_rf_A[$apos] = 1;
    }
  }
  $msa->set_rf($rfseq);
}


# determine CONS and CONSFRACT annotation 
my @cons_seq_A = ();
my @cons_fract_A = ();
my @cons_fract_code_A = ();
my $cons_seq = $msa->consensus_iupac_sequence($cons_thr1, $cons_thr2, $use_rf, $opt_weights, \@cons_fract_A);
@cons_seq_A = split("", $cons_seq);
# create the @cons_fract_A
for($apos = 0; $apos < $msa->alen; $apos++) {
  if($cons_seq_A[$apos] eq "-") {
    $cons_seq_A[$apos]        = "."; # use '.' for RF gaps
    $cons_fract_code_A[$apos] = "."; # use '.' for RF gaps
  }
  else { 
    #printf("apos: $apos do_skip: $do_skip gap_fract_A[$apos] $gap_fract_A[$apos] opt_skip_thr $skip_thr\n");
    if(($do_skip) && ($gap_fract_A[$apos] >= $skip_thr)) {
      $cons_seq_A[$apos] = $skip_char;
      # use gap fraction not cons_fract for determining cons_fract code
      $cons_fract_code_A[$apos] = frequency_to_annotation_code($gap_fract_A[$apos]);
    }
    else { 
      $cons_fract_code_A[$apos] = frequency_to_annotation_code($cons_fract_A[$apos]);
    }
  }
}

# determine INFO annotation
my @info_code_A = ();
my @info_A = ();
if($opt_info) { 
  @info_A = $msa->pos_infocontent($opt_weights);
  for($apos = 0; $apos < $alen; $apos++) {
    if(($use_rf) && (! $i_am_rf_A[$apos])) {
      $info_code_A[$apos] = ".";
    }
    elsif(($do_skip) && ($gap_fract_A[$apos] >= $skip_thr)) {
      $info_code_A[$apos] = $skip_char;
    }
    else { 
      $info_code_A[$apos] = frequency_to_annotation_code($info_A[$apos] / 2.);
    }
  }
}

# determine RELENT annotation
my @relent_code_A = ();
my @relent_A = ();
if($opt_relent) { 
  @relent_A = $msa->pos_relentropy($opt_weights, 1, 0, undef);
  my @squashed_relent_A = (); # relative entropy values converted to a value between 0 and 1
  for($apos = 0; $apos < $alen; $apos++) {
    $squashed_relent_A[$apos] = 1 - exp(-1 * $relent_A[$apos]);
    if(($use_rf) && (! $i_am_rf_A[$apos])) {
      $relent_code_A[$apos] = ".";
    }
    elsif(($do_skip) && ($gap_fract_A[$apos] >= $skip_thr)) {
      $relent_code_A[$apos] = $skip_char;
    }
    else { 
      $relent_code_A[$apos] = frequency_to_annotation_code($squashed_relent_A[$apos]);
    }
  }
}

# determine MIS annotation
my @mis_A = ();
if($opt_mis || $opt_rf_mis) { 
  my $mis = $msa->most_informative_sequence(0., $opt_weights);
  @mis_A = split("", $mis);
  for($apos = 0; $apos < $alen; $apos++) {
    if(($use_rf) && (! $i_am_rf_A[$apos])) {
      $mis_A[$apos] = ".";
    }
    elsif(($do_skip) && ($gap_fract_A[$apos] >= $skip_thr)) {
      $mis_A[$apos] = $skip_char;
    }
    # else leave $mis_A[$apos] as is
  }
}

# add all annotation, including possibly updating RF
my @gc_added_A = ();
my $added_x_as_rf = 0;
my $added_mis_as_rf = 0;
my $added_cons_as_rf = 0;
# add RF annotation if --rf_no not used and
# if --rf_cons, --rf_mis or --rf_x were used
if((! $opt_rf_no) && (! $opt_rf_ignore)) {
  if((! $in_has_rf) || $opt_rf_cons || $opt_rf_mis || $opt_rf_x) {
    if($opt_rf_x) {
      ; # do nothing, we already added the 'x' based RF at the beginning of the script
      $added_x_as_rf = 1;
    }
    elsif($opt_rf_mis) {
      $msa->set_rf(join("", @mis_A));
      $added_mis_as_rf = 1;
    }
    else { # either $opt_rf_cons is true of (! $in_has_rf) (or both), set RF as CONS
      $msa->set_rf(join("", @cons_seq_A));
      $added_cons_as_rf = 1;
    }
  }
}

# add CONS annotation unless --cons_no or it was added as RF above
if((! $opt_cons_no) && (! $added_cons_as_rf)) {
    $msa->addGC("CONS", \@cons_seq_A);
    push(@gc_added_A, "CONS");
}  
# add CONSFRACT if --cons_fract
if($opt_cons_fract) { 
  $msa->addGC("CONSFRACT", \@cons_fract_code_A);
  push(@gc_added_A, "CONSFRACT");
}  
# add INFO annotation if --info
if($opt_info) {
    $msa->addGC("INFO", \@info_code_A);
    push(@gc_added_A, "INFO");
}  
# add RELENT annotation if --relent
if($opt_relent) {
    $msa->addGC("RELENT", \@relent_code_A);
    push(@gc_added_A, "RELENT");
}  
# add GAPFRACT annotation if --gapfract
if($opt_gapfract) {
    $msa->addGC("GAPFRACT", \@gap_fract_code_A);
    push(@gc_added_A, "GAPFRACT");
}  
# add MIS annotation if --mis
if($opt_mis) {
    $msa->addGC("MIS", \@mis_A);
    push(@gc_added_A, "MIS");
}  

# write comments explaining what annotation was added and with what cmdline
my $comment_line1 = "";
my $comment_line2 = "";
my $added_rf = ($added_x_as_rf || $added_mis_as_rf || $added_cons_as_rf) ? 1 : 0;
if((scalar(@gc_added_A == 0)) && (! $added_rf)) {
  die "ERROR no GC annotation added, probably due to strange choice of options.";
}
   
if(scalar(@gc_added_A) > 0) { 
  for(my $g = 0; $g < scalar(@gc_added_A) - 1; $g++) {
    $comment_line1 .= $gc_added_A[$g] . ", ";
  }
  $comment_line1 .= $gc_added_A[(scalar(@gc_added_A)-1)] . " GC annotation added";
  if($added_rf) { 
    $comment_line1 .= " and ";
  }
}
if($added_rf) { 
  if($in_has_rf) { 
    $comment_line1 .= "RF annotation redefined as ";
  }
  else {
    $comment_line1 .= "RF annotation defined as ";
  }
  if($added_x_as_rf) {
    $comment_line1 .= "\'x/.\'.";
    if(! $in_has_rf) { 
      $comment_line1 .= ", with gap positions defined based on fraction of gaps";
    }
  }
  if($added_mis_as_rf) {
    $comment_line1 .= "most-informative-sequence";
    if(! $in_has_rf) { 
      $comment_line1 .= ", with gap positions defined based on fraction of gaps";
    }
  }
  if($added_cons_as_rf) {
    $comment_line1 .= "CONS (consensus) annotation";
    if(! $in_has_rf) { 
      $comment_line1 .= ", with gap positions defined based on fraction of gaps";
    }
  }
}
# add commets to msa with command and list of GC tags added
$comment_line1 .= " with command:";
$comment_line2 = "'$cmdline' [Bio-Easel v$version]";
if(! $opt_nocomment) { 
  $msa->addGF("CC", $comment_line1);
  $msa->addGF("CC", $comment_line2);
}

# output data
if(defined $opt_data) {
  # GAPFRACT
  if((scalar(@gap_fract_A) > 0) && (scalar(@gap_fract_code_A) > 0)) {
    output_to_file($OUT_FH, "GAPFRACT", \@gap_fract_A, \@gap_fract_code_A, \@i_am_rf_A, $do_skip, $skip_thr, \@gap_fract_A);
  }
  # CONS 
  if(scalar(@cons_seq_A)) { 
    output_to_file($OUT_FH, "CONS", undef, \@cons_seq_A, \@i_am_rf_A, $do_skip, $skip_thr, \@gap_fract_A);
  }
  # CONSFRACT
  if((scalar(@cons_fract_A) > 0) && (scalar(@cons_fract_code_A) > 0)) {
    output_to_file($OUT_FH, "CONSFRACT", \@cons_fract_A, \@cons_fract_code_A, \@i_am_rf_A, $do_skip, $skip_thr, \@gap_fract_A);
  }
  # INFO
  if((scalar(@info_A) > 0) && (scalar(@info_code_A) > 0)) {
    output_to_file($OUT_FH, "INFO", \@info_A, \@info_code_A, \@i_am_rf_A, $do_skip, $skip_thr, \@gap_fract_A);
  }
  # RELENT
  if((scalar(@relent_A) > 0) && (scalar(@relent_code_A) > 0)) {
    output_to_file($OUT_FH, "RELENT", \@relent_A, \@relent_code_A, \@i_am_rf_A, $do_skip, $skip_thr, \@gap_fract_A);
  }
  # MIS
  if(scalar(@mis_A) > 0) {
    output_to_file($OUT_FH, "MIS", undef, \@mis_A, \@i_am_rf_A, $do_skip, $skip_thr, \@gap_fract_A);
  }
  close(OUT);
}

# output the msa
$msa->write_msa("STDOUT", "stockholm", 0);


########################
# subroutines
sub frequency_to_annotation_code {
  my ($frequency) = (@_);
  if   ($frequency >= 0.95) { return "*"; }
  elsif($frequency >= 0.85) { return "9"; }
  elsif($frequency >= 0.75) { return "8"; }
  elsif($frequency >= 0.65) { return "7"; }
  elsif($frequency >= 0.55) { return "6"; }
  elsif($frequency >= 0.45) { return "5"; }
  elsif($frequency >= 0.35) { return "4"; }
  elsif($frequency >= 0.25) { return "3"; }
  elsif($frequency >= 0.15) { return "2"; }
  elsif($frequency >= 0.05) { return "1"; }
  else                      { return "0"; }
}

########################
# subroutines
sub output_to_file { 
  my ($FH, $code, $num_AR, $val_AR, $i_am_rf_AR, $do_skip, $skip_thr, $gap_AR) = (@_);

  my $rfpos = 0;
  my @alen = scalar(@{$val_AR});
  my $i_am_rf_valid = ((defined $i_am_rf_AR) && (scalar(@{$i_am_rf_AR}) > 0)) ? 1 : 0;
  my $num_valid     = ((defined $num_AR)     && (scalar(@{$num_AR})     > 0)) ? 1 : 0;
  for(my $apos = 0; $apos < $alen; $apos++) { 
    my $rfcolidx = "-";
    if($i_am_rf_valid) {
      if($i_am_rf_AR->[$apos]) {
        $rfpos++;
        $rfcolidx = $rfpos;
      }
    }
    my $skipped  = "-";
    if($do_skip) {
      $skipped = ($gap_AR->[$apos] >= $skip_thr) ? "yes" : "no";
    }
    printf $FH ("%d\t\%s\t%s\t%s\t%s\t%s\n",
                ($apos+1),
                $rfcolidx,
                $code,
                ($num_valid) ? sprintf("%.6f", $num_AR->[$apos]) : "-",
                $val_AR->[$apos],
                $skipped);
  }
}
