#!/usr/bin/env perl
# 
# esl-alitransfer.pl:  transfer GC or GF annotation from one alignment to another.
#                       
# EPN, Fri Sep 13 14:01:37 2024
# 
# This script uses BioEasel's MSA module.

use strict;
use Getopt::Long;
use Bio::Easel::MSA;

my $version     = "0.16";
my $date        = "Dec 2022";

my $in_gc = undef;
my $in_gf = undef;
my $do_force = 0;
my $gapstr = ".-~";

my $usage;
$usage  = "# esl-alitransfer.pl :: tranfer GC and/or GF annotation from one alignment to another\n";
$usage .= "# Bio-Easel $version ($date)\n";
$usage .= "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
$usage .= "\n";
$usage .= "Usage: esl-alitransfer.pl [OPTIONS] <source alignment> <destination alignment>\n";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t--gc <s>     : transfer GC annotation keys listed in comma separated <s>\n";
$usage .= "\t\t--gf <s>     : transfer GF annotation keys listed in comma separated <s>\n";
$usage .= "\t\t--force      : allow nongap characters in gap RF positions in source GC annotation\n";
$usage .= "\t\t--gapstr <s> : define gaps as any characters in string <s> [df:$gapstr]\n";

&GetOptions( "gc=s"     => \$in_gc, 
             "gf=s"     => \$in_gf,
             "force"    => \$do_force,
             "gapstr=s"); 

if(scalar(@ARGV) != 2) { die $usage; }
my ($src_alifile, $dst_alifile) = @ARGV;

if((! defined $in_gc) && (! defined $in_gf)) {
  die "ERROR, at least one of --gc or --gf must be used\n";
}

# validate input args
if(! -e $src_alifile) { die "ERROR $src_alifile does not exist"; }
if(! -e $dst_alifile) { die "ERROR $dst_alifile does not exist"; }

# open file 
my $src_msa = Bio::Easel::MSA->new({ fileLocation => $src_alifile });
my $dst_msa = Bio::Easel::MSA->new({ fileLocation => $dst_alifile });

# if --gc: make sure we have RF annotation in which case alignments need to be the same
if(defined $in_gc) { 
  if(! $src_msa->has_rf) { die "ERROR $src_alifile does not have RF annotation, it must"; }
  if(! $dst_msa->has_rf) { die "ERROR $dst_alifile does not have RF annotation, it must"; }
}

# get the list of tags we want to transfer
my @transfer_gc_tag_A = ();
my @transfer_gf_tag_A = ();
if(defined $in_gc) {
  @transfer_gc_tag_A = split(",", $in_gc);
}
if(defined $in_gf) {
  @transfer_gf_tag_A = split(",", $in_gf);
  # check that GA, NC, TC don't exist in $in_gf, these are the only 3 we can't handle
  foreach my $transfer_gf_tag (@transfer_gf_tag_A) {
    if(($transfer_gf_tag eq "GA") ||
       ($transfer_gf_tag eq "NC") ||
       ($transfer_gf_tag eq "TC")) {
      die "ERROR, this script is unable to transfer GF tags GA, NC or TC annotation\n";
    }
  }
}

if(defined $in_gc) { 
  # make maps of RF positions to alignment positions in both alignments             
  my @src_rf2a_map_A = ();
  my @src_a2rf_map_A = ();
  my @dst_rf2a_map_A = ();
  my @dst_a2rf_map_A = ();
  $src_msa->get_rf_map(\@src_rf2a_map_A, \@src_a2rf_map_A, undef);
  $dst_msa->get_rf_map(\@dst_rf2a_map_A, \@dst_a2rf_map_A, undef);
  my $src_rflen = scalar(@src_rf2a_map_A) - 1;
  my $dst_rflen = scalar(@dst_rf2a_map_A) - 1;
  if($src_rflen != $dst_rflen) {
    die "ERROR nongap RF lengths are not identical; $src_alifile: $src_rflen, $dst_alifile: $dst_rflen";
  }

  # get all the GC tags in the src_msa
  my @src_gc_tag_A = (); # array of GC tags in src_msa
  my $src_ngc = $src_msa->getGC_number();
  my $gc_idx = 0;
  for($gc_idx = 0; $gc_idx < $src_ngc; $gc_idx++) {
    push(@src_gc_tag_A, $src_msa->getGC_tag($gc_idx));
  }

  # foreach tag we want to transfer, make sure it exists in src_msa,
  # and tranfser it based on the RF annotation
  foreach my $transfer_gc_tag (@transfer_gc_tag_A) {
    my $src_aln_gc   = undef; # the GC annotation, length src_msa->alen
    my $src_unaln_gc = undef; # the GC annotation, gap RF positions removed, length src_msa->rflen
    if($transfer_gc_tag eq "RF") {
      if(! $src_msa->has_rf()) {
        die "ERROR, trying to transfer RF annotation, but MSA read from $src_alifile does not have RF annotation\n";
      }
      $src_aln_gc = $src_msa->get_rf();
    }
    elsif($transfer_gc_tag eq "SS_cons") {
      if(! $src_msa->has_ss_cons()) {
        die "ERROR, trying to transfer SS_cons annotation, but MSA read from $src_alifile does not have SS_cons annotation\n";
      }
      $src_aln_gc = $src_msa->get_ss_cons();
    }
    elsif($transfer_gc_tag eq "SA_cons") {
      if (! $src_msa->has_sa_cons()) {
        die "ERROR, trying to transfer SA_cons annotation, but MSA read from $src_alifile does not have SA_cons annotation\n";
      }
      $src_aln_gc = $src_msa->get_sa_cons();
    }
    elsif($transfer_gc_tag eq "PP_cons") {
      if(! $src_msa->has_pp_cons()) {
        die "ERROR, trying to transfer PP_cons annotation, but MSA read from $src_alifile does not have PP_cons annotation\n";
      }
      $src_aln_gc = $src_msa->get_pp_cons();
    }
    elsif($transfer_gc_tag eq "MM") {
      if(! $src_msa->has_mm()) {
        die "ERROR, trying to transfer MM annotation, but MSA read from $src_alifile does not have MM annotation\n";
      }
      $src_aln_gc = $src_msa->get_mm();
    }
    else { # not RF, SS_cons, SA_cons, PP_cons, MM
      my $found_it = 0;
      for($gc_idx = 0; $gc_idx < $src_ngc; $gc_idx++) {
        if($src_gc_tag_A[$gc_idx] eq $transfer_gc_tag) {
          $found_it = 1;
          $src_aln_gc = $src_msa->getGC_given_idx($gc_idx);
          $gc_idx = $src_ngc; # breaks loop
        }
      }
      if(! $found_it) {
        die "ERROR, trying to transfer GC $transfer_gc_tag annotation, but MSA read from $src_alifile does not have GC $transfer_gc_tag annotation\n";
      }
    }
    # if we get here, $src_unaln_gc is GC annotation we want from src_msa including gap RF positions
    my $src_unaln_gc = "";
    my @src_aln_gc_A = split("", $src_aln_gc);

    # create the unaligned GC annotation for $transfer_gc_tag
    # check to make sure there are no non-gap characters in gap RF positions, if so, exit (unless --force used)
    for(my $src_apos = 1; $src_apos <= $src_msa->alen; $src_apos++) {
      if($src_a2rf_map_A[$src_apos] == -1) { # gap in RF
        if(($src_aln_gc_A[($src_apos-1)] !~ m/[\Q$gapstr\E]/) && # not a gap in $src_aln_gc
           (! $do_force)) { # --force not used
          die "ERROR, in $transfer_gc_tag annotation from MSA read from $src_alifile position $src_apos is a gap in RF, but not in $transfer_gc_tag (use --force to ignore this)\n";
        }
      }
      else { # not a gap in RF
        $src_unaln_gc .= $src_aln_gc_A[($src_apos-1)];
      }
    }

    for(my $src_rfpos = 1; $src_rfpos <= $src_rflen; $src_rfpos++) {
      $src_unaln_gc .= $src_aln_gc_A[($src_rf2a_map_A[$src_rfpos]-1)]
    }
    my @src_unaln_gc_A = split("", $src_unaln_gc);
    my @dst_gc_A = ();
    for(my $dst_apos = 1; $dst_apos <= $dst_msa->alen; $dst_apos++) {
      if($dst_a2rf_map_A[$dst_apos] != -1) { 
        push(@dst_gc_A, $src_unaln_gc_A[($dst_a2rf_map_A[$dst_apos]-1)]);
      }
      else {
        push(@dst_gc_A, ".");
      }
    }
    $dst_msa->addGC($transfer_gc_tag, \@dst_gc_A);
  } # end of 'foreach my $transfer_gc_tag (@transfer_gc_tag_A)'
} # end of 'if(defined $in_gc)'

if(defined $in_gf) {
  # get all the GF tags in the src_msa
  my @src_gf_tag_A = (); 
  my @src_gf_value_A = (); 
  $src_msa->get_all_GF(\@src_gf_tag_A, \@src_gf_value_A);
  my $ngf = scalar(@src_gf_tag_A);
  
  # foreach tag we want to transfer, make sure it exists in src_msa, and transfer it
  # the same tag can have multiple entries
  foreach my $transfer_gf_tag (@transfer_gf_tag_A) {
    my $found_it = 0;
    for(my $f = 0; $f < $ngf; $f++) {
      if($src_gf_tag_A[$f] eq $transfer_gf_tag) {
        $dst_msa->addGF($transfer_gf_tag, $src_gf_value_A[$f]);
        $found_it = 1;
      }
    }
    if(! $found_it) { 
      die "ERROR, trying to transfer GF $transfer_gf_tag annotation, but MSA read from $src_alifile does not have GF $transfer_gf_tag annotation\n";
    }
  }
}

$dst_msa->write_msa("STDOUT", "stockholm", 0);
