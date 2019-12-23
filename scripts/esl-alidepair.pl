#!/usr/bin/env perl
# 
# esl-alidepair.pl: remove some consensus basepairs from an alignment based 
#                   on posterior probabilities of aligned residues in the 
#                   paired positions.
#
# EPN, Mon Jul  7 09:03:48 2014
# 
# This script uses BioEasel's MSA module and creates a new alignment
# with >= 0 of the consensus basepairs from the original alignment
# removed. Basepairs with average posterior probabilities of less
# than 0.9 by default (changeable to <x> with -a <x>) are removed.

use strict;
use Getopt::Long;
use Bio::Easel::MSA;

my $in_alifile  = "";    # name of input MSA file
my $outfile     = "";    # name of output alignment file
my $min_avgpp   = 0.9;   # minimum average posterior probability for keeping a consensus basepair
my $use_weights = 0;     # TRUE to use sequence weights
my $do_pp       = 1;     # TRUE by default, set to 0 if --nc and --dg
my $do_nc       = 0;     # TRUE to ignore PPs and filter based on fraction of non-canonicals
my $min_fractnc = undef; # defined if --nc used
my $do_dg       = 0;     # TRUE to ignore PPs and filter based on fraction of double-gap basepairs
my $min_fractdg = undef; # defined if --dg used

my $usage;
$usage  = "# esl-alidpair.pl :: remove consensus basepairs based on alignment posterior probabilities\n";
$usage .= "# Bio-Easel 0.11 (December 2019)\n";
$usage .= "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
$usage .= "\n";
$usage .= "Usage: esl-alidepair.pl [OPTIONS] <alignment file to work on> <name of output alignment file>\n";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t-a <f>   : change minimum average posterior probability to keep to <f> [df: $min_avgpp]\n";
$usage .= "\t\t-w       : use sequence weights in the alignment file to weight counts [df: do not]\n";
$usage .= "\t\t--nc <f> : ignore PPs, remove consensus basepairs for which > <f> fraction are non-canonical\n";
$usage .= "\t\t--dg <f> : ignore PPs, remove consensus basepairs for which > <f> fraction of seqs are double gaps\n";

&GetOptions( "a=s"  => \$min_avgpp,
             "w"    => \$use_weights, 
             "nc=s" => \$min_fractnc,
             "dg=s" => \$min_fractdg);

if(scalar(@ARGV) != 2) { die $usage; }
($in_alifile, $outfile) = @ARGV;

if(defined $min_fractnc) { 
  $do_nc = 1;
  $do_pp = 0;
}
if(defined $min_fractdg) { 
  $do_dg = 1;
  $do_pp = 0;
}
if($do_nc && $do_dg) { 
  die "ERROR, --nc and --dg are incompatible, choose one. Exectute twice in succession if necessary.";
}
if((! $do_nc) && (! $do_dg)) { 
  $do_pp = 1;
}

# validate input args
if(! -e $in_alifile) { die "ERROR $in_alifile does not exist"; }

# open file 
my $msa = Bio::Easel::MSA->new({ fileLocation => $in_alifile });

# check if we have sequence weights (if we need them)
if($use_weights) { 
  if(! $msa->has_sqwgts) { die "ERROR, with -w the alignment must have sequence weights, but $in_alifile does not"; }
}

# get SS_cons and convert to a CT array.
my @ctA = $msa->get_ss_cons_ct();

my $nseq = $msa->nseq;
my $alen = $msa->alen;
my ($apos, $lpos, $rpos, $pp_lpos, $pp_rpos, $char_lpos, $char_rpos);

my @avgppA   = (); # if  do_pp: [0..apos-1] if apos is left half of bp: summed/average posterior probability at left and right half of bp
my @ncfractA = (); # if  do_nc: [0..apos-1] if apos is left half of bp: summed/fraction of seqs with noncanonical (including half gaps) at bp
my @dgfractA = (); # if  do_nc: [0..apos-1] if apos is left half of bp: fraction of seqs that are double gaps for this bp
my @ngapA    = (); # [0..apos-1] if do_pp: if apos is left half of bp: number of seqs that double gaps at left and right half of bp
                   #             else    : if apos is left half of bp: number of seqs that are a gap at both left and right half of bp
                   # NOTE: ngapA is filled the same way if --nc or --dg
# initialize
for($apos = 0; $apos < $alen; $apos++) { 
  $avgppA[$apos]   = 0.; 
  $ncfractA[$apos] = 0.;
  $dgfractA[$apos] = 0.;
  $ngapA[$apos]    = 0.;
}

my $tot_nseq = 0.; # this will probably be equal to $nseq, but maybe not if $use_weights is TRUE and we have funky weights

for(my $i = 0; $i < $nseq; $i++) { 
  my $ppstr = ($do_pp) ? $msa->get_ppstring_aligned($i) : undef;
  my @ppA   = ($do_pp) ? split("", $ppstr)              : ();
  my $sqstr = ($do_pp) ? undef : $msa->get_sqstring_aligned($i);
  my @sqA   = ($do_pp) ? ()    : split("", $sqstr);
  my $seqwt = $use_weights ? $msa->get_sqwgt($i) : 1.0;
  $tot_nseq += $seqwt;
  for($apos = 0; $apos < $alen; $apos++) { 
    $lpos = $apos+1;
    $rpos = $ctA[$lpos];
    if($rpos > $lpos) { # lpos and rpos make a basepair ($lpos < $rpos)
      if(($do_nc) || ($do_dg)) { # same steps for --nc or --dg, even though for --dg we don't care about noncanonicals
        $char_lpos = $sqA[$apos];
        $char_rpos = $sqA[($rpos-1)];
        if(($char_lpos !~ m/\w/) && ($char_rpos !~ m/\w/)) { 
          $ngapA[$apos] += $seqwt;
        }
        else { # not a double gap
          $ncfractA[$apos] += (bp_is_noncanonical($char_lpos, $char_rpos) * $seqwt); 
        }
      }
      else { 
        $pp_lpos = $ppA[$apos];
        $pp_rpos = $ppA[($rpos-1)];
        if($pp_lpos ne ".") { $avgppA[$apos] += (pp_to_fraction($pp_lpos) * $seqwt); }
        else                { $ngapA[$apos] += $seqwt; }
        if($pp_rpos ne ".") { $avgppA[$apos] += (pp_to_fraction($pp_rpos) * $seqwt); }
        else                { $ngapA[$apos] += $seqwt; }
      }
    }
  }
}

# normalize to get averages, and remove those under the min avg pp
if($do_nc) { 
  printf("#%4s  %5s  %5s  %11s  %11s  %7s\n", "lpos",  "rpos",  "fnc",   "nnondblgap",  "ndblgap",   "remove?");
  printf("#%4s  %5s  %5s  %11s  %11s  %7s\n", "----", "-----", "-----", "-----------", "-----------", "-------");
}
elsif($do_dg) { 
  printf("#%4s  %5s  %5s  %11s  %11s  %7s\n", "lpos",  "rpos",  "fdg",   "nnondblgap",  "ndblgap",   "remove?");
  printf("#%4s  %5s  %5s  %11s  %11s  %7s\n", "----", "-----", "-----", "-----------", "-----------", "-------");
}
else { 
  printf("#%4s  %5s  %5s  %9s  %9s  %7s\n", "lpos", "rpos", "avgpp", "nnongap", "ngap", "remove?");
  printf("#%4s  %5s  %5s  %9s  %9s  %7s\n", "----", "-----", "-----", "---------", "---------", "-------");
}
my $remove_str;
my @new_ssconsA = split("", $msa->get_ss_cons());
for($apos = 0; $apos < $alen; $apos++) { 
  $lpos = $apos+1;
  $rpos = $ctA[$lpos];
  if($rpos > $lpos) { # lpos and rpos make a basepair ($lpos < $rpos)
    if($do_nc) { 
      $ncfractA[$apos] /= ($tot_nseq - $ngapA[$apos]);
      if($ncfractA[$apos] > $min_fractnc) { 
        $remove_str = "yes";
        $new_ssconsA[$apos] = ".";
        $new_ssconsA[$rpos-1] = ".";
      }
      else { 
        $remove_str = "no";
      }
      printf("%5d  %5d  %5.3f  %11.1f  %11.1f  %7s\n", $lpos, $rpos, $ncfractA[$apos], $tot_nseq - $ngapA[$apos], $ngapA[$apos], $remove_str);
    }
    elsif($do_dg) { 
      $dgfractA[$apos] = $ngapA[$apos] / $tot_nseq;
      if($dgfractA[$apos] > $min_fractdg) { 
        $remove_str = "yes";
        $new_ssconsA[$apos] = ".";
        $new_ssconsA[$rpos-1] = ".";
      }
      else { 
        $remove_str = "no";
      }
      printf("%5d  %5d  %5.3f  %11.1f  %11.1f  %7s\n", $lpos, $rpos, $dgfractA[$apos], $tot_nseq - $ngapA[$apos], $ngapA[$apos], $remove_str);
    }
    else { 
      $avgppA[$apos] /= (($tot_nseq * 2.) - $ngapA[$apos]);
      if($avgppA[$apos] < $min_avgpp) { 
        $remove_str = "yes";
        $new_ssconsA[$apos] = ".";
        $new_ssconsA[$rpos-1] = ".";
      }
      else { 
        $remove_str = "no";
      }
      printf("%5d  %5d  %5.3f  %9.1f  %9.1f  %7s\n", $lpos, $rpos, $avgppA[$apos], $tot_nseq - ($ngapA[$apos] / 2.), $ngapA[$apos] / 2., $remove_str);
    }
  }
}
my $new_sscons = "";
for($apos = 0; $apos < $alen; $apos++) { $new_sscons .= $new_ssconsA[$apos]; }
$msa->set_ss_cons_wuss($new_sscons);

$msa->write_msa($outfile);

exit 0;

sub pp_to_fraction {
  my ($pp) = @_;
  if($pp eq '*') { return 0.975; }
  if($pp eq '9') { return 0.9; }
  if($pp eq '8') { return 0.8; }
  if($pp eq '7') { return 0.7; }
  if($pp eq '6') { return 0.6; }
  if($pp eq '5') { return 0.5; }
  if($pp eq '4') { return 0.4; }
  if($pp eq '3') { return 0.3; }
  if($pp eq '2') { return 0.2; }
  if($pp eq '1') { return 0.1; }
  if($pp eq '0') { return 0.025; }
  die "ERROR unexpected value $pp in pp_to_fraction"; 
}

sub bp_is_noncanonical { 
  my ($lchar, $rchar) = (@_);
  # convert to uppercase
  $lchar =~ tr/a-z/A-Z/;
  $rchar =~ tr/a-z/A-Z/;
  if(($lchar eq 'A') && ($rchar eq 'U')) { return 0.; }
  if(($lchar eq 'C') && ($rchar eq 'G')) { return 0.; }
  if(($lchar eq 'G') && ($rchar eq 'C')) { return 0.; }
  if(($lchar eq 'G') && ($rchar eq 'U')) { return 0.; }
  if(($lchar eq 'U') && ($rchar eq 'A')) { return 0.; }
  if(($lchar eq 'U') && ($rchar eq 'G')) { return 0.; }
  if(($lchar eq 'A') && ($rchar eq 'T')) { return 0.; }
  if(($lchar eq 'G') && ($rchar eq 'T')) { return 0.; }
  if(($lchar eq 'T') && ($rchar eq 'A')) { return 0.; }
  if(($lchar eq 'T') && ($rchar eq 'G')) { return 0.; }
  if(($lchar !~ m/\w/) && ($rchar !~ m/\w/)) { die "ERROR unexpectedly got double gap in bp_is_noncanonical()"; }
  return 1.0; # non-canonical, including half gaps
}
