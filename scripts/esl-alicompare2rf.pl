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

my $version      = "0.16";
my $date        = "Dec 2022";

my $in_alifile  = "";    # name of input MSA file

my $usage;
$usage  = "# esl-alicompare2rf.pl :: output differences between aligned sequences and RF\n";
$usage .= "# Bio-Easel $version ($date)\n";
$usage .= "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
$usage .= "\n";
$usage .= "Usage: esl-alicompare2rf.pl <alignment file with RF anntotation>\n";
#$usage .= "Usage: esl-alicompare2rf.pl [OPTIONS] <alignment file with RF anntotation>\n";
#$usage .= "\tOPTIONS:\n";
#$usage .= "\t\t-a <f>   : change minimum average posterior probability to keep to <f> [df: $min_avgpp]\n";
#$usage .= "\t\t-w       : use sequence weights in the alignment file to weight counts [df: do not]\n";
#$usage .= "\t\t--nc <f> : ignore PPs, remove consensus basepairs for which > <f> fraction are non-canonical\n";
#$usage .= "\t\t--dg <f> : ignore PPs, remove consensus basepairs for which > <f> fraction of seqs are double gaps\n";

#&GetOptions( "a=s"  => \$min_avgpp,
#             "w"    => \$use_weights, 
#             "nc=s" => \$min_fractnc,
#             "dg=s" => \$min_fractdg);

if(scalar(@ARGV) != 1) { die $usage; }
($in_alifile) = @ARGV;

# validate input args
if(! -e $in_alifile) { die "ERROR $in_alifile does not exist"; }

# open file 
my $msa = Bio::Easel::MSA->new({ fileLocation => $in_alifile });

# check if we have RF
if(! $msa->has_rf) { die "ERROR, alignment must have RF annotation, it does not"; }

my $alen = $msa->alen;

# get RF
my $rf_str = $msa->get_rf;
my @rf_A = split("", $rf_str);
if(scalar(@rf_A) != $alen) { 
  die "ERROR unexpected alignment length mismatch $alen != %d\n";
}

printf("%-30s  %5s  %5s  %5s  %6s  %6s  description\n", 
       "#seqname", "rfpos", "sqpos", "apos", "rfchar", "sqchar");

my $nseq = $msa->nseq; 
# for each sequence, go through each position and output differences with RF
for(my $i = 0; $i < $nseq; $i++) { 
  my $seq_name = $msa->get_sqname($i);
  my $asqstring = $msa->get_sqstring_aligned($i);
  my @asqstring_A = split("", $asqstring);
  my $rfpos = 0;
  my $sqpos = 0;
  for(my $apos = 1; $apos <= $alen; $apos++) { 
    my $sqchar = $asqstring_A[($apos-1)];
    my $rfchar = $rf_A[($apos-1)];
    my $orig_sqchar = $sqchar;
    my $orig_rfchar = $rfchar;
    $sqchar =~ tr/a-z/A-Z/; # uppercase-ize
    $rfchar =~ tr/a-z/A-Z/; # uppercase-ize
    my $sq_is_gap = ($sqchar =~ m/[A-Z]/) ? 0 : 1;
    my $rf_is_gap = ($rfchar =~ m/[A-Z]/) ? 0 : 1;
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
        $desc = "deletion";
      }
      elsif($rfchar ne $sqchar) { 
        $desc = "substitution";
      }
    }
    if(defined $desc) { 
      printf("%-30s  %5d  %5d  %5d  %6s  %6s  $desc\n", 
             $seq_name, $rfpos, $sqpos, $apos, $orig_rfchar, $orig_sqchar);
    }
  }
}
