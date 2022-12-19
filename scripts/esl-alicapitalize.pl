#!/usr/bin/env perl
# 
# esl-alicapitalize.pl: set case of aligned sequences based on RF: 
#                       nongap RF positions -> uppercase
#                          gap RF positions -> lowercase
#                       also modify SS_cons based on RF:
#                          gap RF positions -> change to '.'
#                       remaining SS_cons gets WUSSified
#                       any BPs that with left or right half that are gaps in RF
#                       are removed.
#                       Do not modify SS_cons with --keepsscons.
#                       Do not modify seqs    with --keepseqs.
#                       Do not output MSA with --checkonly, instead output
#                       '1' if alignment is unchanged (already follows RF gap conventions)
#                       '0' if alignment would be changed, and output details of changes
#                       
# EPN, Wed Dec  7 14:16:43 2022
# 
# This script uses BioEasel's MSA module and creates a new alignment
# following the rules above.

use strict;
use Getopt::Long;
use Bio::Easel::MSA;

my $version      = "0.16";
my $date        = "Dec 2022";

my $in_alifile    = "";    # name of input MSA file
my $outfile       = "";    # name of output alignment file
my $do_checkonly  = 0;     # set to '1' if --checkonly 
my $do_keepseqs   = 0;     # set to '1' if --keepseqs
my $do_keepsscons = 0;     # set to '1' if --keepsscons
my $do_perposn    = 0;     # set to '1' if --perposn

my $usage;
$usage  = "# esl-alicapitalize.pl :: set case of residues (upper/lower) and SS_cons based on RF annotation\n";
$usage .= "# Bio-Easel $version ($date)\n";
$usage .= "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
$usage .= "\n";
$usage .= "Usage: esl-alicapitalize.pl [OPTIONS] <input alignment>\n";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t--checkonly  : do not output new alignment, only '1' if alignment would have changed (with details), '0' if not\n";
$usage .= "\t\t--keepseqs   : do not modify case of residues in aligned sequences [df: do modify]\n";
$usage .= "\t\t--keepsscons : do not modify SS_cons [df: do modify]\n";
$usage .= "\t\t--perposn    : with --checkonly, output per-position comparison\n";

&GetOptions( "checkonly"  => \$do_checkonly,
             "keepseqs"   => \$do_keepseqs,
             "keepsscons" => \$do_keepsscons, 
             "perposn"    => \$do_perposn);

if(scalar(@ARGV) != 1) { die $usage; }
($in_alifile) = @ARGV;

if($do_keepseqs && $do_keepsscons) { 
  die "ERROR, --keepseqs and --keepsscons are incompatible, choose one.";
}
if($do_perposn && (! $do_checkonly)) { 
  die "ERROR, --perposn only makes sense in combination with --checkonly";
}

# validate input args
if(! -e $in_alifile) { die "ERROR $in_alifile does not exist"; }

# open file 
my $msa = Bio::Easel::MSA->new({ fileLocation => $in_alifile,
                                 forceText    => 1,
});

my $orig_msa = undef;
if($do_checkonly) { 
  $orig_msa = $msa->clone_msa();
}

if(! $do_keepseqs) { 
  $msa->capitalize_based_on_rf();
}

if(! $do_keepsscons) { 
  $msa->remove_gap_rf_basepairs(1); # 1: wussify SS_cons and change gap RF columns to '.' in SS_cons
}

if(! $do_checkonly) { 
  # default: output new MSA
  $msa->write_msa("STDOUT", "stockholm", 0);
}
else { 
  # --checkonly used
  # do not output new MSA
  # only compare new MSA with original MSA to see if it was modified
  # if it was not modified print '1' and exit

  my $was_modified = 0;
  my $details = "";
  for(my $i = 0; $i < $msa->nseq; $i++) { 
    my $orig_sqstring = $orig_msa->get_sqstring_aligned($i);
    my $new_sqstring  = $msa->get_sqstring_aligned($i);
    if($orig_sqstring ne $new_sqstring) { 
      $was_modified = 1;
      $details .= sprintf("aligned sequence %d (%s) changed:\n\tinput:  $orig_sqstring\n\toutput: $new_sqstring\n", ($i+1), $msa->get_sqname($i));
      if($do_perposn) {
        $details .= sprintf("per-position differences:\n");
        my @orig_A = split("", $orig_sqstring);
        my @new_A  = split("", $new_sqstring);
        my $len = scalar(@orig_A);
        if($len != (scalar(@new_A))) { die sprintf("ERROR length of sequence %d unexpectedly changed from %d to %d", ($i+1), $len, scalar(@new_A)); }
        for(my $j = 1; $j <= $len; $j++) { 
          if($orig_A[$j] ne $new_A[$j]) { 
            $details .= sprintf("\tposition %5d changed from %s to %s\n", $j, $orig_A[$j], $new_A[$j]);
          }
        }
        $details .= "#\n";
      }
    }
  }
  my $orig_sscons = $orig_msa->get_ss_cons();
  my $new_sscons  = $msa->get_ss_cons();
  if($orig_sscons ne $new_sscons) { 
    $was_modified = 1;
    $details .= "SS_cons changed:\n\tinput:  $orig_sscons\n\toutput: $new_sscons\n";
    if($do_perposn) { 
      $details .= "per-position differences:\n";
      my @orig_A = split("", $orig_sscons);
      my @new_A  = split("", $new_sscons);
      my $len = scalar(@orig_A);
      if($len != (scalar(@new_A))) { die sprintf("ERROR length of SS_cons unexpectedly changed from %d to %d", $len, scalar(@new_A)); }
      for(my $j = 1; $j <= $len; $j++) { 
        if($orig_A[$j] ne $new_A[$j]) { 
          $details .= sprintf("\tposition %5d changed from %s to %s\n", $j, $orig_A[$j], $new_A[$j]);
        }
      }
      $details .= "#\n";
    }
  }
  if(! $do_perposn) { 
    $details .= "# Use --perposn to see per-position changes\n";
  }
  if($was_modified) { 
    print "1\n$details";
  }
  else { 
    print "0\n";
  }
}
