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
#                       (Do not modify SS_cons with --keepsscons)
#                       (Do not modify seqs    with --keepseqs)
#
# EPN, Wed Dec  7 14:16:43 2022
# 
# This script uses BioEasel's MSA module and creates a new alignment
# following the rules above.

use strict;
use Getopt::Long;
use Bio::Easel::MSA;

my $version      = "0.15";
my $date        = "June 2021";

my $in_alifile    = "";    # name of input MSA file
my $outfile       = "";    # name of output alignment file
my $do_checkonly  = 0;     # set to '1' if --checkonly 
my $do_keepseqs   = 0;     # set to '1' if --keepseqs
my $do_keepsscons = 0;     # set to '1' if --keepsscons
my $do_notwussify = 0;     # set to '1' if --notwussify

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
$usage .= "\t\t--notwussify : from SS_cons only remove gap RF basepairs, do not convert to full WUSS after\n";

&GetOptions( "checkonly"  => \$do_checkonly,
             "keepseqs"   => \$do_keepseqs,
             "keepsscons" => \$do_keepsscons, 
             "notwussify" => \$do_notwussify);

if(scalar(@ARGV) != 1) { die $usage; }
($in_alifile) = @ARGV;

if($do_keepseqs && $do_keepsscons) { 
  die "ERROR, --keepseqs and --keepsscons are incompatible, choose one.";
}
if($do_keepsscons && $do_notwussify) { 
  die "ERROR, --keepsscons and --notwussify are incompatible, choose one.";
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
  $msa->remove_gap_rf_basepairs(! $do_notwussify);
}

if(! $do_checkonly) { 
  # output new MSA
  $msa->write_msa("STDOUT", "stockholm", 0);
}
else { 
  # compare possibly modified with original MSA to see if it was modified
  my $was_modified = 0;
  my $details = "";
  for(my $i = 0; $i < $msa->nseq; $i++) { 
    my $orig_sqstring = $orig_msa->get_sqstring_aligned($i);
    my $new_sqstring  = $msa->get_sqstring_aligned($i);
    if($orig_sqstring ne $new_sqstring) { 
      $was_modified = 1;
      $details .= sprintf("aligned sequence %d (%s) changed:\n\tinput:  $orig_sqstring\n\toutput: $new_sqstring\n", ($i+1), $msa->get_sqname($i));
    }
  }
  my $orig_sscons = $orig_msa->get_ss_cons();
  my $new_sscons  = $msa->get_ss_cons();
  if($orig_sscons ne $new_sscons) { 
    $was_modified = 1;
    $details .= "SS_cons changed:\n\tinput:  $orig_sscons\n\toutput: $new_sscons\n";
  }
  if($was_modified) { 
    print "1\n$details";
  }
  else { 
    print "0\n";
  }
}
