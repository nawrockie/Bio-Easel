#!/usr/bin/env perl
# 
# esl-aliconsensus.pl:  add GC annotation to an alignment summarizing the
#                       conservation in the alignment.
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
$usage .= "Usage: esl-aliconsensus.pl [OPTIONS] <input alignment>\n";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t--consthr1 <x> : threshold for fraction of seqs that must be covered by consensus iupac nt [df: 0.5]\n";
$usage .= "\t\t--consthr2 <x> : threshold for making consensus iupac nt uppercase [df: 0.75]\n";
$usage .= "\t\t--consrf       : for CONS annotation, mark gap RF positions as gaps\n";
$usage .= "\t\t--nocons       : do not add GC CONS consensus sequence annotation [df: do add it]\n";
$usage .= "\t\t--noconsfract  : add GC CONSFRACT consensus sequence fraction annotation\n";
$usage .= "\t\t--gapfract     : add GC GAPFRACT fraction of seqs that are gaps annotation\n";
$usage .= "\t\t--gaprf        : for GAPFRACT annotation, mark gap RF positions as gaps\n";
$usage .= "\t\t--mis          : add GC MIS 'most informative sequence' annotation\n";
$usage .= "\t\t--misgaprf <x> : with --mis, set threshold for a gap in MIS as <x> [df: 0.5]\n"; 
$usage .= "\t\t--weights      : use sequence weights in the alignment\n";
$usage .= "\t\t--describe     : output descriptions of possible annotation and exit\n";

# set defaults
my $consthr1       = 0.5;
my $consthr2       = 0.75; 
my $do_consrf      = 0;
my $do_nocons      = 0;
my $do_noconsfract = 0;
my $do_gapfract    = 0;
my $do_gaprf       = 0;
my $do_mis         = 0;
my $misgaprf       = 0.5;
my $do_weights     = 0;
my $do_describe    = 0;

my $cmdline = "esl-aliconsensus.pl ". join(" ", @ARGV);

&GetOptions( "consthr1=s"  => \$consthr1, 
             "consthr2=s"  => \$consthr2, 
             "consrf"      => \$do_consrf, 
             "nocons"      => \$do_nocons,
             "noconsfract" => \$do_nocons,
             "gapfract"    => \$do_gapfract,
             "gaprf"       => \$do_gaprf,
             "mis"         => \$do_mis,
             "misgaprf=s"  => \$misgaprf,
             "weights"     => \$do_weights,
             "describe"    => \$do_describe );

if($do_describe) {
  # describe annotations and exit
  print("Per-column annotation descriptions:\n");
  print("CONS: most specific IUPAC nt that explains <y> > $consthr1 (changeable to <x> with \"--consthr1 <x>\")\n");
  print("      if lower case, then <y> < $consthr2 (changeable to <y> with \"--consthr2 <y>\")\n");
  print("      If --consrf and alignment has RF annotation, gap RF positions will always be gaps\n");
  print("\n");
  print("CONSFRACT: <y> value for CONS annotation, encoded as explained below\n");
  print("\n");
  print("GAPFRACT: fraction of seqs that are a gap, encoded as explained below\n");
  print("          if --gaprf and alignment has RF annotation, gap RF positions will always be gaps\n");
  print("\n");
  print("MIS: most-informative-sequence, IUPAC code that corresponds to all nt above background\n");
  print("     background calculated as the frequency of nt in alignment, across all columns.\n");
  print("     Columns in which <z> > $misgaprf (changeable to <z> with \"--misgaprf\ <z>\")\n");
  print("     sequences are gaps are annotated as '.'\n");
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

if(scalar(@ARGV) != 1) { die $usage; }
($in_alifile) = @ARGV;

if($do_nocons) { # doesn't make sense to add CONSFRACT without CONS
  $do_noconsfract = 1;
} 

# validate input args
if(! -e $in_alifile) { die "ERROR $in_alifile does not exist"; }

# open file 
my $msa = Bio::Easel::MSA->new({ fileLocation => $in_alifile,
});

# make sure MSA has RF if we need it
if($do_consrf || $do_gaprf) {
  if(! $msa->has_rf) {
    die "ERROR with --consrf and --gaprf, input alignment must have #=GC RF annotation";
  }
}

my @gc_added_A = ();
# determine and add CONS and CONSFRACT annotation 
if(! $do_nocons) { 
  my @cons_fract_A = ();
  my @cons_fract_code_A = ();
  my $cons_seq = $msa->consensus_iupac_sequence($consthr1, $consthr2, $do_consrf, $do_weights, \@cons_fract_A);
  my @cons_seq_A = split("", $cons_seq);
  #printf("$cons_seq\n");
  
  # create the @cons_fract_A
  for(my $i = 0; $i < $msa->alen; $i++) {
    if($cons_seq_A[$i] eq "-") {
      $cons_seq_A[$i]        = "."; # use '.' for RF gaps
      $cons_fract_code_A[$i] = "."; # use '.' for RF gaps
    }
    else { 
      $cons_fract_code_A[$i] = frequency_to_annotation_code($cons_fract_A[$i]);
    }
  }
  $msa->addGC("CONS", \@cons_seq_A);
  push(@gc_added_A, "CONS");
  if(! $do_noconsfract) { 
    $msa->addGC("CONSFRACT", \@cons_fract_code_A);
    push(@gc_added_A, "CONSFRACT");
  }
}

# determine and add GAPFRACT annotation
if($do_gapfract) {
  my @gap_fract_A = ();
  @gap_fract_A = $msa->pos_gap($do_weights);
  my @gap_fract_code_A = ();
  for(my $i = 0; $i < $msa->alen; $i++) {
    $gap_fract_code_A[$i] = frequency_to_annotation_code($gap_fract_A[$i]);
  }

  if($do_gaprf) {
    my $rfstr = $msa->get_rf();
    my @rf_A = split("", $rfstr);
    for(my $i = 0; $i < $msa->alen; $i++) {
      if($rf_A[$i] =~ m/[\.\-\~]/) {
        $gap_fract_code_A[$i] = ".";
      }
    }
  }
  $msa->addGC("GAPFRACT", \@gap_fract_code_A);
  push(@gc_added_A, "GAPFRACT");
}

# determine and add MIS annotation
if($do_mis) { 
  my $mis = $msa->most_informative_sequence($misgaprf, $do_weights);
  my @mis_A = split("", $mis);
  $msa->addGC("MIS", \@mis_A);
  push(@gc_added_A, "MIS");
}

# write comment explaining what annotation was added and with what cmdline
my $comment = "";
if(scalar(@gc_added_A) > 0) {
  for(my $g = 0; $g < scalar(@gc_added_A) - 1; $g++) {
    $comment .= $gc_added_A[$g] . ", ";
  }
  $comment .= $gc_added_A[(scalar(@gc_added_A)-1)] . " GC annotation added with command '$cmdline' [Bio-Easel v$version]";
}
$msa->addGF("CC", $comment);

$msa->write_msa("STDOUT", "stockholm", 0);

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
