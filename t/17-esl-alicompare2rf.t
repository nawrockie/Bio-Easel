#! /usr/bin/perl
#
# Tests for esl-alicompare2rf.pl: a script that uses BioEasel to compare sequences in an alignment  to RF
#
# EPN, Wed Nov 20 14:04:16 2024
use strict;
use warnings FATAL => 'all';
use Test::More tests => 5;

BEGIN {
  use_ok( 'Bio::Easel::MSA')     || print "Bail out!\n";
}

my $datadir   = "./t/data/esl-alicompare2rf";
my $scriptdir = "./scripts";               

my $in_df_file   = "$datadir/RF00006.lc.stk";
my $in_norf_file = "$datadir/RF00006.norf.stk";
my @unlinkA = (); # array of files to unlink at end

my $diff = undef;
# test default parameters 
run_command("$scriptdir/esl-alicompare2rf.pl $in_df_file > out.df.out");
$diff = diff("$datadir/exp.RF00006.df.out", "out.df.out", "df.out.diff");
is($diff, "", "esl-alicompare2rf generated correct output in default mode");
push(@unlinkA, ("out.df.out", "df.out.diff"));

# test --alldel
run_command("$scriptdir/esl-alicompare2rf.pl --alldel $in_df_file > out.alldel.out");
$diff = diff("$datadir/exp.RF00006.alldel.out", "out.alldel.out", "df.alldel.diff");
is($diff, "", "esl-alicompare2rf generated correct output with --alldel");
push(@unlinkA, ("out.alldel.out", "df.alldel.diff"));
     
# test --seqrf
run_command("$scriptdir/esl-alicompare2rf.pl --seqrf AF045145.1/1-88 $in_df_file > out.seqrf.out");
$diff = diff("$datadir/exp.RF00006.seqrf.out", "out.seqrf.out", "df.seqrf.diff");
is($diff, "", "esl-alicompare2rf generated correct output with --seqrf");
push(@unlinkA, ("out.seqrf.out", "df.seqrf.diff"));
     
# test --seqrf without RF
run_command("$scriptdir/esl-alicompare2rf.pl --seqrf AF045145.1/1-88 $in_norf_file > out.seqrf.norf.out");
$diff = diff("$datadir/exp.RF00006.seqrf.norf.out", "out.seqrf.norf.out", "df.seqrf.norf.diff");
is($diff, "", "esl-alicompare2rf generated correct output with --seqrf without RF");
push(@unlinkA, ("out.seqrf.norf.out", "df.seqrf.norf.diff"));
     
     
clean_up(\@unlinkA);

exit 0;

###############
# SUBROUTINES #
###############

###############
sub run_command {
  if(scalar(@_) != 1) { die "ERROR run_command entered with wrong number of input args"; }
  my ($cmd) = (@_);
  #printf("running $cmd\n");
  system($cmd);
  if($? != 0) { die "ERROR command $cmd failed (\$? = $?)"; }
  return;
}
###############
sub clean_up {
  if(scalar(@_) != 1) { die "ERROR clean_up entered with wrong number of input args"; }
  my ($unlinkAR) = (@_);
  foreach my $file (@{$unlinkAR}) { 
#    if(-e $file) { printf("unlinking $file\n"); }
    if(-e $file) { unlink $file; }
    if(-e $file) { die "ERROR, unable to unlink $file"; }
  }
  return;
}
###############
sub diff { 
  if(scalar(@_) != 3) { die "ERROR diff() entered with wrong number of input args"; }
  my ($file1, $file2, $diff_file) = (@_);

  # compare with diff
  if(-e $diff_file) { unlink $diff_file; }
  my $cmd = "diff $file1 $file2 > $diff_file";
  run_command($cmd);
  if(! -e $diff_file) { die "ERROR diff output file $diff_file was not created"; }
  if(  -s $diff_file) { die "ERROR diff output file is not empty, script failed to correctly split file"; }
  # read in $diff_file and return it's text (should be "")
  open(IN, $diff_file) || die "ERROR unable to open $diff_file";
  my $diff_output = "";
  while(my $line = <IN>) { 
    $diff_output .= $line;
  }
  close(IN);
  
  return $diff_output;
}
###############
