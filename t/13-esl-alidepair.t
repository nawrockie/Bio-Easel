#! /usr/bin/perl
#
# Tests for esl-alidepair.pl: a script that uses BioEasel to remove basepairs
#
# EPN, Wed Dec 18 05:33:16 2019
use strict;
use warnings FATAL => 'all';
use Test::More tests => 4;

BEGIN {
  use_ok( 'Bio::Easel::SqFile' ) || print "Bail out!\n";
  use_ok( 'Bio::Easel::MSA')     || print "Bail out!\n";
}

my $datadir   = "./t/data/esl-alidepair";      
my $scriptdir = "./scripts";               

my $in_file = "$datadir/in1.stk";
my @unlinkA = (); # array of files to unlink at end

my $diff = undef;
# test default parameters 
run_command("$scriptdir/esl-alidepair.pl $in_file out1.df.stk > out1.df.stdout");
$diff = diff("$datadir/exp.out1.df.stk",    "out1.df.stk",    "out1.df.stk.diff");
is($diff, "", "esl-alidepair $in_file depaired correctly (1) with default parameters");
$diff = diff("$datadir/exp.out1.df.stdout", "out1.df.stdout", "out1.df.stdout.diff");
is($diff, "", "esl-alidepair $in_file depaired correctly (2) with default parameters");
push(@unlinkA, ("out1.df.stk", "out1.df.stdout"));

clean_up(\@unlinkA);

exit 0;

###############
# SUBROUTINES #
###############
sub run_command {
  if(scalar(@_) != 1) { die "ERROR run_command entered with wrong number of input args"; }
  my ($cmd) = (@_);
  printf("running $cmd\n");
  system($cmd);
  if($? != 0) { die "ERROR command $cmd failed (\$? = $?)"; }
  return;
}
###############
sub clean_up {
  if(scalar(@_) != 1) { die "ERROR clean_up entered with wrong number of input args"; }
  my ($unlinkAR) = (@_);
  foreach my $file (@{$unlinkAR}) { 
    if(-e $file) { printf("unlinking $file\n"); }
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

  clean_up(\@unlinkA);
  
  return $diff_output;
}
###############
