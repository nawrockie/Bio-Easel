#! /usr/bin/perl
#
# Tests for esl-aliconsensus.pl: a script that uses BioEasel to add GC annotation to a MSA
#
# EPN, Fri Oct 25 13:38:54 2024
use strict;
use warnings FATAL => 'all';
use Test::More tests => 25;

BEGIN {
  use_ok( 'Bio::Easel::MSA')     || print "Bail out!\n";
}

my $datadir   = "./t/data/esl-aliconsensus";      
my $scriptdir = "./scripts";               

my $in_df_file   = "$datadir/RF00006.stk";
my $in_norf_file = "$datadir/RF00006.norf.stk";
my @unlinkA = (); # array of files to unlink at end

# run all tests twice, once without --weights, and once with --weights
#for(my $i = 0; $i < 2; $i++) { 

my $diff = undef;
# test default parameters 
run_command("$scriptdir/esl-aliconsensus.pl --nocomment $in_df_file > out.df.stk");
$diff = diff("$datadir/exp.df.stk", "out.df.stk", "df.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with default parameters");
push(@unlinkA, ("out.df.stk", "df.stk.diff"));
     
# test --rf_no
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --rf_no $in_norf_file > out.rfno.stk");
$diff = diff("$datadir/exp.rfno.stk", "out.rfno.stk", "rfno.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --rf_no");
push(@unlinkA, ("out.rfno.stk", "rfno.stk.diff"));

# test --rf_ignore
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --rf_ignore $in_df_file > out.rfignore.stk");
$diff = diff("$datadir/exp.rfignore.stk", "out.rfignore.stk", "rfignore.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --rf_ignore");
push(@unlinkA, ("out.rfignore.stk", "rfignore.stk.diff"));

# test --rf_cons --rf_gapthr 0.1
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --gapfract --rf_cons --rf_gapthr 0.1 $in_norf_file > out.rfcons.stk");
$diff = diff("$datadir/exp.rfcons.stk", "out.rfcons.stk", "rfcons.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --rf_cons --rf_gapthr");
push(@unlinkA, ("out.rfcons.stk", "rfcons.stk.diff"));

# test --rf_mis
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --rf_mis $in_df_file > out.rfmis.stk");
$diff = diff("$datadir/exp.rfmis.stk", "out.rfmis.stk", "rfmis.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --rf_mis");
push(@unlinkA, ("out.rfmis.stk", "rfmis.stk.diff"));

# test --rf_x
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --rf_x $in_norf_file > out.rfx.stk");
$diff = diff("$datadir/exp.rfx.stk", "out.rfx.stk", "rfx.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --rf_x");
push(@unlinkA, ("out.rfx.stk", "rfx.stk.diff"));

# test --cons_thr1 
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --cons_thr1 0.75 $in_df_file > out.consthr1.stk");
$diff = diff("$datadir/exp.consthr1.stk", "out.consthr1.stk", "consthr1.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --cons_thr1");
push(@unlinkA, ("out.consthr1.stk", "consthr1.stk.diff"));
     
# test --cons_thr2
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --cons_thr2 0.95 $in_df_file > out.consthr2.stk");
$diff = diff("$datadir/exp.consthr2.stk", "out.consthr2.stk", "consthr2.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --cons_thr2");
push(@unlinkA, ("out.consthr2.stk", "consthr2.stk.diff"));
     
# test --cons_no
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --cons_no $in_norf_file > out.consno.stk");
$diff = diff("$datadir/exp.consno.stk", "out.consno.stk", "consno.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --cons_no");
push(@unlinkA, ("out.consno.stk", "consno.stk.diff"));
     
# test --cons_fract
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --cons_fract $in_df_file > out.consfract.stk");
$diff = diff("$datadir/exp.consfract.stk", "out.consfract.stk", "consfract.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --cons_fract");
push(@unlinkA, ("out.consfract.stk", "consfract.stk.diff"));
     
# test --info
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --info $in_df_file > out.info.stk");
$diff = diff("$datadir/exp.info.stk", "out.info.stk", "info.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --info");
push(@unlinkA, ("out.info.stk", "info.stk.diff"));
     
# test --relent
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --relent $in_df_file > out.relent.stk");
$diff = diff("$datadir/exp.relent.stk", "out.relent.stk", "relent.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --relent");
push(@unlinkA, ("out.relent.stk", "relent.stk.diff"));
     
# test --gapfract
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --gapfract $in_df_file > out.gapfract.stk");
$diff = diff("$datadir/exp.gapfract.stk", "out.gapfract.stk", "gapfract.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --gapfract");
push(@unlinkA, ("out.gapfract.stk", "gapfract.stk.diff"));
     
# test --mis
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --mis $in_df_file > out.mis.stk");
$diff = diff("$datadir/exp.mis.stk", "out.mis.stk", "mis.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --mis");
push(@unlinkA, ("out.mis.stk", "mis.stk.diff"));
     
# test --skip
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --skip $in_norf_file > out.skip.stk");
$diff = diff("$datadir/exp.skip.stk", "out.skip.stk", "skip.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --skip");
push(@unlinkA, ("out.skip.stk", "skip.stk.diff"));

# test --skip_thr
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --skip --skip_thr 0.001 $in_df_file > out.skipthr.stk");
$diff = diff("$datadir/exp.skipthr.stk", "out.skipthr.stk", "skipthr.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --skip_thr");
push(@unlinkA, ("out.skipthr.stk", "skipthr.stk.diff"));
     
# test --skip_char
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --skip --skip_char ! $in_df_file > out.skipchar.stk");
$diff = diff("$datadir/exp.skipchar.stk", "out.skipchar.stk", "skipchar.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --skip_char");
push(@unlinkA, ("out.skipchar.stk", "skipchar.stk.diff"));
     
# test --weights
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --weights $in_norf_file > out.weights.stk");
$diff = diff("$datadir/exp.weights.stk", "out.weights.stk", "weights.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --weights");
push(@unlinkA, ("out.weights.stk", "weights.stk.diff"));
     
# test --describe
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --describe $in_norf_file > out.describe.out");
$diff = diff("$datadir/exp.describe.out", "out.describe.out", "describe.out.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --describe");
push(@unlinkA, ("out.describe.out", "describe.out.diff"));
     
# test --data
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --data out.data $in_df_file > out.data.stk");
$diff = diff("$datadir/exp.data", "out.data", "data.diff");
is($diff, "", "esl-aliconsensus output data correct with --data");
push(@unlinkA, ("out.data.stk", "data.diff", "out.data"));
     
# test a bunch of options with default model
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --info --relent --mis --gapfract --cons_fract --data out.alldf.data $in_df_file > out.alldf.stk");
$diff = diff("$datadir/exp.alldf.stk", "out.alldf.stk", "alldf.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --info --relent --mis --gapfract --cons_fract on aln with RF annotation");
$diff = diff("$datadir/exp.alldf.data", "out.alldf.data", "alldf.data.diff");
is($diff, "", "esl-aliconsensus output data correct with --data and --info --relent --mis --gapfract --cons_fract on aln with RF annotation");
push(@unlinkA, ("out.alldf.data", "alldf.stk.diff", "alldf.data.diff", "out.alldf.stk"));
     
# test a bunch of options with non-rf model
run_command("$scriptdir/esl-aliconsensus.pl --nocomment --info --relent --mis --gapfract --cons_fract --rf_no --data out.allnorf.data $in_norf_file > out.allnorf.stk");
$diff = diff("$datadir/exp.allnorf.stk", "out.allnorf.stk", "allnorf.stk.diff");
is($diff, "", "esl-aliconsensus correctly added annotation with --info --relent --mis --gapfract --cons_fract on aln without RF annotation");
$diff = diff("$datadir/exp.allnorf.data", "out.allnorf.data", "allnorf.data.diff");
is($diff, "", "esl-aliconsensus output data correct with --data and --info --relent --mis --gapfract --cons_fract on aln without RF annotation");
push(@unlinkA, ("out.allnorf.data", "allnorf.stk.diff", "allnorf.data.diff", "out.allnorf.stk"));
     
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
