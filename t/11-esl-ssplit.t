#! /usr/bin/perl
#
# Tests for esl-ssplit.pl: a script that uses BioEasel to split up a fasta file
# into smaller files.
#
# EPN, Thu Jan 16 09:49:03 2014
use strict;
use warnings FATAL => 'all';
use Test::More tests => 17;

BEGIN {
  use_ok( 'Bio::Easel::SqFile' ) || print "Bail out!\n";
}

my $datadir   = "./t/data/esl-ssplit";      
my $scriptdir = "./scripts";               
my $miniappdir = "./src/easel/miniapps"; 

# first run three main modes with three different input files.
my @arg1A     = ("rna-10Kb.fa",         "rna-1Mb.fa",         "aa-10k.fa");         # the three test files, command line argument for test 1
my @zarg1A    = ("sort.z.rna-10Kb.seq", "sort.z.rna-1Mb.seq", "sort.z.aa-10k.seq"); # the three test files, with seqs in random order for testing -z
my @arg2A     = ("3",                   "7",                  "100");               # command line argument for each test file
my @nfiles1A  = ("7",                   "22",                 "101");               # number of files generated for each test file for default options
my @nfiles2A  = ("3",                   "7",                  "100");               # number of files generated for each test file for -n option
my @nfiles3A  = ("3",                   "7",                  "100");               # number of files generated for each test file for -n and -r option
my @nfiles4A  = ("3",                   "7",                  "100");               # number of files generated for each test file for -n and -r and -z option
my $ntestfiles = 3;

my @unlinkA    = ();     # array of files to unlink after each test
my @reqdfilesA = (@arg1A, @zarg1A); # list of files to copy to current dir 
my $tmpdir = "tmp-esl-ssplit-dir";
if(! -e $tmpdir) { 
  run_command("mkdir $tmpdir");
}
copy_orig_files($datadir, \@reqdfilesA, $tmpdir);
push(@unlinkA, @reqdfilesA);

for(my $f = 0; $f < $ntestfiles; $f++) { 
  my $arg1  = $arg1A[$f];
  my $arg2  = $arg2A[$f];
  my $zarg1 = $zarg1A[$f];

  # test default parameters 
  run_command("$scriptdir/esl-ssplit.pl $tmpdir/$arg1 $arg2");
  # test the output by concatenating, reformatting and diff'ing against original file
  my $diff = concatenate_reformat_maybe_sort_and_diff($miniappdir, "$tmpdir/$arg1", "$tmpdir/$arg1", $nfiles1A[$f], 0); # 0: don't sort before diff
  is($diff, "", "esl-ssplit $arg1 split correctly with default parameters");

  # test -n
  run_command($scriptdir . "/esl-ssplit.pl -n $tmpdir/$arg1 $arg2");
  $diff = concatenate_reformat_maybe_sort_and_diff($miniappdir, "$tmpdir/$arg1", "$tmpdir/$arg1", $nfiles2A[$f], 0); # 0: don't sort before diff
  is($diff, "", "esl-ssplit $arg1 split correctly with -n option");

  # test -n and -r 
  run_command($scriptdir . "/esl-ssplit.pl -n -r $tmpdir/$arg1 $arg2");
  $diff = concatenate_reformat_maybe_sort_and_diff($miniappdir, "$tmpdir/$arg1", "$tmpdir/$arg1", $nfiles3A[$f], 0); # 0: don't sort before diff 
  is($diff, "", "esl-ssplit $arg1 split correctly with -n and -r options");

  # test -n and -r and -z, compare against randomly constructed file, only on first 
  run_command($scriptdir . "/esl-ssplit.pl -n -r -z $tmpdir/$arg1 $arg2");
  # note we pass in $zarg1, this is the version of the file with sequences in random order
  $diff = concatenate_reformat_maybe_sort_and_diff($miniappdir, "$tmpdir/$zarg1", "$tmpdir/$arg1", $nfiles4A[$f], 1); # 1: do sort before diff
  is($diff, "", "esl-ssplit $arg1 split correctly with -n and -r and -z options");
}

# now test other options: -d -oroot and -odir
my $arg1 = $arg1A[0];
my $arg2 = $arg2A[0];

# test -d
run_command($scriptdir . "/esl-ssplit.pl -d $tmpdir/$arg1 $arg2");
my $diff = concatenate_reformat_maybe_sort_and_diff($miniappdir, "$tmpdir/$arg1", "$tmpdir/$arg1", $nfiles1A[0], 0); # 0: don't sort before diff 
is($diff, "", "esl-ssplit $arg1 split correctly with -d");
my $ssi_exists = (-e "$tmpdir/$arg1.ssi") ? 1 : 0;
is($ssi_exists, 1, "esl-ssplit -d correctly leaves .ssi file with -d option.");
push(@unlinkA, "$tmpdir/$arg1.ssi");

# test -oroot
my $oroot = "root";
run_command($scriptdir . "/esl-ssplit.pl -oroot $oroot $tmpdir/$arg1 $arg2");
$diff = concatenate_reformat_maybe_sort_and_diff($miniappdir, "$tmpdir/$arg1", $oroot, $nfiles1A[0], 0); # 0: don't sort before diff 
is($diff, "", "esl-ssplit $arg1 split correctly with -oroot $oroot option");

# test -odir
my $odir = $tmpdir . "/odir-test";
if(! -e $odir) { 
  run_command("mkdir $odir");
}
run_command($scriptdir . "/esl-ssplit.pl -odir $odir $tmpdir/$arg1 $arg2");
$diff = concatenate_reformat_maybe_sort_and_diff($miniappdir, "$tmpdir/$arg1", "$odir/$arg1", $nfiles1A[0], 0); # 0: don't sort before diff 
is($diff, "", "esl-ssplit $arg1 split correctly with -odir $odir option");

clean_up(\@unlinkA);
rm_orig_files(\@reqdfilesA, $tmpdir);
system("rmdir $odir");
system("rmdir $tmpdir");

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
sub copy_orig_files {
  if(scalar(@_) != 3) { die "ERROR copy_orig_files entered with wrong number of input args"; }
  my ($datadir, $fileAR, $destdir) = (@_);
  foreach my $file (@{$fileAR}) { 
    if(! -e $datadir . "/" . $file) { die "ERROR, unable to copy required file $file from $datadir"; }
    run_command("cp $datadir/$file $destdir");
  }
  return;
}
###############
sub rm_orig_files {
  if(scalar(@_) != 2) { die "ERROR copy_orig_files entered with wrong number of input args"; }
  my ($fileAR, $destdir) = (@_);
  foreach my $file (@{$fileAR}) { 
    if(-e "$destdir/$file") { unlink "$destdir/$file"; }
  }
  return;
}
###############
sub concatenate_reformat_maybe_sort_and_diff { 
  if(scalar(@_) != 5) { die "ERROR concatenate_reformat_maybe_sort_and_diff entered with wrong number of input args"; }
  my ($miniappdir, $origfile, $smallfileroot, $nfiles, $do_sort) = (@_);

  my $testfile    = "test.fa";
  my $rf_testfile = "test.rf.fa";
  
  my @unlinkA = ();

  my $cmd = "cat ";
  for(my $i = 1; $i <= $nfiles; $i++) { 
    $cmd .= " $smallfileroot.$i";
    push(@unlinkA, "$smallfileroot.$i");
  }
  $cmd .= " > $testfile";
  run_command($cmd);

  # reformat with esl-reformat
  if($do_sort) { 
    $cmd = $miniappdir . "/esl-reformat fasta $testfile | grep -v \^\\> | sort > $rf_testfile"; # remove FASTA header lines, '>' character seems to get sorted differently on different systems
  }
  else { 
    $cmd = $miniappdir . "/esl-reformat fasta $testfile > $rf_testfile";
  }
  run_command($cmd);

  # compare with diff
  my $diff_file = "$testfile.diff";
  $cmd = "diff $origfile $rf_testfile > $diff_file";
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

  push(@unlinkA, ($diff_file, $testfile, $rf_testfile));
  clean_up(\@unlinkA);
  
  return $diff_output;
}
###############
