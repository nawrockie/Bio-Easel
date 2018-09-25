use strict;
use warnings FATAL => 'all';
use Test::More tests => 16;


BEGIN {
    use_ok( 'Bio::Easel::SqFile' ) || print "Bail out!\n";
}

##################################################################
# We do all tests twice, once reading the sqfile read in text    #
# mode and again reading the sqfile in digital mode - that's     # 
# what the big for loop is for.                                  #
##################################################################
my $infile = "./t/data/trna-100.fa";
my ($sqfile, $sqfile2);
my $mode;
my $path;
my $seqstring;
my @AA = ();
my $tmpfile;
my $tmpsqfile;

# first test new without a forceDigital value
$sqfile = Bio::Easel::SqFile->new({
   fileLocation => $infile, 
   forceIndex   => 1, 
});
isa_ok($sqfile, "Bio::Easel::SqFile");
undef $sqfile;

# now do all tests with in both text (mode == 0) and digital (mode == 1) modes
for($mode = 0; $mode <= 1; $mode++) { 
  undef $sqfile;

  $sqfile = Bio::Easel::SqFile->new({
      fileLocation => $infile, 
      forceDigital => $mode,
      forceIndex   => 1, 
  });
  isa_ok($sqfile, "Bio::Easel::SqFile");

  # test sqfile
  $sqfile2 = $sqfile->sqfile;
  isa_ok($sqfile2, "ESL_SQFILE");

  # test path
  $path = $sqfile->path;
  is($path, "./t/data/trna-100.fa");

  # test fetch_subseq_to_fasta_string with no line length
  $seqstring = $sqfile->fetch_subseq_to_fasta_string("tRNA5-sample33", 13, 31);
  is ($seqstring, ">tRNA5-sample33/13-31\nAAGUGGCAUCGCACUUGAC\n");

  # test fetch_subseqs
  @AA = ();
  @{$AA[0]} = ("tRNA5-sample33/13-31", 13, 31, "tRNA5-sample33");
  @{$AA[1]} = ("revcomp", 31, 13, "tRNA5-sample33");
  $seqstring = $sqfile->fetch_subseqs(\@AA, -1);
  if($mode == 0) { # text mode, should be DNA 
    is ($seqstring, ">tRNA5-sample33/13-31\nAAGUGGCAUCGCACUUGAC\n>revcomp\nGTCAAGTGCGATGCCACTT\n");
  } 
  else { # digital mode, should be RNA
    is ($seqstring, ">tRNA5-sample33/13-31\nAAGUGGCAUCGCACUUGAC\n>revcomp\nGUCAAGUGCGAUGCCACUU\n");
  }

  # fetch same subseqs to a temporary file 
  $tmpfile = "t/data/tmp.fa";
  $sqfile->fetch_subseqs(\@AA, 60, $tmpfile);
  # open it 
  $tmpsqfile = Bio::Easel::SqFile->new({
     fileLocation => $tmpfile, 
     forceDigital => $mode,
     forceIndex   => 1, 
  });
  isa_ok($tmpsqfile, "Bio::Easel::SqFile");

  # now fetch these seq again, revcomp it
  @AA = ();
  @{$AA[0]} = ("revrevcomp", 19, 1, "revcomp"); 
  $seqstring = $tmpsqfile->fetch_subseqs(\@AA, -1);
  if($mode == 0) { # text mode, should be DNA 
    is ($seqstring, ">revrevcomp\nAAGTGGCATCGCACTTGAC\n");
  } 
  else { # digital mode, should be RNA
    is ($seqstring, ">revrevcomp\nAAGUGGCAUCGCACUUGAC\n");
  }

  # clean up files we just created
  unlink ($tmpfile);
  unlink ($tmpfile . ".ssi");
}

