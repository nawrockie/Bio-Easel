use strict;
use warnings FATAL => 'all';
use Test::More tests => 48;

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
my $exists;
my $mode; 
my $path;
my $seqstring;
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

  # test check_seq_exists
  $exists = $sqfile->check_seq_exists("tRNA5-sample33");
  is ($exists, "1");

  $exists = $sqfile->check_seq_exists("tRNA6-sample33");
  is ($exists, "0");

  # test check_subseq_exists
  $exists = $sqfile->check_subseq_exists("tRNA5-sample33", 1, 10);
  is ($exists, "1");

  $exists = $sqfile->check_subseq_exists("tRNA5-sample33", 10, 1);
  is ($exists, "1");

  $exists = $sqfile->check_subseq_exists("tRNA5-sample33", 1, 71);
  is ($exists, "1");

  $exists = $sqfile->check_subseq_exists("tRNA6-sample33", 1, 10);
  is ($exists, "0");

  $exists = $sqfile->check_subseq_exists("tRNA5-sample33", 0, 10);
  is ($exists, "0");

  $exists = $sqfile->check_subseq_exists("tRNA5-sample33", 10, 0);
  is ($exists, "0");

  $exists = $sqfile->check_subseq_exists("tRNA5-sample33", 1, 72);
  is ($exists, "0");

  # test fetch_seq_to_fasta_string 
  $seqstring = $sqfile->fetch_seq_to_fasta_string("tRNA5-sample33");
  is ($seqstring, ">tRNA5-sample33\nAUAACCACAGCGAAGUGGCAUCGCACUUGACUUCCGAUCAAGAGACCGCGGUUCGAUUCC\nGCUUGGUGAUA\n");

  # test fetch_seq_to_sqstring
  $seqstring = $sqfile->fetch_seq_to_sqstring("tRNA5-sample33");
  is ($seqstring, "AUAACCACAGCGAAGUGGCAUCGCACUUGACUUCCGAUCAAGAGACCGCGGUUCGAUUCCGCUUGGUGAUA");

  # test fetch_seq_to_fasta_string_given_ssi_number with no line length
  $seqstring = $sqfile->fetch_seq_to_fasta_string_given_ssi_number(33);
  is ($seqstring, ">tRNA5-sample39\nUCCUCCUUAACCGAAUGGUAUGGUUCCCGCCAUUCAAGCGGGCGAUCAUAUGUUCAAUCC\nAUAUAGGAGGCA\n");

  # test fetch_seq_to_fasta_string_given_ssi_number with line length of 42
  $seqstring = $sqfile->fetch_seq_to_fasta_string_given_ssi_number(33, 42);
  is ($seqstring, ">tRNA5-sample39\nUCCUCCUUAACCGAAUGGUAUGGUUCCCGCCAUUCAAGCGGG\nCGAUCAUAUGUUCAAUCCAUAUAGGAGGCA\n");

  # test fetch_seq_to_fasta_string_given_ssi_number with unlimited line length
  $seqstring = $sqfile->fetch_seq_to_fasta_string_given_ssi_number(33, -1);
  is ($seqstring, ">tRNA5-sample39\nUCCUCCUUAACCGAAUGGUAUGGUUCCCGCCAUUCAAGCGGGCGAUCAUAUGUUCAAUCCAUAUAGGAGGCA\n");

  # test fetch_seq_to_fasta_string with unlimited line length
  $seqstring = $sqfile->fetch_seq_to_fasta_string("tRNA5-sample33", -1);
  is ($seqstring, ">tRNA5-sample33\nAUAACCACAGCGAAGUGGCAUCGCACUUGACUUCCGAUCAAGAGACCGCGGUUCGAUUCCGCUUGGUGAUA\n");

  # test fetch_consecutive_seqs
  $seqstring = $sqfile->fetch_consecutive_seqs(3, "", 60);
  is ($seqstring,">tRNA5-sample34\nGUCCACAAAGCGUAAUGGUCAGCGUAGCCAACCUCAAGUUGGCAGGUCUUUGUUCGAUUC\nACAGUGUGGAC\n>tRNA5-sample35\nUCAAGGGGCGUAACUUCGGUAGCGUACCUUUCUGGCAAGAGGGAGAUUUGGGGUUCAACU\nCCCUACUUGAU\n>tRNA5-sample36\nUUGCCGAUGCGCCAGUGGGGAGGCGGACGUUCUGUCACUACGUAGGUCCGUUGUUCAAUA\nCAGUGUCGGCAAC\n");

  # create seq file
  my $tmpfile = "t/data/tmp.trna-30.fa";
  $sqfile->fetch_consecutive_seqs(30, "tRNA5-sample31", 60, $tmpfile);
  # open it 
  my $tmpsqfile = Bio::Easel::SqFile->new({
     fileLocation => $tmpfile, 
     forceDigital => $mode,
     forceIndex   => 1, 
     });
  isa_ok($tmpsqfile, "Bio::Easel::SqFile");

  # test fetch_seq_to_fasta_string with no line length
  $seqstring = $tmpsqfile->fetch_seq_to_fasta_string("tRNA5-sample33");
  is ($seqstring, ">tRNA5-sample33\nAUAACCACAGCGAAGUGGCAUCGCACUUGACUUCCGAUCAAGAGACCGCGGUUCGAUUCC\nGCUUGGUGAUA\n");

  # test fetch_consecutive_seqs
  $seqstring = $tmpsqfile->fetch_consecutive_seqs(3, "", 60);
  is ($seqstring,">tRNA5-sample34\nGUCCACAAAGCGUAAUGGUCAGCGUAGCCAACCUCAAGUUGGCAGGUCUUUGUUCGAUUC\nACAGUGUGGAC\n>tRNA5-sample35\nUCAAGGGGCGUAACUUCGGUAGCGUACCUUUCUGGCAAGAGGGAGAUUUGGGGUUCAACU\nCCCUACUUGAU\n>tRNA5-sample36\nUUGCCGAUGCGCCAGUGGGGAGGCGGACGUUCUGUCACUACGUAGGUCCGUUGUUCAAUA\nCAGUGUCGGCAAC\n");

  # close file, and then fetch again, which should open it 
  # note we get first three seqs this time
  $tmpsqfile->close_sqfile();
  $seqstring = $tmpsqfile->fetch_consecutive_seqs(3, "", 60);
  is ($seqstring, ">tRNA5-sample31\nGCUGACUUAUCGGAGAAGGCCACUAGGGGAGCUUGCCAUGCUUUCUACUCGAGCGCGAUC\nCUCGAAGUCAGCG\n>tRNA5-sample32\nUCGGCCUUGGUGUAAUGGUGUAUCACGGGAGGUUGCCGUCCUCCUAGGACCGGUUGGAUC\nCCGGUAGGCUGAC\n>tRNA5-sample33\nAUAACCACAGCGAAGUGGCAUCGCACUUGACUUCCGAUCAAGAGACCGCGGUUCGAUUCC\nGCUUGGUGAUA\n");

  # clean up files we just created
  unlink ($tmpfile);
  unlink ($tmpfile . ".ssi");
}

