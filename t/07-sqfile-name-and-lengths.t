use strict;
use warnings FATAL => 'all';
use Test::More tests => 18;


BEGIN {
    use_ok( 'Bio::Easel::SqFile' ) || print "Bail out!\n";
}

##################################################################
# We do all tests twice, once reading the sqfile read in text    #
# mode and again reading the sqfile in digital mode - that's     # 
# what the big for loop is for.                                  #
##################################################################
my $infile = "./t/data/trna-100.fa";
my $sqfile;
my $nseq;
my $nres;
my $L;
my $seqname;
my $mode;

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

  # test new 
  $infile = "./t/data/trna-100.fa";
  $sqfile = Bio::Easel::SqFile->new({
       fileLocation => $infile, 
       forceDigital => $mode,
       forceIndex   => 1, 
 });
  isa_ok($sqfile, "Bio::Easel::SqFile");

  # test nseq_ssi
  $nseq = $sqfile->nseq_ssi();
  is ($nseq, 100);

  # test nres_ssi
  $nres = $sqfile->nres_ssi();
  is ($nres, 7087);

  # test fetch_seq_name_and_length_given_ssi_number
  ($seqname, $L) = $sqfile->fetch_seq_name_and_length_given_ssi_number(33);
  is ($seqname, "tRNA5-sample39");
  is ($L, 72);

  # test fetch_seq_name_given_ssi_number
  $seqname = $sqfile->fetch_seq_name_given_ssi_number(40);
  is ($seqname, "tRNA5-sample45");

  # test fetch_seq_name_given_ssi_number
  $L = $sqfile->fetch_seq_length_given_ssi_number(40);
  is ($L, 73);

  # test fetch_seq_length_given_name
  my $L = $sqfile->fetch_seq_length_given_name("tRNA5-sample31");
  is ($L, 73, "fetch_seq_length_given_name() failed");
}
