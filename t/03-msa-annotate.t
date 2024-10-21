use strict;
use warnings FATAL => 'all';
use Test::More tests => 107;


BEGIN {
    use_ok( 'Bio::Easel::MSA' ) || print "Bail out!\n";
}

my ($alnfile, $line, $msa, $outfile, $has_rf, $alnfile2, $has_gc, $test_gc, $foo_gc, $tagidx, $tag, $tagnum);
my ($alnfile3, $has_sa_cons, $has_pp_cons, $has_mm, $sa_cons, $pp_cons, $mm, $set_sa_cons, $set_pp_cons, $set_mm);
my ($pp_alnfile, $grstr1, $grstr2, $grstr3, $has_gr, $test_gr);
my @gcA  = ();
my @gcA2 = ();
$alnfile = "./t/data/test.sto";
$alnfile2 = "./t/data/test.rf.sto";
$alnfile3 = "./t/data/test.gc.sto";
$pp_alnfile = "./t/data/test-pp.sto";

# do all tests twice, once in digital and once in text mode
for(my $mode = 0; $mode <= 1; $mode++) { 
   my $msa = Bio::Easel::MSA->new({
     fileLocation => $alnfile, 
     forceText    => $mode,
   });
   isa_ok($msa, "Bio::Easel::MSA");

   my $pp_msa = Bio::Easel::MSA->new({
     fileLocation => $pp_alnfile, 
     forceText    => $mode,
   });
   isa_ok($msa, "Bio::Easel::MSA");

   $outfile = "./t/data/test-msa-annot.out";

   # test getAccession and getDesc, these should return "" before we
   # set accession and desc
   my $desc = $msa->getDesc;
   is($desc, "", "getDesc correctly noted absence of desc annotation (mode $mode)");
   my $acc = $msa->getAccession;
   is($acc, "", "getAccession correctly noted absence of accession (mode $mode)");

   # test setAccession, setDesc and addGF
   $msa->setAccession("RF99999");
   $msa->setDesc("Bogus RNA family");
   $msa->addGF("BM", "cmbuild CM SEED");

   # test getAccession and getDesc again
   $desc = $msa->getDesc;
   is($desc, "Bogus RNA family", "addDesc/getDesc correctly set/obtained valid description (mode $mode)");
   $acc = $msa->getAccession;
   is($acc, "RF99999", "addAccession/getAccession correctly set/obtained valid accession (mode $mode)");

   # write out MSA and check they get printed properly
   $msa->write_msa($outfile);

   open(IN, $outfile) || die "ERROR unable to open $outfile";
   $line = <IN>;
   $line = <IN>;
   chomp $line;
   is($line, "#=GF AC RF99999", "addAccession properly added accession annotation (mode $mode)");
   $line = <IN>;
   chomp $line;
   is($line, "#=GF DE Bogus RNA family", "addAccession properly added accession annotation (mode $mode)");
   $line = <IN>;
   chomp $line;
   is($line, "#=GF BM cmbuild CM SEED", "addGF properly added GF annotation (mode $mode)");
   unlink $outfile;

   # test has_rf
   $has_rf = $msa->has_rf;
   is($has_rf, "0", "has_rf correctly noted absence of RF annotation (mode $mode)");
   undef $msa;

   $msa = Bio::Easel::MSA->new({
     fileLocation => $alnfile2, 
     forceText    => $mode,
   });

   $has_rf = $msa->has_rf;
   is($has_rf, "1", "has_rf correctly detected RF annotation (mode $mode)");

  # test addGC, getGC_given_tag and hasGC
  @gcA  = ("1","1","1","1","1", "2","2","2","2","2",  "3","3","3","3","3",  "4","4","4","4","4",  "5","5","5","5","5",   "6","6","6");
  @gcA2 = ("7","7","7","7","7", "2","2","2","2","2",  "3","3","3","3","3",  "4","4","4","4","4",  "J","J","J","J","J",   "6","6","6");

  $has_gc = $msa->hasGC("SS_cons");
  is($has_gc, "1", "hasGC correctly notes presence of SS_cons annotation");

  $has_gc = $msa->hasGC("test");
  is($has_gc, "0", "hasGC correctly notes absence of GC annotation");

  $msa->addGC("test", \@gcA);

  $has_gc = $msa->hasGC("test");
  is($has_gc, "1", "hasGC correctly notes presence of GC annotation");

  $msa->addGC("foo", \@gcA2);

  $has_gc = $msa->hasGC("foo");
  is($has_gc, "1", "hasGC correctly notes presence of GC annotation");

  $test_gc = $msa->getGC_given_tag("test");
  is($test_gc, "1111122222333334444455555666", "addGC and getGC_given_tag correctly add and get GC annotation");

  $foo_gc = $msa->getGC_given_tag("foo");
  is($foo_gc, "77777222223333344444JJJJJ666", "addGC and getGC_given_tag correctly add and get GC annotation");

  # test getGC_tagidx, getGC_tag, getGC_given_idx
  $tagidx = $msa->getGC_tagidx("test");
  is($tagidx, "0", "getGC_tagidx seems to work");

  $tagidx = $msa->getGC_tagidx("foo");
  is($tagidx, "1", "getGC_tagidx seems to work");

  $tagnum = $msa->getGC_number();  
  is($tagnum, "2", "getGC_number seems to work");

  $tag = $msa->getGC_tag(1);  
  is($tag, "foo", "getGC_tag seems to work");

  $foo_gc = $msa->getGC_given_idx(1);
  is($foo_gc, "77777222223333344444JJJJJ666", "getGC_given_idx correctly gets GC annotation");

  ################################  
  # test addGR and hasGR

  $grstr1 = ".123123123123123123.123.123.123";
  $grstr2 = "1231231231231231231231231231231";
  $grstr3 = "this-is-a-test-of-GR-annotation";

  $has_gr = $pp_msa->hasGR_given_tag_sqidx("PP", 0);
  is($has_gr, "1", "hasGR_given_tag_sqidx correctly notes presence of PP annotation for seq 1");

  $has_gr = $pp_msa->hasGR_given_tag_sqidx("SS", 0);
  is($has_gr, "0", "hasGR_given_tag_sqidx correctly notes absence of SS annotation for seq 1");

  $has_gr = $pp_msa->hasGR_given_tag_sqidx("test", 0);
  is($has_gr, "0", "hasGR_given_tag_sqidx correctly notes absence of GR annotation for seq 1");

  $pp_msa->addGR("test", 0, $grstr1);

  $has_gr = $pp_msa->hasGR_given_tag_sqidx("test", 0);
  is($has_gr, "1", "hasGR_given_tag_sqidx correctly notes presence of GR annotation for seq 1");

  $has_gr = $pp_msa->hasGR_given_tag_sqidx("test", 2);
  is($has_gr, "0", "hasGR_given_tag_sqidx correctly notes absence of GR annotation for seq 3");

  $pp_msa->addGR("foo", 0, $grstr2);

  $has_gr = $pp_msa->hasGR_given_tag_sqidx("foo", 0);
  is($has_gr, "1", "hasGR_given_tag_sqidx correctly notes presence of second GR annotation for seq 1");

  $pp_msa->addGR("fooey", 2, $grstr3);

  $has_gr = $pp_msa->hasGR_given_tag_sqidx("fooey", 2);
  is($has_gr, "1", "hasGR_given_tag_sqidx correctly notes presence of GR annotation for seq 3");

  $has_gr = $pp_msa->hasGR_given_tag_sqidx("fooey", 0);
  is($has_gr, "0", "hasGR_given_tag_sqidx correctly notes absence of GR annotation for seq 1");

  $has_gr = $pp_msa->hasGR_given_tag_sqidx("fooey", 1);
  is($has_gr, "0", "hasGR_given_tag_sqidx correctly notes absence of GR annotation for seq 2");

  $test_gr = $pp_msa->getGR_given_tag_sqidx("test", 0);
  is($test_gr, ".123123123123123123.123.123.123", "addGR and getGR_given_tag_sqidx correctly add and get GR annotation for seq 1 (test)");

  $test_gr = $pp_msa->getGR_given_tag_sqidx("foo", 0);
  is($test_gr, "1231231231231231231231231231231", "addGR and getGR_given_tag_sqidx correctly add and get GR annotation for seq 1 (foo)");

  $test_gr = $pp_msa->getGR_given_tag_sqidx("fooey", 2);
  is($test_gr, "this-is-a-test-of-GR-annotation", "addGR and getGR_given_tag_sqidx correctly add and get GR annotation for seq 3 (fooey)");

  ########
  # test getGR_tagidx, getGR_tag, hasGR_given_idx, getGR_given_idx, 
  $tagidx = $pp_msa->getGR_tagidx("test");
  is($tagidx, "0", "getGR_tagidx seems to work");

  $tagidx = $pp_msa->getGR_tagidx("foo");
  is($tagidx, "1", "getGR_tagidx seems to work");

  $tagnum = $pp_msa->getGR_number();  
  is($tagnum, "3", "getGR_number seems to work");

  $tag = $pp_msa->getGR_tag(2);  
  is($tag, "fooey", "getGR_tag seems to work");

  $has_gr = $pp_msa->hasGR_given_tagidx_sqidx(1, 0);
  is($has_gr, "1", "hasGR_given_tagidx_sqidx correctly detects presence of tag for seq 1 tag 2");

  $has_gr = $pp_msa->hasGR_given_tagidx_sqidx(2, 0);
  is($has_gr, "0", "hasGR_given_tagidx_sqidx correctly detects absence of tag for seq 1 tag 3");

  $has_gr = $pp_msa->hasGR_given_tagidx_sqidx(0, 2);
  is($has_gr, "0", "hasGR_given_tagidx_sqidx correctly detects absence of tag for seq 3 tag 1");

  $has_gr = $pp_msa->hasGR_given_tagidx_sqidx(2, 2);
  is($has_gr, "1", "hasGR_given_tagidx_sqidx correctly detects presence of tag for seq 3 tag 3");

  $test_gr = $pp_msa->getGR_given_tagidx_sqidx(1, 0);
  is($test_gr, "1231231231231231231231231231231", "getGR_given_tagidx_sqidx correctly gets GR annotation for seq 1 tag 2");

  $test_gr = $pp_msa->getGR_given_tagidx_sqidx(2, 2);
  is($test_gr, "this-is-a-test-of-GR-annotation", "getGR_given_tagidx_sqidx correctly gets GR annotation for seq 3 tag 3");

   # test {has,get,set}_{sa_cons,pp_cons,mm}
   $alnfile3 = "./t/data/test.gc.sto";
   $msa = Bio::Easel::MSA->new({
     fileLocation => $alnfile3, 
     forceText    => $mode,
   });

   # has_
   $has_sa_cons = $msa->has_sa_cons;
   is($has_sa_cons, "1", "has_sa_cons correctly noted presence of SA_cons annotation (mode $mode)");

   $has_pp_cons = $msa->has_pp_cons;
   is($has_pp_cons, "1", "has_pp_cons correctly noted presence of PP_cons annotation (mode $mode)");

   $has_mm = $msa->has_mm;
   is($has_mm, "1", "has_mm correctly noted presence of MM annotation (mode $mode)");

   # get_
   $sa_cons = $msa->get_sa_cons;
   is($sa_cons, ".abcdefghijklmonpqr.stu.vwx.", "get_sa_cons correctly fetched SA_cons");

   $pp_cons = $msa->get_pp_cons;
   is($pp_cons, ".**************8988.***.333.", "get_pp_cons correctly fetched PP_cons");

   $mm = $msa->get_mm;
   is($mm, ".111111111111111111.000.111.", "get_mm correctly fetched MM");

   # set_
   $set_sa_cons = ".zyxwvutsrqpnomlkji.hgf.edc.";
   $msa->set_sa_cons($set_sa_cons);
   $sa_cons = $msa->get_sa_cons();
   is($sa_cons, $set_sa_cons, "set_sa_cons correctly set SA_cons");

   $set_pp_cons = "0000000000000000000000000000";

   $msa->set_pp_cons($set_pp_cons);
   $pp_cons = $msa->get_pp_cons();
   is($pp_cons, $set_pp_cons, "set_pp_cons correctly set PP_cons");

   $set_mm = "0000000000010000000000000000";
   $msa->set_mm($set_mm);
   $mm = $msa->get_mm();
   is($mm, $set_mm, "set_mm_correctly set MM");

   undef $msa;
}
