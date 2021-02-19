use strict;
use warnings FATAL => 'all';
use Test::More tests => 31;

BEGIN {
    use_ok( 'Bio::Easel::MSA' ) || print "Bail out!\n";
}

#####################################################################
# We do all tests twice, once reading the MSA read in digital mode, #
# and again reading the MSA in text mode - that's what the big for  #
# loop is for.                                                      #
#####################################################################
my $alnfile1  = "./t/data/test.sto";
my $alnfile2  = "./t/data/test-pp-sa-ss.sto";
my $ret_val  = undef;
my $human_o  = undef;
my $human_s1 = undef;
my $sample2_sq_o = undef;
my $sample2_pp_o = undef;
my $sample2_sa_o = undef;
my $sample2_ss_o = undef;
my $sample2_sq_s1 = undef;
my $sample2_pp_s1 = undef;
my $sample2_sa_s1 = undef;
my $sample2_ss_s1 = undef;

for(my $mode = 0; $mode <= 1; $mode++) { 
  # test new 
  my $msa1 = Bio::Easel::MSA->new({
      fileLocation => $alnfile1, 
      forceText    => $mode,
  });
  isa_ok($msa1, "Bio::Easel::MSA");

  #human              .AAGACUUCGGAUCUGGCG.ACA.CCC.
  #mouse              aUACACUUCGGAUG-CACC.AAA.GUGa
  #orc                .AGGUCUUC-GCACGGGCAgCCAcUUC.
  ##=GC SS_cons       .::<<<____>->>:<<-<.___.>>>.
  #//
  #                   1234567890123456789012345678
  #                            1         2
  $human_o = ".AAGACUUCGGAUCUGGCG.ACA.CCC.";
  if($mode == 0) { $human_o =~ tr/a-z/A-Z/; $human_o =~ s/\./\-/g; } 

  is($msa1->get_sqstring_aligned(0), $human_o, "get_sqstring_aligned fetched aligned seq correctly seq 1 (mode $mode)");

  $human_s1 = ".AAGACUUCGGAUCUGGC.GACA.CCC.";
  if($mode == 0) { $human_s1 =~ tr/a-z/A-Z/; $human_s1 =~ s/\./\-/g; } 
  $ret_val = $msa1->swap_gap_and_closest_residue(0, 20, 1);
  is($ret_val, "", "swap_gap_and_closest_residue, returned without error (mode $mode)");
  is($msa1->get_sqstring_aligned(0), $human_s1, "swap_gap_and_closest_residue, corrected swapped residue/gap, test 1 (before) (mode $mode)");

  # swap back
  $ret_val = $msa1->swap_gap_and_closest_residue(0, 19, 0);
  is($ret_val, "", "swap_gap_and_closest_residue, returned without error (mode $mode)");
  is($msa1->get_sqstring_aligned(0), $human_o, "swap_gap_and_closest_residue, corrected swapped residue/gap, test 1 (before) (mode $mode)");

  #-----------------
  # test new 
  my $msa2 = Bio::Easel::MSA->new({
      fileLocation => $alnfile2, 
      forceText    => $mode,
  });
  isa_ok($msa2, "Bio::Easel::MSA");
  

  #se-sample1         UACACUUC-.....GAUG-GACC.AAA.GUC
  ##=GR se-sample1 PP ********......****.****.***.***
  ##=GR se-sample1 SS ::<<<____.....>->>-<<-<.___.>>>
  ##=GR se-sample1 SA AAABBBCCC.....DDDD.EEEF.FFF.FFF
  #se-sample2         UUAACUUUG.....G---GCUCCaAAA.---
  ##=GR se-sample2 PP *********.....*...999776444....
  ##=GR se-sample2 SS ::::<____.....>----<<-<.___.>>>
  ##=GR se-sample2 SA AAABBBCCC.....D...EEEFFFFFF....
  #se-sample3         ---------aaaccGGUGGGUCCgUUAgGAC
  ##=GR se-sample3 PP .........3333323379****99999***
  ##=GR se-sample3 SS ...................<<.<.....>>>
  ##=GR se-sample3 SA .........AAABBBCCCDDDDEEEFFFFFF
  ##=GC SS_cons       ::<<<____.....>->>-<<-<.___.>>>
  ##=GC RF            aAgaCUUCG~~~~~GAucgggCg.AcA.ccc
  #                   1234567890123456789012345678901
  #                            1         2         3
  $sample2_sq_o  = "UUAACUUUG.....G---GCUCCaAAA.---";
  if($mode == 0) { $sample2_sq_o =~ tr/a-z/A-Z/; $sample2_sq_o =~ s/\./\-/g; } 
  $sample2_pp_o  = "*********.....*...999776444....";
  $sample2_sa_o  = "AAABBBCCC.....D...EEEFFFFFF....";
  $sample2_ss_o  = "::::<____.....>----<<-<.___.>>>";

  is($msa2->get_sqstring_aligned(1), $sample2_sq_o, "get_sqstring_aligned fetched aligned seq correctly seq 2 (mode $mode)");
  is($msa2->get_ppstring_aligned(1), $sample2_pp_o, "get_ppstring_aligned fetched aligned seq correctly seq 2 (mode $mode)");
  is($msa2->get_sastring_aligned(1), $sample2_sa_o, "get_sastring_aligned fetched aligned seq correctly seq 2 (mode $mode)");
  is($msa2->get_ssstring_aligned(1), $sample2_ss_o, "get_ssstring_aligned fetched aligned seq correctly seq 2 (mode $mode)");

  $sample2_sq_s1 = "UUAACUUUG.....---GGCUCCaAAA.---";
  if($mode == 0) { $sample2_sq_s1 =~ tr/a-z/A-Z/; $sample2_sq_s1 =~ s/\./\-/g; } 
  $sample2_pp_s1 = "*********........*999776444....";
  $sample2_sa_s1 = "AAABBBCCC........DEEFFFFFFF....";
  $sample2_ss_s1 = "::::<____.....--->-<<-<.___.>>>";

  $ret_val = $msa2->swap_gap_and_closest_residue(1, 18, 1);
  is($ret_val, "", "swap_gap_and_closest_residue, returned without error (mode $mode)");
  is($msa2->get_sqstring_aligned(1), $sample2_sq_s1, "swap_gap_and_closest_residue, corrected swapped residue/gap, test 1 (before) (mode $mode)");

  # swap back
  $ret_val = $msa2->swap_gap_and_closest_residue(1, 15, 0);
  is($ret_val, "", "swap_gap_and_closest_residue, returned without error (mode $mode)");
  is($msa2->get_sqstring_aligned(1), $sample2_sq_o, "swap_gap_and_closest_residue, corrected swapped residue/gap, test 1 (before) (mode $mode)");

}  
