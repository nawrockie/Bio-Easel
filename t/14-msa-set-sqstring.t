use strict;
use warnings FATAL => 'all';
use Test::More tests => 7;

BEGIN {
    use_ok( 'Bio::Easel::MSA' ) || print "Bail out!\n";
}

#####################################################################
# We do all tests twice, once reading the MSA read in digital mode, #
# and again reading the MSA in text mode - that's what the big for  #
# loop is for.                                                      #
#####################################################################
my $alnfile = "./t/data/test.sto";
my $exp_sqstring = undef;
my $obs_sqstring = undef;

for(my $mode = 0; $mode <= 1; $mode++) { 
  # test new 
  my $msa = Bio::Easel::MSA->new({
      fileLocation => $alnfile, 
      forceText    => $mode,
  });
  isa_ok($msa, "Bio::Easel::MSA");

  #human              .AAGACUUCGGAUCUGGCG.ACA.CCC.
  #mouse              aUACACUUCGGAUG-CACC.AAA.GUGa
  #orc                .AGGUCUUC-GCACGGGCAgCCAcUUC.
  ##=GC SS_cons       .::<<<____>->>:<<-<.___.>>>.
  #//

  my $new_human = "A.AGACUUCGGAUCUGGCG.ACA.CCC.";
  my $new_mouse = "aUACACUUCGGAUG-CACCA.AA.GUGa";
  my $new_human_dig = $new_human;
  my $new_mouse_dig = $new_mouse;
  $new_human_dig =~ tr/a-z/A-Z/;
  $new_human_dig =~ s/\./\-/g;
  $new_mouse_dig =~ tr/a-z/A-Z/;
  $new_mouse_dig =~ s/\./\-/g;

  $msa->set_sqstring_aligned($new_human, 1);
  $msa->set_sqstring_aligned($new_mouse, 2);
  $exp_sqstring = ($mode == 0) ? $new_human_dig : $new_human;
  $obs_sqstring = $msa->get_sqstring_aligned(1);
  is($obs_sqstring, $exp_sqstring, "set_sqstring_aligned set aligned seq correctly seq 1 (mode $mode)");

  $exp_sqstring = ($mode == 0) ? $new_mouse_dig : $new_mouse;
  $obs_sqstring = $msa->get_sqstring_aligned(2);
  is($obs_sqstring, $exp_sqstring, "set_sqstring_aligned set aligned seq correctly seq 2 (mode $mode)");
}  
