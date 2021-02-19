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

  $msa->set_sqstring_aligned($new_human, 1);
  if($mode == 0) { 
    $new_human =~ tr/a-z/A-Z/;
    $new_human =~ s/\./\-/g;
  }
  is($msa->get_sqstring_aligned(1), $new_human, "set_sqstring_aligned set aligned seq correctly seq 1 (mode $mode)");

  $msa->set_sqstring_aligned($new_mouse, 2);
  if($mode == 0) { 
    $new_mouse =~ tr/a-z/A-Z/;
    $new_mouse =~ s/\./\-/g;
  }
  is($msa->get_sqstring_aligned(2), $new_mouse, "set_sqstring_aligned set aligned seq correctly seq 2 (mode $mode)");
}  
