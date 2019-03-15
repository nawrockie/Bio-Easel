use strict;
use warnings FATAL => 'all';
use Test::More tests => 6;

BEGIN {
    use_ok( 'Bio::Easel::MSA' ) || print "Bail out!\n";
}

# test new 
my $alnfile = "./t/data/test.sto";
my $msa = Bio::Easel::MSA->new({
   fileLocation => $alnfile, 
});
isa_ok($msa, "Bio::Easel::MSA");

# fetch 2nd sequence, aligned
my $aln = $msa->get_sqstring_aligned(1);
is($aln, "AUACACUUCGGAUG-CACC-AAA-GUGA");

# fetch 2nd sequence, unaligned
my $unaln = $msa->get_sqstring_unaligned(1);
is($unaln, "AUACACUUCGGAUGCACCAAAGUGA");

# fetch 2nd sequence, aligned and truncated
$aln = $msa->get_sqstring_aligned_and_truncated(1, 5, 24);
is($aln, "ACUUCGGAUG-CACC-AAA-");

# fetch 2nd sequence, unaligned and truncated
$aln = $msa->get_sqstring_unaligned_and_truncated(1, 5, 24);
is($aln, "ACUUCGGAUGCACCAAA");

