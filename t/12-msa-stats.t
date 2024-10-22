use strict;
use warnings FATAL => 'all';
use Test::More tests => 65;

BEGIN {
    use_ok( 'Bio::Easel::MSA' ) || print "Bail out!\n";
}

###########################################################
# The 'stats' functions only work on digitized alignments #
# so we force digital mode for these tests.               #
###########################################################
my $alnfile         = "./t/data/test.sto";
my $rfamfile_allgap = "./t/data/RF00014-seed-allgap.sto";

my $msa1 = Bio::Easel::MSA->new({
    fileLocation => $alnfile, 
});
isa_ok($msa1, "Bio::Easel::MSA");

my $mis = $msa1->most_informative_sequence(0.5, 0);
is($mis, "-WRSWCUUCGGMWSKSRCV-MMA-BYS-", "calculate_most_informative_sequence() worked.");

my @cons_fract_A = ();
my $cons_seq = $msa1->consensus_iupac_sequence(0.5, 0.75, 0, 0, \@cons_fract_A);
is($cons_seq, "AaagaCUUCGGaucgggCmGacACsucA", "consensus_iupac_sequence() sequence construction worked.");
my $cons_fract_val = sprintf("%.3f", $cons_fract_A[0]);
is($cons_fract_val, "1.000", "consensus_iupac_sequence() sequence fraction worked (pos 1).");
$cons_fract_val = sprintf("%.3f", $cons_fract_A[1]);
is($cons_fract_val, "0.667", "consensus_iupac_sequence() sequence fraction worked (pos 2).");

my @gapA = $msa1->pos_gap();
is(int(($gapA[0] * 100) + 0.5), 67,  "calculate_pos_gap() seems to work (pos 1)");
is(int(($gapA[1] * 100) + 0.5), 0,   "calculate_pos_gap() seems to work (pos 2)");
is(int(($gapA[9] * 100) + 0.5), 33,  "calculate_pos_gap() seems to work (pos 10)");

my @fcbpA = $msa1->pos_fcbp();
is(int(($fcbpA[2] * 100) + 0.5), 0,   "calculate_pos_fcbp() seems to work (pos 3)");
is(int(($fcbpA[3] * 100) + 0.5), 100, "calculate_pos_fcbp() seems to work (pos 4)");
is(int(($fcbpA[4] * 100) + 0.5), 100, "calculate_pos_fcbp() seems to work (pos 5)");
is(int(($fcbpA[5] * 100) + 0.5), 100, "calculate_pos_fcbp() seems to work (pos 6)");
is(int(($fcbpA[6] * 100) + 0.5), 0,   "calculate_pos_fcbp() seems to work (pos 7)");

my @covA = $msa1->pos_covariation();
is(int(($covA[2] * 100) + 0.5), 0,   "calculate_pos_covariation() seems to work (pos 3)");
is(int(($covA[3] * 100) + 0.5), 133, "calculate_pos_covariation() seems to work (pos 4)");
is(int(($covA[4] * 100) + 0.5), 133, "calculate_pos_covariation() seems to work (pos 5)");
is(int(($covA[5] * 100) + 0.5), 0,   "calculate_pos_covariation() seems to work (pos 6)");
is(int(($covA[6] * 100) + 0.5), 0,   "calculate_pos_covariation() seems to work (pos 7)");

undef $msa1;
  
$msa1 = Bio::Easel::MSA->new({
    fileLocation => $rfamfile_allgap, 
});
isa_ok($msa1, "Bio::Easel::MSA");
my @entA = $msa1->pos_entropy();
is(int(($entA[1] * 100) + 0.5), 200, "calculate_pos_ent() seems to work (pos 2)");
is(int(($entA[2] * 100) + 0.5), 137, "calculate_pos_ent() seems to work (pos 3)");
is(int(($entA[3] * 100) + 0.5), 0,   "calculate_pos_ent() seems to work (pos 4)");
is(int(($entA[4] * 100) + 0.5), 72,  "calculate_pos_ent() seems to work (pos 5)");
is(int(($entA[5] * 100) + 0.5), 0,   "calculate_pos_ent() seems to work (pos 6)");

# $use_weights, $gaps_as_miss, $use_uniform_bg, $bgcounts_AR
my @relentA = $msa1->pos_relentropy(0, 1, 0, undef);
is(int(($relentA[1] * 100) + 0.5), 0,   "calculate_pos_relent() 1 seems to work (pos 2)");
is(int(($relentA[2] * 100) + 0.5), 64,  "calculate_pos_relent() 1 seems to work (pos 3)");
is(int(($relentA[3] * 100) + 0.5), 203, "calculate_pos_relent() 1 seems to work (pos 4)");
is(int(($relentA[4] * 100) + 0.5), 143, "calculate_pos_relent() 1 seems to work (pos 5)");
is(int(($relentA[5] * 100) + 0.5), 203, "calculate_pos_relent() 1 seems to work (pos 6)");

# if use_uniform_bg is 1, gaps_as_miss value is irrelevant
@relentA = $msa1->pos_relentropy(0, 1, 1, undef);
is(int(($relentA[1] * 100) + 0.5), 0,   "calculate_pos_relent() 2 seems to work (pos 2)");
is(int(($relentA[2] * 100) + 0.5), 63,  "calculate_pos_relent() 2 seems to work (pos 3)");
is(int(($relentA[3] * 100) + 0.5), 200, "calculate_pos_relent() 2 seems to work (pos 4)");
is(int(($relentA[4] * 100) + 0.5), 128, "calculate_pos_relent() 2 seems to work (pos 5)");
is(int(($relentA[5] * 100) + 0.5), 200, "calculate_pos_relent() 2 seems to work (pos 6)");

@relentA = $msa1->pos_relentropy(0, 0, 1, undef);
is(int(($relentA[1] * 100) + 0.5), 0,    "calculate_pos_relent() 3 seems to work (pos 2)");
is(int(($relentA[2] * 100) + 0.5), 63,   "calculate_pos_relent() 3 seems to work (pos 3)");
is(int(($relentA[3] * 100) + 0.5), 200,  "calculate_pos_relent() 3 seems to work (pos 4)");
is(int(($relentA[4] * 100) + 0.5), 128,  "calculate_pos_relent() 3 seems to work (pos 5)");
is(int(($relentA[5] * 100) + 0.5), 200,  "calculate_pos_relent() 3 seems to work (pos 6)");

my @infoctA = $msa1->pos_infocontent(0);
is(int(($infoctA[1] * 100) + 0.5), 0,   "calculate_pos_infocontent() seems to work (pos 2)");
is(int(($infoctA[2] * 100) + 0.5), 63,  "calculate_pos_infocontent() seems to work (pos 3)");
is(int(($infoctA[3] * 100) + 0.5), 200, "calculate_pos_infocontent() seems to work (pos 4)");
is(int(($infoctA[4] * 100) + 0.5), 128, "calculate_pos_infocontent() seems to work (pos 5)");
is(int(($infoctA[5] * 100) + 0.5), 200, "calculate_pos_infocontent() seems to work (pos 6)");

@relentA = $msa1->pos_relentropy(0, 0, 0, undef);
is(int(($relentA[1] * 100) + 0.5), 2,   "calculate_pos_relent() 4 seems to work (pos 2)");
is(int(($relentA[2] * 100) + 0.5), 64,  "calculate_pos_relent() 4 seems to work (pos 3)");
is(int(($relentA[3] * 100) + 0.5), 203, "calculate_pos_relent() 4 seems to work (pos 4)");
is(int(($relentA[4] * 100) + 0.5), 143, "calculate_pos_relent() 4 seems to work (pos 5)");
is(int(($relentA[5] * 100) + 0.5), 203, "calculate_pos_relent() 4 seems to work (pos 6)");

my @bgcounts_A = (30, 10, 10, 30);
@relentA = $msa1->pos_relentropy(0, 0, 0, \@bgcounts_A);
is(int(($relentA[1] * 100) + 0.5), 21,  "calculate_pos_relent() 5 seems to work (pos 2)");
is(int(($relentA[2] * 100) + 0.5), 36,  "calculate_pos_relent() 5 seems to work (pos 3)");
is(int(($relentA[3] * 100) + 0.5), 300, "calculate_pos_relent() 5 seems to work (pos 4)");
is(int(($relentA[4] * 100) + 0.5), 101, "calculate_pos_relent() 5 seems to work (pos 5)");
is(int(($relentA[5] * 100) + 0.5), 300, "calculate_pos_relent() 5 seems to work (pos 6)");

# try same thing with "" around ints
@bgcounts_A = ("30", "10", "10", "30");
@relentA = $msa1->pos_relentropy(0, 0, 0, \@bgcounts_A);
is(int(($relentA[1] * 100) + 0.5), 21,  "calculate_pos_relent() 6 seems to work (pos 2)");
is(int(($relentA[2] * 100) + 0.5), 36,  "calculate_pos_relent() 6 seems to work (pos 3)");
is(int(($relentA[3] * 100) + 0.5), 300, "calculate_pos_relent() 6 seems to work (pos 4)");
is(int(($relentA[4] * 100) + 0.5), 101, "calculate_pos_relent() 6 seems to work (pos 5)");
is(int(($relentA[5] * 100) + 0.5), 300, "calculate_pos_relent() 6 seems to work (pos 6)");

my @consA = $msa1->pos_conservation();
is(int(($consA[1] * 100) + 0.5),   0, "calculate_pos_conservation() seems to work (pos 2)");
is(int(($consA[2] * 100) + 0.5),  60, "calculate_pos_conservation() seems to work (pos 3)");
is(int(($consA[3] * 100) + 0.5), 100, "calculate_pos_conservation() seems to work (pos 4)");
is(int(($consA[4] * 100) + 0.5),  80, "calculate_pos_conservation() seems to work (pos 5)");
is(int(($consA[31] * 100) + 0.5), 40, "calculate_pos_conservation() seems to work (pos 32)");

