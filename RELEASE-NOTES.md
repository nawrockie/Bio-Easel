# Bio-Easel 0.x release notes 

### Bio-Easel 0.18 release (September 2026): Minor update
  * Adds addGS and a family of getGS_*/hasGS_* subroutines to MSA.pm
    for reading and testing for #=GS per-sequence annotation, as
    analogs of the existing #=GR per-residue annotation subroutines.

---

### Bio-Easel 0.17 release (September 2025): Minor update
  * Adds esl-aliconsensus.pl script for adding GC annotation
    that summarizes the per-column conservation to a stockholm
    alignment file.
  * Adds esl-alitransfer.pl script for transferring GC and GF
    annotation from one stockholm alignment file to another.
  * In esl-compare2rf.pl, adds --seqrf <s> option to specify
    reference sequence is <s>

---

### Bio-Easel 0.16 release (December 2022): Minor update
  * Adds esl-alicapitalize.pl script for enforcing Infernal/HMMER
    conventions related to gap/nongap RF columns on a Stockholm
    alignment file. 

---

### Bio-Easel 0.15 release (June 2021): Minor update
  * Adds support for numbering RF columns MSA.pm:capitalize_based_on_rf().
    Implemented specifically for vadr 1.2.2.

---

For more information, see the [git log for the develop
branch](https://github.com/nawrockie/Bio-Easel/commits/develop).

