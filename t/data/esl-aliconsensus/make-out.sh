#perl ../../../scripts/esl-aliconsensus.pl --nocomment  RF00006.stk > exp.df.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --rf_no RF00006.norf.stk > exp.rfno.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --rf_ignore RF00006.stk > exp.rfignore.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --gapfract --rf_cons --rf_gapthr 0.1 RF00006.norf.stk > exp.rfcons.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --rf_mis RF00006.stk > exp.rfmis.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --rf_x RF00006.norf.stk > exp.rfx.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --skip RF00006.norf.stk > exp.skip.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --skip --skip_thr 0.001 RF00006.stk > exp.skipthr.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --skip --skip_char ! RF00006.stk > exp.skipchar.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --cons_thr1 0.75 RF00006.stk > exp.consthr1.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --cons_thr2 0.95 RF00006.stk > exp.consthr2.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --cons_no RF00006.norf.stk > exp.consno.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --cons_fract RF00006.stk > exp.consfract.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --info RF00006.stk > exp.info.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --relent RF00006.stk > exp.relent.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --gapfract RF00006.stk > exp.gapfract.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --mis RF00006.stk > exp.mis.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --weights RF00006.norf.stk > exp.weights.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --data exp.data RF00006.stk > tmp.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --info --relent --mis --gapfract --cons_fract --data exp.alldf.data RF00006.stk > exp.alldf.stk
#perl ../../../scripts/esl-aliconsensus.pl --nocomment  --info --relent --mis --gapfract --cons_fract --rf_no --data exp.allnorf.data RF00006.norf.stk > exp.allnorf.stk
perl ../../../scripts/esl-aliconsensus.pl --nocomment  --describe RF00006.norf.stk > exp.describe.out
