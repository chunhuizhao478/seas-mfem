#!/bin/bash
# R-402 regression: dead-rank detector for the tpv102 RESULT.txt post-run
# check.  Matches "effectively zero" ||Q||_inf values under the driver's
# fixed %.3e output format (M.MMMe[+-]NN).  Used by the 50-rank and 400-rank
# dev-queue sbatch scripts.
#
# Must be kept in sync with the `grep -qE` expression in
#   jobs/tpv102/tpv102_200m_p1_1.5s_50rank_dev.sbatch
#   jobs/tpv102/tpv102_200m_p1_1.5s_400rank_dev.sbatch
#
# Run manually: bash miniapps/seas/jobs/tpv102/test_result_regex.sh

set -eu
PATTERN='r[0-9]+=(0\.0{3}e\+[0-9]{2}|[0-9]\.[0-9]{3}e-[2-9][0-9])([[:space:]]|$)'

# Cases that MUST match (dead rank / effectively-zero value).
for line in \
   "[qnorm:watch] r10=0.000e+00 r11=1.000e+05 (hypo_rank=10)" \
   "[qnorm:watch] r2=2.470e-61 r11=1.000e+05 (hypo_rank=10)"  \
   "[qnorm:watch] r2=5.000e-21 r11=1.000e+05 (hypo_rank=10)"  \
   "[qnorm:watch] r2=9.821e-33 r11=1.000e+05 (hypo_rank=10)"  \
; do
   echo "$line" | grep -qE "$PATTERN" \
      || { echo "REGRESS: should match: $line"; exit 1; }
done

# Cases that MUST NOT match (alive ranks, including small-but-alive).
for line in \
   "[qnorm:watch] r10=1.000e+06 r11=5.000e+05 (hypo_rank=10)" \
   "[qnorm:watch] r10=1.234e+03 r11=5.000e+05 (hypo_rank=10)" \
   "[qnorm:watch] r10=1.234e-03 r11=5.000e+05 (hypo_rank=10)" \
   "[qnorm:watch] r10=5.000e-10 r11=5.000e+05 (hypo_rank=10)" \
   "[qnorm:watch] r10=9.999e-19 r11=5.000e+05 (hypo_rank=10)" \
; do
   if echo "$line" | grep -qE "$PATTERN"; then
      echo "REGRESS: should NOT match: $line"; exit 1
   fi
done

echo "R-402 RESULT.txt regex: OK"
