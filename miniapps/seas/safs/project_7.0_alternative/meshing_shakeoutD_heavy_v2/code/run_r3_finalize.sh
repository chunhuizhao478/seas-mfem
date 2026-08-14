#!/bin/zsh
# run_r3_finalize.sh -- the unattended chain to the ALT deliverable:
#   surgery round 3 (thr 0.1, refills now POLISHED with box-local optim)
#   -> alt_finalize (PUML, gate census, LEB close to an ASSERTED zero, checks,
#      promote to results/safalt_shakeoutD_heavy_1Hz_p5_v2.puml.h5)
set -u
cd "${0:a:h}/.."
LOG=build_tmp/r3_finalize.log
echo "=== chain start $(date) ===" >> $LOG
zsh code/surgery_driver2.sh build_tmp/v2_final2.mesh 0.1 build_tmp/v2_final3.mesh >> $LOG 2>&1
if [ ! -f build_tmp/v2_final3.mesh ]; then
  echo "CHAIN-ABORT: round 3 produced no output" | tee -a $LOG
  exit 1
fi
zsh code/alt_finalize.sh build_tmp/v2_final3.mesh >> $LOG 2>&1
tail -1 build_tmp/alt_finalize.log | tee -a $LOG
