#!/bin/zsh
# Bounded watch: does PREF's mmg pass-1 analysis complete within 30 min?
LOG=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_preferred/meshing_shakeoutD_heavy/build_tmp/pref_optim_campaign.log
PID=${1:?mmg pid}
n=0
until grep -q "PHASE 1 COMPLETED" $LOG 2>/dev/null; do
  kill -0 $PID 2>/dev/null || { echo PREF-MMG-EXITED-EARLY; exit 0; }
  n=$((n+1))
  [ $n -ge 30 ] && { echo PREF-PHASE1-TIMEOUT-30MIN; exit 0; }
  sleep 60
done
echo PREF-PHASE1-COMPLETED
