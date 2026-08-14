#!/bin/zsh
# guard_mem2.sh -- combined-memory guard while the PREF optim campaign runs
# beside the ALT surgery chain. Same policy as guard_mem.sh: PREF is the
# sacrificial job (each optim pass is ~20-40 min and restartable from the last
# accepted pass); ALT must not die mid-stitch.
#
# A rewrite of the inline v2 guard whose nested quoting broke its awk -- a guard
# with a broken trigger is worse than no guard, because it looks armed.
LOG=/tmp/safs_guard_mem.log
MARK=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_preferred/meshing_shakeoutD_heavy/build_tmp/pref_optim.running
echo "[guard3] start $(date)" >> $LOG
low=0
while [ -f $MARK ]; do
  swfree=$(sysctl -n vm.swapusage | sed 's/.*free = \([0-9.]*\)M.*/\1/')
  rss=$(ps -Ao rss=,command= | grep -E 'mmg3d|stitch_patch|extract_patch|refill_box|sliver_localize|leb_gate|relax_rv' | grep -v grep | awk '{s+=$1} END {printf "%d", s/1048576}')
  [ -z "$rss" ] && rss=0
  if (( $(echo "$swfree < 700" | bc -l) )); then low=$((low+1)); else low=0; fi
  if [ $low -ge 2 ] || [ $rss -gt 29 ]; then
    echo "[guard3] TRIGGER $(date): swapfree=${swfree}M heavyRSS=${rss}G -- killing PREF optim" >> $LOG
    pkill -f iterate_optim.sh 2>/dev/null
    sleep 1
    pkill -x mmg3d_O3 2>/dev/null
    rm -f $MARK
    echo "[guard3] PREF optim killed; ALT continues. Resume from the last accepted pass." >> $LOG
    exit 0
  fi
  sleep 30
done
echo "[guard3] campaign finished cleanly, exiting $(date)" >> $LOG
