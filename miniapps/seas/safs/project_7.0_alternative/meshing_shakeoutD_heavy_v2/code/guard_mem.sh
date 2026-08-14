#!/bin/zsh
# guard_mem.sh -- keep the PREF finalize + ALT surgery under the 36 GB box.
#
# Two builds share the machine at the user's direction. The measured hazard: two
# ~20 GB jobs die at ~23 GB combined because SWAP runs out before RAM does, and
# macOS kills processes rather than growing the swapfile fast enough.
#
# Policy: PREF is the SACRIFICIAL job. Its LEB stage checkpoints every 10 rounds
# and restarts losslessly from the checkpoint; ALT's surgery chain must not die
# mid-stitch. So on pressure, kill the whole PREF pipeline (the finalize.sh
# wrapper FIRST -- it runs `set -uo pipefail` without -e, so killing only the
# python would let the wrapper roll on into rv repair on a stale input).
#
# Triggers (checked every 30 s):
#   * swap free < 700 MB on TWO consecutive checks   (the historical killer)
#   * combined RSS of heavy jobs > 29 GB             (headroom under 36)
LOG=${1:-/tmp/safs_guard_mem.log}
echo "[guard] start $(date)" >> $LOG
low=0
while true; do
  if ! pgrep -f "finalize.sh" >/dev/null 2>&1; then
    echo "[guard] PREF finalize gone -- exiting $(date)" >> $LOG
    exit 0
  fi
  swfree=$(sysctl -n vm.swapusage | sed 's/.*free = \([0-9.]*\)M.*/\1/')
  rss=$(ps -Ao rss=,command= | grep -E "leb_gate_close|relax_rv|rv_census|gate_census2|stitch_patch|extract_patch|refill_box|sliver_localize|tetgen" | grep -v grep | awk '{s+=$1} END {printf "%.0f", s/1048576}')
  if (( $(echo "$swfree < 700" | bc -l) )); then
    low=$((low+1))
  else
    low=0
  fi
  if [ $low -ge 2 ] || [ "${rss:-0}" -gt 29 ]; then
    echo "[guard] TRIGGER $(date): swapfree=${swfree}M rss=${rss}G -- killing PREF pipeline" >> $LOG
    pkill -f "finalize.sh" 2>/dev/null
    sleep 1
    pkill -f "shakeoutD_pref_heavy" 2>/dev/null
    echo "[guard] PREF killed; ALT continues. Restart PREF from its last checkpoint." >> $LOG
    exit 0
  fi
  sleep 30
done
