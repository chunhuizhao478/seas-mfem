#!/bin/zsh
# run_leb.sh -- LEB refinement of the stage-1 mesh against the g=0.15 target field.
#
# This is the stage that makes the mesh the deliverable: refining against a
# PRE-GRADED field closes the 1 Hz p5 gate AND delivers the specified g=0.15
# gradation in ONE pass, because both terms live in the field
# (h = min(h_fault(d) via PD/PH, Vs/0.8, hmax)).
#
# Measured expectations, from the PREFERRED run of the same stage:
#   * round 0 marks ~44 % of cells -- budget by MARKED cells, not mesh size
#   * residency scales with the FINAL count, ~10 GB at 70.6 M, so ~18 GB at the
#     ~130 M ALT target -- and the peak is at the END, so a comfortable reading
#     mid-run means nothing
#   * multi-hour
#
# --include-parent is what protects the fault: with it, NPAR is zeroed and the
# ONLY frozen edges are the fault's own, so the DR triangulation (and every deck
# asset keyed to its facet count) is untouched while everything else may refine.
# NOT --allow-fault-split, which would change the DR facet count.
#
# Detached via nohup + start_new_session for the same reason as the mmg pass: a
# multi-hour job must not die with whatever launched it. setsid does not exist on
# macOS -- an earlier version used it and silently launched nothing.
set -u
cd "${0:a:h}/.."
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
B=build_tmp
LOG=$B/leb.log
MESH=$B/shakeoutD_alt_heavy_s1.puml.h5
FIELD=$B/target_g015.npz
OUT=$B/shakeoutD_alt_heavy_s1_leb.puml.h5

for f in $MESH $FIELD; do
  [ -f $f ] || { echo "MISSING $f"; exit 1; }
done

if [ "${1:-}" != "--force" ]; then
  if pgrep -x mmg3d_O3 >/dev/null 2>&1 || pgrep -f "leb_gate_close.py" >/dev/null 2>&1; then
    echo "REFUSING: a heavy job is already running."
    ps -Ao pid,etime,rss,command | grep -E "mmg3d_O3|leb_gate_close" | grep -v grep \
      | awk '{printf "  pid %s  %s  %.2f GB\n", $1, $2, $3/1048576}'
    exit 1
  fi
fi

echo "[start] $(date)  swap $(sysctl -n vm.swapusage | sed 's/.*used = //;s/ .*//')"
nohup $PY -c "
import subprocess, sys
subprocess.Popen(
    [sys.argv[1], '-u', 'code/leb_gate_close.py',
     '--mesh', '$MESH', '--parent', '$MESH',
     '--target-field', '$FIELD', '--gate', '0.8', '--include-parent',
     '--drive', 'pooled', '--pool-frac', '0.25', '--max-rounds', '200',
     '--checkpoint-every', '10', '--stats', '$B/leb_stats.json',
     '--out', '$OUT'],
    stdout=open('$LOG', 'ab'), stderr=subprocess.STDOUT,
    start_new_session=True)
" $PY >> $LOG 2>&1 &!
sleep 4
pid=$(pgrep -f "code/leb_gate_close.py" | head -1)
if [ -z "$pid" ]; then
  echo "[FAILED] nothing launched -- see $LOG"; tail -15 $LOG; exit 1
fi
echo "[launched] pid $pid -> $LOG"
