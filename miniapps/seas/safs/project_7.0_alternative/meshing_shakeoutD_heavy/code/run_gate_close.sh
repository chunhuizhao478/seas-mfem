#!/bin/zsh
# run_gate_close.sh -- close the 1 Hz p5 gate on the rv-repaired stage-1 mesh.
#
# NOT the same job as the g=0.15 refinement that diverged. That one refined
# against the GRADED field, where the mesh already sat at median dx/h = 0.94, so
# 43.7 % of cells missed marginally and LEB's factor-of-2 bisection doubled the
# mesh to fix 5 % misses. This refines against the GATE ALONE, which the mesh
# already satisfies at 99.48 %:
#
#   safety  cells refined        added tets      final
#     1.00     327,558  0.52%        +0.3 M      63.4 M
#     1.15     707,964  1.12%        +0.7 M      63.8 M
#     1.30   1,169,083  1.85%        +1.5 M      64.5 M
#
# Under 1 M added cells, because most failures miss by only 1.08x and cost
# 1.08^3 - 1 ~ 0.26 cells each. LEB's overshoot is irrelevant at this population
# size, and LEB is local and memory-bounded where a second mmg pass would rebuild
# all 63 M tets over ~9 h to change 1 %.
#
# Why the two populations fail (measured):
#   A  82.9 %  need h >= 1400 m -- mmg treats -hmax as a bound on the METRIC, not
#              a hard cap on a tet's longest edge, so max edges run ~15 % over.
#              Median shortfall 1.08x.
#   B  17.1 %  need h <  1400 m -- the stage-1 metric floor. mmg was never asked.
#              Median shortfall 1.42x, worst 8.4x (Vs 194 m/s, needs dx <= 242 m).
#
# --gate 0.8 with --vs-source muscal: 1 Hz at p5 is (5/4)*Vs/dx >= 1, i.e.
# Vs/dx >= 0.8. Scored on NATIVE MUSCAL, never the deck cube, whose 250 m
# z-binning manufactures shallow failures (~6x more demanding).
#
# --include-parent so the whole mesh may refine and the ONLY frozen edges are the
# fault's own. NOT --allow-fault-split, which would change the DR facet count.
set -u
cd "${0:a:h}/.."
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
B=build_tmp
LOG=$B/gate_close.log
IN=${1:-$B/shakeoutD_alt_heavy_s1_rv.puml.h5}
OUT=$B/shakeoutD_alt_heavy_s1_gate.puml.h5
CVM=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/seisol_quakeworx/safs_seisol_v3_2_0_RSSRW_ALTv2_freesurf/safs_material_cvm.nc
MUSCAL=/Users/chunhuizhao/Downloads/muscal_nc/MUSCAL.nc

[ -f $IN ] || { echo "MISSING $IN -- rv repair not finished?"; exit 1; }
if pgrep -f "leb_gate_close.py|relax_rv.py" >/dev/null 2>&1 || pgrep -x mmg3d_O3 >/dev/null 2>&1; then
  echo "REFUSING: a heavy job is already running."
  ps -Ao pid,etime,rss,command | grep -E "mmg3d_O3|leb_gate_close|relax_rv" | grep -v grep \
    | awk '{printf "  pid %s  %s  %.2f GB\n", $1, $2, $3/1048576}'
  exit 1
fi

echo "[start] $(date)   in: $IN"
nohup $PY -c "
import subprocess, sys
subprocess.Popen(
    [sys.argv[1], '-u', 'code/leb_gate_close.py',
     '--mesh', '$IN', '--parent', '$IN', '--include-parent',
     '--cvm', '$CVM', '--muscal', '$MUSCAL',
     '--vs-source', 'muscal', '--gate', '0.8',
     '--drive', 'pooled', '--pool-frac', '0.25',
     '--max-rounds', '40', '--checkpoint-every', '5',
     '--stats', '$B/gate_close_stats.json', '--out', '$OUT'],
    stdout=open('$LOG', 'ab'), stderr=subprocess.STDOUT,
    start_new_session=True)
" $PY >> $LOG 2>&1 &!
sleep 4
pid=$(pgrep -f "code/leb_gate_close.py" | head -1)
[ -z "$pid" ] && { echo "[FAILED] nothing launched"; tail -15 $LOG; exit 1; }
echo "[launched] pid $pid -> $LOG"
echo "WATCH FOR DIVERGENCE: if 'gate-fail' rises for 3 consecutive rounds while"
echo "collar grows >1.05x/round, stop it -- that is the pattern the graded run showed."
