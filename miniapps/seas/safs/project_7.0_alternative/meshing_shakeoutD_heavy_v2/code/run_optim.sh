#!/bin/zsh
# run_optim.sh -- mmg OPTIMISATION-ONLY pass on the v1 mesh. Repair, not rebuild.
#
# Flags mirror mmg_cli_cleanup.py, which is the repo's canonical de-sliver driver
# and whose documented flow is exactly `mmg3d_O3 -opnbdy -optim -noinsert -nosurf`.
# That MODE has never been run in this repo at any scale -- the recorded "mmg
# passes don't help" measurement used the SIZEMAP route on a 4 km box, which is a
# different mode and a reproducer too small to stress the real problem.
#
# Why it should work here: the slivers are INTERIOR, not fault-locked. Measured on
# all 87,240,911 tets -- ZERO tets with 2+ fault faces; the worst 1000 (eta
# 0.000000 .. 0.000698) have NO fault face; fault-face tets are 5.88 % of the mesh
# but only 4.51 % of the sub-0.05 population. And the exact eta ceiling over all
# 2,564,480 fault triangles is 0.106, so eta_min ~0.076 is attainable.
#
#   -optim     improve shape using the mesh's own sizes; no user metric
#   -noinsert  mmg may swap / collapse / move but NOT add points -- so it cannot
#              run away in cell count the way the gate LEB did
#   -nosurf    do not touch any surface: the flat top stays exactly flat and the
#              fault cannot move. Belt and braces with RequiredTriangles.
#   -hmin 50   collapse floor. v1's min edge is 0.83 m; the gate needs ~194 m even
#              in the slowest MUSCAL material, so clearing sub-50 m interior edges
#              cannot break the gate, and those edges ARE the sliver signature.
#
# RISKS, both watched afterwards rather than assumed away:
#   1. -optim can COLLAPSE, so the 1 Hz gate (v1: 35,353 below, 0.041 %) may
#      regress. Re-score on native MUSCAL after. If it regresses, leb_gate_close
#      re-runs in ~40 min.
#   2. NO memory or runtime figure exists anywhere in this repo for mmg with a
#      50 M+ tet input -- every measured run grew a ~9.6 M input. -m caps it, and
#      the watcher kills it rather than let swap take the machine down.
set -u
cd "${0:a:h}/.."
MMG=/Users/chunhuizhao/miniforge/envs/mmg/bin/mmg3d_O3
B=build_tmp
IN=$B/v1.mesh
OUT=$B/v1_optim.mesh
LOG=$B/optim.log
MEM=${1:-26000}

[ -f $IN ] || { echo "MISSING $IN -- conversion not finished?"; exit 1; }
if pgrep -x mmg3d_O3 >/dev/null 2>&1; then echo "REFUSING: an mmg3d is already running."; exit 1; fi

echo "[start] $(date)  mmg -optim -noinsert -nosurf  (-m $MEM)"
nohup /Users/chunhuizhao/miniforge/envs/pythonenv/bin/python -c "
import subprocess, sys
subprocess.Popen(['$MMG', '-in', '$IN', '-out', '$OUT',
                  '-opnbdy', '-optim', '-noinsert', '-nosurf',
                  '-hmin', '50', '-hmax', '200000', '-hgrad', '3',
                  '-m', '$MEM', '-v', '5'],
                 stdout=open('$LOG','ab'), stderr=subprocess.STDOUT,
                 start_new_session=True)
" >> $LOG 2>&1 &!
sleep 5
pid=$(pgrep -x mmg3d_O3 | head -1)
[ -z "$pid" ] && { echo "[FAILED] nothing launched"; tail -20 $LOG; exit 1; }
echo "[launched] pid $pid -> $LOG"
echo "ACCEPT IF: eta_min >= 0.05 (stretch 0.076), 0 cells below 0.05,"
echo "           fault triangle count unchanged, gate not regressed."
