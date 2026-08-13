#!/bin/zsh
# run_qfill.sh -- QUALITY base fill: the one flag that fixes the sliver problem.
#
# The v1 mesh's base was `tetgen -pY`: a constrained Delaunay with NO interior
# Steiner points. It CONNECTS the PLC, it does not SHAPE it, and mmg could never
# recover from that input. Measured on two independent boxes of real ALT fault
# geometry (code/box_switch_sweep.py):
#
#   deep junction box       tets    eta<0.05   median eta   fault lost
#     pY                 270,352   160,909 (59.5 %)  0.0289          0
#     pq1.4/10Y          794,242     6,986 ( 0.88 %) 0.8128          0
#   SE tip box (where recoversubfaces failed in v1)
#     pq1.4/10Y          460,526     9,230 ( 2.0  %) 0.8087          0
#
# Sliver rate falls 67x, median eta rises 28x, at ~3x the tets, and ZERO fault
# facets are lost. The fault was never the obstacle: the 87.2 M v1 mesh has ZERO
# tets with 2+ fault faces, its worst 1000 tets have NO fault face, and the exact
# eta ceiling over all 2,564,480 fault triangles is 0.106 -- the target is
# geometrically attainable with margin.
#
# MANDATORY FLAGS, each for a recorded reason:
#   --minratio 1.4 --mindihedral 10   the whole point: allows interior Steiner points
#   --hmax 0                          drops `-a`. The one recorded tetgen death was
#                                     the `-a` variant: 134,481,852 tets queued, died
#                                     in locate_point_walk. `-q` alone is bounded in
#                                     shape only, measured ~3x growth.
#   --shuffle 1                       NOT cosmetic. In input order this PLC dies in
#                                     recoversubfaces. Seed 1 is the recorded best
#                                     (279 Steiner, 6 facets lost); seed 2 gave 8.
#   (never drop -Y)                   recorded: without it, 215,780 fault facets are
#                                     silently re-diagonalised (8.4 % of the fault).
#
# Budget: 608 B/tet measured x ~35 M expected = ~21 GB; 12-30 h.
set -u
cd "${0:a:h}/.."
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
LOG=build_tmp/qfill.log

for f in build_tmp/hull.npz build_tmp/fault_surface.npz; do
  [ -f $f ] || { echo "MISSING $f"; exit 1; }
done
if pgrep -x mmg3d_O3 >/dev/null 2>&1 || pgrep -f "fill_domain.py|leb_gate_close.py|relax_rv.py" >/dev/null 2>&1; then
  echo "REFUSING: a heavy job is already running."
  ps -Ao pid,etime,rss,command | grep -E "mmg3d_O3|fill_domain|leb_gate_close|relax_rv" | grep -v grep \
    | awk '{printf "  pid %s  %s  %.2f GB\n", $1, $2, $3/1048576}'
  exit 1
fi

echo "[start] $(date)  quality fill -pq1.4/10Y (shuffle 1)"
nohup $PY -c "
import subprocess, sys
subprocess.Popen(
    [sys.argv[1], '-u', 'code/fill_domain.py', '--out', 'build_tmp/fill.npz',
     '--minratio', '1.4', '--mindihedral', '10', '--hmax', '0', '--shuffle', '1'],
    stdout=open('$LOG', 'ab'), stderr=subprocess.STDOUT, start_new_session=True)
" $PY >> $LOG 2>&1 &!
sleep 5
pid=$(pgrep -f "code/fill_domain.py" | head -1)
[ -z "$pid" ] && { echo "[FAILED] nothing launched"; tail -20 $LOG; exit 1; }
echo "[launched] pid $pid -> $LOG"
echo
echo "GO/NO-GO on completion:"
echo "  volume must be 12,480,379 km3; PLC vertices preserved < 1e-6 m;"
echo "  fault facets lost <= 6 (v1 lost 6 here); median eta > 0.75 (v1 base was 0.0054)."
echo "  If tetgen exceeds ~30 GB, kill it -- the box says 21 GB, and swap is the"
echo "  binding constraint on this machine, not RAM."
