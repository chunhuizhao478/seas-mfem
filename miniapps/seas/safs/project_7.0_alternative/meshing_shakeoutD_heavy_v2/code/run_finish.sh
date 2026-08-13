#!/bin/zsh
# run_finish.sh -- rv repair then gate closing, chained, detached, unattended.
#
# Chained rather than run separately because the two stages have no decision
# between them and a round-trip would idle the box for the length of a
# notification cycle.
#
# ITERS is capped deliberately. relax_rv defaults to 120, but the returns decay
# hard and it writes ONLY at the end -- no checkpoint, no convergence break -- so
# a long run is all-or-nothing. Measured on this mesh:
#
#   initial   rv>2 559,382   rv>3 261,973   med 1.4209
#   iter  0   rv>2 434,611   rv>3 230,165   med 1.4064   (139 s)
#   iter 20   rv>2 388,614   rv>3 219,728   med 1.3963   (3831 s, 185 s/iter)
#
# Iteration 0 alone delivered 22 % of the reduction; the next 20 added 10.6 %.
# That is the structural limit, not a tuning failure: every sliver in this mesh
# is fault-adjacent (97.4 % within 250 m, zero beyond 2 km), so three of each
# tet's four vertices are frozen and only the apex can move.
#
# The rv>3 tail is what actually matters -- it is the empirical tensile-flip risk
# cliff (48.8x above it) -- and relax_rv already refuses any iteration that grows
# it, so the tail is protected regardless of iteration count.
set -u
cd "${0:a:h}/.."
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
B=build_tmp
ITERS=${1:-20}
LOG=$B/finish.log
MESH=$B/shakeoutD_alt_heavy_s1.puml.h5
RV=$B/shakeoutD_alt_heavy_s1_rv.puml.h5
GATE=$B/shakeoutD_alt_heavy_s1_gate.puml.h5
CVM=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/seisol_quakeworx/safs_seisol_v3_2_0_RSSRW_ALTv2_freesurf/safs_material_cvm.nc
MUSCAL=/Users/chunhuizhao/Downloads/muscal_nc/MUSCAL.nc

[ -f $MESH ] || { echo "MISSING $MESH"; exit 1; }
if pgrep -f "leb_gate_close.py|relax_rv.py" >/dev/null 2>&1 || pgrep -x mmg3d_O3 >/dev/null 2>&1; then
  echo "REFUSING: a heavy job is already running."; exit 1
fi

cat > $B/_finish_inner.sh <<INNER
#!/bin/zsh
set -uo pipefail          # NOT -e: a stage exiting non-zero on a fault-count
                          # guard must not kill the chain before it can be judged
echo "=== 1/2  rv repair (--iters $ITERS, per-tet guards) \$(date) ==="
$PY -u code/relax_rv.py --mesh $MESH --out $RV \\
    --cvm $CVM --muscal $MUSCAL --gate 0.8 --qual-tol 0.05 --iters $ITERS
echo "rv exit=\$?"
if [ ! -f $RV ]; then echo "ABORT: rv produced no output"; exit 1; fi
echo
echo "=== 2/2  gate close (1 Hz p5, native MUSCAL) \$(date) ==="
$PY -u code/leb_gate_close.py --mesh $RV --parent $RV --include-parent \\
    --cvm $CVM --muscal $MUSCAL --vs-source muscal --gate 0.8 \\
    --drive pooled --pool-frac 0.25 --max-rounds 40 --checkpoint-every 5 \\
    --stats $B/gate_close_stats.json --out $GATE
echo "gate exit=\$?"
echo "=== finished \$(date) ==="
INNER
chmod +x $B/_finish_inner.sh

echo "[start] $(date)  rv --iters $ITERS, then gate close"
nohup $PY -c "
import subprocess, sys
subprocess.Popen(['/bin/zsh', '$B/_finish_inner.sh'],
                 stdout=open('$LOG','ab'), stderr=subprocess.STDOUT,
                 start_new_session=True)
" >> $LOG 2>&1 &!
sleep 5
pid=$(pgrep -f "code/relax_rv.py" | head -1)
[ -z "$pid" ] && { echo "[FAILED] nothing launched"; tail -20 $LOG; exit 1; }
echo "[launched] rv pid $pid -> $LOG"
