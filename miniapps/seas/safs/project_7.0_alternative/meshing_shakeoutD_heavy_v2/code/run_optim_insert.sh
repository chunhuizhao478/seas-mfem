#!/bin/zsh
# run_optim_insert.sh -- whole-mesh mmg -optim WITH insertion. The validated repair.
#
# Why this differs from the rejected pass: the earlier whole-mesh run used
# -noinsert and plateaued at -15 % (864,825 residual slivers) -- reshaping needs
# somewhere to PUT vertices, and -noinsert forbids exactly that. Measured on a
# 553,690-tet patch cut from the real v1 mesh (its actual slivers, its actual
# frozen fault, code/extract_patch.py):
#
#                          tets      eta<0.05          skin/fault
#   patch baseline       553,690   74,672 (13.5 %)
#   -optim, insertion    451,575   15,279 ( 3.4 %)    preserved exactly
#   second pass          399,880   14,777             plateau -- one pass is it
#
# And of the 15,279 residuals, 94.2 % touch the patch's ARTIFICIAL frozen box
# wall -- a constraint that does not exist at whole-mesh scale. Only 794 (5.2 %)
# fail free of both skin and fault. So the whole-mesh single-pass expectation is
# a >90 % sliver reduction, not 80 %.
#
# Gate interplay: the sliver band is DEEP (z median -14.6 km) where cells are
# 115-250 m against a gate demand of ~4-6 km -- collapses there cannot break the
# 1 Hz gate. Far-field collapses might; the LEB gate close (max-rounds 250) runs
# AFTER this as the finisher either way, and measured LEB leaves eta bit-identical.
#
# -nosurf + RequiredTriangles freeze every surface: fault count/area cannot change.
set -u
cd "${0:a:h}/.."
MMG=/Users/chunhuizhao/miniforge/envs/mmg/bin/mmg3d_O3
B=build_tmp
IN=$B/v1.mesh
OUT=$B/v1_optimI.mesh
LOG=$B/optimI.log
MEM=${1:-26000}

[ -f $IN ] || { echo "MISSING $IN"; exit 1; }
if pgrep -x mmg3d_O3 >/dev/null 2>&1 || pgrep -f "fill_domain.py" >/dev/null 2>&1; then
  echo "REFUSING: a heavy job is already running."; exit 1
fi

echo "[start] $(date)  mmg -optim (insertion ALLOWED) -nosurf  (-m $MEM)"
nohup /Users/chunhuizhao/miniforge/envs/pythonenv/bin/python -c "
import subprocess
subprocess.Popen(['$MMG', '-in', '$IN', '-out', '$OUT',
                  '-opnbdy', '-optim', '-nosurf',
                  '-hmin', '20', '-hmax', '200000', '-hgrad', '3',
                  '-m', '$MEM', '-v', '4'],
                 stdout=open('$LOG','ab'), stderr=subprocess.STDOUT,
                 start_new_session=True)
" >> $LOG 2>&1 &!
sleep 5
pid=$(pgrep -x mmg3d_O3 | head -1)
[ -z "$pid" ] && { echo "[FAILED] nothing launched"; tail -20 $LOG; exit 1; }
echo "[launched] pid $pid -> $LOG"
echo "ACCEPT IF: surface triangles unchanged (6,154,447 expected in = out),"
echo "           eta<0.05 falls ~10x or better, THEN re-census the 1 Hz gate."
