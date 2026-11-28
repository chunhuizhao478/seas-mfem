#!/bin/bash
#SBATCH -J bp5fd_e1
#SBATCH -o bp5fd_e1_%j.out
#SBATCH -e bp5fd_e1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 06:00:00
#SBATCH -A EAR24034
# BP5-FD debug plan, stage B (Chunhui Zhao, 2026-11-27: "proceed with the
# plan, test cases"): event 1 alone at p = 2 on the 1000 m mesh, one change
# per run, to 2.7e6 s (the first dynamic phase and 31 days of afterslip).
#   $1 = the run: b0 (baseline), b1 (Kelvin-Voigt off), b1p (Prakash-Clifton
#        off), b2 (cfl 0.25), b3 (absorbing outer walls), b4 (V_off 1e-4,
#        T_dwell 60 s), b5 (exact nucleation-patch edges)
# Deck config/fd_cycle/bp5_fd_p2_event1_<run>.toml, committed; a fresh run in
# $SCRATCH/scec_bp5/fd_p2_event1_<run>_<job id>.  No restart, no export.
# From miniapps/seas/jobs/fd_cycle, one job per run:
#   for r in b0 b1 b1p b2 b3 b4 b5; do sbatch vista_fd_p2_event1.sbatch $r; done
# The H-matrix is the p = 2 run's file $SCRATCH/scec_bp5/bp5_fd_p2_h1e-10.bin
# (same mesh, order and elastic operator in every run); the job refuses to
# start without it, so seven jobs do not each build it for 80 min.
# Same 4 x 144 layout as the p = 2 run (peak committed memory there
# 109.5-120.2 GiB per node).
# Expected (from the p = 2 run 1094927: 0.0036 s per dynamic step, 0.23 s per
# QD step): the dynamic phase about 25 min (b2 twice that; b4 longer, and it
# stops at T_max = 1200 s of phase time if max V stays above V_off); the QD
# part is not measured.  Cost: 4 nodes x <= 6 h = <= 8 SU per run (EAR24034).
# Binaries: vista_build.sbatch at this code (git diff <built> HEAD may touch
# only config/fd_cycle and jobs/fd_cycle).
# Submit only on Chunhui Zhao's typed go for these jobs.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NVPL_BLAS_NUM_THREADS=1 NVPL_LAPACK_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
MODULES="gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5"
VAR=${1:-}
case "$VAR" in
   b0|b1|b1p|b2|b3|b4|b5) ;;
   *) echo "ERROR: \$1 must be one of b0 b1 b1p b2 b3 b4 b5 (got '$VAR')"; exit 2 ;;
esac
DECK=config/fd_cycle/bp5_fd_p2_event1_$VAR.toml
PREFIX=e1_$VAR
HCACHE=$SCRATCH/scec_bp5/bp5_fd_p2_h1e-10.bin
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
cd "$ROOT" || exit 2

# ---- GitHub-first guard (CLAUDE.md), and the binaries match this code ----
# No fetch: the compute nodes may have no internet (vista_build.sbatch).
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
HEAD=$(git rev-parse HEAD)
BUILT=$(cat miniapps/seas/BUILD_COMMIT_VISTA 2>/dev/null)
git cat-file -e "$BUILT^{commit}" 2>/dev/null \
   || { echo "ERROR: binaries built at '$BUILT', not a commit here; run vista_build.sbatch"; exit 2; }
CODE=$(git diff --name-only "$BUILT" HEAD | grep -v -E '^miniapps/seas/(config|jobs)/fd_cycle/')
[ -z "$CODE" ] || { echo "ERROR: code changed since the build at $BUILT:"; echo "$CODE"; exit 2; }
echo "commit $HEAD, binaries built at $BUILT"

[ -f "miniapps/seas/$DECK" ] || { echo "ERROR: no deck miniapps/seas/$DECK"; exit 2; }
grep -q "^prefix       = \"$PREFIX\"$" "miniapps/seas/$DECK" \
   || { echo "ERROR: $DECK does not have prefix $PREFIX"; exit 2; }
[ -f "$HCACHE" ] || { echo "ERROR: H-matrix $HCACHE not found (the p = 2 run's file)"; exit 2; }
OUT=$SCRATCH/scec_bp5/fd_p2_event1_${VAR}_$SLURM_JOB_ID
[ ! -e "$OUT" ] || { echo "ERROR: $OUT exists"; exit 2; }
mkdir -p "$OUT" || exit 2
echo "run $VAR: deck $DECK, a fresh run in $OUT"
echo "--- $DECK against config/fd_cycle/bp5_fd_p2_1800yr.toml (comment lines left out)"
diff miniapps/seas/config/fd_cycle/bp5_fd_p2_1800yr.toml "miniapps/seas/$DECK" | grep -v -E '^[<>] *#|^[0-9,]+[acd][0-9,]+$|^---$'
echo "H-matrix: reading $HCACHE ($(stat -c %s "$HCACHE") B)"

module purge
for m in $MODULES; do module load $m || exit 2; done
cd miniapps/seas
echo "host $(hostname)  nodes $SLURM_JOB_NODELIST  ranks $SLURM_NTASKS  output $OUT"
LOG=$OUT/run_$SLURM_JOB_ID.log
srun --overlap -N "$SLURM_NNODES" -n "$SLURM_NNODES" --ntasks-per-node=1 \
   "$ROOT/miniapps/seas/jobs/fd_cycle/vista_memlog.sh" 60 > "$OUT/memlog_$SLURM_JOB_ID.txt" 2>&1 &
MEMLOG=$!
trap 'kill $MEMLOG 2>/dev/null' EXIT
s=$(date +%s)
ibrun ./seas_fd_cycle_driver \
   --deck "$DECK" --output-dir "$OUT" --qd-green-hcache "$HCACHE" --log-every 1000 > "$LOG" 2>&1
rc=$?
echo "seas_fd_cycle_driver exit $rc after $(( $(date +%s) - s )) s"
grep -E "QD stages|Green|LTS:|absorbing boundary:|benchmark output|time per dynamic step" "$LOG"
grep -A12 "dynamic phase [0-9]* summary" "$LOG"
grep -E "QD -> dynamic at|back to QD at|STOP|\[QD-GREEN\]" "$LOG"
for f in switch catalog events; do echo "--- ${PREFIX}_$f.txt"; cat "$OUT/${PREFIX}_$f.txt"; done
grep -A10 "=== fd_cycle summary" "$LOG"
echo "--- station rows"
wc -l "$OUT/${PREFIX}"_fltst_*.txt
if [ $rc -eq 0 ]; then
   grep -q "^  status: done$" "$LOG" \
      || { echo "ERROR: exit 0 but the summary does not say 'status: done'"; rc=9; }
fi
[ $rc -eq 0 ] && echo "run $VAR reached tfinal" || echo "ERROR: run $VAR ended with exit $rc"
exit $rc
