#!/bin/bash
#SBATCH -J fdc_p2_1800
#SBATCH -o fdc_p2_1800_%j.out
#SBATCH -e fdc_p2_1800_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 48:00:00
#SBATCH -A EAR20006
# Fully dynamic BP5 cycle on Vista: SCEC BP5-FD to 1800 yr at p = 2 on the
# 1000 m mesh with the QD and dynamic speed-ups and every SCEC output
# (config/fd_cycle/bp5_fd_p2_1800yr.toml), in a chain of 48 h jobs that all
# write one run folder (as discovery_1500yr.sbatch):
#   first link (no argument):  a fresh run from t = 0 in
#                              <checkout>-runs/p2_1800_<job id>
#   later links ($1 = that folder): --restart from its last checkpoint, in
#                                   place (the files are trimmed to the
#                                   checkpoint and continued)
# The driver stops itself at the first QD step after FDC_WALL_S of wall time
# (45 h, --wall-limit): a checkpoint, then exit 8, and this job exits 0 so
# the next link (afterok) continues.  The 3 h left cover a dynamic phase
# that starts just before the limit (the p = 2 event 1 took 5.6 h on 80
# Discovery ranks without LTS, so well under 3 h on 576 ranks with it,
# inferred).  Exit 0 of the driver (t_final reached) writes DONE and the
# SCEC upload files (cycle/scec_bp5_export.py into <run>/scec); a link
# whose folder has DONE exits 0 at once, so spare links cost nothing.  Any
# other driver exit fails the job and the links after it never start.
# The H-matrix of the Green's function: FDC_HCACHE (default
# $SCRATCH/seas-qd-green/bp5_fd_p2_h1e-10.bin).  The first run without it
# builds and writes it (about 80 min, from 1086777); every later link reads
# it.  vista_p2_check.sbatch builds it, so the chain should start after it.
# Chain, from miniapps/seas/jobs/fd_cycle:
#   R=$(cd ../../../.. && pwd)-runs
#   j1=$(sbatch --parsable vista_p2_1800yr.sbatch)
#   j2=$(sbatch --parsable --dependency=afterok:$j1 vista_p2_1800yr.sbatch $R/p2_1800_$j1)
#   j3=$(sbatch --parsable --dependency=afterok:$j2 vista_p2_1800yr.sbatch $R/p2_1800_$j1)
# Expected: QD about 24-34 h for 1800 yr (the QD thread's estimate for
# 1500 yr at 0.222 s per step, scaled), plus about seven dynamic phases;
# two links and a spare (inferred).
# The binaries must be built from a commit whose code is HEAD's: git diff
# <built> HEAD may touch only config/fd_cycle and jobs/fd_cycle.
# FDC_DECK, FDC_PREFIX, FDC_WALL_S and FDC_HCACHE override the deck, the
# file prefix, the wall limit and the H-matrix file (for tests only).
# Memory: Vista compute nodes refuse allocations past 118.8 GiB committed
# per node (vm.overcommit_memory = 2); 1086777 ran the p = 2 QD regime
# with this G on these 4 nodes.
# Cost: 4 nodes x <= 48 h per link on gg = <= 64 SU per link (EAR20006).
# Submit only on Chunhui Zhao's typed go for these jobs.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NVPL_BLAS_NUM_THREADS=1 NVPL_LAPACK_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
MODULES="gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5"
RUN=${1:-}
DECK=${FDC_DECK:-config/fd_cycle/bp5_fd_p2_1800yr.toml}
PREFIX=${FDC_PREFIX:-bp5_fd_p2_1800yr}
WALL=${FDC_WALL_S:-162000}                                # 45 h
HCACHE=${FDC_HCACHE:-$SCRATCH/seas-qd-green/bp5_fd_p2_h1e-10.bin}
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

# ---- this link: fresh, or a restart in place ----
RESTART=()
if [ -z "$RUN" ]; then
   OUT=$ROOT-runs/p2_1800_$SLURM_JOB_ID
   [ ! -e "$OUT" ] || { echo "ERROR: $OUT exists"; exit 2; }
   mkdir -p "$OUT" || exit 2
   echo "first link: a fresh run in $OUT"
else
   OUT=$RUN
   [ ! -e "$OUT/DONE" ] || { echo "chain finished: $OUT reached tfinal; nothing to do"; exit 0; }
   [ -f "$OUT/${PREFIX}_checkpoint_r0.txt" ] \
      || { echo "ERROR: no checkpoint $OUT/${PREFIX}_checkpoint_r0.txt"; exit 2; }
   NCK=$(ls "$OUT/${PREFIX}"_checkpoint_r*.txt | wc -l)
   [ "$NCK" = "$SLURM_NTASKS" ] \
      || { echo "ERROR: $NCK rank checkpoints in $OUT for $SLURM_NTASKS ranks"; exit 2; }
   STEPS=$(grep -m1 -h '^step ' "$OUT/${PREFIX}"_checkpoint_r*.txt | sort -u)
   [ "$(echo "$STEPS" | wc -l)" = 1 ] \
      || { echo "ERROR: the rank checkpoints in $OUT disagree:"; echo "$STEPS"; exit 2; }
   echo "continuing $OUT from its checkpoint at QD $STEPS"
   RESTART=(--restart "$OUT/$PREFIX")
fi
mkdir -p "$(dirname "$HCACHE")" || exit 2
[ -f "$HCACHE" ] && echo "H-matrix: reading $HCACHE ($(stat -c %s "$HCACHE") B)" \
   || echo "H-matrix: $HCACHE not found; this run builds and writes it"

module purge
for m in $MODULES; do module load $m || exit 2; done
cd miniapps/seas
echo "host $(hostname)  nodes $SLURM_JOB_NODELIST  ranks $SLURM_NTASKS  output $OUT  wall limit $WALL s"
LOG=$OUT/run_$SLURM_JOB_ID.log
s=$(date +%s)
mpirun -np "$SLURM_NTASKS" --map-by ppr:144:node --bind-to core ./seas_fd_cycle_driver \
   --deck "$DECK" --output-dir "$OUT" "${RESTART[@]}" --qd-green-hcache "$HCACHE" \
   --wall-limit "$WALL" --log-every 1000 > "$LOG" 2>&1
rc=$?
echo "seas_fd_cycle_driver exit $rc after $(( $(date +%s) - s )) s"
grep -E "restarted from|QD stages|Green|LTS:|benchmark output|snapshots:|wall limit .* reached" "$LOG"
grep -A12 "dynamic phase [0-9]* summary" "$LOG"
grep -E "QD -> dynamic at|back to QD at|STOP|\[QD-GREEN\]" "$LOG"
for f in switch catalog events; do echo "--- ${PREFIX}_$f.txt"; cat "$OUT/${PREFIX}_$f.txt"; done
grep -A10 "=== fd_cycle summary" "$LOG"
ls -la "$OUT"/*.vtkhdf 2>/dev/null
if [ $rc -eq 0 ]; then
   touch "$OUT/DONE"; echo "reached tfinal: the chain is finished"
   python3 cycle/scec_bp5_export.py "$OUT" "$PREFIX" "$OUT/scec" --modeler "Chunhui Zhao" \
      --element-size "1000 m" --order 2 --version "$HEAD" \
      || echo "WARNING: cycle/scec_bp5_export.py failed; run it again by hand"
elif [ $rc -eq 8 ]; then
   echo "stopped at the wall limit; the next link continues from its checkpoint"
   rc=0
else
   echo "ERROR: driver exit $rc; the chain stops here"
fi
exit $rc
