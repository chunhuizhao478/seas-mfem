#!/bin/bash
#SBATCH --job-name=fdc_1500
#SBATCH --clusters=condo
#SBATCH --partition=scec
#SBATCH --account=scec_608
#SBATCH --nodes=2 --ntasks-per-node=20 --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=48:00:00
#SBATCH --output=fdc_1500_%j.out
#SBATCH --error=fdc_1500_%j.err
# Fully dynamic BP5 cycle on Discovery: SCEC BP5-FD to 1500 yr with fault
# and free-surface snapshots (config/fd_cycle/bp5_fd_1500yr.toml), p = 1 on
# the 1000 m mesh, in a chain of 48 h jobs that all write one run folder.
#   first link (no argument):  a fresh run from t = 0 in runs/fd1500_<job id>
#   later links ($1 = that folder): --restart from its last checkpoint, in
#                                   place (the files are trimmed to the
#                                   checkpoint and continued)
# The driver stops itself at the first QD step after FDC_WALL_S of wall time
# (45 h, --wall-limit): a checkpoint, then exit 8, and this job exits 0 so
# the next link (afterok) continues.  The 3 h left cover a dynamic phase
# (about 1.5 h at p = 1) that starts just before the limit.  Exit 0 of the
# driver (t_final reached) writes DONE, and a link whose folder has DONE
# exits 0 at once, so spare links cost nothing.  Any other driver exit (a
# NaN, the guard, a phase past T_max) fails the job and the links after it
# never start.  If Slurm kills a job instead (TIMEOUT), the snapshot files
# may be left unreadable; the station files are text and are trimmed on the
# next restart.
# Chain, from miniapps/seas/jobs/fd_cycle (condo prints "<id>;condo"):
#   R=$(cd ../../../../.. && pwd)/runs
#   j1=$(sbatch --parsable discovery_1500yr.sbatch | cut -d';' -f1)
#   j2=$(sbatch --parsable --dependency=afterok:$j1 discovery_1500yr.sbatch $R/fd1500_$j1 | cut -d';' -f1)
#   j3=$(sbatch --parsable --dependency=afterok:$j2 discovery_1500yr.sbatch $R/fd1500_$j1 | cut -d';' -f1)
#   j4=$(sbatch --parsable --dependency=afterok:$j3 discovery_1500yr.sbatch $R/fd1500_$j1 | cut -d';' -f1)
# Expected: about 23 h per 275 yr cycle on 40 ranks (P4), so 1500 yr in
# about 127 h, three links and a spare (inferred).
# The binaries must be built from a commit whose code is HEAD's: git diff
# <built> HEAD may touch only config/fd_cycle and jobs/fd_cycle.
# FDC_DECK, FDC_PREFIX and FDC_WALL_S override the deck, the file prefix and
# the wall limit (for tests only).
# Cost: 2 nodes x <= 48 h per link on the SCEC condo (scec_608).
# Submit only on Chunhui Zhao's typed go for these jobs.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
MODULES="gcc/13.3.0 openmpi/5.0.5 openblas/0.3.28 hdf5/1.14.3 cmake/3.29.4 python/3.11.9"
RUN=${1:-}
DECK=${FDC_DECK:-config/fd_cycle/bp5_fd_1500yr.toml}
PREFIX=${FDC_PREFIX:-bp5_fd_1500yr}
WALL=${FDC_WALL_S:-162000}                                # 45 h
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
cd "$ROOT" || exit 2

# ---- GitHub-first guard (CLAUDE.md), and the binaries match this code ----
git fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
HEAD=$(git rev-parse HEAD)
BUILT=$(cat miniapps/seas/BUILD_COMMIT_DISCOVERY 2>/dev/null)
git cat-file -e "$BUILT^{commit}" 2>/dev/null \
   || { echo "ERROR: binaries built at '$BUILT', not a commit here; run discovery_build.sbatch"; exit 2; }
CODE=$(git diff --name-only "$BUILT" HEAD | grep -v -E '^miniapps/seas/(config|jobs)/fd_cycle/')
[ -z "$CODE" ] || { echo "ERROR: code changed since the build at $BUILT:"; echo "$CODE"; exit 2; }
echo "commit $HEAD, binaries built at $BUILT"

# ---- this link: fresh, or a restart in place ----
RESTART=()
if [ -z "$RUN" ]; then
   OUT=$(dirname "$ROOT")/runs/fd1500_$SLURM_JOB_ID
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

module purge
for m in $MODULES; do module load $m || exit 2; done
cd miniapps/seas
echo "host $(hostname)  nodes $SLURM_JOB_NODELIST  ranks $SLURM_NTASKS  output $OUT  wall limit $WALL s"
LOG=$OUT/run_$SLURM_JOB_ID.log
s=$(date +%s)
mpirun -np "$SLURM_NTASKS" --map-by ppr:20:node --bind-to core ./seas_fd_cycle_driver \
   --deck "$DECK" --output-dir "$OUT" "${RESTART[@]}" --wall-limit "$WALL" \
   --log-every 1000 > "$LOG" 2>&1
rc=$?
echo "seas_fd_cycle_driver exit $rc after $(( $(date +%s) - s )) s"
grep -E "restarted from|snapshots:|wall limit .* reached" "$LOG"
grep -A12 "dynamic phase [0-9]* summary" "$LOG"
grep -E "QD -> dynamic at|back to QD at|STOP" "$LOG"
for f in switch catalog; do echo "--- ${PREFIX}_$f.txt"; cat "$OUT/${PREFIX}_$f.txt"; done
grep -A10 "=== fd_cycle summary" "$LOG"
ls -la "$OUT"/*.vtkhdf 2>/dev/null
for f in "$OUT"/*.vtkhdf; do [ -e "$f" ] && h5ls "$f/VTKHDF/Steps/Values"; done   # snapshot counts
if [ $rc -eq 0 ]; then
   touch "$OUT/DONE"; echo "reached tfinal: the chain is finished"
elif [ $rc -eq 8 ]; then
   echo "stopped at the wall limit; the next link continues from its checkpoint"
   rc=0
else
   echo "ERROR: driver exit $rc; the chain stops here"
fi
exit $rc
