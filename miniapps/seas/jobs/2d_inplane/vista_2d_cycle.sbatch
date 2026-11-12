#!/bin/bash
#SBATCH -J cycle2d
#SBATCH -o cycle2d_%j.out
#SBATCH -e cycle2d_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 24:00:00
#SBATCH -A EAR20006
# A 2D map-view cycle deck on Vista (Chunhui Zhao, 2026-11-02: "apply to
# UCERF 2D plane strain, with 10MPa cap, and rerun the problem, ... let's
# also explore cluster to run this test because 7hr is too long"): any
# seas_fd_cycle_driver deck of this branch outside BP3-FD (the UCERF decks in
# config/fd_cycle, the simple cases in simple_cases/configs), named by
# CYCLE_DECK (required, relative to miniapps/seas).  One node, in a chain of
# 24 h jobs that all write one run folder (as vista_bp3_fd.sbatch):
#   first link (no argument):  a fresh run from t = 0 in
#                              <run root>/<prefix>_<job id>; the run root is
#                              <checkout>-runs, or CYCLE_RUN_ROOT (below)
#   later links ($1 = that folder): --restart from its last checkpoint, in
#                                   place (the files are trimmed to the
#                                   checkpoint and continued)
# The driver stops itself at the first QD step after CYCLE_WALL_S of wall
# time (default 21 h, --wall-limit): a checkpoint, then exit 8, and this job
# exits 0 so the next link (afterok) continues.  The 3 h left cover a
# dynamic phase that starts just before the limit (UCERF R70, Mac 8 ranks:
# phases of up to 3,600 s of model time at 0.011-0.0125 s per step of
# 3.47 ms, about 3.2 h each there; inferred 3-6 times faster here).
# The run is finished when the driver exits 0: "status: done" (t_final), or
# "stopped in dynamic phase N (exit 0)" at [cycle.switch]
# max_dynamic_phases; this job writes DONE and a link whose folder has DONE
# exits 0 at once, so spare links cost nothing.  Any other driver exit (3
# non-finite, 4 T_max, 7 the hand-back guard) fails the job and the links
# after it never start; the log and the outputs stay for the diagnosis.
# Ranks: SLURM_NTASKS (144 by default; sbatch -n N for fewer, e.g. the
# simple cases' 30 km fault, N_f 240, need far fewer).  A restart needs the
# rank count the run started with (the checkpoints and the snapshot files).
# The binaries must be built (vista_build.sbatch) from a commit whose code is
# HEAD's: git diff <built> HEAD may touch only config/fd_cycle, jobs/fd_cycle,
# jobs/2d_inplane, simple_cases/configs and simple_cases/scripts.
# Memory: a node refuses allocations past 118.8 GiB committed; the dense
# Green's function of the UCERF deck (N_f 3,290) took 12.4 GB on 4 ranks in
# the cloud (block 8); every minute the node's committed memory goes to
# <run>/memlog_<job id>.txt.
# Disk: the UCERF decks' volume snapshots, about 35 MB each every 10 s of a
# dynamic phase: about 6.5 GB per 1,850 s phase (the deck's note).
# CYCLE_RUN_ROOT (optional; Chunhui Zhao, 2026-11-11: "i still need volume
# output, so put results on scratch"): an absolute path under which the first
# link makes its run folder instead of <checkout>-runs, e.g.
# CYCLE_RUN_ROOT=$SCRATCH/seas-mfem-fully-dynamic-fault-opening-runs for a
# run with volume snapshots (40 events: about 190 GB).  The folder is made if
# it is missing.  Later links take the run folder as $1 as before.  TACC
# purges scratch files that were not accessed for some days (the user guide's
# purge policy): copy what must be kept out of scratch after the run.
# Cost: 1 node x <= 24 h per link on gg = <= 8 SU per link (EAR20006).
# CYCLE_DRIVER_ARGS (optional): more driver options for this link, e.g.
# --allow-other-commit for a restart from a checkpoint another commit wrote
# (only on Chunhui Zhao's approval of that change); PETSC_OPTIONS reaches
# PETSc as usual (e.g. -ts_adapt_monitor).  Both go to RUN_INFO.txt.
# Submit only on Chunhui Zhao's typed go for these jobs.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NVPL_BLAS_NUM_THREADS=1 NVPL_LAPACK_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
MODULES="gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5"
RUN=${1:-}
DECK=${CYCLE_DECK:?set CYCLE_DECK to the deck, relative to miniapps/seas}
WALL=${CYCLE_WALL_S:-75600}                               # 21 h
read -r -a EXTRA <<< "${CYCLE_DRIVER_ARGS:-}"
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
cd "$ROOT" || exit 2

# ---- GitHub-first guard (CLAUDE.md), and the binaries match this code ----
# No fetch: the compute nodes may have no internet (vista_build.sbatch);
# submit_vista_2d_cycle.sh fetches on the login node.
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
HEAD=$(git rev-parse HEAD)
BUILT=$(cat miniapps/seas/BUILD_COMMIT_VISTA 2>/dev/null)
git cat-file -e "$BUILT^{commit}" 2>/dev/null \
   || { echo "ERROR: binaries built at '$BUILT', not a commit here; run vista_build.sbatch"; exit 2; }
CODE=$(git diff --name-only "$BUILT" HEAD \
       | grep -v -E '^miniapps/seas/(config/fd_cycle|jobs/fd_cycle|jobs/2d_inplane|simple_cases/configs|simple_cases/scripts)/')
[ -z "$CODE" ] || { echo "ERROR: code changed since the build at $BUILT:"; echo "$CODE"; exit 2; }
echo "commit $HEAD, binaries built at $BUILT"

cd miniapps/seas
[ -f "$DECK" ] || { echo "ERROR: no deck $DECK"; exit 2; }
git ls-files --error-unmatch "$DECK" > /dev/null 2>&1 \
   || { echo "ERROR: $DECK is not tracked (decks run as committed)"; exit 2; }
case "$(basename "$DECK")" in
   bp3_fd_*) echo "ERROR: $DECK is a BP3-FD deck: vista_bp3_fd.sbatch runs those"; exit 2 ;;
esac
PREFIX=$(sed -n 's/^prefix *= *"\([^"]*\)"$/\1/p' "$DECK" | head -1)
[ -n "$PREFIX" ] || { echo "ERROR: no [cycle.output] prefix in $DECK"; exit 2; }
MESH=$(sed -n 's/^path *= *"\([^"]*\)"$/\1/p' "$DECK" | head -1)
[ -f "$MESH" ] || { echo "ERROR: the deck's mesh '$MESH' is missing"; exit 2; }
git ls-files --error-unmatch "$MESH" > /dev/null 2>&1 \
   || { echo "ERROR: the deck's mesh $MESH is not tracked"; exit 2; }

# ---- this link: fresh, or a restart in place ----
RESTART=()
if [ -z "$RUN" ]; then
   RUNROOT=${CYCLE_RUN_ROOT:-$ROOT-runs}
   case "$RUNROOT" in
      /*) ;;
      *) echo "ERROR: CYCLE_RUN_ROOT '$RUNROOT' is not an absolute path"; exit 2 ;;
   esac
   OUT=$RUNROOT/${PREFIX}_$SLURM_JOB_ID
   [ ! -e "$OUT" ] || { echo "ERROR: $OUT exists"; exit 2; }
   mkdir -p "$OUT" || exit 2
   echo "first link: a fresh run in $OUT"
else
   OUT=$RUN
   [ ! -e "$OUT/DONE" ] || { echo "chain finished: $OUT is done; nothing to do"; exit 0; }
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
{ echo "commit $HEAD"; echo "built $BUILT"; echo "deck $DECK"; echo "deck sha256 $(sha256sum "$DECK" | cut -d' ' -f1)"
  echo "mesh $MESH"; echo "mesh sha256 $(sha256sum "$MESH" | cut -d' ' -f1)"
  echo "job $SLURM_JOB_ID ranks $SLURM_NTASKS nodes $SLURM_JOB_NODELIST wall limit $WALL s"
  echo "link $([ -z "$RUN" ] && echo first || echo "restart from QD $STEPS")"
  echo "driver args ${CYCLE_DRIVER_ARGS:-none}  PETSC_OPTIONS ${PETSC_OPTIONS:-none}"; } >> "$OUT/RUN_INFO.txt"

module purge
for m in $MODULES; do module load $m || exit 2; done
echo "host $(hostname)  nodes $SLURM_JOB_NODELIST  ranks $SLURM_NTASKS  deck $DECK  output $OUT  wall limit $WALL s"
LOG=$OUT/run_$SLURM_JOB_ID.log
srun --overlap -N "$SLURM_NNODES" -n "$SLURM_NNODES" --ntasks-per-node=1 \
   "$ROOT/miniapps/seas/jobs/fd_cycle/vista_memlog.sh" 60 > "$OUT/memlog_$SLURM_JOB_ID.txt" 2>&1 &
MEMLOG=$!
trap 'kill $MEMLOG 2>/dev/null' EXIT
s=$(date +%s)
ibrun ./seas_fd_cycle_driver \
   --deck "$DECK" --output-dir "$OUT" "${RESTART[@]}" "${EXTRA[@]}" \
   --wall-limit "$WALL" --log-every 1000 > "$LOG" 2>&1
rc=$?
echo "seas_fd_cycle_driver exit $rc after $(( $(date +%s) - s )) s"
grep -E "restarted from|QD stages|Green|snapshots:|absorbing boundary:|dynamic stepping|sigma_n' (floor|relaxation)|wall limit .* reached" "$LOG"
grep -A12 "dynamic phase [0-9]* summary" "$LOG"
grep -E "QD -> dynamic at|back to QD at|STOP" "$LOG"
for f in switch catalog events; do echo "--- ${PREFIX}_$f.txt"; cat "$OUT/${PREFIX}_$f.txt" 2>/dev/null; done
grep -A12 "=== fd_cycle summary" "$LOG"
ls -la "$OUT"/*.vtkhdf 2>/dev/null
du -sh "$OUT"
# The outcome from the log as well as the exit code, which comes through
# ibrun and mpirun (as vista_bp3_fd.sbatch): a wall-limit stop is exit 8
# whatever code arrives, and exit 0 counts as the end only when the summary
# says "status: done" or "stopped in dynamic phase N (exit 0)".
if [ $rc -ne 8 ] && grep -q "s reached at step .*checkpoint written, stopping" "$LOG"; then
   echo "the log shows the wall-limit stop (exit code $rc): treated as exit 8"; rc=8
fi
if [ $rc -eq 0 ]; then
   grep -E "^  (final time|dynamic phases|events|status): " "$LOG"
   grep -q -E "^  status: (done|stopped in dynamic phase [0-9]+ \(exit 0\))$" "$LOG" \
      || { echo "ERROR: exit 0 but the summary says neither 'done' nor the last dynamic phase"; rc=9; }
fi
if [ $rc -eq 0 ]; then
   touch "$OUT/DONE"; echo "finished: the chain ends here"
elif [ $rc -eq 8 ]; then
   echo "stopped at the wall limit; the next link continues from its checkpoint"
   rc=0
else
   echo "ERROR: driver exit $rc; the chain stops here"
fi
exit $rc
