#!/bin/bash
#SBATCH --job-name=cycle2d
#SBATCH --clusters=condo
#SBATCH --partition=scec
#SBATCH --account=scec_608
#SBATCH --nodes=1 --ntasks-per-node=20 --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=24:00:00
#SBATCH --output=cycle2d_%j.out
#SBATCH --error=cycle2d_%j.err
# A 2D map-view cycle deck on Discovery (USC CARC, SCEC condo), the condo
# twin of jobs/2d_inplane/vista_2d_cycle.sbatch (Chunhui Zhao, 2026-11-05:
# "also let's try discovery, which you should have access"): any
# seas_fd_cycle_driver deck of this branch outside BP3-FD, named by
# CYCLE_DECK (required, relative to miniapps/seas).  One node, in a chain of
# 24 h jobs that all write one run folder:
#   first link (no argument):  a fresh run from t = 0 in
#                              <checkout>-runs/<prefix>_<job id>
#   later links ($1 = that folder): --restart from its last checkpoint, in
#                                   place (the files are trimmed to the
#                                   checkpoint and continued)
# The driver stops itself at the first QD step after CYCLE_WALL_S of wall
# time (default 21 h, --wall-limit): a checkpoint, then exit 8, and this job
# exits 0 so the next link (afterok) continues.  The run is finished when
# the driver exits 0: "status: done" (t_final), or "stopped in dynamic phase
# N (exit 0)" at [cycle.switch] max_dynamic_phases; this job writes DONE and
# a link whose folder has DONE exits 0 at once, so spare links cost nothing.
# Any other driver exit (3 non-finite, 4 T_max, 7 the hand-back guard) fails
# the job and the links after it never start; the log and the outputs stay
# for the diagnosis.
# Ranks: SLURM_NTASKS (20, one condo node; sbatch --ntasks-per-node=N for
# fewer).  A restart needs the rank count the run started with (the
# checkpoints and the snapshot files).
# The binaries must be built (discovery_2d_build.sbatch) from a commit whose
# code is HEAD's: git diff <built> HEAD may touch only config/fd_cycle,
# jobs/fd_cycle, jobs/2d_inplane, simple_cases/configs and
# simple_cases/scripts.
# No memlog twin here: the condo nodes enforce --mem through cgroups, where
# Vista's node-wide overcommit needed its own watcher.
# Cost: 1 node x <= 24 h per link on the SCEC condo (scec_608).
# Submit only on Chunhui Zhao's typed go for these jobs.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
MODULES="gcc/13.3.0 openmpi/5.0.5 openblas/0.3.28 hdf5/1.14.3 cmake/3.29.4 python/3.11.9"
RUN=${1:-}
DECK=${CYCLE_DECK:?set CYCLE_DECK to the deck, relative to miniapps/seas}
WALL=${CYCLE_WALL_S:-75600}                               # 21 h
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
cd "$ROOT" || exit 2

# ---- GitHub-first guard (CLAUDE.md), and the binaries match this code ----
# No fetch: run as checked out (the login-node submit fetched, as the Vista
# twin).
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
HEAD=$(git rev-parse HEAD)
BUILT=$(cat miniapps/seas/BUILD_COMMIT_DISCOVERY 2>/dev/null)
git cat-file -e "$BUILT^{commit}" 2>/dev/null \
   || { echo "ERROR: binaries built at '$BUILT', not a commit here; run discovery_2d_build.sbatch"; exit 2; }
CODE=$(git diff --name-only "$BUILT" HEAD \
       | grep -v -E '^miniapps/seas/(config/fd_cycle|jobs/fd_cycle|jobs/2d_inplane|simple_cases/configs|simple_cases/scripts)/')
[ -z "$CODE" ] || { echo "ERROR: code changed since the build at $BUILT:"; echo "$CODE"; exit 2; }
echo "commit $HEAD, binaries built at $BUILT"

cd miniapps/seas
[ -f "$DECK" ] || { echo "ERROR: no deck $DECK"; exit 2; }
git ls-files --error-unmatch "$DECK" > /dev/null 2>&1 \
   || { echo "ERROR: $DECK is not tracked (decks run as committed)"; exit 2; }
case "$(basename "$DECK")" in
   bp3_fd_*) echo "ERROR: $DECK is a BP3-FD deck: the BP3-FD jobs run those"; exit 2 ;;
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
   OUT=$ROOT-runs/${PREFIX}_$SLURM_JOB_ID
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
  echo "link $([ -z "$RUN" ] && echo first || echo "restart from QD $STEPS")"; } >> "$OUT/RUN_INFO.txt"

module purge
for m in $MODULES; do module load $m || exit 2; done
echo "host $(hostname)  nodes $SLURM_JOB_NODELIST  ranks $SLURM_NTASKS  deck $DECK  output $OUT  wall limit $WALL s"
LOG=$OUT/run_$SLURM_JOB_ID.log
s=$(date +%s)
mpirun -np "$SLURM_NTASKS" --bind-to core ./seas_fd_cycle_driver \
   --deck "$DECK" --output-dir "$OUT" "${RESTART[@]}" \
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
# The outcome from the log as well as the exit code (as the Vista twin): a
# wall-limit stop is exit 8 whatever code arrives, and exit 0 counts as the
# end only when the summary says "status: done" or "stopped in dynamic phase
# N (exit 0)".
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
