#!/bin/bash
#SBATCH -J bp3fd
#SBATCH -o bp3fd_%j.out
#SBATCH -e bp3fd_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 24:00:00
#SBATCH -A EAR20006
# SCEC BP3-FD on Vista (Chunhui Zhao, 2026-10-27: "cancel this run, add
# snapshots, submit it on cluster"): config/fd_cycle/bp3_fd_d60_thrust.toml,
# the 25 m mesh at p = 1 to 1500 yr, the QD and dynamic regimes switched on
# max V, the dynamic regime's bottom and sides absorbing, fault and surface
# snapshots.  144 ranks on one gg node, in a chain of 24 h jobs that all
# write one run folder (as jobs/fd_cycle/vista_p2_1800yr.sbatch):
#   first link (no argument):  a fresh run from t = 0 in
#                              <checkout>-runs/bp3fd_<deck>_<job id>
#   later links ($1 = that folder): --restart from its last checkpoint, in
#                                   place (the files are trimmed to the
#                                   checkpoint and continued)
# The driver stops itself at the first QD step after BP3FD_WALL_S of wall
# time (22 h, --wall-limit): a checkpoint, then exit 8, and this job exits 0
# so the next link (afterok) continues.  The 2 h left cover a dynamic phase
# that starts just before the limit: on the Mac (4 ranks) event 1 ran 0.014 s
# per dynamic step at dt 1.757e-4 s, and its phase lasts about 110 s (626,000
# steps, 2.4 h there); 144 ranks should be 15 to 30 times faster (inferred,
# not measured).  The Mac reached event 1 (176.6 yr) after 14 min of QD.
# Exit 0 of the driver ("status: done", t_final reached) writes DONE and the
# SCEC upload files (cycle/scec_bp3_export.py into <run>/scec); a link whose
# folder has DONE exits 0 at once, so spare links cost nothing.  Any other
# driver exit fails the job and the links after it never start.
# Snapshots (scaled from the 1 km files, inferred): about 0.26 MB per fault
# snapshot (3,200 points) and 0.02 MB per surface one (380 points), about
# 1,100 per dynamic phase of 110 s at dt_dyn 0.1 s: about 0.3 GB per event,
# about 5 GB for 1500 yr.
# The binaries must be built (jobs/2d_inplane/vista_build.sbatch) from a
# commit whose code is HEAD's: git diff <built> HEAD may touch only
# config/fd_cycle, jobs/fd_cycle and jobs/2d_inplane.
# BP3FD_DECK picks another BP3-FD deck (its [cycle.output] prefix names the
# files); BP3FD_WALL_S the wall limit (tests only).  A restart needs the rank
# count the run started with (the checkpoints and the snapshot files).
# Memory: the Mac run peaked at 7.7 GB in all (MUMPS and the Green's
# function at set-up); a node refuses allocations past 118.8 GiB committed.
# Every minute the node's committed memory goes to <run>/memlog_<job id>.txt.
# Cost: 1 node x <= 24 h per link on gg = <= 8 SU per link (EAR20006).
# Submit only on Chunhui Zhao's typed go for these jobs.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NVPL_BLAS_NUM_THREADS=1 NVPL_LAPACK_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
MODULES="gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5"
RUN=${1:-}
DECK=${BP3FD_DECK:-config/fd_cycle/bp3_fd_d60_thrust.toml}
WALL=${BP3FD_WALL_S:-79200}                               # 22 h
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
CODE=$(git diff --name-only "$BUILT" HEAD | grep -v -E '^miniapps/seas/(config/fd_cycle|jobs/fd_cycle|jobs/2d_inplane)/')
[ -z "$CODE" ] || { echo "ERROR: code changed since the build at $BUILT:"; echo "$CODE"; exit 2; }
echo "commit $HEAD, binaries built at $BUILT"

cd miniapps/seas
[ -f "$DECK" ] || { echo "ERROR: no deck $DECK"; exit 2; }
PREFIX=$(sed -n 's/^prefix *= *"\([^"]*\)"$/\1/p' "$DECK" | head -1)
[ -n "$PREFIX" ] || { echo "ERROR: no [cycle.output] prefix in $DECK"; exit 2; }
NAME=${PREFIX#bp3_fd_}
case "$NAME" in
   d30_thrust) SCEN="dip 30, thrust" ;;  d30_normal) SCEN="dip 30, normal" ;;
   d60_thrust) SCEN="dip 60, thrust" ;;  d60_normal) SCEN="dip 60, normal" ;;
   d90) SCEN="dip 90" ;;
   *) echo "ERROR: $DECK is not one of the five BP3-FD decks (prefix $PREFIX)"; exit 2 ;;
esac

# ---- this link: fresh, or a restart in place ----
RESTART=()
if [ -z "$RUN" ]; then
   OUT=$ROOT-runs/bp3fd_${NAME}_$SLURM_JOB_ID
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
echo "host $(hostname)  nodes $SLURM_JOB_NODELIST  ranks $SLURM_NTASKS  deck $DECK  output $OUT  wall limit $WALL s"
LOG=$OUT/run_$SLURM_JOB_ID.log
srun --overlap -N "$SLURM_NNODES" -n "$SLURM_NNODES" --ntasks-per-node=1 \
   "$ROOT/miniapps/seas/jobs/fd_cycle/vista_memlog.sh" 60 > "$OUT/memlog_$SLURM_JOB_ID.txt" 2>&1 &
MEMLOG=$!
trap 'kill $MEMLOG 2>/dev/null' EXIT
s=$(date +%s)
ibrun ./seas_fd_cycle_driver \
   --deck "$DECK" --output-dir "$OUT" "${RESTART[@]}" \
   --wall-limit "$WALL" --log-every 1000 > "$LOG" 2>&1
rc=$?
echo "seas_fd_cycle_driver exit $rc after $(( $(date +%s) - s )) s"
grep -E "restarted from|QD stages|Green|benchmark output|snapshots:|absorbing boundary:|dynamic stepping|wall limit .* reached" "$LOG"
grep -A12 "dynamic phase [0-9]* summary" "$LOG"
grep -E "QD -> dynamic at|back to QD at|STOP" "$LOG"
for f in switch catalog events; do echo "--- ${PREFIX}_$f.txt"; cat "$OUT/${PREFIX}_$f.txt"; done
grep -A10 "=== fd_cycle summary" "$LOG"
ls -la "$OUT"/*.vtkhdf 2>/dev/null
# The outcome from the log as well as the exit code, which comes through
# ibrun and mpirun (as vista_p2_1800yr.sbatch): a wall-limit stop is exit 8
# whatever code arrives, and exit 0 counts as the end only when the driver's
# summary says "status: done".
if [ $rc -ne 8 ] && grep -q "s reached at step .*checkpoint written, stopping" "$LOG"; then
   echo "the log shows the wall-limit stop (exit code $rc): treated as exit 8"; rc=8
fi
if [ $rc -eq 0 ]; then
   grep -E "^  (t_final|final time|status): " "$LOG"
   grep -q "^  status: done$" "$LOG" \
      || { echo "ERROR: exit 0 but the summary does not say 'status: done'"; rc=9; }
fi
if [ $rc -eq 0 ]; then
   touch "$OUT/DONE"; echo "reached tfinal: the chain is finished"
   python3 cycle/scec_bp3_export.py "$OUT" "$PREFIX" "$OUT/scec" --modeler "Chunhui Zhao" \
      --element-size "25 m" --order 1 --scenario "$SCEN" --problem "SEAS Benchmark BP3-FD" \
      --version "$HEAD" \
      || echo "WARNING: cycle/scec_bp3_export.py failed; run it again by hand"
elif [ $rc -eq 8 ]; then
   echo "stopped at the wall limit; the next link continues from its checkpoint"
   rc=0
else
   echo "ERROR: driver exit $rc; the chain stops here"
fi
exit $rc
