#!/bin/bash
#SBATCH -J scec_bp5fd_p2
#SBATCH -o scec_bp5fd_p2_%j.out
#SBATCH -e scec_bp5fd_p2_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 48:00:00
#SBATCH -A EAR24034
# SCEC BP5-FD to 1800 yr at p = 2 on the 1000 m mesh, rerun from the
# scec-benchmarks branch for the SCEC SEAS upload (Chunhui Zhao,
# 2026-11-21: "let's submit the full BP5-FD from this branch and directly
# run and test it").  The deck is the BP5-FD p = 2 gold's,
# config/fd_cycle/bp5_fd_p2_1800yr.toml, unchanged (gold: Vista 1094927 +
# 1094928 at 3926e9c8, bp5/benchmark_data/golden_fd_vista_576r_p2_1000m);
# the run tests that this branch reproduces the gold.
# vista_qd_p2_1800yr.sbatch (the BP5-QD run) with three changes: the FD
# deck and prefix; the run folder $SCRATCH/scec_bp5/fd_p2_1800_<job id> and
# the H-matrix file $SCRATCH/scec_bp5/bp5_fd_p2_h1e-10.bin (results only in
# scratch, project rule 2026-11-15; the gold's vista_p2_1800yr.sbatch wrote
# to <checkout>-runs); the export with --problem "SEAS Benchmark BP5-FD".
# Same 4 x 144 layout, 45 h wall limit per link, restart in place:
#   first link (no argument):  a fresh run in $SCRATCH/scec_bp5/fd_p2_1800_<job id>
#   later links ($1 = that folder): --restart from its last checkpoint
# Chain, from miniapps/seas/jobs/fd_cycle:
#   j1=$(sbatch --parsable vista_fd_p2_1800yr_scec.sbatch)
#   j2=$(sbatch --parsable --dependency=afterok:$j1 vista_fd_p2_1800yr_scec.sbatch $SCRATCH/scec_bp5/fd_p2_1800_$j1)
#   j3=$(sbatch --parsable --dependency=afterok:$j2 vista_fd_p2_1800yr_scec.sbatch $SCRATCH/scec_bp5/fd_p2_1800_$j1)
# (TACC's sbatch prints a banner before the id: take its last line.)
# The H-matrix (eps 1e-10, about 38 GiB): the first link builds and writes
# it (about 80 min, from 1086777); a later link reads it.
# Expected: the gold took 1094927 (45 h, to 1575 yr) + 1094928 (to 1800 yr),
# so two links; the third is a spare (a link whose folder has DONE exits at
# once).  Peak committed memory in the gold: 109.5-120.2 GiB per node.
# Binaries: vista_build.sbatch at this code (git diff <built> HEAD may touch
# only config/fd_cycle and jobs/fd_cycle).
# Cost: 4 nodes x <= 48 h per link on gg = <= 64 SU per link (EAR24034).
# Submit only on Chunhui Zhao's typed go for these jobs.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NVPL_BLAS_NUM_THREADS=1 NVPL_LAPACK_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
MODULES="gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5"
RUN=${1:-}
DECK=${FDC_DECK:-config/fd_cycle/bp5_fd_p2_1800yr.toml}
PREFIX=${FDC_PREFIX:-bp5_fd_p2_1800yr}
WALL=${FDC_WALL_S:-162000}                                # 45 h
HCACHE=${FDC_HCACHE:-$SCRATCH/scec_bp5/bp5_fd_p2_h1e-10.bin}
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
   OUT=$SCRATCH/scec_bp5/fd_p2_1800_$SLURM_JOB_ID
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
srun --overlap -N "$SLURM_NNODES" -n "$SLURM_NNODES" --ntasks-per-node=1 \
   "$ROOT/miniapps/seas/jobs/fd_cycle/vista_memlog.sh" 60 > "$OUT/memlog_$SLURM_JOB_ID.txt" 2>&1 &
MEMLOG=$!
trap 'kill $MEMLOG 2>/dev/null' EXIT
s=$(date +%s)
ibrun ./seas_fd_cycle_driver \
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
# The outcome from the log as well as the exit code, which comes through
# ibrun and mpirun: a wall-limit stop is exit 8 whatever code arrives, and
# exit 0 counts as the end only when the driver's summary says "status:
# done", which it prints when the QD run reaches t_final (this job passes
# no --max-steps).  The summary's final time has 3 digits only, so it is
# shown, not compared with t_final (1094928 reached t_final and failed
# that comparison: exit 9, no DONE, no SCEC files).
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
   python3 cycle/scec_bp5_export.py "$OUT" "$PREFIX" "$OUT/scec" --modeler "Chunhui Zhao" \
      --problem "SEAS Benchmark BP5-FD" \
      --element-size "1000 m" --order 2 --version "$HEAD" \
      || echo "WARNING: cycle/scec_bp5_export.py failed; run it again by hand"
elif [ $rc -eq 8 ]; then
   echo "stopped at the wall limit; the next link continues from its checkpoint"
   rc=0
else
   echo "ERROR: driver exit $rc; the chain stops here"
fi
exit $rc
