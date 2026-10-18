#!/bin/bash
#SBATCH --job-name=fdc_p4c
#SBATCH --clusters=condo
#SBATCH --partition=scec
#SBATCH --account=scec_608
#SBATCH --nodes=2 --ntasks-per-node=20 --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=48:00:00
#SBATCH --output=fdc_p4c_%j.out
#SBATCH --error=fdc_p4c_%j.err
# Fully dynamic BP5 cycle on Discovery, P4 (plan rev 3): continue a P4 run
# (discovery_p4.sbatch) from its last QD checkpoint to 400 yr
# (config/fd_cycle/bp5_fd_p4_400yr.toml, bp5_fd_p4.toml with tfinal 400 yr),
# so that its second event runs dynamic.  The run folder given as $1 is
# copied to runs/p4c_<job id> and the restart writes there; the original
# folder is not touched.  The restart truncates the copied files to the
# checkpoint and appends to them, so the copy holds the whole run.
# Usage, from miniapps/seas/jobs/fd_cycle, once the P4 job has ended (it
# writes its last checkpoint at its tfinal):
#   sbatch [--dependency=afterok:<P4 job id>] discovery_p4_continue.sbatch \
#          <runs>/p4_<P4 job id>
# Same 2 x 20 ranks as discovery_p4.sbatch (a restart needs the same rank
# count).  The binaries must be built from a commit whose code is HEAD's:
# git diff <built> HEAD may touch only config/fd_cycle and jobs/fd_cycle.
# Cost: 2 nodes x <= 48 h on the SCEC condo (scec_608).
# Submit only on Chunhui Zhao's typed go for this job.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MODULES="gcc/13.3.0 openmpi/5.0.5 openblas/0.3.28 hdf5/1.14.3 cmake/3.29.4 python/3.11.9"
SRC=${1:?usage: discovery_p4_continue.sbatch <P4 run folder>}
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

[ -f "$SRC/bp5_fd_p4_checkpoint_r0.txt" ] \
   || { echo "ERROR: no checkpoint $SRC/bp5_fd_p4_checkpoint_r0.txt"; exit 2; }
module purge
for m in $MODULES; do module load $m || exit 2; done
cd miniapps/seas
OUT=$(dirname "$ROOT")/runs/p4c_$SLURM_JOB_ID
[ ! -e "$OUT" ] || { echo "ERROR: $OUT exists"; exit 2; }
cp -a "$SRC" "$OUT" || exit 2
echo "copied $SRC to $OUT ($(du -sh "$OUT" | cut -f1))"
echo "host $(hostname)  nodes $SLURM_JOB_NODELIST  ranks $SLURM_NTASKS  output $OUT"
s=$(date +%s)
mpirun -np "$SLURM_NTASKS" --map-by ppr:20:node --bind-to core ./seas_fd_cycle_driver \
   --deck config/fd_cycle/bp5_fd_p4_400yr.toml --output-dir "$OUT" \
   --restart "$OUT/bp5_fd_p4" --log-every 1000 > "$OUT/run_continue.log" 2>&1
rc=$?
echo "seas_fd_cycle_driver exit $rc after $(( $(date +%s) - s )) s"
grep -E "restarted from" "$OUT/run_continue.log"
grep -A12 "dynamic phase [0-9]* summary" "$OUT/run_continue.log"
grep -E "QD -> dynamic at|back to QD at|STOP" "$OUT/run_continue.log"
for f in switch catalog; do echo "--- bp5_fd_p4_$f.txt"; cat "$OUT/bp5_fd_p4_$f.txt"; done
grep -A8 "=== fd_cycle summary" "$OUT/run_continue.log"
exit $rc
