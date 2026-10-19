#!/bin/bash
#SBATCH --job-name=fdc_e1p2
#SBATCH --clusters=condo
#SBATCH --partition=scec
#SBATCH --account=scec_608
#SBATCH --nodes=4 --ntasks-per-node=20 --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=48:00:00
#SBATCH --output=fdc_e1p2_%j.out
#SBATCH --error=fdc_e1p2_%j.err
# Fully dynamic BP5 cycle on Discovery: P4's first event at p = 2 on the same
# mesh (config/fd_cycle/bp5_fd_event1_p2.toml, bp5_fd_eta.toml with order 2
# and max_dynamic_phases = 1).  Dynamic from t = 0 at V_i = 0.03 m/s, damped
# (Kelvin-Voigt 0.1 + Prakash-Clifton), until max V has stayed below V_off =
# 1 mm/s for 10 s; then the run stops with no QD step.  Compare slip, tau and
# the sigma_n_eff range with the p = 1 run 6500504 and with BICyclE and SBEM.
# Expected: one dynamic phase in bp5_fd_event1_p2_switch.txt, one closed
# event in the catalog, fault table bp5_fd_event1_p2_fault.txt; exit 0.
# 4 x 20 ranks, twice P4's, for memory: p = 2 has 2.5x the unknowns, and a
# 2-rank cloud setup passed 12.5 GB before assembly ended.  p = 2 has not run
# in this driver before, so the first log lines are its smoke test.
# Pace (estimate, not measured): dt 0.81 ms (cfl 0.5 x C(2) 0.212, against
# 1.72 ms at p = 1), 3-5x the p = 1 cost per step (0.0335 s on 40 ranks), so
# 270 s of event takes 4.5-8 h on 80 ranks.  The first log lines give dt
# and s/step.
# The binaries must be built from a commit whose code is HEAD's:
# git diff <built> HEAD may touch only config/fd_cycle and jobs/fd_cycle.
# Usage, from miniapps/seas/jobs/fd_cycle:
#   sbatch discovery_event1_p2.sbatch
# Cost: 4 nodes x <= 48 h on the SCEC condo (scec_608), expected 5-8 h.
# Submit only on Chunhui Zhao's typed go for this job.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MODULES="gcc/13.3.0 openmpi/5.0.5 openblas/0.3.28 hdf5/1.14.3 cmake/3.29.4 python/3.11.9"
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

module purge
for m in $MODULES; do module load $m || exit 2; done
cd miniapps/seas
OUT=$(dirname "$ROOT")/runs/event1_p2_$SLURM_JOB_ID
mkdir -p "$OUT" || exit 2
echo "host $(hostname)  nodes $SLURM_JOB_NODELIST  ranks $SLURM_NTASKS  output $OUT"
s=$(date +%s)
mpirun -np "$SLURM_NTASKS" --map-by ppr:20:node --bind-to core ./seas_fd_cycle_driver \
   --deck config/fd_cycle/bp5_fd_event1_p2.toml --output-dir "$OUT" --log-every 1000 > "$OUT/run.log" 2>&1
rc=$?
echo "seas_fd_cycle_driver exit $rc after $(( $(date +%s) - s )) s"
grep -A12 "dynamic phase [0-9]* summary" "$OUT/run.log"
grep -E "IC set|initial max V|dt = |QD -> dynamic at|back to QD at|STOP" "$OUT/run.log"
for f in switch catalog; do echo "--- bp5_fd_event1_p2_$f.txt"; cat "$OUT/bp5_fd_event1_p2_$f.txt"; done
grep -A8 "=== fd_cycle summary" "$OUT/run.log"
exit $rc
