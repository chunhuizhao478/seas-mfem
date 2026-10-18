#!/bin/bash
#SBATCH --job-name=fdc_eta
#SBATCH --clusters=condo
#SBATCH --partition=scec
#SBATCH --account=scec_608
#SBATCH --nodes=2 --ntasks-per-node=20 --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=24:00:00
#SBATCH --output=fdc_eta_%j.out
#SBATCH --error=fdc_eta_%j.err
# Fully dynamic BP5 cycle on Discovery, the P3 run started at 0.03 m/s
# (Chunhui Zhao, 2026-10-17: "let's submit a run with 0.03 m/s"): the first
# event in the dynamic regime from the scec_fd IC set with eta V_i in the
# prestress, so the nucleation patch starts at V_i = 0.03 m/s
# (config/fd_cycle/bp5_fd_eta.toml as committed; damped: Kelvin-Voigt 0.1 +
# Prakash-Clifton), the switch back to QD once max V has stayed below
# V_off = 1 mm/s for 10 s, then QD for 50 yr.  No switch probe.  It compares
# with the P3 run 6494469 (the patch from 0.015 m/s) and with the SBEM and
# BICyclE uploads.
# Expected (the cloud test, 2 ranks, to 8 s): at -24dp+10 the patch near
# 0.1 m/s at 1.7 s as BICyclE, then about 0.09 m/s to 6.3 s (BICyclE: 1 m/s
# at 8.3 s), so an earlier event than 6494469's but later than BICyclE's.
# Exit 0 at t_final, 7 if the QD balance right after the switch is at or
# above V_on (the guard), 6 at a later QD -> dynamic switch (P4).  A QD
# checkpoint is written at the switch and every 500 QD steps; if the job
# runs out of time in QD, continue with --restart <out>/bp5_fd_eta.
# Submit after discovery_build.sbatch, e.g.
#   sbatch --dependency=afterok:<build job id> discovery_eta.sbatch
# Cost: 2 nodes x <= 24 h on the SCEC condo (scec_608).
# Submit only on Chunhui Zhao's typed go for this job.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MODULES="gcc/13.3.0 openmpi/5.0.5 openblas/0.3.28 hdf5/1.14.3 cmake/3.29.4 python/3.11.9"
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
cd "$ROOT" || exit 2

# ---- GitHub-first guard (CLAUDE.md), and the binaries are from this commit ----
git fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
HEAD=$(git rev-parse HEAD)
BUILT=$(cat miniapps/seas/BUILD_COMMIT_DISCOVERY 2>/dev/null)
[ "$BUILT" = "$HEAD" ] || { echo "ERROR: binaries built at '$BUILT', HEAD is $HEAD; run discovery_build.sbatch"; exit 2; }
echo "commit $HEAD"

module purge
for m in $MODULES; do module load $m || exit 2; done
cd miniapps/seas
OUT=$(dirname "$ROOT")/runs/eta_$SLURM_JOB_ID
mkdir -p "$OUT" || exit 2
echo "host $(hostname)  nodes $SLURM_JOB_NODELIST  ranks $SLURM_NTASKS  output $OUT"
s=$(date +%s)
mpirun -np "$SLURM_NTASKS" --map-by ppr:20:node --bind-to core ./seas_fd_cycle_driver \
   --deck config/fd_cycle/bp5_fd_eta.toml --output-dir "$OUT" --log-every 1000 > "$OUT/run.log" 2>&1
rc=$?
echo "seas_fd_cycle_driver exit $rc after $(( $(date +%s) - s )) s"
grep -A12 "dynamic phase 1 summary" "$OUT/run.log"
grep -E "IC set|initial max V|back to QD at|STOP" "$OUT/run.log"
for f in switch catalog; do echo "--- bp5_fd_eta_$f.txt"; cat "$OUT/bp5_fd_eta_$f.txt"; done
grep -A8 "=== fd_cycle summary" "$OUT/run.log"
exit $rc
