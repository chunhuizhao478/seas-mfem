#!/bin/bash
#SBATCH --job-name=fdc_p4
#SBATCH --clusters=condo
#SBATCH --partition=scec
#SBATCH --account=scec_608
#SBATCH --nodes=2 --ntasks-per-node=20 --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=48:00:00
#SBATCH --output=fdc_p4_%j.out
#SBATCH --error=fdc_p4_%j.err
# Fully dynamic BP5 cycle on Discovery, P4 (plan rev 3): the first two
# switched cycles, DYN -> QD -> DYN -> QD (config/fd_cycle/bp5_fd_p4.toml as
# committed).  The first event dynamic from t = 0 at V_i = 0.03 m/s (eta V_i
# in the prestress, as the 0.03 m/s run), the switch back to QD once max V
# has stayed below V_off = 1 mm/s for 10 s, QD until max V passes V_on =
# 1 cm/s again (the second event, about 225 yr in the QD reference), that
# event dynamic, back to QD, and QD to 300 yr.  Damped dynamic regime
# (Kelvin-Voigt 0.1 + Prakash-Clifton).
# Expected: two dynamic phases in bp5_fd_p4_switch.txt, two closed events in
# the catalog, fault tables bp5_fd_p4_fault.txt and _fault_phase2.txt.
# Exit 0 at t_final; 7 if the QD balance right after a switch back is at or
# above V_on (the guard); 4 if a dynamic phase passes T_max = 1200 s.  A QD
# checkpoint is written at each switch back and every 500 QD steps; if the
# job runs out of time in QD, continue with --restart <out>/bp5_fd_p4.
# Pace (the P3 run 6494469, 40 ranks): 0.036 s per dynamic step, 1.3-1.4 s
# per QD step; the QD reference takes about 48,500 steps per cycle, so the
# 300 yr fit in 48 h (inferred).
# discovery_build.sbatch replaces the binaries, so it runs only when no
# other job uses this checkout; then
#   sbatch --dependency=afterok:<build job id> discovery_p4.sbatch
# Cost: 2 nodes x <= 48 h on the SCEC condo (scec_608).
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
OUT=$(dirname "$ROOT")/runs/p4_$SLURM_JOB_ID
mkdir -p "$OUT" || exit 2
echo "host $(hostname)  nodes $SLURM_JOB_NODELIST  ranks $SLURM_NTASKS  output $OUT"
s=$(date +%s)
mpirun -np "$SLURM_NTASKS" --map-by ppr:20:node --bind-to core ./seas_fd_cycle_driver \
   --deck config/fd_cycle/bp5_fd_p4.toml --output-dir "$OUT" --log-every 1000 > "$OUT/run.log" 2>&1
rc=$?
echo "seas_fd_cycle_driver exit $rc after $(( $(date +%s) - s )) s"
grep -A12 "dynamic phase [0-9]* summary" "$OUT/run.log"
grep -E "IC set|initial max V|QD -> dynamic at|back to QD at|STOP" "$OUT/run.log"
for f in switch catalog; do echo "--- bp5_fd_p4_$f.txt"; cat "$OUT/bp5_fd_p4_$f.txt"; done
grep -A8 "=== fd_cycle summary" "$OUT/run.log"
exit $rc
