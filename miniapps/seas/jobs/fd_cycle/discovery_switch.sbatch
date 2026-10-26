#!/bin/bash
#SBATCH --job-name=fdc_p3_switch
#SBATCH --clusters=condo
#SBATCH --partition=scec
#SBATCH --account=scec_608
#SBATCH --nodes=2 --ntasks-per-node=20 --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=24:00:00
#SBATCH --output=fdc_p3_switch_%j.out
#SBATCH --error=fdc_p3_switch_%j.err
# Fully dynamic BP5 cycle, plan phase P3 on Discovery (Chunhui Zhao,
# 2026-10-16): the first event in the dynamic regime from the scec_fd IC set
# (damped: Kelvin-Voigt 0.1 + Prakash-Clifton), the switch back to QD once max
# V has stayed below V_off = 1 mm/s for 10 s, then QD for 50 yr
# (config/fd_cycle/bp5_fd_p3.toml as committed).  The switch probe measures
# the tau and V jump at the stations for V_off {1, 3} mm/s x T_dwell
# {0.1, 1, 10} s; the deck's own pair (1 mm/s, 10 s) switches last, so the
# dynamic phase ends at the deck's switch.
# Expected: the switch near t = 288.6 s (the Mac P2 run of the same event);
# exit 0 at t_final, 7 if the QD balance right after the switch is at or above
# V_on (the guard); 6 at a later QD -> dynamic switch in builds before P4
# (not before ~225 yr; from P4 on that event runs dynamic).
# A QD checkpoint is written at the switch and every 500 QD steps; if the job
# runs out of time in QD, continue with --restart <out>/bp5_fd_p3.
# Submit after discovery_build.sbatch, e.g.
#   sbatch --dependency=afterok:<build job id> discovery_switch.sbatch
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
OUT=$(dirname "$ROOT")/runs/p3_switch_$SLURM_JOB_ID
mkdir -p "$OUT" || exit 2
echo "host $(hostname)  nodes $SLURM_JOB_NODELIST  ranks $SLURM_NTASKS  output $OUT"
s=$(date +%s)
mpirun -np "$SLURM_NTASKS" --map-by ppr:20:node --bind-to core ./seas_fd_cycle_driver \
   --deck config/fd_cycle/bp5_fd_p3.toml --output-dir "$OUT" --log-every 1000 > "$OUT/run.log" 2>&1
rc=$?
echo "seas_fd_cycle_driver exit $rc after $(( $(date +%s) - s )) s"
grep -A12 "dynamic phase 1 summary" "$OUT/run.log"
grep -E "back to QD at|STOP" "$OUT/run.log"
for f in switch switch_probe_summary catalog; do echo "--- bp5_fd_p3_$f.txt"; cat "$OUT/bp5_fd_p3_$f.txt"; done
grep -A8 "=== fd_cycle summary" "$OUT/run.log"
exit $rc
