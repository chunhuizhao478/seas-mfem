#!/bin/bash
#SBATCH -J h_gpu_gates
#SBATCH -o h_gpu_gates_%j.out
#SBATCH -e h_gpu_gates_%j.err
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 08:00:00
#SBATCH -A ear26024

# GPU plan on Horizon: the G3a, G4 and G3c gates of jobs/vista
# (g3_gpu_gate, g4_gpu_gate, g3c_gpu_gate), one after the other on one
# GB200 node (1- and 2-GPU runs on GPUs 0 and 1), with the settings of
# jobs/horizon/h_env.sh.  Each gate rebuilds seas_bp5_full, compares its
# 50-step runs with Vista's G0 goldens at the tolerances in their
# floor.json, and writes bp5/goldens/<gate>/job<id>/.  The G3c gate needs
# GREEN_CACHE (h_green_cache.sbatch).  The job fails if any gate fails.
#
#   sbatch jobs/horizon/h_gpu_gates.sbatch   (from miniapps/seas of the Horizon clone)

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/horizon/h_env.sh" || exit 1
overall=0
for GATE in g3_gpu_gate g4_gpu_gate g3c_gpu_gate; do
   echo "===== ${GATE}: start $(date -u +%T) ====="
   bash "${SLURM_SUBMIT_DIR}/jobs/vista/${GATE}.sbatch"
   rc=$?
   echo "===== ${GATE}: exit ${rc} at $(date -u +%T) ====="
   [ "${rc}" = 0 ] || overall=1
done
echo "Horizon GPU gates overall: $([ "${overall}" = 0 ] && echo PASS || echo FAIL)"
exit "${overall}"
