#!/bin/bash
#SBATCH -J bp5_cycle
#SBATCH -o bp5_cycle_%j.out
#SBATCH -e bp5_cycle_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# Full BP5-QD earthquake cycle (default 1500 yr) on one Vista gg node, with
# the G0 build (jobs/vista/g0_build.sbatch) and the v63 production flags.
#
#   sbatch [--export=ALL,ORDER=1,TFINAL_YR=1500] jobs/vista/bp5_full_cycle.sbatch
#   (submit from miniapps/seas)
#
# Station output uses the adaptive SCEC schedule (no --write-every-step).
# Checkpoints every 5,000 accepted steps, as in production.  Compare the
# result with bp5/benchmark_data/golden_tandem_p1_1000m (Tandem p=1, same
# mesh) and with golden_parallel_400r_p1_1000m (Frontera, written before the
# station sigma_n fix ef70517, so its station V and tau differ near
# nucleation; slip, state and event times should agree).
#
# No --verify / --ref-dir: the built-in comparison expects SCEC station file
# names, which the Tandem golden does not use.  No --paraview-dt: this build
# has no HDF5.

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g0_env.sh"
g0_check_binary

ORDER="${ORDER:-1}"
TFINAL_YR="${TFINAL_YR:-1500}"
TFINAL_S="$(awk -v y="${TFINAL_YR}" 'BEGIN { printf "%.6e", y * 365.25 * 24 * 3600 }')"
NRANKS="${SLURM_NTASKS:-144}"
cd "${SEAS_DIR}"
OUT="bp5/results_vista_cycle/p${ORDER}_${TFINAL_YR}yr_r${NRANKS}_job${SLURM_JOB_ID}"
PFX="bp5_cycle_p${ORDER}"
mkdir -p "${OUT}"
g0_provenance "${OUT}"
echo "BP5 full cycle: p=${ORDER}, t_final ${TFINAL_YR} yr (${TFINAL_S} s), ${NRANKS} ranks -> ${OUT}"

ibrun ./seas_bp5_full --mesh "${G0_MESH_PROD}" "${G0_V63_FLAGS[@]}" --qd-operators loops \
      --order "${ORDER}" --tfinal "${TFINAL_S}" \
      --checkpoint-interval 5000 \
      --output-dir "${OUT}" --output-prefix "${PFX}" > "${OUT}/run.log" 2>&1
rc=$?
echo "${rc}" > "${OUT}/exit_code"
echo "seas_bp5_full exit ${rc}"
tail -n 20 "${OUT}/run.log"
exit "${rc}"
