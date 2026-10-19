#!/bin/bash
#SBATCH -J g0_accst
#SBATCH -o g0_accst_%j.out
#SBATCH -e g0_accst_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 04:00:00
#SBATCH -A EAR20006

# GPU plan G0: do the values the monitor and the writers take from the last
# operator call describe the accepted state?
#
#   sbatch [-n ranks] --export=ALL,ORDER=1 jobs/vista/g0_accepted_state.sbatch
#   sbatch -N n -n ranks --export=ALL,ORDER=1,RESTART_PREFIX=<dir/prefix>,MAX_STEPS=<step> \
#          jobs/vista/g0_accepted_state.sbatch
#   (submit from miniapps/seas)
#
# --check-accepted-state: after every accepted step, one fresh operator
# evaluation at the accepted (t, state) recomputes displacement, shear and
# normal traction, local and global V_max and the station slip rates, and
# prints [ACCEPTED-STATE] lines with the largest relative differences from
# the cached values (limit 1e-12) and a summary at exit.
#
# Without RESTART_PREFIX: 50 steps from t = 0 on the production mesh.
# With RESTART_PREFIX: continues a v63 checkpoint (written by the same rank
# count; -n must match it) up to accepted step MAX_STEPS (the TS step count
# is absolute after a restart), for the stretch where V_max crosses the
# output and event thresholds.

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g0_env.sh"
g0_check_binary

ORDER="${ORDER:?set ORDER=1 or 2}"
NRANKS="${SLURM_NTASKS}"
cd "${SEAS_DIR}"
if [ -n "${RESTART_PREFIX:-}" ]; then
   : "${MAX_STEPS:?set MAX_STEPS (absolute accepted-step count) with RESTART_PREFIX}"
   TAG="restart_p${ORDER}_r${NRANKS}_job${SLURM_JOB_ID}"
   EXTRA=(--restart "${RESTART_PREFIX}" --max-steps "${MAX_STEPS}")
else
   TAG="fresh_p${ORDER}_r${NRANKS}_job${SLURM_JOB_ID}"
   EXTRA=(--max-steps 50 --write-every-step)
fi
OUT="bp5/results_g0_vista/accepted_state/${TAG}"
mkdir -p "${OUT}"
g0_provenance "${OUT}"

ibrun ./seas_bp5_full --mesh "${G0_MESH_PROD}" "${G0_V63_FLAGS[@]}" --qd-operators loops \
      --order "${ORDER}" --check-accepted-state --checkpoint-interval 0 \
      "${EXTRA[@]}" --output-dir "${OUT}" --output-prefix "g0_accst" \
      > "${OUT}/run.log" 2>&1
rc=$?
echo "${rc}" > "${OUT}/exit_code"
echo "seas_bp5_full exit ${rc}"
grep '^\[ACCEPTED-STATE\]' "${OUT}/run.log" > "${OUT}/accepted_state.txt"
grep -c '^\[ACCEPTED-STATE\] step .*ABOVE' "${OUT}/accepted_state.txt" | sed 's/^/steps above 1e-12: /'
tail -n 3 "${OUT}/accepted_state.txt"
exit "${rc}"
