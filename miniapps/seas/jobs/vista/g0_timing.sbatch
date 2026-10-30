#!/bin/bash
#SBATCH -J g0_timing
#SBATCH -o g0_timing_%j.out
#SBATCH -e g0_timing_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 12:00:00
#SBATCH -A EAR20006

# GPU plan G0: per-stage CPU time of the v63 deck on one Vista gg node.
#
#   sbatch --export=ALL,ORDER=P jobs/vista/g0_timing.sbatch
#   (submit from miniapps/seas; jobs/vista/g0_submit.sh does this)
#
# 1,000 accepted steps on the production mesh, 144 ranks, v63 flags, with
# --stage-timers 100: every 100 accepted steps and at exit the log has a
# [STAGE-TIMERS] table (slip RHS, Dirichlet RHS, solve, traction, friction,
# the monitor's every-step work, output and checkpoints, and the TS's own
# work).  Checkpoints every 250 steps so the table has checkpoint cost
# (production writes one every 5,000 steps).
#
# Differences from the v63 production job, on purpose:
#   * no --matrix-stats: it raises MUMPS's print level to 2, which prints
#     statistics on every solve and would be timed as solve; the factor
#     statistics for this layout come from g0_factor_size.sbatch CASE=prod_pP;
#   * no --paraview-dt: this build has no HDF5 (build_horizon.sh), so ParaView
#     output would be VTU rather than production's VTKHDF and its time would
#     not be production's;
#   * no --verify and no --ref-dir (start-up checks and an end-of-run
#     comparison, outside the step loop).

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g0_env.sh"
g0_check_binary

ORDER="${ORDER:?set ORDER=1 or 2}"
NRANKS="${SLURM_NTASKS:-144}"
cd "${SEAS_DIR}"
OUT="bp5/results_g0_vista/timing_p${ORDER}_r${NRANKS}_job${SLURM_JOB_ID}"
PFX="g0_timing_p${ORDER}"
mkdir -p "${OUT}"
g0_provenance "${OUT}"

ibrun ./seas_bp5_full --mesh "${G0_MESH_PROD}" "${G0_V63_FLAGS[@]}" --qd-operators loops \
      --order "${ORDER}" --max-steps 1000 --stage-timers 100 \
      --checkpoint-interval 250 \
      --output-dir "${OUT}" --output-prefix "${PFX}" > "${OUT}/run.log" 2>&1
rc=$?
echo "${rc}" > "${OUT}/exit_code"
echo "seas_bp5_full exit ${rc}"

# The end-of-run table, for the report.
sed -n '/^\[STAGE-TIMERS\] end of run/,$p' "${OUT}/run.log" \
   | grep '^\[STAGE-TIMERS\]' > "${OUT}/stage_timers_final.txt"
cat "${OUT}/stage_timers_final.txt"
exit "${rc}"
