#!/bin/bash
#SBATCH -J g5_gpu_cycle
#SBATCH -o g5_gpu_cycle_%j.out
#SBATCH -e g5_gpu_cycle_%j.err
#SBATCH -p gh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# GPU plan G5: the full BP5-QD earthquake cycle (default 1500 yr) on one
# Vista Grace-Hopper (gh) node, one GH200, in the CUDA clone
# ($WORK/seas-mfem-gpu-cuda, built by g1_cuda_build.sbatch).  Same deck as
# the CPU gold 1074099 (jobs/vista/bp5_full_cycle.sbatch): v63 flags,
# production mesh, p=1, t_final 1500 yr, checkpoints every 5,000 accepted
# steps, station output on the adaptive SCEC schedule.  The stage timers
# are on, so the run also measures the GPU speed.
#
#   sbatch [--export=ALL,QD_PATH=cudss,ORDER=1,TFINAL_YR=1500] \
#          jobs/vista/g5_gpu_full_cycle.sbatch      (from miniapps/seas)
#
# QD_PATH picks the stage solve:
#   cudss  (default) K u = b solved per stage with cuDSS on the GPU (G3a),
#          stage operators (--qd-operators matrix), friction and time
#          stepping on the GPU (G4);
#   green  opt-in: the dense Green's-function stage on the GPU (G3c,
#          --qd-green).  G is read from GREEN_CACHE (default: the QD speed
#          thread's file $SCRATCH/seas-qd-green/bp5_prod_p1.bin, production
#          mesh, p=1); cuDSS still factors K once for the startup check.
#          Linear elastic bulk only.
# Submit it only after the 50-step gate of that path passed at this commit
# (g4_gpu_gate.sbatch for cudss, g3c_gpu_gate.sbatch for green).
#
# Speed: speed.txt in the output folder gives the simulated years, accepted
# steps, wall time and s/step, next to the CPU gold 1074099 (1 gg node, 144
# ranks, binary a4e6f90: 313,960 steps in 31.6 h, 0.36 s/step).  The
# [STAGE-TIMERS] table, every 5,000 accepted steps and at exit, splits the
# step into its parts.
#
# Compare the result with bp5/benchmark_data/golden_vista_144r_p1_1000m (the
# CPU gold 1074099) and golden_tandem_p1_1000m (Tandem p=1, same mesh), as
# in golden_vista_144r_p1_1000m/compare_tandem_p1.  No --verify / --ref-dir
# (see bp5_full_cycle.sbatch); no --paraview-dt.
#
# GitHub is the version control: the clone must be clean and at the commit
# its upstream branch points to, or the job stops.

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g1_cuda_env.sh"
g1_check_checkout

QD_PATH="${QD_PATH:-cudss}"
ORDER="${ORDER:-1}"
TFINAL_YR="${TFINAL_YR:-1500}"
TFINAL_S="$(awk -v y="${TFINAL_YR}" 'BEGIN { printf "%.6e", y * 365.25 * 24 * 3600 }')"
NP="${SLURM_NTASKS:-1}"

case "${QD_PATH}" in
   cudss)
      PATH_FLAGS=(--qd-operators matrix)
      ;;
   green)
      GREEN_CACHE="${GREEN_CACHE:-${SCRATCH:-${WORK}}/seas-qd-green/bp5_prod_p1.bin}"
      if [ ! -s "${GREEN_CACHE}" ]; then
         echo "ERROR: no Green's-function cache file ${GREEN_CACHE}: run the CPU" >&2
         echo "       build first (jobs/vista_qd_speed, step green, g_n1_r144)." >&2
         exit 1
      fi
      echo "Green cache: ${GREEN_CACHE} ($(stat -c '%s bytes, modified %y' "${GREEN_CACHE}"))"
      PATH_FLAGS=(--qd-operators matrix --qd-green --qd-green-cache "${GREEN_CACHE}")
      ;;
   *)
      echo "ERROR: QD_PATH must be cudss or green, not '${QD_PATH}'." >&2
      exit 1
      ;;
esac

# Rebuilds seas_bp5_full if this commit changed it (MFEM, PETSc and cuDSS
# are the G1 build's, unchanged); after the gate at this commit it is a no-op.
MODE=cuda GPU_ARCH=90 STAGE="seas" JOBS=72 MODULES=none bash "${SEAS_MFEM_ROOT}/build_horizon.sh" || exit 1
g0_check_binary

cd "${SEAS_DIR}"
OUT="bp5/results_vista_cycle/gpu_${QD_PATH}_p${ORDER}_${TFINAL_YR}yr_g${NP}_job${SLURM_JOB_ID}"
PFX="bp5_cycle_p${ORDER}"
mkdir -p "${OUT}"
g0_provenance "${OUT}"
[ "${QD_PATH}" = green ] && echo "green cache: ${GREEN_CACHE}" >> "${OUT}/provenance.txt"
echo "BP5 full cycle on ${NP} GPU(s), QD_PATH ${QD_PATH}: p=${ORDER}, t_final ${TFINAL_YR} yr (${TFINAL_S} s) -> ${OUT}"

T0="$(date +%s)"
ibrun -n "${NP}" -o 0 ./seas_bp5_full --mesh "${G0_MESH_PROD}" "${G0_V63_FLAGS[@]}" \
      --order "${ORDER}" --device cuda --solver cudss "${PATH_FLAGS[@]}" \
      --stage-timers 5000 --tfinal "${TFINAL_S}" \
      --checkpoint-interval 5000 \
      --output-dir "${OUT}" --output-prefix "${PFX}" > "${OUT}/run.log" 2>&1
rc=$?
T1="$(date +%s)"
echo "${rc}" > "${OUT}/exit_code"
echo "seas_bp5_full exit ${rc}"

sed -n '/^\[STAGE-TIMERS\] end of run/,$p' "${OUT}/run.log" \
   | grep '^\[STAGE-TIMERS\]' > "${OUT}/stage_timers_final.txt"
if [ "${QD_PATH}" = green ]; then
   grep '^\[QD-GREEN\]' "${OUT}/run.log" > "${OUT}/qd_green.txt"
fi
STEPS="$(sed -n 's/^ *Total steps: *\([0-9]*\).*/\1/p' "${OUT}/run.log" | tail -n 1)"
YEARS="$(sed -n 's/^ *Final time: *\([0-9.eE+-]*\) years.*/\1/p' "${OUT}/run.log" | tail -n 1)"
EVENTS="$(sed -n 's/^ *Seismic events: *\([0-9]*\).*/\1/p' "${OUT}/run.log" | tail -n 1)"
WALL=$((T1 - T0))
awk -v steps="${STEPS:-0}" -v wall="${WALL}" -v years="${YEARS:-?}" \
    -v events="${EVENTS:-?}" -v path="${QD_PATH}" -v np="${NP}" \
    -v job="${SLURM_JOB_ID:-none}" -v rc="${rc}" 'BEGIN {
   printf "job %s, QD_PATH %s, %s GPU(s), exit %s\n", job, path, np, rc
   printf "simulated: %s yr, %s seismic events\n", years, events
   printf "accepted steps: %d\n", steps
   printf "wall (seas_bp5_full, setup included): %d s = %.2f h\n", wall, wall / 3600
   if (steps > 0) {
      printf "s/step: %.4f\n", wall / steps
      printf "vs CPU gold 1074099 (0.36 s/step, 31.6 h for 1500 yr): %.2fx per step\n", 0.362 / (wall / steps)
   }
}' | tee "${OUT}/speed.txt"
cat "${OUT}/stage_timers_final.txt"
tail -n 20 "${OUT}/run.log"
exit "${rc}"
