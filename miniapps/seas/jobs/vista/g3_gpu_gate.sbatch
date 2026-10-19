#!/bin/bash
#SBATCH -J g3_gpu_gate
#SBATCH -o g3_gpu_gate_%j.out
#SBATCH -e g3_gpu_gate_%j.err
#SBATCH -p gh
#SBATCH -N 2
#SBATCH -n 2
#SBATCH -t 04:00:00
#SBATCH -A EAR20006

# GPU plan G3a on Vista Grace-Hopper (gh) nodes, one GH200 per node, in the
# CUDA clone ($WORK/seas-mfem-gpu-cuda, built by g1_cuda_build.sbatch):
#
#   sbatch jobs/vista/g3_gpu_gate.sbatch   (from miniapps/seas of the CUDA clone)
#
# Rebuilds seas_bp5_full and seas_test_cudss_vs_mumps (MFEM, PETSc and cuDSS
# are the G1 build's, unchanged), then
#   1. seas_test_cudss_vs_mumps on the production mesh at p=1, on 1 GPU and
#      on 2 GPUs (two ranks, cuDSS's row distribution): one factor for 10
#      stage right-hand sides, each cuDSS solution's relative residual and
#      its difference from MUMPS's below 1e-10; it prints the setup times,
#      the device memory cuDSS took, cuDSS's analysis estimates and the
#      time per solve;
#   2. the 50-step gate with --device cuda --solver cudss --qd-operators
#      matrix at p=1 on 1 and 2 GPUs and at p=2 on 1 GPU against G0's
#      goldens (G1_GOLDEN_ROOT/pP_rN/run_a) at TOL = floor.json key
#      G3_to_G5_gpu, judged by check_gate.sh, with --stage-timers 10: the
#      per-stage time table (stage_timers_final.txt, next to G0's) and the
#      host loading share, g(t) callback + copy against the stage time
#      (threshold 5%, a LOADING line);
#   3. for comparison, not a gate: p=1 on 1 GPU with the stage operators
#      on the GPU and MUMPS on the Grace CPU (--solver mumps), the same
#      gate and tables;
#   4. cuDSS's analysis estimates alone (--estimate-only) on 1 GPU for the
#      500 m mesh at p=1 and p=2 and the 250 m mesh at p=1.  250 m at p=2
#      needs at least 4 ranks: g3_cudss_estimate.sbatch.
#
# GitHub is the version control: the clone must be clean and at the commit
# its upstream branch points to, or the job stops.

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g1_cuda_env.sh"
g1_check_checkout

MODE=cuda GPU_ARCH=90 STAGE="seas" JOBS=72 MODULES=none bash "${SEAS_MFEM_ROOT}/build_horizon.sh" || exit 1
g0_check_binary
( cd "${SEAS_DIR}" && rm -f tests/parallel/test_cudss_vs_mumps.o &&
  make -j8 seas_test_cudss_vs_mumps ) || exit 1

cd "${SEAS_DIR}"
OUT_ROOT="bp5/goldens/g3_gpu_gate/job${SLURM_JOB_ID}"
mkdir -p "${OUT_ROOT}"
g0_provenance "${OUT_ROOT}"
overall=0

# ---- 1. cuDSS against MUMPS ----
mkdir -p "${OUT_ROOT}/cudss_vs_mumps"
for NP in 1 2; do
   LOG="${OUT_ROOT}/cudss_vs_mumps/p1_r${NP}.log"
   echo "=== cuDSS vs MUMPS p=1 on ${NP} GPU(s): start $(date -u +%T) ==="
   ibrun -n "${NP}" -o 0 ./seas_test_cudss_vs_mumps --mesh "${G0_MESH_PROD}" \
         --order 1 --device cuda > "${LOG}" 2>&1
   rc=$?
   grep -E "QD-SOLVER|cuDSS analysis|^  [0-9]|^  rank|solve time|PASS|FAIL" "${LOG}"
   if [ "${rc}" = 0 ] && grep -q "ALL PASSED" "${LOG}"; then
      echo "CUDSS_VS_MUMPS p1_r${NP}: PASS"
   else
      echo "CUDSS_VS_MUMPS p1_r${NP}: FAIL (exit ${rc})"
      overall=1
   fi
done

# ---- 2, 3. 50-step gates ----
for CASE in 1:1:cudss 1:2:cudss 2:1:cudss 1:1:mumps; do
   IFS=: read -r ORDER NP SOLVER <<< "${CASE}"
   TAG="p${ORDER}_r${NP}_${SOLVER}"
   GOLD="${G1_GOLDEN_ROOT}/p${ORDER}_r${NP}"
   GPFX="g0_p${ORDER}_r${NP}"
   OUT="${OUT_ROOT}/${TAG}"
   PFX="g3_${TAG}"
   if [ ! -f "${GOLD}/floor.json" ] || [ ! -d "${GOLD}/run_a" ]; then
      echo "ERROR: no G0 golden in ${GOLD}" >&2
      overall=1
      continue
   fi
   TOL="$(python3 -c "import json; print(json.load(open('${GOLD}/floor.json'))['tol']['G3_to_G5_gpu'])")"
   mkdir -p "${OUT}"
   echo "=== ${TAG} on ${NP} GPU(s): TOL ${TOL}, start $(date -u +%T) ==="
   ibrun -n "${NP}" -o 0 ./seas_bp5_full --mesh "${G0_MESH_PROD}" "${G0_V63_FLAGS[@]}" \
         --order "${ORDER}" --device cuda --qd-operators matrix \
         --solver "${SOLVER}" --stage-timers 10 \
         --max-steps 50 --write-every-step --checkpoint-interval 0 \
         --output-dir "${OUT}" --output-prefix "${PFX}" \
         --ref-dir "${GOLD}/run_a" --ref-prefix "${GPFX}" \
         --regression-tolerance "${TOL}" > "${OUT}/run.log" 2>&1
   rc=$?
   echo "${rc}" > "${OUT}/exit_code"
   echo "=== ${TAG}: exit ${rc} at $(date -u +%T) ==="
   sed -n '/^\[STAGE-TIMERS\] end of run/,$p' "${OUT}/run.log" \
      | grep '^\[STAGE-TIMERS\]' > "${OUT}/stage_timers_final.txt"
   { echo "--- ${TAG}"; grep -E "\[QD-SOLVER\]|\[QD-OPERATORS\]" "${OUT}/run.log"; } \
      >> "${OUT_ROOT}/qd_operators.txt"
   cat "${OUT}/stage_timers_final.txt"
   if bash jobs/horizon/check_gate.sh "${OUT}" "${PFX}" "${GOLD}/run_a" "${GPFX}" \
         "${TOL}" "${OUT}/run.log" "${rc}"; then
      echo "GATE ${TAG}: PASS"
   else
      echo "GATE ${TAG}: FAIL"
      [ "${SOLVER}" = cudss ] && overall=1
   fi
   share="$(sed -n 's/.*host loading (g(t) callback + copy): \([0-9.]*\)%.*/\1/p' "${OUT}/run.log" | tail -n 1)"
   if [ -z "${share}" ]; then
      echo "LOADING ${TAG}: no host loading line"
      [ "${SOLVER}" = cudss ] && overall=1
   elif python3 -c "import sys; sys.exit(0 if float('${share}') <= 5.0 else 1)"; then
      echo "LOADING ${TAG}: ${share}% of the stage time (<= 5%)"
   else
      echo "LOADING ${TAG}: ${share}% of the stage time, ABOVE the 5% threshold"
      [ "${SOLVER}" = cudss ] && overall=1
   fi
done

# ---- 4. cuDSS analysis estimates, 500 m and 250 m ----
mkdir -p "${OUT_ROOT}/estimates"
for CASE in 500m:1 500m:2 250m:1; do
   IFS=: read -r RES ORDER <<< "${CASE}"
   case "${RES}" in
      500m) MESH="${G0_MESH_500M}" ;;
      250m) MESH="${G0_MESH_250M}" ;;
   esac
   LOG="${OUT_ROOT}/estimates/${RES}_p${ORDER}_r1.log"
   echo "=== cuDSS estimate ${RES} p=${ORDER} on 1 GPU: start $(date -u +%T) ==="
   ibrun -n 1 -o 0 ./seas_test_cudss_vs_mumps --mesh "${MESH}" --order "${ORDER}" \
         --device cuda --estimate-only > "${LOG}" 2>&1
   rc=$?
   grep -E "^G3a|cuDSS analysis|^  rank|^  [0-9]" "${LOG}"
   echo "ESTIMATE ${RES}_p${ORDER}_r1: exit ${rc} (a measurement, not a gate)"
done

cat "${OUT_ROOT}/qd_operators.txt" 2>/dev/null
echo "G3a GPU gate overall: $([ "${overall}" = 0 ] && echo PASS || echo FAIL)"
exit "${overall}"
