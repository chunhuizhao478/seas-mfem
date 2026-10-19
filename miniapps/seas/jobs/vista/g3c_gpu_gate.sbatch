#!/bin/bash
#SBATCH -J g3c_gpu_gate
#SBATCH -o g3c_gpu_gate_%j.out
#SBATCH -e g3c_gpu_gate_%j.err
#SBATCH -p gh
#SBATCH -N 2
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# GPU plan G3c (opt-in dense Green's-function stage, --qd-green) on Vista
# Grace-Hopper (gh) nodes, one GH200 per node, in the CUDA clone
# ($WORK/seas-mfem-gpu-cuda, built by g1_cuda_build.sbatch):
#
#   sbatch jobs/vista/g3c_gpu_gate.sbatch   (from miniapps/seas of the CUDA clone)
#
# G is not built here: it is read from GREEN_CACHE, the cache file a CPU
# run wrote (default: the QD speed thread's dense-G timing job g_n1_r144,
# $SCRATCH/seas-qd-green/bp5_prod_p1.bin, production mesh, p=1).  The job
# stops if the file is missing.  Rebuilds seas_bp5_full and
# seas_test_qd_green_device (MFEM, PETSc and cuDSS are the G1 build's,
# unchanged), then
#   1. seas_test_qd_green_device --device cuda at the production size
#      (N_f 27,936, G 37.5 GB) on 1 GPU and on 2 GPUs (host staging, and
#      --gpu-mpi), plus a small uneven case on 2 GPUs: the device product
#      against the host product within 1e-13 relative, same bits on a
#      repeat, and the time of one product on each side;
#   2. the CPU reference: the 50-step run on 1 rank without --device, with
#      --qd-green reading GREEN_CACHE (MUMPS for the startup check), judged
#      against G0's p1_r1 golden at TOL = floor.json key G3_to_G5_gpu and
#      reported, not gated (it is QD speed item 3's CPU path);
#   3. the 50-step gate with --device cuda --solver cudss --qd-operators
#      matrix --qd-green on 1 GPU, on 2 GPUs (host staging) and on 2 GPUs
#      with --qd-green-gpu-mpi, against G0's p1 goldens at the same TOL by
#      check_gate.sh, with --stage-timers 10 (stage_timers_final.txt); each
#      run is also compared with the CPU reference of step 2 by
#      bp5_gate_check.py (reported, not gated).  A miss against G0 is
#      reported with both numbers; the tolerance is not changed here;
#   4. the transfer trace: the 1-GPU Green run under Nsight Systems, 50
#      steps with checkpoints at steps 0, 25 and 50, gated by
#      g4_trace_check.py as in G4 (no copy of one value per fault DOF or
#      more inside a stage, between stages or in the monitor outside
#      output).
#   Nsight Systems is NSYS if set, else nsys on PATH, else the CUDA
#   toolkit's; without it the trace step fails.
#
# GitHub is the version control: the clone must be clean and at the commit
# its upstream branch points to, or the job stops.

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g1_cuda_env.sh"
g1_check_checkout

GREEN_CACHE="${GREEN_CACHE:-${SCRATCH:-${WORK}}/seas-qd-green/bp5_prod_p1.bin}"
if [ ! -s "${GREEN_CACHE}" ]; then
   echo "ERROR: no Green's-function cache file ${GREEN_CACHE}: run the CPU" >&2
   echo "       build first (jobs/vista_qd_speed, step green, g_n1_r144)." >&2
   exit 1
fi
echo "Green cache: ${GREEN_CACHE} ($(stat -c '%s bytes, modified %y' "${GREEN_CACHE}"))"

MODE=cuda GPU_ARCH=90 STAGE="seas" JOBS=72 MODULES=none bash "${SEAS_MFEM_ROOT}/build_horizon.sh" || exit 1
g0_check_binary
( cd "${SEAS_DIR}" && rm -f tests/parallel/test_qd_green_device.o &&
  make -j8 seas_test_qd_green_device ) || exit 1

cd "${SEAS_DIR}"
OUT_ROOT="bp5/goldens/g3c_gpu_gate/job${SLURM_JOB_ID}"
mkdir -p "${OUT_ROOT}"
g0_provenance "${OUT_ROOT}"
echo "green cache: ${GREEN_CACHE}" >> "${OUT_ROOT}/provenance.txt"
overall=0

GOLD="${G1_GOLDEN_ROOT}/p1_r1"
GOLD2="${G1_GOLDEN_ROOT}/p1_r2"
for g in "${GOLD}" "${GOLD2}"; do
   if [ ! -f "${g}/floor.json" ] || [ ! -d "${g}/run_a" ]; then
      echo "ERROR: no G0 golden in ${g}" >&2
      exit 1
   fi
done
TOL="$(python3 -c "import json; print(json.load(open('${GOLD}/floor.json'))['tol']['G3_to_G5_gpu'])")"
TOL2="$(python3 -c "import json; print(json.load(open('${GOLD2}/floor.json'))['tol']['G3_to_G5_gpu'])")"
GREEN_FLAGS=(--qd-operators matrix --qd-green --qd-green-cache "${GREEN_CACHE}")

# ---- 1. Device product against the host product ----
for CASE in "1:--nf 27936" "2:--nf 27936" "2:--nf 27936 --gpu-mpi" "2:--nf 2000 --uneven"; do
   NP="${CASE%%:*}"
   ARGS="${CASE#*:}"
   TAG="unit_r${NP}$(echo "${ARGS}" | tr -d ' -' | sed 's/nf/_nf/')"
   LOG="${OUT_ROOT}/${TAG}.log"
   echo "=== unit test ${ARGS} on ${NP} GPU(s): start $(date -u +%T) ==="
   # shellcheck disable=SC2086
   ibrun -n "${NP}" -o 0 ./seas_test_qd_green_device --device cuda ${ARGS} > "${LOG}" 2>&1
   rc=$?
   grep -E "QD Green device test|Device product|PASS|FAIL|Times|Apply" "${LOG}"
   if [ "${rc}" = 0 ] && grep -q "^ALL PASS" "${LOG}"; then
      echo "UNIT ${TAG}: PASS"
   else
      echo "UNIT ${TAG}: FAIL (exit ${rc})"
      overall=1
   fi
done

# run_case TAG NP DEVICE_ARGS... : one 50-step run into OUT_ROOT/TAG
run_case() {
   local TAG="$1" NP="$2"
   shift 2
   local OUT="${OUT_ROOT}/${TAG}" PFX="g3c_${TAG}" G="${GOLD}" T="${TOL}"
   if [ "${NP}" = 2 ]; then G="${GOLD2}"; T="${TOL2}"; fi
   mkdir -p "${OUT}"
   echo "=== ${TAG} on ${NP} rank(s): TOL ${T}, start $(date -u +%T) ==="
   ibrun -n "${NP}" -o 0 ./seas_bp5_full --mesh "${G0_MESH_PROD}" "${G0_V63_FLAGS[@]}" \
         --order 1 "${GREEN_FLAGS[@]}" "$@" --stage-timers 10 \
         --max-steps 50 --write-every-step --checkpoint-interval 0 \
         --output-dir "${OUT}" --output-prefix "${PFX}" \
         --ref-dir "${G}/run_a" --ref-prefix "g0_p1_r${NP}" \
         --regression-tolerance "${T}" > "${OUT}/run.log" 2>&1
   RC=$?
   echo "${RC}" > "${OUT}/exit_code"
   echo "=== ${TAG}: exit ${RC} at $(date -u +%T) ==="
   grep '^\[QD-GREEN\]' "${OUT}/run.log" | tee "${OUT}/qd_green.txt"
   sed -n '/^\[STAGE-TIMERS\] end of run/,$p' "${OUT}/run.log" \
      | grep '^\[STAGE-TIMERS\]' > "${OUT}/stage_timers_final.txt"
   cat "${OUT}/stage_timers_final.txt"
   bash jobs/horizon/check_gate.sh "${OUT}" "${PFX}" "${G}/run_a" "g0_p1_r${NP}" \
        "${T}" "${OUT}/run.log" "${RC}"
}

# ---- 2. CPU reference (QD speed item 3's path; reported) ----
CPU_TAG="cpu_green_r1"
if run_case "${CPU_TAG}" 1 --solver mumps; then
   echo "CPU_GREEN vs G0: PASS (report)"
else
   echo "CPU_GREEN vs G0: FAIL (report; see ${OUT_ROOT}/${CPU_TAG})"
fi
CPU_OUT="${OUT_ROOT}/${CPU_TAG}"

# ---- 3. 50-step gates on the GPU ----
for CASE in "gpu_green_r1:1:" "gpu_green_r2:2:" "gpu_green_r2_gpumpi:2:--qd-green-gpu-mpi"; do
   IFS=: read -r TAG NP EXTRA <<< "${CASE}"
   # shellcheck disable=SC2086
   if run_case "${TAG}" "${NP}" --device cuda --solver cudss ${EXTRA}; then
      echo "GATE ${TAG}: PASS"
   else
      echo "GATE ${TAG}: FAIL"
      overall=1
   fi
   if [ -s "${CPU_OUT}/exit_code" ] && [ "$(cat "${CPU_OUT}/exit_code")" = 0 ]; then
      python3 bp5/bp5_gate_check.py --run-dir "${OUT_ROOT}/${TAG}" \
              --run-prefix "g3c_${TAG}" --golden-dir "${CPU_OUT}" \
              --golden-prefix "g3c_${CPU_TAG}" --tol "${TOL}" \
              --json "${OUT_ROOT}/${TAG}/vs_cpu_green.json" \
         | sed 's/^/  vs CPU green: /'
      echo "VS_CPU_GREEN ${TAG}: exit ${PIPESTATUS[0]} (report)"
   fi
done

# ---- 4. Transfer trace ----
NSYS_BIN="${NSYS:-$(command -v nsys 2>/dev/null || true)}"
if [ -z "${NSYS_BIN}" ]; then
   for root in "${CUDA_HOME:-}" "${TACC_CUDA_DIR:-}" "${CUDA_PATH:-}"; do
      [ -n "${root}" ] || continue
      if [ -x "${root}/bin/nsys" ]; then
         NSYS_BIN="${root}/bin/nsys"
         break
      fi
      cand="$(ls -d "${root}"/nsight-systems-*/bin/nsys 2>/dev/null | sort -V | tail -n 1)"
      if [ -n "${cand}" ]; then
         NSYS_BIN="${cand}"
         break
      fi
   done
fi
if [ -z "${NSYS_BIN}" ] || [ ! -x "${NSYS_BIN}" ]; then
   echo "TRACE: FAIL (no Nsight Systems: set NSYS=/path/to/nsys or load a module that has it)"
   overall=1
else
   echo "Nsight Systems: ${NSYS_BIN} ($("${NSYS_BIN}" --version 2>&1 | head -n 1))"
   TAG="trace_green_r1"
   OUT="${OUT_ROOT}/${TAG}"
   PFX="g3c_${TAG}"
   mkdir -p "${OUT}"
   echo "=== ${TAG} on 1 GPU under Nsight Systems: start $(date -u +%T) ==="
   ibrun -n 1 -o 0 "${NSYS_BIN}" profile --trace=cuda,nvtx --sample=none \
         --cpuctxsw=none --force-overwrite=true -o "${OUT}/trace_rank0" \
         ./seas_bp5_full --mesh "${G0_MESH_PROD}" "${G0_V63_FLAGS[@]}" \
         --order 1 "${GREEN_FLAGS[@]}" --device cuda --solver cudss \
         --stage-timers 10 --max-steps 50 --checkpoint-interval 25 \
         --output-dir "${OUT}" --output-prefix "${PFX}" > "${OUT}/run.log" 2>&1
   rc=$?
   echo "${rc}" > "${OUT}/exit_code"
   echo "=== ${TAG}: exit ${rc} at $(date -u +%T) ==="
   FDOFS="$(sed -n 's/^ *Local fault DOFs: \([0-9]*\) (rank 0).*/\1/p' "${OUT}/run.log" | head -n 1)"
   if [ "${rc}" != 0 ] || [ -z "${FDOFS}" ]; then
      echo "TRACE ${TAG}: FAIL (run exit ${rc}, fault DOFs '${FDOFS}')"
      overall=1
   else
      REP="${OUT}/trace_rank0"
      if ! "${NSYS_BIN}" export --type=sqlite --force-overwrite=true \
            -o "${REP}.sqlite" "${REP}.nsys-rep" > "${REP}_export.log" 2>&1; then
         echo "TRACE ${TAG}: FAIL (nsys export, see ${REP}_export.log)"
         overall=1
      else
         python3 jobs/vista/g4_trace_check.py "${REP}.sqlite" \
                 --fault-dofs "${FDOFS}" --gate > "${REP}_check.txt" 2>&1
         crc=$?
         cat "${REP}_check.txt"
         echo "TRACE ${TAG}: exit ${crc}"
         [ "${crc}" != 0 ] && overall=1
      fi
   fi
fi

echo "G3c GPU gate overall: $([ "${overall}" = 0 ] && echo PASS || echo FAIL)"
exit "${overall}"
