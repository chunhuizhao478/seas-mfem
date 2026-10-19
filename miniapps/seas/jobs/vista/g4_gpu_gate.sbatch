#!/bin/bash
#SBATCH -J g4_gpu_gate
#SBATCH -o g4_gpu_gate_%j.out
#SBATCH -e g4_gpu_gate_%j.err
#SBATCH -p gh
#SBATCH -N 2
#SBATCH -n 2
#SBATCH -t 03:00:00
#SBATCH -A EAR20006

# GPU plan G4 on Vista Grace-Hopper (gh) nodes, one GH200 per node, in the
# CUDA clone ($WORK/seas-mfem-gpu-cuda, built by g1_cuda_build.sbatch):
#
#   sbatch jobs/vista/g4_gpu_gate.sbatch   (from miniapps/seas of the CUDA clone)
#
# Rebuilds seas_bp5_full and seas_test_device_friction (MFEM, PETSc and
# cuDSS are the G1 build's, unchanged), then
#   1. seas_test_device_friction --device cuda on 10^6 random and 10^6
#      BP5-sampled fault nodes plus the edge branches: the GPU kernels
#      against the host friction code within 1e-12 relative in V (not
#      bitwise: GPU exp, pow and log10 differ in the last bit);
#   2. the 50-step gate with --device cuda --solver cudss --qd-operators
#      matrix, friction and time stepping now on the GPU, at p=1 on 1 and 2
#      GPUs and at p=2 on 1 GPU against G0's goldens (G1_GOLDEN_ROOT) at
#      TOL = floor.json key G3_to_G5_gpu (tol_phase 1e-8), judged by
#      check_gate.sh, with --stage-timers 10 (stage_timers_final.txt);
#   3. the transfer trace: p=1 on 1 GPU under Nsight Systems, 50 steps with
#      the adaptive output schedule and checkpoints at steps 0, 25 and 50,
#      so whole accepted steps with and without output are recorded.
#      g4_trace_check.py places every host-device copy in the driver's NVTX
#      ranges and fails on a copy of one value per fault DOF or more inside
#      a stage, between stages or in the monitor outside output (the plan:
#      "no state-sized copies inside a stage").  If it fails, PETSc TS or a
#      host loop moves the state every stage: stop and report before
#      choosing another integrator.  The same trace on 2 GPUs is reported,
#      not gated (halo exchanges may stage through the host).
#   Nsight Systems is NSYS if set, else nsys on PATH, else the CUDA
#   toolkit's; without it the trace step fails.
#
# GitHub is the version control: the clone must be clean and at the commit
# its upstream branch points to, or the job stops.

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g1_cuda_env.sh"
g1_check_checkout

MODE=cuda GPU_ARCH=90 STAGE="seas" JOBS=72 MODULES=none bash "${SEAS_MFEM_ROOT}/build_horizon.sh" || exit 1
g0_check_binary
( cd "${SEAS_DIR}" && rm -f tests/parallel/test_device_friction.o &&
  make -j8 seas_test_device_friction ) || exit 1

cd "${SEAS_DIR}"
OUT_ROOT="bp5/goldens/g4_gpu_gate/job${SLURM_JOB_ID}"
mkdir -p "${OUT_ROOT}"
g0_provenance "${OUT_ROOT}"
overall=0

# ---- 1. Device friction on the GPU ----
LOG="${OUT_ROOT}/device_friction_cuda.log"
echo "=== device friction, --device cuda: start $(date -u +%T) ==="
ibrun -n 1 -o 0 ./seas_test_device_friction --device cuda --n 1000000 > "${LOG}" 2>&1
rc=$?
grep -E "^\[|max rel|bit-identical|status|time:|PASS|FAIL" "${LOG}"
if [ "${rc}" = 0 ] && grep -q "^PASSED" "${LOG}"; then
   echo "DEVICE_FRICTION cuda: PASS"
else
   echo "DEVICE_FRICTION cuda: FAIL (exit ${rc})"
   overall=1
fi

# ---- 2. 50-step gates ----
for CASE in 1:1 1:2 2:1; do
   IFS=: read -r ORDER NP <<< "${CASE}"
   TAG="p${ORDER}_r${NP}"
   GOLD="${G1_GOLDEN_ROOT}/${TAG}"
   GPFX="g0_${TAG}"
   OUT="${OUT_ROOT}/${TAG}"
   PFX="g4_${TAG}"
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
         --solver cudss --stage-timers 10 \
         --max-steps 50 --write-every-step --checkpoint-interval 0 \
         --output-dir "${OUT}" --output-prefix "${PFX}" \
         --ref-dir "${GOLD}/run_a" --ref-prefix "${GPFX}" \
         --regression-tolerance "${TOL}" > "${OUT}/run.log" 2>&1
   rc=$?
   echo "${rc}" > "${OUT}/exit_code"
   echo "=== ${TAG}: exit ${rc} at $(date -u +%T) ==="
   sed -n '/^\[STAGE-TIMERS\] end of run/,$p' "${OUT}/run.log" \
      | grep '^\[STAGE-TIMERS\]' > "${OUT}/stage_timers_final.txt"
   cat "${OUT}/stage_timers_final.txt"
   if bash jobs/horizon/check_gate.sh "${OUT}" "${PFX}" "${GOLD}/run_a" "${GPFX}" \
         "${TOL}" "${OUT}/run.log" "${rc}"; then
      echo "GATE ${TAG}: PASS"
   else
      echo "GATE ${TAG}: FAIL"
      overall=1
   fi
done

# ---- 3. Transfer trace ----
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
   for NP in 1 2; do
      TAG="trace_p1_r${NP}"
      OUT="${OUT_ROOT}/${TAG}"
      PFX="g4_${TAG}"
      mkdir -p "${OUT}"
      cat > "${OUT}/nsys_rank.sh" <<EOF
#!/bin/bash
# One Nsight Systems report per rank.
r=\${OMPI_COMM_WORLD_RANK:-\${PMI_RANK:-\${SLURM_PROCID:-0}}}
exec "${NSYS_BIN}" profile --trace=cuda,nvtx --sample=none --cpuctxsw=none \\
     --force-overwrite=true -o "${OUT}/trace_rank\${r}" "\$@"
EOF
      echo "=== ${TAG} on ${NP} GPU(s) under Nsight Systems: start $(date -u +%T) ==="
      ibrun -n "${NP}" -o 0 bash "${OUT}/nsys_rank.sh" ./seas_bp5_full \
            --mesh "${G0_MESH_PROD}" "${G0_V63_FLAGS[@]}" \
            --order 1 --device cuda --qd-operators matrix \
            --solver cudss --stage-timers 10 \
            --max-steps 50 --checkpoint-interval 25 \
            --output-dir "${OUT}" --output-prefix "${PFX}" > "${OUT}/run.log" 2>&1
      rc=$?
      echo "${rc}" > "${OUT}/exit_code"
      echo "=== ${TAG}: exit ${rc} at $(date -u +%T) ==="
      FDOFS="$(sed -n 's/^ *Local fault DOFs: \([0-9]*\) (rank 0).*/\1/p' "${OUT}/run.log" | head -n 1)"
      if [ "${rc}" != 0 ] || [ -z "${FDOFS}" ]; then
         echo "TRACE ${TAG}: FAIL (run exit ${rc}, fault DOFs '${FDOFS}')"
         [ "${NP}" = 1 ] && overall=1
         continue
      fi
      for (( r = 0; r < NP; r++ )); do
         REP="${OUT}/trace_rank${r}"
         if ! "${NSYS_BIN}" export --type=sqlite --force-overwrite=true \
               -o "${REP}.sqlite" "${REP}.nsys-rep" > "${REP}_export.log" 2>&1; then
            echo "TRACE ${TAG} rank ${r}: FAIL (nsys export, see ${REP}_export.log)"
            [ "${NP}" = 1 ] && overall=1
            continue
         fi
         # Rank 0's fault DOF count is the run log's; on 2 ranks it sets
         # the threshold for both (a report, not a gate).
         GATE_ARG=()
         [ "${NP}" = 1 ] && GATE_ARG=(--gate)
         python3 jobs/vista/g4_trace_check.py "${REP}.sqlite" \
                 --fault-dofs "${FDOFS}" ${GATE_ARG[@]+"${GATE_ARG[@]}"} \
                 > "${REP}_check.txt" 2>&1
         crc=$?
         cat "${REP}_check.txt"
         echo "TRACE ${TAG} rank ${r}: exit ${crc}"
         if [ "${NP}" = 1 ] && [ "${crc}" != 0 ]; then overall=1; fi
      done
   done
fi

echo "G4 GPU gate overall: $([ "${overall}" = 0 ] && echo PASS || echo FAIL)"
exit "${overall}"
