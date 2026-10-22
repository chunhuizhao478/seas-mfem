#!/bin/bash
#SBATCH -J g4_cpu_gate
#SBATCH -o g4_cpu_gate_%j.out
#SBATCH -e g4_cpu_gate_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 06:00:00
#SBATCH -A EAR20006

# GPU plan G4, CPU gate on Vista: friction and the state rate as forall
# kernels, the monitor and output reading the state where it lives.
#
#   sbatch jobs/vista/g4_cpu_gate.sbatch   (from miniapps/seas of a CPU
#   clone, e.g. $WORK/seas-mfem-gpu-g2 built by cpu_clone_build.sbatch)
#
# Rebuilds seas_bp5_full and seas_test_device_friction in this checkout
# (MFEM and its dependencies are the G0 build's, unchanged), then
#   1. seas_test_device_friction on 10^6 random and 10^6 BP5-sampled fault
#      nodes plus the edge branches (tau <= 0, eta = 0, sigma_n <= 0,
#      psi/a > 700), kernel against the host friction code, within 1e-12
#      relative in V: on the host (bit-identical expected) and with
#      --device debug (kernels on the host, host memory guarded);
#   2. the 50-step gate against G0's goldens (bp5/goldens/g0_vista) with
#      --write-every-step:
#        p1_r1_loops       CPU, G1 CPU tolerance; its output files are also
#                          compared byte for byte with the golden (expected
#                          identical: the kernels do the host code's
#                          arithmetic in the same order);
#        p1_r1_loops_debug --device debug and checkpoints at steps 0, 25 and
#                          50 and at the end (the plan's device-memory run:
#                          a host loop that reads device-only memory stops
#                          the run);
#        p1_r1_matrix_debug the same with --qd-operators matrix, G2 tolerance;
#        p1_r8_loops_debug 8 ranks, --device debug, checkpoints;
#        p2_r1_loops       CPU at p=2, G1 CPU tolerance.
#   TOL comes from each golden's floor.json, judged by check_gate.sh.
#
# GitHub is the version control: the checkout must be clean and at the
# commit its upstream branch points to, or the job stops.

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g0_env.sh"

cd "${SEAS_MFEM_ROOT}"
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
   echo "ERROR: checkout has uncommitted changes; build only pushed commits." >&2
   git status --short --untracked-files=no >&2
   exit 1
fi
if [ "$(git rev-parse HEAD)" != "$(git rev-parse '@{u}' 2>/dev/null)" ]; then
   echo "ERROR: HEAD $(git rev-parse --short HEAD) is not its upstream" \
        "($(git rev-parse --short '@{u}' 2>/dev/null || echo none)); pull first." >&2
   exit 1
fi
echo "Checkout: ${SEAS_MFEM_ROOT} at $(git rev-parse HEAD) ($(git rev-parse --abbrev-ref HEAD))"

MODE=cpu STAGE="seas" JOBS=8 MODULES=none bash build_horizon.sh || exit 1
g0_check_binary
( cd "${SEAS_DIR}" && rm -f tests/parallel/test_device_friction.o &&
  make -j8 seas_test_device_friction ) || exit 1

cd "${SEAS_DIR}"
OUT_ROOT="bp5/goldens/g4_cpu_gate/job${SLURM_JOB_ID}"
mkdir -p "${OUT_ROOT}"
g0_provenance "${OUT_ROOT}"

overall=0

# ---- 1. Device friction against the host code ----
mkdir -p "${OUT_ROOT}/device_friction"
for DEVICE in cpu debug; do
   LOG="${OUT_ROOT}/device_friction/${DEVICE}.log"
   echo "=== device friction, --device ${DEVICE}: start $(date -u +%T) ==="
   ibrun -n 1 -o 0 ./seas_test_device_friction --device "${DEVICE}" \
         --n 1000000 > "${LOG}" 2>&1
   rc=$?
   grep -E "^\[|max rel|bit-identical|status|time:|PASS|FAIL" "${LOG}"
   if [ "${rc}" = 0 ] && grep -q "^PASSED" "${LOG}"; then
      echo "DEVICE_FRICTION ${DEVICE}: PASS"
   else
      echo "DEVICE_FRICTION ${DEVICE}: FAIL (exit ${rc})"
      overall=1
   fi
done

# ---- 2. 50-step gates ----
for CASE in 1:1:loops:cpu 1:1:loops:debug 1:1:matrix:debug 1:8:loops:debug \
            2:1:loops:cpu; do
   IFS=: read -r ORDER NP OPS DEVICE <<< "${CASE}"
   GOLD="bp5/goldens/g0_vista/p${ORDER}_r${NP}"
   GPFX="g0_p${ORDER}_r${NP}"
   TAG="p${ORDER}_r${NP}_${OPS}"
   [ "${DEVICE}" = cpu ] || TAG="${TAG}_${DEVICE}"
   OUT="${OUT_ROOT}/${TAG}"
   PFX="g4_${TAG}"
   if [ ! -f "${GOLD}/floor.json" ] || [ ! -d "${GOLD}/run_a" ]; then
      echo "ERROR: no G0 golden in ${GOLD}" >&2
      overall=1
      continue
   fi
   KEY=G1_host_syncs_cpu
   [ "${OPS}" = matrix ] && KEY=G2_operators_cpu
   CKPT=0
   [ "${DEVICE}" = debug ] && CKPT=25
   TOL="$(python3 -c "import json; print(json.load(open('${GOLD}/floor.json'))['tol']['${KEY}'])")"
   mkdir -p "${OUT}"
   echo "=== ${TAG}: TOL ${TOL} (${KEY}), checkpoint interval ${CKPT}, start $(date -u +%T) ==="
   ibrun -n "${NP}" -o 0 ./seas_bp5_full --mesh "${G0_MESH_PROD}" \
         "${G0_V63_FLAGS[@]}" --order "${ORDER}" --qd-operators "${OPS}" \
         --device "${DEVICE}" \
         --max-steps 50 --write-every-step --checkpoint-interval "${CKPT}" \
         --output-dir "${OUT}" --output-prefix "${PFX}" \
         --ref-dir "${GOLD}/run_a" --ref-prefix "${GPFX}" \
         --regression-tolerance "${TOL}" > "${OUT}/run.log" 2>&1
   rc=$?
   echo "${rc}" > "${OUT}/exit_code"
   echo "=== ${TAG}: exit ${rc} at $(date -u +%T) ==="
   if [ "${CKPT}" != 0 ]; then
      n_ckpt="$(grep -c '^Checkpoint written:' "${OUT}/run.log")"
      echo "CHECKPOINT ${TAG}: ${n_ckpt} written (expected 4: steps 0, 25 and 50, and the end of the run)"
      [ "${n_ckpt}" = 4 ] || overall=1
   fi
   if bash jobs/horizon/check_gate.sh "${OUT}" "${PFX}" "${GOLD}/run_a" "${GPFX}" \
         "${TOL}" "${OUT}/run.log" "${rc}"; then
      echo "GATE ${TAG}: PASS"
   else
      echo "GATE ${TAG}: FAIL"
      overall=1
   fi
   if [ "${TAG}" = p1_r1_loops ]; then
      same=1
      for f in "${GOLD}/run_a/${GPFX}"_*; do
         g="${OUT}/${PFX}_${f##*/${GPFX}_}"
         if ! cmp -s "${f}" "${g}"; then
            echo "BITWISE ${TAG}: ${f##*/} differs"
            same=0
         fi
      done
      [ "${same}" = 1 ] && echo "BITWISE ${TAG}: all output files identical"
   fi
done
echo "G4 CPU gate overall: $([ "${overall}" = 0 ] && echo PASS || echo FAIL)"
exit "${overall}"
