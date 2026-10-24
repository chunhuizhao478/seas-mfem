#!/bin/bash
#SBATCH -J g2_cpu_gate
#SBATCH -o g2_cpu_gate_%j.out
#SBATCH -e g2_cpu_gate_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 07:00:00
#SBATCH -A EAR20006

# GPU plan G2, CPU gate on Vista: the per-stage face loops as sparse
# operators (--qd-operators matrix) must give G0's answer.
#
#   sbatch jobs/vista/g2_cpu_gate.sbatch      (submit from miniapps/seas)
#
# Rebuilds seas_bp5_full and seas_test_qd_stage_operators in this checkout
# (MFEM and its dependencies are the G0 build's, unchanged), then
#   1. the operator tests (acceptance 1-4) on the production mesh, p=1 and
#      p=2, on 1, 2, 4 and 8 ranks with METIS, and on 2, 4 and 8 ranks with
#      the partition boundary on the fault (shared fault faces, owners on
#      both sides);
#   2. the 50-step gate with --qd-operators matrix at p=1 on 1, 2, 4 and 8
#      ranks and at p=2 on 1 rank against G0's goldens
#      (bp5/goldens/g0_vista/pP_rN/run_a) at each golden's G2 tolerance,
#      TOL = max(1e-10, 10 floor_N) from its floor.json, judged by
#      check_gate.sh;
#   3. the same 50-step gate with the default face loops at the G1 CPU
#      tolerance, since G2 refactors the loops (p=2 on 1 rank too: G2a
#      reorders the p=2 face nodes, and this run shows that change without
#      the operators); at 1 rank the p=1 output files are also compared
#      byte for byte with the golden (expected identical; a difference is
#      recorded, not a gate failure).  With the 8-rank --device debug run
#      this repeats g1_cpu_gate.sbatch, so one job covers the CPU gates of
#      G1 and G2.
# The matrix runs print the operator sizes, memory and build time, and the
# host loading g(t): evaluations per accepted and rejected step and time
# per evaluation ([QD-OPERATORS] lines, collected in qd_operators.txt).
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
( cd "${SEAS_DIR}" && rm -f tests/parallel/test_qd_stage_operators.o &&
  make -j8 seas_test_qd_stage_operators ) || exit 1

cd "${SEAS_DIR}"
OUT_ROOT="bp5/goldens/g2_cpu_gate/job${SLURM_JOB_ID}"
mkdir -p "${OUT_ROOT}"
g0_provenance "${OUT_ROOT}"

overall=0

# ---- 1. Operator tests ----
mkdir -p "${OUT_ROOT}/operator_tests"
for ORDER in 1 2; do
   for CASE in 1:metis 2:metis 4:metis 8:metis 2:fault 4:fault 8:fault; do
      NP="${CASE%%:*}"
      PART="${CASE##*:}"
      LOG="${OUT_ROOT}/operator_tests/p${ORDER}_r${NP}_${PART}.log"
      echo "=== operator tests p=${ORDER} ${NP} ranks ${PART}: start $(date -u +%T) ==="
      ibrun -n "${NP}" -o 0 ./seas_test_qd_stage_operators \
            --mesh "${G0_MESH_PROD}" --order "${ORDER}" --partition "${PART}" \
            > "${LOG}" 2>&1
      rc=$?
      grep -E "PASS|FAIL|shared fault faces:" "${LOG}"
      if [ "${rc}" = 0 ] && grep -q "ALL PASSED" "${LOG}"; then
         echo "OPERATORS p${ORDER}_r${NP}_${PART}: PASS"
      else
         echo "OPERATORS p${ORDER}_r${NP}_${PART}: FAIL (exit ${rc})"
         overall=1
      fi
   done
done

# ---- 2, 3. 50-step gates ----
for CASE in 1:1:matrix:cpu 1:2:matrix:cpu 1:4:matrix:cpu 1:8:matrix:cpu \
            2:1:matrix:cpu 1:1:loops:cpu 1:2:loops:cpu 1:4:loops:cpu \
            1:8:loops:cpu 2:1:loops:cpu 1:8:loops:debug; do
   IFS=: read -r ORDER NP OPS DEVICE <<< "${CASE}"
   GOLD="bp5/goldens/g0_vista/p${ORDER}_r${NP}"
   GPFX="g0_p${ORDER}_r${NP}"
   TAG="p${ORDER}_r${NP}_${OPS}"
   [ "${DEVICE}" = cpu ] || TAG="${TAG}_${DEVICE}"
   OUT="${OUT_ROOT}/${TAG}"
   PFX="g2_${TAG}"
   if [ ! -f "${GOLD}/floor.json" ] || [ ! -d "${GOLD}/run_a" ]; then
      echo "ERROR: no G0 golden in ${GOLD}" >&2
      overall=1
      continue
   fi
   KEY=G2_operators_cpu
   [ "${OPS}" = loops ] && KEY=G1_host_syncs_cpu
   TOL="$(python3 -c "import json; print(json.load(open('${GOLD}/floor.json'))['tol']['${KEY}'])")"
   mkdir -p "${OUT}"
   echo "=== ${TAG}: TOL ${TOL} (${KEY}), start $(date -u +%T) ==="
   ibrun -n "${NP}" -o 0 ./seas_bp5_full --mesh "${G0_MESH_PROD}" \
         "${G0_V63_FLAGS[@]}" --order "${ORDER}" --qd-operators "${OPS}" \
         --device "${DEVICE}" \
         --max-steps 50 --write-every-step --checkpoint-interval 0 \
         --output-dir "${OUT}" --output-prefix "${PFX}" \
         --ref-dir "${GOLD}/run_a" --ref-prefix "${GPFX}" \
         --regression-tolerance "${TOL}" > "${OUT}/run.log" 2>&1
   rc=$?
   echo "${rc}" > "${OUT}/exit_code"
   echo "=== ${TAG}: exit ${rc} at $(date -u +%T) ==="
   if [ "${OPS}" = matrix ]; then
      { echo "--- ${TAG}"; grep "\[QD-OPERATORS\]" "${OUT}/run.log"; } \
         >> "${OUT_ROOT}/qd_operators.txt"
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
cat "${OUT_ROOT}/qd_operators.txt" 2>/dev/null
echo "G2 CPU gate overall: $([ "${overall}" = 0 ] && echo PASS || echo FAIL)"
exit "${overall}"
