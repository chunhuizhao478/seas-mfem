#!/bin/bash
#SBATCH -J qd_default_gate
#SBATCH -o qd_default_gate_%j.out
#SBATCH -e qd_default_gate_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 05:00:00
#SBATCH -A EAR20006

# QD speed item 1 (plan rev 2): the stage operators are now the default
# (--qd-operators auto: matrix with --dg-method IP, loops with BR2 or
# SEAS_DEBUG_FACE_TRACE).  This gate checks that the default gives G0's
# answer and picks the right path:
#   1. the 50-step gate with no --qd-operators flag at p=1 on 1, 2, 4 and 8
#      ranks and at p=2 on 1 rank, against G0's goldens
#      (bp5/goldens/g0_vista/pP_rN/run_a) at each golden's G2 tolerance
#      (floor.json tol.G2_operators_cpu, as jobs/vista/g2_cpu_gate.sbatch
#      used for --qd-operators matrix), judged by jobs/horizon/check_gate.sh;
#      each log must say "qd-operators: matrix (default)";
#   2. two start-ups on 8 ranks that stop after the factorization
#      (--exit-after-factorize): --dg-method BR2 must say "qd-operators:
#      loops (default; matrix needs --dg-method IP)", and
#      SEAS_DEBUG_FACE_TRACE=1 must say "qd-operators: loops (default;
#      SEAS_DEBUG_FACE_TRACE needs the loops)".
#
#   sbatch jobs/vista_qd_speed/qd_default_gate.sbatch   (from miniapps/seas;
#   jobs/vista_qd_speed/qd_speed_submit.sh does this)
#
# The goldens are G0's, copied into this clone by
# jobs/vista/cpu_clone_build.sbatch.  GitHub is the version control: the
# checkout must be clean and at the commit its upstream branch points to, or
# the job stops.

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g0_env.sh"
g0_check_binary

cd "${SEAS_MFEM_ROOT}"
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
   echo "ERROR: checkout has uncommitted changes; run only pushed commits." >&2
   git status --short --untracked-files=no >&2
   exit 1
fi
if [ "$(git rev-parse HEAD)" != "$(git rev-parse '@{u}' 2>/dev/null)" ]; then
   echo "ERROR: HEAD $(git rev-parse --short HEAD) is not its upstream" \
        "($(git rev-parse --short '@{u}' 2>/dev/null || echo none)); pull first." >&2
   exit 1
fi
echo "Checkout: ${SEAS_MFEM_ROOT} at $(git rev-parse HEAD) ($(git rev-parse --abbrev-ref HEAD))"

cd "${SEAS_DIR}"
OUT_ROOT="bp5/results_qd_speed/default_gate_job${SLURM_JOB_ID}"
mkdir -p "${OUT_ROOT}"
g0_provenance "${OUT_ROOT}"

overall=0

# ---- 1. 50-step gates with the default ----
for CASE in 1:1 1:2 1:4 1:8 2:1; do
   IFS=: read -r ORDER NP <<< "${CASE}"
   GOLD="bp5/goldens/g0_vista/p${ORDER}_r${NP}"
   GPFX="g0_p${ORDER}_r${NP}"
   TAG="p${ORDER}_r${NP}_default"
   OUT="${OUT_ROOT}/${TAG}"
   PFX="qd_${TAG}"
   if [ ! -f "${GOLD}/floor.json" ] || [ ! -d "${GOLD}/run_a" ]; then
      echo "ERROR: no G0 golden in ${GOLD}" >&2
      overall=1
      continue
   fi
   TOL="$(python3 -c "import json; print(json.load(open('${GOLD}/floor.json'))['tol']['G2_operators_cpu'])")"
   mkdir -p "${OUT}"
   echo "=== ${TAG}: TOL ${TOL} (G2_operators_cpu), start $(date -u +%T) ==="
   ibrun -n "${NP}" -o 0 ./seas_bp5_full --mesh "${G0_MESH_PROD}" \
         "${G0_V63_FLAGS[@]}" --order "${ORDER}" \
         --max-steps 50 --write-every-step --checkpoint-interval 0 \
         --output-dir "${OUT}" --output-prefix "${PFX}" \
         --ref-dir "${GOLD}/run_a" --ref-prefix "${GPFX}" \
         --regression-tolerance "${TOL}" > "${OUT}/run.log" 2>&1
   rc=$?
   echo "${rc}" > "${OUT}/exit_code"
   echo "=== ${TAG}: exit ${rc} at $(date -u +%T) ==="
   if ! grep -q "qd-operators: matrix (default)" "${OUT}/run.log"; then
      echo "PATH ${TAG}: FAIL (log does not say 'qd-operators: matrix (default)')"
      overall=1
   fi
   { echo "--- ${TAG}"; grep "\[QD-OPERATORS\]" "${OUT}/run.log"; } \
      >> "${OUT_ROOT}/qd_operators.txt"
   if bash jobs/horizon/check_gate.sh "${OUT}" "${PFX}" "${GOLD}/run_a" "${GPFX}" \
         "${TOL}" "${OUT}/run.log" "${rc}"; then
      echo "GATE ${TAG}: PASS"
   else
      echo "GATE ${TAG}: FAIL"
      overall=1
   fi
done

# ---- 2. Which path the default picks ----
for CASE in br2 face_trace; do
   OUT="${OUT_ROOT}/path_${CASE}"
   mkdir -p "${OUT}"
   if [ "${CASE}" = br2 ]; then
      EXTRA=(--dg-method BR2)
      WANT="qd-operators: loops (default; matrix needs --dg-method IP)"
      ENVSET=()
   else
      EXTRA=()
      WANT="qd-operators: loops (default; SEAS_DEBUG_FACE_TRACE needs the loops)"
      ENVSET=(SEAS_DEBUG_FACE_TRACE=1)
   fi
   echo "=== path_${CASE}: start $(date -u +%T) ==="
   env "${ENVSET[@]}" ibrun -n 8 -o 0 ./seas_bp5_full --mesh "${G0_MESH_PROD}" \
         "${G0_V63_FLAGS[@]}" "${EXTRA[@]}" --order 1 --exit-after-factorize \
         --checkpoint-interval 0 --output-dir "${OUT}" \
         --output-prefix "qd_path_${CASE}" > "${OUT}/run.log" 2>&1
   rc=$?
   echo "${rc}" > "${OUT}/exit_code"
   if [ "${rc}" = 0 ] && grep -qF "${WANT}" "${OUT}/run.log"; then
      echo "PATH ${CASE}: PASS"
   else
      echo "PATH ${CASE}: FAIL (exit ${rc}; want '${WANT}', log says:" \
           "$(grep 'qd-operators:' "${OUT}/run.log" | head -n 1))"
      overall=1
   fi
done

cat "${OUT_ROOT}/qd_operators.txt" 2>/dev/null
echo "QD default gate overall: $([ "${overall}" = 0 ] && echo PASS || echo FAIL)"
exit "${overall}"
