#!/bin/bash
#SBATCH -J g1_cpu_gate
#SBATCH -o g1_cpu_gate_%j.out
#SBATCH -e g1_cpu_gate_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 06:00:00
#SBATCH -A EAR20006

# GPU plan G1, CPU gate on Vista: the host syncs and the --device flag must
# not change the CPU answer.
#
#   sbatch jobs/vista/g1_cpu_gate.sbatch      (submit from miniapps/seas)
#
# Rebuilds seas_bp5_full in this checkout (MFEM and its dependencies are the
# G0 build's, unchanged), then runs the 50-step p=1 gate at 1, 2, 4 and 8
# ranks against G0's goldens (bp5/goldens/g0_vista/p1_rN/run_a) at each
# golden's G1 CPU tolerance, TOL = max(1e-12, 10 floor_N) from its
# floor.json, judged by check_gate.sh.  At 1 rank the output files are
# also compared byte for byte with the golden (expected identical; a
# difference is recorded, not a gate failure).  A fifth run repeats 8 ranks
# with --device debug, MFEM's CPU stand-in for a GPU: it aborts on any host
# access to data that is valid only on the "device", so it finds a missing
# host sync before a GPU run would silently use stale data.
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

cd "${SEAS_DIR}"
OUT_ROOT="bp5/goldens/g1_cpu_gate/job${SLURM_JOB_ID}"
mkdir -p "${OUT_ROOT}"
g0_provenance "${OUT_ROOT}"

overall=0
for CASE in 1:cpu 2:cpu 4:cpu 8:cpu 8:debug; do
   NP="${CASE%%:*}"
   DEVICE="${CASE##*:}"
   GOLD="bp5/goldens/g0_vista/p1_r${NP}"
   GPFX="g0_p1_r${NP}"
   TAG="p1_r${NP}"
   [ "${DEVICE}" = cpu ] || TAG="${TAG}_${DEVICE}"
   OUT="${OUT_ROOT}/${TAG}"
   PFX="g1_${TAG}"
   if [ ! -f "${GOLD}/floor.json" ] || [ ! -d "${GOLD}/run_a" ]; then
      echo "ERROR: no G0 golden in ${GOLD}" >&2
      overall=1
      continue
   fi
   TOL="$(python3 -c "import json; print(json.load(open('${GOLD}/floor.json'))['tol']['G1_host_syncs_cpu'])")"
   mkdir -p "${OUT}"
   echo "=== ${TAG}: TOL ${TOL}, start $(date -u +%T) ==="
   ibrun -n "${NP}" -o 0 ./seas_bp5_full --mesh "${G0_MESH_PROD}" \
         "${G0_V63_FLAGS[@]}" --qd-operators loops --order 1 --device "${DEVICE}" \
         --max-steps 50 --write-every-step --checkpoint-interval 0 \
         --output-dir "${OUT}" --output-prefix "${PFX}" \
         --ref-dir "${GOLD}/run_a" --ref-prefix "${GPFX}" \
         --regression-tolerance "${TOL}" > "${OUT}/run.log" 2>&1
   rc=$?
   echo "${rc}" > "${OUT}/exit_code"
   echo "=== ${TAG}: exit ${rc} at $(date -u +%T) ==="
   if bash jobs/horizon/check_gate.sh "${OUT}" "${PFX}" "${GOLD}/run_a" "${GPFX}" \
         "${TOL}" "${OUT}/run.log" "${rc}"; then
      echo "GATE ${TAG}: PASS"
   else
      echo "GATE ${TAG}: FAIL"
      overall=1
   fi
   if [ "${TAG}" = p1_r1 ]; then
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
echo "G1 CPU gate overall: $([ "${overall}" = 0 ] && echo PASS || echo FAIL)"
exit "${overall}"
