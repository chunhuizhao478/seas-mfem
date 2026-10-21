#!/bin/bash
#SBATCH -J g0_golden
#SBATCH -o g0_golden_%j.out
#SBATCH -e g0_golden_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -A EAR20006

# GPU plan G0: one 50-step golden and its floor at one order and rank count.
#
#   sbatch -n N --export=ALL,ORDER=P jobs/vista/g0_goldens.sbatch
#   (submit from miniapps/seas; jobs/vista/g0_submit.sh does this)
#
# Three runs of the same binary, one after the other, on the production mesh
# with the v63 flags, --max-steps 50 --write-every-step:
#   run_a  the golden;
#   run_b  a repeat, compared by the binary's comparator against run_a;
#   run_p  a repeat with --perturb-initial-state 1 (state scaled by
#          1 + 1e-15 r), to see how far a round-off change moves 50 steps.
# floor_N is the larger of bp5_gate_check.py's max error for run_b and for
# run_p against run_a.  g0_floor_summary.py writes the gate tolerances
# TOL = max(tol_phase, 10 floor_N) and the job ends with check_gate.sh on
# run_b at the G1 CPU tolerance.
#
# p=2 runs on more than one rank wait for the shared fault-face fix
# (domain/elasticity_operator_setup.inl:1139-1143 reorders only the corner
# nodes of a shared fault face); set ALLOW_P2_MULTIRANK=1 once it is in.

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g0_env.sh"
g0_check_binary

ORDER="${ORDER:?set ORDER=1 or 2}"
NRANKS="${SLURM_NTASKS:-1}"
if [ "${ORDER}" -ge 2 ] && [ "${NRANKS}" -gt 1 ] && [ "${ALLOW_P2_MULTIRANK:-0}" != 1 ]; then
   echo "ERROR: p=${ORDER} on ${NRANKS} ranks waits for the shared fault-face fix." >&2
   exit 1
fi

cd "${SEAS_DIR}"
TAG="p${ORDER}_r${NRANKS}"
OUT="bp5/goldens/g0_vista/${TAG}"
PFX="g0_${TAG}"
if [ -e "${OUT}/run_a/${PFX}_global.txt" ]; then
   echo "ERROR: ${OUT} already holds a golden; move it away first." >&2
   exit 1
fi
mkdir -p "${OUT}"/run_a "${OUT}"/run_b "${OUT}"/run_p
g0_provenance "${OUT}"

COMMON=(--mesh "${G0_MESH_PROD}" "${G0_V63_FLAGS[@]}" --order "${ORDER}"
        --max-steps 50 --write-every-step --checkpoint-interval 0)

run() {
   local name="$1"; shift
   echo "=== ${TAG} ${name}: $(date -u +%T) ==="
   ibrun ./seas_bp5_full "${COMMON[@]}" --output-dir "${OUT}/${name}" \
         --output-prefix "${PFX}" "$@" > "${OUT}/${name}/run.log" 2>&1
   local rc=$?
   echo "${rc}" > "${OUT}/${name}/exit_code"
   echo "=== ${TAG} ${name}: exit ${rc} at $(date -u +%T) ==="
   tail -n 25 "${OUT}/${name}/run.log"
   return 0
}

run run_a
run run_b --ref-dir "${OUT}/run_a" --ref-prefix "${PFX}" --regression-tolerance 1e-6
run run_p --perturb-initial-state 1

for r in run_b run_p; do
   python3 bp5/bp5_gate_check.py --run-dir "${OUT}/${r}" --run-prefix "${PFX}" \
      --golden-dir "${OUT}/run_a" --golden-prefix "${PFX}" \
      --json "${OUT}/${r}_vs_a.json" || true
done
python3 jobs/vista/g0_floor_summary.py --order "${ORDER}" --ranks "${NRANKS}" \
   --repeat "${OUT}/run_b_vs_a.json" --perturbed "${OUT}/run_p_vs_a.json" \
   --out "${OUT}/floor.json"
TOL_G1_CPU="$(python3 -c "import json; print(json.load(open('${OUT}/floor.json'))['tol']['G1_host_syncs_cpu'])")"

bash jobs/horizon/check_gate.sh "${OUT}/run_b" "${PFX}" "${OUT}/run_a" "${PFX}" \
   "${TOL_G1_CPU}" "${OUT}/run_b/run.log" "$(cat "${OUT}/run_b/exit_code")"
