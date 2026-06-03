#!/bin/bash
# Submit the TPV102 fault-flux over-integration A/B on Frontera:
#   (1) build seas_spatial_dyn_driver from the current branch (clean-rebuild),
#   (2) FAULT_OVERINT=0 baseline run  (afterok build),
#   (3) FAULT_OVERINT=2 over-int run  (afterok build).
# The two runs differ in EXACTLY the --fault-overint factor.
#
# PREREQ: be on the feat branch with the Phase 1+4 code:
#     git fetch origin && git checkout feat/fault-overint-resample
#
# USAGE (from anywhere in the repo):
#     bash miniapps/seas/jobs/tpv102_spatial/submit_fault_overint_ab.sh          # production: 200 m, normal, tfinal=12
#     CHEAP=1 bash .../submit_fault_overint_ab.sh                                 # cheap 1st pass: 1000 m, flex, tfinal=5
#     NO_BUILD=1 bash .../submit_fault_overint_ab.sh                             # skip the build step (binary already built)
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Walk up to the repo root (dir containing miniapps/seas).
ROOT="${HERE}"
while [[ "${ROOT}" != "/" && ! -d "${ROOT}/miniapps/seas" ]]; do ROOT="$(dirname "${ROOT}")"; done
BUILD_SBATCH="${ROOT}/miniapps/seas/jobs/spatial_dyn_build.sbatch"
RUN_SBATCH="${HERE}/tpv102_p1_aderO2_pureupwind_fault_overint_ab.sbatch"
[[ -f "${BUILD_SBATCH}" ]] || { echo "ERROR: build sbatch not found: ${BUILD_SBATCH}"; exit 1; }
[[ -f "${RUN_SBATCH}"   ]] || { echo "ERROR: run sbatch not found: ${RUN_SBATCH}"; exit 1; }

# CHEAP=1 => coarse 1000 m mesh + flex + small scale + short tfinal (fast first
# answer; the speckle is a per-element artifact, present at 1000 m too).
if [[ "${CHEAP:-0}" == "1" ]]; then
    SBATCH_SCALE=(-p flex -N 2 -n 100 -t 02:00:00)
    RUN_ENV="TPV102_MESH=tpv102/mesh/tpv102_1000m.msh,TPV102_TFINAL=${TPV102_TFINAL:-5}"
    echo "[ab] CHEAP first pass: 1000 m mesh, flex, 2 nodes / 100 ranks, tfinal=${TPV102_TFINAL:-5}."
else
    SBATCH_SCALE=()   # use the run sbatch's own #SBATCH headers (200 m, normal, 10N/500r, 48h)
    RUN_ENV="TPV102_TFINAL=${TPV102_TFINAL:-12}"
    echo "[ab] PRODUCTION A/B: 200 m mesh, normal, 10 nodes / 500 ranks, tfinal=${TPV102_TFINAL:-12}."
fi

# Frontera's `sbatch` is wrapped with a verification banner printed to stdout,
# which pollutes `$(sbatch --parsable ...)`.  Keep ONLY the trailing bare-number
# job id (robust whether the banner lands on stdout, stderr, or the tty).
submit_id() { sbatch --parsable "$@" 2>&1 | grep -xE '[0-9]+' | tail -n1 || true; }

# --- (1) build -------------------------------------------------------------
DEP=()
if [[ "${NO_BUILD:-0}" != "1" ]]; then
    echo "[ab] submitting build (clean-rebuild) ..."
    BID="$(submit_id "${BUILD_SBATCH}" clean-rebuild)"
    [[ -n "${BID}" ]] || { echo "[ab] ERROR: build job submission failed (no job id parsed)."; exit 1; }
    echo "[ab]   build job id = ${BID}"
    DEP=(--dependency="afterok:${BID}")
else
    echo "[ab] NO_BUILD=1 — skipping build (assuming ./seas_spatial_dyn_driver already built)."
fi

# --- (2) baseline (FAULT_OVERINT=0) ---------------------------------------
OFF_ID="$(submit_id ${DEP[@]+"${DEP[@]}"} ${SBATCH_SCALE[@]+"${SBATCH_SCALE[@]}"} \
    --export="ALL,FAULT_OVERINT=0,${RUN_ENV}" "${RUN_SBATCH}")"
[[ -n "${OFF_ID}" ]] || { echo "[ab] ERROR: baseline run submission failed."; exit 1; }
echo "[ab] baseline  (FAULT_OVERINT=0) job id = ${OFF_ID}"

# --- (3) over-int (FAULT_OVERINT=2) ---------------------------------------
ON_ID="$(submit_id ${DEP[@]+"${DEP[@]}"} ${SBATCH_SCALE[@]+"${SBATCH_SCALE[@]}"} \
    --export="ALL,FAULT_OVERINT=2,${RUN_ENV}" "${RUN_SBATCH}")"
[[ -n "${ON_ID}" ]] || { echo "[ab] ERROR: over-int run submission failed."; exit 1; }
echo "[ab] over-int  (FAULT_OVERINT=2) job id = ${ON_ID}"

echo ""
echo "[ab] submitted.  When both finish, compare the on-fault sigma_n station"
echo "     traces (the out dirs auto-encode the factor):"
echo "       tpv102/out_p1_aderO2_pu_overint0_${OFF_ID}/tpv102_station_*.dat   (baseline: drifts)"
echo "       tpv102/out_p1_aderO2_pu_overint2_${ON_ID}/tpv102_station_*.dat    (over-int: should stay ~120 MPa)"
echo "     Confirm the ON run applied it:  grep -h '\\[fault\\] over-integration ON' tpv102_overint_ab_${ON_ID}.out"
