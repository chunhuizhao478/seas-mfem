#!/bin/bash
#SBATCH -J tpv102_ovr0_p2_200m
#SBATCH -o tpv102_ovr0_p2_200m_%j.out
#SBATCH -e tpv102_ovr0_p2_200m_%j.err
#SBATCH -p normal
#SBATCH -N 12
#SBATCH -n 600
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# TPV102 fault-dealiasing A/B/C — p2 PRODUCTION, arm [A of 3]: BASELINE.
#   p2 / ADER-O3 / pure upwind / scalar, 200 m mesh, --fault-overint 0, 0-12 s.
#   Byte-exact with the unmodified driver (K=0 skips SetFaultOverint).  The p2
#   baseline for the order-convergence test vs the p1 baseline.
# Standalone: build first (spatial_dyn_build.sbatch), then `sbatch` this.
# NOTE: p2 dt is smaller than p1; 12 s on the 200 m mesh likely exceeds one 48 h
# segment — it checkpoints every 2000 steps; resubmit with --restart if it times out.
set -euo pipefail

# ---- HARDCODED A/B/C knobs (this file = p2 baseline, 200 m) --------------
FAULT_OVERINT=0
MESH="tpv102/mesh/tpv102_200m.msh"
TFINAL=12
CONFIG_TOML="tpv102/configs/tpv102_spatial_p2.toml"

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8
export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH:-}"

_d="${SEAS_MFEM_ROOT:-${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "${BASH_SOURCE[0]:-$0}")")}}"
while [[ "${_d}" != "/" && ! -d "${_d}/miniapps/seas" ]]; do _d="$(dirname "${_d}")"; done
cd "${_d}/miniapps/seas"

if [[ ! -x ./seas_spatial_dyn_driver ]]; then
    echo "ERROR: ./seas_spatial_dyn_driver not built.  Build it first:"
    echo "       sbatch miniapps/seas/jobs/spatial_dyn_build.sbatch clean-rebuild"
    exit 1
fi
[[ -f "${CONFIG_TOML}" ]] || { echo "ERROR: ${CONFIG_TOML} not found."; exit 1; }
grep -Eq '^[[:space:]]*order[[:space:]]*=[[:space:]]*2([[:space:]]|$)' "${CONFIG_TOML}" \
    || { echo "ERROR: ${CONFIG_TOML} is not [mesh].order=2 (this job is p2)."; exit 4; }
grep -Eq '^[[:space:]]*interior_flux[[:space:]]*=[[:space:]]*"scalar"' "${CONFIG_TOML}" \
    || { echo "ERROR: ${CONFIG_TOML} is not interior_flux=\"scalar\"."; exit 4; }

if [[ ! -f "${MESH}" ]]; then
    echo "[ab] mesh ${MESH} missing — generating from ${MESH%.msh}.geo ..."
    command -v gmsh >/dev/null 2>&1 || { echo "ERROR: gmsh not on PATH; run on a login node: gmsh -format msh22 -3 ${MESH%.msh}.geo -o ${MESH}"; exit 2; }
    gmsh -format msh22 -3 "${MESH%.msh}.geo" -o "${MESH}"
fi

OUT="tpv102/out_p2_aderO3_pu_overint${FAULT_OVERINT}_200m_${SLURM_JOB_ID:-local}"
echo "[ab] p2 BASELINE  FAULT_OVERINT=${FAULT_OVERINT}  mesh=${MESH}  tfinal=${TFINAL}  out=${OUT}"
echo "[ab] Git: $(git rev-parse --abbrev-ref HEAD) $(git rev-parse --short HEAD)"

RESTART_FLAGS=""
if [[ -n "${TPV102_RESTART:-}" ]]; then
    RESTART_FLAGS="--restart ${TPV102_RESTART}"
    echo "[ab] RESTART from ${TPV102_RESTART} (this segment writes ${OUT})."
fi

ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --ader-order 3 \
    --time-integrator ader \
    --mixed-flux none \
    --fault-overint ${FAULT_OVERINT} \
    --tfinal "${TFINAL}" \
    --paraview \
    --paraview-fault-vtu \
    --paraview-fault-dt 0.1 \
    --paraview-max-snapshots 300 \
    --print-derived \
    --checkpoint-every 2000 \
    --output-dir "${OUT}" \
    --deriv-cache --shared-ck-recursion \
    ${RESTART_FLAGS} \
    --verify-dispatch

echo "[ab] done (p2 baseline, 200 m).  sigma_n station traces: ${OUT}/tpv102_station_*.dat"
