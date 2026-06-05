#!/bin/bash
#SBATCH -J tpv102_ovr0_1km
#SBATCH -o tpv102_ovr0_1km_%j.out
#SBATCH -e tpv102_ovr0_1km_%j.err
#SBATCH -p flex
#SBATCH -N 2
#SBATCH -n 100
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV102 fault-flux over-integration A/B — CHEAP pass, run [1 of 4]: BASELINE.
#   p1 / ADER-O2 / pure upwind / scalar, 1000 m mesh, --fault-overint 0.
# Standalone: build first (spatial_dyn_build.sbatch), then `sbatch` this — no
# env vars, no dependency.  Compare its sigma_n station traces against the
# over-int run (ab_overint2_1000m_flex.sbatch).
set -euo pipefail

# ---- HARDCODED A/B knobs (this file = baseline, 1000 m) -------------------
FAULT_OVERINT=0
MESH="tpv102/mesh/tpv102_1000m.msh"
TFINAL=5
CONFIG_TOML="tpv102/configs/tpv102_spatial.toml"

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8
export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH:-}"

# locate the repo root (walk up to the dir containing miniapps/seas)
_d="${SEAS_MFEM_ROOT:-${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "${BASH_SOURCE[0]:-$0}")")}}"
while [[ "${_d}" != "/" && ! -d "${_d}/miniapps/seas" ]]; do _d="$(dirname "${_d}")"; done
cd "${_d}/miniapps/seas"

if [[ ! -x ./seas_spatial_dyn_driver ]]; then
    echo "ERROR: ./seas_spatial_dyn_driver not built.  Build it first:"
    echo "       sbatch miniapps/seas/jobs/spatial_dyn_build.sbatch clean-rebuild"
    exit 1
fi
[[ -f "${CONFIG_TOML}" ]] || { echo "ERROR: ${CONFIG_TOML} not found."; exit 1; }
grep -Eq '^[[:space:]]*order[[:space:]]*=[[:space:]]*1([[:space:]]|$)' "${CONFIG_TOML}" \
    || { echo "ERROR: ${CONFIG_TOML} is not [mesh].order=1."; exit 4; }
grep -Eq '^[[:space:]]*interior_flux[[:space:]]*=[[:space:]]*"scalar"' "${CONFIG_TOML}" \
    || { echo "ERROR: ${CONFIG_TOML} is not interior_flux=\"scalar\"."; exit 4; }

# mesh (Gmsh v2.2 only) — pre-build on a login node if gmsh is not in-job
if [[ ! -f "${MESH}" ]]; then
    echo "[ab] mesh ${MESH} missing — generating from ${MESH%.msh}.geo ..."
    command -v gmsh >/dev/null 2>&1 || { echo "ERROR: gmsh not on PATH; run on a login node: gmsh -format msh22 -3 ${MESH%.msh}.geo -o ${MESH}"; exit 2; }
    gmsh -format msh22 -3 "${MESH%.msh}.geo" -o "${MESH}"
fi

OUT="tpv102/out_p1_aderO2_pu_overint${FAULT_OVERINT}_1000m_${SLURM_JOB_ID:-local}"
echo "[ab] BASELINE  FAULT_OVERINT=${FAULT_OVERINT}  mesh=${MESH}  tfinal=${TFINAL}  out=${OUT}"
echo "[ab] Git: $(git rev-parse --abbrev-ref HEAD) $(git rev-parse --short HEAD)"

ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --ader-order 2 \
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
    --verify-dispatch

echo "[ab] done (baseline).  sigma_n station traces: ${OUT}/tpv102_station_*.dat"
