#!/bin/bash
#SBATCH -J tpv7_p1_ader_sym
#SBATCH -o tpv7_p1_ader_sym_%j.out
#SBATCH -e tpv7_p1_ader_sym_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 12:00:00
#SBATCH -A EAR20006

# SCEC TPV7 (low-contrast bi-material) — ARM 1: BI-MATERIAL FAULT, p1 / ADER-O2 / PURE UPWIND / matrix
#   Riemann, 200 m (first-pass), y-mirror SYMMETRIC mesh.
#
# Arm 1 validates that matrix + upwind + the PER-SIDE bi-material fault Riemann
# (unified-bimaterial plan Part B) runs under ADER on a symmetric mesh.  The
# symmetric near-fault prism strip removes the +y/-y upwind-dissipation asymmetry;
# the per-side fault Riemann carries the (physical) bi-material coupling.
#   * p1            : [mesh].order = 1            (config-only; asserted below).
#   * ADER-O2       : --ader-order 2 --time-integrator ader (the config default).
#   * pure upwind   : --mixed-flux none --fault-overint 0 (no central flux).
#   * matrix Riemann: [numerics].interior_flux = "matrix" (MANDATORY for the
#                     per-side bi-material fault; config-only, asserted below).
#   Config: tpv7/configs/tpv7.toml ([material] halfspace_across_fault; nucleation
#           is the STATIC stress patch tau_strike=81.6 MPa, like TPV205 — no
#           [nucleation] block).  sigma_n is PHYSICALLY non-constant (bi-material)
#           — do NOT treat its on-fault variation as a leak.
# Standalone: build seas_spatial_dyn_driver first (from THIS branch — Part B
# per-side fault + B3 halfspace material).
set -euo pipefail

FAULT_OVERINT=0
MESH="tpv6/mesh/tpv6_200m_symmetric.msh"
TFINAL="${TPV7_TFINAL:-8.0}"                        # TPV6 rupture ~ a few s
CONFIG_TOML="tpv7/configs/tpv7.toml"

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8
# Pin hypre/2.31.0 at the FRONT (petsc/3.15 + mumps/5.3 shadow HYPRE_Initialize).
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
grep -Eq '^[[:space:]]*interior_flux[[:space:]]*=[[:space:]]*"matrix"' "${CONFIG_TOML}" \
    || { echo "ERROR: ${CONFIG_TOML} is not interior_flux=\"matrix\" (TPV6 needs the per-side bi-material fault)."; exit 4; }
grep -Eq '^[[:space:]]*order[[:space:]]*=[[:space:]]*1([[:space:]]|$)' "${CONFIG_TOML}" \
    || { echo "ERROR: ${CONFIG_TOML} is not [mesh].order=1 (this job is p1)."; exit 4; }

# .msh is gitignored — generate from the .geo (Gmsh v2.2) on a login node if absent.
if [[ ! -f "${MESH}" ]]; then
    echo "[tpv7] mesh ${MESH} missing — generating from ${MESH%.msh}.geo (Gmsh v2.2)..."
    command -v gmsh >/dev/null 2>&1 || { echo "ERROR: gmsh not on PATH; pre-build: gmsh -format msh22 -3 ${MESH%.msh}.geo -o ${MESH}"; exit 2; }
    gmsh -format msh22 -3 "${MESH%.msh}.geo" -o "${MESH}"
fi

OUT="tpv7/out_p1_aderO2_upwind_200m_symmetric_${SLURM_JOB_ID:-local}"
echo "[tpv7] ARM 1 (matrix + upwind + per-side fault, ADER)  mesh=${MESH}  tfinal=${TFINAL}  out=${OUT}"
echo "[tpv7] Git: $(git rev-parse --abbrev-ref HEAD) $(git rev-parse --short HEAD)"

RESTART_FLAGS=""
if [[ -n "${TPV7_RESTART:-}" ]]; then
    RESTART_FLAGS="--restart ${TPV7_RESTART}"
    echo "[tpv7] RESTART from ${TPV7_RESTART} (this segment writes ${OUT})."
fi

ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --partition-fault-locality \
    --ader-order 2 \
    --time-integrator ader \
    --mixed-flux none \
    --fault-overint ${FAULT_OVERINT} \
    --tfinal "${TFINAL}" \
    --paraview \
    --paraview-fault-vtu \
    --paraview-fault-dt 0.05 \
    --paraview-max-snapshots 300 \
    --print-derived \
    --checkpoint-every 2000 \
    --output-dir "${OUT}" \
    --deriv-cache --shared-ck-recursion \
    ${RESTART_FLAGS} \
    --verify-dispatch

echo "[tpv7] done (ARM 1: matrix + upwind + per-side fault, ADER, symmetric mesh)."
echo "[tpv7] verify-dispatch banner must show 'BI-MATERIAL fault' (Zp_plus != Zp_minus)."
