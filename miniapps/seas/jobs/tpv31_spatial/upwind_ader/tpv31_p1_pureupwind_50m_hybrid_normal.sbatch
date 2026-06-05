#!/bin/bash
#SBATCH -J tpv31_p1_pu_hybrid
#SBATCH -o tpv31_p1_pu_hybrid_%j.out
#SBATCH -e tpv31_p1_pu_hybrid_%j.err
#SBATCH -p normal
#SBATCH -N 12
#SBATCH -n 600
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# SCEC TPV31 HYBRID-MESH control — p1 / ADER-O2 / pure upwind / matrix Riemann,
#   50 m, NO over-integration, NO resample, 0-15 s.
# Purpose: isolate the MESH effect on the TPV31 LSW blow-up.  tpv31_50m_hybrid.msh
# has a y-mirror structured prism strip near the fault (Pelties 2012 "symmetric
# mesh"), so the +y/-y near-fault halves are exact mirror partners — removing the
# near-fault mesh asymmetry (the upwind-flux / Zhang-2023 channel).  Compare its
# on-fault sigma_n against the same scheme on the UNSTRUCTURED tpv31_50m.msh: if
# the symmetric mesh alone tames the drift/blow-up, it is mesh-asymmetry-driven.
#   * p1            : [mesh].order = 1            (config-only; asserted below).
#   * ADER-O2       : --ader-order 2 --time-integrator ader.
#   * pure upwind   : --mixed-flux none --fault-overint 0  (no central flux, no
#                     over-integration; --fault-resample NOT passed).
#   * matrix Riemann: [numerics].interior_flux = "matrix"  (MANDATORY for TPV31's
#                     depth_profile_1d bi-material; config-only, asserted below).
#   Config: tpv31/configs/tpv31.toml.
# Standalone: build seas_spatial_dyn_driver first (from THIS branch — carries the
# ADER deriv-cache + the over-int/resample code).
# NOTE: even 48 h likely won't reach 15 s — checkpoints every 2000 steps; resubmit
# with TPV31_RESTART=<out>/cp to continue.  Station traces append every step, so
# partial output already shows the early-rupture speckle.
set -euo pipefail

FAULT_OVERINT=0                                   # no over-integration
MESH="tpv31/mesh/tpv31_50m_hybrid.msh"
TFINAL="${TPV31_TFINAL:-15}"                       # spec p. 9 = 15 s
CONFIG_TOML="tpv31/configs/tpv31.toml"             # p1, interior_flux="matrix"

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8
# Pin hypre/2.31.0 at the FRONT (petsc/3.15 + mumps/5.3 prepend a pre-2.31
# libHYPRE that shadows HYPRE_Initialize otherwise — see the noflux TPV31 jobs).
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
    || { echo "ERROR: ${CONFIG_TOML} is not interior_flux=\"matrix\" (TPV31 needs the matrix Riemann path)."; exit 4; }
grep -Eq '^[[:space:]]*order[[:space:]]*=[[:space:]]*1([[:space:]]|$)' "${CONFIG_TOML}" \
    || { echo "ERROR: ${CONFIG_TOML} is not [mesh].order=1 (this job is p1)."; exit 4; }

# The hybrid mesh is large (~350 MB) — prefer to pre-generate it on a login node
# (conda activate pythonenv); it HAS a .geo, so fall back to gmsh if missing.
if [[ ! -f "${MESH}" ]]; then
    echo "[hybrid] mesh ${MESH} missing — generating from ${MESH%.msh}.geo (Gmsh v2.2)..."
    command -v gmsh >/dev/null 2>&1 || { echo "ERROR: gmsh not on PATH; pre-build on a login node: gmsh -format msh22 -3 ${MESH%.msh}.geo -o ${MESH}"; exit 2; }
    gmsh -format msh22 -3 "${MESH%.msh}.geo" -o "${MESH}"
fi

OUT="tpv31/out_p1_aderO2_pu_50m_hybrid_${SLURM_JOB_ID:-local}"
echo "[hybrid] TPV31 p1 pure-upwind on HYBRID mesh  mesh=${MESH}  tfinal=${TFINAL}  out=${OUT}"
echo "[hybrid] Git: $(git rev-parse --abbrev-ref HEAD) $(git rev-parse --short HEAD)"

RESTART_FLAGS=""
if [[ -n "${TPV31_RESTART:-}" ]]; then
    RESTART_FLAGS="--restart ${TPV31_RESTART}"
    echo "[hybrid] RESTART from ${TPV31_RESTART} (this segment writes ${OUT})."
fi

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
    ${RESTART_FLAGS} \
    --verify-dispatch

echo "[hybrid] done (TPV31 p1 pure-upwind, hybrid mesh).  station traces: ${OUT}/"
echo "[hybrid] compare vs unstructured-mesh baseline (tpv31_50m) + scec_{eqdyna,seisol} overlays."
