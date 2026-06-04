#!/bin/bash
#SBATCH -J tpv102_p1_pu_hybrid
#SBATCH -o tpv102_p1_pu_hybrid_%j.out
#SBATCH -e tpv102_p1_pu_hybrid_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 06:00:00
#SBATCH -A EAR20006

# TPV102 HYBRID-MESH control — p1 / ADER-O2 / pure upwind / scalar, 200 m,
#   NO over-integration, NO resample, 0-12 s.
# Purpose: isolate the MESH effect.  tpv102_200m_hybrid.msh has a y-mirror
# structured prism strip near the fault (Pelties 2012 "symmetric mesh"), so the
# +y/-y near-fault halves are exact mirror partners — removing the near-fault
# mesh asymmetry that the plan (§3.4) identifies as the speckle seed.  Compare
# its on-fault sigma_n against the same scheme on the UNSTRUCTURED
# tpv102_200m.msh (out_p1_aderO2_pu_overint0_200m_*): if the symmetric mesh
# alone flattens sigma_n, the speckle is mesh-asymmetry-driven (independent of
# the over-int/resample dealiasing route).
# Standalone: build seas_spatial_dyn_driver first (spatial_dyn_build.sbatch).
# NOTE: 6 h likely will NOT reach 12 s on this (larger) hybrid mesh — it
# checkpoints every 2000 steps; the station traces append every step, so partial
# output is already enough to see the early-rupture speckle.  Resubmit with
# TPV102_RESTART=<out>/cp to continue.
set -euo pipefail

# ---- knobs (pure-upwind baseline on the HYBRID mesh) --------------------
FAULT_OVERINT=0                                   # no over-integration
MESH="tpv102/mesh/tpv102_200m_hybrid.msh"
TFINAL=12
CONFIG_TOML="tpv102/configs/tpv102_spatial.toml"  # p1, interior_flux="scalar"

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
grep -Eq '^[[:space:]]*order[[:space:]]*=[[:space:]]*1([[:space:]]|$)' "${CONFIG_TOML}" \
    || { echo "ERROR: ${CONFIG_TOML} is not [mesh].order=1 (this job is p1)."; exit 4; }
grep -Eq '^[[:space:]]*interior_flux[[:space:]]*=[[:space:]]*"scalar"' "${CONFIG_TOML}" \
    || { echo "ERROR: ${CONFIG_TOML} is not interior_flux=\"scalar\"."; exit 4; }

# The hybrid mesh is large (~170 MB) — prefer to pre-generate it on a login node
# (conda activate pythonenv) rather than here; fall back to gmsh if missing.
if [[ ! -f "${MESH}" ]]; then
    echo "[hybrid] mesh ${MESH} missing — generating from ${MESH%.msh}.geo (Gmsh v2.2)..."
    command -v gmsh >/dev/null 2>&1 || { echo "ERROR: gmsh not on PATH; pre-build on a login node: gmsh -format msh22 -3 ${MESH%.msh}.geo -o ${MESH}"; exit 2; }
    gmsh -format msh22 -3 "${MESH%.msh}.geo" -o "${MESH}"
fi

OUT="tpv102/out_p1_aderO2_pu_200m_hybrid_${SLURM_JOB_ID:-local}"
echo "[hybrid] p1 pure-upwind baseline on HYBRID mesh  mesh=${MESH}  tfinal=${TFINAL}  out=${OUT}"
echo "[hybrid] Git: $(git rev-parse --abbrev-ref HEAD) $(git rev-parse --short HEAD)"

RESTART_FLAGS=""
if [[ -n "${TPV102_RESTART:-}" ]]; then
    RESTART_FLAGS="--restart ${TPV102_RESTART}"
    echo "[hybrid] RESTART from ${TPV102_RESTART} (this segment writes ${OUT})."
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

echo "[hybrid] done (p1 pure-upwind, hybrid mesh).  sigma_n station traces: ${OUT}/tpv102_station_*.dat"
echo "[hybrid] compare vs unstructured-mesh baseline:"
echo "  python3 jobs/tpv102_spatial/compare_sigma_n.py <out_..._overint0_200m_unstructured> ${OUT}"
