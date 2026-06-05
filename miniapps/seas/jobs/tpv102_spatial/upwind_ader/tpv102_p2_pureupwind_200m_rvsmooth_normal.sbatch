#!/bin/bash
#SBATCH -J tpv102_p2_pu_rvsm
#SBATCH -o tpv102_p2_pu_rvsm_%j.out
#SBATCH -e tpv102_p2_pu_rvsm_%j.err
#SBATCH -p normal
#SBATCH -N 12
#SBATCH -n 600
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# TPV102 MESH control on tpv102_200m_rvsmooth.msh — p2 / ADER-O3 / pure upwind /
#   scalar, 200 m, NO over-integration, NO resample, 0-12 s.
# Purpose: isolate the MESH effect at p2.  rvsmooth is an alternative near-fault
# mesh variant (cf. the unstructured tpv102_200m and the y-mirror hybrid); pure
# upwind with no dealiasing knobs, so any on-fault sigma_n change vs the
# unstructured p2 baseline (out_p2_aderO3_pu_overint0_200m_*) is attributable to
# the mesh alone.  Pairs with the p1 rvsmooth control.
# Standalone: build seas_spatial_dyn_driver first (spatial_dyn_build.sbatch).
# NOTE: p2 dt is smaller than p1; even 48 h likely won't reach 12 s — it
# checkpoints every 2000 steps; station traces append every step, so partial
# output already shows the early-rupture speckle.  Resubmit with
# TPV102_RESTART=<out>/cp to continue.
set -euo pipefail

# ---- knobs (pure-upwind baseline on the rvsmooth mesh, p2) -------------
FAULT_OVERINT=0                                      # no over-integration
MESH="tpv102/mesh/tpv102_200m_rvsmooth.msh"
TFINAL=12
CONFIG_TOML="tpv102/configs/tpv102_spatial_p2.toml"  # p2, interior_flux="scalar"

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

# rvsmooth has NO .geo to regenerate from (built offline, ~130 MB) — it must be
# staged on the cluster (scp/rsync) before submitting.
[[ -f "${MESH}" ]] || { echo "ERROR: mesh ${MESH} not found; it has no .geo to regenerate from — stage it (scp/rsync) on the cluster first."; exit 2; }

OUT="tpv102/out_p2_aderO3_pu_200m_rvsmooth_${SLURM_JOB_ID:-local}"
echo "[rvsmooth] p2 pure-upwind on rvsmooth mesh  mesh=${MESH}  tfinal=${TFINAL}  out=${OUT}"
echo "[rvsmooth] Git: $(git rev-parse --abbrev-ref HEAD) $(git rev-parse --short HEAD)"

RESTART_FLAGS=""
if [[ -n "${TPV102_RESTART:-}" ]]; then
    RESTART_FLAGS="--restart ${TPV102_RESTART}"
    echo "[rvsmooth] RESTART from ${TPV102_RESTART} (this segment writes ${OUT})."
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

echo "[rvsmooth] done (p2 pure-upwind, rvsmooth mesh).  sigma_n station traces: ${OUT}/tpv102_station_*.dat"
echo "[rvsmooth] compare vs unstructured-mesh p2 baseline:"
echo "  python3 jobs/tpv102_spatial/compare_sigma_n.py <out_..._p2_overint0_200m_unstructured> ${OUT}"
