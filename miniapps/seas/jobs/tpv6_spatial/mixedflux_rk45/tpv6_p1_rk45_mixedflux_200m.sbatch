#!/bin/bash
#SBATCH -J tpv6_p1_rk_mf
#SBATCH -o tpv6_p1_rk_mf_%j.out
#SBATCH -e tpv6_p1_rk_mf_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 24:00:00
#SBATCH -A EAR20006

# SCEC TPV6 — ARM 2: BI-MATERIAL FAULT, p1 / RK45 / MIXED (central) flux / matrix
#   Riemann, 200 m (first-pass), STANDARD (asymmetric) mesh.
#
# Arm 2 = central flux on the (uniform-per-side) fault-adjacent corridor + upwind
# bulk + the PER-SIDE bi-material fault Riemann (Part B).  Central is non-
# dissipative -> ADER-unstable at p1 -> RK45 required (existing guard).  On the
# asymmetric mesh, central removes the near-fault upwind-dissipation asymmetry.
# The volume contrast guard (Part A, --mixed-flux-contrast-tol 0.05) is ACTIVE but
# nearly INERT for TPV6: only the few corridor faces that straddle the welded
# bi-material interface (y=0) just BEYOND the fault-rectangle edge are reclassified
# to upwind (the half-space material jumps across the whole y=0 plane, not only the
# fault).  The startup histogram reports that small "reclassified N/M" count.
#   * p1 / RK45     : --time-integrator rk45 (the config defaults to ADER; the CLI
#                     flips it — matrix+mixed+ADER is REJECTED, matrix+mixed+RK OK).
#   * mixed flux    : --mixed-flux adjacent.
#   * contrast guard: --mixed-flux-contrast-tol 0.05.
#   * matrix Riemann: [numerics].interior_flux="matrix" (per-side fault; asserted).
#   Config: tpv6/configs/tpv6.toml.  sigma_n is PHYSICALLY non-constant — not a leak.
set -euo pipefail

MESH="tpv6/mesh/tpv6_200m.msh"                      # asymmetric (standard) mesh
TFINAL="${TPV6_TFINAL:-8.0}"
CONFIG_TOML="tpv6/configs/tpv6.toml"

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
grep -Eq '^[[:space:]]*interior_flux[[:space:]]*=[[:space:]]*"matrix"' "${CONFIG_TOML}" \
    || { echo "ERROR: ${CONFIG_TOML} is not interior_flux=\"matrix\" (TPV6 needs the per-side bi-material fault)."; exit 4; }
grep -Eq '^[[:space:]]*order[[:space:]]*=[[:space:]]*1([[:space:]]|$)' "${CONFIG_TOML}" \
    || { echo "ERROR: ${CONFIG_TOML} is not [mesh].order=1 (this job is p1)."; exit 4; }

# .msh is gitignored — generate from the .geo (Gmsh v2.2) on a login node if absent.
if [[ ! -f "${MESH}" ]]; then
    echo "[tpv6] mesh ${MESH} missing — generating from ${MESH%.msh}.geo (Gmsh v2.2)..."
    command -v gmsh >/dev/null 2>&1 || { echo "ERROR: gmsh not on PATH; pre-build: gmsh -format msh22 -3 ${MESH%.msh}.geo -o ${MESH}"; exit 2; }
    gmsh -format msh22 -3 "${MESH%.msh}.geo" -o "${MESH}"
fi

OUT="tpv6/out_p1_rk45_mixedflux_200m_${SLURM_JOB_ID:-local}"
echo "[tpv6] ARM 2 (matrix + mixed + per-side fault, RK45)  mesh=${MESH}  tfinal=${TFINAL}  out=${OUT}"
echo "[tpv6] Git: $(git rev-parse --abbrev-ref HEAD) $(git rev-parse --short HEAD)"

RESTART_FLAGS=""
if [[ -n "${TPV6_RESTART:-}" ]]; then
    RESTART_FLAGS="--restart ${TPV6_RESTART}"
    echo "[tpv6] RESTART from ${TPV6_RESTART} (this segment writes ${OUT})."
fi

ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --partition-fault-locality \
    --time-integrator rk45 \
    --mixed-flux adjacent \
    --mixed-flux-contrast-tol 0.05 \
    --tfinal "${TFINAL}" \
    --paraview \
    --paraview-fault-vtu \
    --paraview-fault-dt 0.05 \
    --paraview-max-snapshots 300 \
    --print-derived \
    --checkpoint-every 2000 \
    --output-dir "${OUT}" \
    ${RESTART_FLAGS} \
    --verify-dispatch

echo "[tpv6] done (ARM 2: matrix + mixed + per-side fault, RK45, asymmetric mesh)."
echo "[tpv6] verify-dispatch banner: 'BI-MATERIAL fault' + a small 'reclassified N/M' (welded-interface corridor faces)."
