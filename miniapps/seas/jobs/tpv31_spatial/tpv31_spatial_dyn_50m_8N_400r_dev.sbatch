#!/bin/bash
#SBATCH -J tpv31_sp_50m_dev
#SBATCH -o tpv31_sp_50m_dev_%j.out
#SBATCH -e tpv31_sp_50m_dev_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# =========================================================================
# SCEC TPV31 — depth-heterogeneous LSW through seas_spatial_dyn_driver
# (Phase 10).  DEV SMOKE: development queue, 8N/400r/2 h.
#
# TPV31 specifics (vs the SAFS rate-and-state runs):
#   * [meta].law = "slip_weakening" (LSW, NOT rate-state).
#   * [material] kind = "depth_profile_1d" — discontinuous 1-D velocity
#     profile (Mode::Coefficient); REQUIRES interior_flux = "matrix"
#     (the bimaterial Riemann path → BimaterialWaveOperator).  NO velocity
#     sidecar (the TOML sets [velocity].use_sidecar = false), so
#     --no-sidecar-material is NOT needed.
#   * [stress] kind = "depth_proportional" — Cauchy tensor scaled per-DOF by
#     mu(depth)/mu_ref.
#   * depth-linear cohesion C0(depth) = 0.000425 MPa/m * max(0, 2400 - depth).
#   * [nucleation] instantaneous_overstress_circular (cosine taper r=1400→2000).
#     NOTE (Phase 10 known limitation): the spatial driver applies the nucleation
#     overstress UNIFORMLY; the spec's per-DOF mu(depth)/mu_0 amplitude scaling
#     is not yet wired (small near the hypocenter depth where mu ~= mu_0).
#
# Reference overlays: tpv31/benchmark_data/scec_{eqdyna,seisol}/.
# Spec: tpv31/benchmark_document/TPV31_32_Description_v03.pdf.
# =========================================================================

set -euo pipefail

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

# --- pin the hypre/2.31.0 runtime libHYPRE.so (REQUIRED) -------------------
# petsc/3.15 + mumps/5.3 (the old PETSc-3.15 tree) prepend a PRE-2.31 libHYPRE
# to LD_LIBRARY_PATH *after* hypre/2.31.0 is loaded.  The driver was compiled
# against hypre/2.31.0 headers, where HYPRE_Init() is a macro for the renamed
# HYPRE_Initialize(); the binary therefore references symbol HYPRE_Initialize,
# which exists only in hypre >= 2.31.0.  Because the binary uses DT_RUNPATH
# (searched AFTER LD_LIBRARY_PATH), the older PETSc/MUMPS libHYPRE shadows the
# 2.31.0 one at run time -> "symbol lookup error: undefined symbol:
# HYPRE_Initialize" and MPI exit code 127.  Re-assert the 2.31.0 lib dir at the
# FRONT of the search path (matches every bp1_*/bp5_* job in this repo).
export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH:-}"

# --- locate the repo root (walk up to the dir containing miniapps/seas) ----
_walk_up_dir="${SEAS_MFEM_ROOT:-${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "${BASH_SOURCE[0]:-$0}")")}}"
while [[ "${_walk_up_dir}" != "/" && ! -d "${_walk_up_dir}/miniapps/seas" ]]; do
    _walk_up_dir="$(dirname "${_walk_up_dir}")"
done
SEAS_MFEM_ROOT="$(cd "${_walk_up_dir}" && pwd)"
cd "${SEAS_MFEM_ROOT}/miniapps/seas"

if [[ ! -x ./seas_spatial_dyn_driver ]]; then
    echo "ERROR: ./seas_spatial_dyn_driver not found/executable in $(pwd)."
    echo "       Build it on a login/idev node:  cd miniapps/seas && make seas_spatial_dyn_driver"
    exit 1
fi

CONFIG_TOML="${TPV31_CONFIG:-tpv31/configs/tpv31.toml}"
MESH="${TPV31_MESH:-tpv31/mesh/tpv31_50m.msh}"

# --- regenerate the Gmsh v2.2 mesh if it is missing -----------------------
# The .msh is NOT committed (only the .geo).  MFEM's reader is Gmsh v2.2 only
# (see miniapps/seas/CLAUDE.md "Known limitation — Gmsh .msh format"), so the
# `-format msh22` flag is MANDATORY.  Run on a login/idev node (gmsh needs the
# `pythonenv` conda env), NOT inside the compute job, if gmsh is unavailable on
# the compute nodes.
if [[ ! -f "${MESH}" ]]; then
    echo "[tpv31] mesh ${MESH} missing — generating from the .geo (Gmsh v2.2)..."
    if command -v gmsh >/dev/null 2>&1; then
        gmsh -format msh22 -3 tpv31/mesh/tpv31_50m.geo -o "${MESH}"
    else
        echo "ERROR: gmsh not on PATH.  Generate the mesh first (conda activate"
        echo "       pythonenv): gmsh -format msh22 -3 tpv31/mesh/tpv31_50m.geo"
        echo "       -o ${MESH}"
        exit 2
    fi
fi

# --- NO separate serial --dry-run pre-flight (intentional; removed) --------
# The spatial driver's --dry-run does NOT short-circuit at CLI-parse time the
# way tpv102/104 do (tpv102_driver.cpp:1000 exits BEFORE the mesh read at 1017):
# it reads the full ~6M-tet 50 m mesh and builds the entire BimaterialWaveOperator
# + fault tables before its dry-run exit (spatial_dyn_driver.cpp:2085).  Launched
# WITHOUT ibrun (a single rank) that construction never finishes inside the wall
# — the p2 sibling (job 7757861) was SLURM-killed at the 2 h TIME LIMIT mid-dry-
# run.  The ibrun run below does the SAME construction + the SAME C++ config
# guards (FaultIteratorSupported, material-vs-interior_flux) in parallel and
# carries --verify-dispatch, so the resolved scheme is still printed.

OUT="${TPV31_OUT:-tpv31/out_${SLURM_JOB_ID:-local}}"
TFINAL="${TPV31_TFINAL:-15s}"   # spec p. 9 (2 h dev wall reaches a fraction; use restart for the rest)

ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --tfinal "${TFINAL}" \
    --paraview \
    --paraview-fault-hdf5 \
    --paraview-fault-zfp-tol 1.0e-12 \
    --paraview-max-snapshots 300 \
    --print-derived \
    --checkpoint-every 2000 \
    --output-dir "${OUT}" \
    --verify-dispatch

echo "[tpv31] done. On-fault station traces in ${OUT}; compare vs"
echo "        tpv31/benchmark_data/scec_{eqdyna,seisol}/ via tpv31/visualize_results.py."
