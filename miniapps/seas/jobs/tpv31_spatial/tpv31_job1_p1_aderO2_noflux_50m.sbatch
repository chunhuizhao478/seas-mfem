#!/bin/bash
#SBATCH -J tpv31_j1_p1_aderO2
#SBATCH -o tpv31_j1_p1_aderO2_%j.out
#SBATCH -e tpv31_j1_p1_aderO2_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# =========================================================================
# SCEC TPV31 — JOB [1]:  p1  ·  no mixed flux  ·  ADER-O2  ·  matrix Riemann
# Driver: seas_spatial_dyn_driver (spatial dynamic-rupture, Phase 10).
#
# NUMERICAL SCHEME (the four knobs that define this job) ------------------
#   * p1            : [mesh].order = 1            (linear elements; config-only,
#                     no CLI override — asserted below).
#   * no mixed flux : --mixed-flux none           (upwind Godunov on every
#                     interior face; no central-flux faces).
#   * ADER-O2       : --ader-order 2 --time-integrator ader
#                     (2nd-order ADER-DG sub-step path; the spatial driver
#                     always sub-steps with O = ader_order, fault_iterator =
#                     "substep").
#   * matrix Riemann: [numerics].interior_flux = "matrix"  (config-only, no
#                     CLI override — asserted below).  This is MANDATORY for
#                     TPV31: its [material].kind = "depth_profile_1d" is a
#                     Mode::Coefficient (bi-material) field, which is mutually
#                     exclusive with the scalar Godunov ctor.  The matrix path
#                     instantiates BimaterialWaveOperator (per-element flux
#                     pool + per-face bi-material Riemann solve).
#
#   The CLI flags re-assert what tpv31.toml already sets, so the resolved
#   scheme is printed in the driver's startup banner (grep the *.out for
#   "ader order:", "time integrator:", "mixed flux:") and cannot silently
#   drift from the config.  interior_flux and [mesh].order have NO CLI
#   override, so they are pre-flight-asserted against the TOML below.
#
# TPV31 physics (vs the SAFS rate-and-state runs) ------------------------
#   * [meta].law = "slip_weakening"  (LSW, NOT rate-state).
#   * [material] kind = "depth_profile_1d" — discontinuous 1-D velocity
#     profile (Mode::Coefficient) ⇒ requires interior_flux = "matrix".
#     NO velocity sidecar ([velocity].use_sidecar = false), so
#     --no-sidecar-material is NOT needed.
#   * [stress] kind = "depth_proportional" — Cauchy tensor scaled per-DOF by
#     mu(depth)/mu_ref.
#   * depth-linear cohesion C0(depth) = 0.000425 MPa/m * max(0, 2400 - depth).
#   * [nucleation] instantaneous_overstress_circular (cosine taper r=1400→2000).
#
# Reference overlays: tpv31/benchmark_data/scec_{eqdyna,seisol}/.
# Spec: tpv31/benchmark_document/TPV31_32_Description_v03.pdf.
#
# This header job is the DEV SMOKE (development queue, 8N/400r/2 h).  To
# promote to a production run see "PRODUCTION PROMOTION" at the bottom.
# =========================================================================

set -euo pipefail

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

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

if [[ ! -f "${CONFIG_TOML}" ]]; then
    echo "ERROR: config ${CONFIG_TOML} not found (cwd $(pwd))."; exit 1
fi

# --- pre-flight: assert the two CONFIG-ONLY knobs (no CLI override) --------
# interior_flux="matrix" is MANDATORY for TPV31 (bi-material depth profile),
# and [mesh].order=1 is the "p1" in this job's name.  Both are config-only,
# so a drifted TOML would silently change the scheme.  Fail loudly here.
if ! grep -Eq '^[[:space:]]*interior_flux[[:space:]]*=[[:space:]]*"matrix"' "${CONFIG_TOML}"; then
    echo "ERROR: ${CONFIG_TOML} does not set [numerics].interior_flux = \"matrix\"."
    echo "       TPV31 (depth_profile_1d bi-material) REQUIRES the matrix Riemann path."
    exit 4
fi
if ! grep -Eq '^[[:space:]]*order[[:space:]]*=[[:space:]]*1([[:space:]]|$)' "${CONFIG_TOML}"; then
    echo "ERROR: ${CONFIG_TOML} does not set [mesh].order = 1 (this job is p1)."
    exit 4
fi

# --- regenerate the Gmsh v2.2 mesh if it is missing -----------------------
# The .msh is NOT committed (only the .geo).  MFEM's reader is Gmsh v2.2 only
# (see miniapps/seas/CLAUDE.md "Known limitation — Gmsh .msh format"), so the
# `-format msh22` flag is MANDATORY.  Run on a login/idev node (gmsh needs the
# `pythonenv` conda env) if gmsh is unavailable on the compute nodes.
if [[ ! -f "${MESH}" ]]; then
    echo "[tpv31-j1] mesh ${MESH} missing — generating from the .geo (Gmsh v2.2)..."
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
# carries --verify-dispatch, so the resolved scheme is still printed.  The TOML
# greps above already fail-fast on a drifted interior_flux / [mesh].order.

OUT="${TPV31_OUT:-tpv31/out_j1_${SLURM_JOB_ID:-local}}"
TFINAL="${TPV31_TFINAL:-15}"   # spec p. 9 = 15 s (2 h dev wall reaches a
                               # fraction; checkpoint/restart for the rest)

# --- run -----------------------------------------------------------------
# Scheme pinned on the CLI (--ader-order 2 / --time-integrator ader /
# --mixed-flux none); matrix Riemann + p1 come from the (asserted) TOML.
ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --ader-order 2 \
    --time-integrator ader \
    --mixed-flux none \
    --tfinal "${TFINAL}" \
    --paraview \
    --paraview-fault-vtu \
    --paraview-fault-dt 0.1 \
    --paraview-max-snapshots 300 \
    --print-derived \
    --checkpoint-every 2000 \
    --output-dir "${OUT}" \
    --verify-dispatch

echo "[tpv31-j1] done. On-fault station traces in ${OUT}; compare vs"
echo "           tpv31/benchmark_data/scec_{eqdyna,seisol}/ via tpv31/visualize_results.py."

# =========================================================================
# PRODUCTION PROMOTION (15 s to completion)
# -------------------------------------------------------------------------
# The header above is a development-queue smoke (2 h wall rarely reaches the
# full 15 s on the few-million-tet 50 m mesh).  To run to completion:
#   1. Swap the queue + wall + node count:
#        #SBATCH -p normal
#        #SBATCH -N 16          # (scale ranks: -n = 56 * N on Frontera CLX)
#        #SBATCH -n 896
#        #SBATCH -t 24:00:00    # normal-queue max is 48 h
#   2. Chain restarts across the wall limit: the run writes checkpoints every
#      2000 steps (--checkpoint-every 2000) to ${OUT}.  Resume a truncated run
#      with:  --restart "${OUT}/cp"   (mirrors tpv104 restart jobs).
#   3. Pre-build the mesh ONCE on a login/idev node (conda activate pythonenv;
#      gmsh -format msh22 -3 tpv31/mesh/tpv31_50m.geo -o tpv31/mesh/tpv31_50m.msh)
#      so the compute job does not spend wall time meshing.
# Sizing the wall accurately needs an empirical steps/hour from the first dev
# smoke — read the cadence off this job's *.out before committing normal-queue
# hours.
# =========================================================================
