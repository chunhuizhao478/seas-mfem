#!/bin/bash
#SBATCH -J tpv102_p1_rk45_mf
#SBATCH -o tpv102_p1_rk45_mf_%j.out
#SBATCH -e tpv102_p1_rk45_mf_%j.err
#SBATCH -p flex
#SBATCH -N 16
#SBATCH -n 800
#SBATCH -t 01:00:00
#SBATCH -A EAR20006

# =========================================================================
# SCEC TPV102 — p1 · SCALAR mixed (central) flux · explicit RK45
# Driver: seas_spatial_dyn_driver (spatial dynamic-rupture).
# DEV SMOKE: development queue, 8N/400r/2 h.
#
# NUMERICAL SCHEME (the four knobs that define this job) ------------------
#   * p1            : [mesh].order = 1            (config-only, asserted below).
#   * mixed flux    : --mixed-flux adjacent       (central flux on the
#                     fault-adjacent faces; scalar Godunov elsewhere).
#   * RK45          : --time-integrator rk45       (Dormand-Prince RK45; NOT
#                     ADER).  At p1 the central flux is UNSTABLE under ADER-O2
#                     (predictor locked to O2 amplifies the imaginary-axis
#                     modes central flux creates); an explicit RK of order >= 3
#                     has imaginary-axis coverage and is stable.  See
#                     document/mixed_flux_dev/
#                     BUILD_mixed_flux_rk_dynamic_rupture_2026-05-28.md.
#   * scalar Riemann: [numerics].interior_flux = "scalar"  (config-only,
#                     asserted below).  MANDATORY for the RK path (mixed/central
#                     flux is scalar-only; matrix/bimaterial is rejected).
#
#   The CLI flags re-assert what tpv102_spatial_rk45_mixedflux.toml already
#   sets, so the resolved scheme is printed in the driver's startup banner
#   (grep the *.out for "time integrator:", "mixed flux:", "law:") and cannot
#   silently drift from the config.  [mesh].order and [numerics].interior_flux
#   have NO CLI override, so they are pre-flight-asserted against the TOML.
#
# TPV102 physics ---------------------------------------------------------
#   * [meta].law = "rate_state" (Dieterich-Ruina, AGING law).  The RK path
#     runs the coupled RK-on-(Q, psi, slip) stepper AdvanceRKCoupled_Spatial,
#     whose wave.Mult fault face uses the instantaneous FaultFaceFlux::Evaluate
#     (same RS Riemann solve as the ADER substep).
#   * [stress] kind = "fault_local_prestress": tau_strike = 75 MPa,
#     sigma_n = 120 MPa.  VW core a=0.008<b=0.012 inside 30x15 km, VS (a=0.016)
#     outside.  Dc = 0.02 m.
#   * Nucleation: compact-bell gradual overstress, delta_tau = 25 MPa, T_nuc = 1 s.
#   * Homogeneous material -> interior_flux = "scalar".
#
# Reference overlays: tpv102 benchmark + drivers/tpv102_driver.cpp.
#
# CFL is a STARTING calibration (RK-aware Adjacent factor 0.6 x cfl 0.25);
# read the empirical stable dt off this job's *.out BEFORE committing
# normal-queue hours.  To promote to a production run see "PRODUCTION
# PROMOTION" below.
#
# LOCAL VERIFICATION SCOPE -----------------------------------------------
#   The RK code is verified locally by compile + unit tests
#   (seas_test_rk_time_stepper, seas_test_fault_face_flux_ader_equivalence,
#   seas_test_wave_operator).  The driver's END-TO-END rk45 path PAST mesh
#   load is NOT exercised locally (no-local-full-mesh rule).  The
#   `--dry-run --verify-dispatch` step BELOW (which passes --mesh) is the
#   first execution that loads a real mesh + builds the operator + resolves
#   the dispatch; it runs HERE on Frontera/idev as the gate.  Confirm it
#   exits 0 before the long run.
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

CONFIG_TOML="${TPV102_CONFIG:-tpv102/configs/tpv102_spatial_rk45_mixedflux.toml}"
# 200 m mesh (real work for 400 ranks; the TOML default is the 1000 m mesh).
MESH="${TPV102_MESH:-tpv102/mesh/tpv102_200m.msh}"
GEO="${TPV102_GEO:-tpv102/mesh/tpv102_200m.geo}"

if [[ ! -f "${CONFIG_TOML}" ]]; then
    echo "ERROR: config ${CONFIG_TOML} not found (cwd $(pwd))."; exit 1
fi

# --- pre-flight: assert the two CONFIG-ONLY knobs (no CLI override) --------
if ! grep -Eq '^[[:space:]]*order[[:space:]]*=[[:space:]]*1([[:space:]]|$)' "${CONFIG_TOML}"; then
    echo "ERROR: ${CONFIG_TOML} does not set [mesh].order = 1 (this job is p1)."
    exit 4
fi
if ! grep -Eq '^[[:space:]]*interior_flux[[:space:]]*=[[:space:]]*"scalar"' "${CONFIG_TOML}"; then
    echo "ERROR: ${CONFIG_TOML} does not set [numerics].interior_flux = \"scalar\"."
    echo "       The RK/mixed-flux path is scalar-only (matrix/bimaterial rejected)."
    exit 4
fi

# --- regenerate the Gmsh v2.2 mesh if it is missing -----------------------
# MFEM's reader is Gmsh v2.2 only (see miniapps/seas/CLAUDE.md "Known
# limitation — Gmsh .msh format"); the `-format msh22` flag is MANDATORY.
if [[ ! -f "${MESH}" ]]; then
    echo "[tpv102-rk45] mesh ${MESH} missing — generating from ${GEO} (Gmsh v2.2)..."
    if command -v gmsh >/dev/null 2>&1; then
        gmsh -format msh22 -3 "${GEO}" -o "${MESH}"
    else
        echo "ERROR: gmsh not on PATH.  Generate the mesh first (conda activate"
        echo "       pythonenv): gmsh -format msh22 -3 ${GEO} -o ${MESH}"
        exit 2
    fi
fi

# --- NO separate serial --dry-run pre-flight (intentional; removed) --------
# The spatial driver's --dry-run reads the full mesh + builds the entire
# WaveOperator + fault tables BEFORE its dry-run exit (spatial_dyn_driver.cpp),
# so launched WITHOUT ibrun (a single rank) it eats the wall and can be
# SLURM-killed before the real run ever starts (this is why the tpv31 jobs
# dropped it, commit 6ef7daf).  The ibrun run below does the SAME construction
# + the SAME C++ dispatch guards in parallel and carries --verify-dispatch,
# which PRINTS the resolved scheme (grep the *.out for "time integrator:",
# "mixed flux:", "law:", "interior flux:") and then CONTINUES into the full
# run — it does NOT exit the way --dry-run does.  The TOML greps above already
# fail-fast on a drifted interior_flux / [mesh].order.

OUT="${TPV102_OUT:-tpv102/out_rk45_mf_${SLURM_JOB_ID:-local}}"
TFINAL="${TPV102_TFINAL:-12}"   # SCEC TPV102 window 0-12 s; the 2 h dev wall
                                # reaches a fraction — restart for the rest.

# --- run -----------------------------------------------------------------
ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --time-integrator rk45 \
    --mixed-flux adjacent \
    --tfinal "${TFINAL}" \
    --paraview \
    --paraview-fault-vtu \
    --paraview-fault-dt 0.1 \
    --paraview-max-snapshots 300 \
    --print-derived \
    --checkpoint-every 2000 \
    --output-dir "${OUT}" \
    --verify-dispatch

echo "[tpv102-rk45] done. On-fault SCEC station traces: ${OUT}/tpv102_station_*.dat (9 stations); compare vs the SCEC TPV102"
echo "              reference.  Confirm V_max is BOUNDED (no ADER-p1-central-flux"
echo "              run-away) and read the empirical stable dt off this *.out"
echo "              before committing normal-queue hours."

# =========================================================================
# PRODUCTION PROMOTION (12 s to completion)
# -------------------------------------------------------------------------
#   1. Swap the queue + wall:  #SBATCH -p normal ;  #SBATCH -t 24:00:00
#      (scale -N/-n as needed; Frontera CLX = 56 cores/node).
#   2. Chain restarts across the wall limit: the run checkpoints every 2000
#      steps to ${OUT}.  Resume with:  --restart "${OUT}/cp".
#   3. Pre-build the mesh ONCE on a login/idev node (conda activate pythonenv;
#      gmsh -format msh22 -3 ${GEO} -o ${MESH}) so the compute job does not mesh.
#   4. Recalibrate dt from the dev smoke's empirical steps/hour (the RK-aware
#      mixed-flux CFL factor 0.6 for Adjacent is an interim placeholder).
# =========================================================================
