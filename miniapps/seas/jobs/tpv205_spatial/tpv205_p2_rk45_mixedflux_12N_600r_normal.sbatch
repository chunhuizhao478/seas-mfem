#!/bin/bash
#SBATCH -J tpv205_p2_rk45_mf
#SBATCH -o tpv205_p2_rk45_mf_%j.out
#SBATCH -e tpv205_p2_rk45_mf_%j.err
#SBATCH -p normal
#SBATCH -N 12
#SBATCH -n 600
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# =========================================================================
# SCEC TPV205 — p2 · SCALAR mixed (central) flux · explicit RK45
# Driver: seas_spatial_dyn_driver (spatial dynamic-rupture).
# PRODUCTION: normal queue, 12 nodes / 600 MPI ranks (50 ranks/node on
# Frontera CLX 56-core nodes), 48 h wall (normal-queue maximum).
#
# NUMERICAL SCHEME (the four knobs that define this job) ------------------
#   * p2            : [mesh].order = 2            (quadratic FE on the same
#                     200 m mesh; config-only, asserted below).  The RK CFL
#                     de-rating auto-scales: RkCflFactor = 3/(2p+1) = 0.6 at
#                     p=2 (vs 1.0 at p=1), so the stable dt shrinks.
#   * mixed flux    : --mixed-flux adjacent       (central flux on the
#                     fault-adjacent faces; scalar Godunov elsewhere).
#   * RK45          : --time-integrator rk45       (Dormand-Prince RK45; NOT
#                     ADER).  The central flux is UNSTABLE under ADER-O2 (any p);
#                     an explicit RK of order >= 3 has imaginary-axis coverage
#                     and is stable.  See document/mixed_flux_dev/
#                     BUILD_mixed_flux_rk_dynamic_rupture_2026-05-28.md.
#   * scalar Riemann: [numerics].interior_flux = "scalar"  (config-only,
#                     asserted below).  MANDATORY for the RK path.
#
#   The CLI flags re-assert what tpv205_spatial_p2_rk45_mixedflux.toml already
#   sets, so the resolved scheme is printed in the driver's startup banner
#   (grep the *.out for "time integrator:", "mixed flux:", "law:").
#   [mesh].order and [numerics].interior_flux have NO CLI override, so they
#   are pre-flight-asserted against the TOML.
#
# TPV205 physics ---------------------------------------------------------
#   * [meta].law = "slip_weakening"  (LSW, NOT rate-state).  The RK path runs
#     the instantaneous FaultFaceFlux::EvaluateLSW on each Mult + the coupled
#     RK-on-(Q,slip) stepper AdvanceRKCoupledLSW_Spatial.
#   * [stress] kind = "fault_local_prestress": background tau_strike = 70 MPa,
#     sigma_n = 120 MPa, + three SCEC patches (nucleation +81.6 MPa).
#   * LSW (mu_s=0.677, mu_d=0.525, d_c=0.40 m), barrier outside 30x15 km.
#   * Homogeneous material -> interior_flux = "scalar".
#
# Reference overlays: tpv205 benchmark + drivers/tpv205_driver.cpp.
#
# CFL is a STARTING calibration (RK-aware Adjacent factor 0.6 x cfl 0.25);
# read the empirical stable dt off this run's *.out and the steps/hour to
# size the restart chain (see RESTART CHAIN below).
#
# LOCAL VERIFICATION SCOPE -----------------------------------------------
#   The RK+LSW code is verified locally by compile + unit tests
#   (seas_test_lsw_rk_mixed_flux L1-L7 + seas_test_lsw_rk_shared_fault_mpi L8).
#   The driver's END-TO-END rk45+LSW path PAST mesh load is NOT exercised
#   locally (no-local-full-mesh rule).  The `--dry-run --verify-dispatch`
#   step BELOW (which passes --mesh) is the first execution that loads a real
#   mesh + builds the operator + resolves the dispatch; it runs HERE on
#   Frontera/idev as the gate.  Confirm it exits 0 with the
#   "[time-integrator] ... [LSW coupled-RK on (Q, slip)]" banner line.
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

CONFIG_TOML="${TPV205_CONFIG:-tpv205/configs/tpv205_spatial_p2_rk45_mixedflux.toml}"
MESH="${TPV205_MESH:-tpv205/mesh/tpv2053d_200m.msh}"
GEO="${TPV205_GEO:-tpv205/mesh/tpv2053d_200m.geo}"

if [[ ! -f "${CONFIG_TOML}" ]]; then
    echo "ERROR: config ${CONFIG_TOML} not found (cwd $(pwd))."; exit 1
fi

# --- pre-flight: assert the two CONFIG-ONLY knobs (no CLI override) --------
if ! grep -Eq '^[[:space:]]*order[[:space:]]*=[[:space:]]*2([[:space:]]|$)' "${CONFIG_TOML}"; then
    echo "ERROR: ${CONFIG_TOML} does not set [mesh].order = 2 (this job is p2)."
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
    echo "[tpv205-rk45] mesh ${MESH} missing — generating from ${GEO} (Gmsh v2.2)..."
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

OUT="${TPV205_OUT:-tpv205/out_p2_rk45_mf_normal_${SLURM_JOB_ID:-local}}"
TFINAL="${TPV205_TFINAL:-12}"   # SCEC TPV205 window 0-12 s; the 48 h wall may
                                # not finish 12 s at p2 — restart for the rest.

# --- restart across the 48 h wall (opt-in via TPV205_RESTART) ----------------------
# The driver ABORTS (return 3) if --output-dir resolves to the SAME directory
# as the PARENT of --restart, so a restart segment MUST write to a DIFFERENT
# --output-dir than the checkpoint it reads.  OUT defaults to a per-job-ID
# directory, so each segment is naturally distinct.  To resume:
#       TPV205_RESTART="<previous-segment-OUT>/cp" sbatch <thisscript>
RESTART_FLAGS=""
if [[ -n "${TPV205_RESTART:-}" ]]; then
    RESTART_FLAGS="--restart ${TPV205_RESTART}"
    echo "[tpv205-rk45] RESTART: resuming from checkpoint prefix '${TPV205_RESTART}'"
    echo "              (this segment writes a fresh output dir: ${OUT})."
fi

# --- run -----------------------------------------------------------------
# OPTIMIZED CODE (ADER hot-path Lever 1 + Lever 3): --deriv-cache precomputes
# the per-element D_d^e = M^-1 K_d and the volume-RHS operator S_d^e, switching
# ComputeVolumeRHS / ApplySpatialDerivative to a dense mat-vec (the dominant
# ~2/3 of step time).  It IS exercised by the RK45 path: the coupled
# AdvanceRKCoupled*_Spatial steppers call wave.Mult -> ComputeVolumeRHS on every
# stage (rk_time_stepper.hpp).  Result changes at round-off only (REVIEW R-002);
# safe on this scalar/homogeneous TPV operator (spatial_dyn_driver.cpp:1146).
# Lever 2 (--shared-ck-recursion) is ADER-predictor-only (consumed solely by
# AdvanceADERWithSubStep_Spatial) and is DELIBERATELY OMITTED — it is a silent
# no-op on the RK45 path.
ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --time-integrator rk45 \
    --mixed-flux adjacent \
    --deriv-cache \
    --tfinal "${TFINAL}" \
    --paraview \
    --paraview-fault-vtu \
    --paraview-fault-dt 0.1 \
    --paraview-max-snapshots 300 \
    --print-derived \
    --checkpoint-every 2000 \
    --output-dir "${OUT}" \
    ${RESTART_FLAGS} \
    --verify-dispatch

echo "[tpv205-rk45] done. On-fault SCEC station traces: ${OUT}/tpv205_station_*.dat (16 stations); compare vs the SCEC TPV205"
echo "              reference.  Confirm V_max is BOUNDED (no ADER-central-flux"
echo "              run-away) and read the steps/hour off this *.out to size"
echo "              the restart chain (see RESTART CHAIN below)."

# =========================================================================
# RESTART CHAIN (full 12 s spans multiple 48 h normal-queue segments)
# -------------------------------------------------------------------------
# The run checkpoints every 2000 steps (--checkpoint-every 2000) into ${OUT}.
# To continue after a wall-limit truncation, re-submit pointing TPV205_RESTART
# at the PREVIOUS segment's cp prefix:
#       TPV205_RESTART="tpv205/out_p2_rk45_mf_normal_<prevjobid>/cp" \
#           sbatch tpv205_p2_rk45_mixedflux_12N_600r_normal.sbatch
# Each segment gets its OWN job-ID output dir (REQUIRED: the driver refuses a
# restart whose --output-dir equals the checkpoint's parent dir).
#
# Pre-build the mesh ONCE on a login/idev node (conda activate pythonenv;
# gmsh -format msh22 -3 ${GEO} -o ${MESH}) so the compute job does not mesh.
# Recalibrate dt from the empirical steps/hour (the RK-aware mixed-flux CFL
# factor 0.6 for Adjacent is a STARTING calibration).
# =========================================================================
