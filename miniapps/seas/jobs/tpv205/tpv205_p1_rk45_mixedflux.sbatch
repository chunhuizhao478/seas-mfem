#!/bin/bash
#SBATCH -J tpv205_p1_rk45_mf
#SBATCH -o tpv205_p1_rk45_mf_%j.out
#SBATCH -e tpv205_p1_rk45_mf_%j.err
#SBATCH -p development
#SBATCH -N 4
#SBATCH -n 224
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# =========================================================================
# SCEC TPV205 — p1 · SCALAR mixed (central) flux · explicit RK45
# Driver: seas_spatial_dyn_driver (spatial dynamic-rupture).
#
# NUMERICAL SCHEME (the four knobs that define this job) ------------------
#   * p1            : [mesh].order = 1            (linear elements; config-only,
#                     asserted below).
#   * mixed flux    : --mixed-flux adjacent       (central flux on the
#                     fault-adjacent faces; scalar Godunov elsewhere).
#   * RK45          : --time-integrator rk45       (Dormand-Prince RK45; NOT
#                     ADER).  At p1 the central flux is UNSTABLE under ADER-O2
#                     (the predictor is locked to O2 and amplifies the
#                     imaginary-axis modes central flux creates); an explicit
#                     RK of order >= 3 has imaginary-axis coverage and is
#                     stable.  See document/mixed_flux_dev/
#                     BUILD_mixed_flux_rk_dynamic_rupture_2026-05-28.md.
#   * scalar Riemann: [numerics].interior_flux = "scalar"  (config-only,
#                     asserted below).  MANDATORY for the RK path (mixed/central
#                     flux is scalar-only; matrix/bimaterial is rejected).
#
#   The CLI flags re-assert what tpv205_spatial_rk45_mixedflux.toml already
#   sets, so the resolved scheme is printed in the driver's startup banner
#   (grep the *.out for "time integrator:", "mixed flux:", "law:") and cannot
#   silently drift from the config.  [mesh].order and [numerics].interior_flux
#   have NO CLI override, so they are pre-flight-asserted against the TOML.
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
# This header job is the DEV SMOKE (development queue, 4N/224r/2 h).  CFL is a
# STARTING calibration (RK-aware Adjacent factor 0.6 x cfl 0.25); read the
# empirical stable dt off this job's *.out BEFORE committing normal-queue
# hours.  To promote to a production run see "PRODUCTION PROMOTION" below.
#
# LOCAL VERIFICATION SCOPE (R-002) ---------------------------------------
#   The RK+LSW *code* is verified locally by compile + unit tests
#   (seas_test_lsw_rk_mixed_flux L1-L7 + seas_test_lsw_rk_shared_fault_mpi L8)
#   and the scheme banner from `--dry-run --verify-dispatch` (which prints
#   law/integrator/flux BEFORE mesh load).  The driver's END-TO-END rk45+LSW
#   path PAST mesh load (wave-operator assembly, iterator construction, the
#   coupled-RK time loop) is NOT exercised locally — per the no-local-full-mesh
#   rule, and because the dry-run aborts at mesh load when the .msh is absent.
#   The `--dry-run --verify-dispatch` step BELOW (which passes --mesh) is the
#   first execution that loads a real mesh + builds the operator + resolves the
#   dispatch; it runs HERE on Frontera/idev as the gate.  Before a long
#   production run, confirm this dry-run exits 0 with the
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

# --- locate the repo root (walk up to the dir containing miniapps/seas) ----
_walk_up_dir="$(dirname "$(readlink -f "${BASH_SOURCE[0]:-$0}")")"
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

CONFIG_TOML="${TPV205_CONFIG:-tpv205/configs/tpv205_spatial_rk45_mixedflux.toml}"
MESH="${TPV205_MESH:-tpv205/mesh/tpv2053d_200m.msh}"

if [[ ! -f "${CONFIG_TOML}" ]]; then
    echo "ERROR: config ${CONFIG_TOML} not found (cwd $(pwd))."; exit 1
fi

# --- pre-flight: assert the two CONFIG-ONLY knobs (no CLI override) --------
# [mesh].order=1 is the "p1" in this job's name and interior_flux="scalar" is
# MANDATORY for the RK/mixed-flux path.  Both are config-only, so a drifted
# TOML would silently change the scheme.  Fail loudly here.
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
# The .msh is NOT committed (only the .geo).  MFEM's reader is Gmsh v2.2 only
# (see miniapps/seas/CLAUDE.md "Known limitation — Gmsh .msh format"), so the
# `-format msh22` flag is MANDATORY.  Run on a login/idev node (gmsh needs the
# `pythonenv` conda env) if gmsh is unavailable on the compute nodes.
if [[ ! -f "${MESH}" ]]; then
    echo "[tpv205-rk45] mesh ${MESH} missing — generating from the .geo (Gmsh v2.2)..."
    if command -v gmsh >/dev/null 2>&1; then
        gmsh -format msh22 -3 tpv205/mesh/tpv2053d_200m.geo -o "${MESH}"
    else
        echo "ERROR: gmsh not on PATH.  Generate the mesh first (conda activate"
        echo "       pythonenv): gmsh -format msh22 -3 tpv205/mesh/tpv2053d_200m.geo"
        echo "       -o ${MESH}"
        exit 2
    fi
fi

# --- dispatch dry-run (verify scalar + slip_weakening + RK45 + mixed-flux) -
echo "[tpv205-rk45] --dry-run --verify-dispatch (scheme: p1 / adjacent / rk45 / scalar):"
./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --time-integrator rk45 \
    --mixed-flux adjacent \
    --dry-run --verify-dispatch || {
    echo "ERROR: TPV205 RK45/mixed-flux config failed the dispatch dry-run."; exit 3; }

OUT="${TPV205_OUT:-tpv205/out_rk45_mf_${SLURM_JOB_ID:-local}}"
TFINAL="${TPV205_TFINAL:-12}"   # SCEC TPV205 request = 0-12 s (config tfinal is
                                # a short smoke; --tfinal overrides for the run).

# --- run -----------------------------------------------------------------
# Scheme pinned on the CLI (--time-integrator rk45 / --mixed-flux adjacent);
# scalar Riemann + p1 come from the (asserted) TOML.
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
    --output-dir "${OUT}"

echo "[tpv205-rk45] done. On-fault traces in ${OUT}; compare vs the SCEC TPV205"
echo "              reference.  Confirm V_max is BOUNDED (the ADER-p1-central"
echo "              flux run-away does NOT occur) and read the empirical stable"
echo "              dt off this *.out before committing normal-queue hours."

# =========================================================================
# PRODUCTION PROMOTION (12 s to completion)
# -------------------------------------------------------------------------
# The header above is a development-queue smoke.  To run to completion:
#   1. Swap the queue + wall + node count:
#        #SBATCH -p normal
#        #SBATCH -N 16          # (scale ranks: -n = 56 * N on Frontera CLX)
#        #SBATCH -n 896
#        #SBATCH -t 24:00:00
#   2. Chain restarts across the wall limit: the run writes checkpoints every
#      2000 steps (--checkpoint-every 2000) to ${OUT}.  Resume a truncated run
#      with:  --restart "${OUT}/cp".
#   3. Pre-build the mesh ONCE on a login/idev node (conda activate pythonenv;
#      gmsh -format msh22 -3 tpv205/mesh/tpv2053d_200m.geo
#      -o tpv205/mesh/tpv2053d_200m.msh) so the compute job does not mesh.
#   4. Recalibrate dt: the RK-aware mixed-flux CFL factors (0.6 for Adjacent)
#      are interim placeholders; size the normal-queue wall from the dev
#      smoke's empirical steps/hour.
# =========================================================================
