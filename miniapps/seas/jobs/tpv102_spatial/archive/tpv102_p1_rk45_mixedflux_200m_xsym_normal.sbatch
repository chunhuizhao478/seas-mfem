#!/bin/bash
#SBATCH -J tpv102_p1_rk45_mf_xsym
#SBATCH -o tpv102_p1_rk45_mf_xsym_%j.out
#SBATCH -e tpv102_p1_rk45_mf_xsym_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# =========================================================================
# SCEC TPV102 — p1 · SCALAR mixed (central) flux · explicit RK45
#   *** STRIKE-SYMMETRIC MESH arm of the dip-slip-drift fix A/B ***
# Driver: seas_spatial_dyn_driver (spatial dynamic-rupture).
# NORMAL QUEUE: 10 nodes / 500 ranks (-N 10 -n 500) / 48 h wall.
#
# WHY THIS JOB ------------------------------------------------------------
#   The on-fault DIP slip drifts (anti-symmetric in along-strike x, ~10
#   microrad spurious block rotation) and is NOT removed by order, flux
#   scheme, over-int/resample, friction, or the y-mirror "symmetric" hybrid
#   mesh.  Compiled symmetry probes prove every pointwise operator is exactly
#   x->-x equivariant, so the only strike asymmetry is the MESH:
#   tpv102_200m_hybrid.msh is y(fault-normal) symmetric but its
#   Frontal-Delaunay fault plane is strike-asymmetric (1.78% of fault nodes
#   have an -x mirror partner).  Root cause + fix:
#   debug_document/dip_drift_local_reproduction_2026-06-06.md.
#
#   THIS JOB runs the same p1 mixed-flux/RK45 scheme on the STRIKE-SYMMETRIC
#   mesh tpv102_200m_hybrid_xsym.msh (transfinite fault, mirrored diagonals:
#   fault triangulation 100% x-mirror verified by gmsh build + MFEM load;
#   y-mirror strip preserved, so sigma_n stays controlled).
#   ACCEPTANCE: the dip slip (slip1 / V_dip channels of tpv102_station_*.dat)
#   should collapse toward ~0 (DRDG3D reference) vs the y-mirror-mesh runs;
#   compare against tpv102_p1_rk45_mixedflux_10N_500r_normal (same scheme,
#   non-strike-symmetric mesh).  Strike slip, strike stress, sigma_n must be
#   unchanged.
#
# NUMERICAL SCHEME (the four knobs that define this job) ------------------
#   * p1            : [mesh].order = 1            (config-only, asserted below).
#   * mixed flux    : --mixed-flux adjacent       (central flux on the
#                     fault-adjacent faces; scalar Godunov elsewhere).
#   * RK45          : --time-integrator rk45       (Dormand-Prince RK45; NOT
#                     ADER).  At p1 the central flux is UNSTABLE under ADER-O2;
#                     an explicit RK of order >= 3 has imaginary-axis coverage
#                     and is stable.  See document/mixed_flux_dev/
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
#   * [meta].law = "rate_state" (Dieterich-Ruina, AGING law).
#   * [stress] kind = "fault_local_prestress": tau_strike = 75 MPa,
#     sigma_n = 120 MPa.  VW core a=0.008<b=0.012 inside 30x15 km, VS (a=0.016)
#     outside.  Dc = 0.02 m.
#   * Nucleation: compact-bell gradual overstress, delta_tau = 25 MPa, T_nuc = 1 s.
#   * Homogeneous material -> interior_flux = "scalar".
#
# CFL is a STARTING calibration (RK-aware Adjacent factor 0.6 x cfl 0.25);
# read the empirical stable dt off this job's *.out BEFORE committing
# normal-queue hours.  To promote to a production run see "PRODUCTION
# PROMOTION" below.
#
# LOCAL VERIFICATION SCOPE -----------------------------------------------
#   The RK code is verified locally by compile + unit tests.  The driver's
#   END-TO-END rk45 path PAST mesh load is NOT exercised locally (no-local-
#   full-mesh rule).  The ibrun run below carries --verify-dispatch, which
#   loads the real mesh + builds the operator + PRINTS the resolved scheme,
#   then CONTINUES into the full run.  Confirm the banner before the long run.
#   The strike-symmetric mesh was verified locally only at the GEOMETRY level
#   (gmsh build, MFEM load of a coarse build, fault 100% x-mirror); its
#   end-to-end driver load happens HERE.
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
# petsc/3.15 + mumps/5.3 prepend a PRE-2.31 libHYPRE to LD_LIBRARY_PATH after
# hypre/2.31.0 is loaded.  The driver was compiled against hypre/2.31.0 headers
# (HYPRE_Init -> HYPRE_Initialize), so the older lib shadowing it gives
# "undefined symbol: HYPRE_Initialize" + MPI exit 127.  Re-assert 2.31.0 at the
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
# STRIKE-SYMMETRIC 200 m hybrid mesh (the dip-drift fix).  Built from the .geo;
# pre-generate ONCE on a login node so the compute job does not mesh.
MESH="${TPV102_MESH:-tpv102/mesh/tpv102_200m_hybrid_xsym.msh}"
GEO="${TPV102_GEO:-tpv102/mesh/tpv102_200m_hybrid_xsym.geo}"

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
# NOTE: the strike-symmetric .geo writes NO $Periodic section (it uses a
# transfinite mirrored-diagonal fault, not a Periodic constraint), so the
# brittle v2.2 reader loads it cleanly.
if [[ ! -f "${MESH}" ]]; then
    echo "[tpv102-rk45-xsym] mesh ${MESH} missing — generating from ${GEO} (Gmsh v2.2)..."
    if command -v gmsh >/dev/null 2>&1; then
        gmsh -format msh22 -3 "${GEO}" -o "${MESH}"
    else
        echo "ERROR: gmsh not on PATH.  Generate the mesh first (conda activate"
        echo "       pythonenv): gmsh -format msh22 -3 ${GEO} -o ${MESH}"
        exit 2
    fi
fi

OUT="${TPV102_OUT:-tpv102/out_rk45_mf_xsym_${SLURM_JOB_ID:-local}}"
TFINAL="${TPV102_TFINAL:-12}"   # SCEC TPV102 window 0-12 s; the 48 h normal wall
                                # reaches a fraction — restart for the rest.

# --- restart from a checkpoint (opt-in via TPV102_RESTART) -----------------
# The run checkpoints every 2000 steps to ${OUT}/cp_checkpoint_r<rank>.txt.
# To resume, point TPV102_RESTART at the PREVIOUS segment's "cp" prefix; resume
# requires the same rank count (-n 500) and a DIFFERENT --output-dir.
RESTART_FLAGS=""
if [[ -n "${TPV102_RESTART:-}" ]]; then
    RESTART_FLAGS="--restart ${TPV102_RESTART}"
    echo "[tpv102-rk45-xsym] RESTART: resuming from checkpoint prefix '${TPV102_RESTART}'"
    echo "              (this segment writes a fresh output dir: ${OUT})."
fi

# --- run -----------------------------------------------------------------
# --deriv-cache precomputes per-element D_d^e = M^-1 K_d + the volume-RHS
# operator (the dominant ~2/3 of step time); exercised by the RK45 path
# (AdvanceRKCoupled*_Spatial -> wave.Mult).  Round-off-only change (REVIEW
# R-002); safe on this scalar/homogeneous TPV operator.
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

echo "[tpv102-rk45-xsym] done. On-fault SCEC station traces: ${OUT}/tpv102_station_*.dat (9 stations)."
echo "    DIP-DRIFT ACCEPTANCE: the dip channel (slip1 / V_dip) should be ~0 (DRDG3D"
echo "    reference) here, vs the secular drift on the y-mirror hybrid mesh.  Confirm"
echo "    strike slip / strike stress / sigma_n are unchanged, and V_max is BOUNDED."

# =========================================================================
# PRODUCTION PROMOTION (12 s to completion)
# -------------------------------------------------------------------------
#   1. The 48 h normal wall is the queue cap; to finish 12 s, chain restarts (2).
#   2. Chain restarts across the wall: checkpoints every 2000 steps to ${OUT};
#      resume with  --restart "${OUT}/cp".
#   3. Pre-build the mesh ONCE on a login/idev node (conda activate pythonenv;
#      gmsh -format msh22 -3 ${GEO} -o ${MESH}) so the compute job does not mesh.
#   4. Recalibrate dt from the dev smoke's empirical steps/hour.
# =========================================================================
