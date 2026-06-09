#!/bin/bash
#SBATCH -J tpv102_p1_rk45_pu
#SBATCH -o tpv102_p1_rk45_pu_%j.out
#SBATCH -e tpv102_p1_rk45_pu_%j.err
#SBATCH -p flex
#SBATCH -N 16
#SBATCH -n 800
#SBATCH -t 01:00:00
#SBATCH -A EAR20006

# =========================================================================
# SCEC TPV102 — p1 · PURE UPWIND (scalar Godunov) · explicit RK45
# Driver: seas_spatial_dyn_driver (spatial dynamic-rupture).
# FLEX SMOKE: development queue, 16N/800r/1 h.
#
# ▶▶ PURPOSE — σ_n ISOLATION EXPERIMENT (the §9 discriminator) ◀◀
#   This job exists to answer ONE question about the on-fault normal-stress
#   (σ_n) drift seen in the pure-upwind ADER runs:
#
#     "Is the σ_n drift caused by the ADER per-sub-step POINTWISE friction
#      application (the element-local Cauchy-Kovalevskaya predictor feeding a
#      noisy [[v_n]] into σ_n_trial = η_p·[[v_n]]+…), or by the BULK upwind
#      dissipation across the (asymmetric) fault-adjacent mesh (Zhang-2023 SSO)?"
#
#   It is the SINGLE-VARIABLE control for the existing TPV102 runs:
#     * job 7761452  PURE UPWIND + ADER-O2   → σ_n COLLAPSES
#                    (120 → 80.7 / −30.3 / 75.2 MPa at flt_0_7.5 / flt_0_3 / flt_9_7.5)
#     * job 7760376  MIXED FLUX  + RK45       → σ_n FLAT (120 ± 0.3 MPa)  [out_rk45_mf]
#     * THIS job     PURE UPWIND + RK45       → ??? (the missing corner)
#
#   Relative to job 7761452 (PURE UPWIND + ADER), the ONLY thing changed here
#   is the time integrator (ADER → RK45).  Flux (pure upwind), interior_flux
#   (scalar), friction law (RS aging), mesh (200 m), stress, nucleation, CFL,
#   and the same production config tpv102/configs/tpv102_spatial.toml are ALL
#   held identical.  So this isolates ADER-vs-RK45 with the flux held at pure
#   upwind — i.e. it removes the bulk Zhang-SSO confound that job 7760376 (which
#   changed BOTH flux AND integrator) could not.
#
#   DECISIVE READ (compare σ_n col 8 of ${OUT}/tpv102_station_*.dat vs 7761452):
#     * σ_n stays FLAT at ~120 MPa  ⇒ the ADER per-sub-step pointwise predictor
#       application is the dominant seed.  PORT SeisSol's single time-integrated
#       imposed-state / common-normal-stress application (this ALSO helps TPV31,
#       which is forced onto pure-upwind+ADER+matrix and can use neither mixed
#       flux nor — today — RK45).
#     * σ_n STILL DRIFTS like 7761452 ⇒ the BULK upwind dissipation on the
#       asymmetric fault-adjacent mesh (Zhang-2023 SSO) is the dominant seed.
#       Then only mixed-flux-on-the-matrix-path or symmetric fault-adjacent
#       meshing can save TPV31; the integrator is not the lever.
#   (Full diagnosis: debug_document/tpv31_debug_document/.)
#
# NUMERICAL SCHEME (the four knobs that define this job) ------------------
#   * p1            : [mesh].order = 1            (config-only, asserted below).
#   * pure upwind   : --mixed-flux none           (upwind Godunov on EVERY
#                     interior face; NO central-flux faces — identical flux to
#                     job 7761452).  MixedFluxMode::None ⇒ CFL factor 1.0
#                     (wave_operator.inl:5883), so RK45 at cfl 0.25 is
#                     conservative/stable (pure upwind is dissipative; DP-RK45
#                     has ample real-axis coverage).
#   * RK45          : --time-integrator rk45       (Dormand-Prince RK45; NOT
#                     ADER).  Supported with rate_state on the scalar path
#                     (spatial_dyn_driver.cpp:784-816; the ONLY RK precondition
#                     is interior_flux="scalar", which TPV102 satisfies — there
#                     is NO requirement that mixed_flux != none).
#   * scalar Riemann: [numerics].interior_flux = "scalar"  (config-only,
#                     asserted below).  TPV102 is a homogeneous half-space.
#
#   The CLI flags (--time-integrator rk45 / --mixed-flux none) OVERRIDE the
#   production config's defaults (which are ADER + mixed_flux="adjacent"); the
#   resolved scheme is printed in the driver's startup banner (grep the *.out
#   for "time integrator:", "mixed flux:", "law:") and cannot silently drift.
#   [mesh].order and [numerics].interior_flux have NO CLI override, so they are
#   pre-flight-asserted against the TOML.
#
# TPV102 physics (IDENTICAL to job 7761452 — same config) -----------------
#   * [meta].law = "rate_state" (Dieterich-Ruina, AGING law).  The RK path runs
#     the coupled RK-on-(Q, psi, slip) stepper AdvanceRKCoupled_Spatial, whose
#     wave.Mult fault face uses the instantaneous FaultFaceFlux::Evaluate (the
#     dt→0 limit of the RS ADER Riemann solve).  [numerics].fault_iterator =
#     "substep" is IGNORED on the RK path (ADER-only; kept for parity).
#   * [stress] kind = "fault_local_prestress": tau_strike = 75 MPa,
#     sigma_n = 120 MPa.  VW core a=0.008<b=0.012 inside 30x15 km, VS (a=0.016)
#     outside.  Dc = 0.02 m.
#   * Nucleation: compact-bell gradual overstress, delta_tau = 25 MPa, T_nuc = 1 s.
#   * Homogeneous material -> interior_flux = "scalar".
#
#   RS-AT-np>1 NOTE: as with 7761452, shared rank-seam fault QPs are treated as
#   preliminary at np>1.  The σ_n onset in 7761452 begins ~1.5-2 s, so even a
#   partial (few-second) run here is enough to read the discriminator.
#
# CFL is the production cfl=0.25 (None factor 1.0); read the empirical stable
# dt and steps/hour off this job's *.out.  If the run does not reach ~2-3 s in
# 1 h, resume via --restart "${OUT}/cp" (checkpoints every 2000 steps).
#
# LOCAL VERIFICATION SCOPE -----------------------------------------------
#   The RK code is verified locally by compile + unit tests
#   (seas_test_rk_time_stepper, seas_test_fault_face_flux_ader_equivalence,
#   seas_test_wave_operator).  The driver's END-TO-END rk45 path PAST mesh load
#   is exercised on Frontera; --verify-dispatch below PRINTS the resolved scheme
#   and CONTINUES into the full run.  Confirm the banner shows
#   "time integrator: rk45" and "mixed flux: none" in the *.out.
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

# Use the SAME production config as job 7761452 (the PURE-UPWIND + ADER run)
# so the ONLY difference is the integrator (overridden to rk45 on the CLI).
CONFIG_TOML="${TPV102_CONFIG:-tpv102/configs/tpv102_spatial.toml}"
# 200 m mesh (real work for 800 ranks; the TOML default is the 1000 m mesh).
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
    echo "       The RK path is scalar-only (matrix/bimaterial rejected)."
    exit 4
fi

# --- regenerate the Gmsh v2.2 mesh if it is missing -----------------------
# MFEM's reader is Gmsh v2.2 only (see miniapps/seas/CLAUDE.md "Known
# limitation — Gmsh .msh format"); the `-format msh22` flag is MANDATORY.
if [[ ! -f "${MESH}" ]]; then
    echo "[tpv102-rk45-pu] mesh ${MESH} missing — generating from ${GEO} (Gmsh v2.2)..."
    if command -v gmsh >/dev/null 2>&1; then
        gmsh -format msh22 -3 "${GEO}" -o "${MESH}"
    else
        echo "ERROR: gmsh not on PATH.  Generate the mesh first (conda activate"
        echo "       pythonenv): gmsh -format msh22 -3 ${GEO} -o ${MESH}"
        exit 2
    fi
fi

OUT="${TPV102_OUT:-tpv102/out_p1_rk45_pureupwind_${SLURM_JOB_ID:-local}}"
TFINAL="${TPV102_TFINAL:-12}"   # SCEC TPV102 window 0-12 s; the 1 h flex wall
                                # reaches a fraction — restart for the rest.

# --- restart from a checkpoint (opt-in via TPV102_RESTART) -----------------
# The run checkpoints every 2000 steps to ${OUT}/cp_checkpoint_r<rank>.txt.
# To resume, point TPV102_RESTART at the PREVIOUS segment's "cp" prefix; the
# driver appends _checkpoint_r<rank>.txt per rank.  Two requirements:
#   * same rank count that wrote the checkpoint (this job = -n 800), and
#   * a DIFFERENT --output-dir than the checkpoint's dir (the driver aborts,
#     return 3, if --output-dir == the parent of --restart).  OUT defaults to
#     a per-job-ID dir, so a fresh submit is naturally distinct.
RESTART_FLAGS=""
if [[ -n "${TPV102_RESTART:-}" ]]; then
    RESTART_FLAGS="--restart ${TPV102_RESTART}"
    echo "[tpv102-rk45-pu] RESTART: resuming from checkpoint prefix '${TPV102_RESTART}'"
    echo "                 (this segment writes a fresh output dir: ${OUT})."
fi

# --- run -----------------------------------------------------------------
# OPTIMIZED CODE (ADER hot-path Lever 1 + Lever 3): --deriv-cache precomputes
# the per-element D_d^e = M^-1 K_d and the volume-RHS operator S_d^e.  It IS
# exercised by the RK45 path (AdvanceRKCoupled_Spatial calls wave.Mult ->
# ComputeVolumeRHS on every stage); result changes at round-off only.
# Lever 2 (--shared-ck-recursion) is ADER-predictor-only and is DELIBERATELY
# OMITTED here (silent no-op on the RK45 path; matches the rk45_mf dev job).
#
# Flux + integrator pinned on the CLI (--time-integrator rk45 / --mixed-flux
# none) — these OVERRIDE the production config's ADER + adjacent defaults so the
# resolved scheme is printed in the banner.  scalar Riemann + p1 come from the
# (asserted) TOML.
ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --time-integrator rk45 \
    --mixed-flux none \
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

echo "[tpv102-rk45-pu] done. On-fault SCEC station traces: ${OUT}/tpv102_station_*.dat (9 stations)."
echo "  σ_n DISCRIMINATOR (col 8): compare against"
echo "    7761452 (PU+ADER, σ_n COLLAPSES 120->80/-30/75 MPa) and"
echo "    7760376 (MF+RK45, σ_n FLAT 120 +-0.3 MPa)."
echo "  σ_n FLAT here  => ADER per-sub-step pointwise application is the seed (port SeisSol-style"
echo "                    time-integrated single-imposed-state application; helps TPV31)."
echo "  σ_n DRIFTS here => bulk Zhang-2023 SSO (asymmetric fault-adjacent upwind) is the seed"
echo "                    (only mixed-flux-on-matrix or symmetric meshing saves TPV31)."

# =========================================================================
# RESTART (resume across the 1 h flex wall)
# -------------------------------------------------------------------------
#   The run checkpoints every 2000 steps to ${OUT}/cp.  The TPV102_RESTART hook
#   above is now WIRED into the ibrun, so resuming is just (each submit gets a
#   fresh per-job-ID OUT, so --output-dir != the checkpoint's parent dir
#   automatically — the driver ABORTS return 3 if they coincide):
#     TPV102_RESTART="tpv102/out_p1_rk45_pureupwind_<PREVJOBID>/cp" \
#       sbatch tpv102_p1_rk45_pureupwind_200m_flex.sbatch
#   Keep -n 800 (restart must match the rank count that wrote the checkpoint).
# =========================================================================
