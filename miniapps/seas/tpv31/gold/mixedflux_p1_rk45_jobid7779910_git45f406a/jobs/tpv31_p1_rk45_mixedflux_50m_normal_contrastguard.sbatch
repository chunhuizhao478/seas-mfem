#!/bin/bash
#SBATCH -J tpv31_p1_rk45_mf_cg
#SBATCH -o tpv31_p1_rk45_mf_n_%j.out
#SBATCH -e tpv31_p1_rk45_mf_n_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# =========================================================================
# SCEC TPV31 — p1 · MATRIX (bi-material) Riemann · mixed (central) flux · RK45
# Driver: seas_spatial_dyn_driver (spatial dynamic-rupture).  Normal-queue,
# 10 nodes / 500 ranks (50 ranks/node on Frontera CLX 56-core nodes) / 48 h wall.
#
# *** FIRST job to combine matrix-bimaterial + central flux + RK45 — the exact
#     path PLAN_mixed_flux_hetero_riemann.md Phases 1-6 build (Phase 7 final). ***
#
# NUMERICAL SCHEME (the four knobs that define this job) ------------------
#   * p1            : [mesh].order = 1            (config-only, asserted below).
#                     DG CFL safety auto-scales 1/(3·(2p+1)) = 1/9 at p=1
#                     (cfl_safety="dg"; spatial_friction.hpp).
#   * matrix Riemann: [numerics].interior_flux = "matrix"  (MANDATORY for
#                     TPV31's depth_profile_1d bi-material; config-only, no CLI
#                     override — asserted below) -> BimaterialWaveOperator.
#   * mixed flux    : --mixed-flux adjacent       (the bi-material CENTRAL flux
#                     on fault-adjacent interior faces; bi-material Godunov
#                     upwind elsewhere — Phases 1-3).  Re-asserts the TOML.
#   * RK45          : --time-integrator rk45       (Dormand-Prince RK45; NOT
#                     ADER).  The central flux is non-dissipative (imaginary-axis
#                     modes) and UNSTABLE under ADER-O2; an explicit RK of order
#                     >= 3 has imaginary-axis coverage and is stable.  The driver
#                     G2 guard + WaveOperator::ComputeMaxDt REJECT matrix+mixed
#                     under ADER, so this scheme MUST be RK (Phases 4-5).  See
#                     document/mixed_flux_dev/
#                     BUILD_mixed_flux_rk_dynamic_rupture_2026-05-28.md.
#   Config: tpv31/configs/tpv31_rk_mixedflux.toml (variant of tpv31.toml: ONLY
#           mixed_flux none->adjacent; [material].seam_continuous=true is RETAINED
#           but is now a DEPRECATED no-op — see PARALLEL SEAM MATERIAL below).
#
#   *** CONTRAST-GUARD VARIANT (unified bi-material plan, Part A) ***
#   Adds --mixed-flux-contrast-tol 0.05: fault-adjacent corridor faces that cross a
#   STRONG impedance contrast (the depth_profile_1d material steps at 2400/5000/
#   10000 m) drop from the central corridor and dispatch the (dissipative) bi-
#   material upwind instead of the non-dissipative central flux.  This targets the
#   sigma_n collapse this session localized to exactly those 3 interface depths
#   (dp024/dp050/dp100: 16/9/64 MPa) while leaving every uniform-material corridor
#   face central.  Confirm engagement: grep the *.out for "CONTRAST GUARD" — the
#   "reclassified central -> upwind" count MUST be > 0 (else the flag was not
#   applied).  Default (no flag / tol<0) is byte-exact vs the gold mixedflux run.
#
#   The CLI re-asserts --time-integrator rk45 / --mixed-flux adjacent so the
#   resolved scheme prints in the startup banner + the --verify-dispatch block
#   (grep the *.out for "time integrator:", "mixed flux:", "interior flux") and
#   cannot silently drift from the config.  [mesh].order and interior_flux have
#   NO CLI override, so they are pre-flight-asserted against the TOML.
#
# CFL IS A STARTING CALIBRATION (read the *.out before committing hours) --
#   The matrix+RK+mixed-flux CFL is an EMPIRICAL starting calibration: the
#   RK-aware Adjacent factor 0.6 (Phase 4, wave_operator.inl §14.4) x cfl=0.5 x
#   the DG-order de-rating (cfl_safety="dg").  NO `--dt` is passed (project
#   memory: explicit-RK CFL — use the driver's auto-CFL; do NOT transfer ADER
#   dt values).  This is NOT a validated stability bound — read the empirical
#   stable dt off `--print-derived` + the *.out, and if V_max runs away or
#   dt->0, recalibrate `cfl` / the Phase-4 factor on Frontera (NOT a local gate).
#
# PARALLEL SEAM MATERIAL (cross-rank true-peer read) ---------------------
#   At np>1 ParMETIS cuts the depth_profile_1d bi-material across partition
#   seams.  The cross-rank material exchange (branch fix-cross-rank-material-
#   exchange) reads the TRUE peer element's material at every seam — for the
#   bi-material CENTRAL flux on fault-adjacent SHARED faces AND the bulk
#   bi-material Godunov UPWIND seam.  This REPLACES the old local-side
#   neighbour-material stub (which reused the LOCAL element's material as the
#   neighbour's — valid only for a depth-only profile on a seam-vertical
#   partition).  Consequences:
#     (1) [material].seam_continuous is now a DEPRECATED no-op (kept only so the
#         config still parses; the guard it cleared is GONE — the startup banner
#         prints "seam_continuous=...; DEPRECATED, no effect").
#     (2) The seam flux at the 2400/5000/10000 m depth discontinuities is
#         RE-BASELINED to the correct 2-sided peer material (== serial); it is
#         byte-exact ONLY for seam-vertical partitions, so the Frontera gold for
#         this job MUST be re-baselined against a fresh np>1 run on this branch.
#   This _contrastguard variant ALSO passes --mixed-flux-contrast-tol 0.05, which
#   drops the strong-contrast corridor faces to upwind (see the CONTRAST-GUARD
#   block above) — that is the strong-contrast safety, NOT seam_continuous.
#   Local guards: tests/parallel/test_bimaterial_seam_material_np2 (135/135) +
#   test_bimaterial_seam_fault_np2 (24/24).
#
# --deriv-cache OMITTED for this FIRST run (deliberate — NOT because unsafe) -
#   The deriv cache works for BOTH the scalar and bimaterial (matrix) operators
#   (geometry-only; round-off-only per R-002; budget-guarded per R-004,
#   ~18 MB/rank at p2 for TPV31) — spatial_dyn_driver.cpp:1140-1158 enables it
#   unconditionally with no bimaterial rejection.  It is omitted here only to
#   (a) isolate the CFL calibration with one fewer variable, and (b) because it
#   accelerates the ADER hot path while THIS job is RK45 (limited benefit).
#   Re-enable `--deriv-cache` once the CFL is calibrated.
#
# LOCAL VERIFICATION SCOPE -----------------------------------------------
#   The merged path is verified locally by unit tests (test_bimaterial_central_
#   flux, _dispatch [incl. Test 3.2 het single-valuedness], _shared np=2,
#   _cfl, _homog_equivalence) + the parse/dry tests (test_tpv_config_parse:
#   matrix+adjacent clears G1+G2; seam_continuous parses).  The driver's
#   END-TO-END matrix+RK+adjacent dispatch PAST mesh load is NOT exercised
#   locally (no-local-full-mesh rule).  The ibrun run below + --verify-dispatch
#   is the FIRST execution that loads the real TPV31 mesh, builds the operator,
#   and resolves the dispatch — it is the Frontera dispatch GATE.  Confirm the
#   banner/verify-dispatch grep (below) before the long run; a short idev slice
#   (small --tfinal) is the recommended first gate.
#
# Reference overlays: tpv31/benchmark_data/scec_{eqdyna,seisol}/.
# Spec: tpv31/benchmark_document/TPV31_32_Description_v03.pdf (15 s, p. 9).
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
    echo "       Build it on a login/idev node from THIS branch (carries the"
    echo "       merged matrix+central-flux+RK code):"
    echo "         cd miniapps/seas && make seas_spatial_dyn_driver"
    exit 1
fi

CONFIG_TOML="${TPV31_CONFIG:-tpv31/configs/tpv31_rk_mixedflux.toml}"
MESH="${TPV31_MESH:-tpv31/mesh/tpv31_50m.msh}"
GEO="${TPV31_GEO:-tpv31/mesh/tpv31_50m.geo}"

if [[ ! -f "${CONFIG_TOML}" ]]; then
    echo "ERROR: config ${CONFIG_TOML} not found (cwd $(pwd))."; exit 1
fi

# --- pre-flight: assert the scheme knobs against the TOML (fail-loud) -------
# interior_flux="matrix" is MANDATORY for TPV31 (bi-material depth profile) and
# is config-only (no CLI override) — a drifted TOML would silently change the
# scheme.  mixed_flux="adjacent" is the feature this job validates (the CLI
# re-asserts it, but assert the config too so a none/all_continuous typo is
# caught here).  [mesh].order=1 is the "p1" in this job's name (config-only).
if ! grep -Eq '^[[:space:]]*interior_flux[[:space:]]*=[[:space:]]*"matrix"' "${CONFIG_TOML}"; then
    echo "ERROR: ${CONFIG_TOML} does not set [numerics].interior_flux = \"matrix\"."
    echo "       TPV31 (depth_profile_1d bi-material) REQUIRES the matrix Riemann path."
    exit 4
fi
if ! grep -Eq '^[[:space:]]*mixed_flux[[:space:]]*=[[:space:]]*"adjacent"' "${CONFIG_TOML}"; then
    echo "ERROR: ${CONFIG_TOML} does not set [numerics].mixed_flux = \"adjacent\"."
    echo "       This job validates the bi-material CENTRAL flux; use"
    echo "       tpv31/configs/tpv31_rk_mixedflux.toml (not tpv31.toml)."
    exit 4
fi
if ! grep -Eq '^[[:space:]]*order[[:space:]]*=[[:space:]]*1([[:space:]]|$)' "${CONFIG_TOML}"; then
    echo "ERROR: ${CONFIG_TOML} does not set [mesh].order = 1 (this job is p1)."
    exit 4
fi

# --- regenerate the Gmsh v2.2 mesh if it is missing -----------------------
# The .msh is NOT committed (only the .geo).  MFEM's reader is Gmsh v2.2 only
# (see miniapps/seas/CLAUDE.md "Known limitation — Gmsh .msh format"), so the
# `-format msh22` flag is MANDATORY.  Pre-build it ONCE on a login/idev node
# (conda activate pythonenv) so the compute job does not spend wall time meshing.
if [[ ! -f "${MESH}" ]]; then
    echo "[tpv31-rk45-mf] mesh ${MESH} missing — generating from ${GEO} (Gmsh v2.2)..."
    if command -v gmsh >/dev/null 2>&1; then
        gmsh -format msh22 -3 "${GEO}" -o "${MESH}"
    else
        echo "ERROR: gmsh not on PATH.  Generate the mesh first (conda activate"
        echo "       pythonenv): gmsh -format msh22 -3 ${GEO} -o ${MESH}"
        exit 2
    fi
fi

OUT="${TPV31_OUT:-tpv31/out_p1_rk45_mf_normal_contrastguard_${SLURM_JOB_ID:-local}}"
TFINAL="${TPV31_TFINAL:-15}"   # spec p. 9 = 15 s (checkpoint/restart across the 48 h wall)

# --- restart across the 48 h wall (opt-in via TPV31_RESTART) ---------------
# The driver ABORTS (return 3) if --output-dir resolves to the SAME directory
# as the PARENT of --restart, so a restart segment MUST write to a DIFFERENT
# --output-dir than the checkpoint it reads.  Because OUT defaults to a
# per-job-ID directory, each segment is naturally distinct.  To resume:
#       TPV31_RESTART="<previous-segment-OUT>/cp" sbatch <thisscript>
# (point it at the cp prefix inside the PREVIOUS job's output dir).  Restart
# requires the SAME rank count that wrote the checkpoint (this job = -n 500).
RESTART_FLAGS=""
if [[ -n "${TPV31_RESTART:-}" ]]; then
    RESTART_FLAGS="--restart ${TPV31_RESTART}"
    echo "[tpv31-rk45-mf] RESTART: resuming from checkpoint prefix '${TPV31_RESTART}'"
    echo "                (this segment writes a fresh output dir: ${OUT})."
fi

# --- run -----------------------------------------------------------------
# Scheme pinned on the CLI (--time-integrator rk45 / --mixed-flux adjacent);
# matrix Riemann + p1 come from the (asserted) TOML.  NO --dt (explicit-RK
# auto-CFL).  NO --ader-order (RK path).  --deriv-cache OMITTED for this first
# run (see header).  --verify-dispatch prints the resolved scheme then CONTINUES
# into the full run (it does NOT exit like --dry-run).
ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --time-integrator rk45 \
    --mixed-flux adjacent \
    --mixed-flux-contrast-tol 0.05 \
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

echo "[tpv31-rk45-mf] done. On-fault station traces in ${OUT}; gate vs the SeisSol"
echo "                reference with the quantitative tolerance exit-code:"
echo "                  python3 tpv31/visualize_results.py --mfem ${OUT} --seisol \\"
echo "                      --tol-peak 0.10 --tol-rms 0.10"
echo "                Confirm V_max is BOUNDED (no central-flux run-away), NaN-free,"
echo "                and read the empirical stable dt off this *.out + --print-derived."

# =========================================================================
# DISPATCH GATE (Phase 7 — the end-to-end gate for BUG-3 + BUG-4)
# -------------------------------------------------------------------------
# Because of the no-local-full-mesh rule, the FIRST execution of the merged
# matrix+RK+adjacent dispatch on the real TPV31 mesh happens HERE on Frontera.
# After the run starts, grep the *.out with WHITESPACE-TOLERANT regex (the
# three scheme lines come from TWO sources — the rank-0 banner and the
# --verify-dispatch block — NOT brittle exact strings):
#       grep -E 'time integrator:'   <out>      # banner -> rk45
#       grep -E 'mixed flux:'        <out>      # banner -> adjacent
#       grep -E 'interior flux'      <out>      # --verify-dispatch -> matrix (...)
# Their presence proves the Phase-5 guard relaxation (G1/G2) + the Phase-2/3
# dispatch resolve on a production mesh.  Run a SHORT idev slice (small
# --tfinal) as the first gate before the 48 h run.
#
# RESTART CHAIN (full 15 s spans multiple 48 h normal-queue segments)
# -------------------------------------------------------------------------
# The run checkpoints every 2000 steps (--checkpoint-every 2000) into ${OUT}.
# To continue after a wall-limit truncation, re-submit pointing TPV31_RESTART
# at the PREVIOUS segment's cp prefix:
#       TPV31_RESTART="tpv31/out_p1_rk45_mf_normal_<prevjobid>/cp" \
#           sbatch tpv31_p1_rk45_mixedflux_50m_normal.sbatch
# Each segment gets its OWN job-ID output dir (REQUIRED: the driver refuses a
# restart whose --output-dir equals the checkpoint's parent dir).  Read the
# steps/hour off the *.out to size the number of segments.
# =========================================================================
