#!/bin/bash
#SBATCH -J tpv104_p1_aderO2_pu
#SBATCH -o tpv104_p1_aderO2_pu_%j.out
#SBATCH -e tpv104_p1_aderO2_pu_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# =========================================================================
# SCEC TPV104 — JOB [1] (PRODUCTION):  p1 · ADER-O2 · pure upwind · scalar Riemann
# Driver: seas_spatial_dyn_driver (spatial dynamic-rupture).  Normal-queue,
# full-12 s sibling of the development/flex job
# (tpv104_p1_rk45_mixedflux_8N_400r_dev.sbatch), but with the ADER + PURE-UPWIND
# scheme (NOT the RK45 + central/mixed flux of the dev job).
#
# NUMERICAL SCHEME (the four knobs that define this job) ------------------
#   * p1            : [mesh].order = 1            (linear FE space on the
#                     straight-tet 200 m mesh; config-only, no CLI override —
#                     asserted below).  DG CFL safety auto-scales:
#                     1/(3·(2p+1)) = 1/9 at p=1 (spatial_friction.hpp:619).
#   * ADER-O2       : --ader-order 2 --time-integrator ader  (2nd-order ADER-DG
#                     sub-step; allowed range {2,3,4}).
#   * pure upwind   : --mixed-flux none           (upwind Godunov on EVERY
#                     interior face; no central-flux faces).  Pure upwind is
#                     stable under ADER (the central/mixed flux is what forced
#                     the dev job onto RK45); CFL factor 1.0.
#   * scalar Riemann: [numerics].interior_flux = "scalar"  (TPV104 is a
#                     HOMOGENEOUS half-space → scalar Godunov; config-only, no
#                     CLI override — asserted below).  This is the key contrast
#                     with TPV31, which is bi-material and REQUIRES "matrix".
#   Config: tpv104/configs/tpv104_spatial.toml (p1 / ADER-O2 baseline).
#
# TPV104 physics ---------------------------------------------------------
#   * [meta].law = "rate_state" (Dieterich-Ruina, SLIP law with strong rate
#     weakening / FVW).  RS uses the ADER sub-step fault iterator
#     ([numerics].fault_iterator = "substep").
#   * [stress] kind = "fault_local_prestress": tau_strike = 40 MPa,
#     sigma_n = 120 MPa.  VW core a=0.01<b=0.014, V_w=0.1 inside 30x15 km;
#     VS (a=0.02, V_w=1.0) outside.  Dc = 0.4 m, f_w = 0.2.
#   * Nucleation: compact-bell gradual overstress, delta_tau = 45 MPa, T_nuc = 1 s.
#   * Homogeneous material → interior_flux = "scalar".
#
# RS-AT-np>1 GATE: the driver prints "Treat parallel RS results as PRELIMINARY"
#   at np>1 (shared rank-seam fault QPs use end-of-step psi, 1st-order at
#   ader_order>=2).  Treat this production run accordingly.
#
# QUEUE / SCALE -----------------------------------------------------------
#   normal queue, 10 nodes, 500 MPI ranks (50 ranks/node on Frontera CLX
#   56-core nodes).  Wall = 48 h (normal-queue maximum).  The 0-12 s SCEC
#   window may span more than one 48 h segment on the 200 m mesh — rely on
#   checkpoint/restart (see RESTART below).
#
# ADER HOT-PATH OPTIMIZATION (ALWAYS ON here) -----------------------------
#   This production job ALWAYS passes the Lever 1/2/3 optimization flags
#   (--deriv-cache --shared-ck-recursion); see "OPTIMIZATION" below.  The
#   --deriv-cache cache works for the SCALAR operator too (driver comment:
#   "Works for both the scalar and bimaterial (matrix) operators").  Result
#   differs from the bit-exact default at ROUND-OFF only (machine-eps).
#   ⇒ The binary MUST be (re)built from this branch (it carries the merged
#     elem_derivative_cache.hpp / DerivMode::Cached code).  A binary built
#     before the merge IGNORES these flags SILENTLY (HasFlag argv-scan) and
#     runs the slow path — grep the *.out for "[deriv] ... enabled" to confirm.
#
# PERFGRAPH (Caliper; ON by default here) ---------------------------------
#   TPV104_PERFGRAPH defaults to 1 so the Caliper region tree is collected.
#   This requires an MFEM_USE_CALIPER=YES build; on a non-Caliper binary the
#   annotations are no-ops and CALI_CONFIG is ignored (the run still proceeds).
#
# Reference overlays: tpv104 benchmark + drivers/tpv104_driver.cpp.
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
    echo "       ADER optimization).  For perfgraph + optimization, do a CLEAN"
    echo "       Caliper build (the WaveOperator templates changed):"
    echo "         cd miniapps/seas && make clean && \\"
    echo "           USE_CALIPER=YES CALIPER_DIR=<prefix> bash build_frontera.sh"
    exit 1
fi

CONFIG_TOML="${TPV104_CONFIG:-tpv104/configs/tpv104_spatial.toml}"
# 200 m mesh (production resolution; the TOML default is the coarse 1000 m mesh).
MESH="${TPV104_MESH:-tpv104/mesh/tpv104_200m.msh}"
GEO="${TPV104_GEO:-tpv104/mesh/tpv104_200m.geo}"

if [[ ! -f "${CONFIG_TOML}" ]]; then
    echo "ERROR: config ${CONFIG_TOML} not found (cwd $(pwd))."; exit 1
fi

# --- pre-flight: assert the two CONFIG-ONLY knobs (no CLI override) --------
# [mesh].order=1 is the "p1" in this job's name; interior_flux="scalar" is
# MANDATORY for homogeneous TPV104 (and for the pure-upwind scalar Godunov).
# Both are config-only, so a drifted TOML would silently change the scheme.
if ! grep -Eq '^[[:space:]]*order[[:space:]]*=[[:space:]]*1([[:space:]]|$)' "${CONFIG_TOML}"; then
    echo "ERROR: ${CONFIG_TOML} does not set [mesh].order = 1 (this job is p1)."
    echo "       Use tpv104/configs/tpv104_spatial.toml (not _p2/_p3)."
    exit 4
fi
if ! grep -Eq '^[[:space:]]*interior_flux[[:space:]]*=[[:space:]]*"scalar"' "${CONFIG_TOML}"; then
    echo "ERROR: ${CONFIG_TOML} does not set [numerics].interior_flux = \"scalar\"."
    echo "       TPV104 is homogeneous; the scalar Godunov path is REQUIRED."
    exit 4
fi

# --- regenerate the Gmsh v2.2 mesh if it is missing -----------------------
# MFEM's reader is Gmsh v2.2 only (see miniapps/seas/CLAUDE.md "Known
# limitation — Gmsh .msh format"); the `-format msh22` flag is MANDATORY.
# Pre-build it ONCE on a login/idev node (conda activate pythonenv) so the
# compute job does not spend wall time meshing.
if [[ ! -f "${MESH}" ]]; then
    echo "[tpv104-p1] mesh ${MESH} missing — generating from ${GEO} (Gmsh v2.2)..."
    if command -v gmsh >/dev/null 2>&1; then
        gmsh -format msh22 -3 "${GEO}" -o "${MESH}"
    else
        echo "ERROR: gmsh not on PATH.  Generate the mesh first (conda activate"
        echo "       pythonenv): gmsh -format msh22 -3 ${GEO} -o ${MESH}"
        exit 2
    fi
fi

# --- positional arg: `perfgraph` (Caliper short slice that flushes) --------
# Submit with the single word `perfgraph` to run a SHORT slice (default
# ~0.004 s) that REACHES a clean exit so the Caliper runtime-report actually
# flushes:   sbatch tpv104_p1_aderO2_pureupwind_200m_normal.sbatch perfgraph
# A wall-clock-killed run (CANCELLED ... DUE TO TIME LIMIT) emits NO report,
# so the full 12 s run's Caliper data only covers whichever restart segment
# reaches t=tfinal cleanly — use this short slice for a representative tree.
# Slice-length override: TPV104_PERFGRAPH_TFINAL (s).
PERFGRAPH_MODE=0
for _arg in "$@"; do
    case "${_arg}" in
        perfgraph) PERFGRAPH_MODE=1 ;;
        *) echo "[tpv104-p1] WARNING: ignoring unknown arg '${_arg}' (expected 'perfgraph')." ;;
    esac
done

OUT="${TPV104_OUT:-tpv104/out_p1_aderO2_pu_normal_${SLURM_JOB_ID:-local}}"

if [[ "${PERFGRAPH_MODE}" == "1" ]]; then
    TFINAL="${TPV104_TFINAL:-${TPV104_PERFGRAPH_TFINAL:-0.004}}"   # short -> finishes -> report flushes
    echo "[tpv104-p1] PERFGRAPH MODE: profiling slice tfinal=${TFINAL} s."
    echo "            Runs to a CLEAN exit so the Caliper report flushes; the full run is tfinal=12 s."
else
    TFINAL="${TPV104_TFINAL:-12}"   # SCEC TPV104 window 0-12 s (checkpoint/restart across the 48 h wall)
fi

# --- ADER hot-path optimization (ALWAYS ON in this production job) ---------
# --deriv-cache         : precompute per-element D_d^e = M^-1 K_d (Lever 1) and the volume-RHS operator S_d^e (Lever 3);
#                         apply as a dense mat-vec instead of per-QP quadrature
#                         (DerivMode::Cached).  Works for the SCALAR operator
#                         (driver: "Works for both the scalar and bimaterial
#                         (matrix) operators").  Changes results at ROUND-OFF
#                         only (machine-eps).  Cache budget (R-004) at p1 is tiny.
# --shared-ck-recursion : run the Cauchy-Kovalevskaya predictor recursion ONCE
#                         per macro-step instead of twice (Lever 2).
#                         BYTE-IDENTICAL at a fixed DerivMode.
# The driver prints "[deriv] ... enabled" on rank 0 when the flags take effect
# — grep the *.out to confirm the binary actually carries the optimization.
DERIV_OPT_FLAGS="--deriv-cache --shared-ck-recursion"

# --- restart across the 48 h wall (opt-in via TPV104_RESTART) --------------
# The driver ABORTS (return 3) if --output-dir resolves to the SAME directory
# as the PARENT of --restart, so a restart segment MUST write to a DIFFERENT
# --output-dir than the checkpoint it reads.  Because OUT defaults to a
# per-job-ID directory, each segment is naturally distinct.  To resume:
#       TPV104_RESTART="<previous-segment-OUT>/cp" sbatch <thisscript>
RESTART_FLAGS=""
if [[ -n "${TPV104_RESTART:-}" ]]; then
    RESTART_FLAGS="--restart ${TPV104_RESTART}"
    echo "[tpv104-p1] RESTART: resuming from checkpoint prefix '${TPV104_RESTART}'"
    echo "            (this segment writes a fresh output dir: ${OUT})."
fi

# --- Caliper "perfgraph" (MFEM profiling facility); ON by default ----------
# MFEM's perfgraph is Caliper: MFEM_PERF_SCOPE annotations (seas:: regions in
# dynamic/wave_operator.inl and drivers/spatial_dyn_driver.cpp) are no-ops
# UNLESS MFEM is built with MFEM_USE_CALIPER=YES + libcaliper.  A region-timing
# TREE is then produced at runtime via CALI_CONFIG.  Recipe:
# document/caliper_perfgraph_dev/README_caliper_build_and_run.md.  This is an
# ADER job, so the ADER-path regions fire (ComputeADER{SubStepStates,
# VolumeUpdate,FaceFluxRHS,SharedFaceFluxRHS}, AdvanceADER, friction_substep,
# step); the Mult / Compute{,Shared}FaceFluxRHS regions are RK-only and are
# correctly ABSENT.  NOTE: with --shared-ck-recursion, ComputeADERSubStepStates
# is replaced by ComputeADERSubStepStatesAndIntegral and the second
# ApplySpatialDerivative subtree under AdvanceADER disappears.
# Set TPV104_PERFGRAPH=0 to disable.
if [[ "${TPV104_PERFGRAPH:-1}" == "1" ]]; then
    # runtime-report         -> human-readable inclusive/exclusive timing tree to stderr (.err)
    # hatchet-region-profile -> JSON for a rendered call-graph via the `hatchet` Python pkg
    export CALI_CONFIG="${CALI_CONFIG:-runtime-report(calc.inclusive=true),hatchet-region-profile(output=${OUT}/caliper_perfgraph.json)}"
    echo "[tpv104-p1] Caliper perfgraph ON (CALI_CONFIG=${CALI_CONFIG})."
    echo "            Produces data ONLY on an MFEM_USE_CALIPER=YES build."
fi

# --- run -----------------------------------------------------------------
# Scheme pinned on the CLI (--ader-order 2 / --time-integrator ader /
# --mixed-flux none); scalar Riemann + p1 come from the (asserted) TOML.
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
    ${DERIV_OPT_FLAGS} \
    ${RESTART_FLAGS} \
    --verify-dispatch

echo "[tpv104-p1] done. On-fault SCEC station traces: ${OUT}/tpv104_station_*.dat"
echo "            (9 stations); compare vs the SCEC TPV104 reference.  Confirm"
echo "            V_max is BOUNDED and read the empirical steps/hour off this"
echo "            *.out to size the restart chain."

# =========================================================================
# RESTART CHAIN (full 12 s may span multiple 48 h normal-queue segments)
# -------------------------------------------------------------------------
# The run checkpoints every 2000 steps (--checkpoint-every 2000) into ${OUT}.
# To continue after a wall-limit truncation, re-submit pointing TPV104_RESTART
# at the PREVIOUS segment's cp prefix:
#       TPV104_RESTART="tpv104/out_p1_aderO2_pu_normal_<prevjobid>/cp" \
#           sbatch tpv104_p1_aderO2_pureupwind_200m_normal.sbatch
# Each segment gets its OWN job-ID output dir (REQUIRED: the driver refuses a
# restart whose --output-dir equals the checkpoint's parent dir).
# =========================================================================
