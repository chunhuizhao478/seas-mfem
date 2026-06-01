#!/bin/bash
#SBATCH -J tpv31_p3_aderO4_n
#SBATCH -o tpv31_p3_aderO4_n_%j.out
#SBATCH -e tpv31_p3_aderO4_n_%j.err
#SBATCH -p normal
#SBATCH -N 16
#SBATCH -n 800
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# =========================================================================
# SCEC TPV31 — JOB [3] (PRODUCTION):  p3 · ADER-O4 · pure upwind · matrix Riemann
# Driver: seas_spatial_dyn_driver (spatial dynamic-rupture).  Normal-queue,
# full-15 s, highest-resolution sibling of the p1/O2 and p2/O3 normal jobs
# (tpv31_p1_aderO2_noflux_50m_normal.sbatch, tpv31_p2_aderO3_noflux_50m_normal.sbatch).
#
# NUMERICAL SCHEME (the four knobs that define this job) ------------------
#   * p3            : [mesh].order = 3            (cubic FE space on the same
#                     straight-tet 50 m mesh; config-only, no CLI override —
#                     asserted below).  DG CFL safety auto-scales:
#                     1/(3·(2p+1)) = 1/21 at p=3 (vs 1/15 at p=2, 1/9 at p=1),
#                     so dt shrinks accordingly (spatial_friction.hpp:619).
#   * ADER-O4       : --ader-order 4 --time-integrator ader  (4th-order ADER-DG
#                     sub-step; O4 is the MAXIMUM of the allowed range {2,3,4}).
#   * pure upwind   : --mixed-flux none           (upwind Godunov on EVERY
#                     interior face; no central-flux faces).
#   * matrix Riemann: [numerics].interior_flux = "matrix"  (MANDATORY for
#                     TPV31's depth_profile_1d bi-material; config-only, no CLI
#                     override — asserted below) → BimaterialWaveOperator.
#   Config: tpv31/configs/tpv31_p3.toml (p3 / ADER-O4 variant of tpv31.toml).
#
# QUEUE / SCALE -----------------------------------------------------------
#   normal queue, 16 nodes, 800 MPI ranks (50 ranks/node on Frontera CLX
#   56-core nodes).  Wall = 48 h (normal-queue maximum).  p3/O4 is the heaviest
#   of the three (cubic DOFs/elem + tightest CFL + O4 sub-step work), so the
#   15 s spec run will span SEVERAL 48 h segments — rely on checkpoint/restart
#   (see RESTART below).
#
# ADER HOT-PATH OPTIMIZATION (ALWAYS ON here) -----------------------------
#   This production job ALWAYS passes the Lever 1/2/3 optimization flags
#   (--deriv-cache --shared-ck-recursion); see "OPTIMIZATION" below.  Result
#   differs from the bit-exact default at ROUND-OFF only (machine-eps; verified
#   by the equivalence unit test, NOT the byte-exact regression contract).
#   ⇒ The binary MUST be (re)built from this branch (it carries the merged
#     elem_derivative_cache.hpp / DerivMode::Cached code).  A binary built
#     before the merge IGNORES these flags SILENTLY (HasFlag argv-scan) and
#     runs the slow path — grep the *.out for "[deriv] ... enabled" to confirm.
#
# PERFGRAPH (Caliper; ON by default here) ---------------------------------
#   TPV31_PERFGRAPH defaults to 1 so the Caliper region tree is collected.
#   This requires an MFEM_USE_CALIPER=YES build; on a non-Caliper binary the
#   annotations are no-ops and CALI_CONFIG is ignored (the run still proceeds).
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
    echo "       ADER optimization).  For perfgraph + optimization, do a CLEAN"
    echo "       Caliper build (the WaveOperator templates changed):"
    echo "         cd miniapps/seas && make clean && \\"
    echo "           USE_CALIPER=YES CALIPER_DIR=<prefix> bash build_frontera.sh"
    exit 1
fi

CONFIG_TOML="${TPV31_CONFIG:-tpv31/configs/tpv31_p3.toml}"
MESH="${TPV31_MESH:-tpv31/mesh/tpv31_50m.msh}"

if [[ ! -f "${CONFIG_TOML}" ]]; then
    echo "ERROR: config ${CONFIG_TOML} not found (cwd $(pwd))."; exit 1
fi

# --- pre-flight: assert the two CONFIG-ONLY knobs (no CLI override) --------
# interior_flux="matrix" is MANDATORY for TPV31; [mesh].order=3 is the "p3"
# in this job's name.  Both are config-only, so a drifted TOML would silently
# change the scheme.  Fail loudly here.
if ! grep -Eq '^[[:space:]]*interior_flux[[:space:]]*=[[:space:]]*"matrix"' "${CONFIG_TOML}"; then
    echo "ERROR: ${CONFIG_TOML} does not set [numerics].interior_flux = \"matrix\"."
    echo "       TPV31 (depth_profile_1d bi-material) REQUIRES the matrix Riemann path."
    exit 4
fi
if ! grep -Eq '^[[:space:]]*order[[:space:]]*=[[:space:]]*3([[:space:]]|$)' "${CONFIG_TOML}"; then
    echo "ERROR: ${CONFIG_TOML} does not set [mesh].order = 3 (this job is p3)."
    echo "       Use tpv31/configs/tpv31_p3.toml (not tpv31.toml/tpv31_p2.toml)."
    exit 4
fi

# --- regenerate the Gmsh v2.2 mesh if it is missing -----------------------
# The .msh is NOT committed (only the .geo).  MFEM's reader is Gmsh v2.2 only,
# so the `-format msh22` flag is MANDATORY.  Pre-build it ONCE on a login/idev
# node (conda activate pythonenv) so the compute job does not spend wall time
# meshing.
if [[ ! -f "${MESH}" ]]; then
    echo "[tpv31-p3] mesh ${MESH} missing — generating from the .geo (Gmsh v2.2)..."
    if command -v gmsh >/dev/null 2>&1; then
        gmsh -format msh22 -3 tpv31/mesh/tpv31_50m.geo -o "${MESH}"
    else
        echo "ERROR: gmsh not on PATH.  Generate the mesh first (conda activate"
        echo "       pythonenv): gmsh -format msh22 -3 tpv31/mesh/tpv31_50m.geo"
        echo "       -o ${MESH}"
        exit 2
    fi
fi

# --- positional arg: `perfgraph` (Caliper short slice that flushes) --------
# Submit with the single word `perfgraph` to run a SHORT slice (default
# ~0.004 s) that REACHES a clean exit so the Caliper runtime-report actually
# flushes:   sbatch tpv31_p3_aderO4_noflux_50m_normal.sbatch perfgraph
# A wall-clock-killed run (CANCELLED ... DUE TO TIME LIMIT) emits NO report,
# so the full 15 s run's Caliper data only covers whichever restart segment
# reaches t=tfinal cleanly — use this short slice for a representative tree.
# Slice-length override: TPV31_PERFGRAPH_TFINAL (s).  (p3 dt is smaller than
# p2's because CFL safety is 1/21 at p=3, so 0.004 s is MORE steps.)
PERFGRAPH_MODE=0
for _arg in "$@"; do
    case "${_arg}" in
        perfgraph) PERFGRAPH_MODE=1 ;;
        *) echo "[tpv31-p3] WARNING: ignoring unknown arg '${_arg}' (expected 'perfgraph')." ;;
    esac
done

OUT="${TPV31_OUT:-tpv31/out_p3_aderO4_normal_${SLURM_JOB_ID:-local}}"

if [[ "${PERFGRAPH_MODE}" == "1" ]]; then
    TFINAL="${TPV31_TFINAL:-${TPV31_PERFGRAPH_TFINAL:-0.004}}"   # short -> finishes -> report flushes
    echo "[tpv31-p3] PERFGRAPH MODE: profiling slice tfinal=${TFINAL} s."
    echo "           Runs to a CLEAN exit so the Caliper report flushes; the full run is tfinal=15 s."
else
    TFINAL="${TPV31_TFINAL:-15}"   # spec p. 9 = 15 s (checkpoint/restart across the 48 h wall)
fi

# --- ADER hot-path optimization (ALWAYS ON in this production job) ---------
# --deriv-cache         : precompute per-element D_d^e = M^-1 K_d (Lever 1) and
#                         the volume-RHS operator S_d^e (Lever 3); apply as a
#                         dense mat-vec instead of per-QP quadrature
#                         (DerivMode::Cached).  Changes results at ROUND-OFF
#                         only (machine-eps; equivalence unit test, NOT the
#                         byte-exact contract).  Cache budget (R-004) at p3 is
#                         ~72 MB/rank — still under the 1 GiB/rank default
#                         (≈3.6 GB/node at 50 ranks/node; CLX nodes have 192 GB).
# --shared-ck-recursion : run the Cauchy-Kovalevskaya predictor recursion ONCE
#                         per macro-step instead of twice (Lever 2).
#                         BYTE-IDENTICAL at a fixed DerivMode.
# The driver prints "[deriv] ... enabled" on rank 0 when the flags take effect
# — grep the *.out to confirm the binary actually carries the optimization.
DERIV_OPT_FLAGS="--deriv-cache --shared-ck-recursion"

# --- restart across the 48 h wall (opt-in via TPV31_RESTART) ---------------
# The driver ABORTS (return 3) if --output-dir resolves to the SAME directory
# as the PARENT of --restart, so a restart segment MUST write to a DIFFERENT
# --output-dir than the checkpoint it reads.  Because OUT defaults to a
# per-job-ID directory, each segment is naturally distinct.  To resume:
#       TPV31_RESTART="<previous-segment-OUT>/cp" sbatch <thisscript>
# (point it at the cp prefix inside the PREVIOUS job's output dir).  This
# segment then checkpoints into its own fresh OUT for the next link.
RESTART_FLAGS=""
if [[ -n "${TPV31_RESTART:-}" ]]; then
    RESTART_FLAGS="--restart ${TPV31_RESTART}"
    echo "[tpv31-p3] RESTART: resuming from checkpoint prefix '${TPV31_RESTART}'"
    echo "           (this segment writes a fresh output dir: ${OUT})."
fi

# --- Caliper "perfgraph" (MFEM profiling facility); ON by default ----------
# MFEM's perfgraph is Caliper: MFEM_PERF_SCOPE annotations (seas:: regions in
# dynamic/wave_operator.inl, bimaterial_wave_operator.inl, and
# drivers/spatial_dyn_driver.cpp) are no-ops UNLESS MFEM is built with
# MFEM_USE_CALIPER=YES + libcaliper.  A region-timing TREE is then produced at
# runtime via CALI_CONFIG.  Recipe: document/caliper_perfgraph_dev/
# README_caliper_build_and_run.md.  This is an ADER job, so the ADER-path
# regions fire (ComputeADER{SubStepStates,VolumeUpdate,FaceFluxRHS,
# SharedFaceFluxRHS}, AdvanceADER, friction_substep, step); the Mult /
# Compute{,Shared}FaceFluxRHS regions are RK-only and are correctly ABSENT.
# NOTE: with --shared-ck-recursion, ComputeADERSubStepStates is replaced by
# ComputeADERSubStepStatesAndIntegral and the second ApplySpatialDerivative
# subtree under AdvanceADER disappears.
# Set TPV31_PERFGRAPH=0 to disable.
if [[ "${TPV31_PERFGRAPH:-1}" == "1" ]]; then
    # runtime-report         -> human-readable inclusive/exclusive timing tree to stderr (.err)
    # hatchet-region-profile -> JSON for a rendered call-graph via the `hatchet` Python pkg
    export CALI_CONFIG="${CALI_CONFIG:-runtime-report(calc.inclusive=true),hatchet-region-profile(output=${OUT}/caliper_perfgraph.json)}"
    echo "[tpv31-p3] Caliper perfgraph ON (CALI_CONFIG=${CALI_CONFIG})."
    echo "           Produces data ONLY on an MFEM_USE_CALIPER=YES build."
fi

# --- run -----------------------------------------------------------------
# Scheme pinned on the CLI (--ader-order 4 / --time-integrator ader /
# --mixed-flux none); matrix Riemann + p3 come from the (asserted) TOML.
ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --ader-order 4 \
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

echo "[tpv31-p3] done. On-fault station traces in ${OUT}; compare vs"
echo "           tpv31/benchmark_data/scec_{eqdyna,seisol}/ via tpv31/visualize_results.py."

# =========================================================================
# RESTART CHAIN (full 15 s spans multiple 48 h normal-queue segments)
# -------------------------------------------------------------------------
# The run checkpoints every 2000 steps (--checkpoint-every 2000) into ${OUT}
# as cp_checkpoint_r<rank>.txt.  To continue after a wall-limit truncation,
# re-submit pointing TPV31_RESTART at the PREVIOUS segment's cp prefix:
#       TPV31_RESTART="tpv31/out_p3_aderO4_normal_<prevjobid>/cp" \
#           sbatch tpv31_p3_aderO4_noflux_50m_normal.sbatch
# Each segment gets its OWN job-ID output dir (REQUIRED: the driver refuses a
# restart whose --output-dir equals the checkpoint's parent dir).  Read the
# steps/hour off the *.out to size the number of segments.
# =========================================================================
