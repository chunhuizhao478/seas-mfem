#!/bin/bash
#SBATCH -J tpv31_p2_aderO3
#SBATCH -o tpv31_p2_aderO3_%j.out
#SBATCH -e tpv31_p2_aderO3_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# =========================================================================
# SCEC TPV31 — p2 · ADER-O3 · no mixed flux · matrix Riemann
# Driver: seas_spatial_dyn_driver.  Sibling of the p1/ADER-O2 job
# (tpv31_job1_p1_aderO2_noflux_50m.sbatch); same physics, higher resolution.
#
# NUMERICAL SCHEME --------------------------------------------------------
#   * p2            : [mesh].order = 2            (quadratic FE space on the
#                     same straight-tet 50 m mesh; config-only, no CLI
#                     override — asserted below).  The DG CFL safety factor
#                     auto-scales: 1/(3·(2p+1)) = 1/15 at p=2 (vs 1/9 at p=1),
#                     so dt shrinks accordingly (spatial_friction.hpp:619).
#   * ADER-O3       : --ader-order 3 --time-integrator ader  (3rd-order ADER-DG
#                     sub-step; allowed range {2,3,4}).
#   * no mixed flux : --mixed-flux none.
#   * matrix Riemann: [numerics].interior_flux = "matrix"  (MANDATORY for
#                     TPV31's depth_profile_1d bi-material; config-only, no CLI
#                     override — asserted below) → BimaterialWaveOperator.
#   Config: tpv31/configs/tpv31_p2.toml (p2/O3 variant of tpv31.toml).
#
# QUEUE / SCALE -----------------------------------------------------------
#   flex queue, 16 nodes, 800 MPI ranks (50 ranks/node on Frontera CLX 56-core
#   nodes).  flex = discounted/opportunistic scheduling.  NOTE: there is no
#   in-repo flex-queue precedent — CONFIRM the current flex limits (max wall,
#   node range, SU rate, any preemption/hold behaviour) against the TACC
#   Frontera user guide before submitting.  p2/O3 is MUCH heavier than p1/O2
#   (more DOFs/elem + ~1.7× more steps from the tighter CFL + O3 work), so the
#   15 s run will span the wall limit — rely on checkpoint/restart (below).
#
# Reference overlays: tpv31/benchmark_data/scec_{eqdyna,seisol}/.
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

CONFIG_TOML="${TPV31_CONFIG:-tpv31/configs/tpv31_p2.toml}"
MESH="${TPV31_MESH:-tpv31/mesh/tpv31_50m.msh}"

if [[ ! -f "${CONFIG_TOML}" ]]; then
    echo "ERROR: config ${CONFIG_TOML} not found (cwd $(pwd))."; exit 1
fi

# --- pre-flight: assert the two CONFIG-ONLY knobs (no CLI override) --------
# interior_flux="matrix" is MANDATORY for TPV31; [mesh].order=2 is the "p2"
# in this job's name.  Both are config-only, so a drifted TOML would silently
# change the scheme.  Fail loudly here.
if ! grep -Eq '^[[:space:]]*interior_flux[[:space:]]*=[[:space:]]*"matrix"' "${CONFIG_TOML}"; then
    echo "ERROR: ${CONFIG_TOML} does not set [numerics].interior_flux = \"matrix\"."
    echo "       TPV31 (depth_profile_1d bi-material) REQUIRES the matrix Riemann path."
    exit 4
fi
if ! grep -Eq '^[[:space:]]*order[[:space:]]*=[[:space:]]*2([[:space:]]|$)' "${CONFIG_TOML}"; then
    echo "ERROR: ${CONFIG_TOML} does not set [mesh].order = 2 (this job is p2)."
    echo "       Use tpv31/configs/tpv31_p2.toml (not tpv31.toml, which is p1)."
    exit 4
fi

# --- regenerate the Gmsh v2.2 mesh if it is missing -----------------------
# The .msh is NOT committed (only the .geo).  MFEM's reader is Gmsh v2.2 only,
# so the `-format msh22` flag is MANDATORY.  Pre-build it ONCE on a login/idev
# node (conda activate pythonenv) so the compute job does not spend flex wall
# time meshing.
if [[ ! -f "${MESH}" ]]; then
    echo "[tpv31-p2] mesh ${MESH} missing — generating from the .geo (Gmsh v2.2)..."
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
# — this job (7757861) was SLURM-killed at the 2 h TIME LIMIT mid-dry-run.  The
# ibrun run below does the SAME construction + the SAME C++ config guards
# (FaultIteratorSupported, material-vs-interior_flux) in parallel and carries
# --verify-dispatch, so the resolved scheme is still printed.  The TOML greps
# above already fail-fast on a drifted interior_flux / [mesh].order pre-launch.

OUT="${TPV31_OUT:-tpv31/out_p2_${SLURM_JOB_ID:-local}}"
TFINAL="${TPV31_TFINAL:-15}"   # spec p. 9 = 15 s (checkpoint/restart across
                               # the flex wall limit)

# --- optional Caliper "perfgraph" (MFEM's profiling facility) -------------
# MFEM's performance-graph facility is Caliper: general/annotation.hpp defines
# MFEM_PERF_FUNCTION / MFEM_PERF_SCOPE / MFEM_PERF_BEGIN/END, which are no-ops
# UNLESS MFEM is built with MFEM_USE_CALIPER=YES and linked against libcaliper.
# A region-timing TREE (the "perfgraph") is then produced at runtime via the
# CALI_CONFIG env var.
#
# PREREQUISITES:
#   (1) Build with MFEM_USE_CALIPER=YES (+ libcaliper).  NOT the default yet
#       (config/config.mk has MFEM_USE_CALIPER=NO).  Build it with:
#         USE_CALIPER=YES CALIPER_DIR=<prefix> bash build_frontera.sh
#       (force-cleans the seas build so the annotations below compile in).
#   (2) MFEM_PERF_SCOPE annotations in the dynamic-rupture hot loops —
#       DONE as of the caliper-perfgraph merge (seas:: regions in
#       dynamic/wave_operator.inl, bimaterial_wave_operator.inl, and
#       drivers/spatial_dyn_driver.cpp).  Recipe + full region list:
#       document/caliper_perfgraph_dev/README_caliper_build_and_run.md.
#
# So on a default (MFEM_USE_CALIPER=NO) binary this block is INERT (the
# annotations are no-ops, CALI_CONFIG is ignored); on a Caliper build it emits
# the perfgraph.  NOTE: this is an ADER job, so the ADER-path regions fire
# (ComputeADER{SubStepStates,VolumeUpdate,FaceFluxRHS,SharedFaceFluxRHS},
# AdvanceADER, friction_substep, step); the Mult / Compute{,Shared}FaceFluxRHS
# regions are RK-only and are correctly ABSENT here.
#
# Run with TPV31_PERFGRAPH=1:
if [[ "${TPV31_PERFGRAPH:-0}" == "1" ]]; then
    # runtime-report      -> human-readable inclusive/exclusive timing tree to stderr (.err)
    # hatchet-region-profile -> JSON for a rendered call-graph via the `hatchet` Python pkg
    export CALI_CONFIG="${CALI_CONFIG:-runtime-report(calc.inclusive=true),hatchet-region-profile(output=${OUT}/caliper_perfgraph.json)}"
    echo "[tpv31-p2] Caliper perfgraph ON (CALI_CONFIG=${CALI_CONFIG})."
    echo "           Produces data ONLY on an MFEM_USE_CALIPER=YES build with"
    echo "           MFEM_PERF_*-annotated hot loops (see note above)."
fi

# --- run -----------------------------------------------------------------
# Scheme pinned on the CLI (--ader-order 3 / --time-integrator ader /
# --mixed-flux none); matrix Riemann + p2 come from the (asserted) TOML.
ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --ader-order 3 \
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

echo "[tpv31-p2] done. On-fault station traces in ${OUT}; compare vs"
echo "           tpv31/benchmark_data/scec_{eqdyna,seisol}/ via tpv31/visualize_results.py."

# =========================================================================
# RESTART (flex wall limit)
# -------------------------------------------------------------------------
# The run checkpoints every 2000 steps (--checkpoint-every 2000) to ${OUT}.
# Resume a truncated run by re-submitting with the same ${OUT} and adding:
#       --restart "${OUT}/cp"
# (mirrors the tpv104 restart jobs).  p2/O3 at 15 s will need several flex
# segments; read the steps/hour off this job's *.out to size the chain.
# =========================================================================
