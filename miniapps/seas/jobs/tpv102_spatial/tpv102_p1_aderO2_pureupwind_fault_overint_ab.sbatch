#!/bin/bash
#SBATCH -J tpv102_overint_ab
#SBATCH -o tpv102_overint_ab_%j.out
#SBATCH -e tpv102_overint_ab_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# =========================================================================
# TPV102 — FAULT-FLUX OVER-INTEGRATION A/B  (Phase 1+4 of the fault-dealiasing
# plan: document/fault_dealiasing_dev/seisol_overintegration_resample_speckle_2026-06-02.md)
#
# ONE parameterized job, submitted TWICE (OFF then ON) so the two runs differ in
# EXACTLY one knob — the fault-flux over-integration factor:
#     FAULT_OVERINT=0  -> baseline (degree 2*order fault quadrature; byte-exact
#                         with pre-Phase-4 = the on-fault sigma_n SPECKLE/DRIFT
#                         the plan targets).
#     FAULT_OVERINT=2  -> --fault-overint 2 (degree 2*(order+2)=6 fault
#                         quadrature; the per-step "speckle lever", plan §4.1).
# Everything else is identical: p1 · ADER-O2 · pure upwind · scalar Riemann,
# same TOML, same mesh, same optimization flags.  Submit via
# submit_fault_overint_ab.sh (chains build -> OFF -> ON), or by hand:
#     FAULT_OVERINT=0 sbatch tpv102_p1_aderO2_pureupwind_fault_overint_ab.sbatch
#     FAULT_OVERINT=2 sbatch tpv102_p1_aderO2_pureupwind_fault_overint_ab.sbatch
# The output dir auto-encodes the factor, so the two runs never collide.
#
# WHAT TO COMPARE: the on-fault SCEC station traces ${OUT}/tpv102_station_*.dat
# (sigma_n column).  Baseline drifts/speckles (TPV102: 120 -> 80, -30, 75 MPa);
# the prediction is that FAULT_OVERINT=2 stays much closer to the constant
# 120 MPa — i.e. over-integration removes the per-step flux aliasing.
#
# NUMERICAL SCHEME (pinned; matches tpv102_p1_aderO2_pureupwind_200m_normal):
#   p1 (config), ADER-O2 (--ader-order 2 --time-integrator ader), pure upwind
#   (--mixed-flux none), scalar Riemann (config interior_flux="scalar").
#   --fault-overint requires pure-upwind + no precomputed flux (both hold here);
#   it is honoured on the scalar AND matrix paths (the fault-flux assembly lives
#   in base WaveOperator).
#
# CHEAP FIRST PASS (recommended before the full 200 m / 12 s production A/B):
#   run on the coarse 1000 m mesh, flex queue, short tfinal — the speckle is a
#   per-element artifact and shows up at 1000 m too (the Phase-0 harness uses
#   ~1 km elements).  See README_fault_overint_ab.md.
#       FAULT_OVERINT=0 TPV102_MESH=tpv102/mesh/tpv102_1000m.msh TPV102_TFINAL=5 \
#           sbatch -p flex -N 2 -n 100 -t 02:00:00 \
#           tpv102_p1_aderO2_pureupwind_fault_overint_ab.sbatch
#       (and again with FAULT_OVERINT=2)
#
# RS-AT-np>1 GATE: the driver prints "Treat parallel RS results as PRELIMINARY"
#   at np>1 (shared rank-seam fault QPs use end-of-step psi).  Both A/B runs see
#   the SAME gate, so the COMPARISON is still valid.
# =========================================================================

set -euo pipefail

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

# Pin the hypre/2.31.0 runtime libHYPRE.so at the FRONT of the search path
# (petsc/3.15 + mumps/5.3 prepend a pre-2.31 libHYPRE; the binary needs the
# 2.31.0 symbol HYPRE_Initialize).  Matches every spatial run job in this repo.
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
    echo "       Build it FIRST from the feat/fault-overint-resample branch:"
    echo "         git fetch origin && git checkout feat/fault-overint-resample"
    echo "         sbatch miniapps/seas/jobs/spatial_dyn_build.sbatch clean-rebuild"
    echo "       (the run job does NOT build; see README_fault_overint_ab.md)."
    exit 1
fi

# --- the ONE A/B knob -----------------------------------------------------
FAULT_OVERINT="${FAULT_OVERINT:-0}"
if ! [[ "${FAULT_OVERINT}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: FAULT_OVERINT must be a non-negative integer, got '${FAULT_OVERINT}'."
    exit 4
fi

CONFIG_TOML="${TPV102_CONFIG:-tpv102/configs/tpv102_spatial.toml}"
# Production resolution.  The TOML default is the coarse 1000 m mesh; the full
# production A/B overrides to 200 m (set TPV102_MESH for the cheap 1000 m pass).
MESH="${TPV102_MESH:-tpv102/mesh/tpv102_200m.msh}"

if [[ ! -f "${CONFIG_TOML}" ]]; then
    echo "ERROR: config ${CONFIG_TOML} not found (cwd $(pwd))."; exit 1
fi

# --- pre-flight: assert the two CONFIG-ONLY knobs (no CLI override) --------
# p1 (=[mesh].order=1) and interior_flux="scalar" are config-only; a drifted
# TOML would silently change the scheme and break the A/B equivalence.
if ! grep -Eq '^[[:space:]]*order[[:space:]]*=[[:space:]]*1([[:space:]]|$)' "${CONFIG_TOML}"; then
    echo "ERROR: ${CONFIG_TOML} does not set [mesh].order = 1 (this job is p1)."
    exit 4
fi
if ! grep -Eq '^[[:space:]]*interior_flux[[:space:]]*=[[:space:]]*"scalar"' "${CONFIG_TOML}"; then
    echo "ERROR: ${CONFIG_TOML} does not set [numerics].interior_flux = \"scalar\"."
    exit 4
fi

# --- regenerate the Gmsh v2.2 mesh if missing (MFEM reader is v2.2-only) ---
# Derive the .geo from the .msh name (tpv102_<res>.msh -> tpv102_<res>.geo).
GEO="${TPV102_GEO:-${MESH%.msh}.geo}"
if [[ ! -f "${MESH}" ]]; then
    echo "[overint-ab] mesh ${MESH} missing — generating from ${GEO} (Gmsh v2.2)..."
    if command -v gmsh >/dev/null 2>&1; then
        gmsh -format msh22 -3 "${GEO}" -o "${MESH}"
    else
        echo "ERROR: gmsh not on PATH.  Pre-build the mesh on a login node"
        echo "       (conda activate pythonenv): gmsh -format msh22 -3 ${GEO} -o ${MESH}"
        exit 2
    fi
fi

# Output dir auto-encodes the over-int factor so OFF/ON never collide.
OUT="${TPV102_OUT:-tpv102/out_p1_aderO2_pu_overint${FAULT_OVERINT}_${SLURM_JOB_ID:-local}}"
TFINAL="${TPV102_TFINAL:-12}"   # SCEC TPV102 window 0-12 s

# --- the over-int CLI flag (the A/B knob) ---------------------------------
# FAULT_OVERINT=0 -> pass --fault-overint 0 (explicit baseline; the driver
# skips SetFaultOverint at K=0, so this is byte-exact with omitting the flag).
FAULT_OVERINT_FLAG="--fault-overint ${FAULT_OVERINT}"

# --- ADER hot-path optimization (ON in BOTH A/B runs => not a confound) ----
DERIV_OPT_FLAGS="--deriv-cache --shared-ck-recursion"

# --- restart across the wall (opt-in via TPV102_RESTART) -------------------
RESTART_FLAGS=""
if [[ -n "${TPV102_RESTART:-}" ]]; then
    RESTART_FLAGS="--restart ${TPV102_RESTART}"
    echo "[overint-ab] RESTART: resuming from '${TPV102_RESTART}' (this segment writes ${OUT})."
fi

echo "[overint-ab] FAULT_OVERINT=${FAULT_OVERINT}  config=${CONFIG_TOML}  mesh=${MESH}"
echo "[overint-ab] tfinal=${TFINAL}  output=${OUT}"
echo "[overint-ab] Git: $(git rev-parse --abbrev-ref HEAD 2>/dev/null) $(git rev-parse --short HEAD 2>/dev/null)"

# --- run ------------------------------------------------------------------
ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --ader-order 2 \
    --time-integrator ader \
    --mixed-flux none \
    ${FAULT_OVERINT_FLAG} \
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

echo "[overint-ab] done (FAULT_OVERINT=${FAULT_OVERINT}).  Station traces: ${OUT}/tpv102_station_*.dat"
echo "[overint-ab] grep the *.out for '[fault] over-integration ON' to confirm the"
echo "             ON run actually applied it (absent => FAULT_OVERINT=0 baseline)."
echo "             Compare the sigma_n column OFF vs ON at each station."
