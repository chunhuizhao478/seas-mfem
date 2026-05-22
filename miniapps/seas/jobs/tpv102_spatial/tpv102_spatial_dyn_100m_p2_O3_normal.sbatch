#!/bin/bash
#SBATCH -J tpv102_spdyn_100m_p2_O3_norm
#SBATCH -o tpv102_spdyn_100m_p2_O3_norm_%j.out
#SBATCH -e tpv102_spdyn_100m_p2_O3_norm_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# TPV102 (SCEC FL=101 rate-and-state) 100 m mesh, P=2 + ADER-O3,
# NORMAL-QUEUE production run, heterogeneous Riemann.
# =======================================================================
#
# Runs TPV102 through the TOML-driven seas_spatial_dyn_driver (NOT the
# native drivers/tpv102_driver.cpp) at polynomial order 2 (P2) with ADER
# time order 3 (Pelties O3).  This is the heterogeneous-Riemann path
# verification for TPV102 (R-001): on a homogeneous Mode::Constant
# material the bimaterial interior flux is byte-identical to the scalar
# ctor, so a clean match against the native TPV102 reference is expected.
#
# How P2/O3 is selected
# ---------------------
# No --order CLI flag exists; the mesh polynomial order is set ONLY via
# the TOML [mesh].order field.  We build a sed-override config that flips
# `order = 1` -> `order = 2` (override file under RESULT_DIR for
# provenance).  ADER order IS a CLI flag → `--ader-order 3` (validated in
# {2,3,4}).  The sed-override rewrites BOTH the mesh path (1000m -> 100m) AND
# the order, exactly like the tpv205_spatial jobs.
#
# Why the heterogeneous Riemann path
# ----------------------------------
# tpv102.toml omits [numerics].interior_flux → defaults to "bimaterial",
# so the heterogeneous (MaterialField, BoundaryConfig) ctor is used and
# the banner must print `[wave_operator] BimaterialFlux precomputation: ...`.
# The material is homogeneous (rho=2670, cs=3464, cp=6000), so the cached
# per-face GodunovFlux is built from the exact (λ,μ,ρ) and the result is
# byte-identical to the scalar path (CLAUDE.md "Phase H Stage 2").
#
# Spec encoding (verbatim per tpv102/configs/tpv102.toml):
#   * Material:   homogeneous (μ=32.038119520e9, λ=32.043760960e9, ρ=2670).
#   * Friction:   regularised rate-and-state, AGEING law (FL=101).
#                 f0=0.6, V0=1e-6, b=0.012, Dc=0.02; `a` via C∞ tanh
#                 boxcar (a_vw=0.008 inside |x|≤15 km ∧ |z+7.5 km|≤7.5 km,
#                 a_vs=0.016 outside, 3 km transition).
#   * Loading:    σ_n=120 MPa (compression POSITIVE), τ_ini=75 MPa
#                 right-lateral along-strike, V_init=1e-12 m/s.
#   * Nucleation: SCEC compact bell F(r)=exp(r²/(r²−R²)) (R=3 km) ×
#                 smoothStep G(t) (T=1 s); peak Δτ₀=25 MPa, strike-only.
#   * Fault:      vertical y=0 plane; hypocenter (0,0,−7500).
#
# Rank count / wall time
# ----------------------
# Fixed at 10 N x 500 r.  The BimaterialFlux precomputation footprint is
# per-FACE (order-independent), so memory headroom is unchanged from P1.
# At P2 the per-step cost is ~3-5x P1 and the CFL tightens, so the 48 h
# budget reaches less simulated time than a P1/O2 run; chain 48 h blocks
# via `--restart <prefix>` (checkpoint every 10000 steps per tpv102.toml
# [output].checkpoint_every_steps) to advance further.

export LC_ALL=C
export LANG=C

# Fail fast on module-load / build / pipe failures.
set -euo pipefail

SBATCH_LOG_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

# Project checkout root.  Auto-detect by walking up from $SLURM_SUBMIT_DIR
# (or $PWD if not running under Slurm) until we find build_frontera.sh.
# Override at submit time when running from a different clone:
#   SEAS_MFEM_ROOT=$(pwd) sbatch ...
if [ -z "${SEAS_MFEM_ROOT:-}" ]; then
    _probe="${SLURM_SUBMIT_DIR:-$PWD}"
    while [ "${_probe}" != "/" ] && [ -n "${_probe}" ]; do
        if [ -f "${_probe}/build_frontera.sh" ]; then
            SEAS_MFEM_ROOT="${_probe}"
            break
        fi
        _probe="$(dirname "${_probe}")"
    done
    unset _probe
fi
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-/scratch2/10024/zhaochun/seas-project/seas-mfem}"
echo "  SEAS_MFEM_ROOT = ${SEAS_MFEM_ROOT}"

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8
# NOTE: do NOT load phdf5 on Frontera.  The phdf5/1.14.x module requires
# intel/23.1.0 + impi/21.9.0 (incompatible with the intel/19 stack);
# build_frontera.sh builds HDF5 1.14.x from source into
# ${SEAS_MFEM_ROOT}/extern/hdf5/install instead; we link/run against that.

HDF5_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/hdf5/install"
H5Z_ZFP_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install"

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${HDF5_BUILT_DIR}/lib:${LD_LIBRARY_PATH:-}"

# REQUIRED for tpv102.toml's `paraview_fault_zfp_tol = 1e-12`: the HDF5
# ZFP filter plugin must be loadable at run time.
if [ -d "${H5Z_ZFP_BUILT_DIR}/plugin" ]; then
   export HDF5_PLUGIN_PATH="${H5Z_ZFP_BUILT_DIR}/plugin"
else
   echo "WARNING: ${H5Z_ZFP_BUILT_DIR}/plugin not found — run build_frontera.sh first."
fi

cd "${SEAS_MFEM_ROOT}"

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES'    config/config.mk || { echo "ERROR: no PETSc";   exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'    config/config.mk || { echo "ERROR: no MUMPS";   exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'     config/config.mk || { echo "ERROR: no HDF5 (paraview_fault=hdf5 needs it)"; exit 1; }
# tpv102.toml's `paraview_*_zfp_tol` fields make the driver abort at parse
# time without H5Z-ZFP — fail before ibrun so we don't waste the slot.
grep -q '^MFEM_USE_H5Z_ZFP *= YES'  config/config.mk || { echo "ERROR: no H5Z-ZFP (tpv102.toml paraview_*_zfp_tol requires it)"; exit 1; }

cd miniapps/seas

# Regenerate the 100 m .msh from .geo on first run or when FORCE_MESHGEN=1.
MESH_SRC="tpv102/mesh/tpv102_100m.geo"
MESH_OUT="tpv102/mesh/tpv102_100m.msh"
GMSH_BIN="${SEAS_MFEM_ROOT}/extern/gmsh/bin/gmsh"
if [ ! -s "${MESH_OUT}" ] || [ "${FORCE_MESHGEN:-0}" = "1" ] || [ "${FORCE_MESHGEN:-}" = "YES" ]; then
    if [ ! -x "${GMSH_BIN}" ]; then
        echo "ERROR: ${GMSH_BIN} not found.  Run build_frontera.sh (or"
        echo "       USE_GMSH=YES bash build_frontera.sh) on the login node first."
        exit 1
    fi
    if [ ! -s "${MESH_SRC}" ]; then
        echo "ERROR: ${MESH_SRC} missing from the checkout — cannot regenerate."
        exit 1
    fi
    echo "=== Regenerating ${MESH_OUT} from ${MESH_SRC} (FORCE_MESHGEN=${FORCE_MESHGEN:-0}) ==="
    # MFEM requires Gmsh v2.2 format (CLAUDE.md "Known limitation").
    "${GMSH_BIN}" -format msh22 -3 "${MESH_SRC}" -o "${MESH_OUT}" \
        || { echo "MESH GENERATION FAILED"; exit 1; }
fi

if [ ! -s "${MESH_OUT}" ]; then
   echo "ERROR: ${MESH_OUT} missing after mesh-regen step."
   exit 1
fi

# Build a 100m+P2 override config from the canonical tpv102.toml by sed-
# replacing the mesh path (1000m -> 100m) AND `order = 1` -> `order = 2`
# (no --order CLI flag).  Override file lives under RESULT_DIR.
RESULT_DIR="tpv102/results_spatial_dyn_100m_p2_O3_normal_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"
CONFIG_OVERRIDE="${RESULT_DIR}/tpv102_100m_p2.toml"
sed -e 's|tpv102_1000m\.msh|tpv102_100m.msh|' \
    -e 's|^order = 1|order = 2|' tpv102/configs/tpv102.toml \
   > "${CONFIG_OVERRIDE}"

# The seas_spatial_dyn_driver binary is built ONCE by the separate
# one-shot build job (jobs/build_spatial_dyn_driver_dev.sbatch).  Submit
# that FIRST and let it reach COMPLETED before submitting this run job.
# Run jobs do NOT build, so any number of them may run concurrently
# without racing on the shared object files / binary (the cause of the
# "drivers/spatial_dyn_driver.o: file truncated" link error).
if [ ! -x seas_spatial_dyn_driver ]; then
   echo "ERROR: seas_spatial_dyn_driver not found in $(pwd)."
   echo "       Build it first (one-shot, development queue):"
   echo "         sbatch ${SEAS_MFEM_ROOT}/miniapps/seas/jobs/build_spatial_dyn_driver_dev.sbatch"
   echo "       Wait for COMPLETED, then resubmit this run job."
   exit 1
fi

echo "========================================"
echo "Git commit:    $(git rev-parse HEAD)"
echo "Git branch:    $(git rev-parse --abbrev-ref HEAD)"
echo "Driver:        seas_spatial_dyn_driver (TOML-driven)"
echo "Config:        ${CONFIG_OVERRIDE} (mesh 100m + order=2 override of tpv102.toml)"
echo "Mesh:          tpv102/mesh/tpv102_100m.msh"
echo "Polynomial:    order 2 (TOML mesh.order override) + ader-order 3 (Pelties O3)"
echo "Interior flux: bimaterial (heterogeneous Riemann; byte-identical on homogeneous material)"
echo "Friction law:  rate_state (regularised R&S, ageing law, FL=101)"
echo "Iterator:      one-shot (O=1 quadrature)"
echo "Queue:         normal (10N x 500r, 48h)"
echo "tfinal:        12.0 s (48 h budget may reach less; resume via --restart)"
echo "Checkpoint:    every 10000 steps (TOML default; resume via --restart)"
echo "Stations:      9 SCEC on-fault traces (prefix tpv102_)"
echo "ParaView:      fault=hdf5 every 0.001 s (TOML default)"
echo "========================================"

ibrun ./seas_spatial_dyn_driver \
      --config "${CONFIG_OVERRIDE}" \
      --output-dir "${RESULT_DIR}" \
      --ader-order 3 \
      --verify-dispatch

# -----------------------------------------------------------------------
# Post-run summary.  TPV102 writes SCEC on-fault traces as
# tpv102_station_<name>.dat (TPV102StationWriter, prefix tpv102).  The
# hypocenter station is flt_0_7.5 (0 km strike, 7.5 km down-dip).
# -----------------------------------------------------------------------
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv102_station_flt_0_7.5.dat"

PEAK_V_HYPO="N/A"
PEAK_T_HYPO="N/A"
SIGMA_N_PERTURB_HYPO="N/A"
SLIP_STRIKE_END="N/A"
T_END="N/A"
N_ROWS="N/A"

if [ -s "${HYPO_FILE}" ]; then
   read PEAK_V_HYPO PEAK_T_HYPO < <(awk '!/^#/ && NF>=3 && $1+0==$1 {
      v=$3; if (v<0) v=-v;
      if (v>m) { m=v; tm=$1 }
   } END { printf "%g %g\n", m, tm }' "${HYPO_FILE}")
   # σ_n perturbation relative to the 120 MPa background.
   SIGMA_N_PERTURB_HYPO=$(awk '!/^#/ && NF>=8 && $1+0==$1 {
      d=$8 - 120e6; if (d<0) d=-d;
      if (d>m) m=d
   } END { printf "%.4e", m+0 }' "${HYPO_FILE}")
   SLIP_STRIKE_END=$(awk '!/^#/ && NF>=2 && $1+0==$1 {a=$2} END{printf "%.4e", a+0}' "${HYPO_FILE}")
   T_END=$(awk '!/^#/ && NF>=1 && $1+0==$1 {tt=$1} END{printf "%.4f", tt+0}' "${HYPO_FILE}")
   N_ROWS=$(awk '!/^#/ && NF>=1 && $1+0==$1' "${HYPO_FILE}" | wc -l | tr -d ' ')
fi

DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv102_spdyn_100m_p2_O3_norm_${SLURM_JOB_ID}.out"

# Gate 1 — heterogeneous Riemann precomputation banner.
BIMAT_VERDICT="UNKNOWN"
if grep -q '\[wave_operator\] BimaterialFlux precomputation' "${DISPATCH_LOG}" 2>/dev/null; then
   BIMAT_VERDICT="ENGAGED — heterogeneous Riemann path active (R-001)"
elif grep -q '\[wave\] interior_flux=scalar' "${DISPATCH_LOG}" 2>/dev/null; then
   BIMAT_VERDICT="REGRESSION — scalar ctor selected"
else
   BIMAT_VERDICT="UNKNOWN — no interior-flux banner in log"
fi

# Gate 2 — ParaView fault output present.
PARAVIEW_VERDICT="UNKNOWN"
if [ -f "${RESULT_DIR}/fault.vtkhdf" ] || ls "${RESULT_DIR}"/*.vtkhdf 1>/dev/null 2>&1 || ls "${RESULT_DIR}"/*.pvd 1>/dev/null 2>&1; then
   PARAVIEW_VERDICT="ENGAGED — fault output written"
elif grep -q 'ParaView output: OFF' "${DISPATCH_LOG}" 2>/dev/null; then
   PARAVIEW_VERDICT="REGRESSION — driver banner reports 'OFF'; check tpv102.toml [output].paraview_enabled"
else
   PARAVIEW_VERDICT="UNKNOWN — no fault file and no 'OFF' banner"
fi

# Gate 3 — SCEC station traces written.
STATION_VERDICT="UNKNOWN"
if grep -q '\[stations\] TPV102 station writer active' "${DISPATCH_LOG}" 2>/dev/null; then
   STATION_VERDICT="ENGAGED — 9-station SCEC trace writer active"
else
   STATION_VERDICT="REGRESSION — no [stations] TPV102 banner; check problem.tag in TOML"
fi

{
   echo "TPV102 spatial_dyn_driver (heterogeneous Riemann, 100 m mesh, P=2 + ADER-O3, normal queue)"
   echo "Driver:                  seas_spatial_dyn_driver (TOML-driven)"
   echo "Config:                  ${CONFIG_OVERRIDE}"
   echo "Mesh:                    tpv102_100m.msh"
   echo "Polynomial order:        order 2 + ader-order 3"
   echo "Interior flux:           bimaterial (heterogeneous Riemann; byte-identical on homogeneous)"
   echo "Iterator:                one-shot (TOML default)"
   echo "Friction law:            rate_state (regularised R&S, ageing law, FL=101)"
   echo "tfinal target:           12.0 s (48 h budget may reach less; resume via --restart)"
   echo "tfinal achieved (hypo):  ${T_END} s"
   echo "hypo station rows:       ${N_ROWS}  (flt_0_7.5)"
   echo ""
   echo "peak |V_strike| at hypo  : ${PEAK_V_HYPO} m/s at t = ${PEAK_T_HYPO} s"
   echo "peak |σ_n − 120 MPa| @ h : ${SIGMA_N_PERTURB_HYPO} Pa"
   echo "slip_strike at hypo (end): ${SLIP_STRIKE_END} m"
   echo ""
   echo "Heterogeneous Riemann    : ${BIMAT_VERDICT}"
   echo "ParaView output gate     : ${PARAVIEW_VERDICT}"
   echo "SCEC station traces      : ${STATION_VERDICT}"
   echo ""
   if [ "${PEAK_V_HYPO}" != "N/A" ] \
      && echo "${BIMAT_VERDICT}"   | grep -q "^ENGAGED" \
      && echo "${PARAVIEW_VERDICT}" | grep -q "^ENGAGED" \
      && echo "${STATION_VERDICT}"  | grep -q "^ENGAGED"; then
      echo "STATUS: PASS — heterogeneous Riemann + ParaView + SCEC traces all active"
   else
      echo "STATUS: FAIL — see verdicts above"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
