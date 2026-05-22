#!/bin/bash
#SBATCH -J tpv104_spdyn_200m_p3_O4_norm
#SBATCH -o tpv104_spdyn_200m_p3_O4_norm_%j.out
#SBATCH -e tpv104_spdyn_200m_p3_O4_norm_%j.err
#SBATCH -p normal
#SBATCH -N 16
#SBATCH -n 800
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# TPV104 (SCEC FL=103 slip-law strong-rate-weakening) 200 m mesh,
# P=3 + ADER-O4, NORMAL-QUEUE production run, heterogeneous Riemann.
# =======================================================================
#
# Runs TPV104 through the TOML-driven seas_spatial_dyn_driver (NOT the
# native drivers/tpv104_driver.cpp) at polynomial order 3 (P3) with ADER
# time order 4 (Pelties O4).  This is the heterogeneous-Riemann path
# verification for TPV104 (R-001): on a homogeneous Mode::Constant
# material the bimaterial interior flux is byte-identical to the scalar
# ctor, so a clean match against the native TPV104 reference is expected.
#
# How P3/O4 is selected
# ---------------------
# No --order CLI flag exists; the mesh polynomial order is set ONLY via
# the TOML [mesh].order field.  We build a sed-override config that flips
# `order = 1` -> `order = 3` (override file under RESULT_DIR for
# provenance).  ADER order IS a CLI flag → `--ader-order 4` (validated in
# {2,3,4}).  The sed-override rewrites BOTH the mesh path (1000m -> 200m) AND
# the order, exactly like the tpv205_spatial jobs.
#
# Why the heterogeneous Riemann path
# ----------------------------------
# tpv104.toml omits [numerics].interior_flux → defaults to "bimaterial",
# so the heterogeneous (MaterialField, BoundaryConfig) ctor is used and
# the banner must print `[wave_operator] BimaterialFlux precomputation: ...`.
# The material is homogeneous (rho=2670, cs=3464, cp=6000), so the cached
# per-face GodunovFlux is built from the exact (λ,μ,ρ) and the result is
# byte-identical to the scalar path (CLAUDE.md "Phase H Stage 2").
#
# Spec encoding (verbatim per tpv104/configs/tpv104.toml):
#   * Material:   homogeneous (μ=32.038119520e9, λ=32.043760960e9, ρ=2670).
#   * Friction:   slip-law strong-rate-weakening (Noda & Lapusta 2013,
#                 FL=103).  f0=0.6, V0=1e-6, b=0.014, Dc=L=0.4 m, f_w=0.2;
#                 `a` AND `V_w` via C∞ tanh boxcar (a_in=0.01/a_out=0.02,
#                 V_w_in=0.1/V_w_out=1.0 m/s) over |x|≤15 km ∧
#                 |z+7.5 km|≤7.5 km with 3 km transition.
#   * Loading:    σ_n=120 MPa (compression POSITIVE), τ_ini=40 MPa
#                 right-lateral along-strike, V_init=1e-16 m/s.
#   * Nucleation: SCEC compact bell F(r)=exp(r²/(r²−R²)) (R=3 km) ×
#                 smoothStep G(t) (T=1 s); peak Δτ₀=45 MPa, strike-only.
#   * Fault:      vertical y=0 plane; hypocenter (0,0,−7500).
#
# Rank count / wall time
# ----------------------
# Fixed at 16 N x 800 r (50 ranks/node).  The BimaterialFlux precomputation
# footprint is per-FACE (order-independent), so memory headroom is unchanged
# from P1.  At P3 the per-step cost is ~10-20x P1 and the CFL tightens
# sharply, so the 48 h budget reaches less simulated time than a P2/O3 run;
# chain 48 h blocks via `--restart <prefix>` (checkpoint every 10000 steps
# per tpv104.toml [output].checkpoint_every_steps) to advance further.

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

# REQUIRED for tpv104.toml's `paraview_fault_zfp_tol = 1e-12`: the HDF5
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
# tpv104.toml's `paraview_*_zfp_tol` fields make the driver abort at parse
# time without H5Z-ZFP — fail before ibrun so we don't waste the slot.
grep -q '^MFEM_USE_H5Z_ZFP *= YES'  config/config.mk || { echo "ERROR: no H5Z-ZFP (tpv104.toml paraview_*_zfp_tol requires it)"; exit 1; }

cd miniapps/seas

# Regenerate the 200 m .msh from .geo on first run or when FORCE_MESHGEN=1.
MESH_SRC="tpv104/mesh/tpv104_200m.geo"
MESH_OUT="tpv104/mesh/tpv104_200m.msh"
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

# Build a 200m+P3 override config from the canonical tpv104.toml by sed-
# replacing the mesh path (1000m -> 200m) AND `order = 1` -> `order = 3`
# (no --order CLI flag).  Override file lives under RESULT_DIR.
RESULT_DIR="tpv104/results_spatial_dyn_200m_p3_O4_normal_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"
CONFIG_OVERRIDE="${RESULT_DIR}/tpv104_200m_p3.toml"
sed -e 's|tpv104_1000m\.msh|tpv104_200m.msh|' \
    -e 's|^order = 1|order = 3|' tpv104/configs/tpv104.toml \
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
echo "Config:        ${CONFIG_OVERRIDE} (mesh 200m + order=3 override of tpv104.toml)"
echo "Mesh:          tpv104/mesh/tpv104_200m.msh"
echo "Polynomial:    order 3 (TOML mesh.order override) + ader-order 4 (Pelties O4)"
echo "Interior flux: bimaterial (heterogeneous Riemann; byte-identical on homogeneous material)"
echo "Friction law:  rate_state (slip-law strong-rate-weakening, FL=103)"
echo "Iterator:      one-shot (O=1 quadrature)"
echo "Queue:         normal (16N x 800r, 48h)"
echo "tfinal:        12.0 s (48 h budget may reach less; resume via --restart)"
echo "Checkpoint:    every 10000 steps (TOML default; resume via --restart)"
echo "Stations:      9 SCEC on-fault traces (prefix tpv104_)"
echo "ParaView:      fault=hdf5 every 0.001 s (TOML default)"
echo "========================================"

ibrun ./seas_spatial_dyn_driver \
      --config "${CONFIG_OVERRIDE}" \
      --output-dir "${RESULT_DIR}" \
      --ader-order 4 \
      --verify-dispatch

# -----------------------------------------------------------------------
# Post-run summary.  TPV104 writes SCEC on-fault traces as
# tpv104_station_<name>.dat (TPV104StationWriter, prefix tpv104).  The
# hypocenter station is x2_0_x3_7.5 (0 km strike, 7.5 km down-dip).
# -----------------------------------------------------------------------
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv104_station_x2_0_x3_7.5.dat"

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

DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv104_spdyn_200m_p3_O4_norm_${SLURM_JOB_ID}.out"

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
   PARAVIEW_VERDICT="REGRESSION — driver banner reports 'OFF'; check tpv104.toml [output].paraview_enabled"
else
   PARAVIEW_VERDICT="UNKNOWN — no fault file and no 'OFF' banner"
fi

# Gate 3 — SCEC station traces written.
STATION_VERDICT="UNKNOWN"
if grep -q '\[stations\] TPV104 station writer active' "${DISPATCH_LOG}" 2>/dev/null; then
   STATION_VERDICT="ENGAGED — 9-station SCEC trace writer active"
else
   STATION_VERDICT="REGRESSION — no [stations] TPV104 banner; check problem.tag in TOML"
fi

{
   echo "TPV104 spatial_dyn_driver (heterogeneous Riemann, 200 m mesh, P=3 + ADER-O4, normal queue)"
   echo "Driver:                  seas_spatial_dyn_driver (TOML-driven)"
   echo "Config:                  ${CONFIG_OVERRIDE}"
   echo "Mesh:                    tpv104_200m.msh"
   echo "Polynomial order:        order 3 + ader-order 4"
   echo "Interior flux:           bimaterial (heterogeneous Riemann; byte-identical on homogeneous)"
   echo "Iterator:                one-shot (TOML default)"
   echo "Friction law:            rate_state (slip-law strong-rate-weakening, FL=103)"
   echo "tfinal target:           12.0 s (48 h budget may reach less; resume via --restart)"
   echo "tfinal achieved (hypo):  ${T_END} s"
   echo "hypo station rows:       ${N_ROWS}  (x2_0_x3_7.5)"
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
