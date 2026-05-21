#!/bin/bash
#SBATCH -J tpv104_spdyn_1000m_p1_O2_dev
#SBATCH -o tpv104_spdyn_1000m_p1_O2_dev_%j.out
#SBATCH -e tpv104_spdyn_1000m_p1_O2_dev_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV104 (SCEC FL=103 slip-law strong-rate-weakening) 1000 m mesh,
# P=1 + ADER-O2, DEV-QUEUE smoke run through the heterogeneous Riemann
# solver.
# =======================================================================
#
# Same driver/mesh/config as tpv104_spatial_dyn_1000m_p2_O3_normal.sbatch
# but at polynomial order 1 (P1) with ADER time order 2 (Pelties O2), on
# the development queue with a 2 h envelope.  Use this to confirm the run
# starts cleanly (build OK, mesh OK, bimaterial banner, station + ParaView
# output) BEFORE committing a 48 h normal-queue allocation.
#
# How P1/O2 is selected
# ---------------------
# No --order CLI flag exists; the mesh polynomial order is set ONLY via
# the TOML [mesh].order field, which is already `order = 1` in the
# canonical tpv104.toml — so this run consumes the config VERBATIM (no
# sed override, unlike the P2/O3 sbatch).  ADER order IS a CLI flag →
# `--ader-order 2` (matches tpv104.toml [numerics].ader_order; stated
# explicitly for provenance).
#
# Why the heterogeneous Riemann path
# ----------------------------------
# tpv104.toml omits [numerics].interior_flux → defaults to "bimaterial",
# so the heterogeneous ctor is used and the banner must print
# `[wave_operator] BimaterialFlux precomputation: ...`.  The material is
# homogeneous (rho=2670, cs=3464, cp=6000), so the cached per-face
# GodunovFlux is built from the exact (λ,μ,ρ) and the result is byte-
# identical to the scalar path (CLAUDE.md "Phase H Stage 2").
#
# Spec encoding: identical to the p2/O3 normal sbatch (see its header).
# In brief: homogeneous material; slip-law strong-rate-weakening (FL=103,
# f0=0.6/V0=1e-6/b=0.014/Dc=0.4/f_w=0.2, C∞ tanh boxcar a_in=0.01/a_out=
# 0.02 + V_w_in=0.1/V_w_out=1.0); σ_n=120 MPa, τ_ini=40 MPa right-lateral;
# SCEC compact-bell nucleation (Δτ₀=45 MPa, R=3 km, T=1 s); vertical y=0
# fault, hypocenter (0,0,−7500).
#
# Wall-time budget
# ----------------
# Dev queue 2 h envelope.  At P1/O2 on the 1000 m mesh the auto dt is
# ~9e-3 s (cfl_safety=dg) and tfinal = 12 s needs O(10^3) macro-steps; a
# 2 h dev slot at 400 ranks should reach a meaningful fraction of (or all
# of) the run.  The full-fidelity run is the companion
# `_p2_O3_normal.sbatch`.

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

# Regenerate the 1000 m .msh from .geo on first run or when FORCE_MESHGEN=1.
MESH_SRC="tpv104/mesh/tpv104_1000m.geo"
MESH_OUT="tpv104/mesh/tpv104_1000m.msh"
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

RESULT_DIR="tpv104/results_spatial_dyn_1000m_p1_O2_dev_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

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
echo "Config:        tpv104/configs/tpv104.toml"
echo "Mesh:          tpv104/mesh/tpv104_1000m.msh"
echo "Polynomial:    --order 1 (TOML mesh.order) + --ader-order 2 (Pelties O2)"
echo "Interior flux: bimaterial (heterogeneous Riemann; byte-identical on homogeneous material)"
echo "Friction law:  rate_state (slip-law strong-rate-weakening, FL=103)"
echo "Iterator:      one-shot (O=1 quadrature)"
echo "Queue:         development (8N x 400r, 2h)"
echo "tfinal:        12.0 s (dev-queue 2 h SLURM-kill may occur before this)"
echo "Stations:      9 SCEC on-fault traces (prefix tpv104_)"
echo "ParaView:      fault=hdf5 every 0.001 s (TOML default)"
echo "========================================"

# CLI overrides: output_dir (avoid clobbering tpv104/out) + explicit
# --ader-order 2 (matches tpv104.toml [numerics].ader_order; stated for
# provenance/symmetry with the p2/O3 run).  All physics from the TOML;
# P1 needs no config override (order = 1 is the canonical TOML default).
ibrun ./seas_spatial_dyn_driver \
      --config tpv104/configs/tpv104.toml \
      --output-dir "${RESULT_DIR}" \
      --ader-order 2 \
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

DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv104_spdyn_1000m_p1_O2_dev_${SLURM_JOB_ID}.out"

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
   echo "TPV104 spatial_dyn_driver (heterogeneous Riemann, 1000 m mesh, P=1 + ADER-O2, dev queue)"
   echo "Driver:                  seas_spatial_dyn_driver (TOML-driven)"
   echo "Config:                  tpv104/configs/tpv104.toml"
   echo "Mesh:                    tpv104_1000m.msh"
   echo "Polynomial order:        --order 1 (TOML mesh.order) + --ader-order 2"
   echo "Interior flux:           bimaterial (heterogeneous Riemann; byte-identical on homogeneous)"
   echo "Iterator:                one-shot (TOML default)"
   echo "Friction law:            rate_state (slip-law strong-rate-weakening, FL=103)"
   echo "tfinal target:           12.0 s (dev-queue 2 h SLURM-kill may occur before)"
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
