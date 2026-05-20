#!/bin/bash
#SBATCH -J tpv205_spdyn_100m_p1_O2_dev
#SBATCH -o tpv205_spdyn_100m_p1_O2_dev_%j.out
#SBATCH -e tpv205_spdyn_100m_p1_O2_dev_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV205 (SCEC TPV5) 100 m mesh, DEV-QUEUE smoke run, heterogeneous Riemann.
# =======================================================================
#
# Same driver/config as the 200 m dev sibling
# (tpv205_spatial_dyn_200m_p1_O2_dev.sbatch) but on the refined
# tpv205/mesh/tpv2053d_100m.geo.  Dev-queue 8 N x 400 r / 2 h envelope.
#
# Wall time
# ---------
# 100 m mesh + P=1 + cfl_safety=dg + ADER-O2 → auto dt ≈ 1.17e-4 s (CFL
# halves vs 200 m).  At 400 ranks the per-step cost (~6 s) means the 2 h
# budget reaches only ~0.1 s simulated — enough to confirm the
# heterogeneous-Riemann path engages and the nucleation patch starts to
# unzip at the finer resolution.  Full-run output is the companion
# `100m_p1_O2_normal.sbatch`.

export LC_ALL=C
export LANG=C

# R-510: fail fast on module-load / build / pipe failures.
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
# NOTE: do NOT load phdf5 on Frontera.  See 200m dev sbatch for rationale.

HDF5_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/hdf5/install"
H5Z_ZFP_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install"

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${HDF5_BUILT_DIR}/lib:${LD_LIBRARY_PATH:-}"

if [ -d "${H5Z_ZFP_BUILT_DIR}/plugin" ]; then
   export HDF5_PLUGIN_PATH="${H5Z_ZFP_BUILT_DIR}/plugin"
else
   echo "WARNING: ${H5Z_ZFP_BUILT_DIR}/plugin not found — run build_frontera.sh first."
fi

cd "${SEAS_MFEM_ROOT}"

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES'    config/config.mk || { echo "ERROR: no PETSc";   exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'    config/config.mk || { echo "ERROR: no MUMPS";   exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'     config/config.mk || { echo "ERROR: no HDF5";    exit 1; }
grep -q '^MFEM_USE_H5Z_ZFP *= YES'  config/config.mk || { echo "ERROR: no H5Z-ZFP (tpv205.toml paraview_*_zfp_tol requires it)"; exit 1; }

cd miniapps/seas

# Regenerate the 100 m .msh from .geo on first run or when FORCE_MESHGEN=1.
MESH_SRC="tpv205/mesh/tpv2053d_100m.geo"
MESH_OUT="tpv205/mesh/tpv2053d_100m.msh"
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

# Build a 100m-override config from the canonical tpv205.toml by sed-
# replacing the mesh path.  Keeps the canonical file authoritative; the
# override file lives in the run's RESULT_DIR for provenance.
RESULT_DIR="tpv205/results_spatial_dyn_100m_p1_O2_dev_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"
CONFIG_OVERRIDE="${RESULT_DIR}/tpv205_100m.toml"
sed 's|tpv2053d_200m\.msh|tpv2053d_100m.msh|g' tpv205/configs/tpv205.toml \
   > "${CONFIG_OVERRIDE}"

rm -f seas_spatial_dyn_driver drivers/spatial_dyn_driver.o \
      dynamic/tpv205_substep_iterator.o dynamic/tpv102_substep_iterator.o \
      dynamic/tpv104_substep_iterator.o dynamic/wave_operator.o \
      dynamic/godunov_flux_bimaterial.o dynamic/fault_face_flux.o \
      dynamic/heterogeneous_material.o spatial/code/spatial_friction.o \
      spatial/code/spatial_stress.o dynamic/spatial_nucleation.o
make seas_spatial_dyn_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

echo "========================================"
echo "Git commit:    $(git rev-parse HEAD)"
echo "Git branch:    $(git rev-parse --abbrev-ref HEAD)"
echo "Driver:        seas_spatial_dyn_driver (TOML-driven)"
echo "Config:        ${CONFIG_OVERRIDE} (mesh-path-overridden tpv205.toml)"
echo "Mesh:          tpv205/mesh/tpv2053d_100m.msh"
echo "Polynomial:    order 1 (TOML mesh.order) + ader-order 2 (Pelties O2)"
echo "Interior flux: bimaterial (heterogeneous Riemann — branch purpose)"
echo "Friction law:  slip_weakening (SCEC TPV5 LSW closed-form)"
echo "Iterator:      one-shot (O=1 quadrature)"
echo "Queue:         development (8N x 400r, 2h)"
echo "tfinal:        12.0 s (dev-queue 2 h SLURM-kill expected long before this)"
echo "ParaView:      fault=hdf5 every 0.001 s (TOML default; SCEC trace cadence)"
echo "========================================"

ibrun ./seas_spatial_dyn_driver \
      --config "${CONFIG_OVERRIDE}" \
      --output-dir "${RESULT_DIR}" \
      --verify-dispatch

# -----------------------------------------------------------------------
# Post-run summary — same template as the other tpv205_spatial scripts.
# -----------------------------------------------------------------------
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv205_x2_0_x3_7.5.dat"

PEAK_V_HYPO="N/A"
PEAK_T_HYPO="N/A"
SIGMA_N_PERTURB_HYPO="N/A"
SLIP_STRIKE_END="N/A"
T_END="N/A"
N_ROWS="N/A"

if [ -s "${HYPO_FILE}" ]; then
   read PEAK_V_HYPO PEAK_T_HYPO < <(awk '!/^#/ && NF>=3 {
      v=$3; if (v<0) v=-v;
      if (v>m) { m=v; tm=$1 }
   } END { printf "%g %g\n", m, tm }' "${HYPO_FILE}")
   SIGMA_N_PERTURB_HYPO=$(awk '!/^#/ && NF>=8 && $1+0==$1 {
      d=$8 - 120e6; if (d<0) d=-d;
      if (d>m) m=d
   } END { printf "%.4e", m+0 }' "${HYPO_FILE}")
   SLIP_STRIKE_END=$(awk '!/^#/ && NF>=2 && $1+0==$1 {a=$2} END{printf "%.4e", a+0}' "${HYPO_FILE}")
   T_END=$(awk '!/^#/ && NF>=1 && $1+0==$1 {tt=$1} END{printf "%.4f", tt+0}' "${HYPO_FILE}")
   N_ROWS=$(awk '!/^#/ && NF>=1 && $1+0==$1' "${HYPO_FILE}" | wc -l | tr -d ' ')
fi

DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv205_spdyn_100m_p1_O2_dev_${SLURM_JOB_ID}.out"

BIMAT_VERDICT="UNKNOWN"
if grep -q '\[wave_operator\] BimaterialFlux precomputation' "${DISPATCH_LOG}" 2>/dev/null; then
   BIMAT_VERDICT="ENGAGED — heterogeneous Riemann path active (R-001)"
elif grep -q '\[wave\] interior_flux=scalar' "${DISPATCH_LOG}" 2>/dev/null; then
   BIMAT_VERDICT="REGRESSION — scalar ctor selected"
else
   BIMAT_VERDICT="UNKNOWN — no interior-flux banner"
fi

PARAVIEW_VERDICT="UNKNOWN"
if [ -f "${RESULT_DIR}/fault.vtkhdf" ] || ls "${RESULT_DIR}"/*.pvd 1>/dev/null 2>&1; then
   PARAVIEW_VERDICT="ENGAGED — fault output written"
elif grep -q 'ParaView output: OFF' "${DISPATCH_LOG}" 2>/dev/null; then
   PARAVIEW_VERDICT="REGRESSION — driver banner reports 'OFF'"
else
   PARAVIEW_VERDICT="UNKNOWN"
fi

STATION_VERDICT="UNKNOWN"
if grep -q '\[stations\] TPV205 station writer active' "${DISPATCH_LOG}" 2>/dev/null; then
   STATION_VERDICT="ENGAGED — SCEC trace writer active"
else
   STATION_VERDICT="REGRESSION — no station banner"
fi

{
   echo "TPV205 spatial_dyn_driver (heterogeneous Riemann, 100 m mesh, P=1 + ADER-O2, dev queue)"
   echo "Driver:                  seas_spatial_dyn_driver"
   echo "Config:                  ${CONFIG_OVERRIDE}"
   echo "Mesh:                    tpv2053d_100m.msh"
   echo "Polynomial order:        order 1 + ader-order 2"
   echo "Interior flux:           bimaterial (heterogeneous Riemann)"
   echo "Iterator:                one-shot"
   echo "Friction law:            slip_weakening (SCEC TPV5 LSW)"
   echo "tfinal target:           12.0 s (dev-queue 2 h SLURM-kill expected long before)"
   echo "tfinal achieved (hypo):  ${T_END} s"
   echo "hypo station rows:       ${N_ROWS}"
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
      && echo "${BIMAT_VERDICT}"  | grep -q "^ENGAGED" \
      && echo "${PARAVIEW_VERDICT}" | grep -q "^ENGAGED" \
      && echo "${STATION_VERDICT}"  | grep -q "^ENGAGED"; then
      echo "STATUS: PASS"
   else
      echo "STATUS: FAIL — see verdicts above"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
