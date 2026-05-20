#!/bin/bash
#SBATCH -J tpv31_spdyn_50m_p2_O3_norm
#SBATCH -o tpv31_spdyn_50m_p2_O3_norm_%j.out
#SBATCH -e tpv31_spdyn_50m_p2_O3_norm_%j.err
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 500
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# SCEC TPV31 — 50 m mesh, P=2 + ADER-O3, NORMAL-QUEUE production run,
# heterogeneous Riemann.
# =======================================================================
#
# Same driver/mesh/config as tpv31_spatial_dyn_50m_p1_O2_dev.sbatch but
# at polynomial order 2 (P2) with ADER time order 3 (Pelties O3), on the
# normal queue with a 48 h budget.
#
# How P2/O3 is selected
# ---------------------
# No --order CLI flag exists; the mesh polynomial order is set ONLY via
# the TOML [mesh].order field.  We build a sed-override config that flips
# `order = 1` -> `order = 2` (override file under RESULT_DIR for
# provenance).  ADER order IS a CLI flag → `--ader-order 3` (validated in
# {2,3,4}).
#
# Why TPV31 needs the heterogeneous Riemann path
# ----------------------------------------------
# TPV31's discontinuous 1D velocity structure (depth_profile_1d material,
# MaterialField::Mode::Coefficient) is mutually exclusive with the scalar
# WaveOperator ctor.  tpv31.toml omits [numerics].interior_flux → defaults
# to "bimaterial"; the banner must print
# `[wave_operator] BimaterialFlux precomputation: ...`.
#
# Spec encoding: identical to the dev sbatch (see its header).  In brief:
#   discontinuous 1D bi-material (depth_profile_1d), depth-proportional
#   pre-stress (σ_xx=60, σ_yy=60, σ_zz=0, σ_xy=30 MPa-per-µ, compression
#   POSITIVE), instantaneous circular cosine-tapered nucleation (Δτ=4.95
#   MPa·µ/µ0, r≤1400→taper→0 at 2000 m), LSW µ_s=0.580/µ_d=0.450/d_c=0.18
#   with depth-linear cohesion (1.02 MPa at surface → 0 below 2.4 km),
#   vertical y=0 right-lateral fault reaching the free surface.
#
# Rank count / wall time
# ----------------------
# 10 N x 500 r.  The BimaterialFlux precomputation footprint is per-FACE
# (order-independent), so memory headroom is unchanged from P1.  P2 raises
# per-step cost ~3-5x and tightens the CFL, so the 48 h budget reaches
# noticeably less simulated time than the P1/O2 run.  Resume via
# `--restart <prefix>` for follow-on blocks (checkpoint every 10000 steps
# per tpv31.toml [output].checkpoint_every_steps), or scale nodes/walltime
# if the full 15 s tfinal is required in one shot.

export LC_ALL=C
export LANG=C

set -euo pipefail

SBATCH_LOG_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

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
# NOTE: do NOT load phdf5 on Frontera.  See the dev sbatch for rationale.

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
grep -q '^MFEM_USE_H5Z_ZFP *= YES'  config/config.mk || { echo "ERROR: no H5Z-ZFP (tpv31.toml paraview_*_zfp_tol requires it)"; exit 1; }

cd miniapps/seas

# Regenerate the .msh from .geo on first run or when FORCE_MESHGEN=1.
MESH_SRC="tpv31/mesh/tpv31_50m.geo"
MESH_OUT="tpv31/mesh/tpv31_50m.msh"
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
    "${GMSH_BIN}" -format msh22 -3 "${MESH_SRC}" -o "${MESH_OUT}" \
        || { echo "MESH GENERATION FAILED"; exit 1; }
fi

if [ ! -s "${MESH_OUT}" ]; then
   echo "ERROR: ${MESH_OUT} missing after mesh-regen step."
   exit 1
fi

# Build a P2-override config from the canonical tpv31.toml by sed-
# replacing `order = 1` -> `order = 2` (mesh polynomial order; there is
# no --order CLI flag).  Override file lives under RESULT_DIR.
RESULT_DIR="tpv31/results_spatial_dyn_50m_p2_O3_normal_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"
CONFIG_OVERRIDE="${RESULT_DIR}/tpv31_50m_p2.toml"
sed 's|^order = 1|order = 2|' tpv31/configs/tpv31.toml \
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
echo "Config:        ${CONFIG_OVERRIDE} (order=2-overridden tpv31.toml)"
echo "Mesh:          tpv31/mesh/tpv31_50m.msh"
echo "Polynomial:    order 2 (TOML mesh.order override) + ader-order 3 (Pelties O3)"
echo "Interior flux: bimaterial (heterogeneous Riemann — REQUIRED for TPV31)"
echo "Friction law:  slip_weakening (LSW + depth-linear cohesion)"
echo "Iterator:      one-shot (O=1 quadrature)"
echo "Queue:         normal (10N x 500r, 48h)"
echo "tfinal:        15.0 s (48 h budget may reach less; resume via --restart)"
echo "Checkpoint:    every 10000 steps (TOML default; resume via --restart)"
echo "Stations:      30 SCEC on-fault traces (prefix tpv31_; spec Part 5)"
echo "ParaView:      fault=hdf5 every 0.001 s (TOML default)"
echo "========================================"

ibrun ./seas_spatial_dyn_driver \
      --config "${CONFIG_OVERRIDE}" \
      --output-dir "${RESULT_DIR}" \
      --ader-order 3 \
      --verify-dispatch

# -----------------------------------------------------------------------
# Post-run summary — same template as the dev sbatch.  Hypocenter station
# = faultst000dp075 (0 km strike, 7.5 km down-dip).  Columns per
# dynamic/tpv31_stations.hpp::WriteStep:
#   1:t 2:h-slip 3:h-slip-rate 4:h-shear(Pa) 5:v-slip 6:v-slip-rate
#   7:v-shear 8:n-stress(Pa,compression+) 9:mu_eff
# -----------------------------------------------------------------------
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv31_station_faultst000dp075.dat"

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
   SIGMA_N_PERTURB_HYPO=$(awk '!/^#/ && NF>=8 && $1+0==$1 {
      if (n==0) { base=$8 } n++;
      d=$8 - base; if (d<0) d=-d;
      if (d>m) m=d
   } END { printf "%.4e", m+0 }' "${HYPO_FILE}")
   SLIP_STRIKE_END=$(awk '!/^#/ && NF>=2 && $1+0==$1 {a=$2} END{printf "%.4e", a+0}' "${HYPO_FILE}")
   T_END=$(awk '!/^#/ && NF>=1 && $1+0==$1 {tt=$1} END{printf "%.4f", tt+0}' "${HYPO_FILE}")
   N_ROWS=$(awk '!/^#/ && NF>=1 && $1+0==$1' "${HYPO_FILE}" | wc -l | tr -d ' ')
fi

DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv31_spdyn_50m_p2_O3_norm_${SLURM_JOB_ID}.out"

BIMAT_VERDICT="UNKNOWN"
if grep -q '\[wave_operator\] BimaterialFlux precomputation' "${DISPATCH_LOG}" 2>/dev/null; then
   BIMAT_VERDICT="ENGAGED — heterogeneous Riemann path active (REQUIRED for TPV31)"
elif grep -q '\[wave\] interior_flux=scalar' "${DISPATCH_LOG}" 2>/dev/null; then
   BIMAT_VERDICT="REGRESSION — scalar ctor selected; TPV31 needs bimaterial"
else
   BIMAT_VERDICT="UNKNOWN — no interior-flux banner in log"
fi

PARAVIEW_VERDICT="UNKNOWN"
if [ -f "${RESULT_DIR}/fault.vtkhdf" ] || ls "${RESULT_DIR}"/*.vtkhdf 1>/dev/null 2>&1 || ls "${RESULT_DIR}"/*.pvd 1>/dev/null 2>&1; then
   PARAVIEW_VERDICT="ENGAGED — fault output written"
elif grep -q 'ParaView output: OFF' "${DISPATCH_LOG}" 2>/dev/null; then
   PARAVIEW_VERDICT="REGRESSION — driver banner reports 'OFF'"
else
   PARAVIEW_VERDICT="UNKNOWN — no fault file and no 'OFF' banner"
fi

STATION_VERDICT="UNKNOWN"
if grep -q '\[stations\] TPV31 station writer active' "${DISPATCH_LOG}" 2>/dev/null; then
   STATION_VERDICT="ENGAGED — 30-station SCEC trace writer active"
else
   STATION_VERDICT="REGRESSION — no [stations] TPV31 banner; check problem.tag in TOML"
fi

{
   echo "TPV31 spatial_dyn_driver (heterogeneous Riemann, 50 m mesh, P=2 + ADER-O3, normal queue)"
   echo "Driver:                  seas_spatial_dyn_driver"
   echo "Config:                  ${CONFIG_OVERRIDE}"
   echo "Mesh:                    tpv31_50m.msh"
   echo "Polynomial order:        order 2 + ader-order 3"
   echo "Interior flux:           bimaterial (heterogeneous Riemann; REQUIRED)"
   echo "Iterator:                one-shot"
   echo "Friction law:            slip_weakening (LSW + depth-linear cohesion)"
   echo "tfinal target:           15.0 s"
   echo "tfinal achieved (hypo):  ${T_END} s"
   echo "hypo station rows:       ${N_ROWS}  (faultst000dp075)"
   echo ""
   echo "peak |V_strike| at hypo  : ${PEAK_V_HYPO} m/s at t = ${PEAK_T_HYPO} s"
   echo "peak |Δσ_n| at hypo      : ${SIGMA_N_PERTURB_HYPO} Pa (vs initial)"
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
