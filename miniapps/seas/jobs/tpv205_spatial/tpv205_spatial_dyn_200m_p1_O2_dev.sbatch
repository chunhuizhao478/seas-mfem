#!/bin/bash
#SBATCH -J tpv205_spdyn_p1_O2_dev
#SBATCH -o tpv205_spdyn_p1_O2_dev_%j.out
#SBATCH -e tpv205_spdyn_p1_O2_dev_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV205 (SCEC TPV5) DEV-QUEUE smoke run through the heterogeneous
# Riemann solver via seas_spatial_dyn_driver.
# =======================================================================
#
# Branch: feature/heterogeneous_riemann_solver
# Driver: drivers/spatial_dyn_driver.cpp (NEW config-driven driver — NOT
#         the byte-parity native drivers/tpv205_driver.cpp).
# Config: tpv205/configs/tpv205.toml (TOML-driven; CLI overrides only for
#         output_dir / mesh path / restart prefix).
#
# Code-fix audit (REVIEW.md 2026-05-19) addressed before this run:
#   R-001  tpv205.toml no longer sets `interior_flux = "scalar"` → defaults
#          to "bimaterial" → WaveOperator(MaterialField, BoundaryConfig)
#          ctor → BimaterialFlux::ApplyPerFaceFlux on every interior face.
#          The banner must print `[wave_operator] BimaterialFlux
#          precomputation: ...` confirming the heterogeneous path engaged.
#   R-004  tpv205.toml sets `paraview_enabled = true` → master gate open
#          → fault.vtkhdf written under <output_dir>/.
#   R-005  spatial_dyn_driver aborts loudly on zero fault faces (instead
#          of silent useless run).
#   R-006  spatial_dyn_driver emits a clear `gmsh -format msh22` command
#          if the .msh path is missing instead of MFEM file-open abort.
#   R-007  TPV205 SCEC station traces are written under <output_dir>/
#          (16 stations from DefaultStations_TPV205, prefix=tpv205_).
#   R-009  TPV205 barrier rules use epsilon-nudged bounds so DOFs at
#          x=±15 km / z=−15 km stay INSIDE the rupture area (matches
#          native InRuptureArea_TPV205 `<=` semantic).
#
# Spec encoding (verbatim per tpv205/configs/tpv205.toml + TPV205Params):
#   ρ = 2670 kg/m³, c_p = 6000 m/s, c_s = 3464 m/s
#   σ_n = 120 MPa, τ_back = 70 MPa, μ_s = 0.677 / 10000 barrier, μ_d = 0.525, d_c = 0.40 m
#   3 static stress patches at depth 7.5 km, half-side 1.5 km:
#     central (0, 0, -7500) σ_xy = 81.6 MPa  (yield-exceeding ⇒ spontaneous nucleation)
#     left    (-7500, 0, -7500) σ_xy = 78.0 MPa
#     right   (+7500, 0, -7500) σ_xy = 62.0 MPa
#   Rupture initiates at t = 0+ in the central patch (μ_s·σ_n = 81.24 MPa < 81.6 MPa).
#
# Wall-time budget
# ----------------
# Dev queue 2 h envelope.  At 200 m mesh + P=1 + cfl_safety=dg + ADER-O2
# the auto dt ≈ 2.34e-4 s and ~51 k macro-steps to reach tfinal = 12 s.
# Per-step cost is dominated by BimaterialFlux precomputed-tensor matvec
# on ~3.7 M interior faces; expect ~5–10 s/step on 400 ranks, so
# 2 h ≈ 720–1440 steps → tfinal_actual ≈ 0.17–0.34 s.  Enough to see the
# nucleation patch unzip and the early rupture front; full-run output is
# the companion `_normal.sbatch` queue.

export LC_ALL=C
export LANG=C

# R-510: fail fast on module-load / build / pipe failures (mirrors
# jobs/tpv205/tpv205_phase6_paraview_zfp_dev_2hr.sbatch).
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
# Last-resort fallback if the probe walked off the top of the filesystem.
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
# intel/23.1.0 + impi/21.9.0 (incompatible with the intel/19 hypre/mumps/
# petsc stack we use), and older modules (1.10/1.12) lack symbols MFEM
# mesh/vtkhdf.cpp needs (H5S_BLOCK, the unsigned-long-long hsize_t
# typedef).  build_frontera.sh builds HDF5 1.14.x from source into
# ${SEAS_MFEM_ROOT}/extern/hdf5/install instead; we link/run against that.

HDF5_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/hdf5/install"
H5Z_ZFP_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install"

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${HDF5_BUILT_DIR}/lib:${LD_LIBRARY_PATH:-}"

# REQUIRED for tpv205.toml's `paraview_fault_zfp_tol = 1e-12`: the HDF5
# ZFP filter plugin must be loadable at run time.  If this path doesn't
# exist, ParaView fault output will fall back to writing uncompressed
# (large) data — keep the warning but don't abort here; the driver
# itself emits a clearer error on the first write attempt.
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
# R-404: fail before ibrun if H5Z-ZFP is missing — tpv205.toml's
# `paraview_*_zfp_tol` fields cause the driver to abort at parse time
# without it, wasting the queue slot.
grep -q '^MFEM_USE_H5Z_ZFP *= YES'  config/config.mk || { echo "ERROR: no H5Z-ZFP (tpv205.toml paraview_*_zfp_tol requires it)"; exit 1; }

cd miniapps/seas

# Regenerate the .msh from .geo on first run or when FORCE_MESHGEN=1.
# The 200 m TPV205 mesh is too large to commit (~125 MB); the .geo
# source is tracked, gmsh is installed by build_frontera.sh.
MESH_SRC="tpv205/mesh/tpv2053d_200m.geo"
MESH_OUT="tpv205/mesh/tpv2053d_200m.msh"
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

RESULT_DIR="tpv205/results_spatial_dyn_p1_O2_dev_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

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
echo "Config:        tpv205/configs/tpv205.toml"
echo "Mesh:          tpv205/mesh/tpv2053d_200m.msh"
echo "Polynomial:    --order 1 + --ader-order 2 (Pelties O2)"
echo "Interior flux: bimaterial (heterogeneous Riemann — branch purpose)"
echo "Friction law:  slip_weakening (SCEC TPV5 LSW closed-form)"
echo "Iterator:      one-shot (O=1 quadrature)"
echo "Queue:         development (8N x 400r, 2h)"
echo "tfinal:        12.0 s (dev-queue 2 h SLURM-kill expected long before this)"
echo "ParaView:      fault=hdf5 every 0.001 s (TOML default; SCEC trace cadence)"
echo "========================================"

# CLI overrides: output_dir is the only field we must override (the TOML
# default 'tpv205/out' would collide with prior runs).  All physics/
# numerics come from the TOML.  The `--paraview` flag is NOT required
# here because tpv205.toml sets `paraview_enabled = true` directly (R-004).
ibrun ./seas_spatial_dyn_driver \
      --config tpv205/configs/tpv205.toml \
      --output-dir "${RESULT_DIR}" \
      --verify-dispatch

# -----------------------------------------------------------------------
# Post-run summary.  The spatial driver writes SCEC trace files under
# tpv205_<station>.dat by default (TPV205StationWriter prefix=tpv205).
# Hypocenter station name follows the SCEC convention `x2_0_x3_7.5`.
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
   # Column layout (tpv205_setup.hpp::TPV205StationWriter::WriteStep):
   #   1: time   2: h-slip(strike=slip2)   3: h-slip-rate(V2)
   #   4: h-shear-stress(tau2_corr)        5: v-slip(dip=slip1)
   #   6: v-slip-rate(V1)                  7: v-shear-stress(tau1_corr)
   #   8: n-stress(sigma_n_corr)           9: mu_eff(δ)
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

DISPATCH_LOG="${SBATCH_LOG_DIR}/tpv205_spdyn_p1_O2_dev_${SLURM_JOB_ID}.out"

# R-001 verdict: must see the heterogeneous Riemann precomputation banner.
BIMAT_VERDICT="UNKNOWN"
if grep -q '\[wave_operator\] BimaterialFlux precomputation' "${DISPATCH_LOG}" 2>/dev/null; then
   BIMAT_VERDICT="ENGAGED — heterogeneous Riemann path active (R-001)"
elif grep -q '\[wave\] interior_flux=scalar' "${DISPATCH_LOG}" 2>/dev/null; then
   BIMAT_VERDICT="REGRESSION — scalar ctor selected; check tpv205.toml [numerics].interior_flux"
else
   BIMAT_VERDICT="UNKNOWN — no interior-flux banner in log"
fi

# R-004 verdict: ParaView fault output must be present.
PARAVIEW_VERDICT="UNKNOWN"
if [ -f "${RESULT_DIR}/fault.vtkhdf" ] || ls "${RESULT_DIR}"/*.pvd 1>/dev/null 2>&1; then
   PARAVIEW_VERDICT="ENGAGED — fault output written"
elif grep -q 'ParaView output: OFF' "${DISPATCH_LOG}" 2>/dev/null; then
   PARAVIEW_VERDICT="REGRESSION — driver banner reports 'OFF'; check tpv205.toml [output].paraview_enabled"
else
   PARAVIEW_VERDICT="UNKNOWN — no fault file and no 'OFF' banner"
fi

# R-007 verdict: TPV205 SCEC traces must be written.
STATION_VERDICT="UNKNOWN"
if grep -q '\[stations\] TPV205 station writer active' "${DISPATCH_LOG}" 2>/dev/null; then
   STATION_VERDICT="ENGAGED — SCEC trace writer active"
else
   STATION_VERDICT="REGRESSION — no [stations] TPV205 banner; check problem.tag in TOML"
fi

{
   echo "TPV205 spatial_dyn_driver (heterogeneous Riemann, P=1 + ADER-O2, dev queue)"
   echo "Driver:                  seas_spatial_dyn_driver (TOML-driven)"
   echo "Config:                  tpv205/configs/tpv205.toml"
   echo "Mesh:                    tpv2053d_200m.msh"
   echo "Polynomial order:        --order 1 (TOML mesh.order) + --ader-order 2"
   echo "Interior flux:           bimaterial (heterogeneous Riemann)"
   echo "Iterator:                one-shot (TOML default)"
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
      echo "STATUS: PASS — heterogeneous Riemann + ParaView + SCEC traces all active"
   else
      echo "STATUS: FAIL — see verdicts above"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
