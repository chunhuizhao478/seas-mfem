#!/bin/bash
#SBATCH -J tpv102_p6_zfp
#SBATCH -o tpv102_p6_zfp_%j.out
#SBATCH -e tpv102_p6_zfp_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 06:00:00
#SBATCH -A EAR20006

# TPV102 — Phase 6 ParaView compaction reference run.
# =======================================================================
# Goal: produce a "post-Phase-6" Frontera reference for the size-estimator
# calibration table (REFERENCE_RUNS.md).  Exercises the new flag triple:
#   --paraview-volume-zfp-tol 1e-3   (primary collection: velocity)
#   --paraview-bulk-zfp-tol   1e-3   (secondary collection: stresses)
#   --paraview-fault-zfp-tol  1e-12  (fault: slip-rate floor preserved)
#
# Phase 6 default-flip note: WITHOUT --paraview-volume-vtu the primary
# volume collection is now ParaViewHDFDataCollection (single .vtkhdf
# per run).  Older sbatch scripts that hard-coded
# <output>/ParaView/<basename>_<rank>_<cycle>.vtu paths will break
# silently; switch to .vtkhdf or pass --paraview-volume-vtu to revert.
#
# Mesh / scale: tpv102_200m.msh, P=1, ADER-O=2.  At 400 ranks ~5.5 s/step
# the 2 h SLURM envelope yields ~1080 steps ≈ 5 s simulated.  R-412:
# `--tfinal 5.0` lets the run complete cleanly and flush HDF5 inside the
# wall-time instead of being SLURM-killed mid-step (VTKHDF restart-mode
# is not yet wired, so a kill could leave a partial volume.vtkhdf).

export LC_ALL=C
export LANG=C

# R-510: fail fast.  `set -e` aborts on any unhandled non-zero exit
# (e.g. failed `module load`); `-u` flags unset variables; `pipefail`
# propagates errors through the post-run `du | sort` pipe.
set -euo pipefail

SBATCH_LOG_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

# Project checkout root.  Override at submit time when running from a
# different clone, e.g. `SEAS_MFEM_ROOT=$(pwd) sbatch ...`.
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-/scratch2/10024/zhaochun/seas-project/seas-mfem}"

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8
# NOTE: do NOT load phdf5 on Frontera.  The phdf5/1.14.x module
# requires intel/23.1.0 + impi/21.9.0 (incompatible with the intel/19
# hypre/mumps/petsc stack we use), and older modules (1.10/1.12) lack
# symbols MFEM mesh/vtkhdf.cpp needs (H5S_BLOCK, the unsigned-long-long
# hsize_t typedef).  build_frontera.sh builds HDF5 1.14.x from source
# into ${SEAS_MFEM_ROOT}/extern/hdf5/install instead; we link/run
# against that.

HDF5_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/hdf5/install"
H5Z_ZFP_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install"

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${HDF5_BUILT_DIR}/lib:${LD_LIBRARY_PATH:-}"

# H5Z-ZFP runtime plugin discovery — required for --paraview-*-zfp-tol.
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
# R-404: fail before ibrun if H5Z-ZFP is missing — the driver aborts at
# parse time on --paraview-*-zfp-tol otherwise, wasting ~30 SU.
grep -q '^MFEM_USE_H5Z_ZFP *= YES'  config/config.mk || { echo "ERROR: no H5Z-ZFP (--paraview-{volume,bulk,fault}-zfp-tol requires it)"; exit 1; }

cd miniapps/seas

# Regenerate the .msh from .geo on first run or when FORCE_MESHGEN=1.
# The 200 m TPV meshes are too large to commit (~125 MB each); the .geo
# source is tracked, gmsh is installed by build_frontera.sh.
MESH_SRC="tpv102/mesh/tpv102_200m.geo"
MESH_OUT="tpv102/mesh/tpv102_200m.msh"
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
    "${GMSH_BIN}" -3 "${MESH_SRC}" -o "${MESH_OUT}" \
        || { echo "MESH GENERATION FAILED"; exit 1; }
fi

if [ ! -s "${MESH_OUT}" ]; then
   echo "ERROR: ${MESH_OUT} missing after mesh-regen step."
   exit 1
fi

RESULT_DIR="tpv102/results_phase6_zfp_dev_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

rm -f seas_tpv102_driver drivers/tpv102_driver.o
make seas_tpv102_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

echo "========================================"
echo "Git commit:    $(git rev-parse HEAD)"
echo "Git branch:    $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:    ${RESULT_DIR}"
echo "Mesh:          tpv102_200m.msh (P=1 + ADER-O2)"
echo "Volume PV:     HDF5 (default), ZFP tol 1e-3 on velocity"
echo "Bulk  PV:      HDF5,           ZFP tol 1e-3 on sigma_yy/xy/xz"
echo "Fault PV:      HDF5,           ZFP tol 1e-12 (slip-rate floor)"
echo "Queue:         development (8N x 400r, 2h)"
echo "========================================"

ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 12.0 \
      --ader-order 2 \
      --friction-solver newton-stable \
      --fric-law aging \
      --fault-iterator one-shot \
      --mixed-flux adjacent \
      --output-dir "${RESULT_DIR}" \
      --output-prefix tpv102_p6_zfp \
      --output-dt 0.05 \
      --paraview \
      --paraview-dt 0.5 \
      --paraview-bulk-dt 0.5 \
      --paraview-volume-zfp-tol 1e-3 \
      --paraview-bulk-zfp-tol 1e-3 \
      --paraview-fault-zfp-tol 1e-12 \
      --paraview-max-snapshots 200 \
      --verify-dispatch

# Post-run size summary — for REFERENCE_RUNS.md.
RESULT_FILE="${RESULT_DIR}/SIZE.txt"
{
   echo "TPV102 Phase 6 ParaView ZFP — size summary"
   echo "Mesh:   tpv102_200m.msh"
   echo "Output: ${RESULT_DIR}"
   echo
   echo "du -sb (total):"
   du -sb "${RESULT_DIR}"
   echo
   echo "du -sh per top-level entry:"
   du -sh "${RESULT_DIR}"/* 2>/dev/null | sort -h
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
