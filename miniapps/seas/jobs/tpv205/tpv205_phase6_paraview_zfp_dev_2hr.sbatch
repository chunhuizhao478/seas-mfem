#!/bin/bash
#SBATCH -J tpv205_p6_zfp
#SBATCH -o tpv205_p6_zfp_%j.out
#SBATCH -e tpv205_p6_zfp_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV205 — Phase 6 ParaView compaction reference run.
# =======================================================================
# Companion to tpv102/tpv104 Phase 6 sbatch.  TPV205 is the SCEC slip-
# weakening (LSW) benchmark with non-planar/branching nucleation patch.

export LC_ALL=C
export LANG=C

# R-510: fail fast on module-load / build / pipe failures.
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

if [ ! -s tpv205/mesh/tpv205_200m.msh ]; then
   echo "ERROR: tpv205/mesh/tpv205_200m.msh missing."
   exit 1
fi

RESULT_DIR="tpv205/results_phase6_zfp_dev_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

rm -f seas_tpv205_driver drivers/tpv205_driver.o
make seas_tpv205_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

echo "========================================"
echo "Git commit:    $(git rev-parse HEAD)"
echo "Git branch:    $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:    ${RESULT_DIR}"
echo "Mesh:          tpv205_200m.msh (P=1 + ADER-O2)"
echo "Volume PV:     HDF5 (default), ZFP tol 1e-3 on velocity"
echo "Bulk  PV:      HDF5,           ZFP tol 1e-3 on sigma_yy/xy/xz"
echo "Fault PV:      HDF5,           ZFP tol 1e-12 (slip-rate floor)"
echo "Queue:         development (8N x 400r, 2h)"
echo "========================================"

ibrun ./seas_tpv205_driver \
      --mesh tpv205/mesh/tpv205_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 5.0 \
      --ader-order 2 \
      --friction-solver newton-stable \
      --fric-law lsw \
      --fault-iterator one-shot \
      --mixed-flux adjacent \
      --output-dir "${RESULT_DIR}" \
      --output-prefix tpv205_p6_zfp \
      --output-dt 0.05 \
      --paraview \
      --paraview-dt 0.5 \
      --paraview-bulk-dt 0.5 \
      --paraview-volume-zfp-tol 1e-3 \
      --paraview-bulk-zfp-tol 1e-3 \
      --paraview-fault-zfp-tol 1e-12 \
      --paraview-max-snapshots 200 \
      --verify-dispatch

RESULT_FILE="${RESULT_DIR}/SIZE.txt"
{
   echo "TPV205 Phase 6 ParaView ZFP — size summary"
   echo "Mesh:   tpv205_200m.msh"
   echo "Output: ${RESULT_DIR}"
   echo
   echo "du -sb (total):"
   du -sb "${RESULT_DIR}"
   echo
   echo "du -sh per top-level entry:"
   du -sh "${RESULT_DIR}"/* 2>/dev/null | sort -h
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
