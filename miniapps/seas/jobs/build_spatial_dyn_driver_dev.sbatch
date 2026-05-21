#!/bin/bash
#SBATCH -J build_spatial_dyn_driver
#SBATCH -o build_spatial_dyn_driver_%j.out
#SBATCH -e build_spatial_dyn_driver_%j.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:40:00
#SBATCH -A EAR20006

# One-shot builder for seas_spatial_dyn_driver (development queue, 1 node).
# =======================================================================
#
# All spatial-dyn RUN jobs — tpv102_spatial/, tpv104_spatial/,
# tpv205_spatial/, tpv31_spatial/, safs/ — run the SAME binary
# `seas_spatial_dyn_driver`.  The benchmark is selected by the `--config`
# TOML at run time, NOT at build time, so this ONE build job produces the
# binary used by every spatial-dyn run job.  The run jobs themselves do
# NOT build anything.
#
# Why a separate build job
# ------------------------
# Previously every run job rebuilt the driver into the shared source tree.
# When two ran at once they raced on the shared object files and one
# truncated drivers/spatial_dyn_driver.o mid-write, giving the linker
# error "file not recognized: file truncated".  Centralising the build
# here removes the race: build once, then submit any number of run jobs
# concurrently (they only read the finished binary).
#
# Workflow
# --------
#   1. sbatch build_spatial_dyn_driver_dev.sbatch        # this job
#   2. wait until it reaches COMPLETED  (squeue -u $USER)
#   3. submit run jobs (dev and/or normal), as many concurrently as you
#      like — none of them rebuild.
#
# Re-run THIS job (only) after any source change; it is the single place
# that compiles the driver.  Do NOT submit two copies of this job at the
# same time — that would reintroduce the concurrent-build race.

export LC_ALL=C
export LANG=C

# Fail fast on module-load / build / pipe failures.
set -euo pipefail

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
# ${SEAS_MFEM_ROOT}/extern/hdf5/install instead; we link against that.

HDF5_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/hdf5/install"

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${HDF5_BUILT_DIR}/lib:${LD_LIBRARY_PATH:-}"

cd "${SEAS_MFEM_ROOT}"

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES'    config/config.mk || { echo "ERROR: no PETSc";   exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'    config/config.mk || { echo "ERROR: no MUMPS";   exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'     config/config.mk || { echo "ERROR: no HDF5 (run-job paraview_fault=hdf5 needs it)"; exit 1; }
# The TPV*/SAFS run-job TOMLs set paraview_*_zfp_tol, which makes the
# driver abort at parse time without H5Z-ZFP — fail here so a run job
# never starts against a driver that can't honour its config.
grep -q '^MFEM_USE_H5Z_ZFP *= YES'  config/config.mk || { echo "ERROR: no H5Z-ZFP (run-job TOML paraview_*_zfp_tol requires it)"; exit 1; }

cd miniapps/seas

# Clean rebuild of the spatial-dyn driver and its objects (the same object
# set the run jobs used to remove before building).  Running this job
# alone means there is no other `make` writing these files concurrently.
rm -f seas_spatial_dyn_driver drivers/spatial_dyn_driver.o \
      dynamic/tpv205_substep_iterator.o dynamic/tpv102_substep_iterator.o \
      dynamic/tpv104_substep_iterator.o dynamic/wave_operator.o \
      dynamic/godunov_flux_bimaterial.o dynamic/fault_face_flux.o \
      dynamic/heterogeneous_material.o spatial/code/spatial_friction.o \
      spatial/code/spatial_stress.o dynamic/spatial_nucleation.o
make seas_spatial_dyn_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

if [ ! -x seas_spatial_dyn_driver ]; then
   echo "ERROR: seas_spatial_dyn_driver missing after make."
   exit 1
fi

echo "========================================"
echo "BUILD OK"
echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
ls -l "${SEAS_MFEM_ROOT}/miniapps/seas/seas_spatial_dyn_driver"
echo "Binary ready — submit the spatial-dyn run jobs now."
echo "========================================"
