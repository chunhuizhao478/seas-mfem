#!/bin/bash
#SBATCH -J safs_build
#SBATCH -o safs_build_%j.out
#SBATCH -e safs_build_%j.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -A EAR20006

# DEDICATED BUILD JOB for seas_spatial_dyn_driver.
# =========================================================================
# Compiles + links ./miniapps/seas/seas_spatial_dyn_driver ONCE on a single
# node.  Submit (or run on a login/idev node) BEFORE the run jobs.
#
# WHY a separate build job: the run sbatch scripts no longer build in-job.  If
# every run job ran `make` on the shared ./seas_spatial_dyn_driver, concurrent
# jobs would race on the same relink — `make` hits ETXTBSY (binary is being
# executed by another job) or a run `ibrun`s a half-written binary.  Building
# exactly once, here, removes that race; the run jobs only check the binary
# exists.
#
# RECOMMENDED submission (chain so run jobs wait for a clean build):
#   bid=$(sbatch --parsable miniapps/seas/jobs/safs/spatial_dyn_build_safs.sbatch)
#   sbatch --dependency=afterok:${bid} miniapps/seas/jobs/safs/spatial_dyn_slipweakening_normalcap_10N_500r_normal_48hr_safs.sbatch
#   sbatch --dependency=afterok:${bid} miniapps/seas/jobs/safs/spatial_dyn_ratestate_depthprofile_normalcap_10N_500r_normal_48hr_safs.sbatch
#   # ... etc. for the dev / triq run jobs.
#
# Override: BUILD_JOBS (parallel make jobs; default 48).

export LC_ALL=C
export LANG=C

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH}"

# Resolve seas-mfem root by walking UP until miniapps/seas/Makefile is found.
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-${SLURM_SUBMIT_DIR}}"
_walk_up_dir="${SEAS_MFEM_ROOT}"
for _ in 0 1 2 3 4 5; do
    if [[ -f "${_walk_up_dir}/miniapps/seas/Makefile" ]]; then
        SEAS_MFEM_ROOT="$(cd "${_walk_up_dir}" && pwd)"
        break
    fi
    _walk_up_dir="${_walk_up_dir}/.."
done
if [[ ! -f "${SEAS_MFEM_ROOT}/miniapps/seas/Makefile" ]]; then
    echo "ERROR: cannot locate seas-mfem root from SLURM_SUBMIT_DIR="
    echo "  ${SLURM_SUBMIT_DIR}.  Set SEAS_MFEM_ROOT explicitly via"
    echo "  sbatch --export=ALL,SEAS_MFEM_ROOT=/path/to/seas-mfem ..."
    exit 1
fi
cd "${SEAS_MFEM_ROOT}"

# config.mk must have the same build options the run jobs depend on (the ZFP
# RUNTIME plugin is checked by the run jobs, not here — it is not needed to
# build).
grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'   config/config.mk || { echo "ERROR: no HDF5"; exit 1; }
grep -q '^MFEM_USE_H5Z_ZFP *= YES' config/config.mk || { echo "ERROR: no H5Z_ZFP"; exit 1; }

BUILD_JOBS="${BUILD_JOBS:-48}"

cd miniapps/seas
echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Building:    seas_spatial_dyn_driver seas_partition_mesh  (make -j${BUILD_JOBS})"
echo "In:          $(pwd)"

# seas_partition_mesh is the pre-split tool (R-OOM): it partitions a large
# serial .msh into per-rank files so production with --mesh-mode presplit avoids
# the serial-mesh-replication startup OOM.  Cheap to build alongside the driver.
make seas_spatial_dyn_driver seas_partition_mesh -j"${BUILD_JOBS}" || { echo "BUILD FAILED"; exit 1; }

if [[ ! -x ./seas_spatial_dyn_driver ]]; then
    echo "ERROR: make returned 0 but ./seas_spatial_dyn_driver is missing/not executable."
    exit 1
fi
if [[ ! -x ./seas_partition_mesh ]]; then
    echo "ERROR: make returned 0 but ./seas_partition_mesh is missing/not executable."
    exit 1
fi

echo "BUILD OK: $(pwd)/seas_spatial_dyn_driver"
ls -la ./seas_spatial_dyn_driver
echo "Now submit the run jobs (they do NOT rebuild)."
