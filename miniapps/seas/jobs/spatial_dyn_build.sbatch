#!/bin/bash
#SBATCH -J spdyn_build
#SBATCH -o spdyn_build_%j.out
#SBATCH -e spdyn_build_%j.err
#SBATCH -p flex
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -A EAR20006

# SINGLE DEDICATED BUILD JOB for seas_spatial_dyn_driver.
# =========================================================================
# Every spatial run-job folder drives the SAME binary
# ./miniapps/seas/seas_spatial_dyn_driver:
#     jobs/tpv31_spatial/   jobs/tpv102_spatial/   jobs/tpv104_spatial/
#     jobs/tpv205_spatial/  jobs/safs/
# so it only needs to be built ONCE.  This is that one build job; compile it
# here (or on a login/idev node) BEFORE submitting any run job.
#
# WHY a single build job: the run sbatch scripts do NOT build in-job (they
# only check the binary exists).  If every run job ran `make` on the shared
# ./seas_spatial_dyn_driver, concurrent jobs would race on the same relink —
# `make` hits ETXTBSY (binary is being executed by another job) or a run
# `ibrun`s a half-written binary.  Building exactly once, here, removes that
# race; the run jobs only check the binary exists.
#
# RECOMMENDED submission (chain so run jobs wait for a clean build):
#   bid=$(sbatch --parsable miniapps/seas/jobs/spatial_dyn_build.sbatch)
#   # then any run job from any spatial folder waits on the build:
#   sbatch --dependency=afterok:${bid} miniapps/seas/jobs/tpv102_spatial/tpv102_p1_rk45_mixedflux_8N_400r_dev.sbatch
#   sbatch --dependency=afterok:${bid} miniapps/seas/jobs/tpv104_spatial/tpv104_p1_rk45_mixedflux_8N_400r_dev.sbatch
#   sbatch --dependency=afterok:${bid} miniapps/seas/jobs/tpv205_spatial/tpv205_p1_rk45_mixedflux_8N_400r_dev.sbatch
#   sbatch --dependency=afterok:${bid} miniapps/seas/jobs/tpv31_spatial/tpv31_spatial_dyn_50m_8N_400r_dev.sbatch
#   # (safs run jobs in jobs/safs/ likewise; they share this binary.)
#
# NOTE: jobs/safs/spatial_dyn_build_safs.sbatch builds the IDENTICAL binary
# (it predates this consolidation and is kept for the safs runbook).  Use
# EITHER this general job OR that one — never two builds concurrently (that
# reintroduces the relink race this job exists to avoid).
#
# Override: BUILD_JOBS (parallel make jobs; default 48).
#           CLEAN=1   force a CLEAN rebuild of the seas objects first (default
#                     0 = incremental).  Submit with:
#                       sbatch --export=ALL,CLEAN=1 .../spatial_dyn_build.sbatch
#                     Use after a `git pull` to be CERTAIN no stale object
#                     survives: the seas Makefile has NO header-dependency
#                     tracking, so an incremental `make` can silently reuse an
#                     out-of-date .o (e.g. a spatial_friction.o predating a
#                     signature/logic change -> ABI/logic mismatch or a wrong
#                     run-time abort).  NOTE: the Makefile's own `make clean` is
#                     INCOMPLETE (clean-build removes only *.o / drivers/ /
#                     tests/, NOT dynamic/ spatial/code/ io/ — where most of the
#                     driver's objects live), so this job removes the driver's
#                     object dirs explicitly.  extern/ (HDF5/ZFP/Caliper/gmsh)
#                     is never touched, so a clean build stays a few minutes.

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
for _ in 0 1 2 3 4 5 6; do
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

# --- config.mk sanity --------------------------------------------------------
# PETSc + MUMPS are the production MFEM config (build_frontera.sh) the spatial
# driver links against; hard-required.
grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || { echo "ERROR: config.mk has no MFEM_USE_PETSC=YES"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || { echo "ERROR: config.mk has no MFEM_USE_MUMPS=YES"; exit 1; }
# HDF5 + H5Z-ZFP are needed only by the run jobs that use --paraview-fault-hdf5
# / --paraview-fault-zfp-tol (tpv31_spatial, safs).  The tpv102/104/205_spatial
# jobs use --paraview-fault-vtu and do NOT need them.  WARN (do not fail): the
# one binary still serves the VTU jobs without HDF5.
grep -q '^MFEM_USE_HDF5 *= YES'    config/config.mk || echo "WARNING: MFEM_USE_HDF5 != YES — HDF5 fault output (tpv31_spatial / safs --paraview-fault-hdf5) will be unavailable; VTU jobs (tpv102/104/205) are fine."
grep -q '^MFEM_USE_H5Z_ZFP *= YES' config/config.mk || echo "WARNING: MFEM_USE_H5Z_ZFP != YES — ZFP-compressed fault output (--paraview-fault-zfp-tol) will be unavailable; VTU jobs are fine."

BUILD_JOBS="${BUILD_JOBS:-48}"
CLEAN="${CLEAN:-0}"

cd miniapps/seas
echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Building:    seas_spatial_dyn_driver  (make -j${BUILD_JOBS})"
echo "In:          $(pwd)"

# --- optional clean step (CLEAN=1) ---------------------------------------
# Remove every object that links into seas_spatial_dyn_driver so each one
# recompiles from the current source.  This works AROUND the incomplete
# Makefile `clean` target (which misses dynamic/ spatial/code/ io/), and is
# scoped to the seas source tree — extern/ (HDF5/ZFP/Caliper/gmsh) and the
# parent libmfem.a are NOT touched, so the clean build is still ~minutes.
case "${CLEAN}" in
    1|YES|yes|true|TRUE|on|ON)
        echo "Clean build (CLEAN=${CLEAN}): removing the driver's object files"
        echo "             (root + drivers/ dynamic/ spatial/code/ io/ config/"
        echo "             friction/) so every translation unit recompiles."
        rm -f ./*.o drivers/*.o dynamic/*.o spatial/code/*.o io/*.o \
              config/*.o friction/*.o
        ;;
    *)
        echo "Incremental build (CLEAN=${CLEAN}).  Re-submit with"
        echo "             'sbatch --export=ALL,CLEAN=1 ...' to force a clean"
        echo "             rebuild (recommended after a git pull)."
        ;;
esac

make seas_spatial_dyn_driver -j"${BUILD_JOBS}" || { echo "BUILD FAILED"; exit 1; }

if [[ ! -x ./seas_spatial_dyn_driver ]]; then
    echo "ERROR: make returned 0 but ./seas_spatial_dyn_driver is missing/not executable."
    exit 1
fi

echo "BUILD OK: $(pwd)/seas_spatial_dyn_driver"
ls -la ./seas_spatial_dyn_driver
echo "Now submit the run jobs from any spatial folder (they do NOT rebuild)."
