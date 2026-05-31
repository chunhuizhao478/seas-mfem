#!/bin/bash
#SBATCH -J safs_build
#SBATCH -o safs_build_%j.out
#SBATCH -e safs_build_%j.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -A EAR20006

# DEDICATED, SELF-VERIFYING BUILD JOB for seas_spatial_dyn_driver (SAFS).
# =========================================================================
# Compiles + links ./miniapps/seas/seas_spatial_dyn_driver ONCE on a single
# node, then SANITY-CHECKS the binary.  Submit (or run on a login/idev node)
# BEFORE the run jobs.  This is the "handle every build issue" build — it is
# deliberately conservative because the SAFS run jobs are expensive
# (10N/500r/48h normal-queue), so a wrong binary is costly.  It addresses the
# three classes of build problem this project keeps hitting:
#
#   (1) STALE OBJECTS  ->  CLEAN BUILD BY DEFAULT.
#       The seas Makefile has NO header-dependency tracking, so an incremental
#       `make` can silently reuse an out-of-date .o after a `git pull` (the
#       V_w side-channel, the ZFP/banner, the print-derived gate, ... — every
#       one of which bit a run job this cycle).  This job removes the driver's
#       objects first so EVERY translation unit recompiles.  Set CLEAN=0 only
#       when you KNOW nothing relevant changed and want a fast relink.
#       NB: we do NOT call `make clean` — its clean-build target is INCOMPLETE
#       (removes only *.o / drivers/ / tests/, NOT dynamic/ spatial/code/ io/,
#       where most of the driver's objects live).  We rm the object dirs
#       explicitly.  extern/ (HDF5/ZFP/Caliper/gmsh) and the parent libmfem.a
#       are NEVER touched, so a clean build is still only a few minutes.
#
#   (2) HYPRE 2.31.0 SHADOWING  ->  PIN + VERIFY.
#       petsc/3.15 + mumps/5.3 prepend a PRE-2.31 libHYPRE to LD_LIBRARY_PATH;
#       the driver references HYPRE_Initialize (renamed in hypre 2.31).  We put
#       TACC_HYPRE_LIB FIRST (exactly as the run jobs do) AND verify post-build
#       that no HYPRE symbol is unresolved — catching the run-time
#       "undefined symbol: HYPRE_Initialize" (MPI exit 127) abort AT BUILD TIME.
#
#   (3) WRONG CONFIG  ->  config.mk sanity + binary existence checks.
#
# WHY a separate build job: the run sbatch scripts do NOT build in-job.  If
# every run job ran `make` on the shared binary, concurrent jobs would race on
# the relink (`make` ETXTBSY / a run `ibrun`s a half-written binary).  Building
# exactly once, here, removes that race; run jobs only check the binary exists.
#
# RECOMMENDED submission (chain so run jobs wait for a clean, verified build):
#   bid=$(sbatch --parsable miniapps/seas/jobs/safs/spatial_dyn_build_safs.sbatch)
#   sbatch --dependency=afterok:${bid} miniapps/seas/jobs/safs/spatial_dyn_slipweakening_normalcap_10N_500r_normal_48hr_safs.sbatch
#   sbatch --dependency=afterok:${bid} miniapps/seas/jobs/safs/spatial_dyn_ratestate_depthprofile_normalcap_10N_500r_normal_48hr_safs.sbatch
#   # ... etc. for the dev / triq run jobs.  --dependency=afterok means a FAILED
#   # build (incl. the symbol check below) cancels the run jobs automatically.
#
# Overrides: BUILD_JOBS  (parallel make jobs; default 48).
#            CLEAN=0     (default 1) -> incremental relink instead of a clean build.

export LC_ALL=C
export LANG=C

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

# (2) HYPRE 2.31.0 FIRST on the runtime search path — see header.  The run
# jobs pin it identically; we do too so the post-build `ldd -r` check below
# resolves symbols the same way the run will.
export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH:-}"

# Resolve seas-mfem root by walking UP until miniapps/seas/Makefile is found.
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
_walk_up_dir="${SEAS_MFEM_ROOT}"
for _ in 0 1 2 3 4 5 6 7; do
    if [[ -f "${_walk_up_dir}/miniapps/seas/Makefile" ]]; then
        SEAS_MFEM_ROOT="$(cd "${_walk_up_dir}" && pwd)"
        break
    fi
    _walk_up_dir="${_walk_up_dir}/.."
done
if [[ ! -f "${SEAS_MFEM_ROOT}/miniapps/seas/Makefile" ]]; then
    echo "ERROR: cannot locate seas-mfem root from SLURM_SUBMIT_DIR="
    echo "  ${SLURM_SUBMIT_DIR:-<unset>}.  Set SEAS_MFEM_ROOT explicitly via"
    echo "  sbatch --export=ALL,SEAS_MFEM_ROOT=/path/to/seas-mfem ..."
    exit 1
fi
cd "${SEAS_MFEM_ROOT}"

# (3) config.mk must carry the build options the SAFS run jobs depend on.
grep -q '^MFEM_USE_PETSC *= YES'   config/config.mk || { echo "ERROR: config.mk has no MFEM_USE_PETSC=YES";  exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'   config/config.mk || { echo "ERROR: config.mk has no MFEM_USE_MUMPS=YES";  exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'    config/config.mk || { echo "ERROR: config.mk has no MFEM_USE_HDF5=YES (SAFS uses --paraview-fault-hdf5)";   exit 1; }
grep -q '^MFEM_USE_H5Z_ZFP *= YES' config/config.mk || { echo "ERROR: config.mk has no MFEM_USE_H5Z_ZFP=YES (SAFS uses --paraview-fault-zfp-tol)"; exit 1; }

BUILD_JOBS="${BUILD_JOBS:-48}"
CLEAN="${CLEAN:-1}"           # default 1 = clean build; CLEAN=0 = incremental
case "${CLEAN}" in 1|YES|yes|true|TRUE|on|ON) CLEAN_ON=1 ;; *) CLEAN_ON=0 ;; esac

cd miniapps/seas
echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Build mode:  $([[ "${CLEAN_ON}" = 1 ]] && echo 'CLEAN (every TU recompiles)' || echo 'incremental (CLEAN=0)')  (make -j${BUILD_JOBS})"
echo "In:          $(pwd)"

# (1) Clean step: remove every object that links into seas_spatial_dyn_driver
# so each recompiles from the current source.  Works around the incomplete
# Makefile clean target; scoped to the seas tree (extern/ + libmfem.a untouched).
if [[ "${CLEAN_ON}" = 1 ]]; then
    echo "Clean build: removing driver objects (root + drivers/ dynamic/ spatial/code/ io/ config/ friction/)."
    rm -f ./*.o drivers/*.o dynamic/*.o spatial/code/*.o io/*.o config/*.o friction/*.o
else
    echo "Incremental build (CLEAN=0): reusing existing objects (use only if nothing relevant changed)."
fi

make seas_spatial_dyn_driver -j"${BUILD_JOBS}" || { echo "BUILD FAILED"; exit 1; }

if [[ ! -x ./seas_spatial_dyn_driver ]]; then
    echo "ERROR: make returned 0 but ./seas_spatial_dyn_driver is missing/not executable."
    exit 1
fi

# (2) Symbol/lib sanity: confirm the binary resolves HYPRE against hypre/2.31.0
# under the SAME LD_LIBRARY_PATH the run jobs use, catching the run-time
# "undefined symbol: HYPRE_Initialize" abort AT BUILD TIME.
echo "=== link sanity (resolved libHYPRE / libpetsc) ==="
ldd ./seas_spatial_dyn_driver 2>/dev/null | grep -iE 'libHYPRE|libpetsc' || echo "  (HYPRE/PETSc not listed by ldd)"
hypre_undef="$(ldd -r ./seas_spatial_dyn_driver 2>&1 | grep -i 'undefined symbol.*HYPRE' || true)"
if [[ -n "${hypre_undef}" ]]; then
    echo "ERROR: HYPRE symbol(s) UNRESOLVED at build time:"
    echo "${hypre_undef}"
    echo "  -> a pre-2.31 petsc/mumps libHYPRE is shadowing hypre/2.31.0 (this is"
    echo "     the run-time 'undefined symbol: HYPRE_Initialize' / MPI exit 127"
    echo "     abort, caught early).  TACC_HYPRE_LIB is pinned FIRST above; check"
    echo "     the loaded module set (module list) and that hypre/2.31.0 loaded."
    exit 1
fi
echo "  HYPRE symbols resolve (HYPRE_Initialize OK)."

echo "BUILD OK: $(pwd)/seas_spatial_dyn_driver"
ls -la ./seas_spatial_dyn_driver
echo "Now submit the run jobs (they do NOT rebuild)."
