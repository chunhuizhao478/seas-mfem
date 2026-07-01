#!/bin/bash
#SBATCH -J safs_build_exp
#SBATCH -o safs_build_exp_%j.out
#SBATCH -e safs_build_exp_%j.err
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH -A lbl107

# DEDICATED build job for seas_spatial_dyn_driver on EXPANSE (SDSC).
# =========================================================================
# Recompiles + links ./miniapps/seas/seas_spatial_dyn_driver from the CURRENT
# source, then sanity-checks the binary.  This is the fast, incremental path:
# you do NOT need to re-run build_expanse.sh after a SEAS source change.
#
# WHEN TO USE WHICH:
#   * bash build_expanse.sh    -> FIRST time on Expanse, or when MFEM itself /
#       the deps (hypre/PETSc/MUMPS/HDF5/ZFP/gmsh) changed.  Builds everything
#       (deps + libmfem.a + config/config.mk + all seas binaries).  Run it in an
#       INTERACTIVE shared session (SDSC OOM-kills gcc on login nodes).
#   * THIS sbatch              -> after editing SEAS code (e.g. the OOBPolicy::
#       Clamp change in io/spatial/code).  Rebuilds ONLY the driver against the
#       already-built libmfem.a + config.mk.  A clean driver rebuild is a few
#       minutes.  Equivalent one-off interactive form:
#         srun --partition=shared --account=<acct> --nodes=1 --ntasks=1 \
#              --cpus-per-task=16 --mem=64G -t 01:00:00 --pty bash
#         module purge && module load slurm cpu/0.17.3b gcc/10.2.0 openmpi/4.1.5
#         cd <repo>/miniapps/seas && rm -f ./*.o drivers/*.o dynamic/*.o \
#              spatial/code/*.o io/*.o config/*.o friction/*.o
#         make seas_spatial_dyn_driver -j16
#
# CLEAN BY DEFAULT (class-1 build bug): the seas Makefile has NO header-
# dependency tracking, so an incremental `make` can silently reuse a stale .o
# after editing a shared header (e.g. data_field_3d.hpp / spatial_friction.hpp).
# This job removes the driver's objects first so EVERY translation unit
# recompiles.  Set CLEAN=0 only when you KNOW nothing relevant changed.
# (We do NOT `make clean` — that target is incomplete; we rm the object dirs.)
# extern/ + the parent libmfem.a are NEVER touched, so it stays fast.
#
# NOTE: unlike Frontera, there is NO HYPRE-2.31 shadowing dance — Expanse's
# hypre/PETSc/MUMPS are STATIC from-source (extern/*/install), linked INTO the
# binary, so no LD_LIBRARY_PATH ordering can break them (build_expanse.sh §hypre).
#
# Submit (chain run jobs after a verified build):
#   bid=$(sbatch --parsable miniapps/seas/jobs/safs/safs_expanse/spatial_dyn_build_safs_expanse.sbatch)
#   sbatch --dependency=afterok:${bid} miniapps/seas/jobs/safs/safs_expanse/spatial_dyn_slipweakening_v3pref_cvm_mixedflux_rk_safv4deep_expanse_2N100r_compute_2hr_safs.sbatch
#   # ... etc.  --dependency=afterok means a FAILED build cancels the run jobs.
#
# Edit -A above or submit with: sbatch --account=<your_expanse_alloc> ...
# Overrides: BUILD_JOBS (default = --cpus-per-task), CLEAN=0 (default 1).

export LC_ALL=C
export LANG=C

# --- Expanse modules (MUST match build_expanse.sh) ---
module purge
module load slurm 2>/dev/null || true
module load cpu/0.17.3b
module load gcc/10.2.0
for _mpi in openmpi/4.1.5 openmpi/4.1.3 openmpi/4.1.1 openmpi; do
    if module load "${_mpi}" 2>/dev/null && module is-loaded "${_mpi}" 2>/dev/null; then
        echo "OpenMPI module: ${_mpi}"; break
    fi
done

# Resolve seas-mfem root by walking UP until miniapps/seas/Makefile is found.
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
_walk_up_dir="${SEAS_MFEM_ROOT}"
for _ in 0 1 2 3 4 5 6 7; do
    if [[ -f "${_walk_up_dir}/miniapps/seas/Makefile" ]]; then
        SEAS_MFEM_ROOT="$(cd "${_walk_up_dir}" && pwd)"; break
    fi
    _walk_up_dir="${_walk_up_dir}/.."
done
if [[ ! -f "${SEAS_MFEM_ROOT}/miniapps/seas/Makefile" ]]; then
    echo "ERROR: cannot locate seas-mfem root from SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-<unset>}."
    echo "  Set SEAS_MFEM_ROOT via sbatch --export=ALL,SEAS_MFEM_ROOT=/path/to/seas-mfem ..."
    exit 1
fi
cd "${SEAS_MFEM_ROOT}"

# Runtime shared libs (so the post-build ldd check resolves the from-source
# HDF5/ZFP; the LINK itself uses the -L/-rpath baked into config.mk, but ldd -r
# needs them on the path).  PETSc/hypre/metis are static (no .so).
export LD_LIBRARY_PATH="${SEAS_MFEM_ROOT}/extern/hdf5/install/lib:${SEAS_MFEM_ROOT}/extern/zfp/install/lib64:${SEAS_MFEM_ROOT}/extern/zfp/install/lib:${SEAS_MFEM_ROOT}/extern/petsc/install/lib:${LD_LIBRARY_PATH:-}"

# Pre-flight: this job rebuilds ONLY the driver; libmfem.a + config.mk must
# already exist (from a prior build_expanse.sh run).
if [[ ! -f "${SEAS_MFEM_ROOT}/libmfem.a" ]]; then
    echo "ERROR: ${SEAS_MFEM_ROOT}/libmfem.a not found."
    echo "  Run the FULL build first (interactive shared session):  bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
    exit 1
fi
# config.mk must carry the build options the SAFS run jobs depend on.
grep -q '^MFEM_USE_PETSC *= YES'   config/config.mk || { echo "ERROR: config.mk has no MFEM_USE_PETSC=YES (re-run build_expanse.sh)";  exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'   config/config.mk || { echo "ERROR: config.mk has no MFEM_USE_MUMPS=YES";  exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'    config/config.mk || { echo "ERROR: config.mk has no MFEM_USE_HDF5=YES (SAFS fault HDF5 output)";   exit 1; }
grep -q '^MFEM_USE_H5Z_ZFP *= YES' config/config.mk || { echo "ERROR: config.mk has no MFEM_USE_H5Z_ZFP=YES (SAFS ZFP output)"; exit 1; }

BUILD_JOBS="${BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-16}}"
CLEAN="${CLEAN:-1}"
case "${CLEAN}" in 1|YES|yes|true|TRUE|on|ON) CLEAN_ON=1 ;; *) CLEAN_ON=0 ;; esac

cd miniapps/seas
echo "Cluster:     Expanse (SDSC)"
echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Build mode:  $([[ "${CLEAN_ON}" = 1 ]] && echo 'CLEAN (every driver TU recompiles)' || echo 'incremental (CLEAN=0)')  (make -j${BUILD_JOBS})"
echo "In:          $(pwd)"

# Clean step: remove every object that links into seas_spatial_dyn_driver so
# each recompiles from current source (works around the incomplete Makefile
# clean target; scoped to the seas tree — extern/ + libmfem.a untouched).
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

# Link sanity: confirm the binary resolves its shared deps (HDF5/ZFP) under the
# same LD_LIBRARY_PATH the run jobs use, and that no symbol is unresolved.
echo "=== link sanity (resolved shared libs) ==="
ldd ./seas_spatial_dyn_driver 2>/dev/null | grep -iE 'libhdf5|libzfp|libmpi|libopen' || echo "  (no matching shared libs listed by ldd)"
undef="$(ldd -r ./seas_spatial_dyn_driver 2>&1 | grep -i 'undefined symbol' || true)"
if [[ -n "${undef}" ]]; then
    echo "ERROR: UNRESOLVED symbol(s) at build time:"
    echo "${undef}"
    echo "  -> check the module set (module list) + that LD_LIBRARY_PATH includes"
    echo "     extern/{hdf5,zfp}/install/lib (set above).  If a dep is missing,"
    echo "     re-run build_expanse.sh."
    exit 1
fi
echo "  no unresolved symbols."

echo "BUILD OK: $(pwd)/seas_spatial_dyn_driver"
ls -la ./seas_spatial_dyn_driver
echo "Now submit the run jobs (they do NOT rebuild)."
