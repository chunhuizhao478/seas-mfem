#!/bin/bash
#SBATCH -J bp5_build
#SBATCH -o jobs/bp5_spatial/rebuild_seas_driver_%j.out
#SBATCH -e jobs/bp5_spatial/rebuild_seas_driver_%j.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -A EAR20006

# DEDICATED, SELF-VERIFYING REBUILD JOB for ./miniapps/seas/seas_driver
# =========================================================================
# seas_driver is the BP5-QD TOML driver that the bp5_spatial jobs depend on.
# It is DISTINCT from seas_spatial_dyn_driver (the SAFS dynamic-rupture driver,
# rebuilt by safs/safs_frontera/spatial_dyn_build_safs_frontera.sbatch) — this
# job builds ONLY seas_driver.
#
# Build the driver ONCE with this job, then submit the run jobs as usual:
#   sbatch jobs/bp5_spatial/rebuild_seas_driver_frontera.sbatch   # wait for it to finish
#   sbatch jobs/bp5_spatial/bp5_verify_iterative_dev_2hr.sbatch
# The run jobs do NOT rebuild — they use the ./seas_driver this job produces.
#
# It addresses the build-issue classes this project keeps hitting:
#   (1) STALE OBJECTS -> CLEAN BUILD BY DEFAULT.  The seas Makefile has NO
#       header-dependency tracking, so an incremental `make` can silently reuse
#       an out-of-date drivers/seas_driver.o after a `git pull`.  seas_driver is
#       a SINGLE translation unit (drivers/seas_driver.o) + the parent
#       libmfem.a — all seas physics is header-only — so removing that object +
#       the binary forces the one TU to recompile against current headers.
#       extern/ and libmfem.a are NEVER touched (a clean build is ~1 min).
#       Set CLEAN=0 for a fast incremental relink when nothing relevant changed.
#   (2) toml11 / SEAS_USE_TOML.  seas_driver REQUIRES the toml11 submodule;
#       without it the binary compiles but MFEM_ABORTs on every TOML config at
#       parse time ("requires SEAS_USE_TOML").  We init the submodule and pass
#       SEAS_USE_TOML=YES explicitly.
#   (3) HYPRE 2.31.0 SHADOWING.  petsc/3.15 + mumps/5.3 prepend a pre-2.31
#       libHYPRE to LD_LIBRARY_PATH; the binary references HYPRE_Initialize
#       (renamed in 2.31).  We pin TACC_HYPRE_LIB FIRST and verify post-build
#       that no HYPRE symbol is unresolved — catching the run-time
#       "undefined symbol: HYPRE_Initialize" (MPI exit 127) abort AT BUILD TIME.
#   (4) WRONG CONFIG.  config.mk sanity (MUMPS for solver_type="mumps"; PETSc
#       for the PETSc-TS RK45 used by the verification configs).
#
# Overrides: BUILD_JOBS (parallel make jobs; default 8).
#            CLEAN=0    (default 1) -> incremental relink instead of clean build.

export LC_ALL=C
export LANG=C

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

# (3) HYPRE 2.31.0 FIRST on the runtime search path (the run jobs pin it
# identically; we do too so the post-build `ldd -r` resolves symbols the same).
export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH:-}"

# Resolve seas-mfem root by walking UP until miniapps/seas/Makefile is found
# (works regardless of the Frontera checkout path).
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
    echo "ERROR: cannot locate seas-mfem root from SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-<unset>}."
    echo "  Set it explicitly: sbatch --export=ALL,SEAS_MFEM_ROOT=/path/to/seas-mfem ..."
    exit 1
fi
cd "${SEAS_MFEM_ROOT}"

# (4) config.mk sanity.
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || { echo "ERROR: config.mk has no MFEM_USE_MUMPS=YES (solver_type=\"mumps\")"; exit 1; }
grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || echo "WARNING: config.mk has no MFEM_USE_PETSC=YES — seas_driver will warn and fall back to its own DormandPrince RK45 (use_petsc_ts ignored)."

BUILD_JOBS="${BUILD_JOBS:-8}"
CLEAN="${CLEAN:-1}"          # default 1 = clean build; CLEAN=0 = incremental relink
case "${CLEAN}" in 1|YES|yes|true|TRUE|on|ON) CLEAN_ON=1 ;; *) CLEAN_ON=0 ;; esac

cd miniapps/seas

# (2) toml11 submodule (SEAS_USE_TOML gating).
if [[ ! -f extern/toml11/toml.hpp && ! -f extern/toml11/include/toml.hpp ]]; then
    echo "toml11 submodule not initialized — initializing."
    git submodule update --init extern/toml11 || { echo "ERROR: git submodule update --init extern/toml11 failed"; exit 1; }
fi
[[ -f extern/toml11/toml.hpp || -f extern/toml11/include/toml.hpp ]] || { echo "ERROR: toml11 still missing after init"; exit 1; }

echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Build mode:  $([[ "${CLEAN_ON}" = 1 ]] && echo 'CLEAN (driver TU recompiles)' || echo 'incremental (CLEAN=0)')  (make -j${BUILD_JOBS})"
echo "In:          $(pwd)"

# (1) Clean step: seas_driver is a single TU + libmfem.a (all physics is in
# headers), so just remove the driver object + binary.
if [[ "${CLEAN_ON}" = 1 ]]; then
    echo "Clean build: removing drivers/seas_driver.o + seas_driver."
    rm -f drivers/seas_driver.o seas_driver
else
    echo "Incremental relink (CLEAN=0)."
fi

make seas_driver SEAS_USE_TOML=YES -j"${BUILD_JOBS}" || { echo "BUILD FAILED"; exit 1; }

if [[ ! -x ./seas_driver ]]; then
    echo "ERROR: make returned 0 but ./seas_driver is missing/not executable."
    exit 1
fi

# (3) Symbol/lib sanity: confirm the binary resolves HYPRE against hypre/2.31.0
# under the SAME LD_LIBRARY_PATH the run jobs use.
echo "=== link sanity (resolved libHYPRE / libpetsc / libdmumps) ==="
ldd ./seas_driver 2>/dev/null | grep -iE 'libHYPRE|libpetsc|libdmumps' || echo "  (HYPRE/PETSc/MUMPS not listed by ldd)"
hypre_undef="$(ldd -r ./seas_driver 2>&1 | grep -i 'undefined symbol.*HYPRE' || true)"
if [[ -n "${hypre_undef}" ]]; then
    echo "ERROR: HYPRE symbol(s) UNRESOLVED at build time:"
    echo "${hypre_undef}"
    echo "  -> a pre-2.31 petsc/mumps libHYPRE is shadowing hypre/2.31.0 (this is the"
    echo "     run-time 'undefined symbol: HYPRE_Initialize' / MPI exit 127 abort,"
    echo "     caught early).  TACC_HYPRE_LIB is pinned FIRST above; check 'module list'."
    exit 1
fi
echo "  HYPRE symbols resolve (HYPRE_Initialize OK)."

echo "BUILD OK: $(pwd)/seas_driver"
ls -la ./seas_driver
echo "Now submit the bp5_spatial verify/run jobs (they do a fast in-job relink, but"
echo "this clean build is the safe pre-flight after a git pull)."
