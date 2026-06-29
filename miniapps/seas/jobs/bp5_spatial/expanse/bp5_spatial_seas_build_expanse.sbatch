#!/bin/bash
#SBATCH -J bp5seas_build_exp
#SBATCH -o bp5seas_build_exp_%j.out
#SBATCH -e bp5seas_build_exp_%j.err
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH -A ddp408

# DEDICATED build job for seas_spatial_seas_driver on EXPANSE (SDSC).
# =========================================================================
# Recompiles + links ./miniapps/seas/seas_spatial_seas_driver (the generic
# QUASI-DYNAMIC SEAS driver — BP5-spatial + SAFS-QD) from the CURRENT source,
# then sanity-checks the binary.  This is the fast, incremental path: you do NOT
# need to re-run build_expanse.sh after a SEAS source change.
#
# It is the jobs/bp5_spatial/expanse sibling of
#   jobs/safs/safs_expanse/spatial_dyn_build_safs_expanse.sbatch
# (which builds the seas_spatial_DYN_driver).  This one builds the QD driver
# that the jobs/bp5_spatial/expanse Caliper run job profiles.
#
# WHEN TO USE WHICH:
#   * bash build_expanse.sh    -> FIRST time on Expanse, or when MFEM itself /
#       the deps (hypre/PETSc/MUMPS/HDF5/ZFP/gmsh/Caliper) changed.  Builds
#       everything (deps + libmfem.a + config/config.mk + all seas binaries,
#       INCLUDING seas_spatial_seas_driver).  Run it in an INTERACTIVE shared
#       session (SDSC OOM-kills gcc on login nodes).  For a CALIPER-profiled
#       binary you MUST do that full build with USE_CALIPER=YES:
#         USE_CALIPER=YES bash build_expanse.sh
#   * THIS sbatch              -> after editing SEAS code (e.g. the QD driver,
#       spatial_friction, or a stress/velocity sidecar).  Rebuilds ONLY the QD
#       driver against the already-built libmfem.a + config.mk.  A clean driver
#       rebuild is a few minutes.  Equivalent one-off interactive form:
#         srun --partition=shared --account=<acct> --nodes=1 --ntasks=1 \
#              --cpus-per-task=16 --mem=64G -t 01:00:00 --pty bash
#         module purge && module load slurm cpu/0.17.3b gcc/10.2.0 openmpi/4.1.5
#         cd <repo>/miniapps/seas && rm -f drivers/spatial_seas_driver.o \
#              spatial/code/*.o io/*.o dynamic/heterogeneous_material.o
#         make seas_spatial_seas_driver SEAS_USE_TOML=YES -j16
#
# CLEAN BY DEFAULT (class-1 build bug): the seas Makefile has NO header-
# dependency tracking, so an incremental `make` can silently reuse a stale .o
# after editing a shared header (e.g. spatial_friction.hpp / data_field_3d.hpp).
# This job removes EVERY object that links into seas_spatial_seas_driver first so
# all of them recompile.  Set CLEAN=0 only when you KNOW nothing relevant
# changed.  (We do NOT `make clean` — that target is incomplete; we rm the
# driver's specific objects.)  extern/ + the parent libmfem.a are NEVER touched.
#
# NOTE: seas_spatial_seas_driver REQUIRES the toml11 submodule — it calls
# LoadSpatialFrictionConfig, gated on SEAS_USE_TOML, which auto-resolves to YES
# from extern/toml11/toml.hpp.  We init the submodule + pass SEAS_USE_TOML=YES
# explicitly so a missing submodule fails loudly instead of building a binary
# that MFEM_ABORTs on every --config at parse time.
#
# Unlike Frontera, there is NO HYPRE-2.31 shadowing dance — Expanse's
# hypre/PETSc/MUMPS are STATIC from-source (extern/*/install), linked INTO the
# binary, so no LD_LIBRARY_PATH ordering can break them (build_expanse.sh §hypre).
#
# Submit (chain the run job after a verified build):
#   bid=$(sbatch --parsable miniapps/seas/jobs/bp5_spatial/expanse/bp5_spatial_seas_build_expanse.sbatch)
#   sbatch --dependency=afterok:${bid} miniapps/seas/jobs/bp5_spatial/expanse/bp5_spatial_seas_cgamg_caliper_expanse_2N80r_100steps.sbatch
#   # --dependency=afterok means a FAILED build cancels the run job.
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
# HDF5/ZFP + Caliper; the LINK itself uses the -L/-rpath baked into config.mk,
# but `ldd -r` needs them on the path).  PETSc/hypre/metis are static (no .so).
export LD_LIBRARY_PATH="${SEAS_MFEM_ROOT}/extern/caliper/install/lib64:${SEAS_MFEM_ROOT}/extern/caliper/install/lib:${SEAS_MFEM_ROOT}/extern/hdf5/install/lib:${SEAS_MFEM_ROOT}/extern/zfp/install/lib64:${SEAS_MFEM_ROOT}/extern/zfp/install/lib:${SEAS_MFEM_ROOT}/extern/petsc/install/lib:${LD_LIBRARY_PATH:-}"

# Pre-flight: this job rebuilds ONLY the driver; libmfem.a + config.mk must
# already exist (from a prior build_expanse.sh run).
if [[ ! -f "${SEAS_MFEM_ROOT}/libmfem.a" ]]; then
    echo "ERROR: ${SEAS_MFEM_ROOT}/libmfem.a not found."
    echo "  Run the FULL build first (interactive shared session):  bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
    echo "  (for a Caliper-profiled binary:  USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh)"
    exit 1
fi
# config.mk must carry the build options the QD driver depends on.  PETSc gives
# the iterative KSP stack; MUMPS is needed only for [solver].type="mumps" (the
# CG-AMG run job does not need it) — warn, don't fail, when it is absent.
grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || { echo "ERROR: config.mk has no MFEM_USE_PETSC=YES (re-run build_expanse.sh)"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || echo "WARNING: config.mk has no MFEM_USE_MUMPS=YES — [solver].type=\"mumps\" configs will fail; the CG-AMG run job is unaffected."
if grep -q '^MFEM_USE_CALIPER *= YES' config/config.mk; then
    echo "config.mk: MFEM_USE_CALIPER=YES — driver will be Caliper-linked (profiling enabled)."
else
    echo "WARNING: config.mk has MFEM_USE_CALIPER=NO — the rebuilt binary will NOT be"
    echo "         Caliper-linked, so the Caliper run job's per-region/sample report is"
    echo "         EMPTY.  For profiling, rebuild libmfem.a first with:"
    echo "             USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
fi

BUILD_JOBS="${BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-16}}"
CLEAN="${CLEAN:-1}"
case "${CLEAN}" in 1|YES|yes|true|TRUE|on|ON) CLEAN_ON=1 ;; *) CLEAN_ON=0 ;; esac

cd miniapps/seas

# toml11 submodule (SEAS_USE_TOML gating) — seas_spatial_seas_driver needs it.
if [[ ! -f extern/toml11/toml.hpp && ! -f extern/toml11/include/toml.hpp ]]; then
    echo "toml11 submodule not initialized — initializing."
    git submodule update --init extern/toml11 || { echo "ERROR: git submodule update --init extern/toml11 failed"; exit 1; }
fi
[[ -f extern/toml11/toml.hpp || -f extern/toml11/include/toml.hpp ]] || { echo "ERROR: toml11 still missing after init"; exit 1; }

echo "Cluster:     Expanse (SDSC)"
echo "Git commit:  $(git rev-parse HEAD)"
echo "Git branch:  $(git rev-parse --abbrev-ref HEAD)"
echo "Build mode:  $([[ "${CLEAN_ON}" = 1 ]] && echo 'CLEAN (every driver TU recompiles)' || echo 'incremental (CLEAN=0)')  (make -j${BUILD_JOBS})"
echo "In:          $(pwd)"

# Clean step: remove every object that links into seas_spatial_seas_driver so
# each recompiles from current source.  The driver's link set (Makefile §2295):
#   drivers/spatial_seas_driver.o            (SPATIAL_SEAS_DRIVER_OBJ)
#   spatial/code/spatial_friction.o          (SPATIAL_FRICTION_OBJ)
#   spatial/code/spatial_velocity.o          (SPATIAL_VELOCITY_OBJ)
#   spatial/code/spatial_stress.o            (SPATIAL_STRESS_OBJ)
#   io/data_field_3d.o io/stress_field_3d.o  (DATA_FIELD_3D_OBJ / STRESS_FIELD_3D_OBJ)
#   io/field_coefficient.o io/material_coefficients.o (FIELD_COEFFICIENT_OBJ / MATERIAL_COEFFICIENTS_OBJ)
#   dynamic/heterogeneous_material.o         (HETEROGENEOUS_MATERIAL_OBJ)
# (rm the whole spatial/code + io dirs; only heterogeneous_material from dynamic/.)
if [[ "${CLEAN_ON}" = 1 ]]; then
    echo "Clean build: removing seas_spatial_seas_driver objects (drivers/ spatial/code/ io/ dynamic/heterogeneous_material)."
    rm -f drivers/spatial_seas_driver.o \
          spatial/code/*.o \
          io/*.o \
          dynamic/heterogeneous_material.o
else
    echo "Incremental build (CLEAN=0): reusing existing objects (use only if nothing relevant changed)."
fi

make seas_spatial_seas_driver SEAS_USE_TOML=YES -j"${BUILD_JOBS}" || { echo "BUILD FAILED"; exit 1; }

if [[ ! -x ./seas_spatial_seas_driver ]]; then
    echo "ERROR: make returned 0 but ./seas_spatial_seas_driver is missing/not executable."
    exit 1
fi

# Link sanity: confirm the binary resolves its shared deps (HDF5/ZFP/Caliper)
# under the same LD_LIBRARY_PATH the run job uses, and that no symbol is unresolved.
echo "=== link sanity (resolved shared libs) ==="
ldd ./seas_spatial_seas_driver 2>/dev/null | grep -iE 'libhdf5|libzfp|libcaliper|libmpi|libopen' || echo "  (no matching shared libs listed by ldd)"
if ldd ./seas_spatial_seas_driver 2>/dev/null | grep -qi 'libcaliper'; then
    echo "  Caliper: LINKED (profiling enabled)."
else
    echo "  Caliper: NOT linked (the run job's Caliper report would be empty — rebuild libmfem.a with USE_CALIPER=YES)."
fi
undef="$(ldd -r ./seas_spatial_seas_driver 2>&1 | grep -i 'undefined symbol' || true)"
if [[ -n "${undef}" ]]; then
    echo "ERROR: UNRESOLVED symbol(s) at build time:"
    echo "${undef}"
    echo "  -> check the module set (module list) + that LD_LIBRARY_PATH includes"
    echo "     extern/{hdf5,zfp,caliper}/install/lib (set above).  If a dep is missing,"
    echo "     re-run build_expanse.sh."
    exit 1
fi
echo "  no unresolved symbols."

echo "BUILD OK: $(pwd)/seas_spatial_seas_driver"
ls -la ./seas_spatial_seas_driver
echo "Now submit the run job (it does NOT rebuild):"
echo "  sbatch jobs/bp5_spatial/expanse/bp5_spatial_seas_cgamg_caliper_expanse_2N80r_100steps.sbatch"
