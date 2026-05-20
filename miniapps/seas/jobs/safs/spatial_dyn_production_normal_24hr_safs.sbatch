#!/bin/bash
#SBATCH -J spatial_dyn_prod
#SBATCH -o spatial_dyn_prod_%j.out
#SBATCH -e spatial_dyn_prod_%j.err
#SBATCH -p normal
#SBATCH -N 32
#SBATCH -n 1600
#SBATCH -t 24:00:00
#SBATCH -A EAR20006

# Phase Z production sbatch for seas_spatial_dyn_driver — first SAFS
# dynamic-rupture production run, normal queue, 24 hr wall, 12 s
# simulated time.  Single-event, gradual_overstress nucleation.
#
# Node count derivation (R-006 fix — done OFFLINE, not from this script):
#   1. Run `seas_spatial_dyn_driver --dry-run --print-derived` on one
#      Frontera dev node with the production TOML + mesh.
#   2. Read num_fault_global and num_global_elements from the log.
#   3. nodes = ceil(num_global_elements / (50 ranks/node · 2500 elem/rank)).
#      For the SAFS 1000m zgraded mesh (~4M tets) this gives ~32 nodes.
#      Bake the result into this header.  The 2500-elem/rank target is
#      empirical (BP5/TPV104 scaling); revisit if the first production
#      run is dominated by MPI time.
#
# Resource budget (SAFS 1000m lcfar3000, constant material,
# regional H&Z stress, mixed-flux "adjacent", O2 ADER, tfinal=12s):
#   Mesh:       ~4M tets, ~120M DOFs (3D velocity-stress, order 1)
#   Per-rank:   ~2.5 GB peak (wave Q + DOFData + GodunovFlux)
#   Wall:       ~18-22 hr at cfl=0.5, dt~6e-5 s, ~200k steps
#   Output:     fault.vtkhdf ~80 GB with ZFP tol 1e-12
#               cp files every 10k steps × ~5 MB/rank = ~8 GB total
#   Scratch:    100 GB recommended

export LC_ALL=C
export LANG=C

SBATCH_LOG_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH}"

# ZFP HDF5 filter plugin path is set below, after SEAS_MFEM_ROOT resolves
# (it lives under ${SEAS_MFEM_ROOT}/extern/h5z-zfp/install/plugin).

# R-006 fix (R-002 round-3 extension): resolve seas-mfem root from the
# sbatch submit dir by walking UP until miniapps/seas/Makefile is found.
# Handles submission from:
#   - seas-mfem checkout root          (0 levels)
#   - miniapps/seas                    (2 levels)
#   - miniapps/seas/jobs/safs          (4 levels)
#   - anywhere else inside seas-mfem (up to 5 levels)
# Set SEAS_MFEM_ROOT explicitly to override.
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

# ZFP HDF5 filter plugin (libh5zzfp.so), required at runtime for
# --paraview-*-zfp-tol.  build_frontera.sh installs it under
# ${SEAS_MFEM_ROOT}/extern/h5z-zfp/install/plugin (same convention as the
# BP5 Phase-6 sbatches).  Fail fast with a build hint if it is missing.
H5Z_ZFP_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install"
if [[ -e "${H5Z_ZFP_BUILT_DIR}/plugin/libh5zzfp.so" ]]; then
    export HDF5_PLUGIN_PATH="${H5Z_ZFP_BUILT_DIR}/plugin:${HDF5_PLUGIN_PATH:-}"
else
    echo "ERROR: ZFP plugin not found at"
    echo "         ${H5Z_ZFP_BUILT_DIR}/plugin/libh5zzfp.so"
    echo "       Build it on Frontera with:  bash build_frontera.sh"
    echo "       (or omit --paraview-*-zfp-tol from the driver invocation)."
    exit 1
fi

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'   config/config.mk || { echo "ERROR: no HDF5"; exit 1; }
grep -q '^MFEM_USE_H5Z_ZFP *= YES' config/config.mk || { echo "ERROR: no H5Z_ZFP"; exit 1; }

cd miniapps/seas
make seas_spatial_dyn_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

SAFS_ROOT=safs/project_7.0_alternative
CONFIG_TOML="${SAFS_ROOT}/config/spatial_friction_slip_weakening_safs_projected_stress.toml"
MESH="${SAFS_ROOT}/meshing/results/msh/safs_fault_box_nwcut_1000m_lcfar3000.msh"

# Fail fast before launching the production allocation: the .msh is
# gitignored and does NOT arrive with `git pull`.  It must be scp'd to
# the cluster (or regenerated via
# ${SAFS_ROOT}/meshing/code/run_nwcut_meshing.py).
if [[ ! -f "${MESH}" ]]; then
    echo "ERROR: SAFS mesh not found:"
    echo "         $(pwd)/${MESH}"
    echo "       The .msh is gitignored and is NOT pulled with the repo."
    echo "       scp it to the path above before submitting this job."
    exit 1
fi

# Determinstic output dir (allows mid-run restart with the same path).
# Override JOBNAME via --export=ALL,JOBNAME=foo to override.
JOBNAME="${JOBNAME:-spatial_dyn_prod}"
OUT="${SCRATCH}/${JOBNAME}"
mkdir -p "${OUT}"

# Automatic restart detection: if any cp_*.h5 from a prior submission
# exists in OUT, resume from the most recent one.  Mirrors the BP5
# restart pattern at jobs/bp5/bp5_v2_restart_*.sbatch.
RESTART_FLAG=""
LATEST_CP=$(ls -1 "${OUT}/cp_t"*.h5 2>/dev/null | sort -V | tail -1)
if [[ -n "${LATEST_CP}" ]]; then
    RESTART_FLAG="--restart ${OUT}/cp"
    echo "Restart detected: resuming from ${LATEST_CP}"
    echo "NOTE: restart MUST use the same --ntasks and --nodes as the "
    echo "      original submission (V1 checkpoint format is per-rank)."
fi

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${OUT}"
echo "Mesh:       ${MESH}"
echo "Config:     ${CONFIG_TOML}"
echo "Restart:    ${RESTART_FLAG:-<none>}"

ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --no-sidecar-material \
    --tfinal 12s \
    --paraview \
    --paraview-fault-hdf5 \
    --paraview-fault-zfp-tol 1.0e-12 \
    --paraview-max-snapshots 5000 \
    --paraview-coseismic-dt 0.005 \
    --paraview-nucleation-dt 0.001 \
    --print-derived \
    --checkpoint-every 10000 \
    --output-dir "${OUT}" \
    ${RESTART_FLAG} \
    > "${OUT}/spatial_dyn_prod_${SLURM_JOB_ID}.log" 2>&1
RC=$?

echo "driver_rc=${RC}"
exit ${RC}
