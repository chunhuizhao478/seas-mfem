#!/bin/bash
#SBATCH -J spatial_dyn_smoke
#SBATCH -o spatial_dyn_smoke_%j.out
#SBATCH -e spatial_dyn_smoke_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# Phase Z smoke for seas_spatial_dyn_driver — first SAFS dynamic-rupture
# production run, dev queue, 2 hr wall, 0.5 s simulated time.
#
# Purpose: end-to-end shakedown of the gradual_overstress nucleation +
# --print-derived pre-flight gate + ParaView fault.vtkhdf + checkpoint
# write pipeline on the 1000 m lcfar3000 SAFS mesh with regional H&Z
# projected stress, constant material, mixed-flux "adjacent", and O2
# ADER.  Reads `safs_smoke_8N_400r_dev.sbatch` for the queue header.
#
# Resource budget (SAFS 1000m lcfar3000, constant material, tfinal=0.5s):
#   Mesh:       ~1.2M tets, ~36M DOFs (3D velocity-stress, order 1)
#   Per-rank:   ~0.6 GB peak (wave Q + DOFData + GodunovFlux)
#   Wall:       ~1.5-2 hr at cfl=0.5, dt~6e-5 s, ~8000 steps
#   Output:     fault.vtkhdf ~5 GB with ZFP tol 1e-12
#               cp files every 2000 steps × ~3 MB/rank
#   Scratch:    50 GB recommended

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

# Fail fast before launching 400 ranks: the .msh is gitignored and does
# NOT arrive with `git pull`.  It must be scp'd to the cluster (or
# regenerated via ${SAFS_ROOT}/meshing/code/run_nwcut_meshing.py).
if [[ ! -f "${MESH}" ]]; then
    echo "ERROR: SAFS mesh not found:"
    echo "         $(pwd)/${MESH}"
    echo "       The .msh is gitignored and is NOT pulled with the repo."
    echo "       scp it to the path above before submitting this job."
    exit 1
fi

# ParaView/checkpoint output → repo safs/ tree (per request); this is a
# shakedown run, so output lives next to the SAFS project for easy
# inspection.  The run dir is gitignored (safs/safs_dyn_smoke_*/).
OUT="${SEAS_MFEM_ROOT}/miniapps/seas/safs/safs_dyn_smoke_${SLURM_JOB_ID}"
# Driver run log + post-run verify summary → jobs/safs (next to this
# sbatch).  The SLURM %j.out/%j.err land in the submit dir, so submit
# from miniapps/seas/jobs/safs to co-locate them here.
LOG_DIR="${SEAS_MFEM_ROOT}/miniapps/seas/jobs/safs"
LOG="${LOG_DIR}/spatial_dyn_smoke_${SLURM_JOB_ID}.log"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${OUT}"
echo "Log file:   ${LOG}"
echo "Mesh:       ${MESH}"
echo "Config:     ${CONFIG_TOML}"

ibrun ./seas_spatial_dyn_driver \
    --config "${CONFIG_TOML}" \
    --mesh   "${MESH}" \
    --no-sidecar-material \
    --tfinal 0.5s \
    --paraview \
    --paraview-fault-hdf5 \
    --paraview-fault-zfp-tol 1.0e-12 \
    --paraview-max-snapshots 200 \
    --print-derived \
    --checkpoint-every 2000 \
    --output-dir "${OUT}" \
    > "${LOG}" 2>&1
RC=$?

if [[ ${RC} -ne 0 ]]; then
    echo "driver_rc=${RC} — see ${LOG} (probably "
    echo "  --print-derived ABORTed; check the [derived] lines)."
    exit ${RC}
fi

# Post-run verifier — OPTIONAL.  It needs h5py, which is NOT in Frontera's
# python3 module and has no conda env on the cluster.  Probe for h5py and
# run the verifier only if it imports; otherwise SKIP with a hint and do
# NOT fail the job — the simulation (driver_rc above) is what matters, and
# a missing post-processing dependency must not mark a good run as failed.
# To verify, run the printed command on a host that has h5py (e.g. scp
# ${OUT}/fault.vtkhdf to your local pythonenv).
# R-001 fix: absolute paths anchored at ${SEAS_MFEM_ROOT} so the verifier
# resolves regardless of which directory we cd to.
module load python3
VERIFIER="${SEAS_MFEM_ROOT}/miniapps/seas/safs/project_7.0_alternative/spatial/code/scripts/verify_spatial_dyn_smoke_safs.py"
CONFIG_TOML_ABS="${SEAS_MFEM_ROOT}/miniapps/seas/${CONFIG_TOML}"
if python3 -c "import h5py" 2>/dev/null; then
    python3 "${VERIFIER}" \
        --fault-vtkhdf "${OUT}/fault.vtkhdf" \
        --toml-config "${CONFIG_TOML_ABS}" \
        --log "${LOG}" \
        > "${LOG_DIR}/verify_summary_${SLURM_JOB_ID}.txt"
    CHECK_RC=$?
    echo "driver_rc=${RC}, verify_rc=${CHECK_RC}"
    exit $(( RC | CHECK_RC ))
else
    echo "driver_rc=${RC}, verify=SKIPPED (no h5py on this host)."
    echo "  Run the verifier on a host that has h5py:"
    echo "    python3 ${VERIFIER} \\"
    echo "        --fault-vtkhdf ${OUT}/fault.vtkhdf \\"
    echo "        --toml-config ${CONFIG_TOML_ABS} \\"
    echo "        --log ${LOG}"
    exit ${RC}
fi
