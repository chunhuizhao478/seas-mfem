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
# write pipeline on the 1000 m zgraded SAFS mesh.  Reads
# `safs_smoke_8N_400r_dev.sbatch` for the queue / allocation header.
#
# Resource budget (SAFS 1000m zgraded, constant material, tfinal=0.5s):
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

# ZFP plugin path — Frontera-specific; adjust if site changes its layout.
# R-007 fix: validate ${WORK} is set and the plugin exists before using
# the default path, so the user gets an actionable error instead of a
# runtime VTKHDF abort at first write.
if [[ -z "${HDF5_PLUGIN_PATH:-}" ]]; then
    if [[ -z "${WORK:-}" ]]; then
        echo "WARNING: \${WORK} is unset; cannot derive HDF5_PLUGIN_PATH."
        echo "  Set HDF5_PLUGIN_PATH explicitly to the directory holding"
        echo "  libh5zzfp.{so,dylib}, or omit --paraview-*-zfp-tol from"
        echo "  the driver invocation."
    else
        export HDF5_PLUGIN_PATH="${WORK}/local/hdf5_plugins:${HDF5_PLUGIN_PATH:-}"
        if [[ ! -e "${WORK}/local/hdf5_plugins/libh5zzfp.so" ]]; then
            echo "WARNING: HDF5_PLUGIN_PATH=${WORK}/local/hdf5_plugins"
            echo "  but libh5zzfp.so not found there.  ZFP-compressed"
            echo "  ParaView output will abort at first write."
        fi
    fi
fi

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

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'   config/config.mk || { echo "ERROR: no HDF5"; exit 1; }
grep -q '^MFEM_USE_H5Z_ZFP *= YES' config/config.mk || { echo "ERROR: no H5Z_ZFP"; exit 1; }

cd miniapps/seas
make seas_spatial_dyn_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

SAFS_ROOT=safs/project_7.0_alternative
CONFIG_TOML="${SAFS_ROOT}/friction/EXAMPLE_spatial_friction_slip_weakening_safs.toml"
MESH="${SAFS_ROOT}/meshing/results/msh/safs_fault_box_nwcut_1000m_zgraded.msh"
OUT="${SCRATCH}/safs_dyn_smoke_${SLURM_JOB_ID}"
mkdir -p "${OUT}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${OUT}"
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
    > "${OUT}/spatial_dyn_smoke.log" 2>&1
RC=$?

if [[ ${RC} -ne 0 ]]; then
    echo "driver_rc=${RC} — see ${OUT}/spatial_dyn_smoke.log (probably "
    echo "  --print-derived ABORTed; check the [derived] lines)."
    exit ${RC}
fi

# Post-run verifier.  Requires conda env `pythonenv` (h5py).
# R-001 fix: use absolute paths anchored at ${SEAS_MFEM_ROOT} so the
# verifier resolves regardless of which directory we cd to.  Both the
# script and the TOML config live under ${SEAS_MFEM_ROOT}/miniapps/seas/.
module load python3
VERIFIER="${SEAS_MFEM_ROOT}/miniapps/seas/safs/project_7.0_alternative/spatial/code/scripts/verify_spatial_dyn_smoke_safs.py"
CONFIG_TOML_ABS="${SEAS_MFEM_ROOT}/miniapps/seas/${CONFIG_TOML}"
python3 "${VERIFIER}" \
    --fault-vtkhdf "${OUT}/fault.vtkhdf" \
    --toml-config "${CONFIG_TOML_ABS}" \
    --log "${OUT}/spatial_dyn_smoke.log" \
    > "${OUT}/verify_summary.txt"
CHECK_RC=$?

echo "driver_rc=${RC}, verify_rc=${CHECK_RC}"
exit $(( RC | CHECK_RC ))
