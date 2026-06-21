#!/bin/bash
#SBATCH -J safs_qd_500m_dev
#SBATCH -o safs_qd_500m_dev_%j.out
#SBATCH -e safs_qd_500m_dev_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# =========================================================================
# Phase 8 — SAF QD 500 m DEV SMOKE (spatial_seas driver).
# =========================================================================
# !!! The config config/safs_qd/safs_qd_500m.toml is a SKELETON with PHYSICS-
# DESIGN PLACEHOLDERS (boundary loading, velocity model, nucleation, sigma_n).
# In particular the boundary LOADING FUNCTION is still BP5 sgn(Y) — WRONG for the
# curved SAF fault — so results are NOT physically meaningful until the SAF
# loading function is implemented + the placeholders are filled from the SCEC
# SAF setup.  This dev job exercises BUILD + INIT + a few steps + AMG convergence
# on the real mesh; it is NOT a physics validation.  See the config header.
#
# Acceptance (mechanical, not physics): builds, SetInitialCondition equilibrium
# < 1e-6, advances >= 100 steps with finite V_max, AMG converging each step.
#
# PREREQUISITE: build the binary ONCE first (this job does NOT build in-job):
#   cd miniapps/seas && make seas_spatial_seas_driver
# For solver scaling (acceptance #2) build with MFEM_USE_CALIPER=YES and set
# CALI_CONFIG=runtime-report (or a perfgraph profile).
#
# Submit (from miniapps/seas):
#   sbatch jobs/safs_qd/safs_qd_500m_dev_2hr.sbatch

set -euo pipefail
export LC_ALL=C
export LANG=C

# --- Frontera module stack (match a working spatial sbatch — [[feedback_sbatch_modules]]) ---
module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8
# HYPRE-shadow fix: TACC's HYPRE must precede any bundled copy on the path.
export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH:-}"

cd "${SLURM_SUBMIT_DIR}"

BINARY=./seas_spatial_seas_driver
if [[ ! -x "${BINARY}" ]]; then
    echo "ERROR: ${BINARY} not found/executable in $(pwd)."
    echo "       Build it first:  cd miniapps/seas && make seas_spatial_seas_driver"
    exit 1
fi

CONFIG_TOML="${SAFS_QD_CONFIG:-config/safs_qd/safs_qd_500m.toml}"
MAX_STEPS="${SAFS_QD_MAX_STEPS:-200}"   # >= 100 for the dev smoke

# Optional solver-scaling capture (no-op unless built with MFEM_USE_CALIPER=YES).
export CALI_CONFIG="${CALI_CONFIG:-runtime-report(output=safs_qd_500m_dev_${SLURM_JOB_ID}_cali.txt)}"

echo "[job] $(date)  ${BINARY}  config=${CONFIG_TOML}  max-steps=${MAX_STEPS}  ranks=${SLURM_NTASKS}"
ibrun "${BINARY}" \
    --config "${CONFIG_TOML}" \
    --max-steps "${MAX_STEPS}"
echo "[job] done $(date)"
