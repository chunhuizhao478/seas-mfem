#!/bin/bash
#SBATCH -J seas_qd_bp5_50step
#SBATCH -o seas_qd_bp5_50step_%j.out
#SBATCH -e seas_qd_bp5_50step_%j.err
#SBATCH -p development
#SBATCH -N 4
#SBATCH -n 200
#SBATCH -t 01:00:00
#SBATCH -A EAR20006

# =========================================================================
# Phase 5 acceptance — quasi-dynamic (QD) SEAS driver, 50-step BP5 run.
# =========================================================================
# Runs the FULL QD coupling on the BP5 mesh (the equilibrium<1e-6 4-phase init +
# the 50-step RK45 advance — the Phase-5 acceptance criteria, which need a
# SAFS-CLEAN fault basis and so cannot run locally on the tiny test fixtures).
#
# PREREQUISITE: build the binary ONCE first (this job does NOT build in-job):
#   cd miniapps/seas && make seas_spatial_seas_driver
# (the QD driver is a DIFFERENT binary from seas_spatial_dyn_driver; the
#  spatial_dyn_build.sbatch jobs do not build it.)
#
# Submit (from miniapps/seas):
#   sbatch jobs/safs_qd/spatial_seas_bp5_50step.sbatch
#
# Expected (PASS): rank-0 stdout shows
#   [spatial_seas] SetInitialCondition done: V_max = <finite> m/s
#   ... and the 4-phase init did NOT abort (equilibrium residual < 1e-6), then
#   50 "[spatial_seas] step N  t=..  dt=..  V_max=.." lines with V_max finite
#   (CLAUDE.md regression list: V_max finite, dt growing toward dt_max).

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
export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH:-}"

# --- Run from miniapps/seas (where the binary + config live) ---
cd "${SLURM_SUBMIT_DIR}"

BINARY=./seas_spatial_seas_driver
if [[ ! -x "${BINARY}" ]]; then
    echo "ERROR: ${BINARY} not found/executable in $(pwd)."
    echo "       Build it first:  cd miniapps/seas && make seas_spatial_seas_driver"
    exit 1
fi

CONFIG_TOML="${SEAS_QD_CONFIG:-config/safs_qd/bp5_phase5_50step.toml}"
if [[ ! -f "${CONFIG_TOML}" ]]; then
    echo "ERROR: config ${CONFIG_TOML} not found."
    exit 1
fi

MAX_STEPS="${SEAS_QD_MAX_STEPS:-50}"

echo "[job] $(date)  ${BINARY}  config=${CONFIG_TOML}  max-steps=${MAX_STEPS}"
ibrun "${BINARY}" \
    --config "${CONFIG_TOML}" \
    --max-steps "${MAX_STEPS}"
echo "[job] done $(date)"
