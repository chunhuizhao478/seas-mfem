#!/bin/bash
#SBATCH -J safs_qd_500m_prod
#SBATCH -o safs_qd_500m_prod_%j.out
#SBATCH -e safs_qd_500m_prod_%j.err
#SBATCH -p normal
#SBATCH -N 16
#SBATCH -n 800
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# =========================================================================
# Phase 8 — SAF QD 500 m PRODUCTION run (spatial_seas driver).
# =========================================================================
# !!! DO NOT SUBMIT until the config config/safs_qd/safs_qd_500m.toml physics-
# design PLACEHOLDERS are resolved from the SCEC SAF setup AND the SAF boundary
# LOADING FUNCTION is implemented (the QD driver currently uses BP5 sgn(Y),
# which is WRONG for the curved SAF fault).  See the config header.  Until then
# this run is physically meaningless.
#
# Acceptance (physics): an inter-seismic -> nucleation -> event sequence with
# plausible recurrence/slip and no CLAUDE.md regression trips.
#
# PREREQUISITE: build once first (cd miniapps/seas && make seas_spatial_seas_driver).
# Pre-submit output-size estimate (avoid scratch-quota overrun):
#   python3 scripts/estimate_output_size.py --driver bp5 ... (adapt for QD)
#
# Submit (from miniapps/seas):
#   sbatch jobs/safs_qd/safs_qd_500m_production.sbatch
# Chained restart (after a 48 h wall-clock limit), with a DIFFERENT --output-dir:
#   sbatch jobs/safs_qd/safs_qd_500m_production.sbatch   # auto-restart wiring is
#   the driver's --restart PREFIX (Phase 6) — set SAFS_QD_RESTART below.

set -euo pipefail
export LC_ALL=C
export LANG=C

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8
export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH:-}"

cd "${SLURM_SUBMIT_DIR}"

BINARY=./seas_spatial_seas_driver
if [[ ! -x "${BINARY}" ]]; then
    echo "ERROR: ${BINARY} not found/executable in $(pwd)."
    exit 1
fi

CONFIG_TOML="${SAFS_QD_CONFIG:-config/safs_qd/safs_qd_500m.toml}"
OUTPUT_DIR="${SAFS_QD_OUTPUT_DIR:-./safs_qd_500m_prod_${SLURM_JOB_ID}}"
RESTART="${SAFS_QD_RESTART:-}"   # set to a checkpoint prefix to resume

mkdir -p "${OUTPUT_DIR}"
ARGS=(--config "${CONFIG_TOML}" --output-dir "${OUTPUT_DIR}" --checkpoint-every 2000)
if [[ -n "${RESTART}" ]]; then
    # Restart safety (driver-enforced): --output-dir must differ from --restart's parent.
    ARGS+=(--restart "${RESTART}")
fi

echo "[job] $(date)  ${BINARY}  config=${CONFIG_TOML}  out=${OUTPUT_DIR}  ranks=${SLURM_NTASKS}"
ibrun "${BINARY}" "${ARGS[@]}"
echo "[job] done $(date)"
