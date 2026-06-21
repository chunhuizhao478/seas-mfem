#!/bin/bash
#SBATCH -J safs_qd_v1_prod
#SBATCH -o safs_qd_v1_prod_%j.out
#SBATCH -e safs_qd_v1_prod_%j.err
#SBATCH -p normal
#SBATCH -N 16
#SBATCH -n 800
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# =========================================================================
# SAF QD 500 m PRODUCTION — STRONG rate weakening (FL=103) cycle model.
#   driver : seas_spatial_seas_driver  (quasi-dynamic / sequence-of-events)
#   config : config/safs_qd/safs_qd_500m_rssrw_v1.toml
# =========================================================================
# !!! REQUIRES THE SRW-IN-QD WIRING (see the dev job + the config header) AND
#     a passing BP5 parity gate (Phase 7) before this is run for science. !!!
# Runs many seismic cycles (interseismic loading -> spontaneous nucleation ->
# event -> postseismic), driven by plate loading Vp/2 = +/-1e-9 m/s on the
# front/back walls.  Checkpoint/restart is enabled so a 48 h block resumes.
#
# Restart (project memory: restart out-dir MUST differ from the checkpoint
# parent): pass  --restart <prev_out>/cp  --output-dir <new_out>.
#
# PREREQUISITE: build the binary first (cd miniapps/seas && make
# seas_spatial_seas_driver); pre-check output size with
# scripts/estimate_output_size.py before submitting.
#
# Submit (from miniapps/seas):
#   sbatch jobs/safs_qd/safs_qd_500m_rssrw_v1_production.sbatch

set -o pipefail
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
    echo "       Build it first:  cd miniapps/seas && make seas_spatial_seas_driver"
    exit 1
fi

CONFIG_TOML="${SAFS_QD_CONFIG:-config/safs_qd/safs_qd_500m_rssrw_v1.toml}"
OUT_DIR="${SAFS_QD_OUT:-output_safs_qd_v1_0_0_RSSRW_${SLURM_JOB_ID}}"
CKPT_EVERY="${SAFS_QD_CKPT_EVERY:-2000}"

# Caliper solver-scaling report (no-op unless built with MFEM_USE_CALIPER=YES).
export CALI_CONFIG="${CALI_CONFIG:-runtime-report(output=safs_qd_v1_prod_${SLURM_JOB_ID}_cali.txt)}"

echo "[job] $(date)  ${BINARY}  config=${CONFIG_TOML}  out=${OUT_DIR}  ranks=${SLURM_NTASKS}"
ibrun "${BINARY}" \
    --config "${CONFIG_TOML}" \
    --output-dir "${OUT_DIR}" \
    --checkpoint-every "${CKPT_EVERY}" \
    "${SAFS_QD_RESTART_ARGS:-}"
echo "[job] done $(date)"
