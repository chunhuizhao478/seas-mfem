#!/bin/bash
#SBATCH -J safs_qd_v1_dev
#SBATCH -o safs_qd_v1_dev_%j.out
#SBATCH -e safs_qd_v1_dev_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# =========================================================================
# SAF QD 500 m DEV SMOKE — STRONG rate weakening (FL=103) from v1.0.0 RSSRW.
#   driver : seas_spatial_seas_driver  (quasi-dynamic / cycle model)
#   config : config/safs_qd/safs_qd_500m_rssrw_v1.toml
# =========================================================================
# !!! REQUIRES THE SRW-IN-QD WIRING (NOT YET IMPLEMENTED) !!!
# The config selects state_evolution="slip_law_srw".  Until the QD fault
# operator is wired to call SlipLawSRWPsi's per-DOF _SRW overloads
# (document/spatial_seas_dev/PLAN_spatial_seas_srw_qd_2026-06-20.md), the driver
# ABORTS on this config (by design — no silent aging-law fall-back).  This job
# is therefore staged-and-ready: it will run the SRW cycle model once the wiring
# lands.  For a standard aging-law cycle-model smoke today, point CONFIG_TOML at
# an aging-law config (state_evolution="aging_law").
#
# Acceptance (mechanical, once SRW is wired): builds, SetInitialCondition
# equilibrium < 1e-6, advances >= 100 steps with finite V_max, AMG converging
# each step.  The loading is the BP5 planar sgn(Y) approximation (Vp/2=+/-1e-9
# m/s on the side walls) the user accepted; physics validation is the BP5 gate
# (Phase 7) + the eventual SAF cycle output, not this smoke.
#
# PREREQUISITE: build the binary ONCE first (this job does NOT build in-job):
#   cd miniapps/seas && make seas_spatial_seas_driver
#
# Submit (from miniapps/seas):
#   sbatch jobs/safs_qd/safs_qd_500m_rssrw_v1_dev_2hr.sbatch

set -o pipefail   # NOT -u: conda/module hooks reference unbound vars
export LC_ALL=C
export LANG=C

# --- Frontera module stack (match a working spatial sbatch) ---
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

CONFIG_TOML="${SAFS_QD_CONFIG:-config/safs_qd/safs_qd_500m_rssrw_v1.toml}"
MAX_STEPS="${SAFS_QD_MAX_STEPS:-200}"   # >= 100 for the dev smoke

# Optional solver-scaling capture (no-op unless built with MFEM_USE_CALIPER=YES).
export CALI_CONFIG="${CALI_CONFIG:-runtime-report(output=safs_qd_v1_dev_${SLURM_JOB_ID}_cali.txt)}"

echo "[job] $(date)  ${BINARY}  config=${CONFIG_TOML}  max-steps=${MAX_STEPS}  ranks=${SLURM_NTASKS}"
ibrun "${BINARY}" \
    --config "${CONFIG_TOML}" \
    --max-steps "${MAX_STEPS}"
echo "[job] done $(date)"
