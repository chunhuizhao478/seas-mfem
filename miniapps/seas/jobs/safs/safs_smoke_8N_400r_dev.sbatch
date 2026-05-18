#!/bin/bash
#SBATCH -J safs_smoke_8N_400r
#SBATCH -o safs_smoke_8N_400r_%j.out
#SBATCH -e safs_smoke_8N_400r_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# SAFS multi-fault SEAS smoke driver (mesh shakedown).
#
# Purpose: run the SAFS production cavity .msh end-to-end through DG
# elasticity + uniform-VW rate-and-state friction for a small number of
# RK45 steps, and check (via check_safs_smoke.py) that the per-step
# diagnostics stay physical. This is NOT a benchmark validation — there
# is no SCEC reference. It is the canary for "does the 6-fault embedded
# triangulation blow up the solver?".
#
# Resource budget: 8 nodes × 50 cores/node = 400 ranks, 2 hr dev queue.

export LC_ALL=C
export LANG=C

# Pin the submit dir before any `cd` so the post-run check can find the
# output CSV from a stable absolute path.
SBATCH_LOG_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH}"

cd /scratch2/10024/zhaochun/seas-project/seas-mfem

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas
make seas_safs_test_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="safs/code/results_8N_400r_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

MESH=safs/mesh/output/newset_6_all6_garnetfirst/output/safs_newset_6_cavity.msh

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Mesh: ${MESH}"

ibrun ./seas_safs_test_driver \
    --mesh "${MESH}" \
    --output-csv "${RESULT_DIR}/safs_smoke.csv" \
    --n-steps 100 \
    --order 1 \
    --solver mumps \
    --print-every 10 \
    > "${RESULT_DIR}/run.log" 2>&1
RC=$?

# Post-run acceptance check
cd "${SBATCH_LOG_DIR}"
python miniapps/seas/safs/code/check_safs_smoke.py \
    "${RESULT_DIR}/safs_smoke.csv" \
    > "${RESULT_DIR}/pass_check.txt"
CHECK_RC=$?

echo "driver_rc=${RC}, check_rc=${CHECK_RC}"
exit $(( RC | CHECK_RC ))
