#!/bin/bash
#SBATCH -J bp5_mstats
#SBATCH -o bp5_mstats_%j.out
#SBATCH -e bp5_mstats_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# BP5 — Stiffness matrix analysis (dev queue, 2h).
# =======================================================================
# One-shot analysis of the assembled DG stiffness matrix used by the BP5
# production sbatch (bp5_phase6_paraview_zfp_normal_48hr.sbatch).
#
# What this job does:
#   1. Build seas_bp5_full identically to the production run.
#   2. Load the production mesh (bp5_tandem_exact.msh, 1000 m, P=1 IP).
#   3. Trigger one stiffness-matrix assembly (--matrix-stats).
#      - Writes <RESULT_DIR>/<OUTPUT_PREFIX>_matrix_stats.json with:
#        global_n, global_nnz (diag + offd), per-row nnz distribution,
#        DG block size, per-rank local "bandwidth", ||A||_F,
#        ||A - A^T||_F (symmetry residual), diag sign/range, K * 1.
#      - Bumps MUMPS print level to 2 — analysis-phase + factorize-phase
#        INFOG / RINFOG land in the .out file.  Look for:
#          " ANALYSIS STEP "
#          " ESTIMATED size in MBYTES for OOC "
#          " REAL space for factors           = "
#          " INTEGER space for factors        = "
#          " Maximum estimated front size     = "
#          " *** Time for factorization "
#   4. Performs one solve (the initial-condition solve), so MUMPS fully
#      factorizes the matrix.
#   5. --exit-after-factorize: bail before the time loop.  No ParaView,
#      no probes.  Total wall time ~ assembly + factorize.
#
# This is the cheapest possible run that exercises the same code path
# the production sbatch uses, so the matrix stats are guaranteed to
# describe exactly the matrix MUMPS factorizes on the long run.

export LC_ALL=C
export LANG=C
set -euo pipefail

SBATCH_LOG_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

if [ -z "${SEAS_MFEM_ROOT:-}" ]; then
    _probe="${SLURM_SUBMIT_DIR:-$PWD}"
    while [ "${_probe}" != "/" ] && [ -n "${_probe}" ]; do
        if [ -f "${_probe}/build_frontera.sh" ]; then
            SEAS_MFEM_ROOT="${_probe}"
            break
        fi
        _probe="$(dirname "${_probe}")"
    done
    unset _probe
fi
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-/scratch2/10024/zhaochun/seas-project/seas-mfem}"
echo "  SEAS_MFEM_ROOT = ${SEAS_MFEM_ROOT}"

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

HDF5_BUILT_DIR="${SEAS_MFEM_ROOT}/extern/hdf5/install"

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${HDF5_BUILT_DIR}/lib:${LD_LIBRARY_PATH:-}"

cd "${SEAS_MFEM_ROOT}"

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES'    config/config.mk || { echo "ERROR: no PETSc";   exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'    config/config.mk || { echo "ERROR: no MUMPS";   exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'     config/config.mk || { echo "ERROR: no HDF5";    exit 1; }

cd miniapps/seas
make clean
make seas_bp5_full -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="bp5/results_matrix_stats_dev_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="bp5_mstats_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR}"

echo "========================================"
echo "Git commit:    $(git rev-parse HEAD)"
echo "Git branch:    $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:    ${RESULT_DIR}"
echo "Mesh:          bp5_tandem_exact.msh (1000 m, P=1 IP) [same as prod sbatch]"
echo "Solver:        MUMPS direct (print-level 2 — fill-in stats logged)"
echo "Mode:          one assembly + one factorize, then exit"
echo "Queue:         development (8N x 400r, 2h)"
echo "========================================"

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver mumps \
      --ref-dir bp5/benchmark_data \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --matrix-stats \
      --exit-after-factorize \
      --tfinal 1

JSON="${RESULT_DIR}/${OUTPUT_PREFIX}_matrix_stats.json"
echo
echo "========================================"
echo "Matrix-stats JSON:"
echo "========================================"
if [ -f "${JSON}" ]; then
   cat "${JSON}"
else
   echo "WARNING: ${JSON} not produced — check rank-0 output above."
fi
echo
echo "========================================"
echo "Greppable MUMPS factor stats (from .out — useful summary):"
echo "========================================"
{
   echo "MUMPS analysis / factor lines:"
   grep -E " (ANALYSIS|FACTORIZATION) STEP " "${SBATCH_LOG_DIR}/bp5_mstats_${SLURM_JOB_ID}.out" || true
   grep -E "(REAL|INTEGER) space for factors" "${SBATCH_LOG_DIR}/bp5_mstats_${SLURM_JOB_ID}.out" || true
   grep -E "Maximum estimated front size"     "${SBATCH_LOG_DIR}/bp5_mstats_${SLURM_JOB_ID}.out" || true
   grep -E "ESTIMATED size in MBYTES"         "${SBATCH_LOG_DIR}/bp5_mstats_${SLURM_JOB_ID}.out" || true
   grep -E "Time for factorization"           "${SBATCH_LOG_DIR}/bp5_mstats_${SLURM_JOB_ID}.out" || true
} | tee "${RESULT_DIR}/mumps_summary.txt"
