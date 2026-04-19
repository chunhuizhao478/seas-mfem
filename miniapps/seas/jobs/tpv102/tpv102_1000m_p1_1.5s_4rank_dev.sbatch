#!/bin/bash
#SBATCH -J tpv102_1000m_4r
#SBATCH -o tpv102_1000m_4r_%j.out
#SBATCH -e tpv102_1000m_4r_%j.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 00:30:00
#SBATCH -A EAR20006

# TPV102 dev-queue shared-fault sanity run (v2 recommendation #1).
#
# Purpose: after the v2 fix round (R-001 swap + R-101 runtime
# consistency diagnostic + R-102 hard-abort peer-rank guard), run a
# cheap 4-rank job on the 1000 m coarse mesh to verify:
#   - `--debug-qnorm` shows non-hypocenter ranks' ||Q||_inf > 0
#     within ~0.2 s (wave transits 1 km at c_s = 3464 m/s ≈ 0.29 ms)
#   - the R-101 runtime diagnostic either passes or SKIPS (if METIS
#     places all fault faces interior to one rank; see test_r101
#     for discussion)
#   - no MFEM_VERIFY abort from R-102 on any shared face
#
# NOTE on rupture physics: at h = 1000 m, the TPV102 cohesive process
# zone Λ_dyn ≈ 160 m gives Λ/h ≪ 1, so this run CANNOT be used to
# verify rupture propagation — that's a 200 m-or-finer question
# (see tpv102_debug_v1.md § "Next-step plan" for discussion).
# This run only verifies the cross-rank-wave-and-consistency path
# introduced by the v1 + v2 fixes.
#
# Expected wall time: ~15 min at 1000 m, 4 ranks, tfinal=1.5 s.
#
# Usage:
#   sbatch jobs/tpv102/tpv102_1000m_p1_1.5s_4rank_dev.sbatch

export LC_ALL=C
export LANG=C

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
make seas_tpv102_driver -j"${JOBS}" || { echo "BUILD FAILED"; exit 1; }

RESULT_DIR="tpv102/results_1000m_p1_1.5s_4r_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="tpv102"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "tfinal=1.5 s, dt=auto (CFL), order=1, debug-qnorm ON"

ibrun ./seas_tpv102_driver \
      --mesh tpv102/mesh/tpv102_1000m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 1.5 \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --output-dt 0.05 \
      --debug-qnorm
