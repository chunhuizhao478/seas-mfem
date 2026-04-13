#!/bin/bash
#SBATCH -J bp5_p250
#SBATCH -o bp5_parametric_250m_%j.out
#SBATCH -e bp5_parametric_250m_%j.err
#SBATCH -p normal
#SBATCH -N 16
#SBATCH -n 896
#SBATCH -t 48:00:00
#SBATCH -A EAR20006

# Parametric study: 250m near-fault resolution (baseline is 1000m).
# Same domain size (400 x 200 x 100 km), 4x finer mesh near fault.
# Expected ~1.1M elements (~15x baseline). Uses 16 nodes (896 ranks)
# to handle the larger problem size.
#
# Compares against 1000m golden (golden_parallel_400r_p1_1000m).
# Expect small differences due to mesh convergence — tolerance 5%.
#
# Before submitting:
#   1. Generate mesh (if not already done):
#      conda activate pythonenv
#      cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas
#      gmsh -3 bp5/parametric_study/case_250m/bp5_250m.geo \
#           -o bp5/parametric_study/case_250m/bp5_250m.msh
#   2. Build driver (once, before submitting any parametric jobs):
#      make clean && make seas_driver -j8
#   3. Submit:
#      cd miniapps/seas && sbatch jobs/bp5/bp5_parametric_250m.sbatch
#
# Post-run comparison:
#   python3 scripts/regression_check.py \
#       --sim-dir $RESULT_DIR --sim-prefix $OUTPUT_PREFIX \
#       --ref-dir bp5/benchmark_data/golden_parallel_400r_p1_1000m \
#       --tolerance 0.05

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

grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas

test -f extern/toml11/include/toml.hpp || test -f extern/toml11/toml.hpp || \
   { echo "ERROR: toml11 not initialized. Run: git submodule update --init miniapps/seas/extern/toml11"; exit 1; }

test -f bp5/parametric_study/case_250m/bp5_250m.msh || \
   { echo "ERROR: Mesh not found. Generate with: gmsh -3 bp5/parametric_study/case_250m/bp5_250m.geo -o bp5/parametric_study/case_250m/bp5_250m.msh"; exit 1; }

test -x ./seas_driver || { echo "ERROR: seas_driver not built. Build first with: make clean && make seas_driver -j8"; exit 1; }

RESULT_DIR="bp5/parametric_study/case_250m/results_job${SLURM_JOB_ID}"
OUTPUT_PREFIX="bp5_parametric_250m_job${SLURM_JOB_ID}"

mkdir -p "${RESULT_DIR}"

echo "Git commit: $(git rev-parse HEAD)"
echo "Git branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir: ${RESULT_DIR}"
echo "Parametric study: 250m near-fault resolution, 896 ranks, 48hr"

ibrun ./seas_driver bp5/parametric_study/case_250m/bp5_250m.toml \
      --override output.output_dir="${RESULT_DIR}" \
      --override output.output_prefix="${OUTPUT_PREFIX}"
