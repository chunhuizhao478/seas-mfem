#!/bin/bash
#SBATCH -J g0_build
#SBATCH -o g0_build_%j.out
#SBATCH -e g0_build_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -c 72
#SBATCH -t 03:00:00
#SBATCH -A EAR20006

# GPU plan G0: CPU-only build on a Vista gg node (login nodes give one core).
# Builds PETSc 3.24.x with MUMPS, ScaLAPACK, METIS, ParMETIS, PT-Scotch and
# hypre from the tarballs that `STAGE=fetch bash build_horizon.sh` put in
# $WORK/seas-gpu-deps/downloads on a login node, then MFEM and seas_bp5_full.
#
# Submit from miniapps/seas:  sbatch jobs/vista/g0_build.sbatch
# Then run the unit tests G0's code touches; the job fails if any fails.

set -euo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/vista/g0_env.sh"

cd "${SEAS_MFEM_ROOT}"
echo "Checkout: ${SEAS_MFEM_ROOT} at $(git rev-parse HEAD) ($(git rev-parse --abbrev-ref HEAD))"

MODE=cpu STAGE="deps mfem seas" JOBS=64 MODULES=none bash build_horizon.sh

cd "${SEAS_DIR}"
echo "=== unit test: seas_test_elasticity_operator (1 rank) ==="
ibrun -n 1 ./seas_test_elasticity_operator > g0_build_unit_elasticity.log 2>&1
tail -n 4 g0_build_unit_elasticity.log
grep -q ', 0 failed out of' g0_build_unit_elasticity.log
for np in 1 2; do
   echo "=== unit test: seas_test_g0_matrix_only (${np} ranks) ==="
   ibrun -n "${np}" ./seas_test_g0_matrix_only > "g0_build_unit_matrix_only_np${np}.log" 2>&1
   tail -n 3 "g0_build_unit_matrix_only_np${np}.log"
done
echo "=== unit test: bp5_gate_check.py ==="
python3 tests/unit/test_bp5_gate_check.py
