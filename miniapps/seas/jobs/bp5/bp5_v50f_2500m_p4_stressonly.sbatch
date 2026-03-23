#!/bin/bash
#SBATCH -J bp5_v50f_p4_so         # Job name
#SBATCH -o bp5_v50f_p4_so_%j.out
#SBATCH -e bp5_v50f_p4_so_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 4                         # 4 nodes
#SBATCH -n 200                       # 200 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v50f: 2500m mesh, p=4, IP, STRESS-ONLY traction (Strategy 1)
#
# Tests whether dropping the penalty correction in ComputeTraction
# fixes the p=4 nucleation failure.
#
# The penalty still appears in K (stiffness matrix) and f (RHS).
# Only the traction RECOVERY changes: T = {sigma(u)·n} only.
#
# Expected: If cancellation hypothesis is correct, V should GROW.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# BUILD job
make clean-build
make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }

mkdir -p bp5/results_v50f_p4_stressonly

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 4 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v50f_p4_stressonly \
      --output-prefix bp5_v50f_p4_so \
      --traction-stress-only \
      --write-every-step
