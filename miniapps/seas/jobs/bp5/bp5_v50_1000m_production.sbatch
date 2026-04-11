#!/bin/bash
#SBATCH -J bp5_v50_prod           # Job name
#SBATCH -o bp5_v50_prod_%j.out
#SBATCH -e bp5_v50_prod_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks — production scale
#SBATCH -t 04:00:00                  # 4 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v50: PRODUCTION RUN — automatic CFL dt + V-guard
#
# First run with NO manual --dt-init or --v-guard flags.
# The code automatically:
#   1. Computes h_min from mesh → dt_CFL = 2*eta*h_min/(4*mu)
#   2. Sets dt_init = min(dt_V, dt_CFL)
#   3. Enables V-guard with factor 100 (default)
#
# Expected: dt_init ≈ 0.072s (CFL-limited for 1000m mesh).
# Should match v49m behavior (stable, healthy nucleation).

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# Build with v50 code changes
make clean-build
make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }

mkdir -p bp5/results_v50_production

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v50_production \
      --output-prefix bp5_v50_prod \
      --write-every-step
