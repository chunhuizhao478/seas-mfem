#!/bin/bash
#SBATCH -J bp5_v48e_blr14        # Job name
#SBATCH -o bp5_v48e_blr14_%j.out
#SBATCH -e bp5_v48e_blr14_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node
#SBATCH -n 8                         # 8 ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v48e: IP p=2, 2500m mesh, 8 ranks, BLR tol=1e-14
#
# If v48d (BLR 1e-12) doesn't work, does 1e-14?
# 1e-14 is near machine epsilon — essentially exact for double precision.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

make clean-build
make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }

mkdir -p bp5/results_v48e_2500m_blr14_8rank

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-14 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v48e_2500m_blr14_8rank \
      --output-prefix bp5_v48e_blr14 \
      --write-every-step
