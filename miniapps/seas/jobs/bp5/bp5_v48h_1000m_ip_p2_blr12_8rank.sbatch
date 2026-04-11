#!/bin/bash
#SBATCH -J bp5_v48h_1k_8r       # Job name
#SBATCH -o bp5_v48h_1k_8r_%j.out
#SBATCH -e bp5_v48h_1k_8r_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 4                         # 4 nodes (memory headroom for BLR on large mesh)
#SBATCH -n 8                         # 8 ranks
#SBATCH -t 04:00:00                  # 4 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v48h: IP p=2, 1000m mesh, 8 ranks, BLR tol=1e-12
#
# Isolate rank count vs mesh resolution:
#   v48d (2500m, 8 ranks, BLR 1e-12): WORKS — V grows
#   v48g (1000m, 400 ranks, BLR 1e-12): CRASH — same blowup as before
#
# If this WORKS → the issue is rank count (shared face density), not mesh
# If this CRASHES → the issue is mesh resolution (1000m has different properties)
#
# 1000m mesh: 63451 elements / 8 ranks = ~8000 elements/rank
# At p=2: ~240K DOFs/rank. BLR should fit in 4 nodes (more memory per rank).

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

mkdir -p bp5/results_v48h_1000m_blr12_8rank

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v48h_1000m_blr12_8rank \
      --output-prefix bp5_v48h_1k_8r \
      --write-every-step
