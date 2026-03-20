#!/bin/bash
#SBATCH -J bp5_v47m_par          # Job name
#SBATCH -o bp5_v47m_par_%j.out
#SBATCH -e bp5_v47m_par_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node
#SBATCH -n 48                        # 48 MPI ranks (shared faces present)
#SBATCH -t 00:30:00                  # 30 min (will crash quickly if bug present)
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v47m: IP p=2, 2500m mesh, 48 RANKS — CONFIRM BLOWUP IN PARALLEL
#
# Counterpart to v47l (serial). Same mesh, same parameters, just multi-rank.
# If this blows up and v47l doesn't → shared face bug confirmed.
# 48 ranks on 13500 tets = ~280 tets/rank — reasonable partitioning.

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

mkdir -p bp5/results_v47m_2500m_parallel

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v47m_2500m_parallel \
      --output-prefix bp5_v47m_par \
      --write-every-step
