#!/bin/bash
#SBATCH -J bp5_v48c_exact8       # Job name
#SBATCH -o bp5_v48c_exact8_%j.out
#SBATCH -e bp5_v48c_exact8_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node
#SBATCH -n 8                         # 8 ranks (small enough for exact MUMPS)
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v48c: IP p=2, 2500m mesh, 8 ranks, EXACT MUMPS (no BLR)
#
# v48b failed: 400 ranks too many for exact MUMPS on this problem size.
# MUMPS INFO(1)=-1, INFO(2)=215 (parallel distribution error).
#
# 8 ranks: ~13500 elements, ~50K DOFs/rank. Exact MUMPS should handle this.
# Still has shared faces to test the parallel physics.
#
# KEY: compare V_max trend with:
#   v47l (serial, BLR):   V GROWS 0.010 → 0.012
#   v48  (400 rank, BLR): V DECAYS 0.010 → 0.004
#
# If V GROWS → BLR parallel accuracy is the issue
# If V DECAYS → something else

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

mkdir -p bp5/results_v48c_2500m_mumps_exact_8rank

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v48c_2500m_mumps_exact_8rank \
      --output-prefix bp5_v48c_exact8 \
      --write-every-step
