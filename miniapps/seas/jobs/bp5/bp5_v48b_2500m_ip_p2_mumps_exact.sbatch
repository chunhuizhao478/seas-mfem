#!/bin/bash
#SBATCH -J bp5_v48b_exact        # Job name
#SBATCH -o bp5_v48b_exact_%j.out
#SBATCH -e bp5_v48b_exact_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v48b: IP p=2, 2500m mesh, 400 ranks, EXACT MUMPS (no BLR)
#
# Test whether MUMPS-BLR distributed factorization is the root cause of
# the parallel vs serial physics difference.
#
# v47l (serial, 2500m, MUMPS-BLR):  V_max GROWS 0.010 → 0.012 (nucleation)
# v48  (parallel, 2500m, MUMPS-BLR): V_max DECAYS 0.010 → 0.004 (no nucleation)
#
# If V_max GROWS with exact MUMPS in parallel → BLR is the problem
# If V_max still DECAYS → MUMPS accuracy is not the issue
#
# 2500m mesh at p=2: ~405K vector DOFs. Distributed across 400 ranks,
# exact MUMPS should fit in memory (~1K DOFs/rank).

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

mkdir -p bp5/results_v48b_2500m_mumps_exact

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v48b_2500m_mumps_exact \
      --output-prefix bp5_v48b_exact \
      --write-every-step
