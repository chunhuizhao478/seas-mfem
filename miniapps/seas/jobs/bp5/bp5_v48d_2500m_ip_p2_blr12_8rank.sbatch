#!/bin/bash
#SBATCH -J bp5_v48d_blr12        # Job name
#SBATCH -o bp5_v48d_blr12_%j.out
#SBATCH -e bp5_v48d_blr12_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node
#SBATCH -n 8                         # 8 ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v48d: IP p=2, 2500m mesh, 8 ranks, BLR tol=1e-12
#
# v48c (exact MUMPS, 8 ranks): V GROWS 0.010→0.019 — correct physics
# v48  (BLR 1e-10, 400 ranks): V DECAYS 0.010→0.004 — wrong physics
#
# Test: is BLR tol=1e-12 tight enough for correct nucleation?
# If V grows → 1e-12 is sufficient; use as default for ×3 penalty
# If V decays → need even tighter (1e-14) or exact MUMPS

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

mkdir -p bp5/results_v48d_2500m_blr12_8rank

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v48d_2500m_blr12_8rank \
      --output-prefix bp5_v48d_blr12 \
      --write-every-step
