#!/bin/bash
#SBATCH -J bp5_v34_p1_full     # Job name
#SBATCH -o bp5_v34_p1_full_%j.out
#SBATCH -e bp5_v34_p1_full_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # MPI tasks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v34: Full run — p=1, BR2, 1000m mesh, 1800 yr
#
# Purpose: Test whether fixing the non-fault Y=0 Dirichlet BC improves
#   the recurrence interval. v33 p=1 baseline: ~435 yr (1.82x Tandem).
#   Tandem reference: ~240 yr.
#
# Fix: Non-fault Y=0 interior faces now get u_D = (Vp*t, 0, 0) instead
#   of u_D = (0, 0, 0). This matches Tandem's bp5.lua boundary().
#   Physically: locked fault (both sides move together at Vp*t).
#
# Expected: Recurrence interval closer to 240 yr if the BC mismatch
#   was artificially constraining the loading rate.

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

mkdir -p bp5/results_v34_1000m_br2_p1_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v34_1000m_br2_p1_full \
      --output-prefix bp5_v34_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
