#!/bin/bash
#SBATCH -J bp5_v34_p1_smoke     # Job name
#SBATCH -o bp5_v34_p1_smoke_%j.out
#SBATCH -e bp5_v34_p1_smoke_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # MPI tasks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v34: Smoke test — p=1, BR2, 1000m mesh, ~50 yr
#
# Purpose: Quick sanity check after v34 fix (non-fault Y=0 Dirichlet BC
#   changed from u=0 to u=Vp*t to match Tandem's bp5.lua boundary()).
#   Verify no blowup and correct initial behavior.
#
# Fix: Non-fault Y=0 interior faces now get u_D = (Vp*t, 0, 0) instead
#   of u_D = (0, 0, 0). This matches Tandem's locked-fault BC.

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

mkdir -p bp5/results_v34_1000m_br2_p1_smoke

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 1.58e9 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v34_1000m_br2_p1_smoke \
      --output-prefix bp5_v34_p1_smoke \
      --checkpoint-interval 500 \
      --write-every-step
