#!/bin/bash
#SBATCH -J bp5_v42b_ip_p4       # Job name
#SBATCH -o bp5_v42b_ip_p4_%j.out
#SBATCH -e bp5_v42b_ip_p4_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes (13.5k tets × 105 DOFs/tet → 1.4M DOFs)
#SBATCH -n 800                       # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v42b: IP p=4, 2500m mesh, full 1800 yr — match Tandem uphoff.2
#
# v42 fixes the c_N_1 hardcoded 1.0 bug. At p=4, c_N_1 = 8.0.
# v41d (same config, bug present) had RHS penalty 8x too weak → total fault locking.
#
# This matches Tandem's SCEC benchmark submission "uphoff.2" (dx=2.5km, P4).
# If p=4 matches Tandem, the formulation is definitively correct and the p=1
# strike slip deficit / dip artifacts are confirmed as pure discretization error.
#
# Expected:
#   - Recurrence: ~240yr (matching Tandem)
#   - Strike slip deficit: <1% (spectral convergence)
#   - Dip quantities: essentially zero
#
# Mesh: bp5_tandem_2500m.msh (13,503 tets)
# DOFs: ~1.4M (3 components × 35 DOFs/tet × 13.5k tets)
# c_N_1: p*(p+D-1)/D = 4*(4+3-1)/3 = 8.0

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

mkdir -p bp5/results_v42b_2500m_ip_p4_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 4 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v42b_2500m_ip_p4_full \
      --output-prefix bp5_v42b_ip_p4 \
      --checkpoint-interval 1000 \
      --write-every-step
