#!/bin/bash
#SBATCH -J bp5_v42a_ip_p2       # Job name
#SBATCH -o bp5_v42a_ip_p2_%j.out
#SBATCH -e bp5_v42a_ip_p2_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes (p=2 → 10 DOFs/tet × 3 comp → 1.9M DOFs)
#SBATCH -n 800                       # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v42a: IP p=2, 1000m mesh, full 1800 yr — c_N_1 fix verification
#
# v42 fixes the c_N_1 hardcoded 1.0 bug in all 5 IP RHS assembly paths.
# At p=2, c_N_1 = 2.667. v41c (same config, bug present) produced ZERO earthquakes
# because the RHS penalty was 2.67x too weak vs the bilinear form.
#
# This is THE critical test: if p=2 nucleates earthquakes with ~240yr recurrence,
# the c_N_1 fix is confirmed and the formulation is verified correct.
#
# Expected improvements over p=1:
#   - Strike slip deficit: ~4% (was ~20% at p=1)
#   - Dip quantities: ~0 (confirmed near-zero in v41c before locking)
#   - Recurrence: ~240yr (matching Tandem)
#
# Mesh: bp5_tandem.msh (63,451 tets, same as v41a baseline)
# DOFs: ~1.9M (3 components × 10 DOFs/tet × 63k tets)
# c_N_1: p*(p+D-1)/D = 2*(2+3-1)/3 = 8/3 ≈ 2.67

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

mkdir -p bp5/results_v42a_1000m_ip_p2_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v42a_1000m_ip_p2_full \
      --output-prefix bp5_v42a_ip_p2 \
      --checkpoint-interval 1000 \
      --write-every-step
