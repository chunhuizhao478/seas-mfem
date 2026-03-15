#!/bin/bash
#SBATCH -J bp5_v33_p2_full     # Job name
#SBATCH -o bp5_v33_p2_full_%j.out
#SBATCH -e bp5_v33_p2_full_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes
#SBATCH -n 800                       # 50 tasks/node × 16 nodes
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v33: Full run — p=2 on 1000m mesh, BR2, full 1800 yr
#
# Purpose: Test whether p=2 fixes the recurrence interval mismatch.
#   v32 result at p=1: ~435 yr recurrence (1.82× Tandem's ~240 yr)
#   Expected at p=2: significantly closer to ~240 yr
#     (BR2 traction correction scales as O(h^p): O(h) at p=1, O(h²) at p=2)
#
# Mesh: bp5_tandem.msh (res_f=1, ~76K tets)
# DOFs: p=2 → 10 DOFs/tet × 3 components × 76K ≈ 2.3M
#
# Changes from v32:
#   - --order 2 (general polynomial order support)
#   - IP trace constant c_N_1 = p*(p+D-1)/D = 8/3 (was 1.0)
#   - BR2 traction correction uses multi-point face quadrature
#
# Success criteria:
#   - Recurrence interval significantly closer to 240 yr than 435 yr
#   - If recurrence ≈ 240-280 yr, the formulation is correct and the
#     mismatch was purely discretization error at p=1
#   - Compare slip profiles, stress drop with Tandem reference

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

mkdir -p bp5/results_v33_1000m_br2_p2_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v33_1000m_br2_p2_full \
      --output-prefix bp5_v33_p2 \
      --checkpoint-interval 1000 \
      --write-every-step
