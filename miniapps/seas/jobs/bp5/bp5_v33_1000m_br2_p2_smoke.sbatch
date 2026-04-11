#!/bin/bash
#SBATCH -J bp5_v33_p2_smoke    # Job name
#SBATCH -o bp5_v33_p2_smoke_%j.out
#SBATCH -e bp5_v33_p2_smoke_%j.err
#SBATCH -p development               # Dev queue (max 8 nodes)
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 50 tasks/node × 8 nodes
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v33: Smoke test — p=2 on 1000m mesh, BR2, short run
#
# Purpose: Verify p=2 works on cluster before committing to full run.
#   - No crash in first ~10 years of simulation
#   - Traction values are O(15 MPa) (reasonable)
#   - Displacement is bounded
#   - MUMPS factorization fits in memory
#
# Mesh: bp5_tandem.msh (res_f=1, ~76K tets)
# DOFs: p=2 → 10 DOFs/tet × 3 components × 76K ≈ 2.3M
#   (vs p=1 → 4 DOFs/tet × 3 × 76K ≈ 912K)
#
# Changes from v32:
#   - --order 2 (general polynomial order support)
#   - IP trace constant c_N_1 = p*(p+D-1)/D = 8/3 (was 1.0)
#   - BR2 traction correction uses multi-point face quadrature
#   - 8 nodes (dev queue max); if MUMPS OOM, move full run to 16 nodes

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

mkdir -p bp5/results_v33_1000m_br2_p2_smoke

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 3.15e8 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v33_1000m_br2_p2_smoke \
      --output-prefix bp5_v33_p2_smoke \
      --checkpoint-interval 500 \
      --write-every-step
