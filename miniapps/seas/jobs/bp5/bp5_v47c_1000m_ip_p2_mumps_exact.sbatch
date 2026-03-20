#!/bin/bash
#SBATCH -J bp5_v47c_ip_p2_exact  # Job name
#SBATCH -o bp5_v47c_ip_p2_exact_%j.out
#SBATCH -e bp5_v47c_ip_p2_exact_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 32                        # 32 nodes (2× v47b for exact MUMPS memory)
#SBATCH -n 800                       # MPI tasks (25/node → 2× memory per task)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v47c: IP p=2, 1000m mesh — EXACT MUMPS (diagnostic)
#
# HYPOTHESIS: v47b (p=2, ×3 penalty, MUMPS-BLR) blows up because:
#   1. BLR approximate factorization introduces local solver errors
#   2. ×3 penalty amplifies local [u]-slip errors: δT = η × ε
#   3. Rate-state friction exponentially amplifies δT: V ∝ exp(Δτ/(a×σ_n))
#   4. Cascade through RK stages → blowup
#
# TEST: Use exact MUMPS (no BLR) to eliminate solver error contribution.
# If this runs successfully, the hypothesis is confirmed.
# Tandem uses exact MUMPS (-pc_type lu -pc_factor_mat_solver_type mumps).
#
# MEMORY NOTE: Exact MUMPS requires ~2-4× more memory than BLR.
# Using 32 nodes (vs 16 for v47b) to provide sufficient memory headroom.
# Each rank gets ~7.7 GB (vs ~3.8 GB with 16 nodes).

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

mkdir -p bp5/results_v47c_1000m_ip_p2_mumps_exact

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v47c_1000m_ip_p2_mumps_exact \
      --output-prefix bp5_v47c_ip_p2_exact \
      --checkpoint-interval 1000 \
      --write-every-step
