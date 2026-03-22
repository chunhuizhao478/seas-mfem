#!/bin/bash
#SBATCH -J bp5_v49l_par_vg       # Job name
#SBATCH -o bp5_v49l_par_vg_%j.out
#SBATCH -e bp5_v49l_par_vg_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 4                         # 4 nodes (memory headroom for 1000m mesh)
#SBATCH -n 8                         # 8 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v49l: 1000m mesh, 8 ranks, V-guard + dt=0.05s
#
# Purpose: The ULTIMATE test — 1000m parallel with both Issue A and Issue B fixes.
# Combines CFL-safe dt (Issue A) + V-guard safety net (Issue A+B) in parallel.
#
# This is the target configuration for production runs.
# 1000m mesh: 63451 tets / 8 ranks ≈ 7930 tets/rank
#
# Both fixes applied simultaneously:
#   --dt-init 0.05  → well below CFL limit (0.072s) for 1000m
#   --v-guard       → stage-level rejection (factor 100×) as safety net
#
# Expected: STABLE if both Issue A and Issue B are cascade-related.
# If crashes: structural parallel bug exists at 1000m scale.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# NO BUILD — v49j builds. Wait for binary to exist.
for i in $(seq 1 60); do
   [ -f seas_bp5_full ] && break
   echo "Waiting for v49j build... ($i/60)"
   sleep 10
done
[ ! -f seas_bp5_full ] && { echo "ERROR: seas_bp5_full not found after 10min"; exit 1; }

mkdir -p bp5/results_v49l_1000m_par_vguard

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v49l_1000m_par_vguard \
      --output-prefix bp5_v49l_par_vg \
      --write-every-step \
      --dt-init 0.05 \
      --v-guard \
      --diag-rk-stages
