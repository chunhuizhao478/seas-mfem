#!/bin/bash
#SBATCH -J bp5_v49m_par400_vg    # Job name
#SBATCH -o bp5_v49m_par400_vg_%j.out
#SBATCH -e bp5_v49m_par400_vg_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks — PRODUCTION SCALE
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v49m: 1000m mesh, 400 ranks, V-guard + dt=0.05s — PRODUCTION SCALE
#
# Purpose: Full production-scale test. 1000m mesh at 400 ranks means:
#   63451 tets / 400 ranks ≈ 159 tets/rank
#   Many more shared faces than v49l (8 ranks)
#   Maximum shared-face exposure — the hardest test for Issue B
#
# Both fixes applied:
#   --dt-init 0.05  → below CFL limit for 1000m
#   --v-guard       → stage-level V rejection (factor 100×)
#
# If v49l (8 ranks) passes but v49m (400 ranks) fails:
#   → More shared faces amplify the parallel issue
#   → Need to investigate shared face assembly at high partition count

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

mkdir -p bp5/results_v49m_1000m_par400_vguard

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
      --output-dir bp5/results_v49m_1000m_par400_vguard \
      --output-prefix bp5_v49m_par400_vg \
      --write-every-step \
      --dt-init 0.05 \
      --v-guard \
      --diag-rk-stages
