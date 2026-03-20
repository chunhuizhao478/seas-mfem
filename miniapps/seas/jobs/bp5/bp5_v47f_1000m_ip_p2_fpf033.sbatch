#!/bin/bash
#SBATCH -J bp5_v47f_ip_p2_fpf033  # Job name
#SBATCH -o bp5_v47f_ip_p2_fpf033_%j.out
#SBATCH -e bp5_v47f_ip_p2_fpf033_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes (p=2 → 10 DOFs/tet, ~1.9M total DOFs)
#SBATCH -n 800                       # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v47f: IP p=2, 1000m mesh — FAULT PENALTY FACTOR = 1/3 (Level 1 test)
#
# Option A from bp5_debug_v47.md Section 13:
#   - Skeleton (non-fault) faces: full ×3 penalty (correct DG accuracy)
#   - Fault faces: 1/3 × penalty (avoids within-face V amplification cascade)
#
# Level 1 implementation: fault_penalty_factor scales penalty in RHS
# (AssembleSlipContributionIP) and traction (ComputeTraction) only.
# Bilinear form K is unchanged (full penalty on all faces).
#
# At α=1/3, the fault-face penalty equals v46's value (which works at
# both p=1 and p=2). Skeleton faces keep the correct ×3. This should be
# strictly better than v46 where ALL faces had the wrong 1/3 penalty.
#
# SUCCESS CRITERIA:
#   1. No blowup at p=2 (no TRACTION BLOWUP messages)
#   2. Nucleation occurs (V_max > 1 m/s)
#   3. Results at least as good as v46b (correct skeleton improves DG accuracy)

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

mkdir -p bp5/results_v47f_1000m_ip_p2_fpf033

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --fault-penalty-factor 0.333333 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v47f_1000m_ip_p2_fpf033 \
      --output-prefix bp5_v47f_ip_p2_fpf033 \
      --checkpoint-interval 1000 \
      --write-every-step
