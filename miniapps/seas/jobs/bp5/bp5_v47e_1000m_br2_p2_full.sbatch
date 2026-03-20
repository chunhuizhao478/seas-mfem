#!/bin/bash
#SBATCH -J bp5_v47e_br2_p2       # Job name
#SBATCH -o bp5_v47e_br2_p2_%j.out
#SBATCH -e bp5_v47e_br2_p2_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes (p=2 → 10 DOFs/tet, ~1.9M total DOFs)
#SBATCH -n 800                       # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v47e: BR2 p=2, 1000m mesh — TEST BR2 AT p=2 WITH ×3 PENALTY
#
# PURPOSE: Diagnose whether the p=2 blowup is specific to IP's penalty
#   mechanism or a general p=2 issue.
#
# BR2 at p=2 differs from IP at p=2 in two critical ways:
#   1. BR2 uses fixed penalty = dim+1 = 4 (no A/V formula, no ×3 issue)
#   2. BR2 uses nbf=1 per fault face (no multi-DOF, no within-face oscillation)
#      This is because face_fe_order is gated: IP && order>=2 only (line 536)
#
# The ×3 penalty fix in the code ONLY affects IP paths (all 15 locations).
# BR2 has its own traction extraction and slip assembly — unaffected.
#
# EXPECTED OUTCOMES:
#   - If BR2 p=2 runs stably: confirms p=2 elastic volume is fine,
#     problem is IP fault coupling (within-face V amplification)
#   - If BR2 p=2 also blows up: deeper issue beyond IP penalty
#   - If BR2 p=2 runs but results poor: nbf=1 loses p=2 fault accuracy
#
# v47 changes retained (only affect IP paths, not BR2):
#   - Penalty ×3 correction at all 15 IP locations
#   - dt_init fix: uses max(V_init, V_nuc) → dt=0.13s
#
# Previous v46 changes (retained, affect both IP and BR2):
#   - Tandem initialization: V_nuc=0.01, delta_tau_factor=0, psi_init="tandem"
#
# See bp5_debug_v47.md Section 11 for context.

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

mkdir -p bp5/results_v47e_1000m_br2_p2_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v47e_1000m_br2_p2_full \
      --output-prefix bp5_v47e_br2_p2 \
      --checkpoint-interval 1000 \
      --write-every-step
