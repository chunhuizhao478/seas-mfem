#!/bin/bash
#SBATCH -J bp5_v47i_diag          # Job name
#SBATCH -o bp5_v47i_diag_%j.out
#SBATCH -e bp5_v47i_diag_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes
#SBATCH -n 800                       # MPI tasks (50/node)
#SBATCH -t 00:30:00                  # 30 min (diagnostic — will crash quickly)
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v47i: IP p=2, 1000m mesh, ×3 penalty — FACE CONSISTENCY DIAGNOSTIC
#
# Purpose: print [SlipRHS] and [Traction] info for the first 3 fault faces
# on each rank, comparing elem1/elem2 ordering, sign, penalty, and per-DOF
# slip values between AssembleSlipContributionIP and ComputeTraction.
#
# This will blowup (same as v47b) but the diagnostic output prints BEFORE
# the crash. Check the .out file for [SlipRHS] and [Traction] lines.
#
# What to look for:
#   - Same face: [SlipRHS] and [Traction] should show same E1, E2, sign, penalty
#   - If they differ → bug in face ordering between bilinear form and RHS/traction
#   - Per-DOF slip values show the within-face heterogeneity pattern
#
# See bp5_debug_v47.md Section 18.5.

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

mkdir -p bp5/results_v47i_diag

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v47i_diag \
      --output-prefix bp5_v47i_diag \
      --checkpoint-interval 1000 \
      --write-every-step
