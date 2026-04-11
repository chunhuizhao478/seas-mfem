#!/bin/bash
#SBATCH -J bp5_v47h_ip_p2_unif   # Job name
#SBATCH -o bp5_v47h_ip_p2_unif_%j.out
#SBATCH -e bp5_v47h_ip_p2_unif_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes
#SBATCH -n 800                       # MPI tasks (50/node)
#SBATCH -t 01:00:00                  # 1 hour (diagnostic only — just need first ~100 steps)
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v47h: IP p=2, 1000m mesh, ×3 penalty — UNIFORM V DIAGNOSTIC
#
# Purpose: isolate whether the multi-DOF blowup is caused by within-face
# V heterogeneity (nucleation boundary cutting across faces) or by a
# fundamental bug in the multi-DOF IP coupling.
#
# Change from v47b: --V-nuc 1e-9 makes V_nuc = V_init = 1e-9 everywhere.
# No nucleation zone boundary → no within-face DOF variation in V.
# All fault DOFs start at V = 1e-9 uniformly.
#
# EXPECTED OUTCOMES (see bp5_debug_v47.md Section 16.6):
#   - If NO blowup: the bug is specifically in how within-face V variation
#     interacts with the multi-DOF coupling. Tandem handles this differently.
#   - If STILL blows up: the bug is more fundamental — in the multi-DOF
#     penalty integrator, mass matrix, or interpolation/projection itself.
#
# Note: with uniform V=1e-9, there will be no earthquake nucleation.
# This is purely a stability diagnostic. We only need ~100 steps to see
# if the simulation is stable or blows up immediately.

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

mkdir -p bp5/results_v47h_1000m_ip_p2_uniform_v

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --V-nuc 1e-9 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v47h_1000m_ip_p2_uniform_v \
      --output-prefix bp5_v47h_ip_p2_unif \
      --checkpoint-interval 1000 \
      --write-every-step
