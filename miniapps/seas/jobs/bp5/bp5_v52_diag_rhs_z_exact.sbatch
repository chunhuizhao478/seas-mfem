#!/bin/bash
#SBATCH -J bp5_v52_rhs_z_exact       # Job name
#SBATCH -o bp5_v52_rhs_z_exact_%j.out
#SBATCH -e bp5_v52_rhs_z_exact_%j.err
#SBATCH -p development               # Dev queue (2hr max, fast turnaround)
#SBATCH -N 4                         # 4 nodes
#SBATCH -n 200                       # 200 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v52: RHS z-component diagnostic — Exact MUMPS solver
#
# Same diagnostic as bp5_v52_diag_rhs_z.sbatch but with exact MUMPS
# (no BLR approximation). This matches what Tandem uses on Frontera:
#   -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps
#
# At p=1 the DOF count is ~8x smaller than p=2, so exact MUMPS
# should fit in memory (v47c OOM was at p=2).
#
# Key comparison with BLR run:
#   If exact MUMPS has u_z ~ 0 but BLR has u_z >> 0:
#     BLR factorization error is the root cause (Hypothesis B)
#   If both have u_z >> 0:
#     The problem is in the RHS assembly (Hypothesis A)
#
# See bp5_debug_v52.md Section 9 for full analysis.

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

mkdir -p bp5/results_v52_rhs_z_exact

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps \
      --check-residual \
      --tfinal 3.15e8 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v52_rhs_z_exact \
      --output-prefix bp5_v52_rhs_z_exact \
      --checkpoint-interval 5000 \
      --write-every-step \
      --diag-rhs-z
