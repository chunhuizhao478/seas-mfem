#!/bin/bash
#SBATCH -J bp5_v52_rhs_z_blr         # Job name
#SBATCH -o bp5_v52_rhs_z_blr_%j.out
#SBATCH -e bp5_v52_rhs_z_blr_%j.err
#SBATCH -p development               # Dev queue (2hr max, fast turnaround)
#SBATCH -N 4                         # 4 nodes
#SBATCH -n 200                       # 200 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v52: RHS z-component diagnostic — MUMPS-BLR solver
#
# Tests whether the RHS vector f has a non-zero z-component (dip forcing).
# Uses MUMPS-BLR (approximate factorization) at p=1.
#
# Companion to bp5_v52_diag_rhs_z_exact.sbatch (exact MUMPS).
# Comparing the two isolates Hypothesis A (RHS bug) vs B (solver error):
#
#   If BOTH show f_z >> 0: assembly bug (Hypothesis A)
#   If BOTH show f_z ~ 0 but exact has u_z ~ 0 and BLR has u_z >> 0:
#     BLR solver error (Hypothesis B)
#   If BOTH show same u_z: neither is BLR-specific
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

mkdir -p bp5/results_v52_rhs_z_blr

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 3.15e8 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v52_rhs_z_blr \
      --output-prefix bp5_v52_rhs_z_blr \
      --checkpoint-interval 5000 \
      --write-every-step \
      --diag-rhs-z
