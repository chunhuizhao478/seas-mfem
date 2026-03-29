#!/bin/bash
#SBATCH -J bp5_v52_rhs_z             # Job name
#SBATCH -o bp5_v52_rhs_z_%j.out
#SBATCH -e bp5_v52_rhs_z_%j.err
#SBATCH -p development               # Dev queue (2hr max, fast turnaround)
#SBATCH -N 4                         # 4 nodes
#SBATCH -n 200                       # 200 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v52: RHS z-component diagnostic — Isolate Hypothesis A vs B
#
# Tests whether the RHS vector f has a non-zero z-component (dip forcing).
# This is the key diagnostic to distinguish:
#
#   Hypothesis A (RHS bug): f_z != 0 => assembly code (EmbedSlip, Dirichlet
#     loading, or tensor coupling at boundaries) creates z-forcing that
#     forces u_z != 0 regardless of solver accuracy.
#
#   Hypothesis B (solver error): f_z = 0 but u_z != 0 => the BLR solver
#     introduces spurious u_z through factorization error at high-penalty
#     fault face blocks.
#
# The diagnostic fires once after the first non-trivial slip step and prints:
#   1. Per-component L2 norms of f_slip (slip contribution to RHS)
#   2. Per-component L2 norms of f_diri (Dirichlet contribution to RHS)
#   3. Per-component L2 norms of total RHS
#   4. Solution u decomposition: ||u_x||, ||u_y||, ||u_z||
#
# If ||f_z||/||f_x|| >> 0: Hypothesis A confirmed (assembly bug)
# If ||f_z||/||f_x|| ~ 0 but ||u_z||/||u_x|| >> 0: Hypothesis B (solver)
#
# Short run: 10 years is enough to trigger the diagnostic and capture
# early interseismic data for comparison.
#
# See bp5_debug_v52.md Section 9.1 for full analysis.

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

mkdir -p bp5/results_v52_rhs_z

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
      --output-dir bp5/results_v52_rhs_z \
      --output-prefix bp5_v52_rhs_z \
      --checkpoint-interval 5000 \
      --write-every-step \
      --diag-rhs-z
