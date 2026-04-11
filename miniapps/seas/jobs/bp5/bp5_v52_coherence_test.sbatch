#!/bin/bash
#SBATCH -J bp5_v52_coherence        # Job name
#SBATCH -o bp5_v52_coherence_%j.out
#SBATCH -e bp5_v52_coherence_%j.err
#SBATCH -p development               # Dev queue (2hr max, fast turnaround)
#SBATCH -N 4                         # 4 nodes
#SBATCH -n 200                       # 200 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v52: Traction Coherence Diagnostic
#
# Tests whether ComputeTraction (post-solve traction extraction) is consistent
# with the elastic solve by decomposing traction into:
#   1. Stress-only: {sigma(u).n}  -- the physical signal
#   2. Penalty correction: eta * ([[u]] - slip)  -- should be ~0
#   3. Solver fault residual: |[[u]] - slip|
#
# Tandem computes traction coherently with the solve (same kernel).
# MFEM uses a separate hand-coded ComputeTraction that re-evaluates stress
# and penalty independently. This test measures the discrepancy.
#
# The diagnostic triggers once when slip is non-trivial (after first step),
# prints a summary with RMS/max statistics, then continues the simulation.
#
# Short run: 10 years is enough to see the diagnostic fire and get early
# interseismic data for quick comparison.
#
# See bp5_debug_v52.md Section 6 for full analysis.

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

mkdir -p bp5/results_v52_coherence

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 3.15e8 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v52_coherence \
      --output-prefix bp5_v52_coherence \
      --checkpoint-interval 5000 \
      --write-every-step \
      --diag-traction-coherence
