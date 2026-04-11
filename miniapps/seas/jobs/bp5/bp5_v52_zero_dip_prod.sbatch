#!/bin/bash
#SBATCH -J bp5_v52_zero_dip        # Job name
#SBATCH -o bp5_v52_zero_dip_%j.out
#SBATCH -e bp5_v52_zero_dip_%j.err
#SBATCH -p normal                    # Queue (48hr max)
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v52: Zero dip traction — PRODUCTION run
#
# The benchmark spec (Eq. 15b) explicitly states V_3 = 0 (zero dip velocity)
# outside Omega_f, and the initial condition (Eq. 16) sets V_3 = V_zero = 1e-20.
# For a pure strike-slip benchmark on a planar vertical fault, tau_dip = 0
# is the analytically correct value everywhere at all times. Even Tandem's
# 0.2-0.7% dip/strike ratio is a discretization artifact.
#
# The --zero-dip-traction flag sets tau_dip = 0 after ComputeTraction,
# eliminating the 21% cross-component DG contamination measured in v51c.
# This is benchmark-correct, not a workaround.
#
# Test: Does removing spurious tau_dip fix both:
#   1. Dip offset (should go to zero)
#   2. Interseismic strike deviations (tau_strike ~1 MPa low, timing drift)
#
# Configuration: p=2, 1000m, IP, BLR 1e-12, v50 defaults, zero-dip-traction
#
# See bp5_debug_v52.md for full analysis and expected outcomes.

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

mkdir -p bp5/results_v52_zero_dip_prod

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
      --output-dir bp5/results_v52_zero_dip_prod \
      --output-prefix bp5_v52_zero_dip \
      --checkpoint-interval 1000 \
      --write-every-step \
      --zero-dip-traction
