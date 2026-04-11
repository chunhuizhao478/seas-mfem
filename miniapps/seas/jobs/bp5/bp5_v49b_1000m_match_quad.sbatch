#!/bin/bash
#SBATCH -J bp5_v49b_quad         # Job name
#SBATCH -o bp5_v49b_quad_%j.out
#SBATCH -e bp5_v49b_quad_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node (192 GB on Frontera)
#SBATCH -n 1                         # 1 rank — SERIAL
#SBATCH -t 04:00:00                  # 4 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v49b: MATCH QUADRATURE ORDER TEST
#
# Tests Candidate 1: quadrature order mismatch between K (2p) and f (2p+1).
#
# Change: --match-quad-order forces RHS and traction to use 2p (same as K)
# instead of 2p+1.
#
# Also includes: --diag-normals, --diag-first-traction, --diag-rk-stages
#
# If STABLE: Quadrature order IS the issue (surprising but possible).
# If CRASH: Candidate 1 eliminated.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# Build already done by v49a (same code). Skip if binary exists.
if [ ! -f seas_bp5_full ]; then
   make clean-build
   make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }
fi

mkdir -p bp5/results_v49b_match_quad

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
      --output-dir bp5/results_v49b_match_quad \
      --output-prefix bp5_v49b_quad \
      --write-every-step \
      --match-quad-order \
      --diag-normals \
      --diag-first-traction \
      --diag-rk-stages
