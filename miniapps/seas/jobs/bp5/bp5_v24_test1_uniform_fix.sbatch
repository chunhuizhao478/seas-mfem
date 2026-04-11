#!/bin/bash
#SBATCH -J bp5_v24_t1_unifix    # Job name
#SBATCH -o bp5_v24_t1_unifix_%j.out
#SBATCH -e bp5_v24_t1_unifix_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours
#SBATCH -A EAR20006             # TACC allocation

# BP5 v24 Test 1: Uniform fault with traction fix (H26)
#
# Fix H26: Remove BR2 lifting from ComputeTraction.
# Traction = {{C:∇u}}·n̂ only (average stress, no penalty correction).
# Matches benchmark document Algorithm line 15 and Tandem.
#
# V_nuc = Vp, delta_tau = 0 (completely uniform, no nucleation).
# Purpose: Verify VS zone maintains V ≈ Vp with the traction fix.
# PASS: V/Vp stays in [0.5, 2.0] at z=22km for 300 years.
# FAIL: VS still locks up — additional issues beyond traction.

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

if [ ! -f bp5/mesh/bp5_1000m.msh ]; then
    gmsh -3 bp5/mesh/bp5.geo -o bp5/mesh/bp5_1000m.msh -setnumber h 1.0
fi

mkdir -p bp5/results_v24_test1_uniform_fix

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --solver mumps-blr \
      --bc-mode far-field \
      --check-residual \
      --delta-tau-factor 0 \
      --V-nuc 1e-9 \
      --tfinal 9.45e9 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v24_test1_uniform_fix \
      --output-prefix bp5_v24_t1_unifix \
      --checkpoint-interval 1000 \
      --write-every-step
