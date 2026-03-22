#!/bin/bash
#SBATCH -J bp5_v49a_smooth       # Job name
#SBATCH -o bp5_v49a_smooth_%j.out
#SBATCH -e bp5_v49a_smooth_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node (192 GB on Frontera)
#SBATCH -n 1                         # 1 rank — SERIAL
#SBATCH -t 04:00:00                  # 4 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v49a: SMOOTH NUCLEATION TEST
#
# Tests Candidate 5: sharp V_nuc/V_init boundary + multi-DOF interpolation
# creates Gibbs-like oscillations that seed the RK cascade.
#
# Change: --smooth-nucleation replaces sharp V boundary with Gaussian taper
# (taper width = 2000m = 2 elements).
#
# Also includes: --diag-normals, --diag-first-traction, --diag-rk-stages
# for baseline diagnostic data.
#
# If STABLE: Root cause is the sharp nucleation boundary. Fix: smooth init.
# If CRASH: Not the sharp boundary. Candidate 5 eliminated.

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

mkdir -p bp5/results_v49a_smooth_nuc

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
      --output-dir bp5/results_v49a_smooth_nuc \
      --output-prefix bp5_v49a_smooth \
      --write-every-step \
      --smooth-nucleation \
      --diag-normals \
      --diag-first-traction \
      --diag-rk-stages
