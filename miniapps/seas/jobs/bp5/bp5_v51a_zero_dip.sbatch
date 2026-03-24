#!/bin/bash
#SBATCH -J bp5_v51a_zero_dip       # Job name
#SBATCH -o bp5_v51a_zero_dip_%j.out
#SBATCH -e bp5_v51a_zero_dip_%j.err
#SBATCH -p development               # Queue (fast turnaround)
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours (confirm fix, then extend)
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v51a: Zero dip traction — production run
#
# Sets tau_dip = 0 after ComputeTraction, before friction law.
# Tests: does eliminating the 21% dip contamination fix the dip offset
# and improve strike-slip match with Tandem?
#
# Also includes coseismic dip diagnostic (dumps ratio during earthquake).

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# BUILD job
make clean-build
make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }

mkdir -p bp5/results_v51a_zero_dip

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 200 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v51a_zero_dip \
      --output-prefix bp5_v51a_zero_dip \
      --write-every-step \
      --zero-dip-traction \
      --diag-coseismic-dip \
      --dt-init 0.01
