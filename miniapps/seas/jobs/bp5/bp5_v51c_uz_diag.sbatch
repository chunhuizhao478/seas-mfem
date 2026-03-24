#!/bin/bash
#SBATCH -J bp5_v51c_uz_diag        # Job name
#SBATCH -o bp5_v51c_uz_diag_%j.out
#SBATCH -e bp5_v51c_uz_diag_%j.err
#SBATCH -p development               # Queue (fast turnaround)
#SBATCH -N 1                         # 1 node — serial for clean output
#SBATCH -n 1                         # 1 rank
#SBATCH -t 00:30:00                  # 30 minutes
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v51c: Displacement u_z diagnostic at fault faces
#
# After the first elastic solve with non-zero slip, dumps u_z at each fault face.
# For pure strike-slip loading, u_z SHOULD be exactly zero.
# If u_z != 0, the K-f assembly has a cross-component mismatch.
#
# Also dumps the dip traction components for cross-reference.
#
# Serial (1 rank) for clean, non-interleaved output.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# Build
make clean-build
make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }

mkdir -p bp5/results_v51c_uz_diag

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 100 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v51c_uz_diag \
      --output-prefix bp5_v51c_uz \
      --write-every-step \
      --diag-uz-fault \
      --diag-dip-traction \
      --dt-init 0.01
