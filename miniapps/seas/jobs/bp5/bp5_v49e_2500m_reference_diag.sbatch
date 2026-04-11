#!/bin/bash
#SBATCH -J bp5_v49e_2500ref      # Job name
#SBATCH -o bp5_v49e_2500ref_%j.out
#SBATCH -e bp5_v49e_2500ref_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node (192 GB on Frontera)
#SBATCH -n 1                         # 1 rank — SERIAL
#SBATCH -t 04:00:00                  # 4 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v49e: REFERENCE DIAGNOSTIC ON 2500m (should be stable)
#
# Same diagnostics as v49cd but on the 2500m mesh, which is KNOWN STABLE
# at p=2 with ×3 penalty in serial (v47l).
#
# Purpose: compare diagnostic output (normals, seed traction, RK stages)
# between the stable 2500m and crashing 1000m to find what's different.
#
# Diagnostics: --diag-normals, --diag-first-traction, --diag-rk-stages

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

mkdir -p bp5/results_v49e_2500m_ref

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 2500 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v49e_2500m_ref \
      --output-prefix bp5_v49e_2500ref \
      --write-every-step \
      --diag-normals \
      --diag-first-traction \
      --diag-rk-stages
