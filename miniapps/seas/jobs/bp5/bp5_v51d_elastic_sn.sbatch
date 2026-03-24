#!/bin/bash
#SBATCH -J bp5_v51d_esn            # Job name
#SBATCH -o bp5_v51d_esn_%j.out
#SBATCH -e bp5_v51d_esn_%j.err
#SBATCH -p development               # Queue (fast turnaround)
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 ranks — production scale
#SBATCH -t 02:00:00                  # 2 hours (through first earthquake)
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v51d: Elastic sigma_n test — 1000m p=2
#
# THE KEY FIX: Use elastic normal stress (sigma_n = sigma_n_pre + T_n_elastic)
# instead of constant sigma_n = 25 MPa. Matches Tandem's DieterichRuinaAgeing.h.
#
# Expected: dip/strike ratio should drop from 0.3-5.7% to ~0.0-0.7% (Tandem level).
# The elastic T_n provides self-consistent feedback suppressing spurious dip motion.

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

mkdir -p bp5/results_v51d_elastic_sn

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
      --output-dir bp5/results_v51d_elastic_sn \
      --output-prefix bp5_v51d_esn \
      --write-every-step \
      --elastic-sigma-n \
      --diag-coseismic-dip-threshold 0.01 \
      --dt-init 0.01
