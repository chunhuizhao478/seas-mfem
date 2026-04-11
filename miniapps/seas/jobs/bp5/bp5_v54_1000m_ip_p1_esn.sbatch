#!/bin/bash
#SBATCH -J bp5_v54_esn_p1            # Job name
#SBATCH -o bp5_v54_esn_p1_%j.out
#SBATCH -e bp5_v54_esn_p1_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v54: IP p=1, 1000m mesh — elastic sigma_n ON by default
#
# v54 cross-verification finding:
#   Tandem ALWAYS uses elastic normal stress feedback:
#     sigma_n = -sn_elastic + SnPre  (DieterichRuinaBase.h:87)
#   MFEM previously used constant sigma_n = 25 MPa by default.
#
#   Unit test (Test 14) showed the DG elastic normal traction is up to
#   9.23 MPa (37% of sigma_n) on test meshes. This creates a 5.54 MPa
#   traction effect that Tandem includes but MFEM was ignoring.
#
# Code change:
#   elastic_sigma_n_ now defaults to true in seas_operator.hpp (v54)
#   Use --no-elastic-sigma-n to revert to old constant sigma_n behavior
#
# Other settings:
#   - Same mesh/solver/resolution as v53
#   - --zero-dip-traction ON (v52 validated)
#   - No --write-every-step (adaptive output for performance)
#   - MUMPS-BLR solver (CG proved identical, v54)

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

mkdir -p bp5/results_v54_1000m_ip_p1_esn

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --zero-dip-traction \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v54_1000m_ip_p1_esn \
      --output-prefix bp5_v54_esn_p1 \
      --checkpoint-interval 1000
