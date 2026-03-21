#!/bin/bash
#SBATCH -J bp5_v48a_sign         # Job name
#SBATCH -o bp5_v48a_sign_%j.out
#SBATCH -e bp5_v48a_sign_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 ranks (same as v47b that crashed)
#SBATCH -t 01:00:00                  # 1 hour
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v48a: IP p=2, 1000m mesh, 400 ranks — SIGN DIAGNOSTIC ON PRODUCTION MESH
#
# The 2500m mesh with 400 ranks (v48 sign diag) did NOT crash — the RK45
# recovered from TRACTION BLOWUP. But V_max decayed (serial grows, parallel decays).
#
# This run uses the 1000m mesh where v47b crashed immediately. With the
# sign diagnostic code (including std::flush), we'll see:
#   1. Whether the 1000m mesh still crashes with current code
#   2. The [SHARED-SIGN] and [SHARED-TRAC-SIGN] output for the crashing faces
#   3. Whether the crash is on shared faces, interior faces, or both
#
# Same parameters as v47b: p=2, IP, ×3 penalty, V_nuc=0.01, MUMPS-BLR.

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

mkdir -p bp5/results_v48a_1000m_sign_diag

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v48a_1000m_sign_diag \
      --output-prefix bp5_v48a_sign \
      --checkpoint-interval 1000 \
      --write-every-step
