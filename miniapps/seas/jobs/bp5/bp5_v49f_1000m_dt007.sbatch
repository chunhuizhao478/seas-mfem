#!/bin/bash
#SBATCH -J bp5_v49f_dt007        # Job name
#SBATCH -o bp5_v49f_dt007_%j.out
#SBATCH -e bp5_v49f_dt007_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node (192 GB on Frontera)
#SBATCH -n 1                         # 1 rank — SERIAL
#SBATCH -t 02:00:00                  # 2 hours (expect stability)
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v49f: CFL FIX — dt_init = 0.07s on 1000m mesh
#
# Tests the CFL hypothesis from v49 Section 4:
#   dt_CFL = 2.0 * eta * h_min / (4*mu) = 0.072s
#   Using dt_init = 0.07s (just below the CFL limit).
#
# Expected: STABLE. Should survive first step and produce healthy nucleation.
# Compare with v49cd (same mesh, dt=0.13s → CRASH).

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# v49f is the BUILD job — always rebuild for Phase 2 code changes
make clean-build
make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }

mkdir -p bp5/results_v49f_dt007

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
      --output-dir bp5/results_v49f_dt007 \
      --output-prefix bp5_v49f_dt007 \
      --write-every-step \
      --dt-init 0.07 \
      --diag-rk-stages
