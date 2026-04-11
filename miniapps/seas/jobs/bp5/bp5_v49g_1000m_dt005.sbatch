#!/bin/bash
#SBATCH -J bp5_v49g_dt005        # Job name
#SBATCH -o bp5_v49g_dt005_%j.out
#SBATCH -e bp5_v49g_dt005_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node (192 GB on Frontera)
#SBATCH -n 1                         # 1 rank — SERIAL
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v49g: CONSERVATIVE CFL — dt_init = 0.05s on 1000m mesh
#
# Tests with extra margin below the CFL limit (0.072s).
# Using dt_init = 0.05s (z = 27.7 * 0.05 = 1.39, well within stability).
#
# Expected: STABLE. Cleanly within stability region.
# Provides confidence that the CFL analysis is correct.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# NO BUILD — v49f builds. Wait for binary to exist.
for i in $(seq 1 60); do
   [ -f seas_bp5_full ] && break
   echo "Waiting for v49f build... ($i/60)"
   sleep 10
done
[ ! -f seas_bp5_full ] && { echo "ERROR: seas_bp5_full not found after 10min"; exit 1; }

mkdir -p bp5/results_v49g_dt005

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
      --output-dir bp5/results_v49g_dt005 \
      --output-prefix bp5_v49g_dt005 \
      --write-every-step \
      --dt-init 0.05 \
      --diag-rk-stages
