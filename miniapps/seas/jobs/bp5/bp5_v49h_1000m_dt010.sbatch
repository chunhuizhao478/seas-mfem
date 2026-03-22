#!/bin/bash
#SBATCH -J bp5_v49h_dt010        # Job name
#SBATCH -o bp5_v49h_dt010_%j.out
#SBATCH -e bp5_v49h_dt010_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node (192 GB on Frontera)
#SBATCH -n 1                         # 1 rank — SERIAL
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v49h: PETSc DEFAULT — dt_init = 0.10s on 1000m mesh
#
# Tests with PETSc's hardcoded default (tscreate.c:45: ts->time_step = 0.1).
# z = 27.7 * 0.10 = 2.77 — marginally within stability region.
#
# Expected: MARGINALLY STABLE. First step may have large amplification
# (similar to 2500m at dt=0.13), but should survive and adaptive controller
# will reduce dt for subsequent steps.
#
# This tells us if Tandem's stability comes from its dt=0.1 default.

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

mkdir -p bp5/results_v49h_dt010

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
      --output-dir bp5/results_v49h_dt010 \
      --output-prefix bp5_v49h_dt010 \
      --write-every-step \
      --dt-init 0.10 \
      --diag-rk-stages
