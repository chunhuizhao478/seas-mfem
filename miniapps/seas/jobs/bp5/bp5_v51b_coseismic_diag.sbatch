#!/bin/bash
#SBATCH -J bp5_v51b_codiag         # Job name
#SBATCH -o bp5_v51b_codiag_%j.out
#SBATCH -e bp5_v51b_codiag_%j.err
#SBATCH -p development               # Queue (fast turnaround)
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours (just need first earthquake)
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v51b: Coseismic dip diagnostic — baseline (NO dip zeroing)
#
# Dumps the dip/strike traction ratio during the first earthquake
# (when V_max > 0.1 m/s). This quantifies the actual cross-component
# contamination at coseismic scale, not just at the first tiny-slip step.
#
# Compare with v51a (zero dip) to see what the contamination ratio IS
# during the earthquake that creates the dip offset.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# NO BUILD — v51a builds
for i in $(seq 1 60); do
   [ -f seas_bp5_full ] && break
   echo "Waiting for v51a build... ($i/60)"
   sleep 10
done
[ ! -f seas_bp5_full ] && { echo "ERROR: seas_bp5_full not found after 10min"; exit 1; }

mkdir -p bp5/results_v51b_codiag

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
      --output-dir bp5/results_v51b_codiag \
      --output-prefix bp5_v51b_codiag \
      --write-every-step \
      --diag-coseismic-dip \
      --dt-init 0.01
