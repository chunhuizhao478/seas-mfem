#!/bin/bash
#SBATCH -J bp5_v51e_esn_p4        # Job name
#SBATCH -o bp5_v51e_esn_p4_%j.out
#SBATCH -e bp5_v51e_esn_p4_%j.err
#SBATCH -p normal                     # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 ranks
#SBATCH -t 04:00:00                  # 4 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v51e: Elastic sigma_n test — 4000m p=4 (Tandem's configuration)
#
# Direct comparison with Tandem p=4 reference (same mesh, same order).
# With elastic sigma_n, the dip/strike ratio should match Tandem's 0.0-0.7%.
#
# Previous v50h (constant sigma_n): dip/strike = 0.4-8.1%
# This run tests: does elastic sigma_n fix the dip offset?

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# NO BUILD — v51d builds. Wait for binary.
for i in $(seq 1 60); do
   [ -f seas_bp5_full ] && break
   echo "Waiting for v51d build... ($i/60)"
   sleep 10
done
[ ! -f seas_bp5_full ] && { echo "ERROR: seas_bp5_full not found after 10min"; exit 1; }

mkdir -p bp5/results_v51e_esn_p4

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_4000m.msh \
      --mesh-scale 1000 \
      --order 4 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 200 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v51e_esn_p4 \
      --output-prefix bp5_v51e_esn_p4 \
      --face-basis-type equi \
      --write-every-step \
      --elastic-sigma-n \
      --dt-init 0.01
