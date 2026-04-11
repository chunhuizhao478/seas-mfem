#!/bin/bash
#SBATCH -J bp5_v50f2_p2_so        # Job name
#SBATCH -o bp5_v50f2_p2_so_%j.out
#SBATCH -e bp5_v50f2_p2_so_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 4                         # 4 nodes
#SBATCH -n 200                       # 200 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v50f2: 2500m mesh, p=2, IP, STRESS-ONLY traction (regression test)
#
# Verifies that stress-only traction does NOT break p=2 nucleation.
# Compare with v50d (same mesh, p=2, standard traction) — should match.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# NO BUILD — v50f builds. Wait for binary.
for i in $(seq 1 60); do
   [ -f seas_bp5_full ] && break
   echo "Waiting for v50f build... ($i/60)"
   sleep 10
done
[ ! -f seas_bp5_full ] && { echo "ERROR: seas_bp5_full not found after 10min"; exit 1; }

mkdir -p bp5/results_v50f2_p2_stressonly

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v50f2_p2_stressonly \
      --output-prefix bp5_v50f2_p2_so \
      --traction-stress-only \
      --write-every-step
