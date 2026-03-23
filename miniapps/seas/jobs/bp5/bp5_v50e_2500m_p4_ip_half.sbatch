#!/bin/bash
#SBATCH -J bp5_v50e_p4_half       # Job name
#SBATCH -o bp5_v50e_p4_half_%j.out
#SBATCH -e bp5_v50e_p4_half_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 4                         # 4 nodes
#SBATCH -n 200                       # 200 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v50e: 2500m mesh, p=4, IP with penalty_factor=0.5
#
# DIAGNOSTIC: Reduces IP penalty by half to test whether over-stiffening
# causes the p=4 nucleation failure.
#
# The penalty_factor scales BOTH the K matrix penalty AND the RHS penalty
# consistently, so K-f consistency is maintained.
#
# If nucleation works with 0.5x penalty → penalty magnitude is the issue.
# If still fails → problem is structural (not just magnitude).

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# NO BUILD — v50c builds. Wait for binary.
for i in $(seq 1 60); do
   [ -f seas_bp5_full ] && break
   echo "Waiting for v50c build... ($i/60)"
   sleep 10
done
[ ! -f seas_bp5_full ] && { echo "ERROR: seas_bp5_full not found after 10min"; exit 1; }

mkdir -p bp5/results_v50e_p4_ip_half

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 4 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v50e_p4_ip_half \
      --output-prefix bp5_v50e_p4_ip_half \
      --penalty-factor 0.5 \
      --write-every-step
