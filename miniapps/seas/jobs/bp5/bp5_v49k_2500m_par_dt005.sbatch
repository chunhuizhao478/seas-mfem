#!/bin/bash
#SBATCH -J bp5_v49k_par_dt005    # Job name
#SBATCH -o bp5_v49k_par_dt005_%j.out
#SBATCH -e bp5_v49k_par_dt005_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node
#SBATCH -n 48                        # 48 MPI ranks (same as v47m)
#SBATCH -t 01:00:00                  # 1 hour
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v49k: 2500m mesh, 48 ranks, small dt (NO V-guard)
#
# Purpose: Test whether Issue B is dt-dependent by using a very
# conservative dt_init=0.05s (well below 2500m CFL limit of 0.18s).
#
# If v49k is STABLE but v47m crashed at dt=0.13s:
#   → Issue B is also a cascade/CFL issue amplified by shared faces
# If v49k ALSO crashes:
#   → Issue B is a structural parallel bug independent of dt

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# NO BUILD — v49j builds. Wait for binary to exist.
for i in $(seq 1 60); do
   [ -f seas_bp5_full ] && break
   echo "Waiting for v49j build... ($i/60)"
   sleep 10
done
[ ! -f seas_bp5_full ] && { echo "ERROR: seas_bp5_full not found after 10min"; exit 1; }

mkdir -p bp5/results_v49k_2500m_par_dt005

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
      --output-dir bp5/results_v49k_2500m_par_dt005 \
      --output-prefix bp5_v49k_par_dt005 \
      --write-every-step \
      --dt-init 0.05 \
      --diag-rk-stages
