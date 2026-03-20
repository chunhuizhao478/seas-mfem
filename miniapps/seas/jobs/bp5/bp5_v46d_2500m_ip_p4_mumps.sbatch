#!/bin/bash
#SBATCH -J bp5_v46d_ip_p4_2500  # Job name
#SBATCH -o bp5_v46d_ip_p4_2500_%j.out
#SBATCH -e bp5_v46d_ip_p4_2500_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes
#SBATCH -n 800                       # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v46d: IP p=4, 2500m mesh — MUMPS-BLR direct solver
#
# Tandem reported configuration: 2500m mesh with p=4.
# Effective resolution ~2500/4 = 625m per direction.
#
# With 2500m mesh: elements scale as (1000/2500)^3 relative to 1000m ≈ 6.4x fewer.
# ~63k/6.4 ≈ 10k elements, 35 DOFs/tet → ~1M total DOFs.
# MUMPS-BLR should handle this.
#
# v46 implements multi-DOF fault discretization (Phases 1-5):
#   face_order = 4, nbf = (4+1)(4+2)/2 = 15 DOFs per fault face
#
# SUCCESS CRITERIA: V_max > 1 m/s, dip slip converges toward Tandem

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

mkdir -p bp5/results_v46d_2500m_ip_p4_mumps

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 2500 \
      --order 4 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v46d_2500m_ip_p4_mumps \
      --output-prefix bp5_v46d_ip_p4_2500 \
      --checkpoint-interval 1000 \
      --write-every-step
