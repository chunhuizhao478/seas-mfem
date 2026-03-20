#!/bin/bash
#SBATCH -J bp5_v46e_ip_p6_4k    # Job name
#SBATCH -o bp5_v46e_ip_p6_4k_%j.out
#SBATCH -e bp5_v46e_ip_p6_4k_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 32                        # 32 nodes (p=6 → 84 DOFs/tet, but coarse mesh)
#SBATCH -n 1600                      # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v46e: IP p=6, 4000m mesh — MUMPS-BLR direct solver
#
# Tandem reported configuration: 4000m mesh with p=6.
# Effective resolution ~4000/6 = 667m per direction.
# nbf = (6+1)(6+2)/2 = 28 DOFs per fault face.
#
# p=6 on tet: (6+1)(6+2)(6+3)/6 = 84 DOFs/tet.
# With 4000m mesh: ~1000 elements, 84 DOFs/tet → ~250k total DOFs.
# MUMPS-BLR should handle this. Using 32 nodes for safety.
#
# v46 implements multi-DOF fault discretization (Phases 1-5):
#   face_order = 6, nbf = 28 DOFs per fault face
#   L2 traction projection, per-quad-point slip interpolation
#   Per-DOF friction at 28 GaussLobatto nodes per face
#
# SUCCESS CRITERIA: V_max > 1 m/s, dip slip matches Tandem reference

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

mkdir -p bp5/results_v46e_4000m_ip_p6_mumps

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 4000 \
      --order 6 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v46e_4000m_ip_p6_mumps \
      --output-prefix bp5_v46e_ip_p6_4k \
      --checkpoint-interval 1000 \
      --write-every-step
