#!/bin/bash
#SBATCH -J bp5_v41d_ip_p4       # Job name
#SBATCH -o bp5_v41d_ip_p4_%j.out
#SBATCH -e bp5_v41d_ip_p4_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes (13.5k tets × 105 DOFs/tet → 1.4M DOFs)
#SBATCH -n 800                       # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v41d: Match Tandem uphoff.2 — IP p=4, 2500m mesh, full 1800 yr
#
# Matches Tandem's SCEC benchmark submission "uphoff.2" (dx=2.5km, P4).
# This is the definitive convergence test — if p=4 matches Tandem, the formulation is correct.
#
# Mesh: bp5_tandem_2500m.msh (13,503 tets)
# DOFs: ~1.4M (3 components × 35 DOFs/tet × 13.5k tets)
# c_N_1: p*(p+D-1)/D = 4*(4+3-1)/3 = 8.0
#
# The combination of higher order + larger elements should give spectral-quality stress
# accuracy with manageable DOF count (similar to 1000m p=1).
#
# Previous p=4/p=6 tests (v34) failed with BR2 (fault locking, fixed in v35-v36).
# This is the first IP p=4 test — expect better behavior since IP has no lifting operator.

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

mkdir -p bp5/results_v41d_2500m_ip_p4_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 4 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v41d_2500m_ip_p4_full \
      --output-prefix bp5_v41d_ip_p4 \
      --checkpoint-interval 1000 \
      --write-every-step
