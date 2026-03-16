#!/bin/bash
#SBATCH -J bp5_v34_p6_full     # Job name
#SBATCH -o bp5_v34_p6_full_%j.out
#SBATCH -e bp5_v34_p6_full_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes
#SBATCH -n 800                       # MPI tasks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v34: Full run — p=6, BR2, 4000m mesh, 1800 yr
#
# Purpose: Match Tandem's SCEC "uphoff" submission (dx=4km, P6).
#   Test whether high polynomial order resolves the recurrence interval
#   and dip stress mismatch. Tandem at P6 gets excellent benchmark results.
#
# DOFs: P6 on 6724 tets → 84 DOFs/tet × 3 components = 1.69M DOFs
#   (comparable to our 1km P2 run at 1.90M DOFs)
#
# Key physics: BR2 traction correction scales as O(h^p).
#   At P6 with h=4km: O(4^6) suppression → DG residual negligible.
#   This is why Tandem's solution is accurate on coarse meshes.

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

mkdir -p bp5/results_v34_4000m_br2_p6_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_4000m.msh \
      --mesh-scale 1000 \
      --order 6 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v34_4000m_br2_p6_full \
      --output-prefix bp5_v34_p6 \
      --checkpoint-interval 1000 \
      --write-every-step
