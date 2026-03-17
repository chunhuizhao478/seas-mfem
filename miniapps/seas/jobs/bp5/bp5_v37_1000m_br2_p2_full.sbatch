#!/bin/bash
#SBATCH -J bp5_v37_p2_full     # Job name
#SBATCH -o bp5_v37_p2_full_%j.out
#SBATCH -e bp5_v37_p2_full_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes
#SBATCH -n 800                       # 50 tasks/node x 16 nodes
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v37: Full run — p=2, BR2, 1000m mesh, 1800 yr
#
# Combined fixes:
#   1. face_int2 sign fix (drives both elements correctly)
#   2. Per-quad-point {sigma.n} traction (correct for p>=2)
#   3. Per-point BR2 correction with cross-element lifting
#
# v34 p=2 showed ~396 yr recurrence (worse than p=1). With all three
# fixes, p=2 should converge toward Tandem's ~240 yr reference.

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

mkdir -p bp5/results_v37_1000m_br2_p2_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v37_1000m_br2_p2_full \
      --output-prefix bp5_v37_p2 \
      --checkpoint-interval 1000 \
      --write-every-step
