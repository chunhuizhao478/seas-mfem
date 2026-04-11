#!/bin/bash
#SBATCH -J bp5_v37_p6_full     # Job name
#SBATCH -o bp5_v37_p6_full_%j.out
#SBATCH -e bp5_v37_p6_full_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes
#SBATCH -n 800                       # MPI tasks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v37: Full run — p=6, BR2, 4000m mesh, 1800 yr
#
# Matches Tandem's SCEC "uphoff" submission (dx=4km, P6).
# v34 p=6: NO earthquakes (fault locked for 1291 yr).
#
# v37 fixes:
#   1. face_int2 sign fix (correct Dirichlet penalty direction)
#   2. Per-quad-point {sigma.n} traction (very large effect at p=6)
#   3. Per-point BR2 correction with cross-element lifting
#
# Expected: Earthquakes should nucleate with recurrence near 240 yr.
#   At P6/4km, DG residual is negligible — this should closely match Tandem.

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

mkdir -p bp5/results_v37_4000m_br2_p6_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_4000m.msh \
      --mesh-scale 1000 \
      --order 6 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v37_4000m_br2_p6_full \
      --output-prefix bp5_v37_p6 \
      --checkpoint-interval 1000 \
      --write-every-step
