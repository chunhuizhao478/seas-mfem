#!/bin/bash
#SBATCH -J bp5_v50h_4000m_p4      # Job name
#SBATCH -o bp5_v50h_4000m_p4_%j.out
#SBATCH -e bp5_v50h_4000m_p4_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks — same as p=6 production
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v50h: 4000m mesh, p=4, IP, ClosedUniform nodes, 400 ranks
#
# DIAGNOSTIC: Same mesh as the working p=6 run (v50g p=6 prod).
# Tests whether p=4 works on 4000m mesh at 400 ranks.
#
# Tandem's p=4 reference also uses the 4000m mesh (element_size=4000m).
# This gives a direct apples-to-apples comparison with bp5qd_tandem_p4_*.txt.
#
# If WORKS: p=4 failure was 2500m mesh + 400 ranks specific (too few elements/rank)
# If FAILS: p=4 has a deeper parallel bug even with fewer shared faces

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# Build
make clean-build
make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }

mkdir -p bp5/results_v50h_4000m_p4

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_4000m.msh \
      --mesh-scale 1000 \
      --order 4 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v50h_4000m_p4 \
      --output-prefix bp5_v50h_4000m_p4 \
      --face-basis-type equi \
      --write-every-step
