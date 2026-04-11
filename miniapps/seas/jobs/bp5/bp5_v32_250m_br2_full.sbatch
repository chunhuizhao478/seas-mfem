#!/bin/bash
#SBATCH -J bp5_v32_250m        # Job name
#SBATCH -o bp5_v32_250m_%j.out
#SBATCH -e bp5_v32_250m_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 64                        # nodes
#SBATCH -n 3200                      # MPI tasks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v32: Convergence study — 250m fault resolution, BR2, full 1800 yr
#
# Mesh: bp5_tandem_250m.msh (res_f=0.25, ~1.1M elements est.)
# Includes: SCEC output sign fix (slip/stress negated for right-lateral)
# Purpose: Fine resolution for h-convergence study

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

mkdir -p bp5/results_v32_250m_br2_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_250m.msh \
      --mesh-scale 1000 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v32_250m_br2_full \
      --output-prefix bp5_v32_250m \
      --checkpoint-interval 1000 \
      --write-every-step
