#!/bin/bash
#SBATCH -J bp5_v34_p4_smoke    # Job name
#SBATCH -o bp5_v34_p4_smoke_%j.out
#SBATCH -e bp5_v34_p4_smoke_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 2                         # 2 nodes
#SBATCH -n 100                       # MPI tasks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v34: Smoke test — p=4, BR2, 2500m mesh, ~50 yr
#
# Purpose: Quick sanity check that P4 runs without blowup on the 2.5km mesh.
#
# DOFs: P4 on 13503 tets → 1.42M DOFs

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

mkdir -p bp5/results_v34_2500m_br2_p4_smoke

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 4 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 1.58e9 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v34_2500m_br2_p4_smoke \
      --output-prefix bp5_v34_p4_smoke \
      --checkpoint-interval 500 \
      --write-every-step
