#!/bin/bash
#SBATCH -J bp5_v34_p6_smoke    # Job name
#SBATCH -o bp5_v34_p6_smoke_%j.out
#SBATCH -e bp5_v34_p6_smoke_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 2                         # 2 nodes
#SBATCH -n 100                       # MPI tasks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v34: Smoke test — p=6, BR2, 4000m mesh, ~50 yr
#
# Purpose: Quick sanity check that P6 runs without blowup on the 4km mesh.
#   Verify correct initial behavior and no NaN/Inf before committing
#   to a full 1800 yr run.
#
# DOFs: P6 on 6724 tets → 1.69M DOFs

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

mkdir -p bp5/results_v34_4000m_br2_p6_smoke

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_4000m.msh \
      --mesh-scale 1000 \
      --order 6 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --tfinal 1.58e9 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v34_4000m_br2_p6_smoke \
      --output-prefix bp5_v34_p6_smoke \
      --checkpoint-interval 500 \
      --write-every-step
