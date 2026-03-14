#!/bin/bash
#SBATCH -J bp5_v21_med_mblr     # Job name
#SBATCH -o bp5_v21_med_mblr_%j.out
#SBATCH -e bp5_v21_med_mblr_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours (medium run)
#SBATCH -A EAR20006             # TACC allocation

# BP5 v21: MUMPS BLR medium run with FarField BC fix (H25)
# Fix: attrs 1-4 Dirichlet, attrs 5-6 Natural (matching Tandem)
# Previous runs had all-Dirichlet BCs which over-constrained top/bottom

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

if [ ! -f bp5/mesh/bp5_1000m.msh ]; then
    gmsh -3 bp5/mesh/bp5.geo -o bp5/mesh/bp5_1000m.msh -setnumber h 1.0
fi

mkdir -p bp5/results_v21_med_mumps_blr

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --solver mumps-blr \
      --bc-mode far-field \
      --check-residual \
      --tfinal 9.45e9 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v21_med_mumps_blr \
      --output-prefix bp5_v21_med_mblr \
      --checkpoint-interval 1000 \
      --write-every-step
