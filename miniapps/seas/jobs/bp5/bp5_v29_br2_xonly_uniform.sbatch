#!/bin/bash
#SBATCH -J bp5_v29_xonly        # Job name
#SBATCH -o bp5_v29_xonly_%j.out
#SBATCH -e bp5_v29_xonly_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours
#SBATCH -A EAR20006             # TACC allocation

# BP5 v29b: Original BR2 + new mesh + x-only loading
#
# Dirichlet only on x=+-Lx (fault-normal faces), Natural everywhere else.
# This is antiplane-style loading — simpler stress field, may reduce
# artifacts at y-boundary intersections.

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

mkdir -p bp5/results_v29_br2_xonly_uniform

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --solver mumps-blr \
      --bc-mode x-only \
      --check-residual \
      --delta-tau-factor 0 \
      --V-nuc 1e-9 \
      --tfinal 9.45e9 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v29_br2_xonly_uniform \
      --output-prefix bp5_v29_xonly \
      --checkpoint-interval 1000 \
      --write-every-step
