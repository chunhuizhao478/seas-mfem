#!/bin/bash
#SBATCH -J bp5_v47k_4rank        # Job name
#SBATCH -o bp5_v47k_4rank_%j.out
#SBATCH -e bp5_v47k_4rank_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node
#SBATCH -n 4                         # 4 MPI ranks (creates shared faces)
#SBATCH -t 01:00:00                  # 1 hour
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v47k: IP p=2, COARSE mesh, 4 RANKS — CONFIRM SHARED FACE BLOWUP
#
# v47j (1 rank, same mesh) was STABLE — 390 steps, no blowup.
# This run uses 4 ranks to introduce shared faces. If it blows up,
# confirms the bug is in shared face processing.
# If it's still stable, the coarse mesh doesn't place nucleation boundary
# faces at partition boundaries — try more ranks.
#
# Same mesh: bp5_tandem_coarse.msh (res=40km, res_f=10km, ~2480 tets)
# Same parameters: p=2, IP, ×3 penalty, V_nuc=0.01

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

mkdir -p bp5/results_v47k_coarse_4rank

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_coarse.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v47k_coarse_4rank \
      --output-prefix bp5_v47k_4rank \
      --write-every-step
