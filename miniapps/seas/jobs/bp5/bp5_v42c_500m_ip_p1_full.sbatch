#!/bin/bash
#SBATCH -J bp5_v42c_500m        # Job name
#SBATCH -o bp5_v42c_500m_%j.out
#SBATCH -e bp5_v42c_500m_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 32                        # 32 nodes (~500K tets × 12 DOFs → 6M DOFs, need more memory)
#SBATCH -n 1600                      # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v42c: h-refinement — IP p=1, 500m mesh, full 1800 yr
#
# Tests mesh resolution effect on the 20% strike slip deficit at p=1.
# Halving h should reduce O(h) stress error by 2x → deficit ~10%.
#
# v41b (same config, pre-c_N_1 fix) only reached 199yr due to walltime.
# c_N_1 fix is no-op at p=1, so physics is identical to v41b but with
# longer walltime allocation.
#
# Compare with v42a (p=2, 1000m) — both ~2x DOF increase, different convergence path.
# If p=2 gives better results with fewer DOFs, p-refinement is preferred over h-refinement.
#
# Mesh: bp5_tandem_500m.msh (~500K tets)
# DOFs: ~6M (3 components × 4 DOFs/tet × 500k tets)
# c_N_1: 1.0 (p=1, same as v41b — no change from c_N_1 fix)

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

mkdir -p bp5/results_v42c_500m_ip_p1_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_500m.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v42c_500m_ip_p1_full \
      --output-prefix bp5_v42c_500m_ip_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
