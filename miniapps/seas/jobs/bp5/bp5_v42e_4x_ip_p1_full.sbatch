#!/bin/bash
#SBATCH -J bp5_v42e_4x_ip       # Job name
#SBATCH -o bp5_v42e_4x_ip_%j.out
#SBATCH -e bp5_v42e_4x_ip_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes (147k tets, ~2.3x baseline)
#SBATCH -n 800                       # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v42e: Domain truncation — 4x domain, IP p=1, 1000m mesh, full 1800 yr
#
# Tests boundary truncation effects by quadrupling domain size:
#   Standard: X=[-200,200], Y=[-100,100], Z=[-100,0] km
#   4x:       X=[-800,800], Y=[-400,400], Z=[-400,0] km
#
# v41f (same config) reached 378yr with 2 events, recurrence ~377yr.
# c_N_1 fix is no-op at p=1.
#
# Compare with v42d (2x) to quantify convergence of boundary effects:
#   If 4x ≈ 2x → boundary effects converged at 2x, no need for larger domain
#   If 4x < 2x → standard domain (1x) has significant truncation error
#
# Mesh: bp5_tandem_4x_1000m.msh (147K tets)
# DOFs: ~1.8M

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

mkdir -p bp5/results_v42e_4x_1000m_ip_p1_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_4x_1000m.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v42e_4x_1000m_ip_p1_full \
      --output-prefix bp5_v42e_4x_ip_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
