#!/bin/bash
#SBATCH -J bp5_v42d_2x_ip       # Job name
#SBATCH -o bp5_v42d_2x_ip_%j.out
#SBATCH -e bp5_v42d_2x_ip_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 12                        # 12 nodes (104k tets, ~1.6x baseline)
#SBATCH -n 600                       # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v42d: Domain truncation — 2x domain, IP p=1, 1000m mesh, full 1800 yr
#
# Tests boundary truncation effects by doubling domain size in all directions:
#   Standard: X=[-200,200], Y=[-100,100], Z=[-100,0] km
#   2x:       X=[-400,400], Y=[-200,200], Z=[-200,0] km
#
# v41e (same config, pre-c_N_1 fix) reached 926yr with 3 events, recurrence ~333yr.
# c_N_1 fix is no-op at p=1, so physics is identical to v41e.
# Longer walltime may allow reaching steady-state cycling.
#
# Key questions:
#   1. Does longer recurrence persist (boundary stiffening at 1x)?
#   2. Are dip quantities reduced vs 1x? (v41e showed ~50% reduction)
#   3. Is the 20% strike deficit affected by domain size?
#
# Mesh: bp5_tandem_2x_1000m.msh (104K tets)
# DOFs: ~1.2M

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

mkdir -p bp5/results_v42d_2x_1000m_ip_p1_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2x_1000m.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v42d_2x_1000m_ip_p1_full \
      --output-prefix bp5_v42d_2x_ip_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
