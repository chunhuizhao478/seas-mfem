#!/bin/bash
#SBATCH -J bp5_v44a_ip_p1       # Job name
#SBATCH -o bp5_v44a_ip_p1_%j.out
#SBATCH -e bp5_v44a_ip_p1_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # MPI tasks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v44a: IP p=1, 1000m mesh — REGRESSION TEST
#
# v44 implements multi-DOF fault discretization (Tandem approach).
# At p=1: face_order=0, nbf=1 → mathematically identical to v42 face averaging.
# This run MUST reproduce v41a/v42 results: ~8 events, ~250yr recurrence.
#
# If results differ from v41a, the v44 changes broke backward compatibility.

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

mkdir -p bp5/results_v44a_1000m_ip_p1_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v44a_1000m_ip_p1_full \
      --output-prefix bp5_v44a_ip_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
