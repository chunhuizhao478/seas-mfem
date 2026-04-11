#!/bin/bash
#SBATCH -J bp5_v47d_ip_p2_blr14 # Job name
#SBATCH -o bp5_v47d_ip_p2_blr14_%j.out
#SBATCH -e bp5_v47d_ip_p2_blr14_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes (same as v47b)
#SBATCH -n 800                       # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v47d: IP p=2, 1000m mesh — TIGHTER BLR TOLERANCE (diagnostic)
#
# Same as v47b but with BLR tolerance 1e-14 (vs default 1e-10).
# Tests whether tighter BLR factorization reduces local solver errors
# enough to prevent the penalty amplification cascade.
#
# If this works:
#   - BLR local accuracy at 1e-10 is insufficient for ×3 penalty at p=2
#   - 1e-14 adds ~50% more BLR memory and ~30% more factorization time
#   - This is a viable production fix (cheaper than exact MUMPS)
#
# If this also blows up:
#   - Need exact MUMPS (v47c) or the issue isn't BLR accuracy

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

mkdir -p bp5/results_v47d_1000m_ip_p2_blr14

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-14 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v47d_1000m_ip_p2_blr14 \
      --output-prefix bp5_v47d_ip_p2_blr14 \
      --checkpoint-interval 1000 \
      --write-every-step
