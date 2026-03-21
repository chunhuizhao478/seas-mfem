#!/bin/bash
#SBATCH -J bp5_v48_sign          # Job name
#SBATCH -o bp5_v48_sign_%j.out
#SBATCH -e bp5_v48_sign_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node
#SBATCH -n 4                         # 4 ranks (small count to identify shared faces)
#SBATCH -t 00:30:00                  # 30 min
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v48 sign diagnostic: verify whether CalcOrtho produces the same or
# opposite normal on shared fault faces across MPI ranks.
#
# Output lines [SHARED-SIGN] and [SHARED-TRAC-SIGN] show:
#   rank, shared face index, Elem1No, Elem2No, normal vector, sign, face center
#
# KEY CHECK: find two [SHARED-SIGN] lines with the SAME face_center from
# DIFFERENT ranks. Compare their nor and sign values:
#   - If nor and sign are THE SAME on both ranks → v48 analysis correct (sign bug)
#   - If nor and sign are OPPOSITE on both ranks → sign flips, code may be correct,
#     bug is elsewhere
#
# Using 2500m mesh (blows up in parallel), 4 ranks for clarity.

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

mkdir -p bp5/results_v48_sign_diag

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v48_sign_diag \
      --output-prefix bp5_v48_sign \
      --write-every-step
