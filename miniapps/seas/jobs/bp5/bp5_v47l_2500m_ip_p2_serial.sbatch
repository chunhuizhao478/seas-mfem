#!/bin/bash
#SBATCH -J bp5_v47l_serial       # Job name
#SBATCH -o bp5_v47l_serial_%j.out
#SBATCH -e bp5_v47l_serial_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node
#SBATCH -n 1                         # 1 MPI rank (SERIAL — no shared faces)
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v47l: IP p=2, 2500m mesh, SINGLE RANK — SHARED FACE BUG DIAGNOSTIC
#
# The 2500m mesh blows up at 400 ranks (v47g at p=4 was immediate blowup).
# At p=2: ~13500 tets × 10 DOFs × 3 components = ~405K DOFs.
# Fits in single-rank MUMPS memory.
#
# If STABLE at 1 rank: shared face bug confirmed.
# If BLOWS UP at 1 rank: bug is in interior face path (not shared-specific).
#
# Note: the coarse mesh (v47j/k) didn't reproduce the blowup at any rank
# count because the nucleation zone (12x12 km) is barely resolved at 10km.
# The 2500m mesh has ~5 fault faces across the nucleation zone boundary,
# sufficient to trigger the heterogeneous-DOF instability.

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

mkdir -p bp5/results_v47l_2500m_serial

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2500m.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v47l_2500m_serial \
      --output-prefix bp5_v47l_serial \
      --write-every-step
