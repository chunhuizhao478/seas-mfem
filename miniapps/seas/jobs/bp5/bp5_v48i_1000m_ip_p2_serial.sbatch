#!/bin/bash
#SBATCH -J bp5_v48i_1k_ser       # Job name
#SBATCH -o bp5_v48i_1k_ser_%j.out
#SBATCH -e bp5_v48i_1k_ser_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node (192 GB on Frontera)
#SBATCH -n 1                         # 1 rank — SERIAL, no shared faces
#SBATCH -t 04:00:00                  # 4 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v48i: IP p=2, 1000m mesh, SERIAL (1 rank), BLR tol=1e-12
#
# CRITICAL TEST: We have NEVER tested the 1000m mesh at p=2 in serial.
#
# What we know:
#   - 2500m serial (v47l): STABLE, V grows (nucleation) — p=2 works on this mesh
#   - 1000m serial: NEVER TESTED at p=2
#   - 1000m 8-rank (v48h): CRASH at DOF 1166, x=(-500,0,-2342) — interior DOF,
#     NOT in nucleation zone, NOT a shared face issue
#   - Same crash with exact MUMPS (v48f) and BLR 1e-12 (v48h)
#   - p=1 at 400 ranks on 1000m: stable (no crash, v47a)
#
# If STABLE: the 1000m mesh works in serial at p=2, crash is parallel-specific
# If CRASHES at same DOF 1166: the 1000m mesh has a mesh/element issue at p=2
#   that exists regardless of MPI (not a shared face or solver issue)
#
# Memory: 63451 elements × 10 DOFs × 3 = ~1.9M vector DOFs.
# MUMPS-BLR on 1 rank should fit in 192 GB.

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

mkdir -p bp5/results_v48i_1000m_serial

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v48i_1000m_serial \
      --output-prefix bp5_v48i_1k_ser \
      --write-every-step
