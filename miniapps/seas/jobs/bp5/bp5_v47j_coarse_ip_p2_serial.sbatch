#!/bin/bash
#SBATCH -J bp5_v47j_serial       # Job name
#SBATCH -o bp5_v47j_serial_%j.out
#SBATCH -e bp5_v47j_serial_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node
#SBATCH -n 1                         # 1 MPI rank (SERIAL — no shared faces)
#SBATCH -t 01:00:00                  # 1 hour
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v47j: IP p=2, COARSE mesh, SINGLE RANK — SHARED FACE BUG DIAGNOSTIC
#
# Purpose: test whether the multi-DOF blowup occurs when ALL fault faces
# are interior (no shared faces). With 1 MPI rank, there are no shared faces.
#
# All previous blowup DOFs were "(shared)" — never "(interior)".
# If this run is STABLE: the bug is in the shared face processing.
# If this run BLOWS UP: the bug is in the interior face path too.
#
# Mesh: bp5_tandem_coarse.msh (res=40km, res_f=10km, ~2480 tets)
# At p=2: ~74K vector DOFs — fits in single-rank MUMPS easily.
# The nucleation zone (12x12 km) is barely resolved (~1-2 fault faces),
# but we just need faces with heterogeneous V_nuc/V_init DOFs.
#
# See bp5_debug_v47.md Section 16-17.

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

mkdir -p bp5/results_v47j_coarse_serial

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_coarse.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v47j_coarse_serial \
      --output-prefix bp5_v47j_serial \
      --write-every-step
