#!/bin/bash
#SBATCH -J bp5_v23_t1_uniform   # Job name
#SBATCH -o bp5_v23_t1_uniform_%j.out
#SBATCH -e bp5_v23_t1_uniform_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours
#SBATCH -A EAR20006             # TACC allocation

# BP5 v23 Test 1-corrected: Completely uniform fault (no nucleation)
#
# V_nuc = Vp = 1e-9 everywhere, delta_tau = 0.
# Every fault node starts at V = Vp with SCEC psi.
# No nucleation zone, no perturbation of any kind.
#
# Purpose: Verify VS zone maintains V ≈ Vp in the absence of any perturbation.
# PASS: V/Vp stays in [0.5, 2.0] everywhere for 300 years.
# FAIL: V drifts — indicates fundamental elasticity/loading problem.

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

if [ ! -f bp5/mesh/bp5_1000m.msh ]; then
    gmsh -3 bp5/mesh/bp5.geo -o bp5/mesh/bp5_1000m.msh -setnumber h 1.0
fi

mkdir -p bp5/results_v23_test1_uniform

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --solver mumps-blr \
      --bc-mode far-field \
      --check-residual \
      --delta-tau-factor 0 \
      --V-nuc 1e-9 \
      --tfinal 9.45e9 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v23_test1_uniform \
      --output-prefix bp5_v23_t1_uniform \
      --checkpoint-interval 1000 \
      --write-every-step
