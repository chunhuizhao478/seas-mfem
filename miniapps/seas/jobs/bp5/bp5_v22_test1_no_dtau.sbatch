#!/bin/bash
#SBATCH -J bp5_v22_t1_nodtau    # Job name
#SBATCH -o bp5_v22_t1_nodtau_%j.out
#SBATCH -e bp5_v22_t1_nodtau_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours
#SBATCH -A EAR20006             # TACC allocation

# BP5 v22 Test 1: No nucleation perturbation (delta_tau = 0)
#
# Purpose: Check if deep VS zone maintains V ≈ Vp without earthquake disruption.
# If VS zone stays at plate rate → premature earthquake is the sole cause of lockup.
# If VS zone drifts below Vp → there is a separate loading problem.
#
# Expected: no earthquake, stress should build at ~0.07 MPa/yr from boundary+deep creep.

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

mkdir -p bp5/results_v22_test1_no_dtau

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --solver mumps-blr \
      --bc-mode far-field \
      --check-residual \
      --delta-tau-factor 0 \
      --tfinal 9.45e9 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v22_test1_no_dtau \
      --output-prefix bp5_v22_t1_nodtau \
      --checkpoint-interval 1000 \
      --write-every-step
