#!/bin/bash
#SBATCH -J bp5_v24_t2_tfix      # Job name
#SBATCH -o bp5_v24_t2_tfix_%j.out
#SBATCH -e bp5_v24_t2_tfix_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours
#SBATCH -A EAR20006             # TACC allocation

# BP5 v24 Test 2: Full BP5 with traction fix (H26) + Tandem-matching params
#
# Fix H26: Remove BR2 lifting from ComputeTraction.
# Uses Tandem-matching initialization: V_nuc=0.01, psi-init-mode=tandem.
#
# Expected: first event ~150yr, deep VS maintains V ≈ Vp, subsequent cycling.
# Run 600 years to capture potential cycling.

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

mkdir -p bp5/results_v24_test2_tandem_fix

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --solver mumps-blr \
      --bc-mode far-field \
      --check-residual \
      --psi-init-mode tandem \
      --V-nuc 0.01 \
      --tfinal 1.89e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v24_test2_tandem_fix \
      --output-prefix bp5_v24_t2_tfix \
      --checkpoint-interval 1000 \
      --write-every-step
