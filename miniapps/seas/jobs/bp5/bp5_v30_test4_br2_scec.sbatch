#!/bin/bash
#SBATCH -J bp5_v30_t4_br2_scec  # Job name
#SBATCH -o bp5_v30_t4_br2_scec_%j.out
#SBATCH -e bp5_v30_t4_br2_scec_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours
#SBATCH -A EAR20006             # TACC allocation

# BP5 v30 Test 4: BR2 method + Tandem coord system + SCEC initialization
#
# All v30 fixes active (see test1 for list).
# BR2 method (default) for comparison with IP.
# Default SCEC init: V_nuc=0.03, delta_tau_factor=1.0, SCEC psi
#
# This test verifies that the coordinate system change + tau_pre fix
# also works with BR2. Previous BR2 runs (v29) showed VS lockup;
# the tau_pre sign fix may resolve this.
# Run 600 years

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

mkdir -p bp5/results_v30_test4_br2_scec

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --dg-method BR2 \
      --solver mumps-blr \
      --check-residual \
      --diag-vtk \
      --tfinal 1.89e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v30_test4_br2_scec \
      --output-prefix bp5_v30_t4_br2_scec \
      --checkpoint-interval 1000 \
      --write-every-step
