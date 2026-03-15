#!/bin/bash
#SBATCH -J bp5_v31_t2_ip_scec   # Job name
#SBATCH -o bp5_v31_t2_ip_scec_%j.out
#SBATCH -e bp5_v31_t2_ip_scec_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours
#SBATCH -A EAR20006             # TACC allocation

# BP5 v31 Test 2: IP method + SCEC initialization (nucleation)
#
# What's new in v31 (on top of v30):
#   - Fixed DG slip sign: (nor(1)>0) ? -1 : 1  -->  ? +1 : -1  (6 locations)
#     Root cause of all v30 blowups: sign*delta_u had wrong direction after
#     v30's tau_pre negation made delta_u < 0 for right-lateral.
#   - Added Dirichlet interior face handling (Y=0 attr-5 faces outside fault)
#
# Default SCEC init: V_nuc=0.03, delta_tau_factor=1.0, SCEC psi
# Expected: immediate earthquake from nucleation zone, then cycling
# PASS criterion: survives t > 1e8 s, traction O(15 MPa), earthquake nucleation
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

mkdir -p bp5/results_v31_test2_ip_scec

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --diag-vtk \
      --tfinal 1.89e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v31_test2_ip_scec \
      --output-prefix bp5_v31_t2_ip_scec \
      --checkpoint-interval 1000 \
      --write-every-step
