#!/bin/bash
#SBATCH -J bp5_v22_t2_bdrload   # Job name
#SBATCH -o bp5_v22_t2_bdrload_%j.out
#SBATCH -e bp5_v22_t2_bdrload_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 00:30:00             # 30 minutes (diagnostic only)
#SBATCH -A EAR20006             # TACC allocation

# BP5 v22 Test 2: Boundary loading traction verification
#
# Purpose: Verify the DG Dirichlet boundary loading produces correct traction.
# Uses --diag-vtk to solve at t=1yr with zero slip and output traction at fault.
#
# Analytical estimate: tau = mu*Vp*1yr/(2*Lx) ≈ 5050 Pa ≈ 0.005 MPa
# Check: Is traction at fault center close to this? Is it purely along-strike?
#
# This is a short diagnostic run — exits after VTK output, no time stepping.

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

mkdir -p bp5/results_v22_test2_bdr_loading

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --solver mumps-blr \
      --bc-mode far-field \
      --check-residual \
      --diag-vtk \
      --tfinal 3.15e7 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v22_test2_bdr_loading \
      --output-prefix bp5_v22_t2_bdrload \
      --write-every-step
