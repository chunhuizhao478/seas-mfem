#!/bin/bash
#SBATCH -J bp5_v41e_2x_ip       # Job name
#SBATCH -o bp5_v41e_2x_ip_%j.out
#SBATCH -e bp5_v41e_2x_ip_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 12                        # 12 nodes (104k tets, ~1.6x baseline)
#SBATCH -n 600                       # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v41e: Domain truncation study — 2x domain, IP p=1, 1000m mesh, full 1800 yr
#
# Tests boundary truncation effects by doubling domain size in all directions:
#   Standard: X=[-200,200], Y=[-100,100], Z=[-100,0] km
#   2x:       X=[-400,400], Y=[-200,200], Z=[-200,0] km
#
# The fault geometry (100x40 km) and fault-region mesh resolution (1000m) are unchanged.
# The far-field uses res=40 km elements with automatic Gmsh grading.
#
# Expected effects:
#   1. tau_dip at center stations: should decrease (less boundary reflection)
#   2. Strike slip rate: should increase slightly (less boundary stiffening)
#   3. Recurrence: should be slightly shorter (more compliant far-field)
#
# PREREQUISITE: Generate the 2x mesh before submitting this job:
#   conda activate pythonenv
#   cd bp5/mesh/reference
#   gmsh -3 bp5_tandem_2x.geo -setnumber res_f 1 -o bp5_tandem_2x_1000m.msh
#
# The Lx, Ly, Lz parameters in the driver are auto-detected from the mesh extents,
# NOT hardcoded. The boundary loading u_D = sgn(Y)*Vp*t/2 works at any domain size
# because it's based on the face centroid Y-coordinate sign.

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

mkdir -p bp5/results_v41e_2x_1000m_ip_p1_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_2x_1000m.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v41e_2x_1000m_ip_p1_full \
      --output-prefix bp5_v41e_2x_ip_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
