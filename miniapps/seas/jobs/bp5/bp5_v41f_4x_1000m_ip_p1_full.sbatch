#!/bin/bash
#SBATCH -J bp5_v41f_4x_ip       # Job name
#SBATCH -o bp5_v41f_4x_ip_%j.out
#SBATCH -e bp5_v41f_4x_ip_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 16                        # 16 nodes (147k tets, ~2.3x baseline)
#SBATCH -n 800                       # MPI tasks (50/node)
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v41f: Domain truncation study — 4x domain, IP p=1, 1000m mesh, full 1800 yr
#
# Tests boundary truncation effects by quadrupling domain size in all directions:
#   Standard: X=[-200,200], Y=[-100,100], Z=[-100,0] km
#   4x:       X=[-800,800], Y=[-400,400], Z=[-400,0] km
#
# The fault geometry (100x40 km) and fault-region mesh resolution (1000m) are unchanged.
# Far-field uses res=60 km elements (larger than 1x/2x to keep element count reasonable).
#
# Expected effects vs 1x domain:
#   1. tau_dip at center: should decrease significantly (boundary reflections ~1/16)
#   2. Strike slip rate: should increase slightly toward Tandem's value
#   3. Essentially eliminates boundary truncation as a confounding factor
#
# Compare with v41e (2x domain) to quantify convergence of boundary effects.
#
# PREREQUISITE: Generate the 4x mesh before submitting this job:
#   conda activate pythonenv
#   cd bp5/mesh/reference
#   gmsh -3 bp5_tandem_4x.geo -setnumber res_f 1 -o bp5_tandem_4x_1000m.msh

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

mkdir -p bp5/results_v41f_4x_1000m_ip_p1_full

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_4x_1000m.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v41f_4x_1000m_ip_p1_full \
      --output-prefix bp5_v41f_4x_ip_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
