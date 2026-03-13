#!/bin/bash
#SBATCH -J bp5_v20_med_mblr     # Job name
#SBATCH -o bp5_v20_med_mblr_%j.out
#SBATCH -e bp5_v20_med_mblr_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours (medium run)
#SBATCH -A EAR20006             # TACC allocation

# BP5 v20: MUMPS BLR medium run with psi_init fix
# Goal: Capture second seismic event with fixed initial state (SCEC Eq.18)
# Fix: use SCEC Eq.18 psi(0) = f0 + b*ln(V0/V_init) instead of InitialStatePsi()

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

mkdir -p bp5/results_v20_med_mumps_blr

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --solver mumps-blr \
      --check-residual \
      --monitor-traction 100 \
      --tfinal 9.45e9 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v20_med_mumps_blr \
      --output-prefix bp5_v20_med_mblr \
      --checkpoint-interval 1000 \
      --write-every-step
