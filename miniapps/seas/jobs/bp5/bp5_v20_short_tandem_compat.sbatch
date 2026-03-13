#!/bin/bash
#SBATCH -J bp5_v20_short_tc     # Job name
#SBATCH -o bp5_v20_short_tc_%j.out
#SBATCH -e bp5_v20_short_tc_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                    # nodes
#SBATCH -n 400                  # MPI tasks
#SBATCH -t 02:00:00             # 2 hours (short run)
#SBATCH -A EAR20006             # TACC allocation

# BP5 v20: Tandem-compatible parameters with psi_init fix
# Goal: Match Tandem nucleation to isolate parameter effects
# Uses V_nuc=0.01 (Tandem default) and delta_tau_factor=0 (no overstress)
# Expected: faster nucleation since psi_init fix + no delta_tau = clean comparison

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

mkdir -p bp5/results_v20_short_tandem_compat

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --solver mumps-blr \
      --V-nuc 0.01 \
      --delta-tau-factor 0 \
      --check-residual \
      --monitor-traction 100 \
      --tfinal 9.45e9 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v20_short_tandem_compat \
      --output-prefix bp5_v20_short_tc \
      --checkpoint-interval 1000 \
      --write-every-step
