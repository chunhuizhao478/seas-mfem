#!/bin/bash
#SBATCH -J bp5_v51_diag_dip        # Job name
#SBATCH -o bp5_v51_diag_dip_%j.out
#SBATCH -e bp5_v51_diag_dip_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 1                         # 1 node — serial for clean output
#SBATCH -n 1                         # 1 rank
#SBATCH -t 00:30:00                  # 30 minutes (only need first few steps)
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v51: Dip-traction component diagnostic
#
# Dumps per-component traction (global xyz + local dip/strike) at each fault DOF
# during the FIRST non-zero-slip ComputeTraction call.
#
# Purpose: identify whether the spurious dip traction comes from:
# (1) T_stress (stress average {sigma.n}) — mesh geometry / DG discretization
# (2) T_correction (penalty correction) — operator-split inconsistency
#
# Uses 1000m mesh, p=2 (same as v50 production, but serial for clean output).
# Short run: only need first earthquake step, not the full cycle.

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module list

export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# Build with v51 diagnostic
make clean-build
make seas_bp5_full || { echo "BUILD FAILED"; exit 1; }

mkdir -p bp5/results_v51_diag_dip

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 2 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 100 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v51_diag_dip \
      --output-prefix bp5_v51_diag \
      --write-every-step \
      --diag-dip-traction
