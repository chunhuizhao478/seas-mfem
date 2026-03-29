#!/bin/bash
#SBATCH -J bp5_v55_d8d4v2_p1         # Job name
#SBATCH -o bp5_v55_d8d4v2_p1_%j.out
#SBATCH -e bp5_v55_d8d4v2_p1_%j.err
#SBATCH -p development               # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks
#SBATCH -t 02:00:00                  # 2 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v55 D8+D4 v2: Corrected Tandem-matched friction solver
#
# D8 (corrected): V_vec negated to match Tandem, with GetSlip/SetSlip
#   negation so the domain solver sees unchanged physics. Previous D8
#   attempt (v55 d8d4) caused blowup because the negated slip flipped
#   the displacement field. Now GetSlip negates S back to physical slip
#   before the domain solver, preserving the same g^F as baseline.
#   Confirmed: D8-only run (7619500) = identical trajectory to baseline.
#
# D4: Brent solver in log10(V) space [-32, log10(tau/eta)], matching
#   Tandem DieterichRuinaBase.h:90-132. Provides uniform relative
#   precision across all slip rates. Fallback bracket included.
#   eta=0 direct inversion matches Tandem Finv formula exactly.
#
# All tests pass: 97+76+32+23 = 228 tests, 0 failures.

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

mkdir -p bp5/results_v55_1000m_ip_p1_d8d4v2

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver mumps-blr \
      --blr-tol 1e-12 \
      --check-residual \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v55_1000m_ip_p1_d8d4v2 \
      --output-prefix bp5_v55_d8d4v2_p1 \
      --checkpoint-interval 1000 \
      --write-every-step
