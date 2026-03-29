#!/bin/bash
#SBATCH -J bp5_v54_cg_p1            # Job name
#SBATCH -o bp5_v54_cg_p1_%j.out
#SBATCH -e bp5_v54_cg_p1_%j.err
#SBATCH -p normal                    # Queue
#SBATCH -N 8                         # 8 nodes
#SBATCH -n 400                       # 400 MPI ranks
#SBATCH -t 48:00:00                  # 48 hours
#SBATCH -A EAR20006                  # TACC allocation

# BP5 v54: IP p=1, 1000m mesh, CG+AMG solver — nucleation fix
#
# Root cause (v54 analysis):
#   v53 showed that MUMPS-BLR solver residual (~1e-7 at fault DOFs) gets
#   amplified ~300x by the IP penalty correction in ComputeTraction, producing
#   a growing traction bias (~0.38 MPa by t=38s) that kills nucleation.
#   MFEM V_strike peaks at t=3.1s then dies, while Tandem (using PETSc CG
#   with rtol=1e-12) recovers and earthquakes at t=82s.
#
# Fix:
#   Switch from MUMPS-BLR (direct, approximate factorization) to CG+AMG
#   (iterative, guarantees ||Ax-b||/||b|| < rtol=1e-10). The existing
#   CG_AMG solver uses HypreBoomerAMG with SetElasticityOptions().
#   Penalty-amplified CG residual: 300 * 1e-13 ≈ 3e-11 Pa (negligible).
#
# Other changes from v53:
#   - --zero-dip-traction: enforce analytically correct tau_dip=0 for BP5
#     (validated in v52, reduces DG cross-component contamination)
#   - Removed --write-every-step (was causing I/O overhead; use adaptive output)
#   - Uncommitted code changes: Tandem-style face-local interpolation for
#     station output (exact FE evaluation instead of nearest-DOF snap)
#
# Configuration:
#   - Mesh: bp5_tandem.msh (1000m reference mesh)
#   - Order: p=1
#   - DG: IP
#   - Solver: CG + BoomerAMG (rtol=1e-10, SetElasticityOptions)
#   - Initialization: Tandem-style (default)
#   - Runtime: 48h wall clock, 8 nodes / 400 MPI ranks
#
# Expected:
#   - Nucleation at strk-24dp+10 should NOT die (V should re-accelerate
#     after initial transient dip, matching Tandem's earthquake at t~82s)
#   - If earthquake occurs, simulation should proceed through post-seismic
#     and into interseismic with growing dt (dt >> 1e6 s)

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

mkdir -p bp5/results_v54_1000m_ip_p1_cg

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --solver cg \
      --check-residual \
      --zero-dip-traction \
      --tfinal 5.68e10 \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_v54_1000m_ip_p1_cg \
      --output-prefix bp5_v54_cg_p1 \
      --checkpoint-interval 1000
