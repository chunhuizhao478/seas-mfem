#!/bin/bash
#SBATCH -J sst_probe
#SBATCH -o sst_probe_%j.out
#SBATCH -e sst_probe_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH --ntasks-per-node=144
#SBATCH -t 01:00:00
#SBATCH -A EAR24034
#SBATCH --exclude=i617-131,i617-101
# SeisSol layout probe (SeisSol and MFEM speed test on Vista, plan rev1 section 6, Chunhui Zhao
# 2026-10-30): TPV205 o3 on the 200 m gold mesh, 4 nodes, 0 to 2.0 s, SeisSol only, one run per entry
# of SS_LAYOUTS (MPI ranks per node, each rank with 144 / ranks - 1 OpenMP threads and one
# communication thread).  The fastest "Simulation time (total)" picks SS_LAYOUTS of the run jobs (a
# follow-up commit).
# First probe 1112916 (2, 4, 8 per node): 22.7 s at 2, 15.4 s at 8; 4 per node crashed in its first
# steps (segfault in UCX, "Local length error on mlx5_0:1/IB", under MPI_Test on SeisSol's
# communication thread).  This repeat (4, 4, 8, 16; Chunhui Zhao 2026-11-01) shows whether the crash
# recurs at 4 and whether 16 per node beats 8; the job prints the MPI and UCX libraries the binary loads.
# Submit from miniapps/seas/jobs/seissol_speed_vista/ AFTER build_speed_test_vista.sbatch passed at this
# commit (1112915 at 1831dfa5 built the binaries and meshes; the build job need not run again).
# Estimate: 1112916 took 2:10 for three runs; four runs about 4 min = 0.3 node-h, 0.1 SU (inferred).
# Requested: 4 nodes x 1 h x 0.33 SU = 1.3 SU at most.
TAG=probe
PROB=tpv205
MESH=tpv205/mesh/tpv2053d_200m.msh
PUML=tpv2053d_200m
ORDER=2
TFINAL=2.0
MFEM_RUN=
SS_LAYOUTS="4 4 8 16"
source "$SLURM_SUBMIT_DIR/run_pair_vista.sh"
