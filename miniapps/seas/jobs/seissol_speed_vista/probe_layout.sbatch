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
# 2026-10-30): TPV205 o3 on the 200 m gold mesh, 4 nodes, 0 to 2.0 s, SeisSol only, three layouts in one
# allocation: 2, 4 and 8 MPI ranks per node, each rank with 144 / ranks - 1 OpenMP threads and one
# communication thread.  The fastest "Simulation time (total)" picks SS_LAYOUTS of the run jobs (a
# follow-up commit).  Also the first run of the decks on this SeisSol commit: a wrong parameter name
# stops here, at about 1 SU, not in the paired runs.
# Submit from miniapps/seas/jobs/seissol_speed_vista/ AFTER build_speed_test_vista.sbatch passed at this
# commit.  Estimate: 3 runs of about 3 to 6 min on 4 nodes = 0.6 to 1.2 node-h, about 0.2 to 0.4 SU (inferred).
# Requested: 4 nodes x 1 h x 0.33 SU = 1.3 SU at most.
TAG=probe
PROB=tpv205
MESH=tpv205/mesh/tpv2053d_200m.msh
PUML=tpv2053d_200m
ORDER=2
TFINAL=2.0
MFEM_RUN=
SS_LAYOUTS="2 4 8"
source "$SLURM_SUBMIT_DIR/run_pair_vista.sh"
