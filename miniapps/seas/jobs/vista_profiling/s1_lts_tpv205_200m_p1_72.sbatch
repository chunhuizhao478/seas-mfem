#!/bin/bash
#SBATCH -J s1_lts_tpv205_200m_p1_72
#SBATCH -o s1_lts_tpv205_200m_p1_72_%j.out
#SBATCH -e s1_lts_tpv205_200m_p1_72_%j.err
#SBATCH -p gg
#SBATCH -N 8
#SBATCH -n 576
#SBATCH --ntasks-per-node=72
#SBATCH -t 00:30:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131,i617-101
# Slow-parts study S1 (Chunhui Zhao, 2026-10-25): the memory-bandwidth test.  The run of
# s1_lts_tpv205_200m_p1.sbatch (the same 576 ranks and cut, 0 to 12.0 s) on 8 nodes x 72 ranks, each rank
# bound to every second core (bind_spread.sh: 36 ranks on each Grace CPU), so every rank has twice the
# memory bandwidth and cache of a full node:
#   spread72  stations meant to be bitwise equal to V3 tets' (1100111); the step time and the phases
#             against s1 timers' say how much of the time, and of the slowest ranks' extra, is memory
# Before the run the job prints where every rank sits (s1_hooks.sh).
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_l2_vista.sbatch passed at this commit.
# Estimate: at most V3's 172 s plus the setup, about 4 min = 0.53 node-h = 0.18 SU.
# Requested: 8 nodes x 0.5 h x 0.33 SU = 1.3 SU at most.
TAG=s1x72
PROB=tpv205
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=1
TFINAL=12.0
GOLD=$SCRATCH/nitsche_profiling/v3_tpv205_200m_p1_1100111/tets
VARIANTS=(
   "spread72|tpv205/configs/tpv205_nitsche.toml||gold|jobs/vista_profiling/bind_spread.sh"
)
source "$SLURM_SUBMIT_DIR/s1_hooks.sh"
source "$SLURM_SUBMIT_DIR/run_m0_vista.sh"
