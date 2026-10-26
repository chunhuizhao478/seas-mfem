#!/bin/bash
#SBATCH -J s1_lts_tpv205_200m_p2
#SBATCH -o s1_lts_tpv205_200m_p2_%j.out
#SBATCH -e s1_lts_tpv205_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2304
#SBATCH -t 00:30:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131,i617-101
# Slow-parts study S1 (Chunhui Zhao, 2026-10-25) at p = 2: TPV205 on the 200 m gold mesh, --lts rate2,
# 0 to 12.0 s, 16 nodes x 144 ranks, the defaults (the tets cut of V3 1100112):
#   timers  each piece of the ticks per rank (tpv205_rank_balance.csv, ph_* columns; the "rank balance
#           phases" log line); stations meant to be bitwise equal to V3 tets'
# Before the run the job prints where every rank sits (s1_hooks.sh).
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_l2_vista.sbatch passed at this commit.
# Estimate: V3's tets took 414 s plus the setup, about 8 min = 2.1 node-h = 0.7 SU.
# Requested: 16 nodes x 0.5 h x 0.33 SU = 2.6 SU at most.
TAG=s1
PROB=tpv205
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=2
TFINAL=12.0
GOLD=$SCRATCH/nitsche_profiling/v3_tpv205_200m_p2_1100112/tets
VARIANTS=(
   "timers|tpv205/configs/tpv205_nitsche.toml||gold"
)
source "$SLURM_SUBMIT_DIR/s1_hooks.sh"
source "$SLURM_SUBMIT_DIR/run_m0_vista.sh"
