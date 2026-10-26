#!/bin/bash
#SBATCH -J s1_lts_tpv205_200m_p1
#SBATCH -o s1_lts_tpv205_200m_p1_%j.out
#SBATCH -e s1_lts_tpv205_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 00:40:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131,i617-101
# Slow-parts study S1 (Chunhui Zhao, 2026-10-25): where the slowest rank's time goes.  TPV205 on the
# 200 m gold mesh, p = 1, --lts rate2, 0 to 12.0 s, 4 nodes x 144 ranks, the defaults (the tets cut of V3
# 1100111), two variants in one allocation:
#   timers  the defaults: each piece of the ticks per rank in tpv205_rank_balance.csv (ph_* columns, the
#           NitscheLtsRun phase timers) and the "rank balance phases" log line; stations meant to be
#           bitwise equal to V3 tets'
#   perf    the same run with four ranks under perf record (perf_wrap.sh, 499 Hz): the slowest by compute
#           time in timers, the slowest on another node, and a typical rank on each of those two nodes
#           (pick_perf_ranks.py); cycles and instructions (and backend memory stalls when the node has
#           the event), so each symbol's share and instructions per cycle compare between slow and typical
#           ranks; stations meant to be bitwise equal too
# Before timers the job prints where every rank sits (s1_hooks.sh, rank_bindings.sh).
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_l2_vista.sbatch passed at this commit.
# Estimate: V3's tets took 172 s plus the setup; two variants about 8 min = 0.53 node-h = 0.18 SU.
# Requested: 4 nodes x 0.67 h x 0.33 SU = 0.9 SU at most.
TAG=s1
PROB=tpv205
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=1
TFINAL=12.0
GOLD=$SCRATCH/nitsche_profiling/v3_tpv205_200m_p1_1100111/tets
VARIANTS=(
   "timers|tpv205/configs/tpv205_nitsche.toml||gold"
   "perf|tpv205/configs/tpv205_nitsche.toml||gold|jobs/vista_profiling/perf_wrap.sh"
)
source "$SLURM_SUBMIT_DIR/s1_hooks.sh"
source "$SLURM_SUBMIT_DIR/run_m0_vista.sh"
