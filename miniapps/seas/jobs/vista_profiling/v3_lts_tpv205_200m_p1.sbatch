#!/bin/bash
#SBATCH -J v3_lts_tpv205_200m_p1
#SBATCH -o v3_lts_tpv205_200m_p1_%j.out
#SBATCH -e v3_lts_tpv205_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 00:40:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131,i617-101
# L2 note rev1 check V3 (Chunhui Zhao, 2026-10-25): the block-work LTS partition
# (--lts-partition-weight blocks) on V2 p1's setup.  TPV205 on the 200 m gold mesh, p = 1, --lts rate2,
# 0 to 12.0 s, 4 nodes x 144 ranks, the defaults otherwise (fault deflate level 1, spread compression),
# two variants in one allocation, since node differences between jobs are as large as the effect:
#   tets    the defaults (--lts-partition-weight tets): stations meant to be bitwise equal to V2's
#           (1098850 spread, kept on scratch), wall time against its 169.3 s
#   blocks  --lts-partition-weight blocks: another cut, so stations meant to agree with V2's at
#           round-off (the 8e-11 of the rank-0 mesh setup is the yardstick); the gain is its wall and
#           ode.Step time against tets in this job
# Both print the block work per cluster and the tick model of their cut on the partition line, and
# write tpv205_rank_balance.csv (L1).
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_l2_vista.sbatch passed at this commit.
# Estimate: V2's spread run took 169 s plus about 60 s of setup; two variants about 8 min = 0.5 node-h
# = 0.18 SU.  Requested: 4 nodes x 0.67 h x 0.33 SU = 0.9 SU at most.
TAG=v3
PROB=tpv205
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=1
TFINAL=12.0
# the station reference: V2 spread (run_m0_vista.sh compares the variants marked gold or roundoff)
GOLD=$SCRATCH/nitsche_profiling/v2_tpv205_200m_p1_1098850/spread
VARIANTS=(
   "tets|tpv205/configs/tpv205_nitsche.toml||gold"
   "blocks|tpv205/configs/tpv205_nitsche.toml|--lts-partition-weight blocks|roundoff"
)
source "$SLURM_SUBMIT_DIR/run_m0_vista.sh"
