#!/bin/bash
#SBATCH -J v3_lts_tpv205_200m_p2
#SBATCH -o v3_lts_tpv205_200m_p2_%j.out
#SBATCH -e v3_lts_tpv205_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2304
#SBATCH -t 00:40:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131,i617-101
# L2 note rev1 check V3 (Chunhui Zhao, 2026-10-25): as v3_lts_tpv205_200m_p1.sbatch at p = 2, the rank
# layout of V2 p2 1098851 (16 nodes x 144 ranks), 0 to 12.0 s, the same two variants:
#   tets    the defaults: stations meant to be bitwise equal to V2's (1098851 spread, kept on
#           scratch), wall time against its 432.8 s
#   blocks  --lts-partition-weight blocks: stations meant to agree with V2's at round-off; the gain is
#           its wall and ode.Step time against tets in this job
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_l2_vista.sbatch passed at this commit.
# Estimate: V2's spread run took 433 s plus about 2 min of setup; two variants about 17 min = 4.5
# node-h = 1.5 SU.  Requested: 16 nodes x 0.67 h x 0.33 SU = 3.5 SU at most.
TAG=v3
PROB=tpv205
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=2
TFINAL=12.0
# the station reference: V2 spread (run_m0_vista.sh compares the variants marked gold or roundoff)
GOLD=$SCRATCH/nitsche_profiling/v2_tpv205_200m_p2_1098851/spread
VARIANTS=(
   "tets|tpv205/configs/tpv205_nitsche.toml||gold"
   "blocks|tpv205/configs/tpv205_nitsche.toml|--lts-partition-weight blocks|roundoff"
)
source "$SLURM_SUBMIT_DIR/run_m0_vista.sh"
