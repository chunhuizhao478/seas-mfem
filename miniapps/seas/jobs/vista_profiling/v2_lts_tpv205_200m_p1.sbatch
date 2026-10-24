#!/bin/bash
#SBATCH -J v2_lts_tpv205_200m_p1
#SBATCH -o v2_lts_tpv205_200m_p1_%j.out
#SBATCH -e v2_lts_tpv205_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 00:30:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131,i617-101
# Fault-writer note rev1 check V2 (Chunhui Zhao, 2026-10-24): W0 (chunks that divide a frame) and W1
# (the frames' deflate on the other ranks, --fault-frame-compress spread, the new default) on V1 p1's
# setup.  TPV205 on the 200 m gold mesh, p = 1, --lts rate2, 0 to 12.0 s, 4 nodes x 144 ranks, the
# deck's fault deflate level 1, two variants:
#   spread  the defaults (W0 + W1)
#   rank0   --fault-frame-compress rank0 (W0 alone: rank 0 compresses, as in V1 o2 1095224)
# Both have station files meant to be bitwise equal to V1 o2's (1095224, kept on scratch; it differs
# from gold 1082586 by the 8.4e-11 round-off of the rank-0 mesh setup); their wall times against V1
# o2's 305.4 s give the gain of W0 and of W0 + W1.  Each run also writes
# tpv205_rank_balance.csv (L1).  Then the two fault files are compared byte for byte in their chunks,
# and the spread file's values against V1 o2's kept file (other chunk lengths, so values only).
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_fault_writer_vista.sbatch passed at this
# commit.  Estimate: V1 o2 took 305 s; two variants at most 10 min plus two comparisons of 15 GB of
# values, about 5 min = 1 node-h = 0.33 SU.  Requested: 4 nodes x 0.5 h x 0.33 SU = 0.7 SU at most.
TAG=v2
PROB=tpv205
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=1
TFINAL=12.0
# the station reference: V1 o2 (run_m0_vista.sh compares the variants marked gold with it)
GOLD=$SCRATCH/nitsche_profiling/v1_tpv205_200m_p1_1095224/o2
VARIANTS=(
   "spread|tpv205/configs/tpv205_nitsche.toml||gold"
   "rank0|tpv205/configs/tpv205_nitsche.toml|--fault-frame-compress rank0|gold"
)
FAULT_COMPARE=(
   "@BASE@/spread/fault.vtkhdf|@BASE@/rank0/fault.vtkhdf|"
   "@BASE@/spread/fault.vtkhdf|$SCRATCH/nitsche_profiling/v1_tpv205_200m_p1_1095224/o2/fault.vtkhdf|--values-only"
)
source "$SLURM_SUBMIT_DIR/run_m0_vista.sh"
