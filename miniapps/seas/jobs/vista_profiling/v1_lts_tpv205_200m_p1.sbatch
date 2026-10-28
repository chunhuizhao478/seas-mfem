#!/bin/bash
#SBATCH -J v1_lts_tpv205_200m_p1
#SBATCH -o v1_lts_tpv205_200m_p1_%j.out
#SBATCH -e v1_lts_tpv205_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 00:30:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131
# Speed-next verification V1 (Chunhui Zhao, 2026-10-20): O2 (stations) and O1 (fault frames) against
# today's fastest setup.  TPV205 on the 200 m gold mesh, p = 1, --lts rate2, 0 to 12.0 s, the rank layout
# of LTS gold 1082586 (4 nodes x 144 ranks), two variants:
#   o2    the defaults (--station-comm sub, --station-flush-every 100; the deck's fault deflate level 1)
#   o1o2  the same with --paraview-fault-deflate-level none (no filter on fault.vtkhdf)
# Both have station files meant to be bitwise equal to 1082586's over the 12 s; their wall times against
# 1082586's 359 s give the speedup of O2 and of O1 + O2.
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_speed_next_vista.sbatch passed at this commit.
# Estimate: 1082586 took 359 s; two variants at most 12 min = 0.8 node-h = 0.3 SU.  Requested: 4 nodes x
# 0.5 h x 0.33 SU = 0.7 SU at most.
TAG=v1
PROB=tpv205
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=1
TFINAL=12.0
GOLD=tpv205/gold/nitsche-kv0p1_p1_cd_pc_lts-rate2-buf3_200m_jobid1082586_gitac3adeb/results
VARIANTS=(
   "o2|tpv205/configs/tpv205_nitsche.toml||gold"
   "o1o2|tpv205/configs/tpv205_nitsche.toml|--paraview-fault-deflate-level none|gold"
)
source "$SLURM_SUBMIT_DIR/run_m0_vista.sh"
