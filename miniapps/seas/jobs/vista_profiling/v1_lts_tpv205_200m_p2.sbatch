#!/bin/bash
#SBATCH -J v1_lts_tpv205_200m_p2
#SBATCH -o v1_lts_tpv205_200m_p2_%j.out
#SBATCH -e v1_lts_tpv205_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2304
#SBATCH -t 00:50:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131
# Speed-next verification V1 (Chunhui Zhao, 2026-10-20): as v1_lts_tpv205_200m_p1.sbatch at p = 2, the
# rank layout of LTS gold 1082384 (16 nodes x 144 ranks), 0 to 12.0 s, variants o2 (the defaults) and
# o1o2 (--paraview-fault-deflate-level none), stations compared with 1082384 over the 12 s (meant to be
# bitwise equal); wall times against 1082384's 739 s.
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_speed_next_vista.sbatch passed at this commit.
# Estimate: two variants at most 2 x 12.5 min = 6.7 node-h = 2.2 SU.  Requested: 16 nodes x 0.83 h x
# 0.33 SU = 4.4 SU at most.
TAG=v1
PROB=tpv205
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=2
TFINAL=12.0
GOLD=tpv205/gold/nitsche-kv0p1_p2_cd_pc_lts-rate2-buf3_200m_jobid1082384_gitac3adeb/results
VARIANTS=(
   "o2|tpv205/configs/tpv205_nitsche.toml||gold"
   "o1o2|tpv205/configs/tpv205_nitsche.toml|--paraview-fault-deflate-level none|gold"
)
source "$SLURM_SUBMIT_DIR/run_m0_vista.sh"
