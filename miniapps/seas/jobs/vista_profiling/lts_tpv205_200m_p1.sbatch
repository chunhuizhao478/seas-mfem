#!/bin/bash
#SBATCH -J lts_tpv205_200m_p1
#SBATCH -o lts_tpv205_200m_p1_%j.out
#SBATCH -e lts_tpv205_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 01:00:00
#SBATCH -A EAR26024
# Local time stepping (plan item 7, P5; Chunhui Zhao, 2026-10-19): tpv205 on the 200 m gold mesh, p = 1,
# --lts rate2 from 0 to 12.0 s (the new gold), after a 1.0 s global-stepping run at this commit for the
# speed reference; deck tpv205/configs/tpv205_nitsche.toml and the rank layout of gold 1076378 (4 nodes x 144 ranks).
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_lts_vista.sbatch passed at this commit.
# Estimate: gold 1076378 stepped 41,614 steps in 742 s (old code); the defaults about 2x faster and LTS 1.5-2x on
# top (inferred): lts 3-5 min plus setup and output; gts 1.0 s about 2 min.  About 10 min = 0.7 node-h = 0.2 SU.
# Requested: 4 nodes x 1 h x 0.33 SU = 1.3 SU at most.
# 144 ranks per node, as the gold.
PROB=tpv205
DECK=tpv205/configs/tpv205_nitsche.toml
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=1
TFINAL=12.0
GOLD=tpv205/gold/nitsche-kv0p1_p1_cd_pc_200m_jobid1076378_git81759da/results
source "$SLURM_SUBMIT_DIR/run_lts_vista.sh"
