#!/bin/bash
#SBATCH -J lts_tpv104_200m_p2
#SBATCH -o lts_tpv104_200m_p2_%j.out
#SBATCH -e lts_tpv104_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2048
#SBATCH -t 01:15:00
#SBATCH -A EAR26024
# Local time stepping (plan item 7, P5; Chunhui Zhao, 2026-10-19): tpv104 on the 200 m gold mesh, p = 2,
# --lts rate2 from 0 to 12.0 s (the new gold), after a 1.0 s global-stepping run at this commit for the
# speed reference; deck tpv104/configs/tpv104_nitsche.toml and the rank layout of gold 1075917 (16 nodes x 128 ranks).
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_lts_vista.sbatch passed at this commit.
# Estimate: gold 1075917 stepped at about 30 ms (old code); the defaults about 2.9x faster (as on TPV205, inferred)
# and LTS 1.5-2x on top: lts 104,388 steps at 5-7 ms plus 160-300 s setup (clusters and the multi-constraint
# cut on rank 0) and the fault output: 12-19 min; gts 1.0 s about 4 min.  About 24 min in all = 6.4 node-h = 2.1 SU.
# Requested: 16 nodes x 1.25 h x 0.33 SU = 6.6 SU at most.
# 128 ranks per node, as the gold (every rank loads the serial mesh, about 1.5 GB per rank here).
PROB=tpv104
DECK=tpv104/configs/tpv104_nitsche.toml
MESH=tpv104/mesh/tpv104_200m.msh
ORDER=2
TFINAL=12.0
GOLD=tpv104/gold/nitsche-kv0p1_p2_cd_pc_200m_jobid1075917_git81759da/results
source "$SLURM_SUBMIT_DIR/run_lts_vista.sh"
