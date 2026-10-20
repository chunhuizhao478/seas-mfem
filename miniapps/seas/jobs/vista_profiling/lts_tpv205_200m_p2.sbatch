#!/bin/bash
#SBATCH -J lts_tpv205_200m_p2
#SBATCH -o lts_tpv205_200m_p2_%j.out
#SBATCH -e lts_tpv205_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2304
#SBATCH -t 01:00:00
#SBATCH -A EAR26024
# Local time stepping (plan item 7, P5; Chunhui Zhao, 2026-10-19): tpv205 on the 200 m gold mesh, p = 2,
# --lts rate2 from 0 to 12.0 s (the new gold), after a 1.0 s global-stepping run at this commit for the
# speed reference; deck tpv205/configs/tpv205_nitsche.toml and the rank layout of gold 1075916 (16 nodes x 144 ranks).
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_lts_vista.sbatch passed at this commit.
# Estimate: gts as job 1079640's dense run (197 s); lts 88,725 steps at 4.7-6.2 ms (9.36 ms of 1079640 over the
# 1.5-2x the cloud measured at one rank, inferred) plus 140-250 s setup (rank 0 also builds the clusters and the
# multi-constraint cut; the P0 report took 6 min for two cuts on one cloud core) and the fault output: 10-19 min.
# About 22 min in all = 5.9 node-h = 1.9 SU.
# Requested: 16 nodes x 1 h x 0.33 SU = 5.3 SU at most.
# 144 ranks per node, as the gold.
PROB=tpv205
DECK=tpv205/configs/tpv205_nitsche.toml
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=2
TFINAL=12.0
GOLD=tpv205/gold/nitsche-kv0p1_p2_cd_pc_200m_jobid1075916_git81759da/results
source "$SLURM_SUBMIT_DIR/run_lts_vista.sh"
