#!/bin/bash
#SBATCH -J map_lts_tpv205_200m_p1
#SBATCH -o map_lts_tpv205_200m_p1_%j.out
#SBATCH -e map_lts_tpv205_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 00:40:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131
# Profile (speed-next, Chunhui Zhao 2026-10-20): the breakdown of each part of the code in today's fastest
# setup at p = 1.  TPV205 on the 200 m gold mesh, the gold deck, --lts rate2, 0 to 3.0 s (150 frames),
# the rank layout of LTS gold 1082586 (4 nodes x 144 ranks), under Linaro MAP (run_map_vista.sh).
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_map_vista.sbatch passed at this commit.
# Estimate: 70 s setup + a quarter of 1082586's 289 s of stepping and output, about 2.5 min, plus MAP's
# start and its merge of the samples, about 10 min = 0.7 node-h = 0.2 SU.  Requested: 4 nodes x 0.67 h
# x 0.33 SU = 0.9 SU at most.
PROB=tpv205
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=1
TFINAL=3.0
DECK=tpv205/configs/tpv205_nitsche.toml
GOLD=tpv205/gold/nitsche-kv0p1_p1_cd_pc_lts-rate2-buf3_200m_jobid1082586_gitac3adeb/results
source "$SLURM_SUBMIT_DIR/run_map_vista.sh"
