#!/bin/bash
#SBATCH -J map_lts_tpv205_200m_p2
#SBATCH -o map_lts_tpv205_200m_p2_%j.out
#SBATCH -e map_lts_tpv205_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2304
#SBATCH -t 01:00:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131
# Profile (speed-next, Chunhui Zhao 2026-10-20): the breakdown of each part of the code in today's fastest
# setup at p = 2.  TPV205 on the 200 m gold mesh, the gold deck, --lts rate2, 0 to 2.0 s (100 frames),
# the rank layout of LTS gold 1082384 (16 nodes x 144 ranks), under Linaro MAP (run_map_vista.sh).
# If MAP's licence on Vista does not cover 2304 processes, MAP exits and the job falls back to perf on
# four ranks; the p = 1 profile (576 ranks) still gives the full breakdown.
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_map_vista.sbatch passed at this commit.
# Estimate: 87 s setup + a sixth of 1082384's 652 s of stepping and output, about 3.3 min, plus MAP's
# start and its merge, about 12 min = 3.2 node-h = 1.1 SU, plus about a minute per failed MAP start.
# Each MAP try is stopped after 25 min (MAP_TRY_S), so two tries and the perf fallback fit.  Requested:
# 16 nodes x 1.0 h x 0.33 SU = 5.3 SU at most.
PROB=tpv205
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=2
TFINAL=2.0
MAP_TRY_S=1500
DECK=tpv205/configs/tpv205_nitsche.toml
GOLD=tpv205/gold/nitsche-kv0p1_p2_cd_pc_lts-rate2-buf3_200m_jobid1082384_gitac3adeb/results
source "$SLURM_SUBMIT_DIR/run_map_vista.sh"
