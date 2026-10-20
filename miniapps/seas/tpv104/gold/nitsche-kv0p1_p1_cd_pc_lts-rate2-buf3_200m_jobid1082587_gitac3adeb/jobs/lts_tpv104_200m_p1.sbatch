#!/bin/bash
#SBATCH -J lts_tpv104_200m_p1
#SBATCH -o lts_tpv104_200m_p1_%j.out
#SBATCH -e lts_tpv104_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 512
#SBATCH -t 01:00:00
#SBATCH -A EAR26024
# Local time stepping (plan item 7, P5; Chunhui Zhao, 2026-10-19): tpv104 on the 200 m gold mesh, p = 1,
# --lts rate2 from 0 to 12.0 s (the new gold), after a 1.0 s global-stepping run at this commit for the
# speed reference; deck tpv104/configs/tpv104_nitsche.toml and the rank layout of gold 1075919 (4 nodes x 128 ranks).
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_lts_vista.sbatch passed at this commit.
# Estimate: as the TPV205 p1 job with 1.2x the steps and a rate-and-state fault (inferred): about 12 min =
# 0.8 node-h = 0.3 SU.
# Requested: 4 nodes x 1 h x 0.33 SU = 1.3 SU at most.
# 128 ranks per node, as the gold (every rank loads the serial mesh, about 1.5 GB per rank here).
PROB=tpv104
DECK=tpv104/configs/tpv104_nitsche.toml
MESH=tpv104/mesh/tpv104_200m.msh
ORDER=1
TFINAL=12.0
GOLD=tpv104/gold/nitsche-kv0p1_p1_cd_pc_200m_jobid1075919_git81759da/results
source "$SLURM_SUBMIT_DIR/run_lts_vista.sh"
