#!/bin/bash
#SBATCH -J b0_tpv104_200m_p1
#SBATCH -o b0_tpv104_200m_p1_%j.out
#SBATCH -e b0_tpv104_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 512
#SBATCH -t 00:30:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131
# B0 of rheology plan rev 1 (Chunhui Zhao, 2026-10-20): tpv104 on the 200 m gold mesh, p = 1, 0 -> 1.0 s, the
# defaults, with the deck, mesh and rank layout of the global-stepping gold (4 nodes x 128 ranks); pass
# when this commit computes what the trunk's code computed in the 1 s run of LTS job 1082587
# (run_b0_vista.sh).  Submit from miniapps/seas/jobs/vista_rheology/ AFTER build_tpv2627_vista.sbatch passed
# at this commit.  Node i617-131 is excluded (out of memory in jobs 1082386-87).
# Estimate from that run: 143 s for 0 -> 1.0 s plus setup, under 10 min.
# Requested: 4 nodes x 0.5 h x 0.33 SU = 0.66 SU at most.
PROB=tpv104
DECK=tpv104/configs/tpv104_nitsche.toml
MESH=tpv104/mesh/tpv104_200m.msh
ORDER=1
GOLD=tpv104/gold/nitsche-kv0p1_p1_cd_pc_200m_jobid1075919_git81759da/results
REF=tpv104/gold/nitsche-kv0p1_p1_cd_pc_lts-rate2-buf3_200m_jobid1082587_gitac3adeb
REF_DIFF="1.000e-10 (tpv104_station_x2_0_x3_3.dat V)"
source "$SLURM_SUBMIT_DIR/run_b0_vista.sh"
