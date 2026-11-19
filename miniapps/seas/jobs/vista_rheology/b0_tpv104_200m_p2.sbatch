#!/bin/bash
#SBATCH -J b0_tpv104_200m_p2
#SBATCH -o b0_tpv104_200m_p2_%j.out
#SBATCH -e b0_tpv104_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2048
#SBATCH -t 00:30:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131
# B0 of rheology plan rev 1 (Chunhui Zhao, 2026-10-20): tpv104 on the 200 m gold mesh, p = 2, 0 -> 1.0 s, the
# defaults, with the deck, mesh and rank layout of the global-stepping gold (16 nodes x 128 ranks); pass
# when this commit computes what the trunk's code computed in the 1 s run of LTS job 1082385
# (run_b0_vista.sh).  Submit from miniapps/seas/jobs/vista_rheology/ AFTER build_tpv2627_vista.sbatch passed
# at this commit.  Node i617-131 is excluded (out of memory in jobs 1082386-87).
# Estimate from that run: 252 s for 0 -> 1.0 s plus setup, under 10 min.
# Requested: 16 nodes x 0.5 h x 0.33 SU = 2.64 SU at most.
PROB=tpv104
DECK=tpv104/configs/tpv104_nitsche.toml
MESH=tpv104/mesh/tpv104_200m.msh
ORDER=2
GOLD=tpv104/gold/nitsche-kv0p1_p2_cd_pc_200m_jobid1075917_git81759da/results
REF=tpv104/gold/nitsche-kv0p1_p2_cd_pc_lts-rate2-buf3_200m_jobid1082385_gitac3adeb
REF_DIFF="1.225e-11 (tpv104_station_x2_0_x3_7.5.dat tau)"
source "$SLURM_SUBMIT_DIR/run_b0_vista.sh"
