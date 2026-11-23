#!/bin/bash
#SBATCH -J b0_tpv205_200m_p1
#SBATCH -o b0_tpv205_200m_p1_%j.out
#SBATCH -e b0_tpv205_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 00:30:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131
# B0 of rheology plan rev 1 (Chunhui Zhao, 2026-10-20): tpv205 on the 200 m gold mesh, p = 1, 0 -> 1.0 s, the
# defaults, with the deck, mesh and rank layout of the global-stepping gold (4 nodes x 144 ranks); pass
# when this commit computes what the trunk's code computed in the 1 s run of LTS job 1082586
# (run_b0_vista.sh).  Submit from miniapps/seas/jobs/vista_rheology/ AFTER build_tpv2627_vista.sbatch passed
# at this commit.  Node i617-131 is excluded (out of memory in jobs 1082386-87).
# Estimate from that run: 100 s for 0 -> 1.0 s plus setup, under 10 min.
# Requested: 4 nodes x 0.5 h x 0.33 SU = 0.66 SU at most.
PROB=tpv205
DECK=tpv205/configs/tpv205_nitsche.toml
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=1
GOLD=tpv205/gold/nitsche-kv0p1_p1_cd_pc_200m_jobid1076378_git81759da/results
REF=tpv205/gold/nitsche-kv0p1_p1_cd_pc_lts-rate2-buf3_200m_jobid1082586_gitac3adeb
REF_DIFF="2.965e-12 (tpv205_station_x2_0_x3_7.5.dat V)"
source "$SLURM_SUBMIT_DIR/run_b0_vista.sh"
