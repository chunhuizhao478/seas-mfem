#!/bin/bash
#SBATCH -J b0_tpv205_200m_p2
#SBATCH -o b0_tpv205_200m_p2_%j.out
#SBATCH -e b0_tpv205_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2304
#SBATCH -t 00:30:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131
# B0 of rheology plan rev 1 (Chunhui Zhao, 2026-10-20): tpv205 on the 200 m gold mesh, p = 2, 0 -> 1.0 s, the
# defaults, with the deck, mesh and rank layout of the global-stepping gold (16 nodes x 144 ranks); pass
# when this commit computes what the trunk's code computed in the 1 s run of LTS job 1082384
# (run_b0_vista.sh).  Submit from miniapps/seas/jobs/vista_rheology/ AFTER build_tpv2627_vista.sbatch passed
# at this commit.  Node i617-131 is excluded (out of memory in jobs 1082386-87).
# Estimate from that run: 199 s for 0 -> 1.0 s plus setup, under 10 min.
# Requested: 16 nodes x 0.5 h x 0.33 SU = 2.64 SU at most.
PROB=tpv205
DECK=tpv205/configs/tpv205_nitsche.toml
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=2
GOLD=tpv205/gold/nitsche-kv0p1_p2_cd_pc_200m_jobid1075916_git81759da/results
REF=tpv205/gold/nitsche-kv0p1_p2_cd_pc_lts-rate2-buf3_200m_jobid1082384_gitac3adeb
REF_DIFF="2.691e-11 (tpv205_station_x2_0_x3_7.5.dat V)"
source "$SLURM_SUBMIT_DIR/run_b0_vista.sh"
