#!/bin/bash
#SBATCH -J m0_lts_tpv205_200m_p1
#SBATCH -o m0_lts_tpv205_200m_p1_%j.out
#SBATCH -e m0_lts_tpv205_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 00:30:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131
# M0 (speed-next plan rev 1.1, §6; Chunhui Zhao, 2026-10-20): the whole-run frame cost in today's fastest
# setup.  TPV205 on the 200 m gold mesh, p = 1, --lts rate2, 0 to 12.0 s (600 frames), the rank layout of
# LTS gold 1082586 (4 nodes x 144 ranks), with the deck copy that sets paraview_fault_deflate_level = 0.
# Its station files are meant to be bitwise equal to 1082586's (the deflate level touches only the fault
# file); its wall time against 1082586's 359 s gives O1's gain at deflate level 0.
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_lts_vista.sbatch passed at this commit.
# Estimate: 1082586 took 359 s, of which 130 s were frames; with half the frame time about 5 min =
# 0.35 node-h = 0.1 SU.  Requested: 4 nodes x 0.5 h x 0.33 SU = 0.7 SU at most.
PROB=tpv205
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=1
TFINAL=12.0
GOLD=tpv205/gold/nitsche-kv0p1_p1_cd_pc_lts-rate2-buf3_200m_jobid1082586_gitac3adeb/results
VARIANTS=(
   "deflate0|tpv205/configs/tpv205_nitsche_m0_deflate0.toml||gold"
)
source "$SLURM_SUBMIT_DIR/run_m0_vista.sh"
