#!/bin/bash
#SBATCH -J m0_lts_tpv205_200m_p2
#SBATCH -o m0_lts_tpv205_200m_p2_%j.out
#SBATCH -e m0_lts_tpv205_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2304
#SBATCH -t 00:45:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131
# M0 (speed-next plan rev 1.1, §6; Chunhui Zhao, 2026-10-20): where the time outside the step goes in
# today's fastest setup.  TPV205 on the 200 m gold mesh, p = 2, --lts rate2, 0 to 1.0 s, the rank layout
# of LTS gold 1082384 (16 nodes x 144 ranks), four variants with existing flags only:
#   base     the defaults (station files bitwise equal to the first 1.0 s of 1082384)
#   every10  --station-every 10: the share of all per-tick station work
#   first    --station-copies first: the share of the copy averaging and its collectives
#   deflate0 deck copy with paraview_fault_deflate_level = 0: the frame cost (stations bitwise as base)
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_lts_vista.sbatch passed at this commit.
# Estimate: each variant about 87 s setup + 7394 ticks x 7.3 ms + 49 frames, about 2.5 min; four about
# 10-12 min = 3.2 node-h = 1.1 SU.  Requested: 16 nodes x 0.75 h x 0.33 SU = 4.0 SU at most.
PROB=tpv205
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=2
TFINAL=1.0
GOLD=tpv205/gold/nitsche-kv0p1_p2_cd_pc_lts-rate2-buf3_200m_jobid1082384_gitac3adeb/results
VARIANTS=(
   "base|tpv205/configs/tpv205_nitsche.toml||gold"
   "every10|tpv205/configs/tpv205_nitsche.toml|--station-every 10|"
   "first|tpv205/configs/tpv205_nitsche.toml|--station-copies first|"
   "deflate0|tpv205/configs/tpv205_nitsche_m0_deflate0.toml||gold"
)
source "$SLURM_SUBMIT_DIR/run_m0_vista.sh"
