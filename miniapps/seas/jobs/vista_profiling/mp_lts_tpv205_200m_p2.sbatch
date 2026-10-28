#!/bin/bash
#SBATCH -J mp_lts_tpv205_200m_p2
#SBATCH -o mp_lts_tpv205_200m_p2_%j.out
#SBATCH -e mp_lts_tpv205_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2304
#SBATCH -t 00:40:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131,i617-101
# Mesh setup check at 16 nodes x 144 ranks (Chunhui Zhao, 2026-10-22), the layout at which TPV31 p = 2
# and p = 3 (1087801, 1087802) stopped in setup at the node commit limit: TPV205 on the 200 m gold mesh,
# p = 2, --lts rate2, the rank layout of LTS gold 1082384.  Four variants through
# run_mesh_parts_vista.sh, 0 to 1.0 s each:
#   every  --mesh-setup every-rank (the path before this change): bitwise equal to the first 1.0 s of
#          1082384; its setup time and memory are the reference
#   rank0  the new default: against 1082384 (round-off differences)
#   split  one rank: --mesh-parts-out <dir> --mesh-parts-n 2304 writes the 2304 parts
#   parts  --mesh-parts <dir>: bitwise equal to rank0
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_mesh_parts_vista.sbatch passed at this commit.
# Estimate: M0's base variant (1085109) took 160 s for the same 1.0 s; three runs and the split (about
# 3 min) about 12 min = 3.2 node-h = 1.1 SU.  Requested: 16 nodes x 0.67 h x 0.33 SU = 3.5 SU at most.
PROB=tpv205
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=2
DECK=tpv205/configs/tpv205_nitsche.toml
LTS_GOLD=tpv205/gold/nitsche-kv0p1_p2_cd_pc_lts-rate2-buf3_200m_jobid1082384_gitac3adeb/results
VARIANTS=(
   "every|all|$DECK|--lts rate2 --mesh-setup every-rank|1.0|$LTS_GOLD|bitwise"
   "rank0|all|$DECK|--lts rate2|1.0|$LTS_GOLD|roundoff"
   "split|1|$DECK|--lts rate2 --mesh-parts-out @BASE/parts --mesh-parts-n @NTASKS|||"
   "parts|all|$DECK|--lts rate2 --mesh-parts @BASE/parts|1.0|@rank0|bitwise"
)
source "$SLURM_SUBMIT_DIR/run_mesh_parts_vista.sh"
