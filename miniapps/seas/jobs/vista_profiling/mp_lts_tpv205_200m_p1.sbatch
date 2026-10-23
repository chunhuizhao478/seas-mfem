#!/bin/bash
#SBATCH -J mp_lts_tpv205_200m_p1
#SBATCH -o mp_lts_tpv205_200m_p1_%j.out
#SBATCH -e mp_lts_tpv205_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 00:45:00
#SBATCH -A EAR26024
#SBATCH --exclude=i617-131,i617-101
# Mesh setup check (partition-memory report rev1 options D and E, speed-next S1; Chunhui Zhao,
# 2026-10-22): TPV205 on the 200 m gold mesh, p = 1, the rank layout of LTS gold 1082586 (4 nodes x 144
# ranks).  Five variants through run_mesh_parts_vista.sh:
#   every  --lts rate2 --mesh-setup every-rank (the path before this change), 0 to 3.0 s: bitwise equal
#          to the first 3.0 s of 1082586; its setup time and memory are the reference
#   rank0  --lts rate2 (rank 0 reads the mesh and sends each rank its part, the new default), 0 to
#          12.0 s against 1082586 (the shared faces come in another order, so round-off differences)
#   split  one rank: --lts rate2 --mesh-parts-out <dir> --mesh-parts-n 576 writes the 576 parts
#   parts  --lts rate2 --mesh-parts <dir> (every rank reads its own part), 0 to 3.0 s: bitwise equal
#          to the first 3.0 s of rank0
#   gts    global stepping with the new default, 0 to 1.0 s against the first 1.0 s of gold 1076378
# Submit from miniapps/seas/jobs/vista_profiling/ AFTER build_mesh_parts_vista.sbatch passed at this commit.
# Estimate: four setups of about 80 s, the split about 2 min, stepping about 6 min (12 s) + 2 x 1.5 min
# (3 s) + 0.5 min (1 s): about 18 min = 1.2 node-h = 0.4 SU.  Requested: 4 nodes x 0.75 h x 0.33 SU =
# 1.0 SU at most.
PROB=tpv205
MESH=tpv205/mesh/tpv2053d_200m.msh
ORDER=1
DECK=tpv205/configs/tpv205_nitsche.toml
LTS_GOLD=tpv205/gold/nitsche-kv0p1_p1_cd_pc_lts-rate2-buf3_200m_jobid1082586_gitac3adeb/results
GTS_GOLD=tpv205/gold/nitsche-kv0p1_p1_cd_pc_200m_jobid1076378_git81759da/results
VARIANTS=(
   "every|all|$DECK|--lts rate2 --mesh-setup every-rank|3.0|$LTS_GOLD|bitwise"
   "rank0|all|$DECK|--lts rate2|12.0|$LTS_GOLD|roundoff"
   "split|1|$DECK|--lts rate2 --mesh-parts-out @BASE/parts --mesh-parts-n @NTASKS|||"
   "parts|all|$DECK|--lts rate2 --mesh-parts @BASE/parts|3.0|@rank0|bitwise"
   "gts|all|$DECK||1.0|$GTS_GOLD|roundoff"
)
source "$SLURM_SUBMIT_DIR/run_mesh_parts_vista.sh"
