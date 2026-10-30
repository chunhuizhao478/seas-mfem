#!/bin/bash
#SBATCH -J sst_tpv102_p3
#SBATCH -o sst_tpv102_p3_%j.out
#SBATCH -e sst_tpv102_p3_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH --ntasks-per-node=144
#SBATCH -t 05:00:00
#SBATCH -A EAR24034
#SBATCH --exclude=i617-131,i617-101
# SeisSol and MFEM speed test on Vista (plan rev1, Chunhui Zhao 2026-10-30): tpv102 on the 200 m mesh
# tpv104/mesh/tpv104_200m.msh, 0 to 12.0 s, 16 nodes, one allocation: MFEM p = 3 (Nitsche, central difference, --lts rate2,
# 128 ranks per node, deck mfem_decks/tpv102_nitsche_sst.toml), then SeisSol o4 (ClusteredLTS 2,
# CFL 0.5, decks/tpv102/, SS_LAYOUTS ranks per node with OpenMP threads).  Both write stations only.
# Submit from miniapps/seas/jobs/seissol_speed_vista/ AFTER build_speed_test_vista.sbatch passed at
# this commit; SS_LAYOUTS is the layout probe's pick.
# Estimate: MFEM 55 to 70 min (gold 1080363, global stepping, 6668 s on 16 x 128; LTS and the speed work 1.7x to 2.1x at p = 2), SeisSol assumed alike: 30 to 37 node-h = 10 to 12 SU; requested 16 x 5 h x 0.33 = 26 SU at most. 128 MFEM ranks per node as the p = 3 gold (memory: 172 to 178 of 238 GB per node).
TAG=sst
PROB=tpv102
MESH=tpv104/mesh/tpv104_200m.msh
PUML=tpv104_200m
ORDER=3
TFINAL=12.0
MFEM_RUN=1
MFEM_RPN=128
SS_LAYOUTS="4"        # set from probe_layout.sbatch's result
SS_SCEC=1
source "$SLURM_SUBMIT_DIR/run_pair_vista.sh"
