#!/bin/bash
#SBATCH -J sst_tpv205_p3
#SBATCH -o sst_tpv205_p3_%j.out
#SBATCH -e sst_tpv205_p3_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH --ntasks-per-node=144
#SBATCH -t 04:00:00
#SBATCH -A EAR24034
#SBATCH --exclude=i617-131,i617-101
# SeisSol and MFEM speed test on Vista (plan rev1, Chunhui Zhao 2026-10-30): tpv205 on the 200 m mesh
# tpv205/mesh/tpv2053d_200m.msh, 0 to 12.0 s, 16 nodes, one allocation: MFEM p = 3 (Nitsche, central difference, --lts rate2,
# 144 ranks per node, deck mfem_decks/tpv205_nitsche_sst.toml), then SeisSol o4 (ClusteredLTS 2,
# CFL 0.5, decks/tpv205/, SS_LAYOUTS ranks per node with OpenMP threads).  Both write stations only.
# Submit from miniapps/seas/jobs/seissol_speed_vista/ AFTER build_speed_test_vista.sbatch passed at
# this commit; SS_LAYOUTS is the layout probe's pick.
# Estimate: MFEM 40 to 60 min (gold 1080361, global stepping, 4663 s; LTS and the speed work 1.7x to 2.1x at p = 2), SeisSol assumed alike: 21 to 32 node-h = 7 to 11 SU; requested 16 x 4 h x 0.33 = 21 SU at most.
TAG=sst
PROB=tpv205
MESH=tpv205/mesh/tpv2053d_200m.msh
PUML=tpv2053d_200m
ORDER=3
TFINAL=12.0
MFEM_RUN=1
MFEM_RPN=144
SS_LAYOUTS="4"        # set from probe_layout.sbatch's result
SS_SCEC=1
source "$SLURM_SUBMIT_DIR/run_pair_vista.sh"
