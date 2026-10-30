#!/bin/bash
#SBATCH -J sst_tpv102_p1
#SBATCH -o sst_tpv102_p1_%j.out
#SBATCH -e sst_tpv102_p1_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH --ntasks-per-node=144
#SBATCH -t 01:00:00
#SBATCH -A EAR24034
#SBATCH --exclude=i617-131,i617-101
# SeisSol and MFEM speed test on Vista (plan rev1, Chunhui Zhao 2026-10-30): tpv102 on the 200 m mesh
# tpv104/mesh/tpv104_200m.msh, 0 to 12.0 s, 4 nodes, one allocation: MFEM p = 1 (Nitsche, central difference, --lts rate2,
# 144 ranks per node, deck mfem_decks/tpv102_nitsche_sst.toml), then SeisSol o2 (ClusteredLTS 2,
# CFL 0.5, decks/tpv102/, SS_LAYOUTS ranks per node with OpenMP threads).  Both write stations only.
# Submit from miniapps/seas/jobs/seissol_speed_vista/ AFTER build_speed_test_vista.sbatch passed at
# this commit; SS_LAYOUTS is the layout probe's pick.
# Estimate: MFEM about 5 min (TPV104 LTS gold 1082587 p1 timing scaled by the V2 gain), SeisSol assumed alike: 0.7 node-h = 0.2 SU; requested 4 x 1 h x 0.33 = 1.3 SU at most.
TAG=sst
PROB=tpv102
MESH=tpv104/mesh/tpv104_200m.msh
PUML=tpv104_200m
ORDER=1
TFINAL=12.0
MFEM_RUN=1
MFEM_RPN=144
SS_LAYOUTS="4"        # set from probe_layout.sbatch's result
SS_SCEC=1
source "$SLURM_SUBMIT_DIR/run_pair_vista.sh"
