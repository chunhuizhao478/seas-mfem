#!/bin/bash
#SBATCH --job-name=tpv205_100m_p1
#SBATCH --clusters=condo
#SBATCH --partition=scec
#SBATCH --account=scec_608
#SBATCH --nodes=19 --ntasks-per-node=14 --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=24:00:00
#SBATCH --output=tpv205_100m_p1_%j.out
#SBATCH --error=tpv205_100m_p1_%j.err
# Nitsche dynamic rupture, tpv205 on the 100 m mesh, p = 1, 0 -> 12.0 s, Kelvin-Voigt eta_bar 0.1 +
# Prakash-Clifton as committed in tpv205/configs/tpv205_nitsche.toml (Chunhui Zhao, 2026-10-16).
# Submit from miniapps/seas/jobs/condo_nitsche/ AFTER build_verify_condo.sbatch passed at this commit.
# Mesh: tpv205/mesh/tpv2053d_100m.msh, gmsh 4.15.0 output of the committed .geo (Mac, 2026-10-16):
# 5,438,053 tets, 98,984 fault triangles, min 6V/A 11.58 m (one sliver next to the fault; 0.01% of
# tets are below 18.7 m) -> dt 1.09e-4 s, about 110,000 steps to 12 s.
# Layout as TPV104 / condo job 6459052: 14 ranks per node (every rank loads the whole serial mesh,
# drivers/nitsche_dyn_driver.cpp:750; about 1.7 GB per rank for this mesh).
# Estimate: 6459052 (TPV102 100 m, same layout, Prakash-Clifton only) took 0.213 s per step on 7.58 M tets
# and 8 s per fault frame.  Scaled by tets, x2 for Kelvin-Voigt (Frontera 200 m damping runs), plus 600
# fault frames scaled by fault size: about 33,600 + 3,100 s = 10 h.  Limit 24 h.
PROB=tpv205
RES=100m
DECK=tpv205/configs/tpv205_nitsche.toml
MESH=tpv205/mesh/tpv2053d_100m.msh
ORDER=1
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_condo.sh"
