#!/bin/bash
#SBATCH -J tpv15_200m_p1
#SBATCH -o tpv15_200m_p1_%j.out
#SBATCH -e tpv15_200m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 6
#SBATCH -n 576
#SBATCH -t 02:00:00
#SBATCH -A EAR26024
# Nitsche dynamic rupture, SCEC TPV15 on the 200 m mesh, p = 1, 0 -> 12.0 s,
# left-lateral main fault with a 30-degree branch 100 m from the junction, per-fault prestress boxes, slip-path weakening (O3),
# Kelvin-Voigt eta_bar 0.1 + Prakash-Clifton, as committed in tpv15/configs/tpv15_nitsche.toml
# (complex-geometries plan rev 2, phase P5; Chunhui Zhao, 2026-10-21).
# Submit from miniapps/seas/jobs/vista_complex/ AFTER build_verify_complex_vista.sbatch passed at this commit.
# Mesh: tpv15/mesh/tpv15_200m.msh, made on Vista by build_verify_complex_vista.sbatch (2.53 M tets in the cloud);
# 96 ranks per node: every rank reads the whole serial mesh; the build-verify setup check measures that memory.
# --lts rate2 (plan Q1): the mesh report gives an ideal LTS speedup of 7.89x (cloud, 200 m); the estimate
# assumes 2.5x (TPV205 reached 1.30x of its ideal 3.20x), so it is 0.5 SU, and 1.2 SU with no LTS gain.
# Estimate (inferred: Vista TPV205 200 m core-s per tet-step at p = 1, 2.39e-06, times 2.53 M tets and
# 202,703 steps, over the assumed LTS gain 2.5): 136 core-h, 0.2 h on 6 nodes = 1.4 node-h = 0.5 SU, plus the fault output.
# Requested: 6 nodes x 2 h x 0.33 SU = 4.0 SU at most.
PROB=tpv15
RES=200m
DECK=tpv15/configs/tpv15_nitsche.toml
MESH=tpv15/mesh/tpv15_200m.msh
ORDER=1
TFINAL=12.0
FLAGS="--lts rate2"
source "$SLURM_SUBMIT_DIR/run_complex_vista.sh"
