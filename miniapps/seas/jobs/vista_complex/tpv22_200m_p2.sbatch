#!/bin/bash
#SBATCH -J tpv22_200m_p2
#SBATCH -o tpv22_200m_p2_%j.out
#SBATCH -e tpv22_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 32
#SBATCH -n 2304
#SBATCH -t 02:30:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, SCEC TPV22 on the 200 m mesh, p = 2, 0 -> 15.0 s,
# two vertical faults with an extensional stepover, depth prestress (O1), slip-path weakening (O3), forced rupture (O4), locked fault borders (O5),
# Kelvin-Voigt eta_bar 0.1 + Prakash-Clifton, as committed in tpv22/configs/tpv22_nitsche.toml
# (complex-geometries plan rev 2, phase P5; Chunhui Zhao, 2026-10-21).
# Submit from miniapps/seas/jobs/vista_complex/ AFTER build_verify_complex_vista.sbatch passed at this commit.
# Mesh: tpv22/mesh/tpv22_200m.msh, made on Vista by build_verify_complex_vista.sbatch (3.63 M tets in the cloud);
# 72 ranks per node: every rank reads the whole serial mesh; the build-verify setup check measures that memory.
# Estimate (inferred: Vista TPV205 200 m core-s per tet-step at p = 2, 1.14e-05, times 3.63 M tets and
# 135,958 steps): 1563 core-h, 0.7 h on 32 nodes = 21.7 node-h = 7.2 SU, plus the fault output.
# Requested: 32 nodes x 2.5 h x 0.33 SU = 26.4 SU at most.
PROB=tpv22
RES=200m
DECK=tpv22/configs/tpv22_nitsche.toml
MESH=tpv22/mesh/tpv22_200m.msh
ORDER=2
TFINAL=15.0
FLAGS=""
source "$SLURM_SUBMIT_DIR/run_complex_vista.sh"
