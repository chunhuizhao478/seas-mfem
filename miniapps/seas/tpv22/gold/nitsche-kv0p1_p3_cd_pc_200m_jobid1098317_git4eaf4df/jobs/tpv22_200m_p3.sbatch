#!/bin/bash
#SBATCH -J tpv22_200m_p3
#SBATCH -o tpv22_200m_p3_%j.out
#SBATCH -e tpv22_200m_p3_%j.err
#SBATCH -p gg
#SBATCH -N 32
#SBATCH -n 2304
#SBATCH -t 08:00:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, SCEC TPV22 on the 200 m mesh, p = 3, 0 -> 15.0 s,
# two vertical faults with an extensional stepover, depth prestress (O1), slip-path weakening (O3), forced rupture (O4), locked fault borders (O5),
# Kelvin-Voigt eta_bar 0.1 + Prakash-Clifton, as committed in tpv22/configs/tpv22_nitsche.toml
# (complex-geometries plan rev 2, phase P5; Chunhui Zhao, 2026-10-21).
# Submit from miniapps/seas/jobs/vista_complex/ AFTER build_verify_complex_vista.sbatch passed at this commit.
# Mesh: tpv22/mesh/tpv22_200m.msh, made on Vista by build_verify_complex_vista.sbatch (3.63 M tets in the cloud);
# 72 ranks per node: every rank reads the whole serial mesh; the build-verify setup check measures that memory.
# Estimate (inferred: Vista TPV205 200 m core-s per tet-step at p = 3, 3.32e-05, times 3.63 M tets and
# 230,426 steps): 7714 core-h, 3.3 h on 32 nodes = 107.1 node-h = 35.4 SU, plus the fault output.
# Requested: 32 nodes x 8 h x 0.33 SU = 84.5 SU at most.
PROB=tpv22
RES=200m
DECK=tpv22/configs/tpv22_nitsche.toml
MESH=tpv22/mesh/tpv22_200m.msh
ORDER=3
TFINAL=15.0
FLAGS=""
source "$SLURM_SUBMIT_DIR/run_complex_vista.sh"
