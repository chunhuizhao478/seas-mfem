#!/bin/bash
#SBATCH -J tpv12_200m_p2
#SBATCH -o tpv12_200m_p2_%j.out
#SBATCH -e tpv12_200m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2304
#SBATCH -t 01:30:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, SCEC TPV12 on the 200 m mesh, p = 2, 0 -> 8.0 s,
# 60-degree normal fault, depth-dependent effective prestress (O1), tension floor (O2), slip-path weakening (O3),
# Kelvin-Voigt eta_bar 0.1 + Prakash-Clifton, as committed in tpv12/configs/tpv12_nitsche.toml
# (complex-geometries plan rev 2, phase P5; Chunhui Zhao, 2026-10-21).
# Submit from miniapps/seas/jobs/vista_complex/ AFTER build_verify_complex_vista.sbatch passed at this commit.
# Mesh: tpv12/mesh/tpv12_200m.msh, made on Vista by build_verify_complex_vista.sbatch (1.88 M tets in the cloud);
# 144 ranks per node: every rank reads the whole serial mesh; the build-verify setup check measures that memory.
# Estimate (inferred: Vista TPV205 200 m core-s per tet-step at p = 2, 1.14e-05, times 1.88 M tets and
# 58,759 steps): 350 core-h, 0.2 h on 16 nodes = 2.4 node-h = 0.8 SU, plus the fault output.
# Requested: 16 nodes x 1.5 h x 0.33 SU = 7.9 SU at most.
PROB=tpv12
RES=200m
DECK=tpv12/configs/tpv12_nitsche.toml
MESH=tpv12/mesh/tpv12_200m.msh
ORDER=2
TFINAL=8.0
FLAGS=""
source "$SLURM_SUBMIT_DIR/run_complex_vista.sh"
