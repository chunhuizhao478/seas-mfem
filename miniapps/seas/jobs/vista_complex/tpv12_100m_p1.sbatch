#!/bin/bash
#SBATCH -J tpv12_100m_p1
#SBATCH -o tpv12_100m_p1_%j.out
#SBATCH -e tpv12_100m_p1_%j.err
#SBATCH -p gg
#SBATCH -N 8
#SBATCH -n 512
#SBATCH -t 04:00:00
#SBATCH -A EAR26024
# Nitsche dynamic rupture, SCEC TPV12 on the 100 m mesh, p = 1, 0 -> 8.0 s,
# 60-degree normal fault, depth-dependent effective prestress (O1), tension floor (O2), slip-path weakening (O3),
# Kelvin-Voigt eta_bar 0.1 + Prakash-Clifton, as committed in tpv12/configs/tpv12_nitsche.toml
# (complex-geometries plan rev 2, phase P5; Chunhui Zhao, 2026-10-21).
# Submit from miniapps/seas/jobs/vista_complex/ AFTER build_verify_complex_vista.sbatch passed at this commit.
# Mesh: made on Vista by make_tpv12_100m_mesh_vista.sbatch (about 5.4 M tets, inferred from TPV205's
# 200 m to 100 m ratio); 64 ranks per node as the TPV205 100 m gold (about 2.3 GB per rank for the serial mesh).
# Estimate (inferred: Vista TPV205 200 m core-s per tet-step at p = 1, 2.39e-06, times 5.40 M tets and
# 55,173 steps): 198 core-h, 0.4 h on 8 nodes = 3.1 node-h = 1.0 SU, plus the fault output.
# Requested: 8 nodes x 4 h x 0.33 SU = 10.6 SU at most.
PROB=tpv12
RES=100m
DECK=tpv12/configs/tpv12_nitsche.toml
MESH=tpv12/mesh/tpv12_100m.msh
ORDER=1
TFINAL=8.0
FLAGS=""
source "$SLURM_SUBMIT_DIR/run_complex_vista.sh"
