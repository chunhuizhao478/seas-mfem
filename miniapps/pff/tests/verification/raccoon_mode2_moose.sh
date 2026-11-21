#!/usr/bin/env bash
# Phase R2 gate: the full mode II Fx-u curve of pff qs_driver (mode2.toml,
# 1000 steps) against the MOOSE run M-R2 (farms_cdms farms-opt, Mac np 2,
# 2026-11-21; provenance in the gold folder README), and the mesh counts
# (2826 elements, 3026 nodes as the M-R2 and RACCOON gold exodus files).
#
#   tests/verification/raccoon_mode2_moose.sh       (from miniapps/pff, after make qs_driver)
#
# MOOSE caps steps 453-1000 at 20 fixed-point passes (REACH_MAX_ITS,
# accepted), so past the peak the curve compares capped iterations. Observed
# (2026-11-21): 6e-12 over all 1000 steps, the same pass count on every step
# (the same 548 capped steps), peak 423.02086 at t = 0.0097; 14.5 min at np 4.
set -euo pipefail
cd "$(dirname "$0")/../.."
NP=${NP:-4}
REL=${REL:-1e-9}
mpirun -np "$NP" ./qs_driver -c raccoon/mode2_brittle_fracture/configs/mode2.toml -o out/verification/mode2 \
   | tee out/verification/mode2.log | grep -E "elements|step (1|500|1000) "
grep -q "elements 2826" out/verification/mode2.log || { echo "mesh: expected 2826 elements"; exit 1; }
python3 raccoon/compare_curve.py out/verification/mode2/mode2_reaction.csv \
   raccoon/mode2_brittle_fracture/gold/moose_mode2_mac_np2_farms46ad2eb/elasticity_out.csv Fx --rel "$REL"
