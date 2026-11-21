#!/usr/bin/env bash
# Phase R1 gate: the full mode I Fy-u curve of pff qs_driver (mode1.toml,
# 175 steps) against the MOOSE run M-R1 (farms_cdms farms-opt, Mac np 4,
# 2026-11-21; provenance in the gold folder README).
#
#   tests/verification/raccoon_mode1_moose.sh       (from miniapps/pff, after make qs_driver)
#
# Gate: relative difference <= 1e-4 on every step with |Fy| > 1, and the
# peak at the same step. Observed: 2.0e-5, peak 691.00227 at t = 2.86e-3 in
# both (equal to 1e-11; the steps where both codes hit the 20-pass cap agree
# to round-off, the rest differ only through the number of passes).
set -euo pipefail
cd "$(dirname "$0")/../.."
NP=${NP:-1}
mpirun -np "$NP" ./qs_driver -c raccoon/mode1_brittle_fracture/configs/mode1.toml -o out/verification/mode1
python3 raccoon/compare_curve.py out/verification/mode1/mode1_reaction.csv \
   raccoon/mode1_brittle_fracture/gold/moose_mode1_mac_np4_farms46ad2eb-dirty/elasticity_out.csv Fy --rel 1e-4
