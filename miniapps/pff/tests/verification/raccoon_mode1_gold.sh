#!/usr/bin/env bash
# Phase R1 gate: pff qs_driver against the RACCOON mode1_brittle_fracture gold
# file (hugary1995/raccoon 9fd3a6f, tutorials/mode1_brittle_fracture/gold).
#
#   tests/verification/raccoon_mode1_gold.sh        (from miniapps/pff, after make qs_driver)
#
# 1. mode1_gold.toml (fixed-point iterations as in the input, with MOOSE's
#    automatic scaling in the fixed-point test, which then stops after one
#    pass at steps 1 and 2 as MOOSE does): rel 1e-8 (the RACCOON test uses
#    exodiff rel 1e-4, abs_zero 1e-9).
# 2. mode1_gold_single_pass.toml (one staggered pass per step): rel 1e-8.
# Needs python3 with numpy and netCDF4.
set -euo pipefail
cd "$(dirname "$0")/../.."
GOLD=raccoon/mode1_brittle_fracture/benchmark_data/raccoon_9fd3a6f/elasticity_out.e
fail=0
run() {
   local cfg=$1 name=$2 rel=$3
   ./qs_driver -c raccoon/mode1_brittle_fracture/configs/$cfg.toml -o out/verification/$name
   for s in 1 2; do
      python3 raccoon/compare_gold.py "$GOLD" out/verification/$name/${name}_nodal_step$s.txt $s \
         --rel "$rel" --abs-zero 1e-9 || fail=1
   done
}
run mode1_gold mode1_gold 1e-8
run mode1_gold_single_pass mode1_gold_single_pass 1e-8
if [ $fail -ne 0 ]; then echo "raccoon_mode1_gold: FAILED"; exit 1; fi
echo "raccoon_mode1_gold: passed"
