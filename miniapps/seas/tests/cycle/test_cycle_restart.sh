#!/bin/bash
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# test_cycle_restart.sh — fully dynamic BP5 cycle, P1 test test_cycle_restart
# (plan rev 3): stop at step K, restart, and get a bit-identical step 2K.
#
#   full/   one run of 2K steps, checkpoint every K
#   split/  K steps (checkpoint at K), then --restart from it to step 2K
# Every station, global and catalog row and the final checkpoint must be
# identical byte for byte (tests/cycle/compare_runs.py --exact --checkpoints).
#
# Usage (from miniapps/seas, after make -f cycle/cycle.mk seas_fd_cycle_driver):
#   tests/cycle/test_cycle_restart.sh DECK WORKDIR [K=25] [NP=1]
# MPIRUN overrides the launcher (default "mpirun").  Three runs of the
# driver: about 4K QD steps in all.  The prefix is read from the deck's
# [cycle.output] prefix.
# Exit code 0 = PASSED.

set -euo pipefail

DECK=${1:?usage: test_cycle_restart.sh DECK WORKDIR [K] [NP]}
WORK=${2:?usage: test_cycle_restart.sh DECK WORKDIR [K] [NP]}
K=${3:-25}
NP=${4:-1}
MPIRUN=${MPIRUN:-mpirun}
DRIVER=./seas_fd_cycle_driver

PREFIX=$(awk '/^\[cycle\.output\]/{s=1;next} /^\[/{s=0} s&&$1=="prefix"{gsub(/"/,"",$3);print $3;exit}' "$DECK")
[ -n "$PREFIX" ] || { echo "FAILED: no [cycle.output] prefix in $DECK"; exit 1; }
[ -x "$DRIVER" ] || { echo "FAILED: $DRIVER not built"; exit 1; }

rm -rf "$WORK/full" "$WORK/split"
mkdir -p "$WORK/full" "$WORK/split"
N=$((2 * K))

echo "test_cycle_restart: deck $DECK, prefix $PREFIX, K = $K, ranks $NP"
echo "[1/3] full: $N steps, checkpoint every $K"
$MPIRUN -np "$NP" $DRIVER --deck "$DECK" --output-dir "$WORK/full" \
   --max-steps "$N" --checkpoint-every "$K" > "$WORK/full/run.log" 2>&1
echo "[2/3] split: $K steps, checkpoint at $K"
$MPIRUN -np "$NP" $DRIVER --deck "$DECK" --output-dir "$WORK/split" \
   --max-steps "$K" --checkpoint-every "$K" > "$WORK/split/run1.log" 2>&1
cp "$WORK/split/${PREFIX}_checkpoint_r0.txt" "$WORK/split/checkpoint_step${K}_r0.keep"
echo "[3/3] split: restart from step $K to $N"
$MPIRUN -np "$NP" $DRIVER --deck "$DECK" --output-dir "$WORK/split" \
   --max-steps "$N" --checkpoint-every "$K" --restart "$WORK/split/$PREFIX" \
   > "$WORK/split/run2.log" 2>&1

python3 tests/cycle/compare_runs.py "$WORK/full" "$PREFIX" "$WORK/split" "$PREFIX" \
   --exact --checkpoints
