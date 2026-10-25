#!/bin/bash
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# tests/2d_inplane/cycle_tests.sh — the fully dynamic cycle's test targets
# (cycle/cycle.mk: test-cycle-api-pin, test-cycle-p2, -p3, -p4, with their
# default rank counts and decks) run with the binaries in BIN, one log per
# run in OUT, so the branch point's and the branch's binaries go through the
# same commands and their printed measurements can be compared
# (parity_3d.sh compare reports differing logs).
#
#   tests/2d_inplane/cycle_tests.sh BIN OUT     (from miniapps/seas)
# Environment: MPIRUN (default "mpirun --oversubscribe"); BLR=1 runs a copy
# of the deck with solver = "mumps_blr" (written to OUT/deck.toml), for a
# host where full-rank MUMPS on the 1000 m BP5 mesh does not fit (it asks for
# 18 GB at 1 rank); every run then reads that copy.
set -u
BIN=$(cd "${1:?usage: cycle_tests.sh BIN OUT}" && pwd); OUT=${2:?}
MPIRUN=${MPIRUN:-mpirun --oversubscribe}
DECK=config/fd_cycle/bp5_fd_first_event.toml
[ -e "$OUT" ] && { echo "$OUT exists: not replacing it"; exit 2; }
mkdir -p "$OUT"
if [ "${BLR:-0}" = 1 ]; then
   python3 - "$DECK" "$OUT/deck.toml" <<'PY' || exit 2
import sys
s = open(sys.argv[1]).read()
a = 'solver        = "mumps"'
assert s.count("\n" + a) == 1, a
open(sys.argv[2], "w").write(s.replace("\n" + a, '\nsolver        = "mumps_blr"', 1))
PY
   DECK=$OUT/deck.toml
fi
export OMP_NUM_THREADS=1
fails=0
run() {   # name, ranks, binary, options...
   local n=$1 np=$2 b=$3; shift 3
   local s=$(date +%s)
   $MPIRUN -np "$np" "$BIN/$b" "$@" < /dev/null > "$OUT/$n.log" 2>&1
   local rc=$?
   echo "exit $rc" >> "$OUT/$n.log"
   echo "$n: exit $rc in $(( $(date +%s) - s )) s on $np rank(s)"
   [ $rc -eq 0 ] || { fails=$((fails + 1)); tail -15 "$OUT/$n.log"; }
}
run api_pin 1 seas_test_cycle_api_pin
for n in 1 2 4 8; do
   run p2_frames_np$n $n seas_test_cycle_p2 --deck $DECK --frames-only
   [ $((n % 2)) -eq 0 ] && run p2_frames_split_np$n $n seas_test_cycle_p2 --deck $DECK --frames-only --split-y
done
run p2_np2 2 seas_test_cycle_p2 --deck $DECK
for n in 1 2; do
   run p3_np$n $n seas_test_cycle_p3 --deck $DECK --work-dir "$OUT/p3_work/np$n"
   run p4_np$n $n seas_test_cycle_p4 --deck $DECK
done
echo "cycle_tests: $fails failed run(s)"
exit $fails
