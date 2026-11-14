#!/bin/bash
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# tests/2d_inplane/relaxed_opening_decks.sh — the deck side of
# [cycle.relaxation] opening = "relaxed" (cycle/normal_stress_relaxation.hpp;
# Chunhui Zhao, 2026-11-13: "bring the relaxation review").  Copies of the
# D9f option I deck through seas_fd_cycle_driver on 1 rank:
#   E1  each deck refusal exits non-zero with its message: the relaxation
#       with integrator = "esdirk" and no opening key; opening = "relaxed"
#       with the DP5 integrator; an unknown opening value
#   E2  opening = "relaxed" with integrator = "esdirk" is accepted and
#       reported, 3 QD steps run, and the checkpoint carries the relaxation
#       record
#
# Usage (from miniapps/seas, after building seas_fd_cycle_driver):
#   MPIEXEC=mpirun MPIEXEC_NP=-np tests/2d_inplane/relaxed_opening_decks.sh
# Exit code 0 = PASSED.

set -u
MPI=${MPIEXEC:-mpirun}
NP=${MPIEXEC_NP:--np}
BASE=simple_cases/configs/D9f_arc08L2p5_h250_avg_floor10_contact_esdirk.toml
RLX=simple_cases/configs/D9f_arc08L2p5_h250_avg_floor10_contact_esdirk_relax5k.toml
P=D9f_arc08L2p5_h250_avg_floor10_contact_esdirk_relax5k
T=${TMPDIR:-/tmp}/seas_test_relaxed_opening
rm -rf "$T"
mkdir -p "$T"
fail=0

# edit NAME FROM OLD NEW: FROM with the one line OLD replaced by NEW (printf %b)
edit()
{
   python3 - "$2" "$T/$1.toml" "$3" "$4" <<'PY'
import sys
src, dst, old, new = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].encode().decode('unicode_escape')
s = open(src).read()
assert s.count(old + "\n") == 1, old
open(dst, "w").write(s.replace(old + "\n", new + "\n", 1))
PY
}
run() { $MPI $NP 1 ./seas_fd_cycle_driver -d "$T/$1.toml" -o "$T/$1_out" "${@:2}" 2>&1; }
refuse()
{
   out=$(run "$1" --max-steps 1)
   rc=$?
   if [ $rc -ne 0 ] && printf '%s' "$out" | grep -qF -- "$2"; then
      echo "  PASSED: E1 deck refusal $1: exit $rc, \"$2\""
   else
      printf '%s\n' "$out" | tail -15
      echo "  FAILED: E1 deck refusal $1: exit $rc, \"$2\" not found"
      fail=1
   fi
}

# ---- E1: refusals
edit nokey "$RLX" 'opening    = "relaxed"' ''
refuse nokey 'or with [cycle.relaxation] opening = "relaxed"'
edit dp5 "$RLX" 'integrator    = "esdirk"' 'integrator    = "dp5"'
refuse dp5 '[cycle.relaxation] opening = "relaxed" needs [cycle.qd] integrator = "esdirk"'
edit bad "$RLX" 'opening    = "relaxed"' 'opening    = "both"'
refuse bad '[cycle.relaxation].opening = "both"'

# ---- E2: accepted, reported, the checkpoint's record
cp "$RLX" "$T/ok.toml"
out=$(run ok --max-steps 3)
rc=$?
ck="$T/ok_out/${P}_checkpoint_r0.txt"
if [ $rc -eq 0 ] && printf '%s' "$out" | grep -qF "sigma_n' relaxation: t_r = 5000 yr ([cycle.relaxation])" \
   && printf '%s' "$out" | grep -qF 'opening = "relaxed"' \
   && printf '%s' "$out" | grep -qF "sigma_n' relaxation at the end (t = " \
   && grep -q "^relax 1\.577880*e+11 " "$ck"; then
   echo "  PASSED: E2 opening = \"relaxed\" with esdirk: accepted, reported, checkpoint with the record"
else
   printf '%s\n' "$out" | tail -15
   echo "  FAILED: E2 opening = \"relaxed\" with esdirk: exit $rc"
   fail=1
fi

[ $fail -eq 0 ] && echo "relaxed_opening_decks: PASSED" || echo "relaxed_opening_decks: FAILED"
exit $fail
