#!/bin/bash
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# tests/2d_inplane/relaxation_decks.sh — the driver side of test-relaxation
# ([cycle.relaxation] t_relax_yr, cycle/normal_stress_relaxation.hpp; Chunhui
# Zhao, 2026-11-02: "implement the relaxation, test and verify it").  Copies
# of the BP3-FD d60 thrust deck on its 1 km mesh, run through
# seas_fd_cycle_driver on 1 rank:
#   D1  each deck refusal exits non-zero with its message: t_relax_yr = 0,
#       with [cycle.probe], with [friction.prakash_clifton]
#   D2  t_relax_yr = 20 is accepted and reported; the checkpoint of 3 QD
#       steps is FD_CYCLE_V5 with a relaxation record; without the section
#       it stays FD_CYCLE_V4 with none
#   D3  restart: a run to 150 yr (t_relax_yr = 20, through event 1 at about
#       130 yr, its switch back and QD after it) against the same run
#       stopped at QD step 1000 (a checkpoint step of both) and restarted
#       from that checkpoint: every station, catalog and switch file the
#       same, byte for byte, and both runs report the same relaxation at
#       the switch and at the end
#
# Usage (from miniapps/seas, after building seas_fd_cycle_driver):
#   MPIEXEC=mpirun MPIEXEC_NP=-np tests/2d_inplane/relaxation_decks.sh
# Exit code 0 = PASSED.

set -u
MPI=${MPIEXEC:-mpirun}
NP=${MPIEXEC_NP:--np}
BP3=config/fd_cycle/bp3_fd_d60_thrust.toml
T=${TMPDIR:-/tmp}/seas_test_relaxation
P=bp3_fd_d60_thrust
rm -rf "$T"
mkdir -p "$T"
fail=0

# base NAME TFINAL: the BP3 deck on the 1 km mesh to TFINAL s
base()
{
   sed -e 's|^path  = "bp3/mesh/bp3_d60_25m.msh"$|path  = "bp3/mesh/bp3_d60_1km.msh"|' \
       -e "s|^tfinal     = 4.7304e10$|tfinal     = $2|" "$BP3" > "$T/$1.toml"
   grep -q '^path  = "bp3/mesh/bp3_d60_1km.msh"$' "$T/$1.toml" \
      && grep -q "^tfinal     = $2$" "$T/$1.toml" \
      || { echo "  FAILED: the BP3 deck lines moved"; exit 1; }
}
# add NAME BEFORE TEXT: TEXT (printf %b) inserted before the line BEFORE
add()
{
   python3 - "$T/$1.toml" "$2" "$3" <<'EOF'
import sys
p, before, text = sys.argv[1], sys.argv[2], sys.argv[3].encode().decode('unicode_escape')
s = open(p).read()
assert s.count(before + "\n") == 1, before
open(p, "w").write(s.replace(before + "\n", text + before + "\n", 1))
EOF
}
run() { $MPI $NP 1 ./seas_fd_cycle_driver -d "$T/$1.toml" -o "$T/$1_out" "${@:2}" 2>&1; }

RELAX='[cycle.relaxation]\nt_relax_yr = 20.0\n\n'

# ---- D1: refusals
refuse()
{
   out=$(run "$1")
   rc=$?
   if [ $rc -ne 0 ] && printf '%s' "$out" | grep -qF -- "$2"; then
      echo "  PASSED: D1 deck refusal $1: exit $rc, \"$2\""
   else
      printf '%s\n' "$out" | tail -15
      echo "  FAILED: D1 deck refusal $1: exit $rc, \"$2\" not found"
      fail=1
   fi
}
base zero 4.7304e10
add zero "[cycle.switch]" '[cycle.relaxation]\nt_relax_yr = 0.0\n\n'
refuse zero "[cycle.relaxation].t_relax_yr = 0.0000000000000000e+00 must be > 0"
base probe 4.7304e10
add probe "[cycle.switch]" "${RELAX}"'[cycle.probe]\nV_off = [1.0e-3]\nT_dwell = [10.0]\n\n'
refuse probe "[cycle.probe] with [cycle.relaxation]"
base pc 4.7304e10
add pc "[cycle.switch]" "${RELAX}"
add pc "[cycle]" '[friction.prakash_clifton]\nv_star = 1.0\nlength = 1.0e-3\n\n'
refuse pc "the relaxation is not combined with Prakash-Clifton"

# ---- D2: accepted, the checkpoint versions
base ok 4.7304e10
add ok "[cycle.switch]" "${RELAX}"
out=$(run ok --max-steps 3)
rc=$?
ck="$T/ok_out/${P}_checkpoint_r0.txt"
if [ $rc -eq 0 ] && printf '%s' "$out" | grep -qF "sigma_n' relaxation: t_r = 20 yr ([cycle.relaxation])" \
   && printf '%s' "$out" | grep -qF "sigma_n' relaxation at the end (t = " \
   && grep -qx "FD_CYCLE_V5" "$ck" && grep -q "^relax 6\.311520*e+08 " "$ck"; then
   echo "  PASSED: D2 t_relax_yr = 20: accepted, reported, checkpoint FD_CYCLE_V5 with the record"
else
   printf '%s\n' "$out" | tail -15
   echo "  FAILED: D2 t_relax_yr = 20: exit $rc"
   fail=1
fi
base off 4.7304e10
out=$(run off --max-steps 3)
rc=$?
ck="$T/off_out/${P}_checkpoint_r0.txt"
if [ $rc -eq 0 ] && grep -qx "FD_CYCLE_V4" "$ck" && ! grep -q "^relax " "$ck" \
   && ! printf '%s' "$out" | grep -qF "sigma_n' relaxation"; then
   echo "  PASSED: D2 without [cycle.relaxation]: checkpoint FD_CYCLE_V4, no record, no report"
else
   printf '%s\n' "$out" | tail -15
   echo "  FAILED: D2 without [cycle.relaxation]: exit $rc"
   fail=1
fi

# ---- D3: restart, bit for bit
base full 4.73e9
add full "[cycle.switch]" "${RELAX}"
cp "$T/full.toml" "$T/part.toml"
out_full=$(run full)
rc_full=$?
out_a=$($MPI $NP 1 ./seas_fd_cycle_driver -d "$T/part.toml" -o "$T/part_out" --max-steps 1000 2>&1)
rc_a=$?
out_b=$($MPI $NP 1 ./seas_fd_cycle_driver -d "$T/part.toml" -o "$T/part_out" \
        --restart "$T/part_out/$P" 2>&1)
rc_b=$?
echo "  D3: full run exit $rc_full, stopped at step 1000 exit $rc_a, restarted exit $rc_b"
nsw=$(printf '%s' "$out_full" | grep -c "back to QD at")
diffs=0
nf=0
for f in "$T/full_out/${P}"_{fltst,srfst}*.txt "$T/full_out/${P}"_{catalog,events,switch,global,fault}.txt; do
   [ -f "$f" ] || continue
   nf=$((nf + 1))
   cmp -s "$f" "$T/part_out/$(basename "$f")" || { echo "    differs: $(basename "$f")"; diffs=$((diffs + 1)); }
done
rl_full=$(printf '%s' "$out_full" | grep "sigma_n' relaxation \(before\|at the end\)")
rl_b=$(printf '%s' "$out_b" | grep "sigma_n' relaxation \(before\|at the end\)")
if [ $rc_full -eq 0 ] && [ $rc_a -eq 0 ] && [ $rc_b -eq 0 ] && [ "$nsw" -ge 1 ] && [ $nf -ge 10 ] \
   && [ $diffs -eq 0 ] && [ -n "$rl_full" ] && [ "$rl_full" = "$rl_b" ]; then
   echo "  PASSED: D3 restart at QD step 1000 through $nsw switch back(s): $nf output files the same," \
        "byte for byte; the relaxation lines the same"
   printf '%s\n' "$rl_full" | sed 's/^/    /'
else
   printf '%s\n' "$out_b" | tail -15
   echo "  FAILED: D3 restart: $diffs of $nf files differ, $nsw switch backs"
   printf 'full:\n%s\nrestarted:\n%s\n' "$rl_full" "$rl_b"
   fail=1
fi

[ $fail -eq 0 ] && rm -rf "$T"
exit $fail
