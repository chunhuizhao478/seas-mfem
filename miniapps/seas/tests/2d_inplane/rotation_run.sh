#!/bin/bash
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# tests/2d_inplane/rotation_run.sh — 2D in-plane design rev 3, phase A: the
# BP5-mini rotation run.  The whole problem rotated rigidly by R (40 degrees
# about (1, 2, 3)) must give the same fault-local and axis-projected
# outputs, the same time-step sequence and the same switch time as the
# unrotated one, to roundoff.
#
#   tests/2d_inplane/rotation_run.sh BIN OUT       (from miniapps/seas)
#
# BIN holds seas_fd_cycle_driver.  The decks are the parity set's cycle decks
# (parity_3d.sh decks: BP5-FD at p = 1 on the 1000 m mesh, mumps_blr), each
# run three times:
#   plate_<run>  unrotated, 1 rank, with [cycle.loading] kind = "plate" and
#                BP5's vectors (equal to the default "bp5" loading bit for
#                bit: test_plate_loading L1, and the parity run when PARITY
#                is set)
#   rot_<run>    1 rank, the mesh's nodes and every deck vector rotated by R
#                (rotate_case.py)
#   ctl_<run>    the plate deck on 2 ranks: the control, a change of
#                summation order only
# Runs: qd      the QD regime only, 10 steps (--max-steps 10)
#       dyn     the first event dynamic from t = 0, 40 dynamic steps
#       switch  V_on = 3.5e-2 m/s, above the patch's initial 3e-2 m/s: QD
#               steps up to the first switch, then 40 dynamic steps.  (In the
#               QD regime max V rises to 4.2e-2 m/s at t = 2.0 s and then
#               falls, cloud run 10-25, so V_on must lie below that; it is
#               crossed after about 40 steps.  --max-steps 300 bounds the QD
#               part, and both runs must log the switch at t > 0.)
# Compared with rotation_compare.py (the vectors' components together,
# global vectors as R times the unrotated run's): first plate against ctl,
# which gives the size C that roundoff alone reaches (the QD regime's
# adaptive steps carry it into the step sizes: about 6e-10 after 10 steps,
# cloud runs 10-25); then plate against rot --rotated at
# the tolerance 100 C (at least 1e-12).  The control must stay below 1e-6.
# With PARITY=<a parity_3d.sh run OUT> the plate runs are also compared
# with the parity runs' cycle_<run>_np1 bit for bit.
# Environment: MPIRUN (default "mpirun --oversubscribe"; add
# --allow-run-as-root as root), ONLY (space-separated run names).  The 1-rank
# mumps_blr QD solve takes about 9 GB; one run at a time, about 4 min each
# (qd, dyn) and 8 min (switch) in the cloud container.
#
# Scope (Chunhui Zhao, 2026-10-25).
set -u
export OMP_NUM_THREADS=1 LC_ALL=C
BIN=$(cd "${1:?usage: rotation_run.sh BIN OUT}" && pwd); OUT=${2:?usage: rotation_run.sh BIN OUT}
MPIRUN=${MPIRUN:-mpirun --oversubscribe}
ONLY=${ONLY:-}
[ -e "$OUT" ] && { echo "$OUT exists: not replacing it"; exit 2; }
HERE=tests/2d_inplane
$HERE/parity_3d.sh decks "$OUT" > /dev/null || exit 2
D=$OUT/decks
python3 - "$D" <<'PY' || exit 2
import sys
d = sys.argv[1]
s = open(f"{d}/cycle_dyn.toml").read()
a = "\nV_on               = 1.0e-2"
assert s.count(a) == 1, a
open(f"{d}/cycle_switch.toml", "w").write(s.replace(a, "\nV_on               = 3.5e-2", 1))
PY
python3 $HERE/rotate_case.py mesh bp5/mesh/reference/bp5_tandem_exact.msh "$OUT/bp5_rotated.msh" || exit 2
for r in qd dyn switch; do
   python3 $HERE/rotate_case.py plate "$D/cycle_$r.toml" "$D/plate_$r.toml" || exit 2
   python3 $HERE/rotate_case.py deck "$D/cycle_$r.toml" "$D/rot_$r.toml" "$OUT/bp5_rotated.msh" || exit 2
done
# DOFs per element of the decks' order (the checkpoint's displacement)
NDOF=$(python3 - "$D/cycle_qd.toml" <<'PY'
import re, sys
s = open(sys.argv[1]).read()
m = re.search(r"(?ms)^\[mesh\]\s*$(.*?)^\[", s)
p = int(re.search(r"(?m)^order\s*=\s*(\d+)", m.group(1)).group(1))
print((p + 1) * (p + 2) * (p + 3) // 6)
PY
) || exit 2

want() { [ -z "$ONLY" ] || [[ " $ONLY " == *" $1 "* ]]; }
fails=0
for r in qd dyn switch; do
   want $r || continue
   case $r in
      qd)     opt="--max-steps 10" ;;
      dyn)    opt="--max-dyn-steps 40" ;;
      switch) opt="--max-steps 300 --max-dyn-steps 40" ;;
   esac
   for v in plate rot ctl; do
      deck=$D/${v}_$r.toml; np=1
      [ $v = ctl ] && { deck=$D/plate_$r.toml; np=2; }
      mkdir -p "$OUT/${v}_$r"
      t0=$(date +%s)
      $MPIRUN -np $np "$BIN/seas_fd_cycle_driver" --deck "$deck" \
         --output-dir "$OUT/${v}_$r" $opt < /dev/null > "$OUT/${v}_$r/run.log" 2>&1
      rc=$?
      echo "${v}_$r: exit $rc in $(( $(date +%s) - t0 )) s on $np rank(s)"
      [ $rc -eq 0 ] || fails=$((fails + 1))
      if [ $r = switch ]; then
         sw=$(awk '!/^#/ && $3 == "QD" && $4 == "dynamic" && $1 > 0 {print $1; exit}' \
              "$OUT/${v}_$r/b_switch.txt" 2>/dev/null)
         if [ -n "$sw" ]; then echo "   switch to dynamic at t = $sw s"
         else echo "   FAILED: no QD-to-dynamic switch at t > 0"; fails=$((fails + 1)); fi
      fi
   done
   echo "control: plate_$r against ctl_$r (1 rank against 2):"
   python3 $HERE/rotation_compare.py "$OUT/plate_$r" "$OUT/ctl_$r" --skip-checkpoint --tol 1e-6 \
      > "$OUT/control_$r.txt"
   crc=$?
   tail -1 "$OUT/control_$r.txt"
   C=$(awk '/^   largest /{print $2}' "$OUT/control_$r.txt")
   if [ $crc -ne 0 ] || [ -z "$C" ]; then
      echo "   FAILED: the control is not below 1e-6 (see $OUT/control_$r.txt)"
      fails=$((fails + 1)); continue
   fi
   tol=$(python3 -c "print('%.3g' % max(100 * float('$C'), 1e-12))")
   echo "compare plate_$r rot_$r (rotated), tolerance 100 x control = $tol:"
   python3 $HERE/rotation_compare.py "$OUT/plate_$r" "$OUT/rot_$r" --rotated \
      --elem-dofs "$NDOF" --tol "$tol" || fails=$((fails + 1))
   if [ -n "${PARITY:-}" ] && [ -d "$PARITY/cycle_${r}_np1" ]; then
      echo "compare $PARITY/cycle_${r}_np1 plate_$r (bit for bit):"
      $HERE/parity_3d.sh compare "$PARITY/cycle_${r}_np1" "$OUT/plate_$r" | sed 's/^/   /'
      [ ${PIPESTATUS[0]} -eq 0 ] || fails=$((fails + 1))
   fi
done
echo "rotation_run: $fails failure(s)"
exit $fails
