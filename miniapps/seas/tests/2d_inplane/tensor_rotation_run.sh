#!/bin/bash
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# tests/2d_inplane/tensor_rotation_run.sh — the cycle's tensor prestress on a
# rotated planar fault line (Chunhui Zhao, 2026-10-29: "test on the rotated
# planar fault line (modified from TPV205-2D) in a certain angle and project
# mesh"; "we should keep the function the same as in dynamic phase, to make
# it consistent and work for 2D and 3D").  The map-view deck
# config/fd_cycle/mapview_rs_tensor_fd.toml (TPV205-2D's layout, rate and
# state, [stress] kind = "constant_tensor") turned rigidly by DEG degrees
# about z, the mesh's nodes and every deck vector with it and the tensor as
# R S R^T, must give the same fault-local and axis-projected outputs, the
# same time-step sequence and the same switch times as the unturned deck, to
# roundoff: the tensor's traction on the turned fault line is the same, by
# the dynamic operator's own projection (domain/fault_stress_projection.hpp).
#
#   tests/2d_inplane/tensor_rotation_run.sh BIN OUT [DEG]   (from miniapps/seas)
#
# BIN holds seas_fd_cycle_driver; DEG defaults to 30.  The deck (without its
# fault snapshots) is run three times:
#   plate_<run>  the deck as committed, 1 rank
#   rot_<run>    1 rank, mapview_rs_500m.msh's nodes and the deck turned by
#                DEG (rotate_case.py mesh2d, deck2d)
#   ctl_<run>    the plate deck on 2 ranks: the control, a change of
#                summation order only
# Runs: qd      the QD regime only, 300 steps (--max-steps 300; the first
#               switch comes after some 3000 steps)
#       event   QD up to the first switch, the whole first dynamic phase,
#               the run ends with it ([cycle.switch] max_dynamic_phases = 1);
#               both runs must log the switch at t > 0
# Each run's log must show the projected prestress of the deck's header on
# every fault DOF, turned or not: sigma_n = 25 MPa and |tau0| = 19.49 MPa
# (the driver's ranges, 6 digits).
# Compared with rotation_compare.py --z-deg DEG (global vectors as R times
# the unturned run's): first plate against ctl, which gives the size C that
# roundoff alone reaches, then plate against rot at the tolerance 100 C (at
# least 1e-12).  The control must stay below 1e-6 (qd) and 1e-4 (event):
# the nucleation amplifies roundoff, and through the whole event the control
# reaches 5.6e-6 (the switch times differ by 2 ms in 2.4e9 s; cloud runs
# 10-29: qd C = 4.4e-11, rot 3.1e-11; event C = 5.6e-6, rot 1.7e-6).
# Environment: MPIRUN (default "mpirun --oversubscribe"; add
# --allow-run-as-root as root), ONLY (space-separated run names).  About
# 30 s (qd) and 14 min (event) on the cloud container, 0.3 GB per rank.
#
# Scope (Chunhui Zhao, 2026-10-29).
set -u
export OMP_NUM_THREADS=1 LC_ALL=C
BIN=$(cd "${1:?usage: tensor_rotation_run.sh BIN OUT [DEG]}" && pwd)
OUT=${2:?usage: tensor_rotation_run.sh BIN OUT [DEG]}
DEG=${3:-30}
MPIRUN=${MPIRUN:-mpirun --oversubscribe}
ONLY=${ONLY:-}
PREFIX=mapview_rs_tensor_fd
[ -e "$OUT" ] && { echo "$OUT exists: not replacing it"; exit 2; }
HERE=tests/2d_inplane
D=$OUT/decks
mkdir -p "$D" || exit 2
python3 - config/fd_cycle/mapview_rs_tensor_fd.toml "$D" <<'PY' || exit 2
import re, sys
src, d = sys.argv[1], sys.argv[2]
s = open(src).read()
# the fault snapshots are not compared: leave them out (disk)
m = re.search(r"(?ms)^# Whole-field snapshots.*?^\[cycle\.snapshots\.fault\]\s*$.*?(?=^\[|\Z)", s)
assert m is not None, "no [cycle.snapshots.fault] section"
s = s[:m.start()] + s[m.end():]
assert "[cycle.snapshots" not in s
open(f"{d}/plate_qd.toml", "w").write(s)
a = "\nmax_dynamic_phases = 0"
assert s.count(a) == 1, a
open(f"{d}/plate_event.toml", "w").write(s.replace(a, "\nmax_dynamic_phases = 1", 1))
PY
python3 $HERE/rotate_case.py mesh2d mapview_rs/mesh/mapview_rs_500m.msh "$OUT/mapview_rs_rot.msh" "$DEG" \
   || exit 2
for r in qd event; do
   python3 $HERE/rotate_case.py deck2d "$D/plate_$r.toml" "$D/rot_$r.toml" "$OUT/mapview_rs_rot.msh" "$DEG" \
      || exit 2
done
NDOF=3   # p = 1 triangles: DOFs per element of the checkpoint's displacement

want() { [ -z "$ONLY" ] || [[ " $ONLY " == *" $1 "* ]]; }
fails=0
for r in qd event; do
   want $r || continue
   case $r in
      qd)    opt="--max-steps 300"; ctl_max=1e-6 ;;
      event) opt="";                ctl_max=1e-4 ;;
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
      # the projected prestress: the header's values on every DOF
      line=$(grep -m1 "^    sigma_n \[" "$OUT/${v}_$r/run.log")
      echo "   $line" | sed 's/; V_init.*//'
      if ! python3 - "$line" <<'PY'
import re, sys
v = [float(u) for pair in re.findall(r"\[([^,\]]+), ([^\]]+)\]", sys.argv[1]) for u in pair]
ok = len(v) == 6 and all(abs(x - 25.0) <= 1e-5 for x in v[0:2]) \
     and all(abs(x - 19.49) <= 1e-4 for x in v[2:4])
sys.exit(0 if ok else 1)
PY
      then echo "   FAILED: the prestress is not sigma_n = 25 MPa, |tau0| = 19.49 MPa on every DOF"
           fails=$((fails + 1)); fi
      if [ $r = event ]; then
         sw=$(awk '!/^#/ && $3 == "QD" && $4 == "dynamic" && $1 > 0 {print $1; exit}' \
              "$OUT/${v}_$r/${PREFIX}_switch.txt" 2>/dev/null)
         if [ -n "$sw" ]; then echo "   switch to dynamic at t = $sw s"
         else echo "   FAILED: no QD-to-dynamic switch at t > 0"; fails=$((fails + 1)); fi
      fi
   done
   echo "control: plate_$r against ctl_$r (1 rank against 2):"
   python3 $HERE/rotation_compare.py "$OUT/plate_$r" "$OUT/ctl_$r" --prefix $PREFIX --dim 2 \
      --skip-checkpoint --tol "$ctl_max" > "$OUT/control_$r.txt"
   crc=$?
   tail -1 "$OUT/control_$r.txt"
   C=$(awk '/^   largest /{print $2}' "$OUT/control_$r.txt")
   if [ $crc -ne 0 ] || [ -z "$C" ]; then
      echo "   FAILED: the control is not below $ctl_max (see $OUT/control_$r.txt)"
      fails=$((fails + 1)); continue
   fi
   tol=$(python3 -c "print('%.3g' % max(100 * float('$C'), 1e-12))")
   echo "compare plate_$r rot_$r (turned by $DEG degrees), tolerance 100 x control = $tol:"
   python3 $HERE/rotation_compare.py "$OUT/plate_$r" "$OUT/rot_$r" --z-deg "$DEG" --prefix $PREFIX \
      --elem-dofs "$NDOF" --tol "$tol" || fails=$((fails + 1))
done
echo "tensor_rotation_run: $fails failure(s)"
exit $fails
