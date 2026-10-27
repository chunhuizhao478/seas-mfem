#!/bin/bash
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# test_cycle_snapshots.sh — the fault and surface snapshots of
# seas_fd_cycle_driver ([cycle.snapshots], cycle/snapshots.hpp) in six short
# runs and a seventh on the BP5 mesh, from config/fd_cycle/bp5_fd_first_event.toml with
# t_final 0.5 s, a station row at every step, a checkpoint every 2 QD steps
# and snapshots every 0.05 s in QD and 0.01 s in a dynamic phase:
#   A  QD from t = 0 (max V 0.015 m/s), dynamic once max V passes V_on =
#      0.016 m/s, ended after 30 steps (--end-phase-after), back to QD, where
#      the guard stops the run (exit 7) after the switch-back checkpoint
#   B  --restart from A's checkpoint with V_on = 100 (QD only) to 0.5 s
#   C  B's folder with A's checkpoint put back, the same restart: the files
#      are trimmed to A's counts and B is redone
#   F  QD only from t = 0 to 0.5 s, straight through
#   D  F with --wall-limit WALL (exit 8), which must stop it after its
#      first QD step and before its last
#   E  --restart from D's checkpoint, in place, to 0.5 s: F again
#   G  as A with V_on = 0.014 m/s: dynamic from t = 0 (the "start"
#      snapshot, as the 1500 yr deck), 30 steps, the guard (exit 7)
# Checks (tests/cycle/check_snapshots.py): layout, ParaView's reader,
# the fault values against the station rows, the schedule and the switch
# and end labels for every run; B's first snapshots are A's, bit for bit;
# C matches B and E matches F to TOL (MUMPS on more than one rank is not bit
# for bit reproducible); D stopped part way; the exit codes.  The read-back
# checks need python3 with numpy and h5py (and vtk for ParaView's reader);
# without them the script says so, checks the exit codes only and prints
# the snapshot counts with h5ls.
#
#   tests/cycle/test_cycle_snapshots.sh WORK [NP]      (from miniapps/seas)
# Environment: MPIRUN (default "mpirun --oversubscribe"), SOLVER
# (mumps_blr, the cloud's memory; or mumps), WALL (D's limit in s; default
# 70% of F's wall time: the setup takes more than 150 s of F's 377 s on 2
# cloud ranks), TOL (default 1e-8).
set -u
WORK=${1:?usage: test_cycle_snapshots.sh WORK [NP]}
NP=${2:-2}
MPIRUN=${MPIRUN:-mpirun --oversubscribe}
SOLVER=${SOLVER:-mumps_blr}
WALL=${WALL:-}
TOL=${TOL:-1e-8}
export OMP_NUM_THREADS=1
BIN=./seas_fd_cycle_driver
CHECK="python3 tests/cycle/check_snapshots.py"
READ=1
NOVTK=
python3 -c "import numpy, h5py" 2>/dev/null || {
   READ=0
   echo "NOTE: $(python3 --version 2>&1) has no numpy or h5py: exit codes only, files not read back"; }
[ $READ -eq 0 ] || python3 -c "import vtk" 2>/dev/null || {
   NOVTK=--no-vtk
   echo "NOTE: no python vtk: ParaView's reader (vtkHDFReader) is not checked"; }
[ -x "$BIN" ] || { echo "no $BIN (make -f cycle/cycle.mk seas_fd_cycle_driver)"; exit 2; }
[ ! -e "$WORK" ] || [ -f "$WORK/snapA.toml" ] \
   || { echo "$WORK exists and is not this test's folder: not replacing it"; exit 2; }
rm -rf "$WORK"
mkdir -p "$WORK"

# ---- the decks ----
python3 - "$WORK" "$SOLVER" <<'EOF'
import sys
work, solver = sys.argv[1], sys.argv[2]
s = open("config/fd_cycle/bp5_fd_first_event.toml").read()
def sub(s, a, b):   # the line that starts with a
    assert s.count("\n" + a) == 1, a
    return s.replace("\n" + a, "\n" + b, 1)
s = sub(s, "tfinal     = 5.680368e10", "tfinal     = 0.5")
s = sub(s, "V_on               = 1.0e-2", "V_on               = 1.6e-2")
s = sub(s, "max_dynamic_phases = 1", "max_dynamic_phases = 0")
s = sub(s, 'solver        = "mumps"', f'solver        = "{solver}"')
s = sub(s, 'prefix       = "bp5_fd_first_event"', 'prefix       = "snap"')
s = sub(s, "fast_V       = 1.0e-3", "fast_V       = 1.0e-12")
s = sub(s, "dt_fast      = 0.1", "dt_fast      = 1.0e-9")
s = sub(s, "dt_slow      = 3155760.0", "dt_slow      = 1.0e-9")
s = sub(s, "every_qd_steps = 0", "every_qd_steps = 2")
s += """
[cycle.snapshots.fault]
dt_qd  = 0.05
dt_dyn = 0.01

[cycle.snapshots.surface]
dt_qd        = 0.05
dt_dyn       = 0.01
attrs        = [1]
plane_normal = [0.0, 0.0, 1.0]
plane_offset = 0.0

[cycle.snapshots.volume]
enabled = false
"""
open(f"{work}/snapA.toml", "w").write(s)
open(f"{work}/snapB.toml", "w").write(sub(s, "V_on               = 1.6e-2",
                                          "V_on               = 100.0"))
open(f"{work}/snapG.toml", "w").write(sub(s, "V_on               = 1.6e-2",
                                          "V_on               = 1.4e-2"))
EOF
[ $? -eq 0 ] || { echo "FAILED: the decks"; exit 1; }

fails=0
run() {   # name, expected exit, log, driver arguments...
   local name=$1 want=$2 log=$3; shift 3
   local s=$(date +%s)
   $MPIRUN -np "$NP" $BIN "$@" > "$log" 2>&1
   local rc=$?
   echo "$name: exit $rc (expected $want), $(( $(date +%s) - s )) s"
   [ "$rc" = "$want" ] || { echo "FAILED: $name"; tail -20 "$log"; fails=$((fails + 1)); }
}
check() {   # dir, deck, extra arguments
   local dir=$1 deck=$2; shift 2
   if [ $READ -eq 0 ]; then
      for k in fault surface; do
         h5ls "$dir/snap_${k}_snapshots.vtkhdf/VTKHDF/Steps/Values" 2>&1 | sed "s|^|$k: |"
      done
      return
   fi
   local dt
   dt=$(grep -m1 -oE "dt = [0-9.e+-]+ s" "$dir"/run*.log | head -1 | awk '{print $3}')
   $CHECK --dir "$dir" --prefix snap --deck "$deck" ${dt:+--dt-dyn-step $dt} $NOVTK "$@" \
      > "$dir/check.log" 2>&1
   local rc=$?
   grep -E "layout|vtkHDFReader|station rows matched| vs |FAIL|ALL OK|FAILURE" "$dir/check.log"
   [ $rc -eq 0 ] || fails=$((fails + 1))
}

W=$WORK
mkdir -p $W/A
run A 7 $W/A/run.log --deck $W/snapA.toml --output-dir $W/A --end-phase-after 30 --log-every 10
cp -a $W/A $W/A_end
cp -a $W/A $W/B
run B 0 $W/B/run_B.log --deck $W/snapB.toml --output-dir $W/B --restart $W/B/snap --log-every 10
cp -a $W/B $W/C
cp $W/A_end/snap_checkpoint_r*.txt $W/C/
run C 0 $W/C/run_C.log --deck $W/snapB.toml --output-dir $W/C --restart $W/C/snap --log-every 10
mkdir -p $W/F
s0=$(date +%s)
run F 0 $W/F/run.log --deck $W/snapB.toml --output-dir $W/F --log-every 10
WALL=${WALL:-$(( ($(date +%s) - s0) * 7 / 10 ))}
echo "D's wall limit: $WALL s"
mkdir -p $W/D
run D 8 $W/D/run.log --deck $W/snapB.toml --output-dir $W/D --wall-limit "$WALL" --log-every 10
cp -a $W/D $W/D_end
mv $W/D $W/E
run E 0 $W/E/run_E.log --deck $W/snapB.toml --output-dir $W/E --restart $W/E/snap \
   --wall-limit 100000 --log-every 10
mkdir -p $W/G
run G 7 $W/G/run.log --deck $W/snapG.toml --output-dir $W/G --end-phase-after 30 --log-every 10

echo "--- A"; check $W/A $W/snapA.toml
echo "--- B"; check $W/B $W/snapB.toml --complete
echo "--- C vs B"; check $W/C $W/snapB.toml --against $W/B --tol "$TOL" --no-vtk --complete
echo "--- D (as it stopped)"; check $W/D_end $W/snapB.toml --no-vtk
echo "--- E vs F"; check $W/E $W/snapB.toml --against $W/F --tol "$TOL" --complete
echo "--- F"; check $W/F $W/snapB.toml --no-vtk --complete
echo "--- G"; check $W/G $W/snapG.toml --no-vtk
[ $READ -eq 0 ] || python3 - $W <<'EOF' || fails=$((fails + 1))
import sys, h5py, numpy as np
f = h5py.File(f"{sys.argv[1]}/G/snap_fault_snapshots.vtkhdf", "r")
r, g, t = f["fd_cycle/reason"][:], f["fd_cycle/regime"][:], f["VTKHDF/Steps/Values"][:]
ok = len(t) > 2 and t[0] == 0 and r[0] == 1 and g[0] == 0 and np.all(g[1:] == 1) and r[-1] == 3
print(f"G: {len(t)} snapshots, 'start' at t = 0, then the dynamic phase to 'end of dynamic': "
      f"{'OK' if ok else 'NO'}")
sys.exit(0 if ok else 1)
EOF
echo "--- B continues A"
[ $READ -eq 0 ] || python3 - $W <<'EOF' || fails=$((fails + 1))
import sys, h5py, numpy as np
w = sys.argv[1]
bad = 0
for kind in ("fault", "surface"):
    a = h5py.File(f"{w}/A_end/snap_{kind}_snapshots.vtkhdf", "r")
    b = h5py.File(f"{w}/B/snap_{kind}_snapshots.vtkhdf", "r")
    na = len(a["VTKHDF/Steps/Values"]); nb = len(b["VTKHDF/Steps/Values"])
    n = int(a["VTKHDF/NumberOfPoints"][0])
    same = all(np.array_equal(a["VTKHDF/PointData"][f][:na * n], b["VTKHDF/PointData"][f][:na * n])
               for f in a["VTKHDF/PointData"])
    same = same and np.array_equal(a["VTKHDF/Steps/Values"][:], b["VTKHDF/Steps/Values"][:na])
    later = b["VTKHDF/Steps/Values"][na:]
    ok = same and nb > na and np.all(later > a["VTKHDF/Steps/Values"][-1]) \
         and np.all(b["fd_cycle/regime"][na:] == 0)
    print(f"{kind}: A has {na} snapshots, B {nb}; B's first {na} are A's bit for bit: {same}; "
          f"later ones QD and after A's last: {ok}")
    bad += 0 if ok else 1
sys.exit(1 if bad else 0)
EOF
echo "--- D stopped at the wall limit and E continued it"
grep -E "wall limit .* reached|status:" $W/D_end/run.log
grep -E "restarted from|status:" $W/E/run_E.log
[ $READ -eq 0 ] || python3 - $W <<'EOF' || fails=$((fails + 1))
import sys, h5py
w = sys.argv[1]
nd = len(h5py.File(f"{w}/D_end/snap_fault_snapshots.vtkhdf", "r")["VTKHDF/Steps/Values"])
nf = len(h5py.File(f"{w}/F/snap_fault_snapshots.vtkhdf", "r")["VTKHDF/Steps/Values"])
ok = 1 < nd < nf
print(f"D stopped with {nd} of F's {nf} fault snapshots: {'OK' if ok else 'not part way (set WALL)'}")
sys.exit(0 if ok else 1)
EOF
ONLY=; [ $READ -eq 1 ] || ONLY=" (exit codes only: no numpy or h5py)"
if [ $fails -eq 0 ]; then echo "test_cycle_snapshots: PASSED$ONLY"; else echo "test_cycle_snapshots: FAILED ($fails)$ONLY"; fi
exit $fails
