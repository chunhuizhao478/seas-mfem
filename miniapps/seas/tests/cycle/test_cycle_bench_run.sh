#!/bin/bash
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# test_cycle_bench_run.sh — fully dynamic BP5 cycle: the p = 2 run's deck
# options in short driver runs, on decks made from
# config/fd_cycle/bp5_fd_p2_1800yr.toml at p = 1 (the cloud's memory) with
# a station row at every step and no snapshots:
#   restart   QD only (V_on = 100 m/s): 2K steps against K + restart + K,
#             every file byte for byte, the benchmark files (source,
#             off-fault, profiles) and the checkpoints included
#             (test_cycle_restart.sh)
#   velocity  the off-fault velocity of those rows (QDPhase::Velocity, the
#             exact quasi-static du/dt) against the central difference of
#             the displacement rows: max |v - (u(k+1) - u(k-1)) / (t(k+1) -
#             t(k-1))| / max |v| per component <= VTOL
#   matrix    QD only, K steps through the stage operators against the face
#             loops (compare_interp.py, MTOL)
#   lts       the first event dynamic from t = 0: LTS rate 2 against global
#             stepping on NP ranks, for DYN_T s of dynamic time
#             (compare_interp.py --scale quantity in their common range
#             after t = 0: each column against the size of its quantity
#             over the stations of its kind, since the dip components stay
#             at round-off and the far off-fault stations see no wave in
#             the first second; the QD row at t = 0 is the same in both,
#             and its quasi-static velocity would set the velocity scale).
#             The on-fault stations and the source rows within LTOL; the
#             off-fault stations within LTOL_OFF, since LTS moves them by
#             up to 3% near a wave front (global stepping at half the step
#             as the reference: 0.4%; p = 1, 2 ranks, 1 s)
# The Green's function is not run here (a dense G at 1000 m is 35 GB; the
# Vista check vista_p2_check.sbatch runs it).
#
#   tests/cycle/test_cycle_bench_run.sh WORK [K=10] [NP=2]   (from miniapps/seas)
# Environment: MPIRUN (default "mpirun --oversubscribe"), SOLVER
# (mumps_blr), VTOL (2e-2), MTOL (1e-8), LTOL (1e-3), LTOL_OFF (5e-2),
# DYN_T (1.0 s).
set -u
WORK=${1:?usage: test_cycle_bench_run.sh WORK [K] [NP]}
K=${2:-10}
NP=${3:-2}
MPIRUN=${MPIRUN:-mpirun --oversubscribe}
SOLVER=${SOLVER:-mumps_blr}
VTOL=${VTOL:-2e-2}
MTOL=${MTOL:-1e-8}
LTOL=${LTOL:-1e-3}
LTOL_OFF=${LTOL_OFF:-5e-2}
DYN_T=${DYN_T:-1.0}
export OMP_NUM_THREADS=1
BIN=./seas_fd_cycle_driver
[ -x "$BIN" ] || { echo "no $BIN (make -f cycle/cycle.mk seas_fd_cycle_driver)"; exit 2; }
[ ! -e "$WORK" ] || [ -f "$WORK/qd_loops.toml" ] \
   || { echo "$WORK exists and is not this test's folder: not replacing it"; exit 2; }
rm -rf "$WORK"
mkdir -p "$WORK"
fails=0

# ---- the decks ----
python3 - "$WORK" "$SOLVER" <<'EOF' || exit 2
import sys
work, solver = sys.argv[1], sys.argv[2]
s = open("config/fd_cycle/bp5_fd_p2_1800yr.toml").read()
def sub(s, a, b):   # the line that starts with a
    assert s.count("\n" + a) == 1, a
    return s.replace("\n" + a, "\n" + b, 1)
def lts_mode(s, m):  # [cycle.lts] mode, whatever the deck sets
    a = [l for l in ('mode         = "off"', 'mode         = "rate2"') if "\n" + l in s]
    assert len(a) == 1, a
    return sub(s, a[0], f'mode         = "{m}"')
s = s[:s.index("\n# Whole-field snapshots")] + "\n"
s = sub(s, "order = 2", "order = 1")
s = sub(s, 'solver        = "mumps"', f'solver        = "{solver}"')
s = sub(s, 'prefix       = "bp5_fd_p2_1800yr"', 'prefix       = "b"')
s = sub(s, "fast_V       = 1.0e-3", "fast_V       = 1.0e-12")
s = sub(s, "dt_fast      = 0.1", "dt_fast      = 1.0e-9")
s = sub(s, "dt_slow      = 3155760.0", "dt_slow      = 1.0e-9")
s = sub(s, "profile_dt_fast = 1.0", "profile_dt_fast = 1.0e-9")
s = sub(s, "[cycle.qd.green]", "[cycle.qd.green]\nenabled          = false")
dyn = s
s = sub(s, "V_on               = 1.0e-2", "V_on               = 1.0e2")
qd_matrix = s
qd_loops = sub(s, 'operators     = "matrix"', 'operators     = "loops"')
dyn = sub(dyn, 'operators     = "matrix"', 'operators     = "loops"')
lts = lts_mode(dyn, "rate2")
gts = lts_mode(dyn, "off")
for name, text in (("qd_loops", qd_loops), ("qd_matrix", qd_matrix), ("lts", lts),
                   ("gts", gts)):
    open(f"{work}/{name}.toml", "w").write(text)
EOF

run() {   # name, ranks, deck, driver options...
   local n=$1 np=$2 deck=$3; shift 3
   mkdir -p "$WORK/$n"
   local s=$(date +%s)
   $MPIRUN -np "$np" $BIN --deck "$deck" --output-dir "$WORK/$n" "$@" > "$WORK/$n/run.log" 2>&1
   local rc=$?
   echo "$n: exit $rc in $(( $(date +%s) - s )) s on $np rank(s)"
   [ $rc -eq 0 ] || { fails=$((fails + 1)); tail -20 "$WORK/$n/run.log"; }
}

echo "--- restart (QD only, $K + $K steps, 1 rank)"
MPIRUN="$MPIRUN" tests/cycle/test_cycle_restart.sh "$WORK/qd_loops.toml" "$WORK/restart" "$K" 1 \
   || fails=$((fails + 1))
grep -E "benchmark output|profile points" "$WORK/restart/full/run.log"

echo "--- velocity: the off-fault rows of the restart test's full run"
python3 - "$WORK/restart/full" "$VTOL" <<'EOF' || fails=$((fails + 1))
import glob, sys
run, tol = sys.argv[1], float(sys.argv[2])
files = sorted(glob.glob(run + "/b_blkst_*.txt"))
ok = len(files) == 9
for p in files:
    r = [list(map(float, l.split())) for l in open(p) if l.strip() and not l.startswith("#")]
    meas = []
    for c in range(3):
        err = vmax = 0.0
        for k in range(1, len(r) - 1):
            fd = (r[k + 1][1 + c] - r[k - 1][1 + c]) / (r[k + 1][0] - r[k - 1][0])
            err = max(err, abs(r[k][4 + c] - fd))
            vmax = max(vmax, abs(r[k][4 + c]))
        meas.append(err / vmax if vmax > 0 else err)
    good = len(r) >= 3 and all(m <= tol for m in meas)
    ok = ok and good
    print(f"  {'ok  ' if good else 'FAIL'} {p.split('/')[-1]}: {len(r)} rows, "
          f"max |v - du/dt| / max |v|: {' '.join(f'{m:.1e}' for m in meas)}")
print(("PASSED" if ok else "FAILED") + f": {len(files)} off-fault files (9 expected)")
sys.exit(0 if ok else 1)
EOF

echo "--- matrix: $K QD steps through the stage operators against the face loops (1 rank)"
run qd_matrix 1 "$WORK/qd_matrix.toml" --max-steps "$K"
grep -E "QD stages|stage operators" "$WORK/qd_matrix/run.log" | head -3
python3 tests/cycle/compare_interp.py "$WORK/qd_matrix" b "$WORK/restart/full" b \
   --tol "$MTOL" --t-max "$(awk '!/^#/{t=$1} END{print t}' "$WORK/qd_matrix/b_source.txt")" \
   || fails=$((fails + 1))

echo "--- lts: the first event dynamic, LTS rate 2 against global stepping, $NP ranks"
NG=$(python3 -c "print(int($DYN_T / 1.72e-3) + 1)")
run lts "$NP" "$WORK/lts.toml" --max-dyn-steps "$(( NG / 4 + 1 ))"
run gts "$NP" "$WORK/gts.toml" --max-dyn-steps "$NG"
grep -E "LTS:|end: " "$WORK/lts/run.log" "$WORK/gts/run.log"
python3 tests/cycle/compare_interp.py "$WORK/lts" b "$WORK/gts" b --tol "$LTOL" --scale quantity \
   --t-min 1e-9 --files "{fltst,source}" || fails=$((fails + 1))
python3 tests/cycle/compare_interp.py "$WORK/lts" b "$WORK/gts" b --tol "$LTOL_OFF" --scale quantity \
   --t-min 1e-9 --files blkst || fails=$((fails + 1))

echo "test_cycle_bench_run: $([ $fails -eq 0 ] && echo PASSED || echo "FAILED ($fails)")"
[ $fails -eq 0 ]
