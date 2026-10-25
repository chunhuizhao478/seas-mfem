#!/bin/bash
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# tests/2d_inplane/parity_3d.sh — 2D in-plane design rev 3: the 3D parity
# set that every phase of the branch must leave bit for bit unchanged.
#
#   tests/2d_inplane/parity_3d.sh run BIN OUT      run the set with the binaries in BIN
#   tests/2d_inplane/parity_3d.sh compare A B      compare two OUT folders file by file
#   tests/2d_inplane/parity_3d.sh decks OUT        only write the run decks, OUT/decks
#                                                  (tests/2d_inplane/rotation_run.sh)
#
# Run from miniapps/seas.  BIN holds seas_nitsche_dyn_driver,
# seas_fd_cycle_driver, seas_bp5_full and the unit-test binaries listed below
# (missing ones are reported and skipped).  Environment: MESHES (folder with
# tpv2053d_1000m.msh and tpv104_1000m.msh, made with
#   gmsh -3 tpv205/mesh/tpv2053d_1000m.geo -format msh22 -o $MESHES/tpv2053d_1000m.msh
#   gmsh -3 tpv104/mesh/tpv104_1000m.geo   -format msh22 -o $MESHES/tpv104_1000m.msh),
# MPIRUN (default "mpirun --oversubscribe"), ONLY (space-separated run names,
# default all), PV_FAULT_OFF=1 (Nitsche decks with paraview_fault = "off",
# for an MFEM built without HDF5; the station files are unchanged).  One run
# at a time.
#
# The runs (each writes OUT/<name>/ with its files and run.log):
#   tpv205_np2     Nitsche driver, TPV205 deck, 1000 m p = 1, 1.0 s, 2 ranks
#   tpv104_np2     Nitsche driver, TPV104 deck, 1000 m p = 1, 1.0 s, 2 ranks
#   cycle_qd_np1   cycle driver, QD only, 10 steps, 1 rank      (decks made from
#   cycle_qd_np2   cycle driver, QD only, 10 steps, 2 ranks      config/fd_cycle/
#   cycle_dyn_np1  cycle driver, the first event dynamic from    bp5_fd_p2_1800yr.toml
#   cycle_dyn_np2  t = 0, 40 dynamic steps, 1 and 2 ranks        at p = 1 with
#                                                                mumps_blr, as
#                                                                tests/cycle/test_cycle_bench_run.sh;
#                                                                full-rank MUMPS asks
#                                                                for 18 GB at 1 rank)
#   bp5_quick_np8  legacy QD driver seas_bp5_full, the Makefile's quick-check:
#                  inline mesh, 20 steps, 8 ranks, against the stored golden
#                  bp5/benchmark_data/golden_inline_8r_20step at tolerance 0
#   unit_<test>    the unit tests' output (exit code in run.log's last line)
#
# compare: every regular file of A except run.log must exist in B with the same
# bytes, after the commit stamp in headers ("commit <hash>[-dirty]") is set
# aside, and in checkpoint files the byte sizes they record for the other
# output files (which change with the stamp's length); such files are counted
# as "same but the commit stamp".  A file that still differs is listed with
# its largest numeric difference, scaled by the column's largest |value|.
# run.log files are compared after dropping lines with wall times, dates,
# commits and paths, and are reported, not failed.  Exit 0 = identical.
#
# QD runs on more than one rank are not reproducible bit for bit from one run
# to the next, at 3926e9c8 as well: MUMPS's parallel solve differs in the
# last bits from one call to the next (tests/cycle/test_cycle_p3.cpp).
# Measured 10-24 in the cloud container with the same binaries twice:
# cycle_qd_np2 up to 3.7e-9 and bp5_quick_np8 up to 9.3e-8 scaled.  For those
# runs, also run the base set a second time (base2) and check that
# compare base B reports differences of the size compare base base2 does.
set -u
export OMP_NUM_THREADS=1 LC_ALL=C

mode=${1:?usage: parity_3d.sh run BIN OUT | compare A B}
if [ "$mode" = compare ]; then
   A=${2:?}; B=${3:?}
   python3 - "$A" "$B" <<'PY'
import os, re, sys
A, B = sys.argv[1], sys.argv[2]
bad = same = stamp = 0
noise = re.compile(r"(wall|time [0-9.]+ s|[0-9.]+ s/step|s/step|seconds|elapsed|date|commit|"
                   r"/tmp/|/Users/|/home/|MB|MiB|GiB|rss|RSS|ms/step|s\)|timing|Timer|took)", re.I)
commit = re.compile(rb"commit [0-9a-f]{7,40}(-dirty)?")
size_line = re.compile(rb"^(\S+\.txt) [0-9]+$", re.M)
def normalized(path, name):
    b = commit.sub(b"commit <hash>", open(path, "rb").read())
    if "checkpoint" in name:
        b = size_line.sub(rb"\1 <size>", b)
    return b
def numbers(b):
    rows = []
    for l in b.decode(errors="replace").splitlines():
        if l.startswith("#"):
            continue
        r = []
        for tok in l.split():
            try: r.append(float(tok))
            except ValueError: r.append(None)
        rows.append(r)
    return rows
def worst_diff(ba, bb):
    ta, tb = numbers(ba), numbers(bb)
    if len(ta) != len(tb):
        return f"{len(ta)} vs {len(tb)} lines"
    w, col = 0.0, -1
    scale = {}
    for r in ta:
        for c, v in enumerate(r):
            if v is not None: scale[c] = max(scale.get(c, 0.0), abs(v))
    for ra, rb in zip(ta, tb):
        for c, (x, y) in enumerate(zip(ra, rb)):
            if x is not None and y is not None and x != y:
                d = abs(x - y) / (scale[c] or 1.0)
                if d > w: w, col = d, c
    return f"largest scaled difference {w:.3e} (column {col})"
for root, _, files in os.walk(A):
    for f in sorted(files):
        pa = os.path.join(root, f)
        rel = os.path.relpath(pa, A)
        pb = os.path.join(B, rel)
        if not os.path.exists(pb):
            print(f"MISSING  {rel}"); bad += 1; continue
        if f == "run.log":
            la = [l for l in open(pa, errors="replace") if not noise.search(l)]
            lb = [l for l in open(pb, errors="replace") if not noise.search(l)]
            if la != lb:
                d = sum(1 for x, y in zip(la, lb) if x != y) + abs(len(la) - len(lb))
                print(f"log-diff {rel}: {d} lines differ (reported, not failed)")
            continue
        ra, rb = open(pa, "rb").read(), open(pb, "rb").read()
        if ra == rb:
            same += 1; continue
        na, nb = normalized(pa, f), normalized(pb, f)
        if na == nb:
            stamp += 1; continue
        print(f"DIFFER   {rel}: {worst_diff(na, nb)}"); bad += 1
for root, _, files in os.walk(B):
    for f in files:
        rel = os.path.relpath(os.path.join(root, f), B)
        if not os.path.exists(os.path.join(A, rel)):
            print(f"EXTRA    {rel}"); bad += 1
print(f"{same} files identical, {stamp} the same but the commit stamp, "
      f"{bad} differ or are missing")
sys.exit(1 if bad else 0)
PY
   exit $?
fi

[ "$mode" = run ] || [ "$mode" = decks ] || { echo "unknown mode $mode"; exit 2; }
if [ "$mode" = decks ]; then
   OUT=${2:?}
else
   BIN=$(cd "${2:?}" && pwd); OUT=${3:?}
   MESHES=${MESHES:?set MESHES to the folder with the 1000 m TPV205/TPV104 meshes}
fi
MPIRUN=${MPIRUN:-mpirun --oversubscribe}
ONLY=${ONLY:-}
[ -e "$OUT" ] && { echo "$OUT exists: not replacing it"; exit 2; }
mkdir -p "$OUT/decks"
fails=0

python3 - "$OUT/decks" "${PV_FAULT_OFF:-0}" <<'PY' || exit 2
import sys
work, pv_off = sys.argv[1], sys.argv[2] == "1"
for name, src in (("tpv205", "tpv205/configs/tpv205_nitsche.toml"),
                  ("tpv104", "tpv104/configs/tpv104_nitsche.toml")):
    t = open(src).read()
    if pv_off:
        a = 'paraview_fault         = "hdf5"'
        assert t.count("\n" + a) == 1, (src, a)
        t = t.replace("\n" + a, '\nparaview_fault         = "off"', 1)
    open(f"{work}/{name}_nitsche.toml", "w").write(t)
s = open("config/fd_cycle/bp5_fd_p2_1800yr.toml").read()
def sub(s, a, b):
    assert s.count("\n" + a) == 1, a
    return s.replace("\n" + a, "\n" + b, 1)
def lts_mode(s, m):
    a = [l for l in ('mode         = "off"', 'mode         = "rate2"') if "\n" + l in s]
    assert len(a) == 1, a
    return sub(s, a[0], f'mode         = "{m}"')
s = s[:s.index("\n# Whole-field snapshots")] + "\n"
s = sub(s, "order = 2", "order = 1")
s = sub(s, 'solver        = "mumps"', 'solver        = "mumps_blr"')
s = sub(s, 'prefix       = "bp5_fd_p2_1800yr"', 'prefix       = "b"')
s = sub(s, "fast_V       = 1.0e-3", "fast_V       = 1.0e-12")
s = sub(s, "dt_fast      = 0.1", "dt_fast      = 1.0e-9")
s = sub(s, "dt_slow      = 3155760.0", "dt_slow      = 1.0e-9")
s = sub(s, "profile_dt_fast = 1.0", "profile_dt_fast = 1.0e-9")
s = sub(s, "[cycle.qd.green]", "[cycle.qd.green]\nenabled          = false")
s = sub(s, 'operators     = "matrix"', 'operators     = "loops"')
dyn = lts_mode(s, "off")
qd = sub(s, "V_on               = 1.0e-2", "V_on               = 1.0e2")
open(f"{work}/cycle_qd.toml", "w").write(qd)
open(f"{work}/cycle_dyn.toml", "w").write(dyn)
PY
[ "$mode" = decks ] && exit 0

want() { [ -z "$ONLY" ] || [[ " $ONLY " == *" $1 "* ]]; }
run() {   # name, ranks, binary, options...
   local n=$1 np=$2 b=$3; shift 3
   want "$n" || return 0
   [ -x "$BIN/$b" ] || { echo "$n: no $BIN/$b, skipped"; return 0; }
   mkdir -p "$OUT/$n"
   local s=$(date +%s)
   $MPIRUN -np "$np" "$BIN/$b" "$@" < /dev/null > "$OUT/$n/run.log" 2>&1
   local rc=$?
   echo "exit $rc" >> "$OUT/$n/run.log"
   echo "$n: exit $rc in $(( $(date +%s) - s )) s on $np rank(s)"
   [ $rc -eq 0 ] || { fails=$((fails + 1)); tail -15 "$OUT/$n/run.log"; }
}

run tpv205_np2 2 seas_nitsche_dyn_driver --config "$OUT/decks/tpv205_nitsche.toml" \
   --mesh "$MESHES/tpv2053d_1000m.msh" --tfinal 1.0 --output-dir "$OUT/tpv205_np2"
run tpv104_np2 2 seas_nitsche_dyn_driver --config "$OUT/decks/tpv104_nitsche.toml" \
   --mesh "$MESHES/tpv104_1000m.msh" --tfinal 1.0 --output-dir "$OUT/tpv104_np2"
for np in 1 2; do
   run cycle_qd_np$np $np seas_fd_cycle_driver --deck "$OUT/decks/cycle_qd.toml" \
      --output-dir "$OUT/cycle_qd_np$np" --max-steps 10
   run cycle_dyn_np$np $np seas_fd_cycle_driver --deck "$OUT/decks/cycle_dyn.toml" \
      --output-dir "$OUT/cycle_dyn_np$np" --max-dyn-steps 40
done
run bp5_quick_np8 8 seas_bp5_full --inline-mesh --max-steps 20 \
   --output-dir "$OUT/bp5_quick_np8" --output-prefix bp5_quick --tandem-time-stepping \
   --regression-tolerance 0 --ref-dir bp5/benchmark_data/golden_inline_8r_20step

# unit tests: name ranks [args]
while read -r t np args; do
   [ -n "$t" ] || continue
   run "unit_${t#seas_test_}_np$np" "$np" "$t" $args
done <<'LIST'
seas_test_elasticity_operator 1
seas_test_elasticity_br2 1
seas_test_fault_basis 1
seas_test_fault_dof_basis 1
seas_test_fault_basis_dip_strike_symmetry 1
seas_test_planar_tpv_basis 1
seas_test_qd_stage_operators 1
seas_test_qd_stage_operators 2
seas_test_nitsche_gates 1
seas_test_nitsche_pf 1
seas_test_nitsche_pf 2
seas_test_nitsche_p1 1
seas_test_nitsche_p2 1
seas_test_nitsche_p2 2
seas_test_parallel_elasticity 2
LIST

echo "parity_3d run: $fails failed run(s)"
exit $fails
