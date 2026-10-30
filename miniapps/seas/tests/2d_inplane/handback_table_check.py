#!/usr/bin/env python3
# tests/2d_inplane/handback_table_check.py — make test-handback-table: the
# decks and the checks of [cycle.output].handback_table (Chunhui Zhao,
# 2026-10-30: the UCERF 2-D run stopped at the switch back, max V dynamic
# 7.8e-5 m/s but QD balance 0.253 m/s; "this sounds like a bug... we need to
# debug"; the card "Hand-back report first, then rest").
#
#   python3 handback_table_check.py decks DECK TMP
#       From the map-view tensor deck: TMP/on.toml (max_dynamic_phases 2, the
#       event started at t = 0 by a V_init patch of 0.02 m/s, |x| <= 3 km,
#       linear to 1e-9 m/s at 11 km, handback_table = true), TMP/off.toml (the
#       same without the table) and TMP/bad.toml (handback_table = 1).
#   python3 handback_table_check.py check TMP PREFIX
#       After seas_fd_cycle_driver ran on.toml into TMP/on (log TMP/on.log)
#       and TMP/on2 (2 ranks, TMP/on2.log), off.toml into TMP/off, each with
#       --end-phase-after N.  Checks:
#        1. TMP/on and TMP/on2 hold PREFIX_handback_phase1.txt, TMP/off none;
#        2. one row of 15 numbers per row of the fault table PREFIX_fault.txt,
#           row by row the same position and the same slip (both are the
#           dynamic regime's slip at its switch, in the same DOF order; a
#           position shows up once per face at it), on 1 and 2 ranks;
#        3. the largest |V| of the dynamic and of the QD columns are the
#           log's "max V dynamic" and "QD balance" (printed with 6 decimals);
#        4. the log's "hand-back table" line: the row count and the largest
#           |tau_qd - tau_dyn|, |sigma_n_qd - sigma_n_dyn| and V_qd of the
#           table;
#        5. every other file of TMP/off is byte for byte the one of TMP/on
#           (the table changes no output); a .vtkhdf file: every dataset and
#           attribute (h5py), since HDF5 stores object times in the file;
#        6. 2 ranks: the rows, sorted by position and slip, are those of 1
#           rank (MUMPS on 2 ranks: about 1e-10 relative).
# Exit 0 when every check passes.

import math
import os
import re
import sys

import h5py
import numpy as np

NC = 15
failed = 0


def check(ok, what):
    global failed
    print(("  PASSED: " if ok else "  FAILED: ") + what)
    if not ok:
        failed += 1


def replace_once(text, old, new, what):
    n = text.count(old)
    if n != 1:
        sys.exit("handback_table_check.py: the deck's %s text occurs %d times (want 1)" % (what, n))
    return text.replace(old, new)


def decks(deck, tmp):
    with open(deck) as f:
        base = f.read()
    base = replace_once(base, "max_dynamic_phases = 0", "max_dynamic_phases = 2",
                        "max_dynamic_phases")
    base = replace_once(base,
                        "[cycle.fields.V_init]\ncenter  = [0.0, 0.0]\naxis1   = [1.0, 0.0]\n"
                        "half1   = 0.0\ninside  = 1.0e-9\n",
                        "[cycle.fields.V_init]\ncenter  = [0.0, 0.0]\naxis1   = [1.0, 0.0]\n"
                        "half1   = 3.0e3\ntaper   = 8.0e3\ninside  = 2.0e-2\n", "V_init")
    pre = re.search(r'^prefix *= *"[^"]*"\n', base, re.M)
    if not pre:
        sys.exit("handback_table_check.py: no [cycle.output] prefix line in the deck")
    line = pre.group(0)
    os.makedirs(tmp, exist_ok=True)
    for name, extra in (("off", ""), ("on", "handback_table  = true\n"),
                        ("bad", "handback_table  = 1\n")):
        with open(os.path.join(tmp, name + ".toml"), "w") as f:
            f.write(replace_once(base, line, line + extra, "prefix"))


def rows(path):
    out = []
    with open(path) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            out.append([float(v) for v in ln.split()])
    return out


def log_numbers(path):
    with open(path) as f:
        log = f.read()
    m = re.search(r"max V dynamic (\S+) m/s, QD balance (\S+) m/s", log)
    h = re.search(r"hand-back table (\S+) \((\d+) DOFs\): largest \|tau_qd - tau_dyn\| (\S+) "
                  r"MPa at .*?; largest \|sigma_n_qd - sigma_n_dyn\| (\S+) MPa at .*?; "
                  r"largest V_qd (\S+) m/s at", log)
    return m, h


def same_hdf5(a, b):
    # Every dataset (shape and values) and every attribute of both files.
    with h5py.File(a, "r") as fa, h5py.File(b, "r") as fb:
        bad = []

        def visit(name, obj):
            if name not in fb:
                bad.append(name)
                return
            other = fb[name]
            if isinstance(obj, h5py.Dataset):
                if not isinstance(other, h5py.Dataset) or not np.array_equal(obj[()], other[()]):
                    bad.append(name)
            if set(obj.attrs.keys()) != set(other.attrs.keys()) or not all(
                    np.array_equal(np.asarray(obj.attrs[k]), np.asarray(other.attrs[k]))
                    for k in obj.attrs.keys()):
                bad.append(name + " (attributes)")

        fa.visititems(visit)
        names = []
        fb.visit(names.append)
        n_a = []
        fa.visit(n_a.append)
        return not bad and sorted(names) == sorted(n_a)


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def run_checks(tmp, prefix):
    name = prefix + "_handback_phase1.txt"
    for d in ("on", "on2"):
        check(os.path.isfile(os.path.join(tmp, d, name)), "%s/%s written" % (d, name))
    stray = [f for f in os.listdir(os.path.join(tmp, "off")) if "_handback_" in f]
    check(not stray, "off/ has no hand-back table (%s)" % (stray or "none"))
    if failed:
        return

    tab = rows(os.path.join(tmp, "on", name))
    for d in ("on", "on2"):
        t = tab if d == "on" else rows(os.path.join(tmp, d, name))
        fault = rows(os.path.join(tmp, d, prefix + "_fault.txt"))
        check(all(len(r) == NC for r in t) and len(t) == len(fault) and len(t) > 0,
              "%s: %d rows of %d numbers, one per fault-table row (%d)"
              % (d, len(t), NC, len(fault)))
        if len(t) != len(fault):
            continue
        # Fault table columns: x y z t_rupture slip_strike slip_dip peak_V psi.
        same_x = all(r[0:3] == q[0:3] for r, q in zip(t, fault))
        dslip = max(max(abs(r[3] - q[4]), abs(r[4] - q[5])) for r, q in zip(t, fault))
        check(same_x and dslip == 0.0,
              "%s: row by row the fault table's position and final slip (largest slip "
              "difference %.3e m)" % (d, dslip))

    vd = max(math.hypot(r[5], r[7]) for r in tab)
    vq = max(math.hypot(r[6], r[8]) for r in tab)
    dtau = max(math.hypot(r[10] - r[9], r[12] - r[11]) for r in tab)
    dsig = max(abs(r[14] - r[13]) for r in tab)
    m, h = log_numbers(os.path.join(tmp, "on.log"))
    check(m is not None, "the log's back-to-QD line")
    if m:
        check(abs(vd - float(m.group(1))) <= 1e-6,
              "largest V_dyn %.6e = the log's max V dynamic %s" % (vd, m.group(1)))
        check(abs(vq - float(m.group(2))) <= 1e-6,
              "largest V_qd %.6e = the log's QD balance %s" % (vq, m.group(2)))
    check(h is not None, "the log's hand-back table line")
    if h:
        check(int(h.group(2)) == len(tab), "the log's row count %s" % h.group(2))
        check(rel(dtau, float(h.group(3))) < 2e-6 and rel(dsig, float(h.group(4))) < 2e-6
              and rel(vq, float(h.group(5))) < 2e-6,
              "the log's largest |dtau| %s, |dsigma_n| %s MPa, V_qd %s m/s = the table's "
              "%.6e, %.6e, %.6e" % (h.group(3), h.group(4), h.group(5), dtau, dsig, vq))

    files = sorted(os.listdir(os.path.join(tmp, "off")))
    diff = []
    for fn in files:
        a = os.path.join(tmp, "off", fn)
        b = os.path.join(tmp, "on", fn)
        if not os.path.isfile(b):
            diff.append(fn + " (missing)")
            continue
        if fn.endswith(".vtkhdf"):
            if not same_hdf5(a, b):
                diff.append(fn)
            continue
        with open(a, "rb") as fa, open(b, "rb") as fb:
            if fa.read() != fb.read():
                diff.append(fn)
    check(not diff and len(files) > 0,
          "the other %d files byte for byte with and without the table%s"
          % (len(files), "" if not diff else ": differ " + ", ".join(diff)))

    t2 = rows(os.path.join(tmp, "on2", name))
    a = sorted(tab, key=lambda r: tuple(r[0:5]))
    b = sorted(t2, key=lambda r: tuple(r[0:5]))
    ok = len(a) == len(b) and all(ra[0:3] == rb[0:3] for ra, rb in zip(a, b))
    check(ok, "2 ranks: the same %d positions" % len(b))
    if ok:
        worst = 0.0
        for ra, rb in zip(a, b):
            for k in range(3, NC):
                scale = max(abs(ra[k]), abs(rb[k]), 1e-12 if k in (5, 6, 7, 8) else 1e-6)
                worst = max(worst, abs(ra[k] - rb[k]) / scale)
        check(worst < 1e-6, "2 ranks: every value within %.3e relative of 1 rank" % worst)


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "decks":
        decks(sys.argv[2], sys.argv[3])
        sys.exit(0)
    if len(sys.argv) == 4 and sys.argv[1] == "check":
        run_checks(sys.argv[2], sys.argv[3])
        sys.exit(1 if failed else 0)
    sys.exit(__doc__ or "usage: handback_table_check.py decks DECK TMP | check TMP PREFIX")
