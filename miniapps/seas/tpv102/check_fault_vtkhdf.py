#!/usr/bin/env python3
"""Check the Nitsche driver's fault ParaView output (<output_dir>/fault.vtkhdf).

usage: check_fault_vtkhdf.py OUTDIR [OTHER_OUTDIR] [--prefix tpv102] [--frames K]

OUTDIR must hold the run's fault.vtkhdf, its CSV fault dumps
<prefix>_fault_r*.csv written at the SAME cadence (--dump-every = the
ParaView interval in steps) and <prefix>_station_conventions.txt.  Checks:
  * exactly one .vtkhdf; the mesh is stored once (static geometry);
  * number of steps == number of CSV dump frames, same times;
  * cells == 4 x fault faces; the 4 sub-areas sum to the face area;
  * no inverted sub-triangle (normal . nhat > 0);
  * every cell value == the mean of its 3 corner values rebuilt from the
    CSV DOF rows (a corner that is not a DOF, the p = 1 edge midpoints, is
    the average of the two vertex DOFs it bisects); tolerance 1e-9 of the
    field's max (the CSV carries 12 digits);
  * with OTHER_OUTDIR (same run at another rank count): same cells (as
    vertex sets), same times, every field equal to 1e-10 of its max.
--frames K checks K evenly spaced frames (first and last included) and
streams the CSV files keeping only those steps, so a full 12 s run (~600
frames, tens of GB of CSV) is checked without loading it; the frame count
is still checked against the CSV's last dump step.
Exits non-zero on the first failed check.  numpy + h5py only (Python 3.7+).
Scope: Chunhui Zhao, 2026-10-14 (P7).
"""
import glob
import os
import re
import sys

import numpy as np
import h5py

FIELDS = ["slip_dip", "slip_strike", "slip_rate_dip", "slip_rate_strike", "slip_rate",
          "traction_dip", "traction_strike", "state_variable", "normal_stress",
          "param_a", "param_Dc"]
# fields rebuilt from the CSV rows (the traction needs sigma_F, not dumped)
CHECKED = ["slip_dip", "slip_strike", "slip_rate_dip", "slip_rate_strike", "slip_rate",
           "state_variable", "normal_stress"]
TOL_CSV, TOL_RANKS = 1e-9, 1e-10


def fail(msg):
    print("FAILED: " + msg)
    sys.exit(1)


def load(d, frames=None):
    """frames: list of step indices to read (None = all)."""
    fs = glob.glob(os.path.join(d, "*.vtkhdf"))
    if len(fs) != 1:
        fail("expected exactly one .vtkhdf in %s, found %s" % (d, fs))
    g = h5py.File(fs[0], "r")["VTKHDF"]
    P = g["Points"][:]
    C = g["Connectivity"][:].reshape(-1, 3)
    times = g["Steps/Values"][:]
    nc = g["NumberOfCells"][:]
    if not (nc[0] == C.shape[0] and np.all(nc[1:] == 0) and np.all(g["Steps/CellOffsets"][:] == 0)):
        fail("the mesh is not stored once (static geometry)")
    sel = list(range(len(times))) if frames is None else list(frames)
    vals = {}
    for k in FIELDS:
        off = g["Steps/PointDataOffsets/" + k][:].astype(np.int64)
        arr = g["PointData/" + k]
        # L2-p0 cell values are written per point: take each cell's first vertex
        rows = []
        for s in sel:
            v = arr[off[s]:off[s] + P.shape[0]]
            if not (np.allclose(v[C[:, 0]], v[C[:, 1]]) and np.allclose(v[C[:, 0]], v[C[:, 2]])):
                fail("field %s is not constant over each cell (step %d)" % (k, s))
            rows.append(v[C[:, 0]])
        vals[k] = np.stack(rows)
    return fs[0], P, C, times, vals


def frame_vectors(d, prefix):
    txt = open(os.path.join(d, prefix + "_station_conventions.txt")).read()
    num = r"(-?[\d.e+-]+)"
    trip = ", ".join([num] * 3)
    n = re.search(r"nhat = pin reference normal \(" + trip + r"\)", txt)
    s = re.search(r"shat = .*?= \(" + trip + r"\)", txt)
    dd = re.search(r"dhat = .*?= \(" + trip + r"\)", txt)
    return tuple(np.array([float(x) for x in m.groups()]) for m in (n, dd, s))


def read_csv_steps(files, steps):
    """CSV rows whose step (first column) is in `steps`, streamed."""
    want = set("%d," % st for st in steps)
    out = []
    for f in files:
        with open(f) as h:
            for line in h:
                c = line.find(",")
                if line[:c + 1] in want:
                    out.append(line)
    return np.array([[float(x) for x in l.split(",")] for l in out]) if out else np.zeros((0, 18))


def last_dump_step(files):
    """(first, last) dump step of a CSV file, from its first and last lines."""
    for f in files:
        if os.path.getsize(f) == 0:
            continue
        with open(f, "rb") as h:
            first = h.readline().decode()
            h.seek(max(0, os.path.getsize(f) - 4096))
            last = h.read().decode().strip().split("\n")[-1]
        return int(float(first.split(",")[0])), int(float(last.split(",")[0]))
    fail("no non-empty CSV fault dump")


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("other", nargs="?")
    ap.add_argument("--prefix", default="tpv102")
    ap.add_argument("--frames", type=int, default=0, help="check K evenly spaced frames (0 = all)")
    a = ap.parse_args()
    d, prefix = a.outdir, a.prefix
    files = sorted(glob.glob(os.path.join(d, prefix + "_fault_r*.csv")))
    every, last = last_dump_step(files)
    with h5py.File(glob.glob(os.path.join(d, "*.vtkhdf"))[0], "r") as h:
        nsteps = len(h["VTKHDF/Steps/Values"])
    if nsteps != last // every:
        fail("steps %d != CSV dump frames %d (dump every %d steps, last %d)"
             % (nsteps, last // every, every, last))
    sel = list(range(nsteps)) if a.frames <= 0 or a.frames >= nsteps else \
        sorted(set(int(round(x)) for x in np.linspace(0, nsteps - 1, a.frames)))
    path, P, C, times, vals = load(d, sel)
    nh, dh, sh = frame_vectors(d, prefix)
    csv_steps = [(i + 1) * every for i in sel]
    rows = read_csv_steps(files, csv_steps)
    nfaces = len(np.unique(np.round(rows[rows[:, 0] == csv_steps[0]][:, 5:8], 3), axis=0))
    print("%s: %.2f MB, %d cells, %d steps (checking %d); CSV dump every %d steps, %d fault faces"
          % (path, os.path.getsize(path) / 1e6, C.shape[0], nsteps, len(sel), every, nfaces))
    csv_t = np.array([rows[rows[:, 0] == s][0, 1] for s in csv_steps])
    if np.max(np.abs(csv_t - times[sel])) > 1e-9 * max(1.0, abs(times[-1])):
        fail("step times differ from the CSV dump times")
    if C.shape[0] != 4 * nfaces:
        fail("cells %d != 4 x fault faces %d" % (C.shape[0], nfaces))
    tri = P[C]
    nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    area = 0.5 * np.linalg.norm(nrm, axis=1)
    if np.any(nrm.dot(nh) <= 0):
        fail("%d inverted sub-triangles (normal . nhat <= 0)" % np.sum(nrm.dot(nh) <= 0))
    pts = tri.reshape(-1, 12, 3)                 # the 4 sub-triangles of a face are consecutive
    face_area, face_cent = [], []
    for b in pts:
        key, cnt = np.unique(np.round(b, 6), axis=0, return_counts=True)
        v = key[cnt == 1]                        # vertices once, edge midpoints 3 times
        if len(v) != 3 or len(key) != 6:
            fail("a face's 4 sub-triangles are not its midpoint subdivision")
        face_area.append(0.5 * np.linalg.norm(np.cross(v[1] - v[0], v[2] - v[0])))
        face_cent.append(v.mean(0))
    face_area = np.array(face_area)
    aerr = np.max(np.abs(area.reshape(-1, 4).sum(1) - face_area) / face_area)
    print("  sub-area sum vs face area: max rel err %.2e; inverted 0 of %d" % (aerr, C.shape[0]))
    if aerr > 1e-6:
        fail("sub-areas do not sum to the face area")

    key_rows = {}
    for i, r in enumerate(rows):
        key_rows.setdefault((r[0],) + tuple(np.round(r[5:8], 3)), []).append(i)

    def nodal(r):
        dv, vv = r[8:11], r[11:14]
        return {"slip_dip": dv.dot(dh), "slip_strike": dv.dot(sh), "slip_rate_dip": vv.dot(dh),
                "slip_rate_strike": vv.dot(sh), "state_variable": r[14], "normal_stress": r[17]}

    err = dict((k, 0.0) for k in CHECKED)
    scale = dict((k, max(np.abs(vals[k]).max(), 1e-300)) for k in CHECKED)
    n_interp = 0
    for si, step in enumerate(csv_steps):
        for fb in range(pts.shape[0]):
            R = rows[key_rows[(step,) + tuple(np.round(face_cent[fb], 3))]]

            def corner(x):
                j = np.where(np.linalg.norm(R[:, 2:5] - x, axis=1) < 1e-6)[0]
                if len(j) == 1:
                    return nodal(R[j[0]]), False
                pairs = [(k, m) for k in range(len(R)) for m in range(k + 1, len(R))
                         if np.linalg.norm(0.5 * (R[k, 2:5] + R[m, 2:5]) - x) < 1e-6]
                if len(pairs) != 1:
                    fail("a corner matches no DOF and no DOF midpoint")
                a, b = nodal(R[pairs[0][0]]), nodal(R[pairs[0][1]])
                return dict((q, 0.5 * (a[q] + b[q])) for q in a), True

            for s in range(4):
                c = 4 * fb + s
                cs = [corner(tri[c, v]) for v in range(3)]
                if si == 0 and any(ip for _, ip in cs):
                    n_interp += 1
                cs = [x for x, _ in cs]
                for q in CHECKED:
                    if q == "slip_rate":
                        exp = np.mean([np.hypot(x["slip_rate_dip"], x["slip_rate_strike"]) for x in cs])
                    else:
                        exp = np.mean([x[q] for x in cs])
                    err[q] = max(err[q], abs(vals[q][si, c] - exp))
    print("  cells with an interpolated (non-DOF) corner: %d of %d" % (n_interp, C.shape[0]))
    print("  cell value vs CSV-rebuilt corner mean, max |err| / max|field|:")
    for q in CHECKED:
        print("    %-18s %.2e" % (q, err[q] / scale[q]))
    if any(err[q] / scale[q] > TOL_CSV for q in CHECKED):
        fail("cell values differ from the CSV corner means (tolerance %.0e)" % TOL_CSV)

    if a.other:
        if len(sel) != nsteps:
            fail("the rank comparison needs all frames (omit --frames)")
        _, P2, C2, t2, v2 = load(a.other)
        if C2.shape != C.shape or len(t2) != len(times):
            fail("the other run has a different cell or step count")

        def order(Pp, Cc):
            cen = np.round(Pp[Cc].mean(1), 4)
            return np.lexsort(cen.T[::-1])
        o1, o2 = order(P, C), order(P2, C2)
        # cells compared as vertex SETS: a partition may list a face's vertices
        # in another local order (orientation is checked above)
        g1 = np.sort(np.round(P[C][o1], 6), axis=1)
        g2 = np.sort(np.round(P2[C2][o2], 6), axis=1)
        s1 = np.sort(np.round(P[C][o1], 6).reshape(len(o1), -1), axis=1)
        s2 = np.sort(np.round(P2[C2][o2], 6).reshape(len(o2), -1), axis=1)
        if not (np.array_equal(s1, s2) and g1.shape == g2.shape):
            fail("the two runs have different cells")
        if not np.array_equal(times, t2):
            fail("the two runs have different step times")
        print("  vs %s: same cells and times" % a.other)
        worst = 0.0
        for q in FIELDS:
            a, b = vals[q][:, o1], v2[q][:, o2]
            e = np.abs(a - b).max() / max(np.abs(a).max(), 1e-300)
            worst = max(worst, e)
            print("    %-18s max |diff| / max|field| %.2e" % (q, e))
        if worst > TOL_RANKS:
            fail("the two rank counts differ (tolerance %.0e)" % TOL_RANKS)
    print("PASSED")


if __name__ == "__main__":
    main()
