#!/usr/bin/env python3
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# check_snapshots.py — reads the snapshot files of seas_fd_cycle_driver
# (cycle/snapshots.hpp) back and checks them:
#   1. layout: the VTKHDF datasets and offsets, /fd_cycle per snapshot,
#      times strictly increasing (h5py);
#   2. ParaView's reader: vtkHDFReader sees every snapshot with its time and
#      the same values (python vtk);
#   3. the fault values against the station files: at every snapshot time
#      that has a station row, the fault point at the station (or the face
#      containing it) holds that row's slip, log10 |V|, shear stress and
#      log10 theta (the station rows are printed with 12 digits);
#   4. the schedule: QD cadence snapshots on new dt_qd marks, dynamic ones on
#      new dt_dyn marks of the phase time; the last QD snapshot before a
#      dynamic phase is "to dynamic" (or "start" at t = 0) and the last one
#      of the phase "end of dynamic"; with --complete (a run that reached
#      t_final) the last snapshot is "end";
#   5. with --against DIR2: the same snapshots in DIR2's files (count,
#      /fd_cycle, times), and each field's largest difference relative to
#      its largest magnitude, at most --tol (default 0: identical).
# Prints the snapshot table.  Exit 1 if a check fails.
#
#   python3 tests/cycle/check_snapshots.py --dir RUN --prefix snap --deck DECK.toml
#          [--dt-dyn-step 1.71963e-3]   (the dynamic dt from run.log, for check 4)
#          [--against RUN2 [--tol 1e-8]] [--no-vtk] [--complete]

import argparse
import math
import sys
import tomllib

import h5py
import numpy as np

REASONS = {0: "cadence", 1: "start", 2: "to dynamic", 3: "end of dynamic", 4: "end"}


def fail(msg):
    print("FAIL:", msg)
    return 1


def check_layout(path, fields):
    bad = 0
    with h5py.File(path, "r") as f:
        v = f["VTKHDF"]
        assert v.attrs["Type"].decode() == "UnstructuredGrid", v.attrs["Type"]
        npts = int(v["NumberOfPoints"][0])
        ncell = int(v["NumberOfCells"][0])
        st = v["Steps"]
        ns = int(np.ravel(st.attrs["NSteps"])[0])
        t = st["Values"][:]
        if len(t) != ns:
            bad += fail(f"{path}: NSteps {ns} but {len(t)} times")
        if ns > 1 and not np.all(np.diff(t) > 0):
            bad += fail(f"{path}: times not strictly increasing")
        for n in ("PartOffsets", "PointOffsets", "CellOffsets", "ConnectivityIdOffsets"):
            if not np.all(st[n][:] == 0) or len(st[n]) != ns:
                bad += fail(f"{path}: Steps/{n} is not {ns} zeros")
        if not np.all(st["NumberOfParts"][:] == 1):
            bad += fail(f"{path}: Steps/NumberOfParts is not all 1")
        if v["Points"].shape != (npts, 3) or len(v["Types"]) != ncell \
                or len(v["Offsets"]) != ncell + 1:
            bad += fail(f"{path}: grid datasets have the wrong shapes")
        for name, nc in fields:
            d = v["PointData"][name]
            want = (ns * npts,) if nc == 1 else (ns * npts, nc)
            if d.shape != want:
                bad += fail(f"{path}: PointData/{name} shape {d.shape}, expected {want}")
            off = st["PointDataOffsets"][name][:]
            if not np.array_equal(off, np.arange(ns) * npts):
                bad += fail(f"{path}: PointDataOffsets/{name} is not k * npts")
        for n in ("regime", "phase", "reason", "step"):
            if len(f["fd_cycle"][n]) != ns:
                bad += fail(f"{path}: /fd_cycle/{n} has {len(f['fd_cycle'][n])} rows, not {ns}")
        meta = {n: f["fd_cycle"][n][:] for n in ("regime", "phase", "reason", "step")}
        print(f"{path}: {npts} points, {ncell} cells (types {sorted(set(v['Types'][:]))}), "
              f"{ns} snapshots, t = {t[0]:.6g} .. {t[-1]:.6g} s; layout {'OK' if not bad else 'BAD'}")
    return bad, t, meta, npts


def check_vtk(path, fields, t):
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    bad = 0
    r = vtk.vtkHDFReader()
    r.SetFileName(path)
    r.UpdateInformation()
    info = r.GetOutputInformation(0)
    key = vtk.vtkStreamingDemandDrivenPipeline.TIME_STEPS()
    ts = np.array([info.Get(key, i) for i in range(info.Length(key))]) if info.Has(key) else np.array([])
    if len(ts) != len(t) or not np.allclose(ts, t, rtol=0, atol=0):
        bad += fail(f"{path}: vtkHDFReader sees {len(ts)} time steps, the file has {len(t)}")
    with h5py.File(path, "r") as f:
        npts = int(f["VTKHDF"]["NumberOfPoints"][0])
        for k in sorted({0, len(t) // 2, len(t) - 1}):
            r.UpdateTimeStep(t[k])
            out = r.GetOutputDataObject(0)
            if out.GetNumberOfPoints() != npts:
                bad += fail(f"{path}: step {k}: VTK has {out.GetNumberOfPoints()} points")
                continue
            for name, nc in fields:
                a = vtk_to_numpy(out.GetPointData().GetArray(name))
                h = f["VTKHDF"]["PointData"][name][k * npts:(k + 1) * npts]
                if not np.array_equal(a.reshape(h.shape), h):
                    bad += fail(f"{path}: step {k}: VTK's {name} differs from the file")
    print(f"{path}: vtkHDFReader {len(ts)} time steps, values at 3 steps {'OK' if not bad else 'BAD'}")
    return bad


def read_station(path):
    rows = np.loadtxt(path, comments="#", ndmin=2)
    return rows


def check_fault_vs_stations(fpath, run_dir, prefix, deck):
    bad = 0
    stations = deck["cycle"]["output"]["station"]
    with h5py.File(fpath, "r") as f:
        v = f["VTKHDF"]
        X = v["Points"][:]
        npts = X.shape[0]
        conn = v["Connectivity"][:]
        types = v["Types"][:]
        t = v["Steps"]["Values"][:]
        pd = v["PointData"]
        names = ["slip_strike", "slip_dip", "slip_rate_strike", "slip_rate_dip",
                 "shear_stress_strike", "shear_stress_dip", "state"]
        data = {n: pd[n][:].reshape(len(t), npts) for n in names}
    tri = conn.reshape(-1, 3) if np.all(types == 5) else None
    nmatch = 0
    worst = 0.0
    for s in stations:
        xs = np.array(s["x"], dtype=float)
        rows = read_station(f"{run_dir}/{prefix}_{s['name']}.txt")
        # Candidate interpolants: points at the station, else containing faces.
        d = np.linalg.norm(X - xs, axis=1)
        cands = [("point", [int(j)], [1.0]) for j in np.nonzero(d < 1e-6 * max(1.0, d.max()))[0]]
        if not cands and tri is not None:
            for c in tri:
                P = X[c]
                # barycentric in the plane (x, z) of the fault y = 0
                A = np.array([[P[0, 0], P[1, 0], P[2, 0]], [P[0, 2], P[1, 2], P[2, 2]],
                              [1.0, 1.0, 1.0]])
                w = np.linalg.solve(A, np.array([xs[0], xs[2], 1.0]))
                if np.all(w >= -1e-9):
                    cands.append(("face", list(c), list(w)))
        if not cands:
            bad += fail(f"station {s['name']}: no fault point or face at {xs}")
            continue
        best = None
        for kind, idx, w in cands:
            err = 0.0
            n = 0
            for k, tk in enumerate(t):
                m = np.nonzero(np.abs(rows[:, 0] - tk) <= 1e-12 * max(1.0, abs(tk)))[0]
                if len(m) == 0:
                    continue
                row = rows[m[-1]]
                e = {nm: float(np.dot(w, data[nm][k, idx])) for nm in names}
                sl = lambda x: math.log10(x) if x > 0 else -300.0
                snap = [e["slip_strike"], e["slip_dip"], sl(abs(e["slip_rate_strike"])),
                        sl(abs(e["slip_rate_dip"])), e["shear_stress_strike"],
                        e["shear_stress_dip"], sl(e["state"])]
                for a, b in zip(snap, row[1:8]):
                    err = max(err, abs(a - b) / (max(abs(a), abs(b)) + 1e-12))
                n += 1
            if best is None or err < best[0]:
                best = (err, n, kind)
        err, n, kind = best
        nmatch += n
        worst = max(worst, err)
        status = "OK" if err <= 1e-10 and n > 0 else "BAD"
        if status == "BAD":
            bad += 1
        print(f"  {s['name']:22s} {kind:5s} {len(cands)} candidate(s): {n} snapshot times "
              f"with a station row, max rel. diff {err:.2e}  {status}")
    print(f"{fpath}: {nmatch} station rows matched, worst {worst:.2e} "
          f"{'OK' if not bad else 'BAD'}")
    return bad


def check_schedule(t, meta, kind_cfg, dt_dyn_step, complete):
    bad = 0
    dt_qd = float(kind_cfg.get("dt_qd", 3.15576e8))
    dt_dyn = float(kind_cfg.get("dt_dyn", 0.1))
    last_qd = -1
    last_dyn = 0
    prev_reg = None
    print(f"   {'#':>4} {'t (s)':>16} {'regime':>8} {'phase':>6} {'step':>6}  reason")
    for k in range(len(t)):
        reg, ph, rs, stp = (int(meta[n][k]) for n in ("regime", "phase", "reason", "step"))
        print(f"   {k:4d} {t[k]:16.9e} {('QD', 'dynamic')[reg]:>8} {ph:6d} {stp:6d}  {REASONS[rs]}")
        if reg == 0:
            mark = math.floor(t[k] / dt_qd * (1 + 1e-12) + 1e-9)
            if rs == 0 and mark <= last_qd:
                bad += fail(f"snapshot {k}: QD cadence at mark {mark}, not after {last_qd}")
            last_qd = max(last_qd, mark)
        else:
            if prev_reg != 1:
                last_dyn = 0
            if dt_dyn_step:
                mark = math.floor(stp * dt_dyn_step / dt_dyn * (1 + 1e-12) + 1e-9)
                if rs == 0 and mark <= last_dyn:
                    bad += fail(f"snapshot {k}: dynamic cadence at mark {mark}, not after {last_dyn}")
                last_dyn = max(last_dyn, mark)
        prev_reg = reg
    reg = [int(r) for r in meta["regime"]]
    rs = [int(r) for r in meta["reason"]]
    for k in range(len(t) - 1):
        if reg[k] == 0 and reg[k + 1] == 1 and rs[k] not in (1, 2):
            bad += fail(f"snapshot {k}: the last QD one before a dynamic phase is "
                        f"'{REASONS[rs[k]]}', not 'start' or 'to dynamic'")
        if reg[k] == 1 and reg[k + 1] == 0 and rs[k] != 3:
            bad += fail(f"snapshot {k}: the last one of a dynamic phase is "
                        f"'{REASONS[rs[k]]}', not 'end of dynamic'")
    if complete and len(t) and rs[-1] != 4:
        bad += fail(f"the last snapshot is '{REASONS[rs[-1]]}', not 'end'")
    return bad


def compare_files(p1, p2, fields, tol):
    bad = 0
    with h5py.File(p1, "r") as f1, h5py.File(p2, "r") as f2:
        v1, v2 = f1["VTKHDF"], f2["VTKHDF"]
        t1, t2 = v1["Steps"]["Values"][:], v2["Steps"]["Values"][:]
        if len(t1) != len(t2):
            return fail(f"{p1}: {len(t1)} snapshots, {p2}: {len(t2)}")
        if not np.array_equal(v1["Points"][:], v2["Points"][:]):
            bad += fail("the grids differ")
        for n in ("regime", "phase", "reason", "step"):
            if not np.array_equal(f1["fd_cycle"][n][:], f2["fd_cycle"][n][:]):
                bad += fail(f"/fd_cycle/{n} differs")
        dt = np.max(np.abs(t1 - t2) / np.maximum(np.abs(t2), 1e-300)) if len(t1) else 0.0
        if dt > max(tol, 0.0) and not (tol == 0 and dt == 0):
            bad += fail(f"times differ by {dt:.2e} (relative)")
        out = [f"times {dt:.1e}"]
        for name, nc in fields:
            a = v1["PointData"][name][:]
            b = v2["PointData"][name][:]
            scale = max(np.max(np.abs(b)), 1e-300)
            d = float(np.max(np.abs(a - b))) / scale
            out.append(f"{name} {d:.1e}")
            if d > tol:
                bad += fail(f"{name}: max difference {d:.2e} of its largest magnitude > {tol}")
        print(f"{p1} vs {p2}: {len(t1)} snapshots; " + ", ".join(out)
              + ("  OK" if not bad else "  BAD"))
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--deck", required=True)
    ap.add_argument("--dt-dyn-step", type=float, default=0.0)
    ap.add_argument("--no-vtk", action="store_true")
    ap.add_argument("--against", default="")
    ap.add_argument("--tol", type=float, default=0.0)
    ap.add_argument("--complete", action="store_true")
    a = ap.parse_args()
    with open(a.deck, "rb") as fh:
        deck = tomllib.load(fh)
    snaps = deck["cycle"].get("snapshots", {})
    bad = 0
    kinds = [("fault", [("slip_strike", 1), ("slip_dip", 1), ("slip_rate_strike", 1),
                        ("slip_rate_dip", 1), ("state", 1), ("normal_stress", 1),
                        ("shear_stress_strike", 1), ("shear_stress_dip", 1)]),
             ("surface", [("displacement", 3), ("velocity", 3)])]
    for kind, fields in kinds:
        if kind not in snaps or not snaps[kind].get("enabled", True):
            continue
        path = f"{a.dir}/{a.prefix}_{kind}_snapshots.vtkhdf"
        b, t, meta, npts = check_layout(path, fields)
        bad += b
        if not a.no_vtk:
            bad += check_vtk(path, fields, t)
        if kind == "fault":
            bad += check_fault_vs_stations(path, a.dir, a.prefix, deck)
        print(f"{path}: schedule (dt_qd {snaps[kind].get('dt_qd')}, dt_dyn {snaps[kind].get('dt_dyn')})")
        bad += check_schedule(t, meta, snaps[kind], a.dt_dyn_step, a.complete)
        if a.against:
            bad += compare_files(path, f"{a.against}/{a.prefix}_{kind}_snapshots.vtkhdf",
                                 fields, a.tol)
    print("ALL OK" if bad == 0 else f"{bad} FAILURE(S)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
