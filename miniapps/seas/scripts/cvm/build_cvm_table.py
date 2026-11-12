#!/usr/bin/env python3
"""Build a 3-D material table (Vp, Vs, density) from a community velocity model.

General tool: nothing here belongs to one benchmark.  A case is a preset
(JSON) with four parts:

  frame    how a local point (X, Y, Z) maps to the model's map coordinates:
               easting  = origin_en_m[0] + x_axis_en[0] X + y_axis_en[0] Y
               northing = origin_en_m[1] + x_axis_en[1] X + y_axis_en[1] Y
               depth    = -Z                      (Z up, Z = 0 the surface)
           x_axis_en / y_axis_en are the map directions (east, north) of the
           local X and Y axes.  Leave "frame" out for a table directly in
           map coordinates (X = easting, Y = northing).
  grid     rectilinear axes in the local frame.  Each axis is a list of
           segments [start, stop, step] (m), joined; "depth_m" is positive
           down; "extra_*_planes_m" adds single planes.
  backend  the model query: a command that reads one "easting northing
           depth" line per point on stdin and prints one line per point, in
           the same order; "output_columns_vp_vs_rho" are its 1-based Vp, Vs,
           rho columns.  {name} fields in the command are filled from
           --set name=value (for example the install paths).
  rules    optional.  "min_query_depth_m": points shallower than this are
           queried at this depth.  "floor": where Vp <= if_vp_le or
           Vs <= or_vs_le, set Vp, Vs, rho to the given values.
           Non-positive Vp, Vs or rho returned by the model is always an
           error and stops the run.

Two steps, so that the machine with the model needs numpy only:

  query  build_cvm_table.py query PRESET.json --out TABLE.npz --set k=v ...
         writes TABLE.npz (axes x, y, z and vp, vs, rho after the rules) and
         TABLE.npz.summary.json (preset, grid, min/max before and after the
         rules, floored points, command, sha256 of the .npz).
  pack   build_cvm_table.py pack TABLE.npz TABLE.h5        (needs h5py)
         writes the HDF5 table in the data_projection_v1 layout
         (document/features_dev/data_projection_schema_v1.md) with the
         preset's "table_crs" label, the frame and the rules as attributes,
         and TABLE.h5.sha256.

Presets: tpv34/configs/tpv34_cvm_table.json (SCEC TPV34, CVM-H 15.1.0).
Tables are not committed (size); they are rebuilt from the preset and
checked with the sha256 files.  Chunhui Zhao, TPV34 plan rev 2, 2026-11-12.
"""

import argparse
import hashlib
import json
import subprocess
import sys

import numpy as np


def build_axis(segments, extra=(), name="axis"):
    """Join [start, stop, step] segments (ends included) and extra planes."""
    parts = [np.asarray(extra, dtype=float)]
    for a, b, h in segments:
        n = round((b - a) / h)
        if n < 1 or abs(a + n * h - b) > 1e-9 * max(1.0, abs(b)):
            sys.exit("%s: segment [%g, %g, %g] is not a whole number of steps" % (name, a, b, h))
        parts.append(a + h * np.arange(n + 1))
    ax = np.unique(np.concatenate(parts))
    if ax.size < 2:
        sys.exit("%s: needs at least two planes" % name)
    return ax


def to_map(frame, X, Y):
    """Local (X, Y) to (easting, northing)."""
    if frame is None:
        return X, Y
    (e0, n0), (xe, xn), (ye, yn) = frame["origin_en_m"], frame["x_axis_en"], frame["y_axis_en"]
    return e0 + xe * X + ye * Y, n0 + xn * X + yn * Y


def run_backend(cmd, be, e, n, d):
    """One model call; returns raw Vp, Vs, rho."""
    fmt = be.get("input_format", "%.3f %.3f %.3f")
    text = "\n".join(fmt % t for t in zip(e, n, d)) + "\n"
    r = subprocess.run(cmd, input=text, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit("model query failed (exit %d): %s" % (r.returncode, r.stderr[-2000:]))
    rows = [ln.split() for ln in r.stdout.splitlines() if ln.strip()]
    if len(rows) != len(e):
        sys.exit("model query returned %d lines for %d points" % (len(rows), len(e)))
    ncol = be.get("output_column_count")
    if ncol is not None and any(len(row) != ncol for row in rows):
        sys.exit("model query: a line does not have %d columns" % ncol)
    cols = [c - 1 for c in be["output_columns_vp_vs_rho"]]
    out = np.array([[row[c] for c in cols] for row in rows], dtype=float)
    return out[:, 0], out[:, 1], out[:, 2]


def apply_rules(rules, vp, vs, rho, where=""):
    bad = (vp <= 0) | (vs <= 0) | (rho <= 0) | ~np.isfinite(vp + vs + rho)
    if bad.any():
        sys.exit("the model returned non-positive or non-finite Vp, Vs or rho at %d points %s"
                 % (bad.sum(), where))
    low = np.zeros(vp.shape, dtype=bool)
    fl = rules.get("floor")
    if fl:
        low = (vp <= fl["if_vp_le"]) | (vs <= fl["or_vs_le"])
        vp, vs, rho = vp.copy(), vs.copy(), rho.copy()
        vp[low], vs[low], rho[low] = fl["set_vp"], fl["set_vs"], fl["set_rho"]
    return vp, vs, rho, low


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def cmd_query(a):
    with open(a.preset) as f:
        ps = json.load(f)
    frame, grid, be, rules = ps.get("frame"), ps["grid"], ps["backend"], ps.get("rules", {})
    subs = dict(kv.split("=", 1) for kv in a.set)
    try:
        cmd = [c.format(**subs) for c in be["command"]]
    except KeyError as k:
        sys.exit("the backend command needs --set %s=..." % k.args[0])
    x = build_axis(grid["x_m"], grid.get("extra_x_planes_m", ()), "x_m")
    y = build_axis(grid["y_m"], grid.get("extra_y_planes_m", ()), "y_m")
    depth = build_axis(grid["depth_m"], grid.get("extra_depth_planes_m", ()), "depth_m")
    if depth[0] != 0.0:
        sys.exit("depth_m must start at 0 (the surface)")
    z = 0.0 - depth[::-1]                       # elevation, increasing, ends at 0
    nx, ny, nz = len(x), len(y), len(z)
    print("grid %d x %d x %d = %d points" % (nx, ny, nz, nx * ny * nz), flush=True)
    dmin = rules.get("min_query_depth_m", 0.0)
    vp = np.empty((nx, ny, nz))
    vs = np.empty_like(vp)
    rho = np.empty_like(vp)
    raw = {k: [np.inf, -np.inf] for k in ("vp", "vs", "rho")}
    floored = 0
    _, Yg, Zg = np.meshgrid(np.zeros(1), y, z, indexing="ij")
    for i0 in range(0, nx, a.planes_per_call):  # a few X planes per model call
        xs = x[i0:i0 + a.planes_per_call]
        e, n = to_map(frame, np.repeat(xs, ny * nz), np.tile(Yg.ravel(), len(xs)))
        d = np.maximum(-np.tile(Zg.ravel(), len(xs)), dmin)
        rvp, rvs, rrho = run_backend(cmd, be, e, n, d)
        for k, v in (("vp", rvp), ("vs", rvs), ("rho", rrho)):
            raw[k] = [min(raw[k][0], float(v.min())), max(raw[k][1], float(v.max()))]
        pvp, pvs, prho, low = apply_rules(rules, rvp, rvs, rrho,
                                          "in the planes X = %g to %g m" % (xs[0], xs[-1]))
        floored += int(low.sum())
        sl = slice(i0, i0 + len(xs))
        vp[sl], vs[sl], rho[sl] = (v.reshape(len(xs), ny, nz) for v in (pvp, pvs, prho))
        print("  planes %d to %d of %d" % (i0 + 1, i0 + len(xs), nx), flush=True)
    meta = {"preset_file": a.preset, "preset": ps, "command": cmd}
    np.savez(a.out, x=x, y=y, z=z, vp=vp, vs=vs, rho=rho, meta=json.dumps(meta))
    summary = dict(meta)
    summary.update({
        "grid": [nx, ny, nz], "points": nx * ny * nz,
        "x_range_m": [x[0], x[-1]], "y_range_m": [y[0], y[-1]], "z_range_m": [z[0], z[-1]],
        "raw_min_max": raw,
        "table_min_max": {k: [float(v.min()), float(v.max())]
                          for k, v in (("vp", vp), ("vs", vs), ("rho", rho))},
        "floored_points": floored, "npz_sha256": sha256(a.out),
    })
    with open(a.out + ".summary.json", "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary, indent=1))


def cmd_pack(a):
    import datetime
    import h5py
    t = np.load(a.npz)
    meta = json.loads(str(t["meta"])) if "meta" in t.files else {}
    ps = meta.get("preset", {})
    with h5py.File(a.h5, "w") as f:
        f.attrs["schema_version"] = "data_projection_v1"
        f.attrs["crs"] = a.crs or ps.get("table_crs", "mesh-local")
        f.attrs["units"] = "m"
        f.attrs["z_positive"] = "elevation"
        f.attrs["created_at"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        f.attrs["source"] = ps.get("description", "")
        f.attrs["frame_json"] = json.dumps(ps.get("frame"))
        f.attrs["rules_json"] = json.dumps(ps.get("rules", {}))
        f.attrs["command_json"] = json.dumps(meta.get("command"))
        f.attrs["npz_sha256"] = sha256(a.npz)
        for k in ("x", "y", "z"):
            f.create_dataset("/grid/" + k, data=t[k].astype("f8"))
        for k, unit in (("vp", "m/s"), ("vs", "m/s"), ("rho", "kg/m3")):
            v = t[k].astype("f8")
            assert v.shape == (len(t["x"]), len(t["y"]), len(t["z"])) and np.isfinite(v).all()
            ds = f.create_dataset("/fields/" + k, data=v)
            ds.attrs["units"] = unit
            ds.attrs["min_value"] = float(v.min())
            ds.attrs["max_value"] = float(v.max())
    digest = sha256(a.h5)
    with open(a.h5 + ".sha256", "w") as f:
        f.write("%s  %s\n" % (digest, a.h5.split("/")[-1]))
    print("wrote %s  sha256 %s" % (a.h5, digest))


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    q = sub.add_parser("query", help="sample the model on the preset's grid, write the .npz")
    q.add_argument("preset", help="preset JSON (frame, grid, backend, rules)")
    q.add_argument("--out", required=True, help="output .npz")
    q.add_argument("--set", action="append", default=[], metavar="NAME=VALUE",
                   help="fill a {NAME} field of the backend command (repeatable)")
    q.add_argument("--planes-per-call", type=int, default=40,
                   help="X planes sent to one model call (memory against start-up cost)")
    q.set_defaults(func=cmd_query)
    k = sub.add_parser("pack", help="write the HDF5 table from the .npz")
    k.add_argument("npz")
    k.add_argument("h5")
    k.add_argument("--crs", default=None, help="override the preset's table_crs label")
    k.set_defaults(func=cmd_pack)
    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
