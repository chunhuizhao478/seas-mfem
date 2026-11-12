#!/usr/bin/env python3
"""Build the SCEC TPV34 material table from the SCEC CVM-H 15.1.0 velocity model.

TPV34 (TPV34_Description_v10.pdf, Part 2) takes Vp, Vs and density from CVM-H
everywhere in the domain.  This script samples CVM-H on one rectilinear grid
given in the MESH frame of the TPV34 decks and writes it as a table the
driver interpolates in.

Frames
    mesh  X along strike, Y fault normal, Z up (m); fault on Y = 0, free
          surface Z = 0 (as TPV205, TPV104, TPV31).
    SCEC  x along strike, y depth, z fault normal:  x = X, y = -Z, z = Y.

Query point of a grid point (description Part 2; UTM zone 11, NAD 1927):
    easting  E = 648446  - 0.5802386 x - 0.8144465 z
    northing N = 3625237 + 0.8144465 x - 0.5802386 z
    depth    d = max(y, 100 m)

Rules applied to every returned point (description Part 2):
    error  Vp <= 0, Vs <= 0 or rho <= 0: stop.
    floor  Vp <= 2984 m/s or Vs <= 1400 m/s:
           Vp = 2984 m/s, Vs = 1400 m/s, rho = 2220.34 kg/m3.

Two steps, so that the cluster step needs numpy only:

    query  (where CVM-H is built: Linux; CVM-H does not run on macOS)
        make_tpv34_material_table.py query --vx-lite $CVMH/bin/vx_lite \
            --model $CVMH/model --out tpv34_cvmh_table.npz
      writes the .npz (axes + vp, vs, rho after the rules) and
      <out>.summary.json (grid, min/max before and after the floor, number of
      floored points, vx_lite path, sha256 of the .npz).

    pack   (anywhere with h5py)
        make_tpv34_material_table.py pack tpv34_cvmh_table.npz tpv34_cvmh_table.h5
      writes the HDF5 table in the data_projection_v1 layout
      (document/features_dev/data_projection_schema_v1.md) with
      crs = "mesh-local", and <h5>.sha256.

Default grid: 200 m for |X|, |Y| <= 25 km and depth <= 15 km, 1000 m outside,
out to |X|, |Y| = 80 km and 80 km depth, plus one plane at Z = -100 m (the
depth clamp d = max(y, 100 m) makes the model constant above it).

The tables are not committed (size); they are rebuilt from this script and
checked with the sha256 files.  Chunhui Zhao, TPV34 plan, 2026-11-12.
"""

import argparse
import hashlib
import json
import subprocess
import sys

import numpy as np

E0, N0 = 648446.0, 3625237.0
C, S = 0.5802386, 0.8144465
DEPTH_MIN = 100.0
VP_FLOOR, VS_FLOOR, RHO_FLOOR = 2984.0, 1400.0, 2220.34


def axis_two_zone(half_fine, h_fine, half_all, h_coarse):
    """Symmetric axis: h_fine inside +-half_fine, h_coarse out to +-half_all."""
    fine = np.arange(-half_fine, half_fine + 0.5 * h_fine, h_fine)
    outer = np.arange(half_fine + h_coarse, half_all + 0.5 * h_coarse, h_coarse)
    ax = np.concatenate([-outer[::-1], fine, outer])
    assert np.all(np.diff(ax) > 0) and ax[0] == -half_all and ax[-1] == half_all
    return ax


def axis_depth(depth_fine, h_fine, depth_all, h_coarse):
    """Elevation axis Z (increasing, ends at 0) with a plane at -DEPTH_MIN."""
    d = np.concatenate([np.arange(0.0, depth_fine + 0.5 * h_fine, h_fine),
                        np.arange(depth_fine + h_coarse, depth_all + 0.5 * h_coarse, h_coarse),
                        [DEPTH_MIN]])
    z = 0.0 - np.unique(d)[::-1]
    assert np.all(np.diff(z) > 0) and z[0] == -depth_all and z[-1] == 0.0
    return z


def query_points(X, Y, Z):
    """Easting, northing, depth (description Part 2) of mesh points."""
    return E0 - C * X - S * Y, N0 + S * X - C * Y, np.maximum(-Z, DEPTH_MIN)


def run_vx_lite(vx_lite, model, e, n, d):
    """One vx_lite call; returns raw Vp, Vs, rho (columns 17, 18, 19)."""
    text = "\n".join("%.3f %.3f %.3f" % t for t in zip(e, n, d)) + "\n"
    r = subprocess.run([vx_lite, "-s", "-z", "dep", "-m", model], input=text,
                       capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit("vx_lite failed (exit %d): %s" % (r.returncode, r.stderr[-2000:]))
    out = np.array([ln.split()[-3:] for ln in r.stdout.splitlines() if ln.strip()], dtype=float)
    if out.shape != (len(e), 3):
        sys.exit("vx_lite returned %d lines for %d points" % (out.shape[0], len(e)))
    return out[:, 0], out[:, 1], out[:, 2]


def apply_rules(vp, vs, rho, where=""):
    bad = (vp <= 0) | (vs <= 0) | (rho <= 0)
    if bad.any():
        sys.exit("CVM-H returned non-positive Vp, Vs or rho at %d points %s (description: "
                 "this is an error, stop)" % (bad.sum(), where))
    low = (vp <= VP_FLOOR) | (vs <= VS_FLOOR)
    vp, vs, rho = vp.copy(), vs.copy(), rho.copy()
    vp[low], vs[low], rho[low] = VP_FLOOR, VS_FLOOR, RHO_FLOOR
    return vp, vs, rho, low


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def cmd_query(a):
    x = axis_two_zone(a.fine_half, a.h_fine, a.box_half, a.h_coarse)
    y = axis_two_zone(a.fine_half, a.h_fine, a.box_half, a.h_coarse)
    z = axis_depth(a.fine_depth, a.h_fine, a.box_depth, a.h_coarse)
    nx, ny, nz = len(x), len(y), len(z)
    print("grid %d x %d x %d = %d points" % (nx, ny, nz, nx * ny * nz), flush=True)
    vp = np.empty((nx, ny, nz))
    vs = np.empty_like(vp)
    rho = np.empty_like(vp)
    raw = {k: [np.inf, -np.inf] for k in ("vp", "vs", "rho")}
    floored = 0
    Xg, Yg, Zg = np.meshgrid(np.zeros(1), y, z, indexing="ij")
    for i0 in range(0, nx, a.planes_per_call):     # a few X planes per vx_lite call
        xs = x[i0:i0 + a.planes_per_call]
        Xc = np.repeat(xs, ny * nz)
        e, n, d = query_points(Xc, np.tile(Yg.ravel(), len(xs)), np.tile(Zg.ravel(), len(xs)))
        rvp, rvs, rrho = run_vx_lite(a.vx_lite, a.model, e, n, d)
        for k, v in (("vp", rvp), ("vs", rvs), ("rho", rrho)):
            raw[k] = [min(raw[k][0], v.min()), max(raw[k][1], v.max())]
        pvp, pvs, prho, low = apply_rules(rvp, rvs, rrho, "in the planes X = %g to %g m"
                                          % (xs[0], xs[-1]))
        floored += int(low.sum())
        sl = slice(i0, i0 + len(xs))
        vp[sl], vs[sl], rho[sl] = (v.reshape(len(xs), ny, nz) for v in (pvp, pvs, prho))
        print("  planes %d to %d of %d" % (i0 + 1, i0 + len(xs), nx), flush=True)
    np.savez(a.out, x=x, y=y, z=z, vp=vp, vs=vs, rho=rho)
    summary = {
        "grid": [nx, ny, nz],
        "x_range_m": [x[0], x[-1]], "y_range_m": [y[0], y[-1]], "z_range_m": [z[0], z[-1]],
        "h_fine_m": a.h_fine, "h_coarse_m": a.h_coarse,
        "fine_half_m": a.fine_half, "fine_depth_m": a.fine_depth,
        "raw_min_max": raw,
        "table_min_max": {"vp": [vp.min(), vp.max()], "vs": [vs.min(), vs.max()],
                          "rho": [rho.min(), rho.max()]},
        "floored_points": floored, "points": nx * ny * nz,
        "vx_lite": a.vx_lite, "model": a.model, "npz_sha256": sha256(a.out),
    }
    with open(a.out + ".summary.json", "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary, indent=1))


def cmd_pack(a):
    import datetime
    import h5py
    t = np.load(a.npz)
    with h5py.File(a.h5, "w") as f:
        f.attrs["schema_version"] = "data_projection_v1"
        f.attrs["crs"] = "mesh-local"
        f.attrs["units"] = "m"
        f.attrs["z_positive"] = "elevation"
        f.attrs["created_at"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        f.attrs["source"] = "SCEC CVM-H 15.1.0 (vx_lite -s -z dep), TPV34 description Part 2"
        f.attrs["mesh_tag"] = "tpv34: X strike, Y fault normal, Z up; fault on Y = 0"
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
    print("wrote %s  sha256 %s  (npz sha256 %s)" % (a.h5, digest, sha256(a.npz)))


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    q = sub.add_parser("query", help="sample CVM-H with vx_lite, write the .npz")
    q.add_argument("--vx-lite", required=True, help="path to CVM-H's vx_lite")
    q.add_argument("--model", required=True, help="CVM-H model directory")
    q.add_argument("--out", required=True, help="output .npz")
    q.add_argument("--h-fine", type=float, default=200.0)
    q.add_argument("--h-coarse", type=float, default=1000.0)
    q.add_argument("--fine-half", type=float, default=25e3, help="fine zone |X|, |Y| (m)")
    q.add_argument("--fine-depth", type=float, default=15e3, help="fine zone depth (m)")
    q.add_argument("--box-half", type=float, default=80e3, help="table |X|, |Y| (m)")
    q.add_argument("--box-depth", type=float, default=80e3, help="table depth (m)")
    q.add_argument("--planes-per-call", type=int, default=40,
                   help="X planes sent to one vx_lite call (memory against start-up cost)")
    q.set_defaults(func=cmd_query)
    k = sub.add_parser("pack", help="write the HDF5 table from the .npz")
    k.add_argument("npz")
    k.add_argument("h5")
    k.set_defaults(func=cmd_pack)
    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
