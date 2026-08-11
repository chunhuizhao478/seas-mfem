#!/usr/bin/env python3
"""census_box.py -- locate and classify frequency-gate failures on a ShakeOut-box mesh.

The gate (locked convention): f = (p/4) * Vs/dx with dx the element MAX edge and
Vs = sqrt(mu/rho) NEAREST-GRID at the element BARYCENTRE.  A mesh "resolves
f_target at order p" when Vs/dx >= 4*f_target/p everywhere.

    p3 @ 0.5 Hz -> 0.6667        p5 @ 1.0 Hz -> 0.8000

What this adds over verify_extension_lean.py's P8:

  * splits parent block vs collar block by tet index (the merge appends the
    collar after the parent, so the split is exact);
  * depth histogram of the failures, and how many rest ON the free surface --
    the two populations need completely different cures;
  * the empirical barycentre-depth / max-edge ratio for surface-resting cells,
    which is the number the self-consistent size field needs and which is
    otherwise guessed at (-h/3 was an estimate, not a measurement);
  * dumps every failing cell (barycentre, dx, Vs, ratio) to npz so the repair
    stage does not have to re-stream a 6 GB file to find its seeds.

Streaming: nothing bigger than one chunk of tets is ever materialised.
"""
import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]
                       / "meshing_shakeoutbox_intermediate" / "code"))
from collar_lib import VsGrid, tet_edge_lengths  # noqa: E402

CH = 1_000_000
CVM = ("/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_"
       "CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--cvm", default=CVM)
    ap.add_argument("--gate", type=float, action="append", required=True,
                    help="repeatable; e.g. --gate 0.6667 --gate 0.8")
    ap.add_argument("--parent-tets", type=int, default=0,
                    help="tets [0, N) are the frozen parent block")
    ap.add_argument("--dump", default=None, help="npz of failing cells (worst gate)")
    ap.add_argument("--chunk", type=int, default=CH)
    ap.add_argument("--muscal", action="store_true",
                    help="judge the gate against MUSCAL (the SOURCE model, 50 m "
                         "depth step near the surface) instead of the deck's "
                         "resampled nc (250 m uniform). The deck file's binning "
                         "hands every barycentre in the top 125 m the z=0 value, "
                         "which overstates the requirement ~2.2x there.")
    a = ap.parse_args()

    gates = sorted(a.gate)
    gmax = gates[-1]                      # strictest gate drives the dump

    if a.muscal:
        from muscal_vs import MuscalVs
        vs = MuscalVs(backup=VsGrid(a.cvm))
    else:
        vs = VsGrid(a.cvm)
    with h5py.File(a.mesh, "r") as f:
        G = f["geometry"][:]
        nt = f["connect"].shape[0]
    print(f"mesh   {a.mesh}")
    print(f"       {nt:,} tets / {len(G):,} verts")
    print(f"parent block [0, {a.parent_tets:,})   collar block [{a.parent_tets:,}, {nt:,})\n",
          flush=True)

    blocks = ("parent", "collar")
    st = {g: {b: dict(n=0, nfail=0, worst=np.inf) for b in blocks} for g in gates}
    # failure dump (strictest gate only)
    fb, fx, fdx, fvs, fr, fsurf = [], [], [], [], [], []
    # empirical surface-cell geometry: |bary z| / max edge, for cells touching z=0
    surf_ratio = []
    surf_dx = []
    ncell = {b: 0 for b in blocks}

    with h5py.File(a.mesh, "r") as f:
        conn = f["connect"]
        for s0 in range(0, nt, a.chunk):
            s1 = min(s0 + a.chunk, nt)
            T = conn[s0:s1].astype(np.int64)
            P = G[T]                                    # (n,4,3)
            dmax = tet_edge_lengths(G, T).max(1)
            bary = P.mean(1)
            v = vs.at(bary)
            r = v / dmax
            is_parent = np.arange(s0, s1) < a.parent_tets
            zmax = P[:, :, 2].max(1)
            touches_top = zmax > -1e-6

            for b, m in (("parent", is_parent), ("collar", ~is_parent)):
                if not m.any():
                    continue
                ncell[b] += int(m.sum())
                for g in gates:
                    d = st[g][b]
                    d["n"] += int(m.sum())
                    d["nfail"] += int((r[m] < g).sum())
                    d["worst"] = min(d["worst"], float(r[m].min()))

            if touches_top.any():
                surf_ratio.append((-bary[touches_top, 2]) / dmax[touches_top])
                surf_dx.append(dmax[touches_top])

            bad = r < gmax
            if bad.any():
                fb.append(bary[bad].astype(np.float32))
                fx.append(np.nonzero(bad)[0].astype(np.int64) + s0)
                fdx.append(dmax[bad].astype(np.float32))
                fvs.append(v[bad].astype(np.float32))
                fr.append(r[bad].astype(np.float32))
                fsurf.append(touches_top[bad])
            del T, P
            if (s0 // a.chunk) % 20 == 0:
                print(f"    ... {s1:,}/{nt:,}", flush=True)

    print()
    out = {"mesh": a.mesh, "tets": int(nt), "parent_tets": int(a.parent_tets),
           "gates": {}}
    for g in gates:
        print(f"  gate {g:.4f}  (= {g*3/4:.3f} Hz at p3, {g*5/4:.3f} Hz at p5)")
        tot = 0
        for b in blocks:
            d = st[g][b]
            if d["n"] == 0:
                continue
            tot += d["nfail"]
            pct = 100.0 * d["nfail"] / d["n"]
            print(f"    {b:<7} {d['nfail']:>10,} of {d['n']:>12,}  ({pct:6.3f} %)   "
                  f"worst {d['worst']:.4f}")
        print(f"    {'TOTAL':<7} {tot:>10,} of {nt:>12,}  "
              f"({100.0*tot/nt:6.3f} %)\n")
        out["gates"][f"{g}"] = {b: {k: (None if not np.isfinite(v) else v)
                                    for k, v in st[g][b].items()} for b in blocks}

    if fb:
        B = np.concatenate(fb)
        IDX = np.concatenate(fx)
        DX = np.concatenate(fdx)
        VS = np.concatenate(fvs)
        R = np.concatenate(fr)
        SURF = np.concatenate(fsurf)
        z = B[:, 2]
        print(f"  failures at the strictest gate {gmax}: {len(B):,}")
        print(f"    barycentre depth   min {z.min():>10,.0f}  p10 {np.percentile(z,10):>9,.0f}  "
              f"med {np.median(z):>9,.0f}  max {z.max():>9,.0f}")
        print(f"    rest on the free surface (a vertex at z=0): "
              f"{int(SURF.sum()):,} ({100.0*SURF.mean():.1f} %)")
        print(f"    dx  med {np.median(DX):,.0f} m   Vs med {np.median(VS):,.0f} m/s")
        print(f"    dx needed (= Vs/gate) med {np.median(VS/gmax):,.0f} m   "
              f"refine factor med {np.median(DX*gmax/VS):.2f}x  p99 {np.percentile(DX*gmax/VS,99):.2f}x")
        print("\n    depth histogram of failures:")
        edges = [0, -125, -250, -500, -1000, -2000, -4000, -8000, -40001]
        for i in range(len(edges) - 1):
            hi, lo = edges[i], edges[i + 1]
            m = (z <= hi) & (z > lo)
            if m.any():
                print(f"      {hi:>7,} .. {lo:>8,} m : {int(m.sum()):>9,}  "
                      f"({100.0*m.mean():5.1f} %)  worst {R[m].min():.4f}")
        out["fail_depth_median"] = float(np.median(z))
        out["fail_on_surface"] = int(SURF.sum())
        if a.dump:
            np.savez_compressed(a.dump, bary=B, idx=IDX, dx=DX, vs=VS, ratio=R,
                                on_surface=SURF, gate=gmax)
            print(f"\n    dumped -> {a.dump}")

    if surf_ratio:
        SR = np.concatenate(surf_ratio)
        SD = np.concatenate(surf_dx)
        print(f"\n  free-surface-resting cells: {len(SR):,}")
        print("    |barycentre depth| / max edge  -- the size field's key constant")
        print(f"      p1 {np.percentile(SR,1):.3f}  p10 {np.percentile(SR,10):.3f}  "
              f"med {np.median(SR):.3f}  p90 {np.percentile(SR,90):.3f}  "
              f"p99 {np.percentile(SR,99):.3f}")
        print(f"    their max edge: med {np.median(SD):,.0f} m  "
              f"p90 {np.percentile(SD,90):,.0f} m  max {SD.max():,.0f} m")
        out["surface_bary_ratio_median"] = float(np.median(SR))
        out["surface_bary_ratio_p10"] = float(np.percentile(SR, 10))

    if hasattr(vs, "report"):
        vs.report()

    print("\n" + json.dumps({k: v for k, v in out.items() if k != "gates"}, indent=None))
    jp = Path(a.mesh).with_suffix("").name + (".muscal" if a.muscal else "") + ".census.json"
    Path(jp).write_text(json.dumps(out, indent=2))
    print(f"json -> {Path(jp).resolve()}")


if __name__ == "__main__":
    main()
