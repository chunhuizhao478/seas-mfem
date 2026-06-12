"""report_fault_quality.py — per-fault surface quality table (Phase 3
acceptance + the end-of-Phase-4 fault quality report).

Two modes:
  --corefined-dir DIR    stats per fault_*_corefined.stl, restricted to
                         below-DEM / above-bottom triangles (the Phase 3
                         gate population: pre-clip, scaffolding excluded);
  --merged STL --markers JSON
                         stats per fault marker of the (clipped) merged
                         soup — the final fault surfaces.

Gates (Phase 3 / G1f): edge median in [425, 600] m, p99 <= 750 m,
min >= 100 m, tri q_min >= 0.30 (faults); report-only for others.

Metric: q = 4*sqrt(3)*A / sum(L^2)  (identical to check_mesh_quality.py).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import meshio
import numpy as np

import check_mesh_quality as cmq


def tri_stats(pts: np.ndarray, tris: np.ndarray) -> dict:
    tri_e, _a, tri_q = cmq._tri_metrics(tris, pts)
    e = tri_e.ravel()
    return {
        "n_tri": int(tris.shape[0]),
        "edge_min_m": float(e.min()),
        "edge_med_m": float(np.median(e)),
        "edge_p99_m": float(np.percentile(e, 99)),
        "edge_max_m": float(e.max()),
        "q_min": float(tri_q.min()),
        "q_med": float(np.median(tri_q)),
        "q_mean": float(tri_q.mean()),
        "n_q_below_0.30": int((tri_q < 0.30).sum()),
        "n_edges_below_100m": int((e < 100.0).sum()),
    }


def gate_eval(s: dict) -> dict:
    return {
        "median_in_[425,600]": bool(425.0 <= s["edge_med_m"] <= 600.0),
        "p99_le_750": bool(s["edge_p99_m"] <= 750.0),
        "min_ge_100": bool(s["edge_min_m"] >= 100.0),
        "q_min_ge_0.30": bool(s["q_min"] >= 0.30),
    }


def below_dem_mask(pts: np.ndarray, tris: np.ndarray,
                   dem_pts: np.ndarray, dem_tris: np.ndarray,
                   z_bottom: float) -> np.ndarray:
    from extend_fault_borders import dem_z_lookup
    z_at = dem_z_lookup(dem_pts, dem_tris)
    cen = pts[tris].mean(axis=1)
    zd = z_at(cen[:, :2])
    below_dem = cen[:, 2] < np.nan_to_num(zd, nan=np.inf)
    above_bot = cen[:, 2] > z_bottom + 1e-6
    return below_dem & above_bot


def fmt_table(rows: dict) -> str:
    hdr = ("| fault | n_tri | edge min | edge med | edge p99 | edge max | "
           "q_min | q_med | q<0.3 | gates |")
    sep = "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|"
    lines = [hdr, sep]
    for name, (s, g) in rows.items():
        gates = " ".join(("PASS" if v else "FAIL") + ":" + k.split("_")[0]
                         for k, v in g.items()) if g else "-"
        lines.append(
            f"| {name[:46]} | {s['n_tri']:,} | {s['edge_min_m']:.1f} | "
            f"{s['edge_med_m']:.1f} | {s['edge_p99_m']:.1f} | "
            f"{s['edge_max_m']:.1f} | {s['q_min']:.3f} | {s['q_med']:.3f} | "
            f"{s['n_q_below_0.30']} | {gates} |")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corefined-dir", type=Path, default=None)
    ap.add_argument("--extracted-dir", type=Path, default=None,
                    help="for dem.stl + bottom plane (below-DEM filter)")
    ap.add_argument("--merged", type=Path, default=None)
    ap.add_argument("--markers", type=Path, default=None)
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args(argv)

    out = {}
    rows = {}

    if args.corefined_dir:
        dem_pts = dem_tris = None
        z_bottom = -np.inf
        if args.extracted_dir:
            dm = meshio.read(str(args.extracted_dir / "dem.stl"))
            dem_pts = np.asarray(dm.points, dtype=np.float64)
            dem_tris = np.asarray(dm.cells_dict["triangle"], dtype=np.int64)
            eman = json.loads(
                (args.extracted_dir / "manifest.json").read_text())
            bb = next(m for m in eman["meshes"]
                      if m["basename"] == "bottom")["bbox"]
            z_bottom = bb["z_lo"]
        for stl in sorted(args.corefined_dir.glob("fault_*_corefined.stl")):
            m = meshio.read(str(stl))
            pts = np.asarray(m.points, dtype=np.float64)
            tris = np.asarray(m.cells_dict["triangle"], dtype=np.int64)
            label = stl.stem.replace("_ext_corefined", "")
            if dem_pts is not None:
                mask = below_dem_mask(pts, tris, dem_pts, dem_tris, z_bottom)
                tris_f = tris[mask]
                n_scaffold = int((~mask).sum())
            else:
                tris_f, n_scaffold = tris, 0
            s = tri_stats(pts, tris_f)
            s["n_scaffold_excluded"] = n_scaffold
            g = gate_eval(s)
            rows[label] = (s, g)
            out[label] = {"stats": s, "gates": g}
        # non-fault participants, report-only
        for stl in sorted(args.corefined_dir.glob("*_corefined.stl")):
            if stl.name.startswith("fault_"):
                continue
            m = meshio.read(str(stl))
            s = tri_stats(np.asarray(m.points, dtype=np.float64),
                          np.asarray(m.cells_dict["triangle"],
                                     dtype=np.int64))
            rows[stl.stem] = (s, None)
            out[stl.stem] = {"stats": s}

    if args.merged and args.markers:
        mk = json.loads(args.markers.read_text())
        markers = np.asarray(mk["markers"], dtype=np.int64)
        m = meshio.read(str(args.merged))
        pts = np.asarray(m.points, dtype=np.float64)
        tris = np.asarray(m.cells_dict["triangle"], dtype=np.int64)
        for fm, bn in enumerate(mk["fault_basenames"], start=1):
            sel = tris[markers == fm]
            if not len(sel):
                continue
            s = tri_stats(pts, sel)
            g = gate_eval(s)
            label = f"[final] {bn.replace('_ext', '')}"
            rows[label] = (s, g)
            out[label] = {"stats": s, "gates": g}

    print(fmt_table(rows))
    if args.json:
        with open(args.json, "w") as f:
            json.dump(out, f, indent=1, sort_keys=True)
            f.write("\n")
        print(f"\nwrote {args.json}")
    all_pass = all(all(g.values()) for s, g in rows.values() if g)
    print(f"\nfault gates: {'ALL PASS' if all_pass else 'FAILURES PRESENT'}")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
