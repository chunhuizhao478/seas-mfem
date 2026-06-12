"""merged_soup_to_vtu.py — ParaView views of the (clipped) merged soup.

Writes, from the Phase 4 merged STL + markers JSON:
  <base>_all.vtu        every triangle, cell data: marker, q (triangle quality)
  <base>_faults.vtu     fault-marked triangles only (marker, q)
  <base>_boundary.vtu   boundary-marked triangles only (marker, q)
  <base>_fault_<k>_<shortname>.vtu   one per fault patch (q)

ParaView recipe: open *_all.vtu, color by `marker` (faults 1..N, boundary
100+); or open *_faults.vtu and color by `q` to inspect triangle quality.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import meshio
import numpy as np

import check_mesh_quality as cmq


def write_vtu(path: Path, pts: np.ndarray, tris: np.ndarray,
              cell_data: dict) -> None:
    used = np.unique(tris.ravel())
    remap = -np.ones(pts.shape[0], dtype=np.int64)
    remap[used] = np.arange(used.size)
    mesh = meshio.Mesh(points=pts[used],
                       cells=[("triangle", remap[tris])],
                       cell_data={k: [v] for k, v in cell_data.items()})
    meshio.write(str(path), mesh)
    print(f"  wrote {path.name}: {len(tris):,} tris")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--merged", type=Path, required=True)
    ap.add_argument("--markers", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)

    mk = json.loads(args.markers.read_text())
    markers = np.asarray(mk["markers"], dtype=np.int64)
    n_faults = int(mk["n_input_faults"])
    box_marker = int(mk["box_marker"])
    mesh = meshio.read(str(args.merged))
    pts = np.asarray(mesh.points, dtype=np.float64)
    tris = np.asarray(mesh.cells_dict["triangle"], dtype=np.int64)
    if len(tris) != len(markers):
        print(f"error: {len(tris)} tris vs {len(markers)} markers",
              file=sys.stderr)
        return 1

    _e, _a, q = cmq._tri_metrics(tris, pts)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    base = args.merged.stem

    write_vtu(args.out_dir / f"{base}_all.vtu", pts, tris,
              {"marker": markers.astype(np.int32), "q": q})
    fsel = (markers >= 1) & (markers <= n_faults)
    write_vtu(args.out_dir / f"{base}_faults.vtu", pts, tris[fsel],
              {"marker": markers[fsel].astype(np.int32), "q": q[fsel]})
    bsel = markers >= box_marker
    write_vtu(args.out_dir / f"{base}_boundary.vtu", pts, tris[bsel],
              {"marker": markers[bsel].astype(np.int32), "q": q[bsel]})
    for fm, bn in enumerate(mk["fault_basenames"], start=1):
        sel = markers == fm
        if not sel.any():
            continue
        short = bn.replace("fault_SAFS-SAFZ-", "").replace("_ext", "")[:40]
        write_vtu(args.out_dir / f"{base}_fault_{fm}_{short}.vtu",
                  pts, tris[sel], {"q": q[sel]})
    return 0


if __name__ == "__main__":
    sys.exit(main())
