"""Per-fault provenance side-car for SAFS unified-fault meshes.

Implements PLAN_domain.md Phase 3.

Reads the .msh, the per-fault STLs (in the local frame, written by
ts_to_stl.py), and assigns every tag-100 fault triangle to the CFM source
fault whose pre-fragment STL has the closest matching triangle centroid
(KDTree).  Outputs output/fault_provenance.json with per-fault triangle
indices.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Iterable

import numpy as np

from safs_origin import FAULT_SHORT_NAMES, read_transform_json

FAULT_TAG = 100


def _read_stl_centroids(path: Path) -> np.ndarray:
    """Read an ASCII STL and return triangle centroids (N, 3) plus the
    triangle vertices themselves (for area-weighted sums in caller).

    Only ASCII STL is expected — ts_to_stl.py writes ASCII exclusively.
    """
    coords: list[tuple[float, float, float]] = []
    with open(path, "r") as fh:
        for line in fh:
            tok = line.strip().split()
            if len(tok) >= 4 and tok[0] == "vertex":
                coords.append((float(tok[1]), float(tok[2]), float(tok[3])))
    if len(coords) % 3 != 0:
        raise ValueError(
            f"{path}: number of vertex lines ({len(coords)}) is not a "
            f"multiple of 3"
        )
    arr = np.asarray(coords, dtype=np.float64).reshape(-1, 3, 3)
    centroids = arr.mean(axis=1)
    return centroids, arr


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def _msh_fault_triangles(msh_path: Path):
    """Read the .msh and return (centroids, triangle_indices_in_msh).

    `triangle_indices_in_msh` is the 0-based index into the *flattened*
    triangle list that matches the order in which meshio yields tag-100
    triangles (callers should not assume any particular global mesh order
    beyond meshio's internal layout).
    """
    import meshio
    m = meshio.read(msh_path)
    gp = m.cell_data_dict.get("gmsh:physical", {})
    if "triangle" not in gp:
        raise ValueError(f"{msh_path} has no triangle elements")

    tris_per_block: list[np.ndarray] = []
    tags_per_block: list[np.ndarray] = []
    for cb, tg in zip(m.cells, gp["triangle"] if "triangle" in gp else []):
        # we'll loop over all cell blocks below
        pass
    # Iterate cell blocks and accumulate triangle vertex indices + tags.
    flat_indices: list[int] = []
    tri_centroids: list[np.ndarray] = []

    block_offset = 0  # not used outside; we just need per-tag triangle list
    # meshio: m.cells is a list of CellBlocks; gp['triangle'] is a flat
    # array aligned with the union of triangle blocks.
    tri_blocks = [(i, cb) for i, cb in enumerate(m.cells)
                  if cb.type == "triangle"]
    if not tri_blocks:
        raise ValueError("no triangle cell blocks")
    tri_tags = gp["triangle"]
    cursor = 0
    flat_idx_global = 0
    for _, cb in tri_blocks:
        n = len(cb.data)
        these_tags = tri_tags[cursor:cursor + n]
        cursor += n
        for j in range(n):
            tag = int(these_tags[j])
            if tag == FAULT_TAG:
                v = cb.data[j]
                centroid = m.points[v].mean(axis=0)
                tri_centroids.append(centroid)
                flat_indices.append(flat_idx_global)
            flat_idx_global += 1
    if not tri_centroids:
        raise ValueError(
            f"{msh_path}: no triangles with physical tag {FAULT_TAG} "
            f"(fault); check the mesh"
        )
    return np.asarray(tri_centroids, dtype=np.float64), \
           np.asarray(flat_indices, dtype=np.int64)


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Assign tag-100 fault triangles to CFM source faults "
                    "via KDTree nearest-centroid matching, write a "
                    "schema-versioned JSON.",
    )
    parser.add_argument("--msh", required=True, type=Path)
    parser.add_argument("--stl-dir", required=True, type=Path)
    parser.add_argument("--transform-json", default=None, type=Path,
                        help="Used only as cross-reference; recorded in "
                             "the output for downstream tools.  "
                             "Defaults to <msh-parent>/../transform.json.")
    parser.add_argument("--include-fault", action="append", default=None,
                        metavar="SHORT_NAME",
                        help="Restrict to this short name; may be "
                             "repeated.  Default: every STL in --stl-dir "
                             "that matches a known short name.")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--ambiguity-distance", type=float, default=None,
                        metavar="METRES",
                        help="Mesh triangles farther than this from any "
                             "STL centroid are flagged 'ambiguous'.  "
                             "Default: 2 * res_f from sizing.json if "
                             "available, else 2000 m.")
    args = parser.parse_args(argv)

    # Resolve transform.json default
    transform_json = args.transform_json or \
        (args.msh.parent.parent / "transform.json")
    if transform_json.exists():
        try:
            tx = read_transform_json(transform_json)
        except Exception as e:
            sys.stderr.write(
                f"warning: failed to load {transform_json}: {e}\n"
            )
            tx = None
    else:
        tx = None

    # Default ambiguity distance: read from sizing.json
    if args.ambiguity_distance is None:
        sizing_path = args.msh.parent / "sizing.json"
        if sizing_path.exists():
            sz = json.loads(sizing_path.read_text())
            ambiguity_dist = 2.0 * float(sz.get("res_f_m", 1000.0))
        else:
            ambiguity_dist = 2000.0
    else:
        ambiguity_dist = float(args.ambiguity_distance)

    # 1. Discover STLs
    candidate_shorts = list(FAULT_SHORT_NAMES.values())
    stls: dict[str, Path] = {}
    for short in candidate_shorts:
        p = args.stl_dir / f"{short}.stl"
        if p.exists():
            stls[short] = p
    if args.include_fault:
        keep = set(args.include_fault)
        unknown = sorted(keep - set(stls))
        if unknown:
            raise SystemExit(
                f"--include-fault asked for {unknown} but those STLs are "
                f"not in {args.stl_dir}"
            )
        stls = {k: v for k, v in stls.items() if k in keep}
    if not stls:
        raise SystemExit(f"no STLs found in {args.stl_dir}")

    # 2. Build KDTree of (centroid, label_idx)
    from scipy.spatial import cKDTree
    short_list = sorted(stls)
    all_centroids: list[np.ndarray] = []
    all_labels: list[int] = []
    cfm_tri_per_fault: dict[str, np.ndarray] = {}
    for label_idx, short in enumerate(short_list):
        c, tri_v = _read_stl_centroids(stls[short])
        all_centroids.append(c)
        all_labels.extend([label_idx] * c.shape[0])
        cfm_tri_per_fault[short] = c
    centroid_arr = np.vstack(all_centroids)
    label_arr = np.asarray(all_labels, dtype=np.int64)
    tree = cKDTree(centroid_arr)

    # 3. Read .msh tag-100 triangles
    msh_centroids, msh_indices = _msh_fault_triangles(args.msh)

    # 4. Nearest-centroid query
    dists, idxs = tree.query(msh_centroids, k=1)
    assigned_label = label_arr[idxs]
    n_ambig = int(np.sum(dists > ambiguity_dist))

    # 5. Build per-fault triangle index lists.
    fault_entries: dict[str, dict] = {}
    for label_idx, short in enumerate(short_list):
        mask = assigned_label == label_idx
        msh_tri_idxs = msh_indices[mask].tolist()
        fault_entries[short] = {
            "cfm_id": next(
                (cid for cid, sn in FAULT_SHORT_NAMES.items() if sn == short),
                None,
            ),
            "stl_file": str(stls[short]),
            "stl_sha256": _sha256(stls[short]),
            "n_triangles_in_msh": len(msh_tri_idxs),
            "triangle_indices_in_msh": msh_tri_idxs,
        }
        if not msh_tri_idxs:
            fault_entries[short]["omitted_reason"] = (
                "no tag-100 triangle in the mesh was nearest to this STL "
                "(possibly fragmented away or filtered out by "
                "--include-fault)"
            )

    # 6. Sanity: every triangle assigned exactly once
    total_assigned = sum(e["n_triangles_in_msh"] for e in fault_entries.values())
    if total_assigned != msh_centroids.shape[0]:
        raise RuntimeError(
            f"provenance partition leak: assigned {total_assigned}, but "
            f"mesh has {msh_centroids.shape[0]} tag-100 triangles"
        )

    out = {
        "schema_version": 1,
        "mesh_file": str(args.msh),
        "fault_tag": FAULT_TAG,
        "transform_json": str(transform_json) if tx is not None else None,
        "ambiguity_distance_m": ambiguity_dist,
        "n_ambiguous_triangles": n_ambig,
        "n_tag100_triangles_total": int(msh_centroids.shape[0]),
        "faults": fault_entries,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")
    print(f"write_fault_provenance OK: {args.out}")
    print(f"  total tag-100 triangles : {out['n_tag100_triangles_total']}")
    print(f"  ambiguous (>{ambiguity_dist:.0f} m): {n_ambig}")
    for short, e in fault_entries.items():
        print(f"  {short}: {e['n_triangles_in_msh']} triangles")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
