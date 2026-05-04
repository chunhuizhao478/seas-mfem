"""orient_fault_surface.py — General fault winding-fixup pass.

Reads a .msh, BFS-propagates fault-triangle winding within each
tag-100 connected component starting from a deterministic seed
triangle.  Writes the corrected .msh.  Idempotent.

This implements R-003 of the 2026-05-02 SAFS topology review: every
mesh-mutating stage in the SAFS pipeline (HXT, mmg3d_post_pass,
mmg3d_local_patch, local_cavity_retet) must be followed by this
canonical winding-normalizer so downstream consumers (paraview_output,
DG flux assembly) see a consistently-oriented fault surface.

The fault surface (tag=100) is unsigned by construction in the
upstream meshers; HXT, mmg3d, and CGAL autorefine all emit
inconsistently-wound triangulations.  Without this pass, ParaView
back-face culling produces visible holes and DG flux signs are
ambiguous.
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import meshio
import numpy as np

FAULT_TAG = 100


def _orient_fault(tris: np.ndarray, ttags: np.ndarray, points: np.ndarray
                   ) -> tuple[np.ndarray, dict]:
    """Return new tris with consistent fault-tri winding + a stats dict.

    Box (tag != 100) tris are passed through unchanged.

    BFS within each connected component; for every propagation edge,
    if the two incident fault tris traverse it in the SAME direction,
    the second tri's winding is flipped (swap last two vertices).

    Edge classification (count = number of fault tris at the edge):
      count==1: fault perimeter — no propagation
      count==2: standard manifold edge — propagate to the unique
                neighbour
      count>=3: non-manifold edge (fault-fault intersection in
                branching SAFS geometry).  Partition the incident
                tris into DIHEDRAL PAIRS — two tris are paired iff
                their third-vertex bearings (perpendicular to the
                edge, projected) are anti-parallel (dot product < 0).
                Each dihedral pair is a "manifold continuation" of
                the same fault sheet across this edge; propagate
                orientation WITHIN each pair only.  Tris with no
                anti-parallel partner are "terminators" (e.g. a
                fault that ends at another fault's interior) and
                receive no propagation through this edge.

    This dihedral-aware rule subsumes the old manifold-only logic
    (count==2 always pairs cleanly with dot ≈ -1) and is essential
    for count==3 (T-junction: 1 pair + 1 terminator) and count==4
    (X-junction: 2 pairs).  Mirrored in
    validate_msh.check_13_fault_orientation.
    """
    fault_idx = np.where(ttags == FAULT_TAG)[0]
    if fault_idx.size == 0:
        return tris.copy(), {"n_fault_tris": 0, "n_flipped": 0,
                              "n_components": 0,
                              "n_nonmanifold_edges": 0}
    fault_tri = tris[fault_idx].copy()
    edge2tri: dict[frozenset, list[int]] = defaultdict(list)
    for ti in range(fault_tri.shape[0]):
        a, b, c = (int(fault_tri[ti, 0]), int(fault_tri[ti, 1]),
                    int(fault_tri[ti, 2]))
        for u, v in ((a, b), (b, c), (c, a)):
            edge2tri[frozenset((u, v))].append(ti)
    n_nonmanifold = sum(1 for v in edge2tri.values() if len(v) > 2)

    # Precompute dihedral partner table for non-manifold edges.
    # partner[(ti, edge)] = nti for the dihedral-paired neighbour
    # through `edge` (or absent if `ti` is a terminator).
    partner: dict[tuple[int, frozenset], int] = {}
    for e, tilist in edge2tri.items():
        if len(tilist) <= 2:
            continue
        uv = list(e)
        u_, v_ = uv[0], uv[1]
        edge_vec = points[v_] - points[u_]
        edge_len = float(np.linalg.norm(edge_vec))
        if edge_len < 1e-12:
            continue
        edge_unit = edge_vec / edge_len
        bearings: list[np.ndarray | None] = []
        for ti in tilist:
            a, b, c = (int(fault_tri[ti, 0]),
                        int(fault_tri[ti, 1]),
                        int(fault_tri[ti, 2]))
            third = (a if a not in (u_, v_)
                     else (b if b not in (u_, v_) else c))
            r = points[third] - points[u_]
            rp = r - float(np.dot(r, edge_unit)) * edge_unit
            n = float(np.linalg.norm(rp))
            bearings.append(rp / n if n > 1e-12 else None)
        used = [False] * len(tilist)
        for i in range(len(tilist)):
            if used[i] or bearings[i] is None:
                continue
            best_j, best_d = -1, 1.0
            for j in range(i + 1, len(tilist)):
                if used[j] or bearings[j] is None:
                    continue
                d = float(np.dot(bearings[i], bearings[j]))
                if d < best_d:
                    best_d = d
                    best_j = j
            # Require strongly anti-parallel bearings (dihedral
            # close to 180°) to count as a manifold pair.  -0.5
            # is cos(120°): a generous threshold that catches
            # smooth fault-bend geometry but excludes ~90°
            # X-/T-junction terminators.
            if best_j < 0 or best_d > -0.5:
                used[i] = True
                continue
            partner[(tilist[i], e)] = tilist[best_j]
            partner[(tilist[best_j], e)] = tilist[i]
            used[i] = True
            used[best_j] = True

    visited = np.zeros(fault_tri.shape[0], dtype=bool)
    n_flipped = 0
    n_components = 0
    for start in range(fault_tri.shape[0]):
        if visited[start]:
            continue
        n_components += 1
        visited[start] = True
        queue = [start]
        while queue:
            ti = queue.pop()
            a, b, c = (int(fault_tri[ti, 0]), int(fault_tri[ti, 1]),
                        int(fault_tri[ti, 2]))
            for u, v in ((a, b), (b, c), (c, a)):
                e = frozenset((u, v))
                e_tris = edge2tri[e]
                cnt = len(e_tris)
                if cnt == 1:
                    continue
                if cnt == 2:
                    nti_candidates = [n for n in e_tris if n != ti]
                else:
                    p = partner.get((ti, e))
                    if p is None:
                        # ti is a terminator at this non-manifold edge
                        continue
                    nti_candidates = [p]
                for nti in nti_candidates:
                    if nti == ti or visited[nti]:
                        continue
                    na, nb, nc = (int(fault_tri[nti, 0]),
                                    int(fault_tri[nti, 1]),
                                    int(fault_tri[nti, 2]))
                    same_dir = ((na == u and nb == v)
                                 or (nb == u and nc == v)
                                 or (nc == u and na == v))
                    if same_dir:
                        # FLIP: swap last two vertices
                        fault_tri[nti, 1], fault_tri[nti, 2] = (
                            fault_tri[nti, 2], fault_tri[nti, 1])
                        n_flipped += 1
                    visited[nti] = True
                    queue.append(nti)
    out = tris.copy()
    out[fault_idx] = fault_tri
    return out, {"n_fault_tris": int(fault_tri.shape[0]),
                 "n_flipped": int(n_flipped),
                 "n_components": int(n_components),
                 "n_nonmanifold_edges": int(n_nonmanifold)}


def _read_mesh(in_msh: Path):
    """Return (points, tris, ttags, tets, tetags) from a gmsh .msh.

    Concatenates ALL triangle/tetra cell-blocks (gmsh22 may split by
    physical tag) so the orient pass operates on a single flat
    array.  Tags are pulled from cell_data['gmsh:physical'] to
    preserve ordering.
    """
    m = meshio.read(in_msh)
    points = np.asarray(m.points, dtype=np.float64)
    gp = m.cell_data_dict.get("gmsh:physical", {})

    tri_blocks = [cb for cb in m.cells if cb.type == "triangle"]
    if tri_blocks:
        tri_data: list[np.ndarray] = []
        tri_tags: list[np.ndarray] = []
        all_tags = gp.get("triangle")
        cursor = 0
        for cb in tri_blocks:
            n = len(cb.data)
            tri_data.append(np.asarray(cb.data, dtype=np.int64))
            if all_tags is not None:
                tri_tags.append(np.asarray(all_tags[cursor:cursor + n],
                                            dtype=np.int32))
                cursor += n
            else:
                tri_tags.append(np.full(n, -1, dtype=np.int32))
        tris = np.concatenate(tri_data, axis=0)
        ttags = np.concatenate(tri_tags, axis=0)
    else:
        tris = np.empty((0, 3), dtype=np.int64)
        ttags = np.empty(0, dtype=np.int32)

    tet_blocks = [cb for cb in m.cells if cb.type == "tetra"]
    if tet_blocks:
        tet_data: list[np.ndarray] = []
        tet_tags: list[np.ndarray] = []
        all_tags = gp.get("tetra")
        cursor = 0
        for cb in tet_blocks:
            n = len(cb.data)
            tet_data.append(np.asarray(cb.data, dtype=np.int64))
            if all_tags is not None:
                tet_tags.append(np.asarray(all_tags[cursor:cursor + n],
                                            dtype=np.int32))
                cursor += n
            else:
                tet_tags.append(np.full(n, -1, dtype=np.int32))
        tets = np.concatenate(tet_data, axis=0)
        tetags = np.concatenate(tet_tags, axis=0)
    else:
        tets = np.empty((0, 4), dtype=np.int64)
        tetags = np.empty(0, dtype=np.int32)

    return points, tris, ttags, tets, tetags


def _write_mesh(out_msh: Path, points, tris, ttags, tets, tetags) -> None:
    """Write a gmsh .msh v2.2 ASCII via meshio, preserving the same
    cell layout as the upstream SAFS writers (`mmg3d_local_patch._write_msh`)."""
    out_msh.parent.mkdir(parents=True, exist_ok=True)
    mesh = meshio.Mesh(
        points=points,
        cells=[("triangle", tris), ("tetra", tets)],
        cell_data={"gmsh:physical": [ttags, tetags],
                   "gmsh:geometrical": [ttags, tetags]},
    )
    meshio.write(out_msh, mesh, file_format="gmsh22", binary=False)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Re-orient fault tri winding consistently per "
                    "connected component (R-003).  Idempotent.")
    p.add_argument("--in-msh", required=True, type=Path)
    p.add_argument("--out-msh", required=True, type=Path)
    args = p.parse_args(argv)
    points, tris, ttags, tets, tetags = _read_mesh(args.in_msh)
    new_tris, stats = _orient_fault(tris, ttags, points)
    print(f"  orient_fault_surface: components={stats['n_components']} "
          f"flipped={stats['n_flipped']} "
          f"nonmanifold={stats['n_nonmanifold_edges']}",
          file=sys.stderr)
    _write_mesh(args.out_msh, points, new_tris, ttags, tets, tetags)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
