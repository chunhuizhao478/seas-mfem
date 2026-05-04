#!/usr/bin/env python
"""Drop coplanar-overlapping triangle pairs from per-fault STLs.

This is the Step-6-unblocking fix for the autorefine R-501 failure
mode (REVIEW_autorefine_mode.md): when two source faults are
physically adjacent (e.g., the "coav_missioncreek" and
"sbmt_missioncreek" CFM branches that meet at the surface trace),
CGAL 6.1's autorefine emits OVERLAPPING triangulated patches of the
same physical surface region.  The per-fault STL output therefore
contains pairs of triangles that:

  - share 2 vertex IDs (i.e., share a common edge), and
  - lie on the same plane (|n_A · n_B| ≈ 1).

HXT's PLC recovery then rejects the input with:

    "Found two exactly self-intersecting facets
     (dihedral angle 0.00000E+00).  1st: [a, b, c] #7  2nd: [a, b, d] #7"

generate_safs_mesh.py's `_combine_stls` already drops triangles with
identical 3-vertex IDs (the trivial-duplicate case), but does NOT
drop coplanar-overlapping pairs whose third vertex differs.  This
script does that pass.

Algorithm:
  1. Read per-fault STLs.
  2. Build a UNIFIED snap-grid vertex map across all faults.
  3. For each fault's triangle, register it in an edge-to-triangle
     index keyed on the (snapped) edge.
  4. Walk the index: for each edge with > 1 incident triangle (across
     all faults), check pairwise normals; coplanar pairs (|n.n|
     > 1 - tol) are flagged.
  5. Drop the triangle that comes from the LATER fault in the
     `--include-fault` list (deterministic).  If both incident
     triangles come from the same fault, drop the one with the higher
     local index.
  6. Write filtered per-fault STLs.

Cross-fault provenance: the chosen-fault rule favours the
`--include-fault` order, so the user can prioritize a specific
fault by listing it first.

Usage:
  python dedup_coplanar_facets.py \\
      --in-stl-dir output/.../stl_conformal \\
      --out-stl-dir output/.../stl_dedup \\
      --include-fault A --include-fault B [...] \\
      [--snap-m 0.1] [--coplanar-tol 1e-3]
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
from typing import Iterable

import numpy as np

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
import break_fault_wedges as bfw  # noqa: E402


def _normal_unit(V: np.ndarray, t: tuple[int, int, int]) -> np.ndarray:
    """Unit normal of triangle t."""
    v0, v1, v2 = V[t[0]], V[t[1]], V[t[2]]
    n = np.cross(v1 - v0, v2 - v0)
    L = np.linalg.norm(n)
    if L == 0.0:
        return np.zeros(3)
    return n / L


def _build_global_index(
        per_fault: list[tuple[str, np.ndarray, np.ndarray, str]],
        snap_m: float
        ) -> tuple[np.ndarray,
                    list[tuple[str, list[tuple[int, int, int]]]]]:
    """Build a global vertex array (deduped by snap_key across all
    faults) and per-fault triangle lists with global vertex indices.

    Returns:
      global_V: (N_unique, 3) float64 — first-seen coords for each
                unique snap_key.
      tagged_tris: list of (fault_name, [(a, b, c), ...]) where
                   indices reference global_V.
    """
    snap_to_global: dict[tuple[int, int, int], int] = {}
    global_V_list: list[tuple[float, float, float]] = []
    tagged_tris: list[tuple[str, list[tuple[int, int, int]]]] = []

    for name, V, T, _ in per_fault:
        local_to_global: list[int] = []
        for vi in range(V.shape[0]):
            key = bfw._snap_key(V[vi], snap_m)
            if key in snap_to_global:
                local_to_global.append(snap_to_global[key])
            else:
                gi = len(global_V_list)
                snap_to_global[key] = gi
                global_V_list.append(tuple(V[vi]))
                local_to_global.append(gi)
        ftris: list[tuple[int, int, int]] = []
        for ti in range(T.shape[0]):
            a = local_to_global[int(T[ti, 0])]
            b = local_to_global[int(T[ti, 1])]
            c = local_to_global[int(T[ti, 2])]
            # Skip degenerate triangles (two coincident snap_keys).
            if a == b or b == c or a == c:
                continue
            ftris.append((a, b, c))
        tagged_tris.append((name, ftris))

    return (np.asarray(global_V_list, dtype=np.float64), tagged_tris)


def _same_side_in_plane(global_V: np.ndarray,
                          a: int, b: int, c: int, d: int,
                          n: np.ndarray
                          ) -> bool:
    """For two triangles (a, b, c) and (a, b, d) sharing edge (a, b)
    in the same plane (normal `n`), return True iff vertices c and
    d lie on the SAME side of line (a, b).

    "Same side" means the two triangles OVERLAP.  "Opposite side"
    means they tile the quad (a, b, c, d) — a legitimate
    triangulation.

    Computed via signed planar cross-product: side(p) = sign(dot(
    cross(b - a, p - a), n)).  Two points on the same side of the
    line in the plane have the same sign.
    """
    pa = global_V[a]
    pb = global_V[b]
    pc = global_V[c]
    pd = global_V[d]
    ab = pb - pa
    s_c = np.dot(np.cross(ab, pc - pa), n)
    s_d = np.dot(np.cross(ab, pd - pa), n)
    # Same-side iff both sides have the same sign.  Treat exact zero
    # as "on the line" — degenerate; treat as same-side
    # (defensively, since a third vertex exactly on the shared edge
    # is itself a degenerate triangle).
    return (s_c * s_d) > 0.0


def _find_coplanar_overlap_losers(
        global_V: np.ndarray,
        tagged_tris: list[tuple[str, list[tuple[int, int, int]]]],
        coplanar_tol: float,
        include_order: list[str]
        ) -> set[tuple[str, int]]:
    """R-004 fix: identify losers using cluster-based resolution.

    For each shared edge, cluster all triangle incidences whose
    third vertices are mutually on the same side of the shared edge
    AND whose normals are coplanar.  For each cluster of size K ≥ 2,
    keep ONE incidence (the one whose fault is earliest in
    `include_order`, breaking ties by lower local index) and mark
    the remaining K-1 as losers.

    The previous implementation enumerated unordered pairs and
    flagged each loser independently.  At a triple-junction with 3
    mutually-coplanar incidences, that approach produced 3 pairs
    each contributing one loss, dropping K-1 = 2 incidences.  But
    when the loser-selection rule is "later-fault loses,"
    pair (A, B) drops B; pair (A, C) drops C; pair (B, C) drops C
    (already dropped).  Net K-1 = 2 dropped.  Equivalent to cluster.

    The cluster-based formulation is clearer AND handles
    asymmetric cases (e.g., (A, C) coplanar, (A, B) coplanar, but
    (B, C) NOT coplanar within tol — graph-cluster vs pairwise
    matters).

    Returns the set of (fault_name, local_tri_idx) losers.
    `include_order` provides the deterministic priority for the
    keeper.
    """
    # Build edge -> [(fault_name, local_tri_idx, third_vertex_id, normal)]
    edge_to_incidence: dict[tuple[int, int],
                             list[tuple[str, int, int, np.ndarray]]] = {}
    for fname, ftris in tagged_tris:
        for li, (a, b, c) in enumerate(ftris):
            for u, v, w in ((a, b, c), (b, c, a), (c, a, b)):
                ek = (u, v) if u < v else (v, u)
                n = _normal_unit(global_V, (a, b, c))
                edge_to_incidence.setdefault(ek, []).append(
                    (fname, li, w, n))

    # R-007 perf fix: build a name → priority lookup once instead of
    # calling include_order.index() per pair.
    pos_index = {name: i for i, name in enumerate(include_order)}

    def _priority(member: tuple[str, int]) -> tuple[int, int]:
        # Tuple comparison: first by fault rank (earlier = lower),
        # then by local index (lower = earlier in input order).
        return (pos_index[member[0]], member[1])

    losers: set[tuple[str, int]] = set()
    for ek, incs in edge_to_incidence.items():
        if len(incs) < 2:
            continue
        a, b = ek

        # Cluster by mutual coplanarity + same-side at this edge.
        # Use a Union-Find over incidence indices.
        parent = list(range(len(incs)))

        def _find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def _union(x, y):
            rx, ry = _find(x), _find(y)
            if rx != ry:
                parent[rx] = ry

        for i in range(len(incs)):
            for j in range(i + 1, len(incs)):
                fA, lA, wA, nA = incs[i]
                fB, lB, wB, nB = incs[j]
                if (fA, lA) == (fB, lB):
                    continue
                d = float(abs(np.dot(nA, nB)))
                if d <= 1.0 - coplanar_tol:
                    continue
                if wA == wB:
                    # Same triangle (vertex-identity duplicate).
                    continue
                if not _same_side_in_plane(global_V,
                                            a, b, wA, wB, nA):
                    continue
                _union(i, j)

        # For each cluster of size ≥ 2, keep the highest-priority
        # member; mark the rest as losers.
        cluster_members: dict[int, list[int]] = {}
        for i in range(len(incs)):
            r = _find(i)
            cluster_members.setdefault(r, []).append(i)
        for r, members in cluster_members.items():
            if len(members) < 2:
                continue
            entries = [(incs[i][0], incs[i][1]) for i in members]
            keeper = min(entries, key=_priority)
            for entry in entries:
                if entry != keeper:
                    losers.add(entry)
    return losers


def _find_coplanar_overlapping_pairs(
        global_V: np.ndarray,
        tagged_tris: list[tuple[str, list[tuple[int, int, int]]]],
        coplanar_tol: float
        ) -> list[tuple[tuple[str, int], tuple[str, int]]]:
    """Legacy pairwise API kept for backward compatibility with
    existing unit tests.  Internally calls
    `_find_coplanar_overlap_losers` with a synthetic include_order
    derived from the tagged_tris fault names.

    The pair list is reconstructed from the loser set by re-walking
    the edge incidences and pairing each loser with the keeper of
    its cluster.  See R-004 fix.
    """
    include_order = [name for name, _ in tagged_tris]
    losers = _find_coplanar_overlap_losers(
        global_V, tagged_tris, coplanar_tol, include_order)
    # For backward compatibility we report each loser paired with
    # the cluster keeper.  This requires re-detecting clusters,
    # which we do by repeating the structure.
    pairs: list[tuple[tuple[str, int], tuple[str, int]]] = []
    edge_to_incidence: dict[tuple[int, int],
                             list[tuple[str, int, int, np.ndarray]]] = {}
    for fname, ftris in tagged_tris:
        for li, (a, b, c) in enumerate(ftris):
            for u, v, w in ((a, b, c), (b, c, a), (c, a, b)):
                ek = (u, v) if u < v else (v, u)
                n = _normal_unit(global_V, (a, b, c))
                edge_to_incidence.setdefault(ek, []).append(
                    (fname, li, w, n))
    pos_index = {name: i for i, name in enumerate(include_order)}
    for ek, incs in edge_to_incidence.items():
        if len(incs) < 2:
            continue
        a, b = ek
        parent = list(range(len(incs)))

        def _find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def _union(x, y):
            rx, ry = _find(x), _find(y)
            if rx != ry:
                parent[rx] = ry

        for i in range(len(incs)):
            for j in range(i + 1, len(incs)):
                fA, lA, wA, nA = incs[i]
                fB, lB, wB, nB = incs[j]
                if (fA, lA) == (fB, lB):
                    continue
                d = float(abs(np.dot(nA, nB)))
                if d <= 1.0 - coplanar_tol:
                    continue
                if wA == wB:
                    continue
                if not _same_side_in_plane(global_V,
                                            a, b, wA, wB, nA):
                    continue
                _union(i, j)

        cluster_members: dict[int, list[int]] = {}
        for i in range(len(incs)):
            r = _find(i)
            cluster_members.setdefault(r, []).append(i)
        for members in cluster_members.values():
            if len(members) < 2:
                continue
            entries = [(incs[i][0], incs[i][1]) for i in members]
            keeper = min(entries, key=lambda m: (pos_index[m[0]], m[1]))
            for entry in entries:
                if entry != keeper:
                    a_pair = (entry, keeper) if entry < keeper \
                             else (keeper, entry)
                    if a_pair not in pairs:
                        pairs.append(a_pair)
    return pairs


def _select_loser(pair: tuple[tuple[str, int], tuple[str, int]],
                   include_order: list[str]
                   ) -> tuple[str, int]:
    """For a coplanar overlap pair, return the (fault_name, idx)
    that should be DROPPED.  Rule: keep the triangle from the fault
    listed earlier in `include_order`; if both come from the same
    fault, drop the higher index (deterministic).
    """
    a, b = pair
    if a[0] == b[0]:
        # Same fault: drop higher local index.
        return a if a[1] > b[1] else b
    # Different faults: drop the one whose fault comes LATER in
    # include_order.
    pos_a = include_order.index(a[0])
    pos_b = include_order.index(b[0])
    return a if pos_a > pos_b else b


def dedup_coplanar(
        in_stl_dir: Path, out_stl_dir: Path,
        include_fault: list[str],
        snap_m: float = 0.1,
        coplanar_tol: float = 1e-3
        ) -> dict:
    """Run the dedup pass.  Returns a per-fault report dict."""
    if not in_stl_dir.is_dir():
        raise SystemExit(f"--in-stl-dir not found: {in_stl_dir}")
    if not include_fault:
        raise SystemExit("--include-fault must appear at least once")
    if snap_m <= 0:
        raise ValueError(f"snap_m must be > 0; got {snap_m}")
    if not (0.0 < coplanar_tol < 1.0):
        raise ValueError(
            f"coplanar_tol must be in (0, 1); got {coplanar_tol}")
    out_stl_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load per-fault STLs.
    per_fault: list[tuple[str, np.ndarray, np.ndarray, str]] = []
    for short in include_fault:
        path = in_stl_dir / f"{short}.stl"
        if not path.exists():
            raise SystemExit(f"missing input STL: {path}")
        V, T, solid = bfw._read_ascii_stl(path)
        per_fault.append((short, V, T, solid))
        print(f"  load {short}: V={V.shape[0]}, T={T.shape[0]}",
              file=sys.stderr)

    # 2. Build global indexing.
    global_V, tagged_tris = _build_global_index(per_fault, snap_m)
    print(f"  global unique vertices (snap={snap_m}): {global_V.shape[0]}",
          file=sys.stderr)

    # 3. Find coplanar overlapping LOSERS directly (R-004 fix).
    # The cluster-based detector resolves K-incidence overlaps
    # correctly: keep one, drop K-1.  No risk of dropping all K
    # at a triple-coplanar junction.
    loser_entries = _find_coplanar_overlap_losers(
        global_V, tagged_tris, coplanar_tol, include_fault)
    print(f"  coplanar overlap losers: {len(loser_entries)}",
          file=sys.stderr)

    # 4. Group losers by fault.
    losers: dict[str, set[int]] = {n: set() for n in include_fault}
    for (fname, li) in loser_entries:
        losers[fname].add(li)

    # 5. Write filtered per-fault STLs.  Translate global indices
    #    back to per-fault local indices for output.
    per_fault_report: dict[str, dict] = {}
    for (name, V, T, solid), (gname, ftris) in zip(per_fault, tagged_tris):
        assert name == gname
        drop_set = losers[name]
        # Build (V_filtered, T_filtered) using ORIGINAL local coords.
        # Map: which local triangles to keep?
        # Note: tagged_tris dropped degenerate triangles in
        # _build_global_index.  We need to keep them aligned.
        # Re-walk the original T and filter both degenerate AND
        # losers.
        kept_tris: list[tuple[int, int, int]] = []
        local_kept_count = 0
        local_li = 0  # Index into ftris (post-degenerate-filter).
        for ti in range(T.shape[0]):
            a, b, c = int(T[ti, 0]), int(T[ti, 1]), int(T[ti, 2])
            # Use local-V snap_keys to detect degenerate.
            ka = bfw._snap_key(V[a], snap_m)
            kb = bfw._snap_key(V[b], snap_m)
            kc = bfw._snap_key(V[c], snap_m)
            if ka == kb or kb == kc or ka == kc:
                # Degenerate.  Already excluded from ftris in
                # _build_global_index.  Drop here too.
                continue
            # This triangle has local_li in ftris.
            if local_li in drop_set:
                local_li += 1
                continue
            kept_tris.append((a, b, c))
            local_kept_count += 1
            local_li += 1
        T_kept = np.asarray(kept_tris, dtype=np.int64)
        out_path = out_stl_dir / f"{name}.stl"
        bfw._write_ascii_stl(out_path, V, T_kept,
                              solid or f"SAFS:{name}:dedup")
        per_fault_report[name] = {
            "T_in":  int(T.shape[0]),
            "T_out": int(T_kept.shape[0]),
            "n_dropped_coplanar": int(len(drop_set)),
        }
        print(f"  write {name}: T={T.shape[0]}->{T_kept.shape[0]} "
              f"(dropped {len(drop_set)} coplanar overlaps)",
              file=sys.stderr)

    return {
        "snap_m": snap_m,
        "coplanar_tol": coplanar_tol,
        # Total triangles dropped across all faults (= number of
        # coplanar overlap losers identified).  Replaces the old
        # `n_coplanar_pairs_detected` field whose semantics were
        # ambiguous under multi-incidence clusters.  Backward-compat
        # alias kept below.
        "n_coplanar_losers": len(loser_entries),
        "n_coplanar_pairs_detected": len(loser_entries),
        "per_fault": per_fault_report,
    }


def main(argv: Iterable[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Drop coplanar-overlapping triangles from "
                    "per-fault STLs.")
    p.add_argument("--in-stl-dir", required=True, type=Path)
    p.add_argument("--out-stl-dir", required=True, type=Path)
    p.add_argument("--include-fault", action="append", default=[],
                   metavar="SHORT_NAME")
    p.add_argument("--snap-m", type=float, default=0.1,
                   help="Coordinate snap precision for cross-fault "
                        "vertex matching.  Default 0.1 m.")
    p.add_argument("--coplanar-tol", type=float, default=1e-3,
                   help="|n_A · n_B| > 1 - coplanar_tol marks a pair "
                        "as coplanar.  Default 1e-3 (= ~2.5° angle "
                        "tolerance).")
    p.add_argument("--report-json", type=Path, default=None)
    args = p.parse_args(list(argv) if argv is not None else None)

    report = dedup_coplanar(
        in_stl_dir=args.in_stl_dir,
        out_stl_dir=args.out_stl_dir,
        include_fault=args.include_fault,
        snap_m=args.snap_m,
        coplanar_tol=args.coplanar_tol,
    )

    report_path = args.report_json or (
        args.out_stl_dir / "dedup_coplanar_report.json")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w") as fh:
        json.dump(report, fh, indent=2)
    print(f"  wrote {report_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
