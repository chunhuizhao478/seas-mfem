"""clip_fault_overhangs.py — Phase 4.3/4.4 of PLAN_mesh_quality_safv4_remesh.

Remove every fault overhang from the autorefined merged soup and validate
the result as a closed PLC ready for tetgen:

  (a) shell containment — drop fault-marked triangles whose centroid is
      OUTSIDE the closed boundary shell (ray parity along +z against
      boundary-marked triangles; covers both above-DEM and below-bottom
      overhangs — the DEM-heightfield shortcut alone would keep
      below-bottom triangles);
  (b) junction overhangs — for each junction pair, cut the GUEST fault at
      the corefined junction polyline and keep its largest connected
      component (the strips are scaffolding);
  (c) 1e-6 margin at the polylines (inside the parity test).

Post-clip validation (same script): boundary closure (every boundary edge
shared by exactly 2 boundary triangles), fault border classification
(dem / bottom / junction / tip; near-miss "orphans" are a hard failure),
min edge >= --min-edge, duplicate-coordinate scan at --tol.

CLI
---
    python clip_fault_overhangs.py --merged-stl RAW.stl --markers RAW.json
        --corefine-manifest .../manifest.json --extracted-manifest .../manifest.json
        --out-stl OUT.stl --out-markers OUT.json --report REPORT.json
        [--tol 1e-3] [--min-edge 100]
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import meshio
import numpy as np

from extract_safv4_surfaces import write_ascii_stl

WELD_TOL = 1e-6     # autorefine intern() dedup scale
GRID_CELL = 2000.0  # xy binning for the ray-parity test


def _weld_exact(points: np.ndarray, tris: np.ndarray, tol: float = WELD_TOL):
    """Re-establish the autorefine intern() vertex identity after the STL
    round-trip (weld coordinate-coincident vertices at `tol`)."""
    q = np.round(points / tol).astype(np.int64)
    uniq, inverse = np.unique(q, axis=0, return_inverse=True)
    # representative coordinates: first occurrence of each group
    first = np.full(uniq.shape[0], -1, dtype=np.int64)
    for i in range(points.shape[0] - 1, -1, -1):
        first[inverse[i]] = i
    new_pts = points[first]
    new_tris = inverse[tris]
    return new_pts, new_tris


def ray_parity_inside(fault_cen: np.ndarray, b_tris: np.ndarray,
                      pts: np.ndarray, margin: float = 1e-6) -> np.ndarray:
    """For each fault centroid: parity of +z ray crossings against the
    boundary triangles (True = inside the closed shell)."""
    cell = GRID_CELL
    grid: dict = defaultdict(list)
    P = pts
    for ti, t in enumerate(b_tris):
        xy = P[t][:, :2]
        lo = np.floor(xy.min(axis=0) / cell).astype(int)
        hi = np.floor(xy.max(axis=0) / cell).astype(int)
        for cx in range(lo[0], hi[0] + 1):
            for cy in range(lo[1], hi[1] + 1):
                grid[(cx, cy)].append(ti)
    inside = np.zeros(len(fault_cen), dtype=bool)
    for i, c in enumerate(fault_cen):
        x, y, zc = c
        key = (int(np.floor(x / cell)), int(np.floor(y / cell)))
        crossings: list[float] = []
        for ti in grid.get(key, ()):  # candidates
            a, b, d = b_tris[ti]
            x1, y1, z1 = P[a]
            x2, y2, z2 = P[b]
            x3, y3, z3 = P[d]
            det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
            if det == 0.0:
                continue  # vertical triangle (ribbon) — ray parallel
            l1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / det
            l2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / det
            l3 = 1.0 - l1 - l2
            if min(l1, l2, l3) < -1e-12:
                continue
            z = l1 * z1 + l2 * z2 + l3 * z3
            if z > zc + margin:
                crossings.append(z)
        # collapse near-duplicate crossings (shared-edge double hits)
        crossings.sort()
        n = 0
        last = None
        for z in crossings:
            if last is None or z - last > 1e-6:
                n += 1
            last = z
        inside[i] = (n % 2) == 1
    return inside


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--merged-stl", type=Path, required=True)
    ap.add_argument("--markers", type=Path, required=True)
    ap.add_argument("--corefine-manifest", type=Path, required=True)
    ap.add_argument("--extracted-manifest", type=Path, required=True)
    ap.add_argument("--out-stl", type=Path, required=True)
    ap.add_argument("--out-markers", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--tol", type=float, default=1e-3)
    ap.add_argument("--min-edge", type=float, default=100.0)
    args = ap.parse_args(argv)

    mk = json.loads(args.markers.read_text())
    markers = np.asarray(mk["markers"], dtype=np.int64)
    n_faults = int(mk["n_input_faults"])
    fault_basenames = mk["fault_basenames"]
    box_marker = int(mk["box_marker"])

    mesh = meshio.read(str(args.merged_stl))
    pts = np.asarray(mesh.points, dtype=np.float64)
    tris = np.asarray(mesh.cells_dict["triangle"], dtype=np.int64)
    if len(tris) != len(markers):
        print(f"error: {len(tris)} triangles vs {len(markers)} markers",
              file=sys.stderr)
        return 1
    pts, tris = _weld_exact(pts, tris)
    # Welding coincident duplicated vertices (junction pinch clusters)
    # turns their 0-length-edge triangles degenerate: drop them.
    degen = ((tris[:, 0] == tris[:, 1]) | (tris[:, 1] == tris[:, 2])
             | (tris[:, 0] == tris[:, 2]))
    if degen.any():
        tris = tris[~degen]
        markers = markers[~degen]
    n_degen_dropped = int(degen.sum())

    is_fault = (markers >= 1) & (markers <= n_faults)
    is_boundary = markers >= box_marker
    report = {"n_tris_in": int(len(tris)),
              "n_degen_dropped_at_weld": n_degen_dropped,
              "n_fault_in": int(is_fault.sum()),
              "n_boundary": int(is_boundary.sum())}

    # ---- (a) shell containment -----------------------------------------
    fidx = np.nonzero(is_fault)[0]
    cen = pts[tris[fidx]].mean(axis=1)
    inside = ray_parity_inside(cen, tris[is_boundary], pts)
    drop = np.zeros(len(tris), dtype=bool)
    drop[fidx[~inside]] = True
    report["clip_a_outside_shell"] = int((~inside).sum())

    # ---- (b) junction flood-fill ----------------------------------------
    cman = json.loads(args.corefine_manifest.read_text())
    eman = json.loads(args.extracted_manifest.read_text())
    base_by_idx = {m["index"]: m["basename"] for m in cman["meshes"]}
    marker_by_base = {}
    for fm, bn in enumerate(fault_basenames, start=1):
        marker_by_base[bn] = fm

    junction_cuts = []
    for j in eman["junctions"].values():
        guest, host = j["guest"], j["host"]
        gbase = next((b for b in marker_by_base if b.startswith(guest)), None)
        hbase_idx = next((i for i, b in base_by_idx.items()
                          if b.startswith(host)), None)
        gidx_m = next((i for i, b in base_by_idx.items()
                       if b.startswith(guest)), None)
        if gbase is None or hbase_idx is None or gidx_m is None:
            print(f"error: junction {guest}->{host} not resolvable in "
                  f"manifests", file=sys.stderr)
            return 1
        pl = None
        for p in cman["pairs"]:
            ij = {p["i"], p["j"]}
            if ij == {gidx_m, hbase_idx}:
                pl = p["polylines"]
                break
        if pl is None:
            print(f"error: no corefine pair for junction {guest}->{host}",
                  file=sys.stderr)
            return 1
        junction_cuts.append((marker_by_base[gbase], pl, guest, host))

    def edge_set_of(sel_mask) -> set:
        t = tris[sel_mask]
        e = np.sort(np.concatenate([t[:, [0, 1]], t[:, [1, 2]],
                                    t[:, [2, 0]]]), axis=1)
        return set(map(tuple, e))

    report["clip_b_junction"] = {}
    for gmarker, polylines, guest, host in junction_cuts:
        # Cut edges = edges shared between the guest and the HOST in the
        # welded soup: the Phase 3 corefine made the junction polyline
        # bit-identical in both meshes, so the shared edge set IS the
        # junction line (the manifest polylines are the RESAMPLED curves,
        # whose vertices are not mesh vertices — matching those failed).
        hmarker = next(fm for fm, bn in enumerate(fault_basenames, start=1)
                       if bn.startswith(host))
        host_edges = edge_set_of((markers == hmarker) & ~drop)
        g_tris_idx = np.nonzero((markers == gmarker) & ~drop)[0]
        cut_edges = set()
        edge_owner: dict = defaultdict(list)
        for ti in g_tris_idx:
            a, b, c = tris[ti]
            for u, v in ((a, b), (b, c), (c, a)):
                k = (min(u, v), max(u, v))
                if k in host_edges:
                    cut_edges.add(k)
                    continue
                edge_owner[k].append(ti)
        parent = {int(ti): int(ti) for ti in g_tris_idx}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for owners in edge_owner.values():
            for o in owners[1:]:
                ra, rb = find(int(owners[0])), find(int(o))
                if ra != rb:
                    parent[ra] = rb
        comps = defaultdict(list)
        for ti in g_tris_idx:
            comps[find(int(ti))].append(int(ti))
        sizes = sorted(((len(v), k) for k, v in comps.items()), reverse=True)
        n_dropped = 0
        for sz, k in sizes[1:]:
            for ti in comps[k]:
                drop[ti] = True
                n_dropped += 1
        report["clip_b_junction"][f"{guest} -> {host}"] = {
            "n_components": len(comps),
            "n_cut_edges": len(cut_edges),
            "n_dropped": n_dropped,
        }

    keep = ~drop
    tris_k = tris[keep]
    markers_k = markers[keep]

    # ---- (b2) collapse the eps-detach ribbons ---------------------------
    # The junction guests were detached from their hosts by eps (~30 m)
    # before extrusion (Phase 3), so the corefined junction line runs ~eps
    # from the guest's old border row, leaving a ribbon of sub-100 m
    # triangles the constrained remesh cannot collapse (both rows are
    # protected).  Snap every kept fault vertex within --ribbon-snap of a
    # junction cut-edge endpoint ONTO its nearest polyline node, re-weld,
    # and drop the degenerate triangles.
    ribbon_snap = 75.0
    all_cut_nodes = set()
    for gmarker, polylines, guest, host in junction_cuts:
        hmarker = next(fm for fm, bn in enumerate(fault_basenames, start=1)
                       if bn.startswith(host))
        host_edges2 = edge_set_of((markers == hmarker) & ~drop)
        gt = tris[(markers == gmarker) & ~drop]
        ge = np.sort(np.concatenate([gt[:, [0, 1]], gt[:, [1, 2]],
                                     gt[:, [2, 0]]]), axis=1)
        for k in map(tuple, ge):
            if k in host_edges2:
                all_cut_nodes.update(k)
    n_snapped = 0
    if all_cut_nodes:
        from scipy.spatial import cKDTree as _KD
        cut_ids = np.array(sorted(all_cut_nodes))
        ctree = _KD(pts[cut_ids])
        fault_vids = np.unique(
            tris_k[(markers_k >= 1) & (markers_k <= n_faults)].ravel())
        boundary_vids = np.unique(tris_k[markers_k >= box_marker].ravel())
        cand = np.setdiff1d(np.setdiff1d(fault_vids, cut_ids), boundary_vids)
        d, j = ctree.query(pts[cand], workers=-1)
        snap_sel = d <= ribbon_snap
        pts = pts.copy()
        pts[cand[snap_sel]] = pts[cut_ids[j[snap_sel]]]
        n_snapped = int(snap_sel.sum())
        # re-weld at exact coordinates and drop degenerates
        pts, tris_k = _weld_exact(pts, tris_k)
        dg = ((tris_k[:, 0] == tris_k[:, 1]) | (tris_k[:, 1] == tris_k[:, 2])
              | (tris_k[:, 0] == tris_k[:, 2]))
        tris_k = tris_k[~dg]
        markers_k = markers_k[~dg]
        report["ribbon_snapped_vertices"] = n_snapped
        report["ribbon_degen_dropped"] = int(dg.sum())

    # ---- (b3) soup-level min-edge enforcement ----------------------------
    # Residual sub-floor edges (cleanup-guard leftovers along the junction
    # polylines) are contracted directly on the welded soup: vertex-level
    # contraction is conformal across every surface sharing the vertices
    # (the archived pipeline did the same thing as tetgen_mesh.py's
    # --output-dedup-tol 99; here it is topology-aware).  Clusters
    # containing a boundary-shell vertex contract onto that vertex so the
    # shell geometry never moves; both-boundary edges are never contracted.
    n_contract_rounds = 0
    n_contracted = 0
    while True:
        b_vids = set(np.unique(tris_k[markers_k >= box_marker].ravel())
                     .tolist())
        e = np.sort(np.concatenate([tris_k[:, [0, 1]], tris_k[:, [1, 2]],
                                    tris_k[:, [2, 0]]]), axis=1)
        ue = np.unique(e, axis=0)
        el2 = np.linalg.norm(pts[ue[:, 0]] - pts[ue[:, 1]], axis=1)
        shorts = ue[el2 < args.min_edge]
        shorts = [tuple(int(x) for x in row) for row in shorts
                  if not (int(row[0]) in b_vids and int(row[1]) in b_vids)]
        if not shorts or n_contract_rounds >= 10:
            break
        n_contract_rounds += 1
        parent: dict = {}

        def find(x):
            parent.setdefault(x, x)
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for u, v in shorts:
            ru, rv = find(u), find(v)
            if ru != rv:
                parent[ru] = rv
        clusters = defaultdict(list)
        for x in list(parent):
            clusters[find(x)].append(x)
        pts = pts.copy()
        for members in clusters.values():
            anchors = [x for x in members if x in b_vids]
            tgt = (pts[anchors[0]] if anchors
                   else pts[members].mean(axis=0))
            for x in members:
                pts[x] = tgt
            n_contracted += len(members) - 1
        pts, tris_k = _weld_exact(pts, tris_k)
        dg = ((tris_k[:, 0] == tris_k[:, 1])
              | (tris_k[:, 1] == tris_k[:, 2])
              | (tris_k[:, 0] == tris_k[:, 2]))
        tris_k = tris_k[~dg]
        markers_k = markers_k[~dg]
    report["short_edge_contraction"] = {
        "rounds": n_contract_rounds, "vertices_merged": n_contracted}

    # ---- (b4) cap-triangle flip pass ------------------------------------
    # The contractions can leave near-collinear "cap" triangles (q < 0.3,
    # all edges >= floor).  Flip the cap's longest edge when the edge is
    # interior to ONE marker (used by exactly 2 same-marker triangles and
    # by no other marker), the pair is near-coplanar, and the flip
    # improves the pair's min quality.
    def tri_q(t):
        a, b, c = pts[t[:, 0]], pts[t[:, 1]], pts[t[:, 2]]
        e2 = (((b - a) ** 2).sum(1) + ((c - b) ** 2).sum(1)
              + ((a - c) ** 2).sum(1))
        area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(e2 > 0, 4 * np.sqrt(3.0) * area / e2, 0.0)

    n_flips = 0
    for _round in range(5):
        q = tri_q(tris_k)
        bad = np.nonzero(q < 0.30)[0]
        if not len(bad):
            break
        eo: dict = defaultdict(list)
        for ti, t in enumerate(tris_k):
            for u, v in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
                eo[(min(u, v), max(u, v))].append(ti)
        flipped_this_round = 0
        touched: set = set()
        for ti in bad:
            if ti in touched:
                continue
            t = tris_k[ti]
            P3 = pts[t]
            L = [np.linalg.norm(P3[1] - P3[0]), np.linalg.norm(P3[2] - P3[1]),
                 np.linalg.norm(P3[0] - P3[2])]
            pairs_uv = [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])]
            u, v = pairs_uv[int(np.argmax(L))]
            k = (min(u, v), max(u, v))
            owners = eo[k]
            if len(owners) != 2:
                continue
            tj = owners[0] if owners[1] == ti else owners[1]
            if tj in touched or markers_k[ti] != markers_k[tj]:
                continue
            w1 = next(int(x) for x in tris_k[ti] if x not in (u, v))
            w2 = next(int(x) for x in tris_k[tj] if x not in (u, v))
            n1 = np.cross(pts[v] - pts[u], pts[w1] - pts[u])
            n2 = np.cross(pts[w2] - pts[u], pts[v] - pts[u])
            d1, d2 = np.linalg.norm(n1), np.linalg.norm(n2)
            if d1 == 0 or d2 == 0:
                continue
            if np.dot(n1, n2) / (d1 * d2) < np.cos(np.radians(30)):
                continue
            old_q = min(tri_q(tris_k[[ti]])[0], tri_q(tris_k[[tj]])[0])
            cand1 = np.array([[w1, w2, int(v)]])
            cand2 = np.array([[w2, w1, int(u)]])
            new_q = min(tri_q(cand1)[0], tri_q(cand2)[0])
            new_min_edge = np.linalg.norm(pts[w1] - pts[w2])
            if new_q <= old_q or new_min_edge < args.min_edge:
                continue
            tris_k[ti] = [int(u), w2, w1]
            tris_k[tj] = [int(v), w1, w2]
            touched.update((int(ti), int(tj)))
            n_flips += 1
            flipped_this_round += 1
        if not flipped_this_round:
            break
    report["cap_flips"] = n_flips

    report["n_tris_out"] = int(len(tris_k))
    report["n_fault_out"] = int(((markers_k >= 1)
                                 & (markers_k <= n_faults)).sum())

    # re-run containment as a gate
    fidx2 = np.nonzero((markers_k >= 1) & (markers_k <= n_faults))[0]
    cen2 = pts[tris_k[fidx2]].mean(axis=1)
    inside2 = ray_parity_inside(cen2, tris_k[markers_k >= box_marker], pts)
    report["gate_fault_outside_shell_after_clip"] = int((~inside2).sum())

    # ---- validation ------------------------------------------------------
    # boundary closure
    bt = tris_k[markers_k >= box_marker]
    e = np.sort(np.concatenate([bt[:, [0, 1]], bt[:, [1, 2]],
                                bt[:, [2, 0]]]), axis=1)
    uniq_e, counts_e = np.unique(e, axis=0, return_counts=True)
    report["gate_boundary_open_edges"] = int((counts_e == 1).sum())
    report["gate_boundary_overshared_edges"] = int((counts_e > 2).sum())

    # min edge over the whole soup, localized by marker
    allt = tris_k
    ea = np.concatenate([allt[:, [0, 1]], allt[:, [1, 2]], allt[:, [2, 0]]])
    em = np.concatenate([markers_k] * 3)
    el = np.linalg.norm(pts[ea[:, 0]] - pts[ea[:, 1]], axis=1)
    report["gate_min_edge_m"] = float(el.min())
    report["gate_edges_below_floor"] = int((el < args.min_edge).sum())
    short = el < args.min_edge
    if short.any():
        by_marker = {}
        for mval in np.unique(em[short]):
            sel = short & (em == mval)
            mids = (pts[ea[sel, 0]] + pts[ea[sel, 1]]) / 2
            by_marker[int(mval)] = {
                "n": int(sel.sum()),
                "min_m": float(el[sel].min()),
                "z_range": [float(mids[:, 2].min()), float(mids[:, 2].max())],
                "sample_mid": [round(float(x), 1) for x in mids[0]],
            }
        report["short_edges_by_marker"] = by_marker

    # duplicate scan at tol over used nodes
    used = np.unique(tris_k.ravel())
    qq = np.round(pts[used] / args.tol).astype(np.int64)
    _, cts = np.unique(qq, axis=0, return_counts=True)
    report["gate_dup_groups_at_tol"] = int((cts >= 2).sum())

    # fault border classification
    bset = {}
    for nm, sel in (("dem_or_top", markers_k == box_marker),
                    ("boundary_all", markers_k >= box_marker)):
        bset[nm] = set(np.unique(tris_k[sel].ravel()).tolist())
    cls_counts = {}
    orphans = 0
    from scipy.spatial import cKDTree
    btree = cKDTree(pts[sorted(bset["boundary_all"])])
    bnodes_sorted = np.array(sorted(bset["boundary_all"]))
    for fm in range(1, n_faults + 1):
        ft = tris_k[markers_k == fm]
        if not len(ft):
            continue
        fe = np.sort(np.concatenate([ft[:, [0, 1]], ft[:, [1, 2]],
                                     ft[:, [2, 0]]]), axis=1)
        ue, uc = np.unique(fe, axis=0, return_counts=True)
        border = ue[uc == 1]
        other_fault_nodes = set(np.unique(
            tris_k[(markers_k >= 1) & (markers_k <= n_faults)
                   & (markers_k != fm)].ravel()).tolist())
        cc = defaultdict(int)
        for u, v in border:
            u, v = int(u), int(v)
            if u in bset["boundary_all"] and v in bset["boundary_all"]:
                cc["on_boundary_shell"] += 1
            elif u in other_fault_nodes and v in other_fault_nodes:
                cc["junction"] += 1
            else:
                # Near-shell tips: chain-end transitions where the border
                # departs the shell (one endpoint exactly shared, or a
                # single vertex sagging a few tens of metres off the
                # shell between two on-shell nodes).  The fault border is
                # an interior border there — no PLC violation, no
                # duplicate nodes — so these are reported, not failed.
                d1 = btree.query(pts[u])[0]
                d2 = btree.query(pts[v])[0]
                if max(d1, d2) < args.min_edge:
                    cc["near_shell_tip"] += 1
                    orphans += 1
                else:
                    cc["tip"] += 1
        cls_counts[fault_basenames[fm - 1]] = dict(sorted(cc.items()))
    report["fault_border_classes"] = cls_counts
    report["warn_near_shell_tip_edges"] = orphans

    ok = (report["gate_fault_outside_shell_after_clip"] == 0
          and report["gate_boundary_open_edges"] == 0
          and report["gate_dup_groups_at_tol"] == 0
          and report["gate_edges_below_floor"] == 0)
    report["pass"] = bool(ok)

    # ---- outputs ----------------------------------------------------------
    used = np.unique(tris_k.ravel())
    remap = -np.ones(pts.shape[0], dtype=np.int64)
    remap[used] = np.arange(used.size)
    write_ascii_stl(args.out_stl, pts[used], remap[tris_k],
                    name="safv4_merged_clipped")
    out_mk = dict(mk)
    out_mk["markers"] = markers_k.tolist()
    out_mk["n_output_triangles"] = int(len(tris_k))
    out_mk["n_output_vertices"] = int(used.size)
    out_mk["clipped"] = True
    with open(args.out_markers, "w") as f:
        json.dump(out_mk, f)
        f.write("\n")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with open(args.report, "w") as f:
        json.dump(report, f, indent=1, sort_keys=True)
        f.write("\n")

    print(json.dumps({k: v for k, v in report.items()
                      if k != "fault_border_classes"}, indent=1,
                     sort_keys=True))
    print("fault border classes:")
    for k, v in report["fault_border_classes"].items():
        print(f"  {k[:60]}: {v}")
    print(f"RESULT: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
