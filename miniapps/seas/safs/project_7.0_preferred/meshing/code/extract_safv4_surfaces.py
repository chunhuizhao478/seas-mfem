"""extract_safv4_surfaces.py — Phase 2 of PLAN_mesh_quality_safv4_remesh.

Extract the SAFv4 model surfaces from the frozen SKUA-GOCAD .inp, weld all
coordinate-coincident duplicate nodes (the GOCAD export is split along the
faults and the trace), and emit canonical per-surface STLs + a manifest:

  fault_<CFM basename>.stl   x3   (canonical = minus side)
  dem.stl, bottom.stl, ribbon_<NE30|NE120|NW60|NW150>.stl

The manifest records, per fault, a border classification: every border
edge is assigned to exactly one of {dem (trace), bottom, ribbon_<name>,
fault_<basename> (junction), interior_tip} by post-weld shared-node test,
emitted as ordered chains per class.  Junction chains carry host/guest
(the fault whose border lies ON the other's interior is the guest).

CLI
---
    python extract_safv4_surfaces.py INP --out-dir DIR [--tol 1e-3]
        [--fault-side minus]

Aborts (exit 2) if: the minus-subset-of-plus assumption breaks (>1%
unmatched), or the boundary shell (dem+ribbons+bottom) has an open rim
edge after the weld (GOCAD export inconsistency).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import meshio
import numpy as np

import check_mesh_quality as cmq

# Classification precedence when a border edge qualifies for more than one
# class (corner nodes at chain meets): trace wins, then junction, bottom,
# ribbon; unmatched edges are interior tips.
CLASS_PRECEDENCE = ("dem", "junction", "bottom", "ribbon")


def write_ascii_stl(path, points: np.ndarray, tris: np.ndarray,
                    name: str = "surface") -> None:
    """ASCII STL at %.15g.  Binary STL is float32 by format definition and
    destroys UTM-scale coordinate precision (~0.06 m at y ~ 4e6); the
    corefine pipeline needs >= 1e-6 m conformality, so we write ASCII at
    full double precision (same convention as the archived
    corefine_faults._off_to_stl)."""
    pts = np.asarray(points, dtype=np.float64)
    with open(path, "w") as f:
        f.write(f"solid {name}\n")
        for t in tris:
            a, b, c = pts[t[0]], pts[t[1]], pts[t[2]]
            n = np.cross(b - a, c - a)
            mag = np.linalg.norm(n)
            if mag > 0:
                n = n / mag
            f.write(f"facet normal {n[0]:.15g} {n[1]:.15g} {n[2]:.15g}\n")
            f.write("  outer loop\n")
            for p in (a, b, c):
                f.write(f"    vertex {p[0]:.15g} {p[1]:.15g} {p[2]:.15g}\n")
            f.write("  endloop\nendfacet\n")
        f.write(f"endsolid {name}\n")


def _qkeys(pts: np.ndarray, tol: float) -> np.ndarray:
    return np.round(np.asarray(pts, dtype=np.float64) / tol).astype(np.int64)


# ---------------------------------------------------------------------------
# Weld
# ---------------------------------------------------------------------------

def build_weld_map(pts: np.ndarray, tol: float):
    """Global quantized-key weld: returns (weld (N,) int64 mapping every
    node ID to its group representative = lowest ID in the group,
    n_groups_total)."""
    q = _qkeys(pts, tol)
    _, inverse, counts = np.unique(q, axis=0, return_inverse=True,
                                   return_counts=True)
    weld = np.arange(pts.shape[0], dtype=np.int64)
    dup_buckets = np.nonzero(counts >= 2)[0]
    if dup_buckets.size:
        order = np.argsort(inverse, kind="stable")
        inv_sorted = inverse[order]
        starts = np.searchsorted(inv_sorted, dup_buckets, side="left")
        ends = np.searchsorted(inv_sorted, dup_buckets, side="right")
        for s, e in zip(starts, ends):
            ids = order[s:e]
            weld[ids] = ids.min()
    return weld, int(dup_buckets.size)


# ---------------------------------------------------------------------------
# Border classification
# ---------------------------------------------------------------------------

def border_edges(tris: np.ndarray) -> np.ndarray:
    """(K,2) sorted node pairs used by exactly one triangle."""
    e = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
    e = np.sort(e, axis=1)
    uniq, counts = np.unique(e, axis=0, return_counts=True)
    return uniq[counts == 1]


def chains_from_edges(edges: np.ndarray) -> list[list[int]]:
    """Order a set of edges into vertex chains (open chains first from
    degree-1 endpoints, then any cycles), deterministically."""
    from collections import defaultdict
    adj = defaultdict(list)
    eset = set()
    for u, v in edges:
        u, v = int(u), int(v)
        adj[u].append(v)
        adj[v].append(u)
        eset.add((min(u, v), max(u, v)))
    for k in adj:
        adj[k].sort()
    used = set()
    chains = []

    def walk(start):
        chain = [start]
        cur = start
        while True:
            nxt = None
            for nb in adj[cur]:
                key = (min(cur, nb), max(cur, nb))
                if key in eset and key not in used:
                    nxt = nb
                    used.add(key)
                    break
            if nxt is None:
                return chain
            chain.append(nxt)
            cur = nxt

    endpoints = sorted(k for k in adj if len(adj[k]) == 1)
    for ep in endpoints:
        if all((min(ep, nb), max(ep, nb)) in used for nb in adj[ep]):
            continue
        chains.append(walk(ep))
    # Remaining cycles.
    for u, v in sorted(eset):
        if (u, v) not in used:
            used.add((u, v))
            chain = [u, v]
            cur = v
            while True:
                nxt = None
                for nb in adj[cur]:
                    key = (min(cur, nb), max(cur, nb))
                    if key in eset and key not in used:
                        nxt = nb
                        used.add(key)
                        break
                if nxt is None:
                    break
                chain.append(nxt)
                cur = nxt
            chains.append(chain)
    return chains


# ---------------------------------------------------------------------------
# Main extraction
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inp", type=Path)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--tol", type=float, default=1e-3)
    ap.add_argument("--fault-side", choices=("minus", "plus"), default="minus")
    args = ap.parse_args(argv)

    if not args.inp.is_file():
        print(f"error: {args.inp} not found", file=sys.stderr)
        return 1
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"reading {args.inp.name} ...", flush=True)
    m = meshio.read(str(args.inp))
    pts = np.asarray(m.points, dtype=np.float64)
    tet_blocks = [cb.data for cb in m.cells if cb.type == "tetra"]
    all_tetra = np.concatenate(tet_blocks, axis=0) if len(tet_blocks) > 1 \
        else tet_blocks[0]
    offsets = cmq._build_block_offsets(m)

    raw: dict[str, np.ndarray] = {}
    for base in cmq._surface_names(m):
        tris = cmq._extract_surface_triangles(m, offsets, all_tetra, base)
        if tris is None or len(tris) == 0:
            continue
        name = base[:-5] if base.endswith("_C3D4") else base
        raw[name] = tris.astype(np.int64)
    print(f"  surfaces: {sorted(raw)}", flush=True)

    # ---- identify roles ------------------------------------------------
    minus_names = sorted(n for n in raw if n.endswith("_minus"))
    fault_bases = [n[:-6] for n in minus_names]          # strip _minus
    dem_name = next(n for n in raw if "dem" in n.lower())
    bottom_name = next(n for n in raw if n.lower() == "bottom")
    ribbon_names = sorted(n for n in raw if n.endswith("_ribbon"))
    side = args.fault_side
    other = "plus" if side == "minus" else "minus"

    # ---- minus/plus centroid match (pre-weld, baseline algorithm) ------
    match_report = {}
    for fb in fault_bases:
        ta, tb = raw[f"{fb}_{side}"], raw[f"{fb}_{other}"]
        ca = set(map(tuple, _qkeys(pts[ta].mean(axis=1), args.tol)))
        cb = set(map(tuple, _qkeys(pts[tb].mean(axis=1), args.tol)))
        unmatched = len(ca - cb)
        match_report[fb] = {
            "n_tri_canonical": int(ta.shape[0]),
            "n_tri_other_side": int(tb.shape[0]),
            "matched": len(ca & cb),
            f"unmatched_{side}": unmatched,
            f"extras_{other}": len(cb - ca),
        }
        frac = unmatched / max(1, ta.shape[0])
        print(f"  {fb}: matched {len(ca & cb)}/{ta.shape[0]} "
              f"({other} extras {len(cb - ca)})", flush=True)
        if frac > 0.01:
            print(f"ABORT: {fb}: {unmatched} canonical-{side} triangles "
                  f"({100*frac:.2f}% > 1%) have no {other}-side match — the "
                  f"subset assumption broke; user decides.", file=sys.stderr)
            return 2

    # ---- global weld ----------------------------------------------------
    print("welding ...", flush=True)
    weld, n_groups_total = build_weld_map(pts, args.tol)

    # Fault-subset group count, computed exactly as the diagnosis baseline:
    # scan restricted to fault-surface node IDs (both sides).  (NOT the
    # "groups containing >=1 fault node" of the global scan — a fault node
    # coincident only with a non-fault node is not a baseline group.)
    fault_ids_raw = np.unique(np.concatenate(
        [raw[f"{fb}_{s}"].ravel() for fb in fault_bases
         for s in ("minus", "plus")]))
    qf = _qkeys(pts[fault_ids_raw], args.tol)
    _, counts_f = np.unique(qf, axis=0, return_counts=True)
    n_groups_fault_subset = int((counts_f >= 2).sum())
    print(f"  weld groups: total={n_groups_total:,} "
          f"fault_subset={n_groups_fault_subset:,}", flush=True)

    # ---- welded canonical surfaces --------------------------------------
    surf: dict[str, np.ndarray] = {}
    out_name: dict[str, str] = {}
    n_degenerate = 0
    members = [(f"fault_{fb}", raw[f"{fb}_{side}"]) for fb in fault_bases]
    members += [("dem", raw[dem_name]), ("bottom", raw[bottom_name])]
    members += [(f"ribbon_{rn.replace('_ribbon', '')}", raw[rn])
                for rn in ribbon_names]
    for oname, tris in members:
        tw = weld[tris]
        degen = ((tw[:, 0] == tw[:, 1]) | (tw[:, 1] == tw[:, 2])
                 | (tw[:, 0] == tw[:, 2]))
        n_degenerate += int(degen.sum())
        surf[oname] = tw[~degen]
        out_name[oname] = f"{oname}.stl"
    if n_degenerate:
        print(f"  WARNING: {n_degenerate} triangles became degenerate "
              f"after weld (dropped)", flush=True)

    # ---- DEM re-triangulation (crack repair) -----------------------------
    # Measured (2026-06-12): the GOCAD export splits the DEM along the fault
    # traces with DIFFERENT trace segmentations per side (e.g. SAF plus 622
    # trace nodes vs minus 380), so the welded DEM has T-junction cracks
    # (977 open edges, all with both endpoints on welded fault nodes).  A
    # DEM is single-valued in (x,y): rebuild it as the 2-D Delaunay
    # triangulation of its own (unchanged) vertex set — the crack vanishes,
    # and rim conformity with the ribbons is preserved because the rim
    # vertices are reused verbatim on the convex footprint boundary.
    dem_tris_old = surf["dem"]
    dem_used = np.unique(dem_tris_old.ravel())
    fault_node_union = set()
    for fb in fault_bases:
        fault_node_union.update(
            np.unique(weld[raw[f"{fb}_minus"]].ravel()).tolist())
        fault_node_union.update(
            np.unique(weld[raw[f"{fb}_plus"]].ravel()).tolist())

    # xy uniqueness (heightfield precondition).
    qxy = _qkeys(pts[dem_used][:, :2], args.tol)
    if np.unique(qxy, axis=0).shape[0] != dem_used.size:
        print("ABORT: DEM vertex set is not single-valued in (x,y) — "
              "cannot re-triangulate as a heightfield.", file=sys.stderr)
        return 2

    # Rim = original DEM border edges that are NOT trace-crack edges.
    dem_border = border_edges(dem_tris_old)
    crack_mask = np.array([(int(u) in fault_node_union and
                            int(v) in fault_node_union)
                           for u, v in dem_border])
    rim_edges = dem_border[~crack_mask]
    n_crack_edges = int(crack_mask.sum())

    from scipy.spatial import Delaunay
    tri2d = Delaunay(pts[dem_used][:, :2])
    new_tris = dem_used[tri2d.simplices]
    # orient consistently: CCW in xy (normal +z)
    v01 = pts[new_tris[:, 1], :2] - pts[new_tris[:, 0], :2]
    v02 = pts[new_tris[:, 2], :2] - pts[new_tris[:, 0], :2]
    cwm = (v01[:, 0] * v02[:, 1] - v01[:, 1] * v02[:, 0]) < 0
    new_tris[cwm] = new_tris[cwm][:, [0, 2, 1]]

    # Drop degenerate (zero-xy-area) Delaunay slivers from collinear rim pts.
    area2 = np.abs(v01[:, 0] * v02[:, 1] - v01[:, 1] * v02[:, 0])
    new_tris = new_tris[area2 > 1e-9]

    # The footprint is CONCAVE (measured: convex-hull chords up to 430 km
    # appear as open edges otherwise): Delaunay triangulates the convex
    # hull, so drop every triangle whose centroid falls outside the rim
    # polygon (the ordered non-crack DEM border cycle).
    rim_chains = chains_from_edges(rim_edges)
    if len(rim_chains) != 1:
        print(f"ABORT: DEM rim is not a single closed chain "
              f"(got {len(rim_chains)} chains of sizes "
              f"{[len(c) for c in rim_chains]})", file=sys.stderr)
        return 2
    poly = pts[rim_chains[0]][:, :2]

    def _in_poly(xy, poly):
        inside = np.zeros(len(xy), dtype=bool)
        x, y = xy[:, 0], xy[:, 1]
        x0, y0 = poly[:, 0], poly[:, 1]
        x1, y1 = np.roll(x0, -1), np.roll(y0, -1)
        for i in range(len(poly)):
            dy = y1[i] - y0[i]
            if dy == 0.0:
                continue
            cond = (y0[i] > y) != (y1[i] > y)
            xin = (x1[i] - x0[i]) * (y - y0[i]) / dy + x0[i]
            inside ^= cond & (x < xin)
        return inside

    cen = pts[new_tris].mean(axis=1)[:, :2]
    inside = _in_poly(cen, poly)
    n_outside = int((~inside).sum())
    new_tris = new_tris[inside]

    # Verify every rim edge survives in the new triangulation.
    ne = np.sort(np.concatenate([new_tris[:, [0, 1]], new_tris[:, [1, 2]],
                                 new_tris[:, [2, 0]]]), axis=1)
    ne_set = set(map(tuple, ne))
    missing_rim = [tuple(int(x) for x in e2) for e2 in np.sort(rim_edges, axis=1)
                   if tuple(e2) not in ne_set]
    if missing_rim:
        print(f"ABORT: DEM re-triangulation lost {len(missing_rim)} rim "
              f"edge(s) (rim not reproducible by Delaunay); first: "
              f"{missing_rim[:3]}", file=sys.stderr)
        return 2
    print(f"  DEM re-triangulated: {dem_tris_old.shape[0]:,} -> "
          f"{new_tris.shape[0]:,} tris ({n_crack_edges} crack edges "
          f"resolved; {n_outside} outside-footprint tris dropped; rim "
          f"{rim_edges.shape[0]} edges preserved)", flush=True)
    surf["dem"] = new_tris.astype(np.int64)
    dem_retri_info = {
        "applied": True,
        "n_tri_before": int(dem_tris_old.shape[0]),
        "n_tri_after": int(new_tris.shape[0]),
        "n_crack_edges_resolved": n_crack_edges,
        "n_outside_footprint_dropped": n_outside,
        "n_rim_edges_preserved": int(rim_edges.shape[0]),
    }

    node_set = {n: set(np.unique(t.ravel()).tolist())
                for n, t in surf.items()}

    # ---- boundary shell watertightness ----------------------------------
    boundary_names = ["dem", "bottom"] + sorted(
        n for n in surf if n.startswith("ribbon_"))
    btris = np.concatenate([surf[n] for n in boundary_names], axis=0)
    e = np.sort(np.concatenate([btris[:, [0, 1]], btris[:, [1, 2]],
                                btris[:, [2, 0]]]), axis=1)
    uniq_e, counts_e = np.unique(e, axis=0, return_counts=True)
    n_open = int((counts_e == 1).sum())
    n_nonmanifold = int((counts_e > 2).sum())
    print(f"  boundary shell: {uniq_e.shape[0]:,} edges, open={n_open}, "
          f"count>2={n_nonmanifold}", flush=True)
    if n_open > 0:
        print(f"ABORT: boundary shell has {n_open} open rim edge(s) after "
              f"weld at tol={args.tol} — GOCAD export inconsistency; user "
              f"decides.", file=sys.stderr)
        return 2

    # ---- per-fault border classification --------------------------------
    fault_names = [f"fault_{fb}" for fb in fault_bases]
    borders = {}
    for fn in fault_names:
        bedges = border_edges(surf[fn])
        classes = {}
        for u, v in bedges:
            u, v = int(u), int(v)
            assigned = None
            # precedence: dem > junction > bottom > ribbon > tip
            if u in node_set["dem"] and v in node_set["dem"]:
                assigned = "dem"
            if assigned is None:
                for other_fn in fault_names:
                    if other_fn == fn:
                        continue
                    if u in node_set[other_fn] and v in node_set[other_fn]:
                        assigned = other_fn
                        break
            if assigned is None and u in node_set["bottom"] \
                    and v in node_set["bottom"]:
                assigned = "bottom"
            if assigned is None:
                for rn in boundary_names[2:]:
                    if u in node_set[rn] and v in node_set[rn]:
                        assigned = rn
                        break
            if assigned is None:
                assigned = "interior_tip"
            classes.setdefault(assigned, []).append((u, v))
        borders[fn] = classes

    # Junction host/guest: the guest's junction-contact nodes lie in the
    # host's INTERIOR (not on the host's own border).
    fault_border_nodes = {fn: set(np.unique(border_edges(surf[fn]).ravel()).tolist())
                          for fn in fault_names}
    junction_info = {}
    for fn in fault_names:
        for cls, edges in borders[fn].items():
            if not cls.startswith("fault_"):
                continue
            nodes = {n for e2 in edges for n in e2}
            on_host_border = sum(1 for n in nodes
                                 if n in fault_border_nodes[cls])
            junction_info[f"{fn} -> {cls}"] = {
                "guest": fn, "host": cls,
                "n_contact_nodes": len(nodes),
                "n_on_host_border": on_host_border,
                "host_interior_contact": bool(
                    on_host_border < len(nodes)),
            }

    # ---- manifest --------------------------------------------------------
    def edge_stats(tris):
        tri_e, _a, _q = cmq._tri_metrics(tris, pts)
        return {"min": float(tri_e.min()), "median": float(np.median(tri_e)),
                "max": float(tri_e.max())}

    manifest = {
        "source_inp": args.inp.name,
        "tol_m": args.tol,
        "fault_side": side,
        "weld": {
            "n_groups_total": n_groups_total,
            "n_groups_fault_subset": n_groups_fault_subset,
            "n_degenerate_dropped": n_degenerate,
        },
        "minus_plus_match": match_report,
        "dem_retriangulation": dem_retri_info,
        "boundary_shell": {
            "names": boundary_names,
            "n_edges": int(uniq_e.shape[0]),
            "n_open": n_open,
            "n_count_gt2": n_nonmanifold,
        },
        "junctions": junction_info,
        "meshes": [],
    }
    for oname in ([f"fault_{fb}" for fb in fault_bases]
                  + ["dem", "bottom"] + boundary_names[2:]):
        tris = surf[oname]
        used = np.unique(tris.ravel())
        entry = {
            "basename": oname,
            "file": out_name[oname],
            "n_tri": int(tris.shape[0]),
            "n_vertices": int(used.size),
            "bbox": {
                "x_lo": float(pts[used, 0].min()), "x_hi": float(pts[used, 0].max()),
                "y_lo": float(pts[used, 1].min()), "y_hi": float(pts[used, 1].max()),
                "z_lo": float(pts[used, 2].min()), "z_hi": float(pts[used, 2].max()),
            },
            "edge_stats": edge_stats(tris),
        }
        if oname.startswith("fault_"):
            cls_out = {}
            for cls, edges in sorted(borders[oname].items()):
                chains = chains_from_edges(np.asarray(edges))
                cls_out[cls] = {
                    "n_edges": len(edges),
                    "n_nodes": len({n for e2 in edges for n in e2}),
                    "chains": [
                        {
                            "n": len(ch),
                            "on_boundary_start": any(
                                ch[0] in node_set[rn]
                                for rn in boundary_names[2:]),
                            "on_boundary_end": any(
                                ch[-1] in node_set[rn]
                                for rn in boundary_names[2:]),
                            "coords": [[float(c) for c in pts[n]]
                                       for n in ch],
                        } for ch in chains],
                }
            entry["border_classes"] = cls_out
        manifest["meshes"].append(entry)

    # ---- STL output -------------------------------------------------------
    for oname, tris in surf.items():
        used = np.unique(tris.ravel())
        remap = -np.ones(pts.shape[0], dtype=np.int64)
        remap[used] = np.arange(used.size)
        out = args.out_dir / out_name[oname]
        write_ascii_stl(out, pts[used], remap[tris], name=oname)
        print(f"  wrote {out.name}: {tris.shape[0]:,} tris, "
              f"{used.size:,} verts", flush=True)

    man_path = args.out_dir / "manifest.json"
    with open(man_path, "w") as f:
        json.dump(manifest, f, indent=1, sort_keys=True)
        f.write("\n")
    print(f"wrote {man_path}", flush=True)

    # ---- summary of border classes ---------------------------------------
    print("\nborder classification (edges per class):")
    for fn in fault_names:
        parts = ", ".join(f"{cls}:{len(edges)}"
                          for cls, edges in sorted(borders[fn].items()))
        print(f"  {fn}: {parts}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
