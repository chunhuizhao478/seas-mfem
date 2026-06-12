"""check_fault_conformity.py — conformity + sliver gate harness (Phase 1 of
PLAN_mesh_quality_safv4_remesh_2026-06-12.md).

Proves/refutes the three SAFv4 mesh goals on any mesh (.inp, .msh, .vtu):

  goal 1 (fault sizing)    — via the sliver/sub-100m localization report;
  goal 2 (slivers)         — Joe-Liu eta census + localization (z-histogram,
                             fault/DEM attribution, distance-to-fault);
  goal 3 (no duplicates)   — coordinate-coincident distinct-ID node groups,
                             plus/minus side conformity, fault embedding
                             (every fault triangle a face of exactly 2 tets).

Exit code 0 iff: zero duplicate groups (mesh-wide) AND, when tets are
present, every fault triangle is interior (n_interior == n_tri per fault).

CLI
---
    python check_fault_conformity.py MESH [--fault-tags 101,102,103]
        [--dem-name NAME | --dem-tag N] [--pairs auto|none] [--tol 1e-3]
        [--json out.json]
    python check_fault_conformity.py --self-test

`--pairs auto` (default for .inp) pairs every `*_minus*` surface with its
`*_plus*` sibling for the conformity table.

The duplicate-pair / conformity quantities reproduce the 2026-06-12
diagnosis baseline algorithm exactly (quantized-key dicts; see
baseline_safv4_2km_topo_conformity.json).
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

CHUNK = 2_000_000
QBINS_Z = [-20000.0, -15000.0, -10000.0, -5000.0, -2500.0, -1000.0, 0.0, 4000.0]


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

@dataclass
class SurfaceModel:
    points: np.ndarray                       # (N,3) float64
    tets: np.ndarray | None                  # (M,4) int64 or None
    surfaces: dict = field(default_factory=dict)   # name -> (K,3) int64
    fault_names: list = field(default_factory=list)
    dem_name: str | None = None


def _strip_c3d4(name: str) -> str:
    return name[:-5] if name.endswith("_C3D4") else name


def load_surfaces(path: Path, fault_tags: list[int] | None = None,
                  dem_name: str | None = None,
                  dem_tag: int | None = None) -> SurfaceModel:
    """Load mesh + named surface triangle sets from .inp / .msh / .vtu.

    .inp: GOCAD Abaqus surfaces via check_mesh_quality's FACE_NODES helpers.
    .msh/.vtu: triangles grouped by gmsh:physical tag (name = 'tag_<N>').
    """
    import meshio
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"mesh not found: {path}")
    suffix = path.suffix.lower()
    m = meshio.read(str(path))
    pts = np.asarray(m.points, dtype=np.float64)

    tet_blocks = [cb.data for cb in m.cells if cb.type == "tetra"]
    tets = (np.concatenate(tet_blocks, axis=0).astype(np.int64)
            if tet_blocks else None)

    surfaces: dict = {}
    if suffix == ".inp":
        import check_mesh_quality as cmq
        offsets = cmq._build_block_offsets(m)
        all_tetra = tets
        if all_tetra is None:
            raise ValueError(f"{path}: .inp has no tetra block")
        for base in cmq._surface_names(m):
            tris = cmq._extract_surface_triangles(m, offsets, all_tetra, base)
            if tris is None or len(tris) == 0:
                continue
            surfaces[_strip_c3d4(base)] = tris.astype(np.int64)
    else:
        # .msh / .vtu: group triangles by gmsh:physical.
        tri_data = []     # list of (tris_block, tags_block or None)
        for bi, cb in enumerate(m.cells):
            if cb.type != "triangle":
                continue
            tags = None
            if "gmsh:physical" in m.cell_data:
                tags = np.asarray(m.cell_data["gmsh:physical"][bi])
            tri_data.append((np.asarray(cb.data, dtype=np.int64), tags))
        by_tag: dict = {}
        for tris, tags in tri_data:
            if tags is None:
                by_tag.setdefault("untagged", []).append(tris)
                continue
            for tag in np.unique(tags):
                sel = tris[tags == tag]
                by_tag.setdefault(int(tag), []).append(sel)
        for tag, blocks in by_tag.items():
            name = f"tag_{tag}" if isinstance(tag, int) else str(tag)
            surfaces[name] = np.concatenate(blocks, axis=0)

    model = SurfaceModel(points=pts, tets=tets, surfaces=surfaces)

    # Fault / DEM identification.
    if suffix == ".inp":
        model.fault_names = sorted(n for n in surfaces if "fault" in n.lower())
        dems = [n for n in surfaces if "dem" in n.lower()]
        model.dem_name = dem_name or (dems[0] if dems else None)
    else:
        wanted = fault_tags if fault_tags else []
        model.fault_names = [f"tag_{t}" for t in wanted if f"tag_{t}" in surfaces]
        if dem_tag is not None and f"tag_{dem_tag}" in surfaces:
            model.dem_name = f"tag_{dem_tag}"
        elif dem_name and dem_name in surfaces:
            model.dem_name = dem_name
    return model


# ---------------------------------------------------------------------------
# Core checks
# ---------------------------------------------------------------------------

def _qkeys(pts: np.ndarray, tol: float) -> np.ndarray:
    """(N,3) int64 quantized keys (same scheme as the diagnosis baseline)."""
    return np.round(np.asarray(pts, dtype=np.float64) / tol).astype(np.int64)


def duplicate_node_groups(pts: np.ndarray, ids: np.ndarray | None = None,
                          tol: float = 1e-3) -> list[list[int]]:
    """Groups of >= 2 distinct node IDs sharing one quantized coordinate key.

    ids=None scans every point (IDs = row indices); otherwise only the given
    node IDs are scanned.
    """
    if ids is None:
        ids = np.arange(pts.shape[0], dtype=np.int64)
    else:
        ids = np.asarray(ids, dtype=np.int64)
    if ids.size == 0:
        return []
    q = _qkeys(pts[ids], tol)
    uniq, inverse, counts = np.unique(q, axis=0, return_inverse=True,
                                      return_counts=True)
    dup_buckets = np.nonzero(counts >= 2)[0]
    if dup_buckets.size == 0:
        return []
    order = np.argsort(inverse, kind="stable")
    inv_sorted = inverse[order]
    ids_sorted = ids[order]
    starts = np.searchsorted(inv_sorted, dup_buckets, side="left")
    ends = np.searchsorted(inv_sorted, dup_buckets, side="right")
    return [ids_sorted[s:e].tolist() for s, e in zip(starts, ends)]


def pair_conformity(model: SurfaceModel, name_a: str, name_b: str,
                    tol: float = 1e-3) -> dict:
    """Conformity quantities for a (minus, plus) surface pair.

    Reproduces the diagnosis-baseline algorithm: per-surface dicts mapping
    quantized coordinate key -> node ID (last writer wins inside one
    surface), coincident = key intersection, duplicated = coincident keys
    whose representative IDs differ; centroid matching at the same tol.
    """
    pts = model.points
    tris_a = model.surfaces[name_a]
    tris_b = model.surfaces[name_b]
    ids_a = np.unique(tris_a.ravel())
    ids_b = np.unique(tris_b.ravel())

    def keymap(ids):
        q = _qkeys(pts[ids], tol)
        return {tuple(row): int(i) for row, i in zip(q, ids)}

    km, kp = keymap(ids_a), keymap(ids_b)
    coincident = set(km).intersection(kp)
    duplicated = sum(1 for k in coincident if km[k] != kp[k])
    shared_ids = int(np.intersect1d(ids_a, ids_b).size)

    def centroid_keyset(tris):
        c = pts[tris].mean(axis=1)
        return set(map(tuple, _qkeys(c, tol)))

    ca, cb = centroid_keyset(tris_a), centroid_keyset(tris_b)
    matches = len(ca & cb)
    return {
        "n_tri_a": int(tris_a.shape[0]),
        "n_tri_b": int(tris_b.shape[0]),
        "n_nodes_a": int(ids_a.size),
        "n_nodes_b": int(ids_b.size),
        "shared_node_ids": shared_ids,
        "coincident_coord_pairs": len(coincident),
        "duplicated_pairs": int(duplicated),
        "centroid_matches": int(matches),
        "extra_a": int(len(ca - cb)),
        "extra_b": int(len(cb - ca)),
    }


def fault_embedding(model: SurfaceModel, fault_name: str) -> dict:
    """Per-fault-triangle tet-face incidence, streaming over tets.

    Builds the key set of FAULT triangles only (sorted node triples), then
    streams tets in CHUNK blocks, pre-filtering candidate faces on the
    minimum vertex ID before set lookup.  Memory O(n_fault_tri + chunk);
    a global face->tet map is forbidden (would exceed 12 GB at 31 M tets).
    """
    tris = model.surfaces[fault_name]
    if model.tets is None:
        raise ValueError("fault_embedding requires a volume mesh")
    fkeys = np.sort(tris, axis=1)
    fault_set = {tuple(row) for row in fkeys}
    min_ids = np.unique(fkeys[:, 0])
    incidence: dict = {k: 0 for k in fault_set}

    FACE_COLS = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
    tets = model.tets
    for s in range(0, tets.shape[0], CHUNK):
        T = tets[s:s + CHUNK]
        for cols in FACE_COLS:
            faces = np.sort(T[:, list(cols)], axis=1)
            mask = np.isin(faces[:, 0], min_ids)
            if not mask.any():
                continue
            for row in faces[mask]:
                k = (int(row[0]), int(row[1]), int(row[2]))
                if k in incidence:
                    incidence[k] += 1
    counts = np.fromiter(incidence.values(), dtype=np.int64)
    return {
        "n_tri": int(tris.shape[0]),
        "n_unique_tri": int(len(fault_set)),
        "n_interior": int((counts == 2).sum()),
        "n_boundary": int((counts == 1).sum()),
        "n_orphan": int((counts == 0).sum()),
        "n_over": int((counts > 2).sum()),
    }


def _tet_eta_chunked(pts: np.ndarray, tets: np.ndarray):
    """Yield (start, eta, min_edge) per CHUNK of tets (Joe-Liu eta)."""
    epairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    for s in range(0, tets.shape[0], CHUNK):
        T = tets[s:s + CHUNK]
        A = pts[T[:, [p[0] for p in epairs]]]
        B = pts[T[:, [p[1] for p in epairs]]]
        e2 = ((A - B) ** 2).sum(axis=2)
        p0, p1, p2, p3 = pts[T[:, 0]], pts[T[:, 1]], pts[T[:, 2]], pts[T[:, 3]]
        V = np.abs(np.einsum("ij,ij->i", p1 - p0,
                             np.cross(p2 - p0, p3 - p0))) / 6.0
        sumL2 = e2.sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            eta = np.where(sumL2 > 0.0,
                           12.0 * np.cbrt((3.0 * V) ** 2) / sumL2, 0.0)
        yield s, eta, np.sqrt(e2.min(axis=1))


def sliver_report(model: SurfaceModel, eta_floor: float = 0.1,
                  min_edge_floor: float = 100.0) -> dict:
    """Sliver census + localization (z-histogram, fault/DEM attribution,
    distance-to-fault percentiles) and the sub-min-edge tet census."""
    if model.tets is None:
        raise ValueError("sliver_report requires a volume mesh")
    pts, tets = model.points, model.tets

    sliver_ids, short_ids = [], []
    for s, eta, emin in _tet_eta_chunked(pts, tets):
        sliver_ids.append(np.nonzero(eta <= eta_floor)[0] + s)
        short_ids.append(np.nonzero(emin < min_edge_floor)[0] + s)
    sliver_ids = np.concatenate(sliver_ids) if sliver_ids else np.array([], dtype=np.int64)
    short_ids = np.concatenate(short_ids) if short_ids else np.array([], dtype=np.int64)

    fault_ids = (np.unique(np.concatenate(
        [model.surfaces[n].ravel() for n in model.fault_names]))
        if model.fault_names else np.array([], dtype=np.int64))
    dem_ids = (np.unique(model.surfaces[model.dem_name].ravel())
               if model.dem_name else np.array([], dtype=np.int64))
    fault_set, dem_set = set(fault_ids.tolist()), set(dem_ids.tolist())

    def touch_counts(tet_rows):
        nf = np.array([sum(1 for v in row if v in fault_set) for row in tet_rows])
        nd = np.array([sum(1 for v in row if v in dem_set) for row in tet_rows])
        return nf, nd

    out = {
        "eta_floor": eta_floor,
        "min_edge_floor_m": min_edge_floor,
        "n_slivers": int(sliver_ids.size),
        "n_short_edge_tets": int(short_ids.size),
    }
    if sliver_ids.size:
        ST = tets[sliver_ids]
        nf, nd = touch_counts(ST)
        cz = pts[ST].mean(axis=1)[:, 2]
        zlo = min(QBINS_Z[0], float(cz.min()))
        zhi = max(QBINS_Z[-1], float(cz.max()) + 1.0)
        bins = [zlo] + [b for b in QBINS_Z[1:-1]] + [zhi]
        hist, edges = np.histogram(cz, bins=bins)
        out["sliver_touch_fault_ge1"] = int((nf >= 1).sum())
        out["sliver_touch_fault_ge3"] = int((nf >= 3).sum())
        out["sliver_touch_dem_ge1"] = int((nd >= 1).sum())
        out["sliver_touch_dem_ge3"] = int((nd >= 3).sum())
        out["sliver_z_hist"] = {
            "edges_m": [float(e) for e in edges],
            "counts": hist.tolist(),
        }
        if fault_ids.size:
            from scipy.spatial import cKDTree
            tree = cKDTree(pts[fault_ids])
            d, _ = tree.query(pts[ST].mean(axis=1), workers=-1)
            out["sliver_dist_to_fault_m"] = {
                f"p{q:02d}": float(np.percentile(d, q))
                for q in (5, 25, 50, 75, 95)}
            out["sliver_frac_within_1km_of_fault"] = float((d < 1000.0).mean())
    if short_ids.size:
        BT = tets[short_ids]
        nf, nd = touch_counts(BT)
        out["short_touch_fault_ge1"] = int((nf >= 1).sum())
        out["short_touch_dem_ge1"] = int((nd >= 1).sum())
    return out


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

def _self_test() -> int:
    """(a) two tets sharing a fault triangle with shared nodes -> PASS;
    (b) the shared triangle's 3 nodes duplicated -> 3 dup groups + an
    orphaned/boundary fault face -> FAIL semantics."""
    # (a) shared nodes
    pts_a = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0],
                      [0, 0, 1], [0, 0, -1]], dtype=np.float64)
    tets_a = np.array([[0, 1, 2, 3], [0, 2, 1, 4]], dtype=np.int64)
    model_a = SurfaceModel(points=pts_a, tets=tets_a,
                           surfaces={"fault": np.array([[0, 1, 2]])},
                           fault_names=["fault"])
    groups_a = duplicate_node_groups(pts_a)
    emb_a = fault_embedding(model_a, "fault")
    ok_a = (len(groups_a) == 0 and emb_a["n_interior"] == 1
            and emb_a["n_boundary"] == 0 and emb_a["n_orphan"] == 0)
    print(f"self-test (a) shared nodes: groups={len(groups_a)} "
          f"emb={emb_a} -> {'PASS' if ok_a else 'FAIL'}")

    # (b) duplicated fault nodes (split mesh)
    pts_b = np.vstack([pts_a, pts_a[:3]])           # ids 5,6,7 duplicate 0,1,2
    tets_b = np.array([[0, 1, 2, 3], [5, 7, 6, 4]], dtype=np.int64)
    model_b = SurfaceModel(points=pts_b, tets=tets_b,
                           surfaces={"fault": np.array([[0, 1, 2]])},
                           fault_names=["fault"])
    groups_b = duplicate_node_groups(pts_b)
    emb_b = fault_embedding(model_b, "fault")
    ok_b = (len(groups_b) == 3 and emb_b["n_interior"] == 0
            and emb_b["n_boundary"] == 1)
    print(f"self-test (b) duplicated nodes: groups={len(groups_b)} "
          f"emb={emb_b} -> {'PASS' if ok_b else 'FAIL'}")
    if ok_a and ok_b:
        print("self-test PASS")
        return 0
    print("self-test FAIL", file=sys.stderr)
    return 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _auto_pairs(model: SurfaceModel) -> list[tuple[str, str]]:
    pairs = []
    for n in sorted(model.surfaces):
        if "_minus" in n:
            sib = n.replace("_minus", "_plus")
            if sib in model.surfaces:
                pairs.append((n, sib))
    return pairs


def run_checks(path: Path, fault_tags: list[int] | None, dem_name: str | None,
               dem_tag: int | None, pairs_mode: str, tol: float) -> tuple[dict, int]:
    model = load_surfaces(path, fault_tags=fault_tags, dem_name=dem_name,
                          dem_tag=dem_tag)
    report: dict = {"mesh": str(path), "tol_m": tol,
                    "n_nodes": int(model.points.shape[0]),
                    "n_tets": int(model.tets.shape[0]) if model.tets is not None else 0,
                    "surfaces": {n: int(t.shape[0])
                                 for n, t in sorted(model.surfaces.items())},
                    "fault_names": model.fault_names,
                    "dem_name": model.dem_name}

    # Duplicates: mesh-wide + fault-node subset.
    groups_all = duplicate_node_groups(model.points, tol=tol)
    report["dup_groups_total"] = len(groups_all)
    report["dup_max_multiplicity"] = (max(len(g) for g in groups_all)
                                      if groups_all else 1)
    if model.fault_names:
        fault_ids = np.unique(np.concatenate(
            [model.surfaces[n].ravel() for n in model.fault_names]))
        groups_fault = duplicate_node_groups(model.points, ids=fault_ids, tol=tol)
        report["dup_groups_fault_subset"] = len(groups_fault)
        report["n_fault_surface_nodes"] = int(fault_ids.size)

    # Pair conformity.
    if pairs_mode == "auto":
        pair_list = _auto_pairs(model)
        report["pairs"] = {f"{a} | {b}": pair_conformity(model, a, b, tol)
                           for a, b in pair_list}

    # Embedding + slivers (volume meshes only).
    embedding_ok = True
    if model.tets is not None and model.fault_names:
        emb = {}
        for n in model.fault_names:
            emb[n] = fault_embedding(model, n)
            if emb[n]["n_interior"] != emb[n]["n_unique_tri"]:
                embedding_ok = False
        report["embedding"] = emb
    if model.tets is not None:
        report["slivers"] = sliver_report(model)

    rc = 0 if (len(groups_all) == 0 and embedding_ok) else 1
    report["pass"] = (rc == 0)
    return report, rc


def _print_report(r: dict) -> None:
    print(f"\n========== {Path(r['mesh']).name} ==========")
    print(f"nodes {r['n_nodes']:,} | tets {r['n_tets']:,} | "
          f"surfaces {len(r['surfaces'])}")
    print(f"duplicate groups (mesh-wide, {r['tol_m']} m): "
          f"{r['dup_groups_total']:,} "
          f"(max multiplicity {r['dup_max_multiplicity']})")
    if "dup_groups_fault_subset" in r:
        print(f"duplicate groups (fault-node subset): "
              f"{r['dup_groups_fault_subset']:,} of "
              f"{r['n_fault_surface_nodes']:,} fault nodes")
    if "pairs" in r and r["pairs"]:
        print("\n| pair | tri -/+ | nodes -/+ | shared IDs | coincident | "
              "dup pairs | centroid match |")
        print("|---|---|---|---:|---:|---:|---|")
        for k, p in r["pairs"].items():
            print(f"| {k[:48]} | {p['n_tri_a']:,}/{p['n_tri_b']:,} | "
                  f"{p['n_nodes_a']:,}/{p['n_nodes_b']:,} | "
                  f"{p['shared_node_ids']:,} | "
                  f"{p['coincident_coord_pairs']:,} | "
                  f"{p['duplicated_pairs']:,} | "
                  f"{p['centroid_matches']:,}/{p['n_tri_a']:,} |")
    if "embedding" in r:
        print("\n| fault | n_tri | interior(2) | boundary(1) | orphan(0) |")
        print("|---|---:|---:|---:|---:|")
        for n, e in r["embedding"].items():
            print(f"| {n[:48]} | {e['n_tri']:,} | {e['n_interior']:,} | "
                  f"{e['n_boundary']:,} | {e['n_orphan']:,} |")
    if "slivers" in r:
        s = r["slivers"]
        print(f"\nslivers (eta<={s['eta_floor']}): {s['n_slivers']:,}; "
              f"sub-{s['min_edge_floor_m']:.0f}m-edge tets: "
              f"{s['n_short_edge_tets']:,}")
    print(f"\nRESULT: {'PASS' if r['pass'] else 'FAIL'}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mesh", nargs="?", type=Path)
    ap.add_argument("--fault-tags", type=str, default=None,
                    help="comma-separated physical tags that are faults (.msh)")
    ap.add_argument("--dem-name", type=str, default=None)
    ap.add_argument("--dem-tag", type=int, default=None)
    ap.add_argument("--pairs", choices=("auto", "none"), default="auto")
    ap.add_argument("--tol", type=float, default=1e-3)
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)

    if args.self_test:
        return _self_test()
    if args.mesh is None:
        ap.error("MESH is required unless --self-test")

    fault_tags = ([int(x) for x in args.fault_tags.split(",")]
                  if args.fault_tags else None)
    report, rc = run_checks(args.mesh, fault_tags, args.dem_name,
                            args.dem_tag, args.pairs, args.tol)
    _print_report(report)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        with open(args.json, "w") as f:
            json.dump(report, f, indent=1, sort_keys=True)
            f.write("\n")
        print(f"wrote {args.json}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
