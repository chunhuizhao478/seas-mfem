"""msh_to_puml.py — convert a Gmsh v2.2 tet .msh (SAFS *_safstags convention:
fault=101, top=102, bottom=103, sides=104) to a SeisSol PUML/HDF5 mesh
(`*.puml.h5`) on macOS, WITHOUT PUMGen.

The PUML conventions here were reverse-engineered (decode_puml.py) from the
production reference `safv4_deep_500m_opt.puml.h5` and validated (0 unmatched
boundary faces):

  datasets   : geometry (Nnode,3) f8  | connect (Ntet,4) u8 (0-based)
               boundary (Ntet,) i4     | group (Ntet,) i4 (=1)
  file attrs : boundary-format='i32', topology-format='geometric'
  BC packing : byte i of `boundary` holds the BC code of tet face i:
                 boundary = sum_i code_i << (8*i)
  face order : SeisSol tet face -> local vertices
                 f0={0,2,1} f1={0,1,3} f2={1,2,3} f3={0,3,2}
  tag -> BC  : 101 fault -> 3 (dynamic rupture)
               102 top   -> 1 (free surface)
               103 bottom-> 5 (absorbing)
               104 sides -> 5 (absorbing)

Nodes/tets are kept in .msh order (SeisSol re-partitions at run time, so the
absolute ordering is irrelevant; PUMGen's own reordering is not reproduced and
is not needed).  HOWEVER each tet's local vertex order IS normalized to a
positive (right-handed) signed volume (orient_tets_positive): Gmsh does not
guarantee tet orientation, and inverted tets make SeisSol blow up
(bulk energy -> Inf/NaN at t~1s).  --validate asserts 0 inverted tets.

Usage:
    conda activate pythonenv
    export PYTHONPATH=.../project_7.0_preferred/meshing/code
    python msh_to_puml.py IN_safstags.msh OUT.puml.h5 [--bc 101:3,102:1,103:5,104:5] [--validate]
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np
import h5py
import meshio

# SeisSol tet face -> local vertex indices (decoded face-map "B")
SEISSOL_FACES = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)]
DEFAULT_BC = {101: 3, 102: 1, 103: 5, 104: 5}


def read_msh(path):
    m = meshio.read(str(path))
    pts = np.asarray(m.points, np.float64)
    phys = m.cell_data.get("gmsh:physical")
    tetb, trib, trir = [], [], []
    for bi, cb in enumerate(m.cells):
        if cb.type == "tetra":
            tetb.append(np.asarray(cb.data, np.int64))
        elif cb.type == "triangle":
            t = np.asarray(cb.data, np.int64)
            r = (np.asarray(phys[bi], np.int64) if phys is not None
                 else np.zeros(len(t), np.int64))
            trib.append(t); trir.append(r)
    tets = np.concatenate(tetb) if len(tetb) > 1 else tetb[0]
    surf = np.concatenate(trib) if trib else np.zeros((0, 3), np.int64)
    sref = np.concatenate(trir) if trir else np.zeros(0, np.int64)
    return pts, tets, surf, sref


def signed_volume(pts, tets):
    """Signed volume per tet (positive = right-handed, the SeisSol convention)."""
    p = pts[tets]                                   # (Ntet,4,3)
    return np.einsum("ij,ij->i",
                     np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                     p[:, 3] - p[:, 0]) / 6.0


def orient_tets_positive(pts, tets):
    """Flip inverted tets so ALL signed volumes are positive (SeisSol/ADER-DG
    requires a consistent right-handed orientation; the reference
    safv4_deep_500m_opt.puml.h5 is 100% positive).

    Gmsh/meshio do NOT normalize tet orientation, so a raw .msh typically has
    ~50% negative-volume tets.  Writing them verbatim makes SeisSol run past init
    and then BLOW UP (bulk elastic energy -> 1e271 -> "Inf/NaN in energies" abort
    at t~1s; observed on safv4_deep_500m_flattop_fixed_opt.puml.h5, 2026-06-25).

    Swapping local vertices 2<->3 negates the signed volume.  This MUST run BEFORE
    build_boundary so the boundary is rebuilt from the corrected connectivity (the
    boundary is keyed by sorted node-set, so each face keeps its BC code).  Returns
    (tets_oriented, n_flipped, n_degenerate)."""
    v = signed_volume(pts, tets)
    neg = v < 0.0
    degen = v == 0.0
    if neg.any():
        tets = tets.copy()
        tets[neg] = tets[neg][:, [0, 1, 3, 2]]
    return tets, int(neg.sum()), int(degen.sum())


def _sview(a):
    """(N,3) int64 -> 1-D structured view, one field per column.

    Field-wise lexicographic ordering -- what np.unique(axis=0) uses internally.
    This REPLACES the former `_tri_key` packed integer (a*nn*nn + b*nn + c),
    which OVERFLOWS int64 once nn > ~2.1e6.  At nn = 9,467,963 the largest key
    is 8.5e20 against an int64 max of 9.2e18, so keys wrapped silently and
    stopped being injective: two distinct faces could collide and swap BC codes.
    No injective 64-bit packing of a sorted triple exists at this node count, so
    the packed key is dropped in favour of a real lexicographic comparison.
    """
    a = np.ascontiguousarray(a)
    return a.view([(f"f{i}", a.dtype) for i in range(a.shape[1])]).ravel()


def _canon(tri_rows):
    """canonical (order-independent) sorted (N,3) int64 face rows."""
    return np.sort(np.asarray(tri_rows, np.int64), axis=1)


def _lookup(sorted_keys, query):
    """(hit mask, insertion positions) of `query` within sorted `sorted_keys`."""
    p = np.searchsorted(sorted_keys, query)
    hit = p < len(sorted_keys)
    if hit.any():
        hit[hit] = sorted_keys[p[hit]] == query[hit]
    return hit, p


def build_boundary(pts, tets, surf, sref, bc_map):
    """(boundary i4 array, stats dict).  byte i = BC code of tet face i."""
    boundary = np.zeros(len(tets), np.int64)
    per_tag = {tg: 0 for tg in set(bc_map.values())}
    n_bfaces = 0
    sv_sorted = codes_sorted = None
    if len(surf):
        sv = _sview(_canon(surf))
        codes_s = np.array([bc_map.get(int(r), 0) for r in sref], np.int64)
        order = np.argsort(sv, kind="stable")
        sv_sorted, codes_sorted = sv[order], codes_s[order]
        n_dup = len(sv_sorted) - len(np.unique(sv_sorted))
        if n_dup:
            print(f"  WARNING: {n_dup:,} duplicate surface triangles in the .msh")
        del order, sv
    for fi, F in enumerate(SEISSOL_FACES):
        codes = np.zeros(len(tets), np.int64)
        if sv_sorted is not None:
            fk = _sview(_canon(tets[:, list(F)]))
            hit, p = _lookup(sv_sorted, fk)
            codes[hit] = codes_sorted[p[hit]]
            del fk, hit, p
        boundary |= (codes << (8 * fi))
        nb = int((codes != 0).sum()); n_bfaces += nb
        for c in np.unique(codes):
            if c: per_tag[int(c)] = per_tag.get(int(c), 0) + int((codes == c).sum())
        del codes
    return boundary.astype(np.int32), {"n_boundary_faces": n_bfaces, "per_bc": per_tag}


def write_puml(path, pts, tets, boundary, group=None):
    if group is None:
        group = np.ones(len(tets), np.int32)
    with h5py.File(str(path), "w") as f:
        f.create_dataset("geometry", data=pts.astype(np.float64))
        f.create_dataset("connect", data=tets.astype(np.uint64))
        f.create_dataset("boundary", data=boundary.astype(np.int32))
        f.create_dataset("group", data=group.astype(np.int32))
        f.attrs["boundary-format"] = np.bytes_(b"i32")
        f.attrs["topology-format"] = np.bytes_(b"geometric")


def validate(puml_path, pts, tets, surf, sref, bc_map):
    """Decode the written boundary via the reference convention and confirm it
    reproduces the .msh surface-tri BCs exactly (per-tag face counts match)."""
    with h5py.File(str(puml_path), "r") as f:
        b = np.asarray(f["boundary"]).astype(np.int64)
        conn = np.asarray(f["connect"]).astype(np.int64)
        geom = np.asarray(f["geometry"])
        grp = np.asarray(f["group"])
        attrs = dict(f.attrs)
    ok = True
    msgs = []
    msgs.append(f"datasets: geometry{geom.shape} connect{conn.shape} "
                f"boundary{b.shape} group(unique={np.unique(grp).tolist()}) attrs={attrs}")
    # decode boundary -> per-BC face count, compare to expected per-tag tri counts
    decoded = {}
    for fi, F in enumerate(SEISSOL_FACES):
        codes = (b >> (8 * fi)) & 0xff
        for c in np.unique(codes):
            if c:
                decoded[int(c)] = decoded.get(int(c), 0) + int((codes == c).sum())
    # Expected coverage: a fault is an INTERNAL surface tri (shared by 2 tets) so
    # BOTH tets tag that face -> 2x; external tris (top/bottom/sides) -> 1x.
    # Multiplicity of each surface tri among ALL tet faces.  The former code used
    # a collections.Counter over every tet face: 4 x 54M = 216M Python int
    # increments into a ~110M-entry dict (>12 GB).  A sorted face array plus
    # searchsorted ranges is exact and memory-bounded.
    nt = len(tets)
    allF = np.empty((4 * nt, 3), np.int64)
    for fi, F in enumerate(SEISSOL_FACES):
        allF[fi * nt:(fi + 1) * nt] = tets[:, list(F)]
    allF.sort(axis=1)
    av = _sview(allF)
    av.sort()                       # in place: also reorders allF's rows
    sv = _sview(_canon(surf)) if len(surf) else None
    expected = {}
    if sv is not None:
        mult = (np.searchsorted(av, sv, side="right")
                - np.searchsorted(av, sv, side="left"))
        codes_s = np.array([bc_map.get(int(r), 0) for r in sref], np.int64)
        for c in np.unique(codes_s):
            if c:
                expected[int(c)] = int(mult[codes_s == c].sum())
        del mult, codes_s
    del av, allF
    msgs.append(f"decoded boundary face counts per BC: {decoded}")
    msgs.append(f"expected (internal faults x2)  per BC: {expected}")
    if decoded != expected:
        ok = False; msgs.append("  !! MISMATCH: boundary does not round-trip")
    else:
        msgs.append("  boundary round-trips exactly (fault internal => x2, "
                    "boundary external => x1)")
    # connectivity in range / 0-based
    nn = len(pts)
    if conn.max() >= nn or conn.min() < 0:
        ok = False; msgs.append(f"  !! connect out of range [0,{nn})")
    # orientation: EVERY tet must have positive signed volume (SeisSol convention);
    # inverted tets make SeisSol blow up (Inf/NaN energies) after init.
    vv = signed_volume(geom, conn)
    n_inv = int((vv <= 0.0).sum())
    msgs.append(f"inverted/zero-volume tets (signed vol <= 0): {n_inv} "
                f"(min vol {vv.min():.3e})")
    if n_inv:
        ok = False
        msgs.append("  !! MESH HAS INVERTED TETS \u2014 SeisSol will blow up "
                    "(bulk energy -> Inf/NaN). orient_tets_positive() should have "
                    "fixed this; do not ship.")
    del vv
    # every boundary face actually coincides with a surface tri (by node set)
    bad = 0
    if sv is not None:
        sv_sorted = np.sort(sv)
        for fi, F in enumerate(SEISSOL_FACES):
            codes = (b >> (8 * fi)) & 0xff
            nz = np.flatnonzero(codes != 0)
            if not len(nz):
                continue
            fk = _sview(_canon(conn[nz][:, list(F)]))
            hit, _ = _lookup(sv_sorted, fk)
            bad += int((~hit).sum())
            del fk, hit, nz
    msgs.append(f"boundary faces NOT matching a surface tri: {bad}")
    if bad: ok = False
    return ok, msgs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("infile", type=Path); ap.add_argument("outfile", type=Path)
    ap.add_argument("--bc", default=None,
                    help="override tag:code map, e.g. 101:3,102:1,103:5,104:5")
    ap.add_argument("--validate", action="store_true",
                    help="decode the output and verify it round-trips the .msh BCs")
    args = ap.parse_args(argv)
    bc_map = dict(DEFAULT_BC)
    if args.bc:
        bc_map = {int(k): int(v) for k, v in (p.split(":") for p in args.bc.split(","))}

    print(f"reading {args.infile} ...", flush=True)
    pts, tets, surf, sref = read_msh(args.infile)
    from collections import Counter
    print(f"  nodes={len(pts):,} tets={len(tets):,} surf={len(surf):,} "
          f"tags={dict(Counter(int(x) for x in sref))}")
    print(f"  tag->BC map: {bc_map}")
    # Normalize tet orientation BEFORE building the boundary (Gmsh does not
    # guarantee a consistent sign; inverted tets -> SeisSol bulk-energy blowup).
    tets, n_flipped, n_degen = orient_tets_positive(pts, tets)
    print(f"  orientation: flipped {n_flipped:,} inverted tets -> positive "
          f"({100.0 * n_flipped / max(len(tets), 1):.1f}%); "
          f"{n_degen:,} zero-volume (degenerate)")
    if n_degen:
        sys.exit(f"ERROR: {n_degen} zero-volume (degenerate) tets cannot be "
                 "oriented — repair the mesh before converting.")
    boundary, stats = build_boundary(pts, tets, surf, sref, bc_map)
    print(f"  boundary faces tagged: {stats['n_boundary_faces']:,}  per BC: {stats['per_bc']}")
    write_puml(args.outfile, pts, tets, boundary)
    print(f"wrote {args.outfile} (PUML/HDF5): {len(pts):,} nodes, {len(tets):,} tets")
    if args.validate:
        ok, msgs = validate(args.outfile, pts, tets, surf, sref, bc_map)
        print("\n=== validation ===")
        for m in msgs: print("  " + m)
        print(f"  RESULT: {'PASS' if ok else 'FAIL'}")
        return 0 if ok else 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
