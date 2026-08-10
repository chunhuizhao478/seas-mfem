#!/usr/bin/env python3
"""merge_collar.py -- weld a lateral collar onto its (frozen) parent PUML mesh.

Every parent vertex, tet and boundary code is carried over unchanged, so the
fault, the free surface, the gate compliance and the min-edge/dt floor come
across bit-for-bit rather than being re-derived.

Three edits happen:
  * shared WALL vertices are NOT duplicated -- collar tets are re-indexed onto
    the parent's own wall vertex ids, which is what makes the weld conformal;
  * the parent's old side-wall faces stop being absorbing and become interior;
  * the collar's outer faces are tagged: z = z_top -> free surface, everything
    else (new outer walls + new bottom) -> absorbing.  Its wall-side faces are
    left interior.

The wall is RELOCATED from the parent being merged, by exact coordinate key --
never by the vertex/tet indices stored in the collar npz.  That is what lets a
single collar serve both the intermediate and the heavy mesh: the two share a
bit-identical side wall (26,704 triangles, verified) but number their vertices
and tets completely differently.

Usage:
    python merge_collar.py --parent <p.puml.h5> --collar <collar.npz> \
        --out <merged.puml.h5> [--full-validate]
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from collar_lib import (BC_ABSORBING, BC_DYNAMIC_RUPTURE, BC_FREE_SURFACE,
                        LOCAL_FACES, extract_wall, face_code, tet_edge_lengths,
                        tet_eta, tet_signed_volume)

CH = 3_000_000


def key3(P):
    """Exact rounded coordinate key (0.1 mm) -- distinct mesh nodes are >= 1 m apart."""
    return np.round(np.asarray(P, float) * 1e4).astype(np.int64)


def boundary_faces(tets):
    """(face_vertices, owner_tet, local_slot) for faces used by exactly one tet."""
    nt = len(tets)
    t = tets.astype(np.int32) if tets.max() < 2 ** 31 - 1 else tets
    faces = np.concatenate([t[:, LOCAL_FACES[s]] for s in range(4)])
    owner = np.tile(np.arange(nt, dtype=np.int64), 4)
    slot = np.repeat(np.arange(4, dtype=np.int8), nt)
    key = np.sort(faces, axis=1)
    order = np.lexsort((key[:, 2], key[:, 1], key[:, 0]))
    ks = key[order]
    dup = np.all(ks[1:] == ks[:-1], axis=1)
    uniq = np.ones(len(ks), bool)
    uniq[:-1] &= ~dup
    uniq[1:] &= ~dup
    sel = order[uniq]
    return faces[sel].astype(np.int64), owner[sel], slot[sel]


def rows_in(a, b):
    """Which rows of a (sorted int triples) appear in b."""
    dt = [("i", a.dtype), ("j", a.dtype), ("k", a.dtype)]
    va = np.ascontiguousarray(a).view(dt).ravel()
    vb = np.ascontiguousarray(np.asarray(b, a.dtype)).view(dt).ravel()
    return np.isin(va, vb)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parent", required=True)
    ap.add_argument("--collar", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tol", type=float, default=1e-3)
    ap.add_argument("--full-validate", action="store_true")
    a = ap.parse_args()

    t0 = time.time()
    import h5py
    with h5py.File(a.parent, "r") as f:
        G = f["geometry"][:]
        C = f["connect"][:].astype(np.int64)
        B = f["boundary"][:].astype(np.int32)
        GR0 = int(f["group"][0]) if "group" in f else 1
    nv0, nt0 = len(G), len(C)
    print(f"[parent] {nt0:,} tets, {nv0:,} verts", flush=True)

    S = np.load(a.collar)
    SP, ST = S["points"], S["tets"]
    z_top, z_bot = float(S["z_top"]), float(S["z_bot"])
    print(f"[collar] {len(ST):,} tets, {len(SP):,} verts, "
          f"z_top {z_top:,.3f} z_bot {z_bot:,.3f}", flush=True)

    # --- relocate the wall on THIS parent, by coordinate --------------------
    PW, TW, wall_global, wall_tet, wall_slot = extract_wall(G, C, B)
    print(f"[wall]   relocated {len(TW):,} tris / {len(PW):,} verts on this parent")

    kp = key3(PW)
    kc = key3(SP[S["wall_collar_idx"]])
    if kp.shape != kc.shape:
        raise RuntimeError(f"wall vertex count differs: parent {len(kp):,} "
                           f"vs collar {len(kc):,} -- collar was built on a "
                           f"different wall")
    # map collar-wall order -> parent-wall order via the shared coordinate key
    op = np.lexsort((kp[:, 2], kp[:, 1], kp[:, 0]))
    oc = np.lexsort((kc[:, 2], kc[:, 1], kc[:, 0]))
    if not np.array_equal(kp[op], kc[oc]):
        raise RuntimeError("wall vertex COORDINATES differ between parent and collar")
    parent_of_collarwall = np.empty(len(kp), np.int64)
    parent_of_collarwall[oc] = wall_global[op]

    c2m = np.full(len(SP), -1, np.int64)
    c2m[S["wall_collar_idx"]] = parent_of_collarwall
    fresh = np.flatnonzero(c2m < 0)
    c2m[fresh] = nv0 + np.arange(len(fresh))
    if not np.allclose(SP[S["wall_collar_idx"]], G[parent_of_collarwall], atol=1e-6):
        raise RuntimeError("wall vertex remap is inconsistent")
    print(f"  shared wall verts {len(kp):,}; new verts {len(fresh):,}")

    # --- boundary word ------------------------------------------------------
    boundary = B.copy()
    bu = boundary.view(np.uint32)
    for s in range(4):
        sel = wall_tet[wall_slot == s]
        if len(sel):
            bu[sel] &= np.uint32(~(0xFF << (8 * s)) & 0xFFFFFFFF)
    print(f"  parent wall faces cleared to interior: {len(wall_tet):,}")

    LF = np.asarray(LOCAL_FACES)
    wall_faces = np.sort(C[wall_tet[:, None], LF[wall_slot]], axis=1)

    sfaces, sowner, sslot = boundary_faces(ST)
    sfaces_m = np.sort(c2m[sfaces], axis=1)
    on_wall = rows_in(sfaces_m, wall_faces)

    x = sfaces_m[on_wall]
    y = wall_faces
    x = x[np.lexsort((x[:, 2], x[:, 1], x[:, 0]))]
    y = y[np.lexsort((y[:, 2], y[:, 1], y[:, 0]))]
    if x.shape != y.shape or not np.array_equal(x, y):
        raise RuntimeError(f"weld mismatch: collar wall faces {x.shape} vs parent {y.shape}")
    print(f"  WELD OK: {len(x):,} interface faces identical on both sides")

    ztri = SP[sfaces][:, :, 2]
    is_top = np.all(np.abs(ztri - z_top) < a.tol, axis=1) & ~on_wall
    is_out = (~on_wall) & (~is_top)
    sword = np.zeros(len(ST), np.int32)
    for s in range(4):
        m = is_top & (sslot == s)
        if m.any():
            sword[sowner[m]] |= np.int32(BC_FREE_SURFACE << (8 * s))
        m = is_out & (sslot == s)
        if m.any():
            sword[sowner[m]] |= np.int32(BC_ABSORBING << (8 * s))
    n_wall_f, n_top_f, n_out_f = int(on_wall.sum()), int(is_top.sum()), int(is_out.sum())
    print(f"  collar faces: wall(interior) {n_wall_f:,}, free surface {n_top_f:,}, "
          f"absorbing {n_out_f:,}")

    # --- stream the merged mesh out (never materialise a 2nd full connect) ---
    ST_m = c2m[ST]
    del c2m
    nvt, ntt = nv0 + len(fresh), nt0 + len(ST)
    with h5py.File(a.out, "w") as f:
        dg = f.create_dataset("geometry", (nvt, 3), np.float64)
        dg[:nv0] = G
        dg[nv0:] = SP[fresh]
        dc = f.create_dataset("connect", (ntt, 4), np.uint64)
        for s0 in range(0, nt0, CH):
            s1 = min(s0 + CH, nt0)          # clamp: the dataset is LONGER than
            dc[s0:s1] = C[s0:s1]            # the parent block, so s0+CH would
                                            # select past it and fail to broadcast
        dc[nt0:] = ST_m
        db = f.create_dataset("boundary", (ntt,), np.int32)
        db[:nt0] = boundary
        db[nt0:] = sword
        dgr = f.create_dataset("group", (ntt,), np.int32)
        dgr[:] = GR0

    # --- checks on the merged result ----------------------------------------
    boundary_all = np.concatenate([boundary, sword])
    ft = sum(int((face_code(boundary_all, s) == BC_DYNAMIC_RUPTURE).sum()) for s in range(4))
    fs = sum(int((face_code(boundary_all, s) == BC_FREE_SURFACE).sum()) for s in range(4))
    ab = sum(int((face_code(boundary_all, s) == BC_ABSORBING).sum()) for s in range(4))
    print(f"  BC census: fault {ft:,} (=2x{ft//2:,}), free surface {fs:,}, absorbing {ab:,}")

    # hull identity: trusting the parent's own validated tagging, the merged
    # hull is (parent hull - wall) + (collar hull - wall).  O(1) instead of a
    # 580M-face global sort on the heavy mesh.
    tagged_p = sum(int((face_code(B, s) != 0).sum()) for s in range(4))
    ft_p = sum(int((face_code(B, s) == BC_DYNAMIC_RUPTURE).sum()) for s in range(4))
    hull_p = tagged_p - ft_p            # fault faces are tagged on BOTH sides
    hull_expect = (hull_p - n_wall_f) + (len(sfaces) - n_wall_f)
    tagged_m = sum(int((face_code(boundary_all, s) != 0).sum()) for s in range(4))
    ok = tagged_m == hull_expect + ft
    print(f"  hull identity: tagged {tagged_m:,} == hull {hull_expect:,} + fault {ft:,} "
          f"-> {'OK' if ok else 'MISMATCH'}")
    if not ok:
        raise RuntimeError("boundary tagging does not match the expected hull")

    # Only the COLLAR block can be inverted -- the parent's tets are copied
    # verbatim and were validated when it was built.
    Gm = np.concatenate([G, SP[fresh]])
    neg = 0
    for s0 in range(0, len(ST_m), CH):
        neg += int((tet_signed_volume(Gm, ST_m[s0:s0 + CH]) < 0).sum())
    print(f"  collar signed volume: {neg} negative")
    if neg:
        raise RuntimeError("merged mesh has inverted tets")

    # Chunked, and only on the collar block: Gm[ST_m] on a 25M-tet collar is a
    # 2.4 GB temporary, and on the heavy parent the int64 connect is still
    # holding 3.9 GB at this point.
    del C
    eta_min, eta_med, n05, n01, e_min, e_max = np.inf, [], 0, 0, np.inf, 0.0
    for s0 in range(0, len(ST_m), CH):
        Tc = ST_m[s0:s0 + CH]
        et = tet_eta(Gm, Tc)
        ed = tet_edge_lengths(Gm, Tc)
        eta_min = min(eta_min, float(et.min()))
        n05 += int((et < 0.05).sum())
        n01 += int((et < 0.1).sum())
        e_min = min(e_min, float(ed.min()))
        e_max = max(e_max, float(ed.max()))
        eta_med.append(et[::13])
    print(f"  merged: {ntt:,} tets, {nvt:,} verts")
    print(f"    x [{Gm[:,0].min():,.1f}, {Gm[:,0].max():,.1f}]  "
          f"y [{Gm[:,1].min():,.1f}, {Gm[:,1].max():,.1f}]  "
          f"z [{Gm[:,2].min():,.1f}, {Gm[:,2].max():,.1f}]")
    print(f"  collar quality: eta min {eta_min:.4f} "
          f"med {np.median(np.concatenate(eta_med)):.3f}, "
          f"eta<0.05 {n05:,}, eta<0.1 {n01:,}")
    print(f"  collar edges: min {e_min:,.0f} max {e_max:,.0f} m")
    print(f"[write] {a.out}  ({time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
