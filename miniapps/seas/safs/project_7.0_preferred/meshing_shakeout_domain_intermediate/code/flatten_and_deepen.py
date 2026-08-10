#!/usr/bin/env python3
"""flatten_and_deepen.py -- snap the free surface to exactly z = 0 and move the
flat bottom to an exact depth, IN PLACE, on a SeisSol PUML mesh.

Only /geometry changes: connectivity, boundary codes and groups are untouched,
so every BC census and the fault TRIANGULATION are invariant by construction.

SAFETY.  Nothing is written until every check has passed on the NEW geometry:
tet orientation is re-evaluated for the whole mesh, the fault triangle multiset
is compared before/after, and the displacement of every moved vertex is
reported.  Only then is /geometry rewritten with h5py mode 'r+', which patches
the existing file without copying it.

THE FAULT-TRACE CAVEAT.  Some vertices belong to BOTH the free surface and the
fault -- they are the surface trace.  Snapping those MOVES THE FAULT SURFACE.
They are therefore excluded unless --include-fault-trace is given explicitly,
and the tool reports exactly which they are and how far they would move.
Measured on safpref_intermediate_shakeout: 818 lid vertices are off z = 0, of
which 808 are ordinary lid nodes (max |z| 14.7 m) and 10 are fault-trace nodes
(max |z| 49.92 m).

Usage:
    python flatten_and_deepen.py --mesh m.puml.h5 --z-bottom -40000 \
        [--include-fault-trace] [--apply]
"""

import argparse
import sys

import numpy as np
import h5py

LOCAL_FACES = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2))
CH = 2_000_000


def face_code(b, s):
    return ((np.ascontiguousarray(b, np.int32).view(np.uint32) >> np.uint32(8 * s))
            & np.uint32(0xFF)).astype(np.int32)


def vert_sets(path, P):
    """Vertex index sets for the free surface, the flat bottom, and the fault."""
    lid, bot, fau = set(), set(), set()
    zbot0 = float(P[:, 2].min())
    with h5py.File(path, "r") as f:
        B = f["boundary"][:].astype(np.int32)
        conn = f["connect"]
        nt = conn.shape[0]
        for s0 in range(0, nt, CH):
            s1 = min(s0 + CH, nt)
            T = conn[s0:s1].astype(np.int64)
            Bc = B[s0:s1]
            for s in range(4):
                fc = face_code(Bc, s)
                m = fc == 1
                if m.any():
                    lid.update(np.unique(T[m][:, list(LOCAL_FACES[s])]).tolist())
                m = fc == 3
                if m.any():
                    fau.update(np.unique(T[m][:, list(LOCAL_FACES[s])]).tolist())
                m = fc == 5
                if m.any():
                    tri = T[m][:, list(LOCAL_FACES[s])]
                    flat = np.all(np.abs(P[tri][:, :, 2] - zbot0) < 1.0, axis=1)
                    if flat.any():
                        bot.update(np.unique(tri[flat]).tolist())
            del T, Bc
    return (np.array(sorted(lid), np.int64), np.array(sorted(bot), np.int64),
            np.array(sorted(fau), np.int64), zbot0)


def fault_key(path, P):
    """Canonical coordinate multiset of the fault triangulation."""
    tri = []
    with h5py.File(path, "r") as f:
        B = f["boundary"][:].astype(np.int32)
        conn = f["connect"]
        nt = conn.shape[0]
        for s0 in range(0, nt, CH):
            s1 = min(s0 + CH, nt)
            T = conn[s0:s1].astype(np.int64)
            Bc = B[s0:s1]
            for s in range(4):
                m = face_code(Bc, s) == 3
                if m.any():
                    tri.append(T[m][:, list(LOCAL_FACES[s])])
            del T, Bc
    tri = np.vstack(tri)
    W = np.round(P[tri] * 1e3).astype(np.int64)
    i = np.lexsort((W[:, :, 2], W[:, :, 1], W[:, :, 0]), axis=1)
    W = np.take_along_axis(W, i[:, :, None], axis=1).reshape(len(W), 9)
    A = np.linalg.norm(np.cross(P[tri][:, 1] - P[tri][:, 0],
                                P[tri][:, 2] - P[tri][:, 0]), axis=1).sum() * 0.5
    return W[np.lexsort([W[:, k] for k in range(8, -1, -1)])], A


def n_inverted(path, P):
    neg = 0
    worst = np.inf
    with h5py.File(path, "r") as f:
        conn = f["connect"]
        nt = conn.shape[0]
        for s0 in range(0, nt, CH):
            T = conn[s0:min(s0 + CH, nt)].astype(np.int64)
            p = P[T]
            v = np.einsum("ij,ij->i", p[:, 1] - p[:, 0],
                          np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0])) / 6.0
            neg += int((v < 0).sum())
            worst = min(worst, float(np.abs(v).min()))
            del T, p, v
    return neg, worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--z-bottom", type=float, default=None)
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--band", type=float, default=10000.0,
                    help="thickness of the graded stretch band above the bottom")
    ap.add_argument("--include-fault-trace", action="store_true")
    ap.add_argument("--trace-safe", action="store_true",
                    help="with --include-fault-trace, snap only those trace vertices "
                         "whose incident tets stay valid; report the rest")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()

    with h5py.File(a.mesh, "r") as f:
        P0 = f["geometry"][:]
        nt = f["connect"].shape[0]
    print(f"[mesh] {a.mesh.split('/')[-1]}   {nt:,} tets, {len(P0):,} verts")
    print(f"       z {P0[:,2].min():,.4f} .. {P0[:,2].max():,.4f}")

    lid, bot, fau, zbot0 = vert_sets(a.mesh, P0)
    print(f"[sets] lid {len(lid):,}   bottom {len(bot):,} (z = {zbot0:,.4f})   "
          f"fault {len(fau):,}")

    P = P0.copy()

    # ---- free surface ------------------------------------------------------
    stray = lid[np.abs(P0[lid][:, 2]) > a.tol]
    on_fault = np.isin(stray, fau)
    move_lid = stray if a.include_fault_trace else stray[~on_fault]
    print(f"[lid]  strays {len(stray):,}: {int((~on_fault).sum()):,} ordinary, "
          f"{int(on_fault.sum()):,} on the FAULT TRACE")
    if len(stray[~on_fault]):
        print(f"       ordinary max |z| {np.abs(P0[stray[~on_fault]][:,2]).max():.4f} m")
    if int(on_fault.sum()):
        zz = P0[stray[on_fault]][:, 2]
        print(f"       fault-trace max |z| {np.abs(zz).max():.4f} m -- "
              f"{'INCLUDED (fault WILL move)' if a.include_fault_trace else 'EXCLUDED (fault preserved)'}")
    if a.include_fault_trace and a.trace_safe and int(on_fault.sum()):
        # Snapping a fault-trace vertex is NOT always free: the trace sits in a
        # thin wedge between the lid and the shallowest fault row, and pulling
        # it onto z = 0 can flip a neighbour.  Measured on
        # safpref_intermediate_shakeout: snapping all 10 inverted 2 tets.
        # So test each candidate against its own incident tets and keep only
        # the ones that are safe; the rest are reported, not forced.
        cand = stray[on_fault]
        inc = {int(v): [] for v in cand}
        with h5py.File(a.mesh, "r") as f:
            conn = f["connect"]
            nt2 = conn.shape[0]
            cs = set(int(v) for v in cand)
            for s0 in range(0, nt2, CH):
                T = conn[s0:min(s0 + CH, nt2)].astype(np.int64)
                hit = np.nonzero(np.isin(T, list(cs)).any(1))[0]
                for r in hit:
                    for v in T[r]:
                        if int(v) in inc:
                            inc[int(v)].append(T[r])
                del T
        ok, bad = [], []
        for v in cand:
            Tl = np.array(inc[int(v)], np.int64)
            if not len(Tl):
                continue
            z_old = P[v, 2]
            P[v, 2] = 0.0
            p = P[Tl]
            vol = np.einsum("ij,ij->i", p[:, 1] - p[:, 0],
                            np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0])) / 6.0
            p0 = P0[Tl]
            vol0 = np.einsum("ij,ij->i", p0[:, 1] - p0[:, 0],
                             np.cross(p0[:, 2] - p0[:, 0], p0[:, 3] - p0[:, 0])) / 6.0
            if np.all(np.sign(vol) == np.sign(vol0)) and np.abs(vol).min() > 0:
                ok.append(int(v))
            else:
                P[v, 2] = z_old
                bad.append((int(v), float(z_old)))
        print(f"       trace-safe: snapped {len(ok)} of {len(cand)}; "
              f"{len(bad)} would invert a neighbour and were LEFT IN PLACE")
        for vi, zz in sorted(bad, key=lambda t: -abs(t[1])):
            print(f"         kept  vert {vi:>9,}  z {zz:+10.4f} m  "
                  f"x {P0[vi,0]:11,.1f}  y {P0[vi,1]:13,.1f}")
        move_lid = stray[~on_fault]
        P[move_lid, 2] = 0.0
        print(f"       snapping {len(move_lid):,} ordinary + {len(ok)} trace vertices")
    else:
        P[move_lid, 2] = 0.0
        print(f"       snapping {len(move_lid):,} vertices to z = 0")

    # ---- bottom: GRADED stretch, not a rigid plane shift --------------------
    #
    # Translating only the bottom plane INVERTS tets: a cell with two vertices
    # on the plane and two above gets sheared through degeneracy.  Measured on
    # safpref_intermediate_shakeout, moving the plane -39,329.4 -> -40,000
    # inverted 15 of 41,985,157 tets.
    #
    # Instead map z through a monotone ramp over a band [z_bot, z_ref]:
    #     f(z) = z + (z_target - z_bot) * (z_ref - z)/(z_ref - z_bot)
    # f(z_bot) = z_target, f(z_ref) = z_ref, and
    #     f'(z) = 1 - (z_target - z_bot)/(z_ref - z_bot) > 0  for a deepening.
    # The map is (x, y, z) -> (x, y, f(z)), whose Jacobian is diag(1, 1, f'),
    # so det = f' > 0 EVERYWHERE: orientation is preserved by construction and
    # no tet can invert, wherever the band edge falls.
    #
    # z_ref is placed strictly BELOW the deepest fault vertex, so the fault
    # cannot move; that is asserted, not assumed.
    nb = 0
    if a.z_bottom is not None:
        if np.isin(bot, fau).any():
            raise RuntimeError("bottom shares vertices with the fault -- refusing")
        zf = float(P0[fau][:, 2].min()) if len(fau) else 0.0
        z_ref = zbot0 + a.band
        if z_ref >= zf - 1.0:
            z_ref = 0.5 * (zbot0 + zf)
            print(f"[bot]  band shortened so it stays below the fault "
                  f"(fault z_min {zf:,.1f})")
        H = z_ref - zbot0
        d = a.z_bottom - zbot0
        slope = 1.0 - d / H
        if slope <= 0:
            raise RuntimeError(f"ramp would fold (f' = {slope:.3f}); widen --band")
        m = P0[:, 2] <= z_ref
        P[m, 2] = P0[m, 2] + d * (z_ref - P0[m, 2]) / H
        nb = int(m.sum())
        print(f"[bot]  graded stretch over [{zbot0:,.1f}, {z_ref:,.1f}] "
              f"(band {H/1e3:.1f} km): {nb:,} vertices, f' = {slope:.4f}")
        print(f"       bottom {zbot0:,.4f} -> {a.z_bottom:,.4f} ({d:+,.4f} m); "
              f"deepest fault vertex z {zf:,.1f} is {zf - z_ref:,.1f} m above the band")
        if np.abs(P[fau, 2] - P0[fau, 2]).max() > 0:
            raise RuntimeError("the ramp moved fault vertices -- ABORT")

    moved = np.nonzero(np.any(P != P0, axis=1))[0]
    print(f"[move] {len(moved):,} vertices changed, "
          f"max |dz| {np.abs(P[moved,2]-P0[moved,2]).max() if len(moved) else 0:,.4f} m")

    # ---- checks on the NEW geometry, before anything is written -------------
    k0, A0 = fault_key(a.mesh, P0)
    k1, A1 = fault_key(a.mesh, P)
    same = k0.shape == k1.shape and np.array_equal(k0, k1)
    print(f"[chk]  fault triangle multiset identical: {same}")
    print(f"       fault area {A0/1e6:,.6f} -> {A1/1e6:,.6f} km2  "
          f"(delta {A1-A0:+.4e} m2)")
    if not a.include_fault_trace and not same:
        raise RuntimeError("fault changed although the trace was excluded -- ABORT")

    neg, vmin = n_inverted(a.mesh, P)
    print(f"[chk]  inverted tets after the move: {neg}   min |volume| {vmin:,.4f} m3")
    if neg:
        raise RuntimeError("the move inverts tets -- ABORT")

    zl = P[lid][:, 2]
    print(f"[chk]  lid now: |z|>{a.tol} on {int((np.abs(zl)>a.tol).sum()):,} of {len(lid):,}"
          f"   z {zl.min():.6f}..{zl.max():.6f}")
    print(f"[chk]  new z range {P[:,2].min():,.4f} .. {P[:,2].max():,.4f}")

    if not a.apply:
        print("\nDRY RUN -- nothing written.  Re-run with --apply.")
        return
    with h5py.File(a.mesh, "r+") as f:
        f["geometry"][:] = P
    print(f"\n[write] /geometry patched in place in {a.mesh}")


if __name__ == "__main__":
    main()
