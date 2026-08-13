#!/usr/bin/env python3
"""relax_rv.py -- bulk across-fault volume-ratio (rv) relaxation, scaled to 10^8 cells.

rv = max(V_A/V_B, V_B/V_A) over a dynamic-rupture face = the ratio of the two APEX
HEIGHTS, so this is a vertex-position problem and the fault triangulation is never
touched (Zhang et al. 2023 JGR Eq.18; high rv -> spurious tensile normal stress ->
strength clamped to 0 -> runaway slip).

Port of `meshing_deep19km/code/fix_rv_apex_relax.py`, which is the pass that did the
bulk of the work there (deep rv>2 cells 20,859 -> 1,131).  The key property, and the
reason a greedy per-vertex polish cannot substitute for it: every movable apex is
displaced SIMULTANEOUSLY toward the geometric mean sqrt(hA*hB), guarded only
GEOMETRICALLY, with a global backtracking line search that halves the step of
vertices sitting in a violating tet.  There is no per-move rv acceptance test --
an apex serves ~7 fault faces with 3 degrees of freedom, so almost no individual
move improves every face it touches, and a greedy rule rejects nearly all of them
(measured on this mesh: 147k accepted moves for a 4 % drop).

Two changes from the original, both required at this scale:
  * every whole-mesh gather is chunked or restricted to the WORKING SET (tets
    incident to a movable apex, ~2 % of the mesh);
  * a FREQUENCY-GATE guard: these meshes sit exactly on the resolved-frequency
    gate, so no tet's MAX EDGE may grow past Vs_pooled/gate (or its current size,
    whichever is larger).  Without it the relaxation silently re-opens the gate.

Pinned (BC 1 free surface / 3 fault / 5 absorbing) vertices never move, so the free
surface and the DR surface come out bit-identical.
"""
import argparse, sys, time
from pathlib import Path
import numpy as np, h5py

sys.path.insert(0, str(Path(__file__).resolve().parent))
from material import Material

FACE = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)]
TET_EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
BC_PINNED = (1, 3, 5)
CH = 3_000_000
BOX = (132000.0, 724000.0, 3490000.0, 4055000.0)


def svol(p):
    return np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                     p[:, 3] - p[:, 0]) / 6.0

def eta_of(p):
    v = np.abs(svol(p)); s = np.zeros(len(p))
    for a_, b_ in TET_EDGES:
        s += np.sum((p[:, a_] - p[:, b_]) ** 2, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(s > 0, 12.0 * np.cbrt((3.0 * v) ** 2) / s, 0.0)

def edges_of(p):
    return np.stack([np.linalg.norm(p[:, a_] - p[:, b_], axis=1) for a_, b_ in TET_EDGES], 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--cvm"); ap.add_argument("--muscal")
    ap.add_argument("--target-field", default=None,
                    help="pre-graded target size field. REQUIRED when the mesh was "
                         "refined against one: the apex-move guard must bound against "
                         "the ACTUAL target, not a raw Vs/gate bound, or it either "
                         "over-constrains or silently under-protects.")
    ap.add_argument("--gate", type=float, required=True)
    ap.add_argument("--rv-target", type=float, default=2.0)
    ap.add_argument("--iters", type=int, default=120)
    ap.add_argument("--omega", type=float, default=0.6)
    ap.add_argument("--hold-weight", type=float, default=1.0)
    ap.add_argument("--disp-cap", type=float, default=0.25, help="x local min edge, per iter")
    ap.add_argument("--cum-cap", type=float, default=3.0, help="x the per-iter cap, cumulative")
    ap.add_argument("--pool-frac", type=float, default=0.25)
    ap.add_argument("--qual-tol", type=float, default=0.05,
                    help="max fractional loss of a tet's OWN eta and inradius")
    a = ap.parse_args()
    t0 = time.time()

    with h5py.File(a.mesh) as f:
        x = f["geometry"][:].astype(np.float64)
        conn = f["connect"][:].astype(np.int64)
        bnd = f["boundary"][:].astype(np.int64)
    nV, nT = len(x), len(conn)
    x0 = x.copy()
    print(f"[mesh] {nT:,} tets  {nV:,} verts", flush=True)

    tri, own = [], []
    for k in range(4):
        m = ((bnd >> (8 * k)) & 0xFF) == 3
        if m.any():
            idx = np.flatnonzero(m)
            tri.append(conn[idx][:, list(FACE[k])]); own.append(idx)
    tri = np.vstack(tri); own = np.concatenate(own)
    key = np.sort(tri, axis=1)
    o = np.lexsort((key[:, 2], key[:, 1], key[:, 0]))
    key, own, tri = key[o], own[o], tri[o]
    assert np.all(key[0::2] == key[1::2])
    pairs = np.stack([own[0::2], own[1::2]], 1); tri_u = tri[0::2]
    nF = len(tri_u)
    apA = conn[pairs[:, 0]].sum(1) - tri_u.sum(1)
    apB = conn[pairs[:, 1]].sum(1) - tri_u.sum(1)

    v0 = x0[tri_u[:, 0]]
    nrm = np.cross(x0[tri_u[:, 1]] - v0, x0[tri_u[:, 2]] - v0)
    nrm /= np.linalg.norm(nrm, axis=1, keepdims=True)
    c0 = np.einsum("ij,ij->i", nrm, v0)

    def sAB(xc):
        return (np.einsum("fj,fj->f", nrm, xc[apA]) - c0,
                np.einsum("fj,fj->f", nrm, xc[apB]) - c0)
    sA, sB = sAB(x)
    rv0 = np.maximum(np.abs(sA) / np.abs(sB), np.abs(sB) / np.abs(sA))
    print(f"[rv] {nF:,} faces  med {np.median(rv0):.4f}  max {rv0.max():.3f}  "
          f">1.5 {int((rv0>1.5).sum()):,}  >2 {int((rv0>2).sum()):,}  >3 {int((rv0>3).sum()):,}", flush=True)

    pinned = np.zeros(nV, bool)
    for k in range(4):
        fc = (bnd >> (8 * k)) & 0xFF
        for bc in BC_PINNED:
            m = fc == bc
            if m.any():
                pinned[np.unique(conn[m][:, list(FACE[k])])] = True
    movA, movB = ~pinned[apA], ~pinned[apB]

    target = rv0 > a.rv_target
    movable = np.unique(np.concatenate([apA[target][movA[target]], apB[target][movB[target]]]))
    vid = -np.ones(nV, np.int64); vid[movable] = np.arange(len(movable))
    # any face touching a movable apex is AFFECTED and must be accounted for
    affected = (vid[apA] >= 0) | (vid[apB] >= 0)
    both = affected & movA & movB
    onlyA = affected & movA & ~movB
    onlyB = affected & movB & ~movA
    used = target & (both | onlyA | onlyB)
    hold = affected & ~target                 # healthy faces: anchor, do not drive
    print(f"[relax] target {int(target.sum()):,}  movable apexes {len(movable):,}  "
          f"affected faces {int(affected.sum()):,}  driven {int(used.sum()):,}  "
          f"hold {int(hold.sum()):,}", flush=True)

    amask = np.zeros(nV, bool); amask[movable] = True
    inc = []
    for s in range(0, nT, CH):
        c = conn[s:s + CH]
        inc.append(np.flatnonzero(amask[c].any(1)) + s)
    inc = np.concatenate(inc)
    print(f"[working] {len(inc):,} tets ({100*len(inc)/nT:.2f} %)", flush=True)

    eta_floor, edge_floor = np.inf, np.inf
    for s in range(0, nT, CH):
        p = x0[conn[s:s + CH]]
        eta_floor = min(eta_floor, float(eta_of(p).min()))
        edge_floor = min(edge_floor, float(edges_of(p).min()))
        del p
    ci = conn[inc]
    p0 = x0[ci]
    vol0 = svol(p0); sign0 = np.sign(vol0); vmin = 1e-6 * np.abs(vol0)
    # PER-TET guards.  The inherited guard compared every tet against the mesh's
    # GLOBAL worst eta, which on this mesh is 0.0307 (a fault-band sliver) -- so a
    # healthy tet could fall from 0.85 to 0.031 and still pass.  Measured cost of
    # that: eta<0.05 8 -> 895, eta<0.1 41 -> 1,260, and r_insphere 0.1649 ->
    # 0.0681 m, i.e. dt 16.129 -> 1.487 us.  Bound each tet against ITSELF.
    eta0 = eta_of(p0)
    def insphere(pp):
        v = np.abs(svol(pp)); A = np.zeros(len(pp))
        for f_ in FACE:
            q = pp[:, list(f_)]
            A += 0.5*np.linalg.norm(np.cross(q[:,1]-q[:,0], q[:,2]-q[:,0]), axis=1)
        return 3.0*v/A
    r0 = insphere(p0)
    eta_min_tet = np.maximum(eta0 * (1.0 - a.qual_tol), eta_floor)
    r_min_tet = r0 * (1.0 - a.qual_tol)
    dx0 = edges_of(p0).max(1)
    if a.target_field:
        from leb_gate_close import FieldTarget
        mat = FieldTarget(a.target_field, a.gate)
    elif a.cvm:
        mat = Material(a.cvm, a.muscal, box=BOX, source="muscal")
    else:
        raise SystemExit("need --cvm or --target-field")
    # (1+1e-6): where allow == dx0 exactly the guard is knife-edge and ANY
    # floating-point growth trips it, so the backtracking never converges and the
    # whole iteration is discarded -- measured as 0 vertices displaced.
    allow = np.maximum(mat.pooled(p0, dx0, a.pool_frac) / a.gate, dx0) * (1.0 + 1e-6)
    print(f"[floors] eta {eta_floor:.6f}  min_edge {edge_floor:.4f} m   "
          f"gate-frozen tets {int((allow <= dx0 + 1e-9).sum()):,}", flush=True)
    del p0

    tmp = np.full(len(movable), np.inf)
    e_min = edges_of(x0[ci]).min(1)
    for j in range(4):
        l = vid[ci[:, j]]
        m = l >= 0
        np.minimum.at(tmp, l[m], e_min[m])
    cap = a.disp_cap * tmp
    cum = a.cum_cap * cap

    nstall = nfail = nrej = 0
    n_severe = int((rv0 > 3.0).sum()); max_rv0 = float(rv0.max())
    for it in range(a.iters):
        sA, sB = sAB(x)
        hA, hB = np.abs(sA), np.abs(sB)
        hsA = np.where(both, np.sqrt(hA * hB), hB)
        hsB = np.where(both, np.sqrt(hA * hB), hA)
        dA = (hsA - hA)[:, None] * (np.sign(sA)[:, None] * nrm)
        dB = (hsB - hB)[:, None] * (np.sign(sB)[:, None] * nrm)
        acc = np.zeros((len(movable), 3)); cnt = np.zeros(len(movable))
        for apx, d, sel in ((apA, dA, used & movA), (apB, dB, used & movB)):
            l = vid[apx[sel]]; m = l >= 0
            np.add.at(acc, l[m], d[sel][m]); np.add.at(cnt, l[m], 1.0)
        for apx, sel in ((apA, hold & movA), (apB, hold & movB)):
            l = vid[apx[sel]]; m = l >= 0
            np.add.at(cnt, l[m], a.hold_weight)
        good = cnt > 0
        disp = np.zeros_like(acc)
        disp[good] = a.omega * acc[good] / cnt[good, None]
        n_ = np.linalg.norm(disp, axis=1)
        disp *= np.where(n_ > cap, cap / np.maximum(n_, 1e-30), 1.0)[:, None]
        new = (x[movable] - x0[movable]) + disp
        nn = np.linalg.norm(new, axis=1)
        ov = nn > cum
        if ov.any():
            disp[ov] = new[ov] / nn[ov, None] * cum[ov, None] - (x[movable] - x0[movable])[ov]

        active = np.ones(len(movable), bool)
        xt = x
        for _ in range(30):
            xt = x.copy(); xt[movable] += disp * active[:, None]
            pt = xt[ci]
            e = edges_of(pt); vt = svol(pt)
            bad = ((np.sign(vt) != sign0) | (np.abs(vt) < vmin)
                   | (eta_of(pt) < eta_min_tet) | (insphere(pt) < r_min_tet)
                   | (e.min(1) < edge_floor) | (e.max(1) > allow))
            if not bad.any():
                break
            bl = vid[np.unique(ci[bad].ravel())]; bl = bl[bl >= 0]
            tiny = np.linalg.norm(disp[bl], axis=1) < 1e-9
            active[bl[tiny]] = False
            disp[bl[~tiny]] *= 0.5
            nstall += 1
        else:
            # Do NOT discard the whole step because a few tets still object --
            # at 163k movable vertices one stubborn tet would veto everything and
            # nothing ever moves.  Freeze only the offenders and let the rest go.
            xt = x.copy(); xt[movable] += disp * active[:, None]
            pt = xt[ci]; e = edges_of(pt); vt = svol(pt)
            bad = ((np.sign(vt) != sign0) | (np.abs(vt) < vmin)
                   | (eta_of(pt) < eta_min_tet) | (insphere(pt) < r_min_tet)
                   | (e.min(1) < edge_floor) | (e.max(1) > allow))
            if bad.any():
                bl = vid[np.unique(ci[bad].ravel())]; bl = bl[bl >= 0]
                active[bl] = False
            xt = x.copy(); xt[movable] += disp * active[:, None]
            nfail += 1
        # GLOBAL tail guard: relaxing the rv>2 bulk can push individual faces UP
        # (measured on the small mesh: rv>3 112 -> 125).  Accept the iteration
        # only if the SEVERE tail does not grow -- rv 3 is the empirical risk
        # cliff (48.8x deep tensile-flip risk above it).  Otherwise halve and retry.
        for _try in range(6):
            sa, sb = sAB(xt)
            rt = np.maximum(np.abs(sa)/np.abs(sb), np.abs(sb)/np.abs(sa))
            if int((rt > 3.0).sum()) <= n_severe and rt.max() <= max_rv0 + 1e-9:
                break
            disp *= 0.5
            xt = x.copy(); xt[movable] += disp * active[:, None]
        else:
            xt = x.copy(); nrej += 1
        x = xt
        if it % 20 == 0 or it == a.iters - 1:
            sA, sB = sAB(x)
            r = np.maximum(np.abs(sA)/np.abs(sB), np.abs(sB)/np.abs(sA))
            print(f"  iter {it:3d}: >2 {int((r>2).sum()):,}  >3 {int((r>3).sum()):,}  "
                  f"max {r.max():.3f}  med {np.median(r):.4f}  "
                  f"[no-step {nfail}, tail-rejected {nrej}]  ({time.time()-t0:.0f} s)", flush=True)

    sA, sB = sAB(x)
    rv1 = np.maximum(np.abs(sA)/np.abs(sB), np.abs(sB)/np.abs(sA))
    mv = np.flatnonzero(np.any(x != x0, axis=1))
    dmax = float(np.linalg.norm(x[mv]-x0[mv],axis=1).max()) if len(mv) else 0.0
    print(f"\n[out] {len(mv):,} vertices displaced (max {dmax:.3f} m), {nfail} iterations took no step")
    for t in (1.5, 2.0, 3.0, 5.0):
        print(f"  rv>{t}: {int((rv0>t).sum()):>9,} -> {int((rv1>t).sum()):>9,}")
    print(f"  max {rv0.max():.3f} -> {rv1.max():.3f}   med {np.median(rv0):.4f} -> {np.median(rv1):.4f}")
    assert not pinned[mv].any(), "a PINNED vertex moved"
    with h5py.File(a.mesh) as fin, h5py.File(a.out, "w") as fo:
        fo.create_dataset("geometry", data=x)
        for k in ("connect", "boundary", "group"):
            if k in fin: fo.create_dataset(k, data=fin[k][:])
    print(f"[write] {a.out}  ({time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
