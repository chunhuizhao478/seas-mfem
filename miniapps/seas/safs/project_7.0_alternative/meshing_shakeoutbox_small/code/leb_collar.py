#!/usr/bin/env python3
"""leb_collar.py -- close the volume frequency gate INSIDE a collar, by
conforming 3-D longest-edge (Rivara) bisection.

Adapted from meshing_deep40km_1Hz_p5_leb/code/leb_refine_1Hz.py, which closed
the 1 Hz @ p5 gate on a 122M-cell mesh.  Two things differ:

  * it runs on the COLLAR ALONE (an npz of points+tets), before the weld.
    Running it on the merged mesh would be wrong: bisecting a collar cell whose
    longest edge lies ON the wall forces the parent's adjacent cell to split
    too, and the parent must stay bit-identical.
  * the FROZEN edge set is therefore the WALL triangulation, not the fault.
    Wall edges are never bisected, so the interface triangulation the weld keys
    on is preserved exactly.  Cells whose longest edge IS a wall edge form a
    structural exempt class, reported and not chased -- the direct analogue of
    the fault-edge exemption.

Everything else is carried over because it was hard-won:
  * the refinement TARGET is separated from the acceptance GATE, and is
    selectable (--target), because no single choice is right at every cell
    size.  The 1 Hz p5 build pooled Vs in z ONLY, valid there because its cells
    were <= ~500 m -- narrower than the CVM's 1500 m lateral lattice, so a
    refined child changed bin essentially only vertically.  Collar cells are
    1-8 km across and change bin LATERALLY too, and z-only pooling drove the
    small collar's gate BACKWARDS: 15 failing cells -> 82, worst 0.1314 ->
    0.0801, for +186 % tets.  Full bounding-box pooling (--target bbox) is a
    strict lower bound on every descendant and therefore provably convergent,
    but on km-scale cells it takes the slowest Vs anywhere inside and
    over-refines hard (+339 % on the same smoke test).  The default "half"
    pools over half the cell's extent -- the region a child's barycentre can
    reach in ONE bisection -- and is re-evaluated every round.
  * the LEPP tie-break is a STRICT TOTAL ORDER (length, then the globally
    unique edge key).  Breaking ties by local slot lets adjacent tets name each
    other's edge, the chain closes into a cycle, and no terminal edge is found.
  * EVERY geometry pass is chunked; P[T[:, EI]] on a 10^8-cell mesh is a
    17.6 GB temporary that swap-thrashes the machine.

Usage:
    python leb_collar.py --collar in.npz --parent p.puml.h5 --cvm cvm.nc \
        --out out.npz --gate 0.6667 [--hops 5] [--max-rounds 40]
"""
import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from collar_lib import (extract_wall, tet_edge_lengths, tet_eta,
                        tet_signed_volume)

PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
EI = np.array([p[0] for p in PAIRS])
EJ = np.array([p[1] for p in PAIRS])
NVMAX = np.int64(1) << np.int64(28)
CH = 3_000_000


def rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1 << 30)


class PoolVs:
    """min Vs over a cell's own BOUNDING BOX -- a true lower bound on any child.

    WHY NOT z-only pooling.  The 1 Hz p5 build pooled in z alone, on the
    argument that the CVM is 1500 m laterally but 250 m vertically, so a
    refined cell's barycentre almost always changes bin VERTICALLY.  That
    argument holds for cells <= ~500 m across.  It fails here: a collar cell is
    1-8 km across, far wider than the 1500 m lateral lattice, so children
    routinely land in a different lateral bin.  Measured with z-only pooling on
    the small collar: the gate went BACKWARDS, 15 failing cells -> 82, worst
    0.1314 -> 0.0801, while the mesh grew +186 %.  That is the treadmill the
    pooled target exists to prevent, just in x/y instead of z.

    Pooling over the cell's own bbox restores the guarantee: a child's bbox is
    contained in its parent's, so the child's pooled Vs is >= the parent's,
    while its dx strictly shrinks.  f is therefore monotone under bisection and
    the loop provably converges.

    Implementation: cells are bucketed by dyadic half-extent in CVM index units
    and each bucket is answered from a cached scipy minimum_filter.  A filter of
    size 2k+1 centred on the cell covers [c-k, c+k], which contains the bbox by
    construction, so the bound is conservative, never optimistic.
    """

    def __init__(self, path, max_cache=8):
        import h5py
        with h5py.File(path, "r") as f:
            self.X, self.Y, self.Z = f["x"][:], f["y"][:], f["z"][:]
            d = f["data"]
            mu = np.asarray(d["mu"]).astype(np.float64)
            rho = np.asarray(d["rho"]).astype(np.float64)
        self.vs = np.sqrt(np.maximum(mu, 0.0) / rho).astype(np.float32)
        self.dx = float(self.X[1] - self.X[0])
        self.dy = float(self.Y[1] - self.Y[0])
        self.dz = float(self.Z[1] - self.Z[0])
        self.n = self.vs.shape                       # (nz, ny, nx)
        self._cache = {}
        self._max = max_cache

    def _idx(self, P):
        ix = np.clip(np.rint((P[:, 0] - self.X[0]) / self.dx), 0, self.n[2] - 1).astype(np.int32)
        iy = np.clip(np.rint((P[:, 1] - self.Y[0]) / self.dy), 0, self.n[1] - 1).astype(np.int32)
        iz = np.clip(np.rint((P[:, 2] - self.Z[0]) / self.dz), 0, self.n[0] - 1).astype(np.int32)
        return ix, iy, iz

    def at(self, P):
        ix, iy, iz = self._idx(P)
        return self.vs[iz, iy, ix].astype(np.float64)

    def _filt(self, kx, ky, kz):
        key = (kx, ky, kz)
        if key not in self._cache:
            from scipy.ndimage import minimum_filter
            if len(self._cache) >= self._max:
                self._cache.pop(next(iter(self._cache)))
            self._cache[key] = (self.vs if (kx == 0 and ky == 0 and kz == 0)
                                else minimum_filter(self.vs,
                                                    size=(2 * kz + 1, 2 * ky + 1, 2 * kx + 1),
                                                    mode="nearest"))
        return self._cache[key]

    @staticmethod
    def _dyadic(k):
        k = np.maximum(k, 0)
        out = np.zeros_like(k)
        nz = k > 0
        out[nz] = 1 << np.ceil(np.log2(k[nz])).astype(np.int32)
        return out

    def bbox_min(self, lo, hi):
        """min Vs over each axis-aligned box [lo, hi].

        PERFORMANCE NOTE (measured, deliberately NOT optimised): fully
        anisotropic dyadic bucketing can generate ~100 distinct windows and each
        is a 127 MB scipy minimum_filter over the CVM cube, so the 8-entry cache
        can thrash.  On the deep collar it does not: a round over 8.8M cells
        costs ~118 s end to end.  Collapsing the buckets (kx = ky, capped) would
        speed it up but CHANGES the pooled value and therefore the mesh, so it
        is left alone -- the shipped collar must be reproducible from this file.
        """
        i0x, i0y, i0z = self._idx(lo)
        i1x, i1y, i1z = self._idx(hi)
        cx, cy, cz = (i0x + i1x) // 2, (i0y + i1y) // 2, (i0z + i1z) // 2
        kx = self._dyadic(np.maximum(cx - i0x, i1x - cx))
        ky = self._dyadic(np.maximum(cy - i0y, i1y - cy))
        kz = self._dyadic(np.maximum(cz - i0z, i1z - cz))
        out = np.empty(len(lo))
        key = (kx.astype(np.int64) << np.int64(40)) + (ky.astype(np.int64) << np.int64(20)) + kz
        for u in np.unique(key):
            m = key == u
            a, b, c = int(kx[m][0]), int(ky[m][0]), int(kz[m][0])
            out[m] = self._filt(a, b, c)[cz[m], cy[m], cx[m]]
        return out


def all_edge_keys(T, chunk=CH):
    n = len(T)
    out = np.empty(n * 6, np.int64)
    for s0 in range(0, n, chunk):
        Tc = T[s0:s0 + chunk]
        a = np.minimum(Tc[:, EI], Tc[:, EJ]).astype(np.int64)
        b = np.maximum(Tc[:, EI], Tc[:, EJ]).astype(np.int64)
        out[s0 * 6:(s0 + len(Tc)) * 6] = (a * NVMAX + b).ravel()
    return out


def longest_edges(P, T, chunk=CH):
    """(key, kslot) of each tet's longest edge under a strict total order."""
    n = len(T)
    LE = np.empty(n, np.int64)
    KS = np.empty(n, np.int8)
    for s0 in range(0, n, chunk):
        Tc = T[s0:s0 + chunk]
        E = np.linalg.norm(P[Tc[:, EI]] - P[Tc[:, EJ]], axis=2)
        a6 = np.minimum(Tc[:, EI], Tc[:, EJ]).astype(np.int64)
        b6 = np.maximum(Tc[:, EI], Tc[:, EJ]).astype(np.int64)
        k6 = a6 * NVMAX + b6
        tie = E == E.max(1, keepdims=True)
        k = np.where(tie, k6, np.int64(-1)).argmax(1)
        r = np.arange(len(Tc))
        LE[s0:s0 + len(Tc)] = k6[r, k]
        KS[s0:s0 + len(Tc)] = k.astype(np.int8)
    return LE, KS


def failing(P, T, vs, gate, chunk=CH, pooled="bary"):
    """Cells with Vs/dx < gate, and the worst Vs/dx.

    pooled=False -- THE GATE.  Vs nearest-grid at the BARYCENTER (locked).
    pooled -- THE REFINEMENT TARGET, one of:
      "bary"  the gate itself; cheapest, but can treadmill if children land in
              slower bins faster than dx shrinks
      "half"  Vs pooled over a box of HALF the cell's extent about the
              barycentre.  After one bisection a child's barycentre lies within
              about half the parent's extent, so this bounds the NEXT
              generation and is re-evaluated every round -- most of the
              protection of a full bound at a fraction of the refinement
      "bbox"  Vs pooled over the whole cell.  A strict lower bound on any
              descendant, so provably convergent, but wildly conservative on
              km-scale cells (it takes the slowest Vs anywhere inside)
    """
    out, fmin = [], np.inf
    for s0 in range(0, len(T), chunk):
        Tc = T[s0:s0 + chunk]
        Pc = P[Tc]
        dx = np.linalg.norm(Pc[:, EI] - Pc[:, EJ], axis=2).max(1)
        if pooled == "bbox":
            v = vs.bbox_min(Pc.min(1), Pc.max(1))
        elif pooled == "half":
            b = Pc.mean(1)
            r = 0.25 * (Pc.max(1) - Pc.min(1))
            v = vs.bbox_min(b - r, b + r)
        else:
            v = vs.at(Pc.mean(1))
        fr = np.where(dx > 0, v / dx, 0.0)
        fmin = min(fmin, float(fr.min()))
        b = np.nonzero(fr < gate)[0]
        if b.size:
            out.append(b + s0)
    return (np.concatenate(out) if out else np.zeros(0, np.int64)), fmin


def ragged(lo, hi):
    cnt = (hi - lo).astype(np.int64)
    if int(cnt.sum()) == 0:
        return np.zeros(0, np.int64)
    off = np.concatenate([[0], np.cumsum(cnt)[:-1]])
    return (np.arange(int(cnt.sum()), dtype=np.int64)
            - np.repeat(off, cnt) + np.repeat(lo.astype(np.int64), cnt))


def split_tets(T, sel, kslot, mid):
    """Bisect tets `sel` at their longest edge; vectorised, orientation-preserving."""
    Ts = T[sel]
    ks = kslot[sel]
    p, q = EI[ks], EJ[ks]
    out = []
    for keep, gone in ((p, q), (q, p)):
        child = Ts.copy()
        child[np.arange(len(Ts)), gone] = mid
        out.append(child)
    return np.concatenate(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--collar", required=True)
    ap.add_argument("--parent", required=True)
    ap.add_argument("--cvm", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--gate", type=float, default=0.6667)
    ap.add_argument("--target", choices=["bary", "half", "bbox"], default="half",
                    help="refinement target; see failing()")
    ap.add_argument("--hops", type=int, default=5)
    ap.add_argument("--all", action="store_true",
                    help="refine the whole collar instead of a k-hop patch; the "
                         "collar is standalone so there is no outer mesh to protect, "
                         "and this removes rim freezing entirely")
    ap.add_argument("--max-rounds", type=int, default=40)
    ap.add_argument("--stats", default=None)
    a = ap.parse_args()
    t0 = time.time()

    vs = PoolVs(a.cvm)
    S = np.load(a.collar)
    P = S["points"].copy()
    T = S["tets"].astype(np.int32)
    wall_collar = S["wall_collar_idx"]
    nt0, nv0 = len(T), len(P)
    print(f"[collar] {nt0:,} tets, {nv0:,} verts   RSS {rss_gb():.1f} GB", flush=True)

    # ---- frozen wall edges -------------------------------------------------
    import h5py
    with h5py.File(a.parent, "r") as f:
        G = f["geometry"][:]
        C = f["connect"][:].astype(np.int64)
        B = f["boundary"][:].astype(np.int32)
    PW, TW, _, _, _ = extract_wall(G, C, B)
    del G, C, B
    if len(PW) != len(wall_collar):
        raise RuntimeError(f"wall vertex count {len(PW):,} != collar's {len(wall_collar):,}")
    if not np.allclose(P[wall_collar], PW, atol=1e-6):
        raise RuntimeError("collar wall vertices do not match the parent's wall")
    wtri = wall_collar[TW].astype(np.int64)
    wk = []
    for i, j in ((0, 1), (1, 2), (0, 2)):
        x = np.minimum(wtri[:, i], wtri[:, j])
        y = np.maximum(wtri[:, i], wtri[:, j])
        wk.append(x * NVMAX + y)
    wall_keys = np.unique(np.concatenate(wk))
    del wk
    print(f"[wall]  {len(TW):,} tris -> {len(wall_keys):,} frozen edges", flush=True)

    bad0, f0 = failing(P, T, vs, a.gate)
    tgt0, t0m = failing(P, T, vs, a.gate, pooled=a.target)
    print(f"[gate]   {len(bad0):,} failing cells, worst {f0:.4f}")
    print(f"[target] ({a.target}) {len(tgt0):,} cells below target, worst {t0m:.4f}   "
          f"RSS {rss_gb():.1f} GB", flush=True)
    if not len(tgt0):
        print("nothing to do")
        return

    # ---- working patch ------------------------------------------------------
    if a.all:
        patch = np.arange(nt0)
        interior = np.ones(nv0, bool)
        touch = None
    else:
        sel_v = np.zeros(nv0, bool)
        sel_v[np.unique(T[tgt0])] = True
        interior = sel_v.copy()
        for h in range(a.hops):
            touch = sel_v[T].any(1)
            if h == a.hops - 1:
                interior = sel_v.copy()
            sel_v[np.unique(T[touch])] = True
        patch = np.nonzero(sel_v[T].any(1))[0]
    print(f"[patch] {len(patch):,} tets ({100*len(patch)/nt0:.2f}%), "
          f"{int(interior.sum()):,} bisectable verts", flush=True)

    keep = np.ones(nt0, bool)
    keep[patch] = False
    Tp = T[patch].copy()
    T_keep = T[keep]
    del T, keep
    bis_ok = interior.copy()
    hist = []

    for rd in range(a.max_rounds):
        marked, fmin_p = failing(P, Tp, vs, a.gate, pooled=a.target)
        if not len(marked):
            print(f"[leb] round {rd}: patch CLEAN on the {a.target} target "
                  f"(worst {fmin_p:.4f})", flush=True)
            break
        LE, kslot = longest_edges(P, Tp)
        allk = all_edge_keys(Tp)
        owner = np.repeat(np.arange(len(Tp), dtype=np.int32), 6)
        o = np.argsort(allk, kind="stable")
        allk, owner = allk[o], owner[o]
        del o
        uk, ustart = np.unique(allk, return_index=True)
        uend = np.concatenate([ustart[1:], [len(allk)]])
        n_shell = uend - ustart
        pos = np.searchsorted(uk, LE)
        n_le = np.bincount(pos, minlength=len(uk))
        terminal = n_le == n_shell

        # LEPP: walk from the marked cells to terminal edges.
        #
        # Membership is a BOOLEAN MASK, not np.isin against the growing front.
        # isin re-sorts the whole front on every hop: with 2.6M marked cells and
        # up to 200 hops that is ~10^10 comparisons and the round never
        # finishes (measured: no round 0 output after 2.5 min of 100 % CPU).
        in_front = np.zeros(len(Tp), bool)
        in_front[marked] = True
        front = marked.copy()
        for _ in range(200):
            kpos = np.searchsorted(uk, LE[front])
            nonterm = ~terminal[kpos]
            if not nonterm.any():
                break
            fi = ragged(ustart[kpos[nonterm]], uend[kpos[nonterm]])
            cand = np.unique(owner[fi])
            new = cand[~in_front[cand]]
            if not len(new):
                break
            in_front[new] = True
            front = new                      # only the NEW frontier needs walking
        front = np.nonzero(in_front)[0]
        kpos = np.searchsorted(uk, LE[front])
        term_keys = np.unique(uk[kpos[terminal[kpos]]])
        if not len(term_keys):
            print(f"[leb] round {rd}: NO terminal edge reachable -- stop", flush=True)
            break

        x = (term_keys // NVMAX).astype(np.int64)
        y = (term_keys % NVMAX).astype(np.int64)
        wp = np.clip(np.searchsorted(wall_keys, term_keys), 0, len(wall_keys) - 1)
        is_wall = wall_keys[wp] == term_keys
        ok = bis_ok[x] & bis_ok[y] & ~is_wall
        n_rim = int((~(bis_ok[x] & bis_ok[y])).sum())
        n_wall = int(is_wall.sum())
        term_keys = term_keys[ok]
        if not len(term_keys):
            print(f"[leb] round {rd}: all {n_rim+n_wall} terminal edges frozen "
                  f"(rim {n_rim}, wall {n_wall}) -- stop", flush=True)
            break

        x = (term_keys // NVMAX).astype(np.int64)
        y = (term_keys % NVMAX).astype(np.int64)
        newP = 0.5 * (P[x] + P[y])
        mid_id = np.arange(len(P), len(P) + len(newP), dtype=np.int64)
        P = np.vstack([P, newP])
        bis_ok = np.concatenate([bis_ok, np.ones(len(newP), bool)])
        lpos = np.clip(np.searchsorted(term_keys, LE), 0, len(term_keys) - 1)
        hit = term_keys[lpos] == LE
        sel = np.nonzero(hit)[0]
        newT = split_tets(Tp, sel, kslot, mid_id[lpos[sel]])
        km = np.ones(len(Tp), bool)
        km[sel] = False
        Tp = np.vstack([Tp[km], newT.astype(np.int32)])
        hist.append(dict(round=rd, marked=int(len(marked)), worst=float(fmin_p),
                         terminal=int(len(term_keys)), frozen_rim=n_rim,
                         frozen_wall=n_wall, split=int(len(sel)),
                         patch_tets=int(len(Tp)), rss_gb=round(rss_gb(), 2)))
        print(f"[leb] round {rd}: {len(marked):,} marked (worst {fmin_p:.4f}) | "
              f"{len(term_keys):,} terminal (froze rim {n_rim} wall {n_wall}) | "
              f"split {len(sel):,} -> patch {len(Tp):,} | RSS {rss_gb():.1f} GB "
              f"| {time.time()-t0:.0f}s", flush=True)
        del allk, owner, uk, ustart, uend, n_shell, n_le, terminal, LE, kslot

    tets = np.vstack([T_keep, Tp]).astype(np.int64)
    del T_keep, Tp
    print(f"\n[out] {len(tets):,} tets (+{len(tets)-nt0:,}, "
          f"{100*(len(tets)-nt0)/nt0:+.2f}%), {len(P):,} verts (+{len(P)-nv0:,})",
          flush=True)

    nb, fnew = failing(P, tets, vs, a.gate)
    print(f"[gate] {len(nb):,} failing cells, worst {fnew:.4f} "
          f"(was {len(bad0):,}, worst {f0:.4f})")

    # the wall must be untouched
    if not np.allclose(P[wall_collar], PW, atol=1e-12):
        raise RuntimeError("LEB moved wall vertices")
    sv_neg = 0
    for s0 in range(0, len(tets), CH):
        sv_neg += int((tet_signed_volume(P, tets[s0:s0 + CH]) < 0).sum())
    emin, etamin, n01 = np.inf, np.inf, 0
    for s0 in range(0, len(tets), CH):
        Tc = tets[s0:s0 + CH]
        e = tet_edge_lengths(P, Tc)
        et = tet_eta(P, Tc)
        emin = min(emin, float(e.min()))
        etamin = min(etamin, float(et.min()))
        n01 += int((et < 0.1).sum())
    print(f"[quality] eta min {etamin:.4f}, eta<0.1 {n01:,}, min edge {emin:,.1f} m, "
          f"inverted {sv_neg}")
    if sv_neg:
        raise RuntimeError("LEB produced inverted tets")

    out = {k: S[k] for k in S.files if k not in ("points", "tets")}
    np.savez_compressed(a.out, points=P, tets=tets, **out)
    print(f"[write] {a.out}  ({time.time()-t0:.0f} s)")
    if a.stats:
        json.dump(dict(rounds=hist, gate_before=len(bad0), gate_after=len(nb),
                       worst_before=f0, worst_after=fnew,
                       tets_before=nt0, tets_after=len(tets)),
                  open(a.stats, "w"), indent=1)


if __name__ == "__main__":
    main()
