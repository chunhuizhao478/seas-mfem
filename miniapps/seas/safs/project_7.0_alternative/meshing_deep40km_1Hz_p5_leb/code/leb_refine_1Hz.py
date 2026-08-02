#!/usr/bin/env python3
"""leb_refine_1Hz.py -- close the 1 Hz @ p5 gate by conforming 3-D longest-edge
(Rivara) bisection, vectorized for a 10^8-cell mesh.

WHY BISECTION AND NOT mmg.  Two mmg size-map regrades of this mesh REGRESSED the
gate (worst Vs/dx 0.6667 -> 0.4104 / 0.4137) because mmg's +-41 % acceptance band
licenses COARSENING: compliant cells drift over tolerance and the pass mints
failures as fast as it fixes them.  Bisection cannot do that --

  * it only ever SPLITS, so no cell's dx can grow;
  * an UNSPLIT cell's barycenter does not move, so its measured Vs is unchanged
    and it cannot newly fail;
  * a split child is strictly smaller, so it can only fail if its barycenter
    crosses into a slower CVM bin -- and as cells shrink their barycenters
    converge inside a single bin, so the outer iteration CONVERGES;
  * midpoints only, so geometry is reproduced exactly (the ALT top is exactly
    flat, so split top edges stay at z = 0 and the free surface is bit-preserved);
  * deterministic -- no mmg draw lottery.

FROZEN: fault edges are never bisected (the fault triangulation must stay
bit-identical for deck compatibility).  Non-fault edges that merely touch a
fault VERTEX are fine -- bisecting them does not move the fault surface.

Usage:
    python leb_refine_1Hz.py --mesh in.puml.h5 --cvm cvm.nc --out out.puml.h5 \
        --gate 0.8 [--hops 5] [--max-rounds 40]
"""
import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
import h5py

sys.path.insert(0, str(Path(__file__).resolve().parent))
from puml_io import (BC_DYNAMIC_RUPTURE, LOCAL_FACES, face_code,  # noqa: E402
                     tet_edge_lengths, tet_eta, tet_signed_volume, write_puml)

PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
EI = np.array([p[0] for p in PAIRS])
EJ = np.array([p[1] for p in PAIRS])
INFACE = np.zeros((4, 4), bool)
for _s, _fv in enumerate(LOCAL_FACES):
    for _v in _fv:
        INFACE[_s, _v] = True
NVMAX = np.int64(1) << np.int64(26)     # > 20.9M verts, and grows during refinement
CH = 4_000_000


def rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1 << 30)


class Vs:
    def __init__(self, path):
        import netCDF4 as ncdf
        d = ncdf.Dataset(str(path))
        self.X, self.Y, self.Z = d["x"][:].data, d["y"][:].data, d["z"][:].data
        raw = d["data"][:]
        self.mu = np.asarray(raw["mu"], np.float64)
        self.rho = np.asarray(raw["rho"], np.float64)
        d.close()

    def at(self, p):
        ix = np.clip(np.rint((p[:, 0] - self.X[0]) / (self.X[1] - self.X[0])).astype(int),
                     0, len(self.X) - 1)
        iy = np.clip(np.rint((p[:, 1] - self.Y[0]) / (self.Y[1] - self.Y[0])).astype(int),
                     0, len(self.Y) - 1)
        iz = np.clip(np.rint((p[:, 2] - self.Z[0]) / (self.Z[1] - self.Z[0])).astype(int),
                     0, len(self.Z) - 1)
        return np.sqrt(np.maximum(self.mu[iz, iy, ix], 0.0) / self.rho[iz, iy, ix])


# ---------------------------------------------------------------------------
# EVERY geometry pass below is CHUNKED.  P[T[:, EI]] on the full mesh is a
# 122e6 x 6 x 3 float64 temporary = 17.6 GB and swap-thrashes the machine; the
# chunked form peaks at ~300 MB per buffer.
# ---------------------------------------------------------------------------

def all_edge_keys(T, chunk=CH):
    """Flat (6n,) sorted-pair edge keys, tet-major."""
    n = len(T)
    out = np.empty(n * 6, np.int64)
    for s0 in range(0, n, chunk):
        Tc = T[s0:s0 + chunk]
        a = np.minimum(Tc[:, EI], Tc[:, EJ]).astype(np.int64)
        b = np.maximum(Tc[:, EI], Tc[:, EJ]).astype(np.int64)
        out[s0 * 6:(s0 + len(Tc)) * 6] = (a * NVMAX + b).ravel()
    return out


def longest_edges(P, T, chunk=CH):
    """(key, kslot) of each tet's longest edge, under a STRICT TOTAL ORDER on edges:
    first by length, ties broken by the (globally unique) edge key.

    *** THIS TIE-BREAK IS LOAD-BEARING, NOT COSMETIC. ***
    Rivara LEPP terminates only because edge length strictly increases along the
    chain.  This mesh was built by EXACT RED (1->4) refinement, whose children are
    geometrically SIMILAR to their parent -- so it contains huge numbers of EXACTLY
    equal edge lengths.  Breaking ties by local slot index (the obvious
    `E.argmax(1)`) is not globally consistent: two adjacent tets can each name a
    different edge of an equal-length pair, the LEPP chain closes into a CYCLE, and
    no terminal edge is ever found.  Measured with the slot tie-break: only 32,557
    terminal edges were reachable from 703,220 marked cells and the count plateaued
    instead of converging.  Ordering ties by edge key makes "longest edge" a strict
    total order, so chains cannot cycle and LEPP always terminates.
    """
    n = len(T)
    LE = np.empty(n, np.int64)
    KS = np.empty(n, np.int8)
    for s0 in range(0, n, chunk):
        Tc = T[s0:s0 + chunk]
        E = np.linalg.norm(P[Tc[:, EI]] - P[Tc[:, EJ]], axis=2)
        a6 = np.minimum(Tc[:, EI], Tc[:, EJ]).astype(np.int64)
        b6 = np.maximum(Tc[:, EI], Tc[:, EJ]).astype(np.int64)
        k6 = a6 * NVMAX + b6                                # (m,6) edge keys
        tie = E == E.max(1, keepdims=True)                  # exact-length ties
        k = np.where(tie, k6, np.int64(-1)).argmax(1)       # largest key among them
        r = np.arange(len(Tc))
        LE[s0:s0 + len(Tc)] = k6[r, k]
        KS[s0:s0 + len(Tc)] = k.astype(np.int8)
    return LE, KS


def failing(P, T, vs, gate, chunk=CH, pooled=False):
    """Indices of cells with Vs/dx < gate, plus the global worst Vs/dx.

    pooled=False  -- THE GATE.  Vs nearest-grid at the BARYCENTER.  This is the
                     locked acceptance rule and is what the census reports.
    pooled=True   -- THE REFINEMENT TARGET.  Vs = MIN over the cell's own VERTICAL
                     extent, sampled at the barycenter's (x, y): the 4 vertex z's
                     and the barycenter z.

    Why Z-ONLY and not a full 3-D pool over the cell: the CVM lattice is 1500 m
    laterally but only 250 m vertically, while these cells are <= ~500 m across, so
    a barycenter that changes bins on refinement almost always does so VERTICALLY.
    A full 3-D min over the cell takes the slowest Vs anywhere inside it, which on
    coarse multi-km cells is wildly over-conservative -- measured +55.8% tets on the
    small ALT mesh versus +0.04% for the raw barycenter rule.  Pooling in z alone
    costs a fraction of that and still removes the treadmill.

    Why the target must be pooled: the CVM is nearest-grid on a 250 m vertical
    lattice and Vs steps ~4x at the z = -125 m bin edge (567 -> 2323 m/s).  Under
    the barycenter rule, bisecting a cell that straddles that edge sends one child
    into the SLOW bin, where it needs a 4x smaller dx than its parent did -- so the
    pass mints failures as fast as it fixes them.  Measured: the count stalled at
    ~65k, falling 1.8%/round while adding 157k tets/round.
    Every child barycenter lies inside the parent tet, so min-over-samples is a
    LOWER bound on any child's measured Vs.  Sizing to it is therefore
    refinement-stable: once a cell complies it STAYS compliant, and because the
    pool is over the cell's OWN extent (not a fixed window) deep fast cells are
    not penalised -- a cell wholly inside one bin pools that bin's value.
    """
    out, fmin = [], np.inf
    for s0 in range(0, len(T), chunk):
        Tc = T[s0:s0 + chunk]
        Pc = P[Tc]
        E = np.linalg.norm(Pc[:, EI] - Pc[:, EJ], axis=2)
        dx = E.max(1)
        bary = Pc.mean(1)
        v = vs.at(bary)
        if pooled:
            q = bary.copy()
            for k in range(4):
                q[:, 2] = Pc[:, k, 2]                  # same (x,y), vertex z
                np.minimum(v, vs.at(q), out=v)
        fr = np.where(dx > 0, v / dx, 0.0)
        fmin = min(fmin, float(fr.min()))
        b = np.nonzero(fr < gate)[0]
        if b.size:
            out.append(b + s0)
    idx = np.concatenate(out) if out else np.zeros(0, np.int64)
    return idx, fmin


def quality(P, T, chunk=CH):
    """Chunked (eta_min, n_eta_lt01, edge_min, n_inverted)."""
    eta_min, emin, n01, ninv = np.inf, np.inf, 0, 0
    for s0 in range(0, len(T), chunk):
        Tc = T[s0:s0 + chunk]
        e = tet_edge_lengths(P, Tc)
        emin = min(emin, float(e.min()))
        et = tet_eta(P, Tc)
        eta_min = min(eta_min, float(et.min()))
        n01 += int((et < 0.1).sum())
        ninv += int((tet_signed_volume(P, Tc) < 0).sum())
    return eta_min, n01, emin, ninv


def ragged(lo, hi):
    """Flat indices of the concatenated ranges [lo_i, hi_i), plus the group id."""
    cnt = (hi - lo).astype(np.int64)
    tot = int(cnt.sum())
    if tot == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    grp = np.repeat(np.arange(len(lo), dtype=np.int64), cnt)
    off = np.concatenate([[0], np.cumsum(cnt)[:-1]])
    idx = np.arange(tot, dtype=np.int64) - np.repeat(off, cnt) + np.repeat(lo.astype(np.int64), cnt)
    return idx, grp


def split_tets(P, T, W, G, sel, kslot, mid_of_tet):
    """Bisect tets `sel` at their longest edge; vectorized, orientation-preserving."""
    Ts, Ws, Gs = T[sel], W[sel].astype(np.int64), G[sel]
    ks = kslot[sel]
    p = EI[ks]                      # local slot of endpoint a
    q = EJ[ks]                      # local slot of endpoint b
    m = mid_of_tet                  # (len(sel),) new midpoint vertex id
    out_T, out_W, out_G = [], [], []
    for keep, gone in ((p, q), (q, p)):
        child = Ts.copy()
        child[np.arange(len(Ts)), gone] = m
        cw = np.zeros(len(Ts), np.int64)
        for s in range(4):
            code = (Ws >> (8 * s)) & 0xFF
            has_p, has_q = INFACE[s, p], INFACE[s, q]
            keep_in, gone_in = INFACE[s, keep], INFACE[s, gone]
            take = (has_p & has_q) | (keep_in & ~gone_in)
            cw |= np.where(take, code, 0) << (8 * s)
        out_T.append(child)
        out_W.append(cw.astype(np.int32))
        out_G.append(Gs)
    return (np.concatenate(out_T), np.concatenate(out_W), np.concatenate(out_G))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--cvm", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--gate", type=float, default=0.8)
    ap.add_argument("--hops", type=int, default=5)
    ap.add_argument("--max-rounds", type=int, default=40)
    ap.add_argument("--stats", default=None)
    args = ap.parse_args()
    t0 = time.time()

    vs = Vs(args.cvm)
    with h5py.File(args.mesh, "r") as f:
        P = f["geometry"][:]
        ds = f["connect"]
        nt0 = ds.shape[0]
        C = np.empty((nt0, 4), np.int32)          # chunked read: never hold the uint64 copy
        for s0 in range(0, nt0, 8_000_000):
            C[s0:s0 + 8_000_000] = ds[s0:s0 + 8_000_000].astype(np.int32)
        B = f["boundary"][:].astype(np.int32)
        GR = f["group"][:].astype(np.int32)
    nv0 = len(P)
    print(f"[mesh] {nt0:,} tets, {nv0:,} verts   RSS {rss_gb():.1f} GB", flush=True)

    # ---- frozen fault edges -------------------------------------------------
    fk = []
    fverts = []
    for s in range(4):
        sel = np.nonzero(face_code(B, s) == BC_DYNAMIC_RUPTURE)[0]
        if not sel.size:
            continue
        tri = C[sel][:, LOCAL_FACES[s]].astype(np.int64)
        fverts.append(tri.ravel())
        for i, j in ((0, 1), (1, 2), (0, 2)):
            a = np.minimum(tri[:, i], tri[:, j])
            b = np.maximum(tri[:, i], tri[:, j])
            fk.append(a * NVMAX + b)
    fault_keys = np.unique(np.concatenate(fk))
    fault_v0 = np.unique(np.concatenate(fverts))
    P_fault0 = P[fault_v0].copy()
    del fk, fverts
    print(f"[fault] {len(fault_keys):,} frozen fault edges  RSS {rss_gb():.1f} GB", flush=True)

    # ---- locate failures and build the working patch ------------------------
    bad0, f0min = failing(P, C, vs, args.gate)                      # THE GATE
    tgt0, t0min = failing(P, C, vs, args.gate, pooled=True)          # refinement target
    print(f"[gate]   {len(bad0):,} failing cells, worst {f0min:.4f}", flush=True)
    print(f"[target] {len(tgt0):,} cells below the POOLED target, worst {t0min:.4f}   "
          f"RSS {rss_gb():.1f} GB", flush=True)
    if not len(tgt0):
        print("nothing to do")
        return

    sel_v = np.zeros(nv0, bool)
    sel_v[np.unique(C[tgt0])] = True
    interior = sel_v.copy()
    for h in range(args.hops):
        touch = sel_v[C].any(1)
        if h == args.hops - 1:
            interior = sel_v.copy()          # verts strictly inside -> bisectable
        sel_v[np.unique(C[touch])] = True
        print(f"   hop {h}: {int(touch.sum()):,} tets touched, "
              f"{int(sel_v.sum()):,} verts   RSS {rss_gb():.1f} GB", flush=True)
    patch = np.nonzero(sel_v[C].any(1))[0]
    print(f"[patch] {len(patch):,} tets ({100*len(patch)/nt0:.2f}%), "
          f"{int(interior.sum()):,} bisectable verts   RSS {rss_gb():.1f} GB", flush=True)
    del touch

    keepmask = np.ones(nt0, bool)
    keepmask[patch] = False
    T = C[patch].copy()
    W = B[patch].copy()
    G = GR[patch].copy()
    C_keep, B_keep, G_keep = C[keepmask], B[keepmask], GR[keepmask]
    del C, B, GR, keepmask
    bis_ok = interior.copy()                 # per-vertex; new midpoints inherit True
    hist = []

    for rd in range(args.max_rounds):
        marked, fmin_p = failing(P, T, vs, args.gate, pooled=True)
        if not len(marked):
            print(f"[leb] round {rd}: patch CLEAN on the pooled target "
                  f"(worst {fmin_p:.4f})", flush=True)
            break
        LE, kslot = longest_edges(P, T)
        allk = all_edge_keys(T)
        owner = np.repeat(np.arange(len(T), dtype=np.int32), 6)   # int32: halves peak RSS
        o = np.argsort(allk, kind="stable")
        allk, owner = allk[o], owner[o]
        del o
        uk, ustart = np.unique(allk, return_index=True)
        uend = np.concatenate([ustart[1:], [len(allk)]])
        n_shell = uend - ustart
        # how many tets have each unique edge as their LONGEST edge
        pos = np.searchsorted(uk, LE)
        n_le = np.bincount(pos, minlength=len(uk))
        terminal = n_le == n_shell           # per unique edge

        # ---- LEPP: walk to terminal edges --------------------------------
        front = marked.copy()
        for _ in range(200):
            kpos = np.searchsorted(uk, LE[front])
            nonterm = ~terminal[kpos]
            if not nonterm.any():
                break
            lo, hi = ustart[kpos[nonterm]], uend[kpos[nonterm]]
            fi, _g = ragged(lo, hi)
            cand = np.unique(owner[fi])
            new = cand[~np.isin(cand, front, assume_unique=False)]
            if not len(new):
                break
            front = np.concatenate([front, new])
        kpos = np.searchsorted(uk, LE[front])
        term_keys = np.unique(uk[kpos[terminal[kpos]]])
        if not len(term_keys):
            print(f"[leb] round {rd}: NO terminal edge reachable -- stop", flush=True)
            break

        # ---- freeze filters ------------------------------------------------
        a = (term_keys // NVMAX).astype(np.int64)
        b = (term_keys % NVMAX).astype(np.int64)
        fp = np.clip(np.searchsorted(fault_keys, term_keys), 0, len(fault_keys) - 1)
        is_fault = fault_keys[fp] == term_keys
        ok = bis_ok[a] & bis_ok[b] & ~is_fault
        n_rim = int((~(bis_ok[a] & bis_ok[b])).sum())
        n_flt = int(is_fault.sum())
        term_keys = term_keys[ok]
        if not len(term_keys):
            print(f"[leb] round {rd}: all {n_rim+n_flt} terminal edges frozen "
                  f"(rim {n_rim}, fault {n_flt}) -- stop", flush=True)
            break

        # ---- create midpoints and split -------------------------------------
        a, b = (term_keys // NVMAX).astype(np.int64), (term_keys % NVMAX).astype(np.int64)
        newP = 0.5 * (P[a] + P[b])
        mid_id = np.arange(len(P), len(P) + len(newP), dtype=np.int64)
        P = np.vstack([P, newP])
        bis_ok = np.concatenate([bis_ok, np.ones(len(newP), bool)])
        lpos = np.searchsorted(term_keys, LE)
        lpos = np.clip(lpos, 0, len(term_keys) - 1)
        hit = term_keys[lpos] == LE
        sel = np.nonzero(hit)[0]
        nt_new, nw_new, ng_new = split_tets(P, T, W, G, sel, kslot, mid_id[lpos[sel]])
        km = np.ones(len(T), bool)
        km[sel] = False
        T = np.vstack([T[km], nt_new.astype(np.int32)])
        W = np.concatenate([W[km], nw_new])
        G = np.concatenate([G[km], ng_new])
        hist.append(dict(round=rd, marked=int(len(marked)), worst=float(fmin_p),
                         terminal=int(len(term_keys)), frozen_rim=n_rim, frozen_fault=n_flt,
                         split=int(len(sel)), patch_tets=int(len(T)),
                         rss_gb=round(rss_gb(), 2)))
        print(f"[leb] round {rd}: {len(marked):,} marked (worst {fmin_p:.4f}) | "
              f"{len(term_keys):,} terminal (froze rim {n_rim} fault {n_flt}) | "
              f"split {len(sel):,} -> patch {len(T):,} | RSS {rss_gb():.1f} GB "
              f"| {time.time()-t0:.0f}s", flush=True)
        del allk, owner, uk, ustart, uend, n_shell, n_le, terminal, LE, kslot

    # ---- reassemble ---------------------------------------------------------
    connect = np.vstack([C_keep, T])
    boundary = np.concatenate([B_keep, W])
    group = np.concatenate([G_keep, G])
    del C_keep, B_keep, G_keep, T, W, G
    print(f"\n[out] {len(connect):,} tets (+{len(connect)-nt0:,}, "
          f"{100*(len(connect)-nt0)/nt0:+.2f}%), {len(P):,} verts (+{len(P)-nv0:,})   "
          f"RSS {rss_gb():.1f} GB", flush=True)

    nb_idx, fmin_new = failing(P, connect, vs, args.gate)
    nb = int(len(nb_idx))
    eta_min, n_eta01, edge_min, ninv = quality(P, connect)
    assert np.array_equal(P[fault_v0], P_fault0), "FAULT VERTICES MOVED"
    topv = np.unique(np.concatenate(
        [connect[face_code(boundary, s) == 1][:, LOCAL_FACES[s]].ravel() for s in range(4)]))
    print(f"  gate: worst Vs/dx {fmin_new:.4f}  ->  p3 {0.75*fmin_new:.4f} Hz, "
          f"p5 {1.25*fmin_new:.4f} Hz;  below {args.gate}: {nb:,} (was {len(bad0):,})")
    print(f"  inverted {ninv};  fault vertices unmoved: PASS")
    print(f"  free surface z [{P[topv,2].min():.6f}, {P[topv,2].max():.6f}]")
    print(f"  eta min {eta_min:.4f}, eta<0.1 {n_eta01:,}, min edge {edge_min:.4f} m")
    assert ninv == 0, f"{ninv} inverted tets"
    if args.stats:
        json.dump(dict(rounds=hist, n_tets=int(len(connect)), n_verts=int(len(P)),
                       n_fail_before=int(len(bad0)), n_fail_after=nb,
                       worst_after=float(fmin_new),
                       resolved_p3=0.75 * float(fmin_new),
                       resolved_p5=1.25 * float(fmin_new),
                       eta_min=float(eta_min), eta_lt01=int(n_eta01),
                       edge_min=float(edge_min), inverted=ninv,
                       peak_rss_gb=round(rss_gb(), 2)),
                  open(args.stats, "w"), indent=2)
    write_puml(args.out, P, connect, boundary, group)
    print(f"[write] {args.out}   peak RSS {rss_gb():.1f} GB   ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
