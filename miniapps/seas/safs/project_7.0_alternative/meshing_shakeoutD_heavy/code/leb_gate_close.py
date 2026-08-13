#!/usr/bin/env python3
"""leb_gate_close.py -- close the resolved-frequency gate on a frozen-parent
collar EXACTLY, by conforming 3-D longest-edge (Rivara) bisection.

Why bisection and not another mmg size-map pass
-----------------------------------------------
mmg accepts any edge within ~[0.71, 1.41]x of its metric, so a regrade LICENSES
coarsening: compliant cells drift back over tolerance and a pass mints failures
as fast as it fixes them.  Bisection can only SPLIT, so compliance is monotone
and the loop provably converges.  It inserts edge MIDPOINTS only, so a flat top
stays exactly flat and every vertex of the frozen parent keeps its coordinates.

The driving metric is NOT the acceptance metric
-----------------------------------------------
Acceptance (locked, `gate_census.py`): f = Vs / dx, dx = element MAX edge, Vs
nearest-grid at the element BARYCENTRE.  Targeting that directly TREADMILLS:
the CVM is a nearest-grid lookup on a coarse vertical lattice, so bisecting a
cell that straddles a bin edge throws one child into the SLOWER bin, where it
needs a smaller dx than its parent did.  The loop then adds cells forever.

This tool drives on `dx <= Vs_pool / gate`, where Vs_pool is the MINIMUM Vs over
the CVM z levels the cell's OWN VERTICAL EXTENT spans (sampled at the
barycentre's column).  Any descendant's barycentre lies inside the cell, so its
measured Vs is >= Vs_pool: the target is a lower bound, it can only RISE under
refinement, and a compliant cell stays compliant.  Pooling is in z ONLY -- a
full 3-D min over a multi-km cell is correct but wildly over-conservative
(measured +55.8 % tets elsewhere), and the CVM's lateral step (1500 m) is far
coarser than its vertical one, so bin changes are essentially always vertical.

What is frozen
--------------
* every tet of the parent region (indices < n_parent_tets) -- never touched;
* every edge whose BOTH endpoints are parent vertices (ids < n_parent_verts).
  That is exactly the welded seam wall (plus the top/bottom rims), where a split
  would leave a hanging node against an unsplit parent tet.  Freezing a few
  collar-interior diagonals as well is conservative and harmless.
Cells whose longest-edge chain terminates on a frozen edge are reported as
`seam-blocked`: a structural class, not a defect.

Usage:
    leb_gate_close.py --mesh merged.puml.h5 --parent parent.puml.h5 \
        --cvm cvm.nc --gate 0.6667 --out closed.puml.h5 [--max-rounds 60]
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from puml_io import (BC_DYNAMIC_RUPTURE, LOCAL_FACES, face_code,
                     read_puml, write_puml)  # noqa: F401
from material import Material

PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
NVMAX = np.int64(1) << np.int64(27)          # vertex-id packing base (134 M)
CH = 3_000_000
import os
DEBUG = bool(os.environ.get('LEB_DEBUG'))


# --------------------------------------------------------------- CVM ------
class Vs:
    """Nearest-grid Vs on the CVM lattice, and the z-pooled minimum."""

    def __init__(self, path):
        import netCDF4 as ncdf
        d = ncdf.Dataset(str(path))
        self.X = d["x"][:].data; self.Y = d["y"][:].data; self.Z = d["z"][:].data
        D = d["data"][:]
        self.V = np.sqrt(np.maximum(D["mu"], 0) / D["rho"]).astype(np.float32)
        del D
        dz = np.diff(self.Z)
        if not np.allclose(dz, dz[0]):
            raise RuntimeError("CVM z axis is not uniform; gate_census assumes it is")
        self.nz = len(self.Z)

    def _ij(self, p):
        i = np.clip(np.rint((p[:, 0] - self.X[0]) / (self.X[1] - self.X[0])).astype(np.int32), 0, len(self.X) - 1)
        j = np.clip(np.rint((p[:, 1] - self.Y[0]) / (self.Y[1] - self.Y[0])).astype(np.int32), 0, len(self.Y) - 1)
        return i, j

    def _k(self, z):
        return np.clip(np.rint((z - self.Z[0]) / (self.Z[1] - self.Z[0])).astype(np.int32), 0, self.nz - 1)

    def at(self, p):
        i, j = self._ij(p)
        return self.V[self._k(p[:, 2]), j, i]

    def pooled(self, v):
        """min Vs over the z levels the cell (v: (n,4,3)) spans, at its column."""
        b = v.mean(1)
        i, j = self._ij(b)
        kb = self._k(b[:, 2])
        k0 = np.minimum(self._k(v[:, :, 2].min(1)), kb)
        k1 = np.maximum(self._k(v[:, :, 2].max(1)), kb)
        out = self.V[k0, j, i].copy()
        span = int((k1 - k0).max()) if len(k0) else 0
        for o in range(1, span + 1):
            # NOTE: `np.minimum(out[m], x, out=out[m])` silently does nothing --
            # `out[m]` is a fancy-index COPY on both sides, so the result lands
            # in a temporary and the pool is never applied.  Index explicitly.
            sel = np.flatnonzero((k0 + o) <= k1)
            if not len(sel):
                continue
            out[sel] = np.minimum(out[sel], self.V[k0[sel] + o, j[sel], i[sel]])
        return out


# ------------------------------------------------------------ geometry ----
def _edge_keys(T):
    ii = np.array([p[0] for p in PAIRS]); jj = np.array([p[1] for p in PAIRS])
    a = np.minimum(T[:, ii], T[:, jj]).astype(np.int64)
    b = np.maximum(T[:, ii], T[:, jj]).astype(np.int64)
    return a * NVMAX + b                                   # (n,6)


def _edge_lengths(P, T):
    return np.stack([np.linalg.norm(P[T[:, a]] - P[T[:, b]], axis=1) for a, b in PAIRS], 1)


def longest_edge(P, T):
    """(key, slot) of each tet's longest edge under a STRICT TOTAL ORDER.

    `L.argmax(1)` breaks ties by LOCAL SLOT.  Meshes carrying red/LEB refinement
    are full of exactly equal edge lengths, so adjacent tets can each name the
    other's edge and the LEPP chain closes into a CYCLE with no terminal edge.
    Ordering ties by the globally unique edge key makes the relation a strict
    total order, so the chain always terminates.
    """
    n = len(T)
    key = np.empty(n, np.int64); slot = np.empty(n, np.int8)
    for s in range(0, n, CH):
        t = T[s:s + CH]
        L = _edge_lengths(P, t); K = _edge_keys(t)
        Kt = np.where(L >= L.max(1, keepdims=True) - 1e-12, K, np.int64(-1))
        k = Kt.max(1)
        key[s:s + CH] = k
        slot[s:s + CH] = (Kt == k[:, None]).argmax(1)
        del L, K, Kt
    return key, slot


def occurrences(T, le_key, cand, nv):
    """All (key, owner) of the candidate edges, plus whether the owner agrees.

    Only tets touching a candidate ENDPOINT can own a candidate edge, and the
    candidates are a localised handful, so the whole-collar key scan is skipped
    for everything else -- this is what keeps a round at seconds rather than
    minutes on a 35 M-cell collar.
    """
    cand = np.sort(cand)
    vm = np.zeros(nv, bool)
    vm[(cand // NVMAX).astype(np.int64)] = True
    vm[(cand % NVMAX).astype(np.int64)] = True
    sub = []
    for s in range(0, len(T), CH):
        sub.append(np.flatnonzero(vm[T[s:s + CH]].any(1)) + s)
    sub = np.concatenate(sub)
    hk, ho = [], []
    for s in range(0, len(sub), CH):
        idx = sub[s:s + CH]
        K = _edge_keys(T[idx])
        pos = np.searchsorted(cand, K)
        np.clip(pos, 0, len(cand) - 1, out=pos)
        hit = cand[pos] == K
        if not hit.any():
            continue
        r, c = np.nonzero(hit)
        hk.append(K[r, c]); ho.append(idx[r])
        del K, pos, hit
    if not hk:
        return (np.zeros(0, np.int64),) * 2 + (np.zeros(0, bool),)
    hk = np.concatenate(hk); ho = np.concatenate(ho)
    o = np.argsort(hk, kind="stable")
    hk, ho = hk[o], ho[o]
    return hk, ho, le_key[ho] == hk


def _write(a, h5py, Pm, T, W, R, NPAR, ntot):
    """Write the product, streaming the frozen parent block straight through."""
    BLK = 8_000_000
    tmp = str(a.out) + ".part"
    with h5py.File(a.mesh) as fin, h5py.File(tmp, "w") as fo:
        fo.create_dataset("geometry", data=np.ascontiguousarray(Pm, np.float64))
        dc = fo.create_dataset("connect", (ntot, 4), np.uint64)
        db = fo.create_dataset("boundary", (ntot,), np.int32)
        dg = fo.create_dataset("group", (ntot,), np.int32)
        for s0 in range(0, NPAR, BLK):                 # frozen parent, verbatim
            s1 = min(s0 + BLK, NPAR)
            dc[s0:s1] = fin["connect"][s0:s1]
            db[s0:s1] = fin["boundary"][s0:s1]
            dg[s0:s1] = (fin["group"][s0:s1] if "group" in fin
                         else np.ones(s1 - s0, np.int32))
        dc[NPAR:] = T.astype(np.uint64)
        db[NPAR:] = W
        dg[NPAR:] = R
    Path(tmp).replace(a.out)          # atomic: a kill never leaves a torn file


# --------------------------------------------------------------- main -----
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--parent", required=True)
    ap.add_argument("--cvm", required=True, help="deck safs_material_cvm.nc")
    ap.add_argument("--muscal", default=None,
                    help="native MUSCAL.nc; the gate is then min(deck, MUSCAL)")
    ap.add_argument("--box", type=float, nargs=4,
                    default=[132000.0, 724000.0, 3490000.0, 4055000.0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--gate", type=float, required=True)
    ap.add_argument("--max-rounds", type=int, default=200)
    ap.add_argument("--allow-fault-split", action="store_true",
                    help="also refine the fault triangulation itself "
                         "(implies --include-parent); changes the DR facet count")
    ap.add_argument("--include-parent", action="store_true",
                    help="also refine the frozen parent block; then the ONLY frozen "
                         "edges are the fault's own, so the fault triangulation "
                         "(and every deck asset built on it) is untouched")
    ap.add_argument("--checkpoint-every", type=int, default=20,
                    help="write the product every N rounds (0 = only at the end)")
    ap.add_argument("--vs-source", choices=("min", "muscal", "deck"), default="min",
                    help="which cube the gate is SCORED on")
    ap.add_argument("--drive", choices=("pooled", "accept"), default="pooled",
                    help="stopping bound for already-split cells")
    ap.add_argument("--safety", type=float, default=0.85,
                    help="--drive accept: required margin on Vs_bary/gate")
    ap.add_argument("--pool-frac", type=float, default=0.25,
                    help="barycentre wander radius as a fraction of dx")
    ap.add_argument("--lepp-depth", type=int, default=60,
                    help="hard cap on LEPP growth iterations per round")
    ap.add_argument("--stats")
    a = ap.parse_args()
    t0 = time.time()

    import h5py
    globals()['h5py'] = h5py
    with h5py.File(a.parent) as f:
        NPAR = f["connect"].shape[0]; NPV = f["geometry"].shape[0]
    # The parent block is frozen, so it is never held in RAM: only the collar is
    # loaded, and the parent rows are copied straight through at write time.
    # At 133.7 M parent tets that block alone is 4.3 GB.
    if a.include_parent or a.allow_fault_split:
        NPAR = 0
    with h5py.File(a.mesh) as f:
        NT = f["connect"].shape[0]
        G = f["geometry"][:]
        T = f["connect"][NPAR:].astype(np.int64)
        W = f["boundary"][NPAR:].astype(np.int32)
        R = (f["group"][NPAR:].astype(np.int32) if "group" in f
             else np.ones(NT - NPAR, np.int32))
    # --- which edges may never be split -----------------------------------
    # default: the welded seam (both ends parent vertices) -- splitting one
    #          would leave a hanging node against an unsplit parent tet.
    # --include-parent: the fault's own edges -- splitting one would move the
    #          dynamic-rupture surface that every deck's stress, friction and
    #          receiver asset is built on.
    FKEY = None
    if a.allow_fault_split:
        # Nothing is frozen.  Splitting a FAULT edge is geometrically safe --
        # the midpoint of an edge of a fault triangle lies ON that triangle, so
        # the dynamic-rupture surface is SUBDIVIDED, not moved, and the two
        # halves inherit BC 3 on both sides (every tet in the edge's shell is
        # split, so the face stays conforming and still interior x2).  What it
        # DOES change is the fault facet count -- and it can shrink inradii in
        # the fault band, which is where dt lives.  Both are measured after.
        FKEY = np.zeros(0, np.int64)
        print("[freeze] NOTHING frozen: fault edges may be split", flush=True)
    elif a.include_parent:
        ftri = []
        for s4 in range(4):
            m = face_code(W, s4) == BC_DYNAMIC_RUPTURE
            if m.any():
                ftri.append(T[m][:, LOCAL_FACES[s4]])
        if ftri:
            ft = np.vstack(ftri)
            e = np.concatenate([ft[:, [0, 1]], ft[:, [1, 2]], ft[:, [0, 2]]])
            lo = np.minimum(e[:, 0], e[:, 1]).astype(np.int64)
            hi = np.maximum(e[:, 0], e[:, 1]).astype(np.int64)
            FKEY = np.unique(lo * NVMAX + hi)
            print(f"[freeze] {len(ft):,} fault facets -> {len(FKEY):,} frozen fault edges",
                  flush=True)
        else:
            FKEY = np.zeros(0, np.int64)

    def frozen_edges(keys):
        if FKEY is None:
            return ((keys // NVMAX) < NPV) & ((keys % NVMAX) < NPV)
        if not len(FKEY):
            return np.zeros(len(keys), bool)
        pos = np.clip(np.searchsorted(FKEY, keys), 0, len(FKEY) - 1)
        return FKEY[pos] == keys
    if NT <= NPAR:
        raise RuntimeError("parent block mismatch")
    print(f"[mesh] {NT:,} tets / {len(G):,} verts   "
          f"parent {NPAR:,} tets / {NPV:,} verts   collar {len(T):,} tets", flush=True)
    vs = Material(a.cvm, a.muscal, box=tuple(a.box), source=a.vs_source)

    P = [G]; nv = len(G)
    NT0 = NT
    le_key = le_slot = None
    blocked = np.zeros(len(T), bool)
    touched = np.zeros(len(T), bool)
    active = None
    hist = []

    for rnd in range(a.max_rounds):
        Pm = P[0] if len(P) == 1 else np.concatenate(P)
        # ---- who to refine -------------------------------------------------
        # ACCEPTANCE (what the locked census scores): dx <= Vs(barycentre)/gate.
        # POOLED (a monotone lower bound): dx <= min Vs over the cell's own
        # vertical extent / gate.  Refining is driven by acceptance failures,
        # but a cell that has ALREADY been split keeps being refined until it is
        # POOLED-compliant -- that is what stops the treadmill, because a
        # pooled-compliant cell cannot produce a child that newly fails, while
        # an acceptance-compliant one can (its child's barycentre may cross a
        # CVM bin edge into slower material).  Cells never touched are left
        # alone: splitting elsewhere cannot change them.
        # Only cells that FAILED at round 0, or that a split has since created,
        # can ever be marked: an untouched, passing cell is not changed by a
        # split elsewhere.  Scanning just those turns a whole-collar sweep per
        # round (35 M cells) into a sweep of the working front.
        fail = np.zeros(len(T), bool); over = np.zeros(len(T), bool)
        worst = np.inf
        idxs = np.arange(len(T)) if active is None else np.flatnonzero(active)
        for s in range(0, len(idxs), CH):
            ii = idxs[s:s + CH]
            t = T[ii]; v = Pm[t]
            dx = _edge_lengths(Pm, t).max(1)
            vb = vs.at(v.mean(1))
            f_ = np.where(vb > 0, vb / dx, 0.0)
            fail[ii] = f_ < a.gate
            if a.drive == "pooled":
                # min Vs over the window a descendant's barycentre can reach
                over[ii] = dx > (vs.pooled(v, dx, a.pool_frac) / a.gate)
            else:
                # A cell only DRAGGED into a chain must not inherit the pooled
                # demand.  Near the free surface the pooled bound is the 116 m/s
                # surface value (~174 m), far stricter than the gate needs
                # there, so every chain reaching the surface layer spawned a
                # cascade and the failure count oscillated instead of falling.
                # Require instead a margin on the ACCEPTANCE metric itself:
                # a child that flips into material up to (1-safety) slower
                # still passes, and dx shrinks geometrically while the demand
                # is bounded below, so the loop converges.
                over[ii] = dx > (a.safety * vb / a.gate)
            worst = min(worst, float(f_.min()) if len(f_) else np.inf)
            del t, v, dx, vb, f_
        if active is None:
            active = fail.copy()
        marked = over & (fail | touched) & ~blocked
        nmark = int(marked.sum())
        hist.append(dict(round=rnd, tets=int(len(T)), fail=int(fail.sum()),
                         over=int(over.sum()), marked=nmark,
                         blocked=int(blocked.sum()), worst_raw_f=round(worst, 5)))
        print(f"[r{rnd:02d}] collar {len(T):,}  gate-fail {int(fail.sum()):,}  "
              f"pooled-over {int(over.sum()):,}  marked {nmark:,}  "
              f"blocked {int(blocked.sum()):,}  worst f {worst:.4f}", flush=True)
        if nmark == 0:
            break

        # ---- LEPP -----------------------------------------------------------
        # le_key is CACHED: bisection only appends vertices, so an unsplit tet's
        # geometry -- and therefore its longest edge -- is untouched.  Only the
        # children need recomputing.  (Recomputing all of them each round was
        # ~40 s of a ~57 s round on the 150 M-cell mesh.)
        if le_key is None:
            le_key, le_slot = longest_edge(Pm, T)
        mk = marked.copy()
        term_keys = []
        prev_nterm = 0
        stagnant = 0
        for it in range(a.lepp_depth):
            if not mk.any():
                break
            cand = np.unique(le_key[mk])
            hk, ho, agree = occurrences(T, le_key, cand, len(Pm))
            edges, first = np.unique(hk, return_index=True)
            cnt = np.diff(np.append(first, len(hk)))
            nagr = np.add.reduceat(agree.astype(np.int64), first)
            isterm = nagr == cnt
            froz = frozen_edges(edges)
            seg = np.repeat(np.arange(len(edges)), cnt)
            nterm = int(isterm.sum())
            if DEBUG:
                print(f"    it{it:03d} mk {int(mk.sum()):,} cand {len(edges):,} "
                      f"term {nterm:,} frozen-term {int((isterm & froz).sum()):,}", flush=True)
            if nterm:
                mk[ho[isterm[seg]]] = False          # resolved: split or blocked
                fb = isterm & froz
                if fb.any():
                    blocked[ho[fb[seg]]] = True
                term_keys.append(edges[isterm & ~froz])
            # Stop when the terminal set stops GROWING, not after a fixed
            # number of links.  A shallow cap sees only the chains that happen
            # to terminate at once, which in a graded collar are exactly the
            # frozen seam edges: measured 82 terminal (all frozen) at depth 2,
            # against 45,527 once the chains were allowed to drain.
            # LEPP growth is UNEVEN -- it plateaus for an iteration or two and
            # then jumps (measured 82, 137, 11933, 21810 ...).  Breaking on the
            # first non-growing iteration truncates the drain and the round
            # splits a handful of edges instead of thousands, which shows up as
            # the terminal count collapsing 139-83-48-26-13-5 while the failure
            # count barely moves.  Require several stagnant iterations.
            if nterm and prev_nterm and nterm < 1.01 * prev_nterm:
                stagnant += 1
                if stagnant >= 3:
                    break
            else:
                stagnant = 0
            prev_nterm = max(prev_nterm, nterm)
            if nterm == len(edges):
                break
            o = ~isterm[seg]
            sh = ho[o]
            dis = sh[le_key[sh] != hk[o]]
            new_ = np.unique(dis[~mk[dis]])
            if not len(new_):
                if nterm:
                    break
                raise RuntimeError(f"round {rnd}: LEPP stalled with no terminal edge")
            mk[new_] = True
        term_keys = np.unique(np.concatenate(term_keys)) if term_keys else np.zeros(0, np.int64)
        if len(term_keys) == 0:
            print("[leb] every remaining chain ends on a frozen edge", flush=True)
            break
        hk, ho, _ = occurrences(T, le_key, term_keys, len(Pm))
        edges, first = np.unique(hk, return_index=True)
        cnt = np.diff(np.append(first, len(hk)))

        own = ho
        eid = np.repeat(np.arange(len(edges)), cnt)
        mid = nv + eid
        P.append(0.5 * (Pm[edges // NVMAX] + Pm[edges % NVMAX])); nv += len(edges)

        childT = np.empty((2 * len(own), 4), np.int64)
        childW = np.empty(2 * len(own), np.int32)
        slots = le_slot[own]
        for k, (pp, qq) in enumerate(PAIRS):
            m = slots == k
            if not m.any():
                continue
            idx = own[m]; mm = mid[m]
            base = T[idx]
            ca = base.copy(); ca[:, qq] = mm            # child that KEEPS pp
            cb = base.copy(); cb[:, pp] = mm            # child that KEEPS qq
            # a face keeps its BC code iff it still contains the vertex kept
            ka = np.int32(0); kb = np.int32(0)
            for s4 in range(4):
                fv = LOCAL_FACES[s4]
                if pp in fv: ka |= np.int32(0xFF) << np.int32(8 * s4)
                if qq in fv: kb |= np.int32(0xFF) << np.int32(8 * s4)
            w = W[idx]
            n0 = np.flatnonzero(m)
            childT[2 * n0] = ca; childT[2 * n0 + 1] = cb
            childW[2 * n0] = w & ka; childW[2 * n0 + 1] = w & kb
        childR = np.repeat(R[own], 2)

        alive = np.ones(len(T), bool); alive[own] = False
        T = np.vstack([T[alive], childT])
        W = np.concatenate([W[alive], childW])
        R = np.concatenate([R[alive], childR])
        blocked = np.concatenate([blocked[alive], np.zeros(2 * len(own), bool)])
        touched = np.concatenate([touched[alive], np.ones(2 * len(own), bool)])
        active = np.concatenate([active[alive], np.ones(2 * len(own), bool)])
        ck, cs = longest_edge(Pm2 := np.concatenate(P) if len(P) > 1 else P[0], childT)
        le_key = np.concatenate([le_key[alive], ck])
        le_slot = np.concatenate([le_slot[alive], cs])
        del Pm2
        print(f"       split {len(edges):,} terminal edges, {len(own):,} tets -> "
              f"{2*len(own):,}   collar now {len(T):,}", flush=True)
        # Checkpoint: this machine is shared and the OOM killer has taken these
        # runs mid-flight.  Writing only at the end throws away every round.
        if a.checkpoint_every and (rnd + 1) % a.checkpoint_every == 0:
            Pw = P[0] if len(P) == 1 else np.concatenate(P)
            _write(a, h5py, Pw, T, W, R, NPAR, NPAR + len(T))
            print(f"       [checkpoint] wrote {NPAR + len(T):,} tets", flush=True)
            del Pw

    Pm = P[0] if len(P) == 1 else np.concatenate(P)
    ntot = NPAR + len(T)
    print(f"\n[out] {ntot:,} tets (+{ntot-NT0:,}), "
          f"{len(Pm):,} verts (+{len(Pm)-len(G):,})", flush=True)
    assert np.array_equal(Pm[:NPV], G[:NPV]), "parent vertices moved"
    _write(a, h5py, Pm, T, W, R, NPAR, ntot)
    if a.stats:

        Path(a.stats).write_text(json.dumps(dict(
            rounds=hist, tets=int(ntot), verts=int(len(Pm)),
            parent_tets=int(NPAR), seconds=round(time.time() - t0, 1)), indent=1))
    print(f"[write] {a.out}  ({time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
