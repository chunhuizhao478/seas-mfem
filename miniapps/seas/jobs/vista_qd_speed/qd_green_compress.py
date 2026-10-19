#!/usr/bin/env python3
"""How far the dense Green's function G of a --qd-green cache compresses.

QD speed item 4 (compressed G) is worth building only if G, compressed so
that the stage product stays at the round-off level, reads clearly fewer
bytes than dense G: on Vista the dense product already streams G at the
node's memory bandwidth.  This script measures that on a real cache file
(SEAS-QD-GREEN-2, domain/qd_green_operator.hpp), with no change to the
solver:

    python3 jobs/vista_qd_speed/qd_green_compress.py G.bin [--leaf 64]
        [--eta 2] [--eps 1e-8 1e-10 1e-12 1e-13] [--procs N]

Method (route 3 of the Mac 2 km study, item 3 notes):
  - cluster tree on the fault DOFs' node positions (the file's keys),
    bisecting the longest box side until a cluster has <= leaf DOFs;
  - a block (target cluster s, source cluster t) is admissible when
    max(diam s, diam t) <= eta * dist(s, t); the others down to the leaves
    stay dense (near field);
  - each admissible block B (3|s| rows x 2|t| columns) is cut to the
    smallest rank k with sqrt(sum_{i>=k} s_i^2) <= eps * ||B||_F, and kept
    dense when k (m + n) >= m n;
  - stored bytes per eps against dense, and the product error
    max|(G - G_H) x| / max|G x| for five slips x: the solver's validation
    slip, uniform strike slip and three random slips.  For scale: G built
    on 4 and on 8 ranks differ by about 1e-15 in that measure (2 km study),
    the round-off floor of G itself.

Memory: G in RAM once (34.9 GiB at p = 1 on the production mesh), shared
copy-on-write by the worker processes (fork).  Blocks are independent, so
the SVDs run in --procs processes with OPENBLAS_NUM_THREADS threads each.
"""

import argparse
import multiprocessing as mp
import os
import struct
import sys
import time

import numpy as np

G = None          # set before the fork; read-only in the workers
X = None
NCOMP = 3
EPS = []


def load(path):
    with open(path, "rb") as f:
        h = f.read(512)
    if h[:15] != b"SEAS-QD-GREEN-2":
        sys.exit(f"ERROR: {path}: not a SEAS-QD-GREEN-2 file ({h[:16]!r})")
    nf, ncomp, ncols = struct.unpack("<3q", h[16:40])
    info = h[80:].split(b"\0")[0].decode()
    off = 512
    keys = np.fromfile(path, np.float64, 6 * nf, offset=off).reshape(nf, 6)
    off += 48 * nf
    g = np.fromfile(path, np.float64, ncomp * nf * ncols, offset=off)
    g = g.reshape(ncomp * nf, ncols)
    off += 8 * g.size
    size = os.path.getsize(path)
    if size != off + 8 * ncomp * nf:
        sys.exit(f"ERROR: {path}: {size} bytes, expected {off + 8 * ncomp * nf}")
    return nf, ncomp, ncols, info, keys, g


def build_tree(P, leaf):
    clusters = []

    def build(ids):
        lo, hi = P[ids].min(0), P[ids].max(0)
        c = dict(ids=ids, kids=[], lo=lo, hi=hi)
        clusters.append(c)
        if len(ids) > leaf:
            ax = int(np.argmax(hi - lo))
            o = ids[np.argsort(P[ids, ax], kind="stable")]
            h = len(o) // 2
            c["kids"] = [build(o[:h]), build(o[h:])]
        return c

    sys.setrecursionlimit(100000)
    return build(np.arange(len(P))), clusters


def diam(c):
    return float(np.linalg.norm(c["hi"] - c["lo"]))


def dist(a, b):
    return float(np.linalg.norm(np.maximum(0, np.maximum(a["lo"] - b["hi"],
                                                         b["lo"] - a["hi"]))))


def partition(root, eta):
    adm, near = [], []

    def walk(s, t):
        d = dist(s, t)
        if d > 0 and max(diam(s), diam(t)) <= eta * d:
            adm.append((s["ids"], t["ids"]))
            return
        if not s["kids"] or not t["kids"]:
            near.append((s["ids"], t["ids"]))
            return
        for a in s["kids"]:
            for b in t["kids"]:
                walk(a, b)

    walk(root, root)
    return adm, near


def rows_of(ids):
    return (NCOMP * ids[:, None] + np.arange(NCOMP)).ravel()


def cols_of(ids):
    return (2 * ids[:, None] + np.arange(2)).ravel()


def compress_block(job):
    """One admissible block: per eps the rank, the bytes, the squared
    Frobenius tail dropped, and the rows' share of (G - G_H) x."""
    sids, tids = job
    ri, ci = rows_of(sids), cols_of(tids)
    B = G[np.ix_(ri, ci)]
    m, n = B.shape
    U, sv, Vt = np.linalg.svd(B, full_matrices=False)
    fro = float(np.sqrt(np.sum(sv ** 2)))
    tail = np.sqrt(np.cumsum((sv ** 2)[::-1])[::-1])
    tail = np.append(tail, 0.0)
    Xc = X[ci]
    out = []
    for e in EPS:
        k = int(np.argmax(tail <= e * fro))
        if k * (m + n) >= m * n:
            out.append((k, 8.0 * m * n, 0.0, None))
        else:
            ex = U[:, k:] @ (sv[k:, None] * (Vt[k:] @ Xc))
            out.append((k, 8.0 * k * (m + n), float(tail[k] ** 2), ex))
    return ri, m, n, out


def main():
    global G, X, NCOMP, EPS
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("file")
    ap.add_argument("--leaf", type=int, default=64)
    ap.add_argument("--eta", type=float, default=2.0)
    ap.add_argument("--eps", type=float, nargs="+",
                    default=[1e-8, 1e-10, 1e-12, 1e-13])
    ap.add_argument("--procs", type=int, default=max(1, (os.cpu_count() or 1) // 4))
    a = ap.parse_args()
    EPS = sorted(a.eps, reverse=True)

    T0 = time.time()
    nf, NCOMP, ncols, info, keys, G = load(a.file)
    G.setflags(write=False)
    print(f"{a.file}: {info}")
    print(f"N_f {nf}, {NCOMP} components, {ncols} columns, dense G "
          f"{G.nbytes / 2**30:.3f} GiB; loaded in {time.time() - T0:.0f} s")

    j = np.arange(ncols, dtype=float)
    rng = np.random.default_rng(12345)
    X = np.column_stack([np.sin(0.37 * j + 0.1) + 0.5 * np.cos(0.011 * j),
                         np.tile([0.0, 1.0], nf)]
                        + [rng.standard_normal(ncols) for _ in range(3)])
    names = ["validation slip", "uniform strike", "random #1", "random #2",
             "random #3"]
    t0 = time.time()
    Y = G @ X
    ymax = np.abs(Y).max(axis=0)
    gnorm = float(np.linalg.norm(G))
    print(f"dense product of the 5 test slips: {time.time() - t0:.1f} s")

    root, clusters = build_tree(keys[:, :3], a.leaf)
    adm, near = partition(root, a.eta)
    near_bytes = sum(8.0 * NCOMP * len(s) * 2 * len(t) for s, t in near)
    print(f"tree: leaf <= {a.leaf} DOFs, {len(clusters)} clusters; eta {a.eta}: "
          f"{len(near)} dense near-field blocks ({near_bytes / 2**30:.3f} GiB), "
          f"{len(adm)} admissible blocks")

    # Biggest blocks first, so the pool does not end on one long SVD.
    adm.sort(key=lambda st: -len(st[0]) * len(st[1]))
    EX = {e: np.zeros_like(Y) for e in EPS}
    byts = {e: 0.0 for e in EPS}
    tail2 = {e: 0.0 for e in EPS}
    ranks = {e: [] for e in EPS}
    t0 = time.time()
    with mp.get_context("fork").Pool(a.procs) as pool:
        for done, (ri, m, n, out) in enumerate(
                pool.imap_unordered(compress_block, adm, chunksize=1), 1):
            for e, (k, b, t2, ex) in zip(EPS, out):
                ranks[e].append(k)
                byts[e] += b
                tail2[e] += t2
                if ex is not None:
                    EX[e][ri] += ex
            if done % max(1, len(adm) // 10) == 0:
                print(f"  {done}/{len(adm)} blocks, {time.time() - t0:.0f} s",
                      flush=True)
    print(f"SVD of the admissible blocks: {time.time() - t0:.0f} s on "
          f"{a.procs} processes")

    dense_bytes = float(G.nbytes)
    print("\neps      stored GiB  of dense  smaller  near share  "
          "||G-G_H||_F/||G||_F  rank median/p90/max  product error per slip")
    for e in EPS:
        tot = byts[e] + near_bytes
        r = np.array(ranks[e])
        ep = np.abs(EX[e]).max(axis=0) / ymax
        print(f"{e:<8g} {tot / 2**30:10.3f}  {tot / dense_bytes:8.3f}  "
              f"x{dense_bytes / tot:6.2f}  {near_bytes / tot:9.1%}  "
              f"{np.sqrt(tail2[e]) / gnorm:19.3e}  "
              f"{int(np.median(r))}/{int(np.percentile(r, 90))}/{r.max():<10d}  "
              + ", ".join(f"{nm} {v:.2e}" for nm, v in zip(names, ep)))
    print(f"\ntotal wall {time.time() - T0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
