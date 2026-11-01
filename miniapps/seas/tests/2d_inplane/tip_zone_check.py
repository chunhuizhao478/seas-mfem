#!/usr/bin/env python3
# tests/2d_inplane/tip_zone_check.py — checks the tip-zone field of
# seas_test_tip_zone's dumps (--dump-mesh, --dump-deck) on a gmsh 2.2 mesh
# against chains built here, independently of cycle/fault_chains.hpp
# (Chunhui Zhao, 2026-10-29: "let's write this in the document, and prepare to
# add this function"; the UCERF 2-D faults, report rev 15, section 3.1).
#
# Here: the fault segments are the mesh's line elements of physical tag 3;
# chains are walked from the end with the smaller (x, y); an end on a line
# element of another tag (the outer boundary) is not a free tip; a point goes
# to its nearest segment (all segments, numpy), s = arc length from the
# chain's first end, d = the distance to the nearer free tip along the chain,
# and the value is written out again:
#   W = min(width, max_fraction L);  tip for d <= W - taper, inside for
#   d >= W, linear in between.
#
#   python3 tip_zone_check.py MESH --dump PREFIX [--deck] [--width 20e3]
#       [--taper 2e3] [--inside 0.004] [--tip 0.04] [--max-fraction F]
#       [--same PREFIX2] [--near-tip X Y --radius R]
# --dump PREFIX reads PREFIX.<rank> (all ranks).  Without --deck the rows are
# x y chain s d a (--dump-mesh: chain, s and d are checked too); with --deck
# they are x y a.  --same PREFIX2: the sorted rows of both dumps must be
# identical text (another rank count).  --near-tip X Y: lists the points of
# OTHER chains within R of the free tip nearest (X, Y), with their own d and
# a (the distance is along the own chain, never a straight line).
# Exit 0 when every check passes.

import argparse
import glob
import math
import sys

import numpy as np


def read_msh(path):
    nodes, segs, other = {}, [], set()
    with open(path) as f:
        it = iter(f)
        for line in it:
            if line.startswith('$MeshFormat'):
                v = next(it).split()[0]
                if not v.startswith('2'):
                    sys.exit('tip_zone_check: gmsh 2.2 ASCII only, got ' + v)
            elif line.startswith('$Nodes'):
                for _ in range(int(next(it))):
                    p = next(it).split()
                    nodes[int(p[0])] = (float(p[1]), float(p[2]))
            elif line.startswith('$Elements'):
                for _ in range(int(next(it))):
                    p = list(map(int, next(it).split()))
                    et, nt = p[1], p[2]
                    vs = p[3 + nt:]
                    if et != 1:
                        continue
                    if p[3] == 3:
                        segs.append((vs[0], vs[1]))
                    else:
                        other.update(vs)
    return nodes, segs, other


def build_chains(nodes, segs, other):
    adj = {}
    for k, (a, b) in enumerate(segs):
        adj.setdefault(a, []).append(k)
        adj.setdefault(b, []).append(k)
    bad = [v for v in adj if len(adj[v]) >= 3]
    if bad:
        sys.exit('tip_zone_check: junction at %s' % (nodes[bad[0]],))
    ends = sorted((v for v in adj if len(adj[v]) == 1), key=lambda v: nodes[v])
    used, chains = set(), []
    for v0 in ends:
        if adj[v0][0] in used:
            continue
        pts, v, k = [v0], v0, adj[v0][0]
        while True:
            used.add(k)
            a, b = segs[k]
            w = b if a == v else a
            pts.append(w)
            v = w
            if len(adj[v]) == 1:
                break
            k = adj[v][1] if adj[v][0] == k else adj[v][0]
        xy = np.array([nodes[p] for p in pts])
        s = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(xy, axis=0).T))])
        chains.append(dict(xy=xy, s=s, L=s[-1], free0=pts[0] not in other,
                           free1=pts[-1] not in other))
    if len(used) != len(segs):
        sys.exit('tip_zone_check: %d segments in closed chains' % (len(segs) - len(used)))
    return chains


def value(d, L, a):
    W = min(a.width, a.max_fraction * L) if a.max_fraction > 0 else a.width
    if d <= W - a.taper:
        return a.tip
    if d >= W:
        return a.inside
    return a.tip + (a.inside - a.tip) * (d - (W - a.taper)) / a.taper


def read_dump(prefix):
    rows = []
    for f in sorted(glob.glob(prefix + '.*')):
        if not f[len(prefix) + 1:].isdigit():   # PREFIX.<rank> only
            continue
        with open(f) as h:
            rows += [line.split() for line in h if line.strip()]
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mesh')
    ap.add_argument('--dump', required=True)
    ap.add_argument('--deck', action='store_true')
    ap.add_argument('--width', type=float, default=20e3)
    ap.add_argument('--taper', type=float, default=2e3)
    ap.add_argument('--inside', type=float, default=0.004)
    ap.add_argument('--tip', type=float, default=0.04)
    ap.add_argument('--max-fraction', type=float, default=0.0)
    ap.add_argument('--same')
    ap.add_argument('--near-tip', type=float, nargs=2)
    ap.add_argument('--radius', type=float, default=10e3)
    a = ap.parse_args()

    nodes, segs, other = read_msh(a.mesh)
    chains = build_chains(nodes, segs, other)
    print('%s: %d chains, %d segments, %d free tips' % (a.mesh, len(chains), len(segs),
          sum(c['free0'] + c['free1'] for c in chains)))
    for i, c in enumerate(chains):
        print('  chain %d: L = %.3f km, (%.2f, %.2f) -> (%.2f, %.2f) km' % (
            i, c['L'] / 1e3, *(c['xy'][0] / 1e3), *(c['xy'][-1] / 1e3)))
    # every segment, in walk order
    A = np.concatenate([c['xy'][:-1] for c in chains])
    B = np.concatenate([c['xy'][1:] for c in chains])
    S0 = np.concatenate([c['s'][:-1] for c in chains])
    CH = np.concatenate([np.full(len(c['xy']) - 1, i) for i, c in enumerate(chains)])
    E = B - A
    LEN = np.hypot(E[:, 0], E[:, 1])

    rows = read_dump(a.dump)
    if not rows:
        sys.exit('tip_zone_check: no rows in %s.*' % a.dump)
    P = np.array([[float(r[0]), float(r[1])] for r in rows])
    ok = True
    worst = dict(off=0.0, a=0.0, s=0.0, d=0.0)
    nchain = 0
    regime = {}
    res = []
    for k in range(0, len(P), 2000):
        p = P[k:k + 2000]
        t = ((p[:, None, 0] - A[None, :, 0]) * E[None, :, 0]
             + (p[:, None, 1] - A[None, :, 1]) * E[None, :, 1]) / (LEN[None, :] ** 2)
        t = np.clip(t, 0.0, 1.0)
        qx = p[:, None, 0] - (A[None, :, 0] + t * E[None, :, 0])
        qy = p[:, None, 1] - (A[None, :, 1] + t * E[None, :, 1])
        dist = np.hypot(qx, qy)
        j = np.argmin(dist, axis=1)
        for m, jj in enumerate(j):
            res.append((jj, t[m, jj], dist[m, jj]))
    for r, (jj, t, dist) in zip(rows, res):
        c = chains[CH[jj]]
        s = S0[jj] + t * LEN[jj]
        d = min(s if c['free0'] else math.inf, c['L'] - s if c['free1'] else math.inf)
        want = value(d, c['L'], a)
        got = float(r[-1])
        worst['off'] = max(worst['off'], dist / LEN[jj])
        worst['a'] = max(worst['a'], abs(got - want))
        if not a.deck:
            nchain += int(int(r[2]) != CH[jj])
            worst['s'] = max(worst['s'], abs(float(r[3]) - s))
            worst['d'] = max(worst['d'], abs(float(r[4]) - d))
        key = (CH[jj], 'tip' if want == a.tip else 'inside' if want == a.inside else 'ramp')
        regime[key] = regime.get(key, 0) + 1
    print('%d rows from %s.*: max distance to the segment / its length %.2e; max |a - a_here| %.3e'
          % (len(rows), a.dump, worst['off'], worst['a']))
    ok &= worst['off'] <= 1e-6 and worst['a'] <= 1e-12
    if not a.deck:
        print('  chain index differs at %d rows; max |s - s_here| %.3e m, max |d - d_here| %.3e m'
              % (nchain, worst['s'], worst['d']))
        ok &= nchain == 0 and worst['s'] <= 1e-6 and worst['d'] <= 1e-6
    for i in range(len(chains)):
        print('  chain %d: tip / ramp / inside points %d / %d / %d' % (
            i, regime.get((i, 'tip'), 0), regime.get((i, 'ramp'), 0), regime.get((i, 'inside'), 0)))

    if a.near_tip:
        q = np.array(a.near_tip)
        tips = [(i, e) for i, c in enumerate(chains) for e, free in
                ((c['xy'][0], c['free0']), (c['xy'][-1], c['free1'])) if free]
        ci, tip = min(tips, key=lambda it: np.hypot(*(it[1] - q)))
        print('near-tip: the free tip of chain %d at (%.3f, %.3f) km; other chains\' points '
              'within %.1f km of it:' % (ci, tip[0] / 1e3, tip[1] / 1e3, a.radius / 1e3))
        n_near = 0
        dmin = math.inf
        for r, (jj, t, dist) in zip(rows, res):
            if CH[jj] == ci:
                continue
            e = math.hypot(float(r[0]) - tip[0], float(r[1]) - tip[1])
            if e > a.radius:
                continue
            c = chains[CH[jj]]
            s = S0[jj] + t * LEN[jj]
            d = min(s, c['L'] - s)
            n_near += 1
            dmin = min(dmin, e)
            if float(r[-1]) != value(d, c['L'], a):
                ok = False
            if d > a.width and float(r[-1]) != a.inside:
                ok = False
        print('  %d points (nearest %.2f km from the tip); every one has its own chain\'s '
              'value%s' % (n_near, dmin / 1e3, ' (a = inside: own d > width)' if n_near else ''))

    if a.same:
        other_rows = read_dump(a.same)
        same = sorted(' '.join(r) for r in rows) == sorted(' '.join(r) for r in other_rows)
        print('same rows as %s.* (sorted text): %s' % (a.same, 'yes' if same else 'NO'))
        ok &= same
    print('RESULT: %s' % ('PASS' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
