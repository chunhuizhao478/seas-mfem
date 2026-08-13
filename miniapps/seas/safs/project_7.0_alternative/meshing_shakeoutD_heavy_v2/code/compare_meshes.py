#!/usr/bin/env python3
"""compare_meshes.py -- identical measurement of two PUML meshes, side by side.

Written because quoting each mesh's numbers from its own build log compares
different definitions: eta conventions, whether the gate was scored on the deck
cube or native MUSCAL, and whether "min edge" is over all tets or only the fault.
Every number here comes from the same code on both files.

Chunked throughout: the small-domain parent is 133.7 M tets and its connect
array alone is 4.3 GB if read whole.
"""
import argparse
import sys

import h5py
import numpy as np

sys.path.insert(0, "code")
from material import Material

FACE = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)]
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
BOX = (132000.0, 724000.0, 3490000.0, 4055000.0)
CH = 4_000_000

ap = argparse.ArgumentParser()
ap.add_argument("--mesh", action="append", required=True, help="label=path")
ap.add_argument("--cvm", required=True)
ap.add_argument("--muscal", required=True, help="REQUIRED -- without it Material "
                                                "silently scores on the deck cube")
ap.add_argument("--gate", type=float, default=0.8)
a = ap.parse_args()

mat = Material(a.cvm, a.muscal, box=BOX, source="muscal", verbose=False)
out = {}

for spec in a.mesh:
    label, path = spec.split("=", 1)
    f = h5py.File(path)
    V = f["geometry"][:]
    C = f["connect"]
    B = f["boundary"]
    NT = C.shape[0]
    vol = 0.0
    neg = 0
    emin = np.inf
    esum = []
    etamin = np.inf
    hist = np.zeros(5, np.int64)          # <0.05 <0.1 <0.3 <0.7 >=0.7
    nfail = 0
    worstf = np.inf
    nf = {1: 0, 3: 0, 5: 0}
    ftri = []
    for s in range(0, NT, CH):
        t = C[s:s + CH][:].astype(np.int64)
        b = B[s:s + CH][:].astype(np.int64)
        p = V[t]
        e = np.stack([np.linalg.norm(p[:, j] - p[:, i], axis=1) for i, j in PAIRS], 1)
        dx = e.max(1)
        emin = min(emin, float(e.min()))
        esum.append(np.median(e, axis=1)[::37])          # subsample for the median
        d = np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                      p[:, 3] - p[:, 0]) / 6.0
        vol += np.abs(d).sum()
        neg += int((d <= 0).sum())
        # Joe-Liu eta, the project convention
        ssum = (e ** 2).sum(1)
        eta = np.where(ssum > 0, 12.0 * np.cbrt((3.0 * np.abs(d)) ** 2) / ssum, 0.0)
        etamin = min(etamin, float(eta.min()))
        hist += np.array([int((eta < 0.05).sum()), int((eta < 0.1).sum()),
                          int((eta < 0.3).sum()), int((eta < 0.7).sum()),
                          int((eta >= 0.7).sum())], np.int64)
        vb = mat.at(p.mean(1))
        fr = np.where(vb > 0, vb / dx, 0.0)
        nfail += int((fr < a.gate).sum())
        worstf = min(worstf, float(fr.min()))
        for k in range(4):
            c = (b >> (8 * k)) & 0xFF
            for u in (1, 3, 5):
                nf[u] += int((c == u).sum())
            m = c == 3
            if m.any():
                ftri.append(t[m][:, list(FACE[k])])
        del t, b, p, e, dx, d, eta, vb, fr
    F = np.vstack(ftri)
    key = np.sort(F, axis=1)
    u, cnt = np.unique(key, axis=0, return_counts=True)
    P = V[u]
    farea = 0.5 * np.linalg.norm(np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]), axis=1).sum()
    out[label] = dict(
        tets=NT, verts=len(V), vol=vol, neg=neg, emin=emin,
        emed=float(np.median(np.concatenate(esum))),
        etamin=etamin, h=hist, nfail=nfail, worstf=worstf,
        ffac=len(u), single=int((cnt == 1).sum()), farea=farea,
        bc1=nf[1], bc3=nf[3], bc5=nf[5],
        ext=(V[:, 0].min(), V[:, 0].max(), V[:, 1].min(), V[:, 1].max(), V[:, 2].min()))
    del V, F, key, u, cnt, P, ftri
    print(f"[done] {label}", flush=True)

ks = list(out)
w = 26
def row(name, fn):
    print(f"{name:<{w}}" + "".join(f"{fn(out[k]):>24}" for k in ks))

print("\n" + "=" * (w + 24 * len(ks)))
print(f"{'':<{w}}" + "".join(f"{k:>24}" for k in ks))
print("=" * (w + 24 * len(ks)))
row("tets",              lambda d: f"{d['tets']:,}")
row("vertices",          lambda d: f"{d['verts']:,}")
row("volume (km3)",      lambda d: f"{d['vol']/1e9:,.0f}")
row("tets per 1000 km3", lambda d: f"{d['tets']/(d['vol']/1e9)*1000:,.0f}")
row("inverted tets",     lambda d: f"{d['neg']:,}")
print("-" * (w + 24 * len(ks)))
row("min edge (m)",      lambda d: f"{d['emin']:.2f}")
row("median edge (m)",   lambda d: f"{d['emed']:.0f}")
print("-" * (w + 24 * len(ks)))
row("eta min",           lambda d: f"{d['etamin']:.4f}")
row("eta < 0.05",        lambda d: f"{d['h'][0]:,}")
row("  as % of tets",    lambda d: f"{100*d['h'][0]/d['tets']:.3f} %")
row("eta < 0.1",         lambda d: f"{d['h'][1]:,}")
row("  as % of tets",    lambda d: f"{100*d['h'][1]/d['tets']:.3f} %")
row("eta < 0.3",         lambda d: f"{d['h'][2]:,}")
row("  as % of tets",    lambda d: f"{100*d['h'][2]/d['tets']:.3f} %")
row("eta >= 0.7",        lambda d: f"{100*d['h'][4]/d['tets']:.1f} %")
print("-" * (w + 24 * len(ks)))
row("fault facets",      lambda d: f"{d['ffac']:,}")
row("fault single-sided",lambda d: f"{d['single']:,}")
row("fault area (km2)",  lambda d: f"{d['farea']/1e6:,.2f}")
row("BC1 free surface",  lambda d: f"{d['bc1']:,}")
row("BC5 absorbing",     lambda d: f"{d['bc5']:,}")
print("-" * (w + 24 * len(ks)))
row("below 1 Hz p5",     lambda d: f"{d['nfail']:,}")
row("  as % of tets",    lambda d: f"{100*d['nfail']/d['tets']:.3f} %")
row("worst Vs/dx",       lambda d: f"{d['worstf']:.4f}")
row("worst resolved Hz", lambda d: f"{1.25*d['worstf']:.4f}")
print("=" * (w + 24 * len(ks)))
