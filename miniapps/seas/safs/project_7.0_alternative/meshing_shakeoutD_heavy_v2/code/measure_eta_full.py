#!/usr/bin/env python3
"""measure_eta_full.py -- one-line chunked Joe-Liu eta summary of a MEDIT mesh.

The per-pass measurement of the iterated -optim campaign, as a script so loop
drivers can call it and parse one stable line:

  ETA <mesh> tets=N min=X lt005=N lt01=N lt03=N edgemin=X tris=N
"""
import sys
import numpy as np, pandas as pd

sys.path.insert(0, sys.path[0] or '.')
from medit_hdr import medit_sections

F = sys.argv[1]
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
CH = 8_000_000

sec = medit_sections(F)
LV, NV = sec['Vertices']
LT, NT = sec['Tetrahedra']
NS = sec['Triangles'][1] if 'Triangles' in sec else 0
P = pd.read_csv(F, sep=r'\s+', header=None, skiprows=LV + 1, nrows=NV,
                usecols=range(3), dtype=np.float64, engine='c').to_numpy()
hist = np.zeros(3, np.int64)
emin = np.inf
etamin = np.inf
for s in range(0, NT, CH):
    n = min(CH, NT - s)
    T = pd.read_csv(F, sep=r'\s+', header=None, skiprows=LT + 1 + s, nrows=n,
                    usecols=range(4), dtype=np.int32, engine='c').to_numpy() - 1
    p = P[T]
    e = np.stack([np.linalg.norm(p[:, j] - p[:, i], axis=1) for i, j in PAIRS], 1)
    emin = min(emin, float(e.min()))
    d6 = np.abs(np.einsum('ij,ij->i', np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                          p[:, 3] - p[:, 0])) / 6.0
    ss = (e ** 2).sum(1)
    eta = np.where(ss > 0, 12.0 * np.cbrt((3.0 * d6) ** 2) / ss, 0.0)
    etamin = min(etamin, float(eta.min()))
    hist += np.array([int((eta < 0.05).sum()), int((eta < 0.1).sum()),
                      int((eta < 0.3).sum())])
    del T, p, e, d6, ss, eta
print(f'ETA {F} tets={NT} min={etamin:.4f} lt005={hist[0]} lt01={hist[1]} '
      f'lt03={hist[2]} edgemin={emin:.2f} tris={NS}', flush=True)
