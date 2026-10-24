#!/usr/bin/env python3
"""Spectral radius of the multirate-leapfrog sync-interval map: profiling plan
item 7 (local time stepping), phase P1, gate (c) of lts-kernels-design-rev1.

For r = r0, r0 + 0.01, ..., r1, `seas_test_nitsche_lts --phi r file` writes the
map of one sync interval on the graded box (cluster 0 at dt0 = r dt_CFL), and
numpy gives its eigenvalues.  The one-cluster map at the same dt0 (--gts,
central difference, stable for r < 1 with every |lambda| = 1) is the baseline:
its max |lambda| - 1 is the eigenvalue round-off.

P2: --kv eta_bar and --absorb give the map Kelvin-Voigt damping and the
absorbing faces of the P2 boxes.  The binary writes the map of the
recursion's own state (u, v after the start kicks), so the spectrum is the
stepping's also with damping; damped modes have |lambda| < 1, and the
one-cluster map is then unstable above its damped limit.

usage: python3 nitsche_lts_eig.py ./seas_test_nitsche_lts [--order p]
           [--r0 0.30] [--r1 1.00] [--jobs 4] [--tmp dir] [--kv eta] [--absorb]
"""
import argparse
import os
import subprocess
import sys
import tempfile

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402


def spectral_radius(binary, order, r, tmp, gts, extra):
    path = os.path.join(tmp, f"phi_p{order}_r{r:.2f}{'_gts' if gts else ''}.bin")
    cmd = [binary, "--phi", f"{r:.2f}", path, "--order", str(order)] + (["--gts"] if gts else []) + extra
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    with open(path, "rb") as f:
        n2 = int(np.frombuffer(f.read(8), dtype=np.int64)[0])
        phi = np.frombuffer(f.read(), dtype=np.float64).reshape(n2, n2).T
    os.remove(path)
    lam = np.abs(np.linalg.eigvals(phi))
    return r, gts, n2, lam.max() - 1.0, int((lam > 1.0 + 1e-10).sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("binary")
    ap.add_argument("--order", type=int, default=1)
    ap.add_argument("--r0", type=float, default=0.30)
    ap.add_argument("--r1", type=float, default=1.00)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--tmp", default=None)
    ap.add_argument("--kv", type=float, default=0.0)
    ap.add_argument("--absorb", action="store_true")
    a = ap.parse_args()
    extra = (["--kv", str(a.kv)] if a.kv > 0 else []) + (["--absorb"] if a.absorb else [])
    tmp = a.tmp or tempfile.mkdtemp(prefix="lts_eig_")
    rs = [round(a.r0 + 0.01 * k, 2) for k in range(int(round((a.r1 - a.r0) / 0.01)) + 1)]
    tasks = [(r, g) for r in rs for g in (False, True)]
    res = {}
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(spectral_radius, os.path.abspath(a.binary), a.order, r, tmp, g, extra)
                for r, g in tasks]
        for f in futs:
            r, g, n2, d, cnt = f.result()
            res[(r, g)] = (d, cnt)
            print(f"  r {r:.2f} {'one cluster' if g else 'clusters   '}: map {n2}x{n2}, "
                  f"max|lambda| - 1 = {d:.3e}, {cnt} above 1 + 1e-10", flush=True)
    print(f"p = {a.order}{', Kelvin-Voigt ' + str(a.kv) if a.kv > 0 else ''}"
          f"{', absorbing faces' if a.absorb else ''}: r, max|lambda| - 1 (clusters), (one cluster)")
    bad = []
    for r in rs:
        d, cnt = res[(r, False)]
        dg, _ = res[(r, True)]
        print(f"  {r:.2f}  {d: .3e}  {dg: .3e}")
        if d > 1e-10:
            bad.append(r)
    print("r with max|lambda| > 1 + 1e-10 (clusters): " + (", ".join(f"{r:.2f}" for r in bad) or "none"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
