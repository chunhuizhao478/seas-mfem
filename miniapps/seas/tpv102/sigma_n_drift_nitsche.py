#!/usr/bin/env python3
"""TPV102 / T-13(b): normal-stress drift of a Nitsche run, two ways.

(a) per station: max over t of |sigma_n - 120 MPa| at the nine SCEC stations,
    ours beside every gold run and DRDG3D / PyLith (the references exist only
    at stations);
(b) fault-wide: over all (t, fault DOF) samples of the fault dumps
    <prefix>_fault_r*.csv, max |sigma_n - 120|/120 and its 99th / 99.9th
    percentiles, so one front sample does not stand for the whole field;
    the percentiles are repeated over the slipping samples (|V| > 1 mm/s),
    since the unruptured fault keeps sigma_n = 120 exactly.  The location of
    the maximum is printed.

Usage
  sigma_n_drift_nitsche.py --ours DIR [--prefix tpv102] [--gold DIR ...]
                           [--bench DIR] [--tmax T]
Dump columns (drivers/nitsche_dyn_driver.cpp): step t x y z cx cy cz d(3) V(3)
psi Xd Xs sig, sig = effective normal stress (compression +, Pa).
Scope: Chunhui Zhao, 2026-10-14 (P7).
"""
import argparse, glob, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_nitsche_scec import STATIONS, SCEC_NAME, load_ours, load_scec  # noqa: E402

SIG0 = 120.0  # MPa
SIG_COL = 17


def station_drift(d, tmax):
    m = (d["t"] <= tmax + 1e-9) & (d["n-stress"] != 0.0)   # skip all-zero t=0 rows
    return np.max(np.abs(np.abs(d["n-stress"][m]) - SIG0)) if m.any() else np.nan


def gold_label(g):
    g = g.rstrip("/")
    name = os.path.basename(os.path.dirname(g)) if os.path.basename(g) == "results" \
        else os.path.basename(g)
    return "_".join(name.split("_")[:2])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ours", required=True)
    ap.add_argument("--prefix", default="tpv102")
    ap.add_argument("--gold", nargs="*", default=[])
    ap.add_argument("--bench", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                    "benchmark_data"))
    ap.add_argument("--tmax", type=float, default=None,
                    help="window end (default: our last station time, capped at 12 s)")
    a = ap.parse_args()

    ours = {s: load_ours(p) for s in STATIONS
            for p in [os.path.join(a.ours, f"{a.prefix}_station_{s}.dat")] if os.path.exists(p)}
    tmax = a.tmax if a.tmax is not None else \
        (min(min(d["t"][-1] for d in ours.values()), 12.0) if ours else 12.0)
    print(f"# window 0 <= t <= {tmax:.4f} s")

    # ---------------- (a) per station ----------------
    cols = [("nitsche", lambda s: ours.get(s))]
    for g in a.gold:
        lab = gold_label(g)
        cols.append((lab, lambda s, g=g: load_ours(os.path.join(g, f"tpv102_station_{s}.dat"))
                     if os.path.exists(os.path.join(g, f"tpv102_station_{s}.dat")) else None))
    for code in ("drdg3d", "pylith"):
        def ref(s, code=code):
            f = os.path.join(a.bench, f"scec_{code}", f"tpv102_{code}_{SCEC_NAME[s]}.txt")
            return load_scec(f) if os.path.exists(f) else None
        cols.append((code, ref))
    print("\n# (a) max_t |sigma_n - 120 MPa| per station (MPa)")
    print(f"{'station':11s} " + " ".join(f"{c[0][:16]:>16s}" for c in cols))
    for s in STATIONS:
        cells = []
        for _, get in cols:
            d = get(s)
            cells.append(f"{station_drift(d, tmax):16.4f}" if d is not None else f"{'-':>16s}")
        print(f"{s:11s} " + " ".join(cells))

    # ---------------- (b) fault-wide ----------------
    files = sorted(glob.glob(os.path.join(a.ours, f"{a.prefix}_fault_r*.csv")))
    parts = []
    for f in files:
        if os.path.getsize(f) == 0:
            continue
        A = np.loadtxt(f, delimiter=",", ndmin=2, usecols=(1, 2, 3, 4, SIG_COL, 11, 12, 13))
        parts.append(A[A[:, 0] <= tmax + 1e-9])
    if not parts:
        print("\n# (b) no fault dumps found")
        return
    A = np.vstack(parts)
    rel = np.abs(A[:, 4] / 1e6 - SIG0) / SIG0
    i = int(np.argmax(rel))
    nt = np.unique(A[:, 0]).size
    print(f"\n# (b) fault-wide |sigma_n - 120|/120 over {A.shape[0]} (t, DOF) samples"
          f" ({nt} dump times)")
    print(f"  max    {rel[i]:.4e}  ({rel[i]*SIG0:.3f} MPa) at t={A[i,0]:.4f} s,"
          f" x={A[i,1]/1e3:.3f} km, depth={-A[i,3]/1e3:.3f} km")
    for q in (99.0, 99.9):
        v = np.percentile(rel, q)
        print(f"  p{q:<5g} {v:.4e}  ({v*SIG0:.3f} MPa)")
    # the unruptured fault keeps sigma_n = 120 exactly and dilutes the
    # percentiles; repeat them over the slipping samples only
    slipping = np.sqrt((A[:, 5:8] ** 2).sum(1)) > 1e-3
    if slipping.any():
        r = rel[slipping]
        print(f"  slipping samples (|V| > 1 mm/s): {int(slipping.sum())},"
              f" max {r.max():.4e} ({r.max()*SIG0:.3f} MPa)")
        for q in (99.0, 99.9):
            v = np.percentile(r, q)
            print(f"  p{q:<5g} {v:.4e}  ({v*SIG0:.3f} MPa)")


if __name__ == "__main__":
    main()
