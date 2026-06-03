#!/usr/bin/env python3
"""Compare on-fault sigma_n drift between two spatial-driver runs (the
fault-flux over-integration A/B).

The spatial driver writes one file per on-fault station:
    tpv102_station_<name>.dat
with a 3-line '#' header and columns:
    time slip1 slip2 V1 V2 tau1 tau2 sigma_n log10_theta
sigma_n is column index 7 (8th column).  For a symmetric planar strike-slip
fault sigma_n is analytically constant at 120 MPa; the drift is
max_t |sigma_n(t) - 120 MPa|.

Usage (from miniapps/seas):
    python3 jobs/tpv102_spatial/compare_sigma_n.py <OFF_dir> <ON_dir> [--plot]

Stdlib-only for the table; --plot additionally needs matplotlib (writes PNGs,
Agg backend, no display required).
"""
import sys, os, glob, argparse

SIGMA0 = 120.0e6      # Pa (analytic on-fault normal stress)
COL_T, COL_SN = 0, 7  # 0-based column indices: time, sigma_n


def load(path):
    t, sn = [], []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            p = s.split()
            if len(p) <= COL_SN:
                continue
            t.append(float(p[COL_T]))
            sn.append(float(p[COL_SN]))
    return t, sn


def drift(sn):
    dev = [abs(s - SIGMA0) for s in sn]
    return (max(dev) if dev else float("nan"),
            dev[-1] if dev else float("nan"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("off_dir", help="baseline (--fault-overint 0) output dir")
    ap.add_argument("on_dir",  help="over-int (--fault-overint 2) output dir")
    ap.add_argument("--plot", action="store_true", help="write per-station PNGs")
    args = ap.parse_args()

    off_files = sorted(glob.glob(os.path.join(args.off_dir,
                                              "tpv102_station_*.dat")))
    if not off_files:
        sys.exit(f"no tpv102_station_*.dat in {args.off_dir}")

    print(f"on-fault sigma_n drift  =  max_t |sigma_n - 120 MPa|   [MPa]")
    print(f"{'station':<16}{'OFF peak':>10}{'ON peak':>10}{'reduction':>11}"
          f"{'OFF final':>11}{'ON final':>10}")
    print("-" * 68)

    rows, worst = [], (None, 0.0)
    for off in off_files:
        name = os.path.basename(off)
        on = os.path.join(args.on_dir, name)
        if not os.path.exists(on):
            print(f"{name:<16}  (no match in ON dir — skipped)")
            continue
        to, sno = load(off)
        tn, snn = load(on)
        op, of_ = drift(sno)
        np_, nf = drift(snn)
        red = (op / np_) if np_ > 1e-9 else float("inf")
        tag = name.replace("tpv102_station_", "").replace(".dat", "")
        print(f"{tag:<16}{op/1e6:>10.3f}{np_/1e6:>10.3f}{red:>10.2f}x"
              f"{of_/1e6:>11.3f}{nf/1e6:>10.3f}")
        rows.append((tag, off, on))
        if op > worst[1]:
            worst = (tag, op)

    print("-" * 68)
    if worst[0]:
        print(f"worst baseline drift: station '{worst[0]}', "
              f"max|dSn| = {worst[1]/1e6:.3f} MPa")
    print("\nreduction = OFF_peak / ON_peak  (>1 means over-integration "
          "shrank the drift; the plan predicts >> 1).")

    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except Exception as e:
            sys.exit(f"--plot needs matplotlib: {e}")
        for tag, off, on in rows:
            to, sno = load(off); tn, snn = load(on)
            plt.figure(figsize=(7, 4))
            plt.axhline(120.0, color="k", lw=0.6, ls="--", label="120 MPa (analytic)")
            plt.plot(to, [s/1e6 for s in sno], label="baseline (overint 0)")
            plt.plot(tn, [s/1e6 for s in snn], label="over-int (overint 2)")
            plt.xlabel("time [s]"); plt.ylabel("sigma_n [MPa]")
            plt.title(f"on-fault sigma_n — station {tag}")
            plt.legend(); plt.tight_layout()
            out = f"sigma_n_{tag}.png"
            plt.savefig(out, dpi=120); plt.close()
            print(f"  wrote {out}")


if __name__ == "__main__":
    main()
