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
    python3 jobs/tpv102_spatial/compare_sigma_n.py <OFF_dir> <ON_dir> [RESAMPLE_dir] [--plot]
The optional 3rd dir is the over-int + --fault-resample arm (C); when given the
table adds the C peak/final and the C/A reduction.

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
    ap.add_argument("resample_dir", nargs="?", default=None,
                    help="optional over-int + --fault-resample dir (arm C)")
    ap.add_argument("--plot", action="store_true", help="write per-station PNGs")
    args = ap.parse_args()

    off_files = sorted(glob.glob(os.path.join(args.off_dir,
                                              "tpv102_station_*.dat")))
    if not off_files:
        sys.exit(f"no tpv102_station_*.dat in {args.off_dir}")

    has_rs = args.resample_dir is not None
    width = 76 if has_rs else 68

    print(f"on-fault sigma_n drift  =  max_t |sigma_n - 120 MPa|   [MPa]")
    if has_rs:
        print(f"{'station':<16}{'A off':>9}{'B ovr':>9}{'C +rs':>9}"
              f"{'B/A':>7}{'C/A':>7}{'A fin':>9}{'C fin':>9}")
    else:
        print(f"{'station':<16}{'OFF peak':>10}{'ON peak':>10}{'reduction':>11}"
              f"{'OFF final':>11}{'ON final':>10}")
    print("-" * width)

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
        tag = name.replace("tpv102_station_", "").replace(".dat", "")
        redB = (op / np_) if np_ > 1e-9 else float("inf")

        rs_path = os.path.join(args.resample_dir, name) if has_rs else None
        if rs_path and os.path.exists(rs_path):
            tr, snr = load(rs_path)
            rp, rf = drift(snr)
            redC = (op / rp) if rp > 1e-9 else float("inf")
            print(f"{tag:<16}{op/1e6:>9.2f}{np_/1e6:>9.2f}{rp/1e6:>9.2f}"
                  f"{redB:>6.1f}x{redC:>6.1f}x{of_/1e6:>9.2f}{rf/1e6:>9.2f}")
            rows.append((tag, off, on, rs_path))
        else:
            if has_rs:
                print(f"{tag:<16}{op/1e6:>9.2f}{np_/1e6:>9.2f}{'--':>9}"
                      f"{redB:>6.1f}x{'--':>7}{of_/1e6:>9.2f}{nf/1e6:>9.2f}"
                      "  (no C match)")
            else:
                print(f"{tag:<16}{op/1e6:>10.3f}{np_/1e6:>10.3f}{redB:>10.2f}x"
                      f"{of_/1e6:>11.3f}{nf/1e6:>10.3f}")
            rows.append((tag, off, on, None))
        if op > worst[1]:
            worst = (tag, op)

    print("-" * width)
    if worst[0]:
        print(f"worst baseline drift: station '{worst[0]}', "
              f"max|dSn| = {worst[1]/1e6:.3f} MPa")
    if has_rs:
        print("\nB/A = over-integration reduction; C/A = over-int+resample "
              "reduction (both > 1 => smaller drift; the plan predicts >> 1).")
    else:
        print("\nreduction = OFF_peak / ON_peak  (>1 means over-integration "
              "shrank the drift; the plan predicts >> 1).")

    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except Exception as e:
            sys.exit(f"--plot needs matplotlib: {e}")
        for tag, off, on, rs_path in rows:
            to, sno = load(off); tn, snn = load(on)
            plt.figure(figsize=(7, 4))
            plt.axhline(120.0, color="k", lw=0.6, ls="--", label="120 MPa (analytic)")
            plt.plot(to, [s/1e6 for s in sno], label="A baseline (overint 0)")
            plt.plot(tn, [s/1e6 for s in snn], label="B over-int (overint 2)")
            if rs_path:
                tr, snr = load(rs_path)
                plt.plot(tr, [s/1e6 for s in snr], label="C +resample")
            plt.xlabel("time [s]"); plt.ylabel("sigma_n [MPa]")
            plt.title(f"on-fault sigma_n — station {tag}")
            plt.legend(); plt.tight_layout()
            out = f"sigma_n_{tag}.png"
            plt.savefig(out, dpi=120); plt.close()
            print(f"  wrote {out}")


if __name__ == "__main__":
    main()
