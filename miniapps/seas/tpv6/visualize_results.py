#!/usr/bin/env python3
"""Plot TPV6/TPV7 on-fault, per-side station traces vs the SCEC DRDG3D reference.

Part C / C2 (unified bi-material plan).  Per the SCEC TPV6/7 spec (Part I), the
on-fault split-node stations report, on EACH side of the bi-material fault, the
absolute DISPLACEMENT and VELOCITY (NOT slip / slip-rate).

near vs far material (CONFIRMED — spec p.3 + the reference data):
    nearside = STRONG/FAST (vp=6000)  ;  farside = WEAK/SLOW (vp=3750 TPV6 / 5000 TPV7)
    -> the WEAK (far) side moves ~3-4x faster (the bi-material signature).

DRDG3D reference file columns (MKS; one file per station per side):
    t  h-disp  h-vel  h-stress  v-disp  v-vel  v-stress  n-disp  n-vel  n-stress
    h = along-strike, v = along-dip, n = fault-normal.
    NOTE: drdg3d n-stress is COMPRESSION-NEGATIVE; this repo is COMPRESSION-POSITIVE.

This run's per-side station traces come from the per-side station writer
`dynamic/tpv6_stations.hpp` (the documented-deferred Part-C remainder).  Until it
lands, this script plots the DRDG3D reference (so the comparison contract is fixed)
and overlays the run output if matching files are present (--out).

Usage:
    python3 visualize_results.py --ref benchmark_data/scec_drdg3d [--out <run_dir>]
                                 [--problem tpv6] [--save plots/]
"""
import argparse
import glob
import os
import sys

# drdg3d column index (0-based) by name.
COL = {"t": 0, "h-disp": 1, "h-vel": 2, "h-stress": 3, "v-disp": 4, "v-vel": 5,
       "v-stress": 6, "n-disp": 7, "n-vel": 8, "n-stress": 9}
STATIONS = ["x2_0_x3_0", "x2_-12_x3_0", "x2_12_x3_0", "x2_-12_x3_7.5", "x2_12_x3_7.5"]


def load_drdg3d(path):
    """Load a '#'-commented whitespace table -> list of float rows."""
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue  # the 'names' header row
    return rows


def col(rows, name):
    j = COL[name]
    return [r[j] for r in rows if len(r) > j]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", default="benchmark_data/scec_drdg3d",
                    help="DRDG3D reference dir")
    ap.add_argument("--out", default="", help="run output dir (overlay if present)")
    ap.add_argument("--problem", default="tpv6", choices=["tpv6", "tpv7"])
    ap.add_argument("--field", default="h-vel",
                    choices=["h-disp", "h-vel", "v-disp", "v-vel", "n-stress"])
    ap.add_argument("--save", default="", help="dir to write PNGs (else show)")
    args = ap.parse_args(argv)

    try:
        import matplotlib
        if args.save:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("ERROR: matplotlib required (conda activate pythonenv).", file=sys.stderr)
        return 2

    if not os.path.isdir(args.ref):
        print(f"ERROR: reference dir {args.ref} not found.", file=sys.stderr)
        return 2
    if args.save:
        os.makedirs(args.save, exist_ok=True)

    for st in STATIONS:
        fig, ax = plt.subplots(figsize=(8, 5))
        for side, c in (("nearside", "C0"), ("farside", "C1")):
            ref = os.path.join(args.ref, f"{args.problem}_drdg3d_{side}_{st}.txt")
            if not os.path.exists(ref):
                print(f"  (missing {os.path.basename(ref)})")
                continue
            rows = load_drdg3d(ref)
            if not rows:
                continue
            t, y = col(rows, "t"), col(rows, args.field)
            # drdg3d normal stress is COMPRESSION-NEGATIVE; this repo is
            # COMPRESSION-POSITIVE -> flip the reference n-stress so a run overlay aligns.
            if args.field == "n-stress":
                y = [-v for v in y]
            label = "near (STRONG)" if side == "nearside" else "far (WEAK)"
            ax.plot(t, y, c, label=f"drdg3d {label}")
            # Overlay this run if a matching per-side CSV exists.
            if args.out:
                hits = glob.glob(os.path.join(args.out, f"*{side}*{st}*")) \
                    or glob.glob(os.path.join(args.out, f"*{st}*{side}*"))
                for h in hits:
                    rr = load_drdg3d(h)
                    if rr:
                        ax.plot(col(rr, "t"), col(rr, args.field), c + "--",
                                label=f"run {label}")
        ax.set_xlabel("t [s]")
        ax.set_ylabel(args.field + (" [MPa]" if "stress" in args.field
                                    else " [m/s]" if "vel" in args.field else " [m]"))
        ax.set_title(f"{args.problem.upper()} on-fault {st} — {args.field} (near=STRONG, far=WEAK)")
        ax.legend(); ax.grid(True, alpha=0.3)
        fig.tight_layout()
        if args.save:
            out_png = os.path.join(args.save, f"{args.problem}_{st}_{args.field}.png")
            fig.savefig(out_png, dpi=120); plt.close(fig)
            print(f"  wrote {out_png}")
        else:
            plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
