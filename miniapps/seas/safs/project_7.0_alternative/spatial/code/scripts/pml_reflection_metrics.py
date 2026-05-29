#!/usr/bin/env python3
"""Phase 12.3 PML-effectiveness metrics for seas_spatial_dyn_driver SAFS runs.

Quantifies whether the convolutional PML (plan §Phase 12) removes the
post-`t_reflect` reflection contamination of the SAFS dynamic-rupture run and
classifies the LSW `V_max` departure at `t ~ 28 s` as reflection-driven vs an
intrinsic friction instability.

Input is the driver's stdout log (the `[DIAG] step ... V_max=...` and
`[DIAG-SIGN] step ... sigma_n_min=...` trace lines that `SEAS_DIAG_BLOWUP=1`
emits).  This tool is STANDARD-LIBRARY ONLY so it runs anywhere (login node,
laptop) without `conda activate pythonenv`; matplotlib is used only if present
and only when `--plot` is requested.

Metrics implemented here (parseable from the log alone):

  1. A/B V_max(t) overlay + departure classification          (plan req 1)
        --ab PML_ON.log PML_OFF.log
     Decision rule: if the `t ~ 28 s` departure disappears / strongly
     delays WITH PML it was reflection-driven; if it persists it is an
     intrinsic friction/LSW instability.  Either is a definitive result.

  2. Reflection-free-window smoothness                         (plan req 2)
        --smooth ONE.log --t-reflect 6.94
     Flags a kink in V_max(t) / sigma_n_min(t) at the reflection arrival.

  3. Big-box reference agreement (gold standard)               (plan req 3)
        --bigbox PML.log BIGBOX.log --window 14
     Accept the PML if max relative |V_max^PML - V_max^bigbox| /
     V_max^bigbox < `--tol` (default 0.05) over [0, window] s.

Metrics 4 (near-boundary seismogram) and 5 (interior energy monitor) require
ParaView volume/receiver fields, NOT the scalar log trace, so they are NOT
implemented here; `--list-metrics` prints the full plan list and marks which
ones this tool covers so coverage is never silently overstated.

Exit code is 0 when every requested metric meets its acceptance gate, 1
otherwise (suitable for an sbatch post-step `&&` chain).
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Log parsing
# ---------------------------------------------------------------------------

# [DIAG] step 1234 t=12.34 V_max=5.6 V_substep_max=... n_rupturing(V>0.5)=... max_slip=... m
_DIAG_RE = re.compile(
    r"^\[DIAG\]\s+step\s+(\d+)\s+t=([-\d.eE+]+)\s+V_max=([-\d.eE+]+)"
)
# [DIAG-SIGN] step 1234 t=12.34 ... | end-of-step: sigma_n_min=1.23e7 ...
_SIGN_RE = re.compile(
    r"^\[DIAG-SIGN\]\s+step\s+(\d+)\s+t=([-\d.eE+]+).*?sigma_n_min=([-\d.eE+]+)"
)


class Series:
    """A parsed V_max(t) (and optional sigma_n_min(t)) time series."""

    def __init__(self, path: Path):
        self.path = path
        self.t: list[float] = []
        self.vmax: list[float] = []
        self.t_sn: list[float] = []
        self.sigma_n_min: list[float] = []
        self._parse()

    def _parse(self) -> None:
        with open(self.path, "r", errors="replace") as fh:
            for line in fh:
                m = _DIAG_RE.match(line)
                if m:
                    self.t.append(float(m.group(2)))
                    self.vmax.append(float(m.group(3)))
                    continue
                s = _SIGN_RE.match(line)
                if s:
                    self.t_sn.append(float(s.group(2)))
                    self.sigma_n_min.append(float(s.group(3)))
        if not self.t:
            raise ValueError(
                f"{self.path}: no '[DIAG] step ... V_max=' lines found. "
                "Run the driver with SEAS_DIAG_BLOWUP=1 (V_max trace ON)."
            )

    @property
    def t_end(self) -> float:
        return self.t[-1]

    @property
    def vmax_peak(self) -> float:
        return max(self.vmax)

    def vmax_at(self, time: float) -> float:
        """Linear interpolation of V_max at `time` (clamped to the range)."""
        return _interp(self.t, self.vmax, time)


def _interp(xs: list[float], ys: list[float], x: float) -> float:
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    # xs is monotincreasing (step loop); binary-search the bracket.
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x:
            lo = mid
        else:
            hi = mid
    x0, x1 = xs[lo], xs[hi]
    y0, y1 = ys[lo], ys[hi]
    if x1 == x0:
        return y0
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


# ---------------------------------------------------------------------------
# Metric 1 — A/B overlay + departure classification
# ---------------------------------------------------------------------------

def _departure_time(s: Series, baseline_end: float = 7.0,
                     factor: float = 2.0) -> float | None:
    """First time V_max exceeds `factor` x its median over the reflection-free
    window [0, baseline_end].  Returns None if it never departs."""
    base = [v for ti, v in zip(s.t, s.vmax) if ti <= baseline_end]
    if not base:
        base = s.vmax[: max(1, len(s.vmax) // 4)]
    base_sorted = sorted(base)
    med = base_sorted[len(base_sorted) // 2]
    thresh = max(med * factor, med + 1e-12)
    for ti, v in zip(s.t, s.vmax):
        if ti > baseline_end and v > thresh:
            return ti
    return None


def metric_ab(pml_on: Path, pml_off: Path, t_reflect: float,
              csv_out: Path | None, plot: bool) -> tuple[bool, str]:
    on = Series(pml_on)
    off = Series(pml_off)
    dep_on = _departure_time(on, baseline_end=t_reflect)
    dep_off = _departure_time(off, baseline_end=t_reflect)

    lines = ["=== Metric 1: A/B V_max(t) overlay (plan req 1) ==="]
    lines.append(f"  PML ON : {pml_on}  (t_end={on.t_end:.3g}s, "
                 f"V_max_peak={on.vmax_peak:.4g} m/s)")
    lines.append(f"  PML OFF: {pml_off}  (t_end={off.t_end:.3g}s, "
                 f"V_max_peak={off.vmax_peak:.4g} m/s)")
    lines.append(f"  t_reflect (reflection-free window end) = {t_reflect:.3g}s")
    lines.append(f"  departure time  PML OFF: "
                 f"{'none' if dep_off is None else f'{dep_off:.3g}s'}")
    lines.append(f"  departure time  PML ON : "
                 f"{'none' if dep_on is None else f'{dep_on:.3g}s'}")

    # Classification (descriptive — both outcomes are definitive results).
    if dep_off is None:
        verdict = ("INCONCLUSIVE: no post-reflection departure in the PML-OFF "
                   "run either — extend tfinal past the t~28s onset.")
        ok = False
    elif dep_on is None:
        verdict = ("REFLECTION-DRIVEN: the departure present without PML "
                   "DISAPPEARS with PML.")
        ok = True
    elif dep_on > dep_off + 1.0:
        verdict = (f"REFLECTION-DRIVEN (delayed): PML pushes the departure "
                   f"from {dep_off:.3g}s to {dep_on:.3g}s.")
        ok = True
    else:
        verdict = ("INTRINSIC: the departure persists at nearly the same time "
                   "with PML -> friction/LSW instability, not reflection.")
        ok = True  # a definitive classification still 'passes' the metric
    lines.append(f"  >> CLASSIFICATION: {verdict}")

    if csv_out is not None:
        _write_overlay_csv(csv_out, on, off)
        lines.append(f"  overlay CSV written: {csv_out}")
    if plot:
        msg = _maybe_plot_overlay(on, off, t_reflect,
                                  (csv_out.with_suffix(".png")
                                   if csv_out else Path("pml_ab_vmax.png")))
        lines.append(f"  {msg}")
    return ok, "\n".join(lines)


def _write_overlay_csv(path: Path, on: Series, off: Series) -> None:
    # Union of sample times; interpolate both series onto it.
    times = sorted(set(on.t) | set(off.t))
    with open(path, "w") as fh:
        fh.write("t,V_max_pml_on,V_max_pml_off\n")
        for ti in times:
            fh.write(f"{ti:.6g},{on.vmax_at(ti):.6g},{off.vmax_at(ti):.6g}\n")


# ---------------------------------------------------------------------------
# Metric 2 — reflection-free-window smoothness (kink detector)
# ---------------------------------------------------------------------------

def _max_rel_jump(t: list[float], y: list[float],
                  t0: float, t1: float) -> float:
    """Largest |Δy|/(|y|+eps) step-to-step over the window [t0,t1]."""
    worst = 0.0
    eps = 1e-30
    for i in range(1, len(t)):
        if t[i] < t0 or t[i] > t1:
            continue
        dy = abs(y[i] - y[i - 1])
        worst = max(worst, dy / (abs(y[i - 1]) + eps))
    return worst


def metric_smooth(log: Path, t_reflect: float, halfwin: float,
                  jump_tol: float) -> tuple[bool, str]:
    s = Series(log)
    lines = ["=== Metric 2: reflection-free-window smoothness (plan req 2) ==="]
    lines.append(f"  log: {log}")
    win_lo, win_hi = t_reflect - halfwin, t_reflect + halfwin
    j_v = _max_rel_jump(s.t, s.vmax, win_lo, win_hi)
    lines.append(f"  V_max max relative step-jump in [{win_lo:.3g},"
                 f"{win_hi:.3g}]s = {j_v:.3g} (tol {jump_tol:.3g})")
    ok = j_v <= jump_tol
    if s.sigma_n_min:
        j_sn = _max_rel_jump(s.t_sn, s.sigma_n_min, win_lo, win_hi)
        lines.append(f"  sigma_n_min max relative step-jump = {j_sn:.3g}")
        ok = ok and (j_sn <= jump_tol)
    else:
        lines.append("  sigma_n_min: no [DIAG-SIGN] lines (skipped; enable "
                     "SEAS_DIAG_BLOWUP / the [DIAG-SIGN] tracker).")
    lines.append(f"  >> {'PASS' if ok else 'FAIL'}: no kink at the reflection "
                 "arrival." if ok else
                 f"  >> FAIL: a kink at t~{t_reflect:.3g}s suggests residual "
                 "reflection contamination.")
    return ok, "\n".join(lines)


# ---------------------------------------------------------------------------
# Metric 3 — big-box reference agreement
# ---------------------------------------------------------------------------

def metric_bigbox(pml: Path, bigbox: Path, window: float,
                  tol: float, nsamp: int) -> tuple[bool, str]:
    p = Series(pml)
    b = Series(bigbox)
    hi = min(window, p.t_end, b.t_end)
    lines = ["=== Metric 3: big-box reference agreement (plan req 3) ==="]
    lines.append(f"  PML run    : {pml}")
    lines.append(f"  big-box ref: {bigbox}")
    lines.append(f"  comparison window: [0, {hi:.3g}]s "
                 f"(reflection-free in the big box)")
    worst = 0.0
    worst_t = 0.0
    for k in range(nsamp + 1):
        ti = hi * k / nsamp
        vp = p.vmax_at(ti)
        vb = b.vmax_at(ti)
        denom = abs(vb)
        if denom < 1e-12:
            continue
        rel = abs(vp - vb) / denom
        if rel > worst:
            worst, worst_t = rel, ti
    lines.append(f"  max relative |V_max^PML - V_max^bigbox| / V_max^bigbox "
                 f"= {worst*100:.3g}%  at t={worst_t:.3g}s")
    lines.append(f"  acceptance tol = {tol*100:.3g}%")
    ok = worst < tol
    lines.append(f"  >> {'PASS' if ok else 'FAIL'}: PML "
                 f"{'matches' if ok else 'DEVIATES FROM'} the big-box truth.")
    return ok, "\n".join(lines)


# ---------------------------------------------------------------------------
# Optional plotting
# ---------------------------------------------------------------------------

def _maybe_plot_overlay(on: Series, off: Series, t_reflect: float,
                        png: Path) -> str:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:  # noqa: BLE001
        return ("matplotlib not available — skipped --plot (CSV is the "
                "portable artifact).")
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(on.t, on.vmax, label="PML ON", lw=1.3)
    ax.plot(off.t, off.vmax, label="PML OFF", lw=1.3, ls="--")
    ax.axvline(t_reflect, color="k", lw=0.8, ls=":",
               label=f"t_reflect={t_reflect:.3g}s")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("V_max [m/s]")
    ax.set_title("SAFS V_max(t): PML vs no-PML (Phase 12.3)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(png, dpi=130)
    return f"overlay plot written: {png}"


# ---------------------------------------------------------------------------
# Metric coverage table (honest scope — plan §12.3 has 5 metrics)
# ---------------------------------------------------------------------------

_METRICS = [
    ("1. A/B V_max(t) overlay + departure classification", True,
     "--ab PML_ON.log PML_OFF.log"),
    ("2. Reflection-free-window smoothness (kink detector)", True,
     "--smooth ONE.log --t-reflect 6.94"),
    ("3. Big-box reference agreement (< tol over window)", True,
     "--bigbox PML.log BIGBOX.log --window 14"),
    ("4. Near-boundary seismogram (reflected/incident peak ratio)", False,
     "requires ParaView receiver traces (volume/bulk fields)"),
    ("5. Interior energy monitor (decay, no t_reflect re-injection)", False,
     "requires ParaView volume-field integration (EnergyDensity)"),
]


def list_metrics() -> str:
    out = ["Plan §12.3 metrics and this tool's coverage:"]
    for name, covered, how in _METRICS:
        tag = "COVERED " if covered else "NOT here"
        out.append(f"  [{tag}] {name}\n             -> {how}")
    out.append("\nMetrics 4 & 5 need ParaView field post-processing (h5py / "
               "the .vtkhdf), out of scope for this log-only tool.")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(
        description="Phase 12.3 PML-effectiveness metrics from driver logs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--ab", nargs=2, metavar=("PML_ON", "PML_OFF"),
                    help="Metric 1: A/B V_max(t) overlay + classification.")
    ap.add_argument("--smooth", metavar="LOG",
                    help="Metric 2: kink check around t_reflect.")
    ap.add_argument("--bigbox", nargs=2, metavar=("PML", "BIGBOX"),
                    help="Metric 3: agreement vs the big-box reference.")
    ap.add_argument("--t-reflect", type=float, default=6.94,
                    help="Reflection arrival time [s] "
                         "(default 6.94 = 41.6km/5996m/s).")
    ap.add_argument("--window", type=float, default=14.0,
                    help="Metric 3 comparison window [s] (default 14).")
    ap.add_argument("--halfwin", type=float, default=2.0,
                    help="Metric 2 half-window around t_reflect [s].")
    ap.add_argument("--tol", type=float, default=0.05,
                    help="Metric 3 acceptance tol (default 0.05 = 5%%).")
    ap.add_argument("--jump-tol", type=float, default=0.5,
                    help="Metric 2 max relative step-jump (default 0.5).")
    ap.add_argument("--nsamp", type=int, default=400,
                    help="Metric 3 sample count over the window.")
    ap.add_argument("--csv", metavar="PATH", type=Path,
                    help="Metric 1: write the V_max overlay CSV here.")
    ap.add_argument("--plot", action="store_true",
                    help="Metric 1: also write a PNG overlay (needs matplotlib).")
    ap.add_argument("--list-metrics", action="store_true",
                    help="Print the plan §12.3 metric list + coverage and exit.")
    args = ap.parse_args(argv)

    if args.list_metrics:
        print(list_metrics())
        return 0

    if not (args.ab or args.smooth or args.bigbox):
        ap.error("nothing to do: pass --ab, --smooth, and/or --bigbox "
                 "(or --list-metrics).")

    reports: list[str] = []
    all_ok = True
    try:
        if args.ab:
            ok, rep = metric_ab(Path(args.ab[0]), Path(args.ab[1]),
                                args.t_reflect, args.csv, args.plot)
            reports.append(rep); all_ok = all_ok and ok
        if args.smooth:
            ok, rep = metric_smooth(Path(args.smooth), args.t_reflect,
                                    args.halfwin, args.jump_tol)
            reports.append(rep); all_ok = all_ok and ok
        if args.bigbox:
            ok, rep = metric_bigbox(Path(args.bigbox[0]), Path(args.bigbox[1]),
                                    args.window, args.tol, args.nsamp)
            reports.append(rep); all_ok = all_ok and ok
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    print("\n\n".join(reports))
    print(f"\n=== OVERALL: {'PASS' if all_ok else 'FAIL'} ===")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
