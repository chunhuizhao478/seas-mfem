#!/usr/bin/env python3
"""v9.2.0 plan Step 2 — peak-value re-analysis vs Poisson bound.

Per the v9.2.0 plan Step 2 (REVIEW round-3 R-V92-E01), re-extract
per-cycle peak |sigma_n - 120 MPa| and peak |slip_dip| from the
existing job-7668434 outlier_report.json and classify whether the
normal-stress deviation SATURATES near the ~45 MPa Poisson bound
(rupture-extent artefact, not amplification) or grows unbounded
(real amplification — H-V92-G primary confirmed).

Why this is distinct from plot_peak_trajectory.py:
- plot_peak_trajectory reports `max_dev = max|value - median|`.
  median for sigma_n stays ~ 120 MPa (no drift), so `max_dev` is
  already effectively `max|sigma_n - 120 MPa|`.  This script makes
  that reference explicit and documents the 45 MPa Poisson bound
  as the classifier threshold so the verdict is reproducible.
- Emits a single JSON verdict + a focused two-panel plot
  (normal_stress vs slip_dip) so the result can be cited directly
  in the plan's Step 2 decision.

Poisson bound derivation (plan §2.4):
    TPV102 has Poisson's ratio ν = λ/(2(λ+μ)) = 0.25.
    Mode-II rupture-tip normal-stress concentration |Δσ_n| has the
    upper envelope |Δσ_n| ≤ |Δτ| · ν / (1 - ν) ≈ 0.33 · |Δτ|.
    Peak shear perturbation |Δτ| ≈ τ_dyn_drop ≈ 2 · (τ_s - τ_d) ≈
    2 · (0.677 - 0.525) · σ_n ≈ 0.30 · 120 MPa = 36 MPa.
    So the Poisson envelope is |Δσ_n| ≲ 12-15 MPa in a
    quasi-static analysis; the 45 MPa figure in the plan is the
    full-dynamic upper bound including radiation.  Either way,
    a SATURATION verdict requires |Δσ_n|_peak ≲ 45 MPa.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


POISSON_BOUND_PA = 45.0e6   # plan §2.4 full-dynamic Mode-II envelope
SIGMA_N0_PA      = 120.0e6
SCEC_SLIP_DIP_SPEC_M = 1.0e-3   # spec ceiling for |slip_dip|


def extract_trajectories(report_path: Path) -> dict:
    with open(report_path) as fh:
        r = json.load(fh)

    times, peak_dsn, peak_sd = [], [], []
    med_sn_min, med_sn_max = None, None
    for cyc in r["cycles"]:
        t = cyc["time"]
        s_sn = cyc["field_stats"].get("normal_stress")
        s_sd = cyc["field_stats"].get("slip_dip")
        if s_sn is None or s_sd is None:
            continue
        # Per-cycle peak |sigma_n - 120 MPa|:
        #   upper bound = |median - 120e6| + max_dev
        # median is ~120 MPa throughout the run (no drift), so this
        # reduces to max_dev to high accuracy.
        med_sn = s_sn.get("median", SIGMA_N0_PA)
        dev_sn = abs(med_sn - SIGMA_N0_PA) + s_sn.get("max_dev", 0.0)
        # Per-cycle peak |slip_dip|:
        #   slip_dip median ~ 0, so max_dev is a tight bound.
        med_sd = s_sd.get("median", 0.0)
        dev_sd = abs(med_sd) + s_sd.get("max_dev", 0.0)

        times.append(t)
        peak_dsn.append(dev_sn)
        peak_sd.append(dev_sd)
        med_sn_min = med_sn if med_sn_min is None else min(med_sn_min, med_sn)
        med_sn_max = med_sn if med_sn_max is None else max(med_sn_max, med_sn)

    return {
        "times": times,
        "peak_dsn": peak_dsn,          # |sigma_n - 120 MPa| per cycle
        "peak_sd":  peak_sd,           # |slip_dip| per cycle
        "median_sn_range": (med_sn_min, med_sn_max),
        "n_cycles": len(times),
        "report_path": str(report_path),
    }


def classify(traj: dict) -> dict:
    times = traj["times"]
    dsn = traj["peak_dsn"]
    sd  = traj["peak_sd"]
    if not times:
        return {"verdict": "NO-DATA"}

    n = len(dsn)
    max_dsn = max(dsn)
    end_dsn = dsn[-1]
    argmax_dsn = dsn.index(max_dsn)
    t_max = times[argmax_dsn]
    t_end = times[-1]

    # Late-run mean over last 25 %.
    cut = max(1, int(0.75 * n))
    late_mean_dsn = sum(dsn[cut:]) / max(1, n - cut)
    late_ratio = late_mean_dsn / max_dsn if max_dsn > 0 else 0.0
    at_end = (t_max / t_end) > 0.80 if t_end > 0 else False

    # Classification per plan Step 2:
    #   SATURATES near 45 MPa  ==> R-V92-E01 correct; rupture-extent misread
    #   GROWS unbounded       ==> real amplification; H-V92-G mechanism open
    verdict = "UNCLEAR"
    reason = ""
    if max_dsn < POISSON_BOUND_PA * 1.1:    # within 10 % of bound
        verdict = "SATURATES-BELOW-POISSON-BOUND"
        reason = (f"max |sigma_n - 120 MPa| = {max_dsn/1e6:.1f} MPa "
                  f"≲ 45 MPa Poisson bound (plan §2.4). "
                  "R-V92-E01 correct: outlier-count growth reflected "
                  "rupture extent, not amplification. §18 H-V92-G "
                  "verdict is not supported by peak-value evidence.")
    elif at_end and late_ratio > 0.85:
        verdict = "MONOTONIC-GROWTH"
        reason = (f"max |sigma_n - 120 MPa| = {max_dsn/1e6:.1f} MPa, "
                  f"{max_dsn / POISSON_BOUND_PA:.1f}× Poisson bound; "
                  "late-run mean tracks max; peak at end. "
                  "Real amplification — H-V92-G remains a viable "
                  "primary mechanism; need §5.1 bulk SYY probe to "
                  "localise the amplifier.")
    elif late_ratio > 0.80:
        verdict = "SATURATES-ABOVE-POISSON-BOUND"
        reason = (f"max |sigma_n - 120 MPa| = {max_dsn/1e6:.1f} MPa "
                  f"≫ 45 MPa Poisson bound BUT plateaus: "
                  f"late_mean/max = {late_ratio:.2f}. "
                  "Something bounds the amplification below runaway; "
                  "may be a nonlinear fault-strength ceiling. "
                  "Investigate whether the cap comes from friction "
                  "saturation or a numerical clipping path.")
    else:
        verdict = "GROWS-THEN-DECAYS"
        reason = (f"max |sigma_n - 120 MPa| = {max_dsn/1e6:.1f} MPa "
                  f"at t = {t_max:.2f} s; late_mean/max = "
                  f"{late_ratio:.2f}. Not consistent with "
                  "rupture-extent-tracking OR pure amplification. "
                  "Inspect plot for reflection-spike signature.")

    return {
        "verdict":         verdict,
        "reason":          reason,
        "max_dsn_pa":      max_dsn,
        "max_dsn_mpa":     max_dsn / 1e6,
        "end_dsn_pa":      end_dsn,
        "end_dsn_mpa":     end_dsn / 1e6,
        "max_slip_dip_m":  max(sd),
        "end_slip_dip_m":  sd[-1],
        "poisson_bound_mpa": POISSON_BOUND_PA / 1e6,
        "sigma_n0_mpa":    SIGMA_N0_PA / 1e6,
        "spec_slip_dip_m": SCEC_SLIP_DIP_SPEC_M,
        "slip_dip_factor_over_spec": max(sd) / SCEC_SLIP_DIP_SPEC_M
                                     if SCEC_SLIP_DIP_SPEC_M > 0 else None,
        "late_ratio":      late_ratio,
        "t_max":           t_max,
        "t_end":           t_end,
        "median_sn_range_mpa": [x / 1e6 for x in traj["median_sn_range"]],
    }


def emit_plot(traj: dict, cls: dict, out_path: Path) -> bool:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False

    t = traj["times"]
    dsn = [x / 1e6 for x in traj["peak_dsn"]]
    sd  = traj["peak_sd"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    ax1.plot(t, dsn, "-", color="C3", lw=1.6, label=r"peak $|\sigma_n - 120\,\mathrm{MPa}|$")
    ax1.axhline(POISSON_BOUND_PA / 1e6, color="k", ls="--", lw=1.0,
                label=r"45 MPa Poisson bound (plan §2.4)")
    ax1.set_xlabel("t (s)")
    ax1.set_ylabel("MPa")
    ax1.set_title(f"normal_stress deviation — {cls['verdict']}")
    ax1.legend(fontsize=9)
    ax1.grid(alpha=0.3)

    ax2.plot(t, sd, "-", color="C0", lw=1.6, label=r"peak $|\mathrm{slip}_{\mathrm{dip}}|$")
    ax2.axhline(SCEC_SLIP_DIP_SPEC_M, color="k", ls="--", lw=1.0,
                label="SCEC spec ceiling = 1 mm")
    ax2.set_xlabel("t (s)")
    ax2.set_ylabel("m")
    ax2.set_yscale("symlog", linthresh=1e-6)
    ax2.set_title(r"slip_dip — grows to %.2f m (%.0f$\times$ spec)" %
                  (cls["max_slip_dip_m"], cls["slip_dip_factor_over_spec"]))
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.3)

    fig.suptitle("v9.2.0 Step 2 re-analysis — peak vs Poisson bound", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    return True


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--report", type=Path,
                   default=Path(__file__).parent / "outlier_report.json")
    p.add_argument("--out-dir", type=Path,
                   default=Path(__file__).parent)
    p.add_argument("--no-plot", action="store_true")
    a = p.parse_args(argv)

    if not a.report.exists():
        print(f"ERROR: report not found: {a.report}", file=sys.stderr)
        return 2

    traj = extract_trajectories(a.report)
    cls = classify(traj)

    a.out_dir.mkdir(parents=True, exist_ok=True)
    json_path = a.out_dir / "step2_peak_vs_poisson_bound.json"
    with open(json_path, "w") as fh:
        json.dump({"trajectory_meta": {"n_cycles": traj["n_cycles"],
                                       "t_end": traj["times"][-1] if traj["times"] else None,
                                       "report_path": traj["report_path"],
                                       "median_sn_range_mpa": [x / 1e6 for x in traj["median_sn_range"]]},
                   "classification": cls},
                  fh, indent=2)

    if not a.no_plot:
        plot_path = a.out_dir / "step2_peak_vs_poisson_bound.png"
        ok = emit_plot(traj, cls, plot_path)
        if ok:
            cls["plot_path"] = str(plot_path)

    # Console verdict.
    print("=" * 72)
    print("TPV102 v9.2.0 Step 2 — peak |sigma_n - 120 MPa| vs Poisson bound")
    print("=" * 72)
    print(f"Report:     {traj['report_path']}")
    print(f"Cycles:     {traj['n_cycles']}  ({traj['times'][-1]:.2f} s end)")
    print(f"sigma_n median range: "
          f"{traj['median_sn_range'][0]/1e6:.3f} - "
          f"{traj['median_sn_range'][1]/1e6:.3f} MPa "
          "(stays at 120 MPa reference -> max_dev = deviation-from-ref)")
    print()
    print(f"Peak |sigma_n - 120 MPa|: {cls['max_dsn_mpa']:.1f} MPa   "
          f"(Poisson bound = {cls['poisson_bound_mpa']:.1f} MPa)")
    print(f"End  |sigma_n - 120 MPa|: {cls['end_dsn_mpa']:.1f} MPa")
    print(f"Peak |slip_dip|         : {cls['max_slip_dip_m']:.4f} m "
          f"(SCEC spec: 1 mm -> {cls['slip_dip_factor_over_spec']:.0f}x spec)")
    print(f"Late-run / max ratio    : {cls['late_ratio']:.3f}")
    print()
    print(f"Verdict: {cls['verdict']}")
    import textwrap
    for line in textwrap.wrap(cls["reason"], 70):
        print(f"  {line}")
    print("=" * 72)
    print(f"JSON:  {json_path}")
    if "plot_path" in cls:
        print(f"Plot:  {cls['plot_path']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
