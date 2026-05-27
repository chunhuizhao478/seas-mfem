#!/usr/bin/env python3
"""
nucleation_patch_strength.py — compute the rate-and-state SHEAR STRENGTH inside
the gradual_overstress nucleation patch from a SAFS fault.vtkhdf, and compare it
to the applied overstress (diagnose "is the overstress strong enough?").

Rate-and-state shear strength (regularized aging law, the form the driver uses):

    tau_str = sigma_n_eff * f,    f = a * asinh( V / (2 V0) * exp(psi / a) )

with sigma_n_eff = normal_stress (already effective; compression positive),
V = |slip_rate|, psi = state_variable, a = param_a, V0 the reference rate.
All of these are PointData fields in the fault.vtkhdf except V0 (a scalar
config constant).

The VTKHDF stores each PointData field as per-step blocks concatenated
([step0 all points][step1 all points]...); step k for point p is
field[offset_k + p] where offset_k = Steps/PointDataOffsets/<field>[k].

Usage:
    conda activate pythonenv   # (python, not python3 — has h5py/numpy)
    python nucleation_patch_strength.py FAULT_VTKHDF \
        [--center-x 606971 --center-y 3707270 --center-z -4965.62] \
        [--radius-m 4000] [--delta-tau-pa 20e6] [--V0 1.0e-6] \
        [--step 0] [--also-step -1]
"""
from __future__ import annotations
import argparse
import numpy as np
import h5py


def _load_step(f, field, step, npts):
    off = int(f[f"VTKHDF/Steps/PointDataOffsets/{field}"][step])
    return f[f"VTKHDF/PointData/{field}"][off:off + npts]


def _stats(name, arr, unit=""):
    a = np.asarray(arr, dtype=np.float64)
    print(f"  {name:22s} min={a.min():.6g}{unit}  mean={a.mean():.6g}{unit}  "
          f"max={a.max():.6g}{unit}")


def analyze(path, cx, cy, cz, radius, delta_tau, V0, step, npts_hint=None):
    f = h5py.File(path, "r")
    pts = f["VTKHDF/Points"][:]                       # (Npt, 3)
    npts = pts.shape[0]
    nsteps = f["VTKHDF/Steps/Values"].shape[0]
    times = f["VTKHDF/Steps/Values"][:]
    step = step % nsteps
    t = float(times[step])

    # Patch = points within `radius` (3-D) of the nucleation centre.
    d = np.sqrt((pts[:, 0] - cx) ** 2 + (pts[:, 1] - cy) ** 2
                + (pts[:, 2] - cz) ** 2)
    patch = d <= radius
    npatch = int(patch.sum())

    print(f"fault.vtkhdf : {path}")
    print(f"points       : {npts}   steps : {nsteps}   "
          f"t-range : [{times.min():.4g}, {times.max():.4g}] s")
    print(f"patch centre : ({cx:.0f}, {cy:.0f}, {cz:.1f})  radius {radius:.0f} m  "
          f"depth ~{abs(cz)/1000:.2f} km")
    print(f"patch points : {npatch}")
    print(f"--- step {step}  (t = {t:.6g} s) ---")
    if npatch == 0:
        print("  (no points within the patch radius — check the centre/radius)")
        return

    sig = _load_step(f, "normal_stress", step, npts)[patch]       # Pa, eff, +comp
    ts = _load_step(f, "traction_strike", step, npts)[patch]      # Pa
    td = _load_step(f, "traction_dip", step, npts)[patch]         # Pa
    vs = _load_step(f, "slip_rate_strike", step, npts)[patch]     # m/s
    vd = _load_step(f, "slip_rate_dip", step, npts)[patch]        # m/s
    psi = _load_step(f, "state_variable", step, npts)[patch]      # -
    a = _load_step(f, "param_a", step, npts)[patch]               # -

    V = np.sqrt(vs ** 2 + vd ** 2)
    tau = np.sqrt(ts ** 2 + td ** 2)
    # Regularized RS friction coefficient and strength.
    f_coef = a * np.arcsinh(V / (2.0 * V0) * np.exp(psi / a))
    tau_str = sig * f_coef                                        # Pa

    _stats("sigma_n_eff", sig / 1e6, " MPa")
    _stats("|slip_rate| V", V, " m/s")
    _stats("state psi", psi)
    _stats("param a", a)
    _stats("RS friction f", f_coef)
    _stats("SHEAR STRENGTH tau_str", tau_str / 1e6, " MPa")
    _stats("|traction| (mobilized)", tau / 1e6, " MPa")

    # Reference peak strength f0*sigma_n (f0 = 0.6) for context.
    f0 = 0.6
    print(f"  reference f0*sigma_n (f0={f0}) : "
          f"{(f0*sig).min()/1e6:.4g} / {(f0*sig).mean()/1e6:.4g} / "
          f"{(f0*sig).max()/1e6:.4g} MPa  (min/mean/max)")

    # Overstress comparison.
    str_mean = float(np.mean(tau_str))
    print(f"--- overstress vs strength (patch mean) ---")
    print(f"  applied overstress  delta_tau     = {delta_tau/1e6:.4g} MPa")
    print(f"  shear strength      tau_str(mean)  = {str_mean/1e6:.4g} MPa")
    if str_mean != 0.0:
        print(f"  ratio delta_tau / tau_str           = "
              f"{delta_tau/str_mean:.4f}")
    print(f"  driving |tau|+delta_tau (patch mean) = "
          f"{(float(np.mean(tau))+delta_tau)/1e6:.4g} MPa")
    print(f"  overshoot (|tau|+delta_tau - tau_str), patch mean = "
          f"{(float(np.mean(tau))+delta_tau-str_mean)/1e6:.4g} MPa")
    f.close()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fault_vtkhdf")
    ap.add_argument("--center-x", type=float, default=606971.0)
    ap.add_argument("--center-y", type=float, default=3707270.0)
    ap.add_argument("--center-z", type=float, default=-4965.62)
    ap.add_argument("--radius-m", type=float, default=4000.0)
    ap.add_argument("--delta-tau-pa", type=float, default=20.0e6)
    ap.add_argument("--V0", type=float, default=1.0e-6)
    ap.add_argument("--step", type=int, default=0)
    ap.add_argument("--also-step", type=int, default=None,
                    help="optionally analyze a second step (e.g. -1 for last)")
    args = ap.parse_args(argv)

    analyze(args.fault_vtkhdf, args.center_x, args.center_y, args.center_z,
            args.radius_m, args.delta_tau_pa, args.V0, args.step)
    if args.also_step is not None:
        print()
        analyze(args.fault_vtkhdf, args.center_x, args.center_y, args.center_z,
                args.radius_m, args.delta_tau_pa, args.V0, args.also_step)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
