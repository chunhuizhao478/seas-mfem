#!/usr/bin/env python3
"""Unit tests for fault_local_overpressure (PLAN Phase 2 acceptance).
Run: /Users/chunhuizhao/miniforge/envs/pythonenv/bin/python test_overpressure.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fault_local_overpressure import (
    build_overpressure_field, save_overpressure_field, load_overpressure_field,
    delta_pp_mpa, apply_overpressure_to_tensor, floor_mu, strike_distance_km)
from project_csm_stress_to_vtu import resolve_tractions

RNG = np.random.default_rng(0)


def test_isotropy_tau_unchanged_sigma_n_drops():
    """Subtracting DeltaPp*I drops sigma_n_eff by exactly DeltaPp and leaves tau
    unchanged (the whole physical basis of approach 2)."""
    n = 200
    # random symmetric comp+ tensors, random unit normals/strikes/dips frame
    A = RNG.normal(size=(n, 3, 3))
    sig = 0.5 * (A + np.transpose(A, (0, 2, 1))) * 20.0 + \
        np.eye(3)[None] * 80.0                       # diag-dominant, comp+
    nrm = RNG.normal(size=(n, 3)); nrm /= np.linalg.norm(nrm, axis=1)[:, None]
    # build an orthonormal (strike, dip) frame per facet
    up = np.array([0.0, 0.0, 1.0])
    s = np.cross(np.broadcast_to(up, nrm.shape), nrm)
    s /= np.linalg.norm(s, axis=1)[:, None]
    d = np.cross(s, nrm)
    dpp = RNG.uniform(1.0, 30.0, size=n)
    r0 = resolve_tractions(sig, s, d, nrm, 0.0)
    sig2 = apply_overpressure_to_tensor(sig, dpp)
    r1 = resolve_tractions(sig2, s, d, nrm, 0.0)
    dsn = (r0["sigma_n_eff"] - r1["sigma_n_eff"])
    assert np.allclose(dsn, dpp, atol=1e-9), \
        f"sigma_n_eff did not drop by DeltaPp: max err {np.abs(dsn-dpp).max():.2e}"
    assert np.allclose(r0["tau_magnitude"], r1["tau_magnitude"], atol=1e-9), \
        "tau changed under isotropic overpressure (must be invariant)"
    print("PASS test_isotropy_tau_unchanged_sigma_n_drops "
          f"(max d(sigma_n)-DPp {np.abs(dsn-dpp).max():.2e}, "
          f"max d(tau) {np.abs(r0['tau_magnitude']-r1['tau_magnitude']).max():.2e})")


def test_off_is_identity():
    """delta_pp_mpa returns 0 outside the field; apply with dpp=0 is identity."""
    sig = np.tile(np.diag([80.0, 70.0, 35.0]), (5, 1, 1))
    out = apply_overpressure_to_tensor(sig, np.zeros(5))
    assert np.array_equal(out, sig), "dpp=0 overpressure is not identity"
    print("PASS test_off_is_identity")


def _toy_field():
    """A toy gate: 1 km of facets at s in [40,71], depth 3..12 km, two
    populations per depth (a clamped low-mu core + a well-oriented high-mu)."""
    az, hypo = 314.0, (0.0, 0.0)
    su = np.array([np.sin(np.radians(az)), np.cos(np.radians(az))])
    s_targets = np.linspace(42.0, 69.0, 60)
    zk = np.linspace(-11.0, -4.0, 40)
    S, Z = np.meshgrid(s_targets, zk, indexing="ij")
    s_km = S.ravel(); depth = -Z.ravel() * 1000.0
    # place x,y so strike_distance recovers s_km (x,y along su)
    xy = np.outer(s_km * 1000.0, su)
    x, y = xy[:, 0], xy[:, 1]
    Sv_total = 26000.0 * 9.81 * depth / 1e6 * 1.0    # ~ rho g z (MPa)
    Pp_h = 1000.0 * 9.81 * depth / 1e6
    sn = 0.6 * Sv_total + 40.0                        # clamped (high sigma_n)
    tau = 0.12 * sn                                   # mu_app ~ 0.12 (gate core)
    # well-oriented (high mu_app) facets only at the band EDGES (a realistic
    # coherent gate: low-mu core, transitioning to higher mu_app at the edges)
    hi = (s_km < 44.0) | (s_km > 67.0)
    tau[hi] = 0.29 * sn[hi]
    seis = (Z.ravel() > -12) & (Z.ravel() < -3)
    return (x, y, s_km, depth, sn, tau, Sv_total, Pp_h, seis, az, hypo)


def test_build_never_supercritical_and_roundtrips(tmp="/tmp/_op_test.npz"):
    x, y, s_km, depth, sn, tau, Sv, Pp, seis, az, hypo = _toy_field()
    mu_s = 0.318
    field = build_overpressure_field(
        s_km, depth, sn, tau, Sv, Pp, seis, [(40.0, 71.0)], mu_s, az, hypo,
        mu_d=0.10, s_target=1.7, lambda_max=0.9)
    # sample DeltaPp at the same facets, apply, recompute mu_app; assert (a) no
    # facet exceeds mu_s (never supercritical) and (b) the LOW cores cross floor.
    dpp = delta_pp_mpa(field, x, y, -depth)           # z = -depth
    sn_new = sn - dpp
    assert np.all(sn_new > 0.4), "overpressure drove sigma_n_eff to tension"
    mu_new = tau / sn_new
    band_seis = seis
    assert np.nanmax(mu_new[band_seis]) <= mu_s + 1e-6, \
        f"supercritical: max mu_app {np.nanmax(mu_new[band_seis]):.3f} > mu_s {mu_s}"
    floor = floor_mu(mu_s, 0.10, 1.7)
    # well inside the homogeneous low-mu core (away from edge transition +
    # neighborhood-min reach), the cores should reach the floor
    core = band_seis & (s_km > 48) & (s_km < 63) & (tau < 0.2 * sn)
    frac = float(np.mean(mu_new[core] >= floor - 1e-3))
    assert frac > 0.9, f"only {100*frac:.0f}% of low cores reached floor {floor:.3f}"
    # round-trip save/load identical
    save_overpressure_field(tmp, field)
    f2 = load_overpressure_field(tmp)
    assert np.array_equal(field["DPp"], f2["DPp"])
    assert np.allclose(delta_pp_mpa(f2, x, y, -depth), dpp)
    print(f"PASS test_build_never_supercritical_and_roundtrips "
          f"(floor {floor:.3f}, low-core cross {100*frac:.0f}%, "
          f"max mu_app after {np.nanmax(mu_new[band_seis]):.3f} <= mu_s {mu_s})")


def test_outside_field_is_zero():
    x, y, s_km, depth, sn, tau, Sv, Pp, seis, az, hypo = _toy_field()
    field = build_overpressure_field(
        s_km, depth, sn, tau, Sv, Pp, seis, [(40.0, 71.0)], 0.318, az, hypo)
    # a point at s=200 km (no gate there) must get DeltaPp = 0
    su = np.array([np.sin(np.radians(az)), np.cos(np.radians(az))])
    far = 200.0 * 1000.0 * su
    val = delta_pp_mpa(field, np.array([far[0]]), np.array([far[1]]),
                       np.array([-6000.0]))
    assert float(val[0]) == 0.0, f"overpressure leaked outside the gate: {val}"
    print("PASS test_outside_field_is_zero")


if __name__ == "__main__":
    test_isotropy_tau_unchanged_sigma_n_drops()
    test_off_is_identity()
    test_build_never_supercritical_and_roundtrips()
    test_outside_field_is_zero()
    print("\nALL OVERPRESSURE UNIT TESTS PASSED")
