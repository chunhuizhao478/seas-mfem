#!/usr/bin/env python3
# Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory. All rights
# reserved.  See files LICENSE and NOTICE for details.  LLNL-CODE-806117.
#
# Pytest gate-metric tests for tpv31/visualize_results.py (Phase 7, BUG-11).
# These run WITHOUT any mesh / MFEM / Frontera run — they synthesize tiny
# .dat / .txt fixtures in a tmp dir and/or call the gate metric functions
# directly.  No matplotlib import is required (the gate is headless).
#
# Coverage:
#   NEW-7.1  offset time grids => exercises np.interp; one >band channel
#            => exit != 0; one within-band => exit 0.
#   NEW-7.4  zero-crossing reference (pre-nucleation zeros + slip-rate pulse)
#            => peak-normalized metric does NOT spuriously fail.
#   NEW-7.5  truncated MFEM trace (covers far less than the reference span)
#            => non-zero exit with an explicit insufficient-coverage failure.
#
# Run:  pytest tpv31/test_visualize_results_tol.py -v
#       (or `python3 -m pytest`)

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np

# Make `visualize_results` importable regardless of pytest's rootdir.
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import visualize_results as vr  # noqa: E402

SCRIPT = HERE / "visualize_results.py"

# Reference 8-column field-name line (begins with a letter => parser skips it).
_FIELD_LINE = ("t h-slip h-slip-rate h-shear-stress "
               "v-slip v-slip-rate v-shear-stress n-stress")


# ---------------------------------------------------------------------------
# Fixture writers (match the formats parsed by load_mfem_file /
# load_reference_file: MFEM 9-col Pa; reference 8-col MPa, n-stress
# compression-NEGATIVE in the file).
# ---------------------------------------------------------------------------

def _write_mfem_dat(path, t, v_strike, slip_strike, tau_strike_mpa,
                    sigma_n_mpa, v_dip=None, slip_dip=None,
                    tau_dip_mpa=None, mu_eff=None):
    """Write a 9-column MFEM station .dat (stresses written in Pa)."""
    n = len(t)
    z = np.zeros(n)
    v_dip = z if v_dip is None else v_dip
    slip_dip = z if slip_dip is None else slip_dip
    tau_dip_mpa = z if tau_dip_mpa is None else tau_dip_mpa
    mu_eff = np.full(n, 0.5) if mu_eff is None else mu_eff
    pa = 1.0e6  # MPa -> Pa for the file
    with open(path, "w") as f:
        f.write("# MFEM TPV31 station (synthetic test fixture)\n")
        f.write("# t h-slip h-slip-rate h-shear v-slip v-slip-rate "
                "v-shear n-stress mu_eff\n")
        for i in range(n):
            f.write(f"{t[i]:.8e} {slip_strike[i]:.8e} {v_strike[i]:.8e} "
                    f"{tau_strike_mpa[i] * pa:.8e} {slip_dip[i]:.8e} "
                    f"{v_dip[i]:.8e} {tau_dip_mpa[i] * pa:.8e} "
                    f"{sigma_n_mpa[i] * pa:.8e} {mu_eff[i]:.8e}\n")


def _write_ref_txt(path, t, v_strike, slip_strike, tau_strike_mpa,
                   sigma_n_mpa_pos, v_dip=None, slip_dip=None,
                   tau_dip_mpa=None):
    """Write an 8-column SCEC reference .txt.

    `sigma_n_mpa_pos` is the desired compression-POSITIVE value (as the
    reader returns it); it is written NEGATED into the file to mirror the
    SCEC compression-negative convention.
    """
    n = len(t)
    z = np.zeros(n)
    v_dip = z if v_dip is None else v_dip
    slip_dip = z if slip_dip is None else slip_dip
    tau_dip_mpa = z if tau_dip_mpa is None else tau_dip_mpa
    with open(path, "w") as f:
        f.write("# problem=TPV31 (synthetic test fixture)\n")
        f.write("# author=test\n")
        f.write("#\n")
        f.write("# The line below lists the names of the data fields:\n")
        f.write(_FIELD_LINE + "\n")
        f.write("#\n")
        for i in range(n):
            f.write(f"  {t[i]:.10e}   {slip_strike[i]:.7e}   "
                    f"{v_strike[i]:.7e}   {tau_strike_mpa[i]:.7e}   "
                    f"{slip_dip[i]:.7e}   {v_dip[i]:.7e}   "
                    f"{tau_dip_mpa[i]:.7e}   {-sigma_n_mpa_pos[i]:.7e}\n")


# ---------------------------------------------------------------------------
# NEW-7.1: offset time grids => np.interp exercised; >band => fail, in-band
# => pass.
# ---------------------------------------------------------------------------

def test_new_7_1_offset_grids_within_band_passes():
    """Offset (different) time grids over a shared span; MFEM ~ reference
    within band => gate_channel reports small error and run_tolerance_gate
    exits 0."""
    ref_t = np.linspace(0.0, 10.0, 101)            # fixed cadence
    mfem_t = np.linspace(0.0, 10.0, 137)           # different (RK45-like) grid
    # Smooth signal so interpolation between offset grids is faithful.
    ref_y = 2.0 + np.sin(ref_t)
    mfem_y = 2.0 + np.sin(mfem_t)                   # same underlying function

    peak_rel, rms_rel, covered = vr.gate_channel(
        mfem_t, mfem_y, ref_t, ref_y)
    assert covered
    # The two sample the SAME smooth function; interpolation residual is tiny.
    assert peak_rel < 0.05
    assert rms_rel < 0.05

    # Full pair through run_tolerance_gate => exit 0.
    sigma = np.full_like(ref_t, 60.0)
    mfem = {
        "time_s": mfem_t, "V_strike": mfem_y,
        "slip_strike": np.cumsum(mfem_y) * 0.0 + 1.0,
        "tau_strike": np.full_like(mfem_t, 35.0),
        "slip_dip": np.zeros_like(mfem_t), "V_dip": np.zeros_like(mfem_t),
        "tau_dip": np.zeros_like(mfem_t),
        "sigma_n": np.full_like(mfem_t, 60.0),
        "mu_eff": np.full_like(mfem_t, 0.5),
    }
    ref = {
        "time_s": ref_t, "V_strike": ref_y,
        "slip_strike": np.full_like(ref_t, 1.0),
        "tau_strike": np.full_like(ref_t, 35.0),
        "slip_dip": np.zeros_like(ref_t), "V_dip": np.zeros_like(ref_t),
        "tau_dip": np.zeros_like(ref_t),
        "sigma_n": sigma,
        "mu_eff": np.full(ref_t.shape, np.nan),
    }
    rc = vr.run_tolerance_gate(
        [("st", mfem, ref)], tol_peak=0.10, tol_rms=0.10)
    assert rc == 0


def test_new_7_1_offset_grids_out_of_band_fails():
    """Same offset grids, but one channel (V_strike) is scaled well beyond
    the band => non-zero exit."""
    ref_t = np.linspace(0.0, 10.0, 101)
    mfem_t = np.linspace(0.0, 10.0, 137)
    ref_y = 2.0 + np.sin(ref_t)
    mfem_y = 1.5 * (2.0 + np.sin(mfem_t))          # 50% larger amplitude

    peak_rel, rms_rel, covered = vr.gate_channel(
        mfem_t, mfem_y, ref_t, ref_y)
    assert covered
    assert peak_rel > 0.10                         # clearly out of a 10% band

    mfem = {
        "time_s": mfem_t, "V_strike": mfem_y,
        "slip_strike": np.full_like(mfem_t, 1.0),
        "tau_strike": np.full_like(mfem_t, 35.0),
        "slip_dip": np.zeros_like(mfem_t), "V_dip": np.zeros_like(mfem_t),
        "tau_dip": np.zeros_like(mfem_t),
        "sigma_n": np.full_like(mfem_t, 60.0),
        "mu_eff": np.full_like(mfem_t, 0.5),
    }
    ref = {
        "time_s": ref_t, "V_strike": ref_y,
        "slip_strike": np.full_like(ref_t, 1.0),
        "tau_strike": np.full_like(ref_t, 35.0),
        "slip_dip": np.zeros_like(ref_t), "V_dip": np.zeros_like(ref_t),
        "tau_dip": np.zeros_like(ref_t),
        "sigma_n": np.full_like(ref_t, 60.0),
        "mu_eff": np.full(ref_t.shape, np.nan),
    }
    rc = vr.run_tolerance_gate(
        [("st", mfem, ref)], tol_peak=0.10, tol_rms=0.10)
    assert rc == 1


def test_new_7_1_end_to_end_files_and_exit_code(tmp_path):
    """Round-trip through the readers + main() with synthesized files on
    offset grids: within-band run exits 0, out-of-band run exits 1."""
    bench = tmp_path / "benchmark_data"
    seisol = bench / "scec_seisol"
    eqdyna = bench / "scec_eqdyna"
    seisol.mkdir(parents=True)
    eqdyna.mkdir(parents=True)

    run_ok = tmp_path / "run_ok"
    run_bad = tmp_path / "run_bad"
    run_ok.mkdir()
    run_bad.mkdir()

    # Use station 1 only (faultst000dp000, ref label x2_0_x3_0).
    mfem_name, ref_label, _s, _d = vr.SCEC_STATIONS[0]

    ref_t = np.linspace(0.0, 5.0, 51)
    mfem_t = np.linspace(0.0, 5.0, 73)             # offset grid
    ref_v = 1.0 + np.sin(ref_t)
    sigma = np.full_like(ref_t, 60.0)
    slip = np.full_like(ref_t, 0.3)
    tau = np.full_like(ref_t, 35.0)

    _write_ref_txt(
        seisol / f"tpv31_seisol_{ref_label}.txt",
        ref_t, ref_v, slip, tau, sigma)

    # OK run: same smooth function on the offset grid.
    _write_mfem_dat(
        run_ok / f"tpv31_station_{mfem_name}.dat",
        mfem_t, 1.0 + np.sin(mfem_t), np.full_like(mfem_t, 0.3),
        np.full_like(mfem_t, 35.0), np.full_like(mfem_t, 60.0))
    # Bad run: V_strike doubled => far out of band.
    _write_mfem_dat(
        run_bad / f"tpv31_station_{mfem_name}.dat",
        mfem_t, 2.0 * (1.0 + np.sin(mfem_t)), np.full_like(mfem_t, 0.3),
        np.full_like(mfem_t, 35.0), np.full_like(mfem_t, 60.0))

    rc_ok = _run_main(
        ["--mfem", str(run_ok), "--seisol",
         "--benchmark-dir", str(bench), "--stations", "1",
         "--tol-peak", "0.10", "--tol-rms", "0.10"])
    assert rc_ok == 0

    rc_bad = _run_main(
        ["--mfem", str(run_bad), "--seisol",
         "--benchmark-dir", str(bench), "--stations", "1",
         "--tol-peak", "0.10", "--tol-rms", "0.10"])
    assert rc_bad == 1


# ---------------------------------------------------------------------------
# NEW-7.4: zero-crossing reference => peak-normalized metric must not
# spuriously fail.
# ---------------------------------------------------------------------------

def test_new_7_4_zero_crossing_no_spurious_fail():
    """A reference that is exactly zero before nucleation, then a slip-rate
    pulse.  A per-sample |mfem-ref|/|ref| metric would blow up over the
    pre-nucleation zeros; the peak-normalized metric must stay small."""
    ref_t = np.linspace(0.0, 10.0, 201)
    mfem_t = np.linspace(0.0, 10.0, 263)           # offset grid

    def pulse(t):
        y = np.zeros_like(t)
        on = (t >= 4.0) & (t <= 6.0)
        # Smooth bump, zero outside [4,6].
        y[on] = np.sin((t[on] - 4.0) * np.pi / 2.0) ** 2 * 3.0
        return y

    ref_y = pulse(ref_t)
    mfem_y = pulse(mfem_t) * 1.02                   # 2% amplitude difference

    # Sanity: the reference genuinely contains exact zeros (the case a
    # per-sample relative metric would choke on).
    assert np.any(ref_y == 0.0)

    peak_rel, rms_rel, covered = vr.gate_channel(
        mfem_t, mfem_y, ref_t, ref_y)
    assert covered
    # Peak-normalized: ~2% diff, NOT a divide-by-zero blowup.
    assert np.isfinite(peak_rel) and np.isfinite(rms_rel)
    assert peak_rel < 0.10
    assert rms_rel < 0.10

    mfem = {
        "time_s": mfem_t, "V_strike": mfem_y,
        "slip_strike": np.cumsum(mfem_y) * (mfem_t[1] - mfem_t[0]),
        "tau_strike": np.full_like(mfem_t, 35.0),
        "slip_dip": np.zeros_like(mfem_t), "V_dip": np.zeros_like(mfem_t),
        "tau_dip": np.zeros_like(mfem_t),
        "sigma_n": np.full_like(mfem_t, 60.0),
        "mu_eff": np.full_like(mfem_t, 0.5),
    }
    ref = {
        "time_s": ref_t, "V_strike": ref_y,
        "slip_strike": np.cumsum(ref_y) * (ref_t[1] - ref_t[0]),
        "tau_strike": np.full_like(ref_t, 35.0),
        "slip_dip": np.zeros_like(ref_t), "V_dip": np.zeros_like(ref_t),
        "tau_dip": np.zeros_like(ref_t),
        "sigma_n": np.full_like(ref_t, 60.0),
        "mu_eff": np.full(ref_t.shape, np.nan),
    }
    rc = vr.run_tolerance_gate(
        [("st", mfem, ref)], tol_peak=0.10, tol_rms=0.10)
    assert rc == 0


# ---------------------------------------------------------------------------
# NEW-7.5 (BUG-23): truncated MFEM trace => insufficient-coverage failure.
# ---------------------------------------------------------------------------

def test_new_7_5_truncated_run_fails_coverage():
    """MFEM run that covers only the first ~40% of the reference span must
    FAIL the coverage guard (not silently pass over the tiny early window),
    EVEN when the overlapping samples agree perfectly."""
    ref_t = np.linspace(0.0, 10.0, 101)
    ref_y = 2.0 + np.sin(ref_t)

    # Truncated: MFEM ends at t=4.0 (40% of the 0..10 span).
    mfem_t = np.linspace(0.0, 4.0, 41)
    mfem_y = 2.0 + np.sin(mfem_t)                   # perfect agreement here

    peak_rel, rms_rel, covered = vr.gate_channel(
        mfem_t, mfem_y, ref_t, ref_y)
    # Over the overlap the data agree, but coverage (40%) < 90% => not covered.
    assert not covered

    mfem = {
        "time_s": mfem_t, "V_strike": mfem_y,
        "slip_strike": np.full_like(mfem_t, 1.0),
        "tau_strike": np.full_like(mfem_t, 35.0),
        "slip_dip": np.zeros_like(mfem_t), "V_dip": np.zeros_like(mfem_t),
        "tau_dip": np.zeros_like(mfem_t),
        "sigma_n": np.full_like(mfem_t, 60.0),
        "mu_eff": np.full_like(mfem_t, 0.5),
    }
    ref = {
        "time_s": ref_t, "V_strike": ref_y,
        "slip_strike": np.full_like(ref_t, 1.0),
        "tau_strike": np.full_like(ref_t, 35.0),
        "slip_dip": np.zeros_like(ref_t), "V_dip": np.zeros_like(ref_t),
        "tau_dip": np.zeros_like(ref_t),
        "sigma_n": np.full_like(ref_t, 60.0),
        "mu_eff": np.full(ref_t.shape, np.nan),
    }

    msgs = []
    rc = vr.run_tolerance_gate(
        [("st", mfem, ref)], tol_peak=0.10, tol_rms=0.10,
        printer=msgs.append)
    assert rc == 1
    joined = "\n".join(msgs)
    assert "insufficient coverage" in joined.lower()


def test_new_7_5_end_to_end_truncated_exit_nonzero(tmp_path):
    """Same truncation but through synthesized files + main()."""
    bench = tmp_path / "benchmark_data"
    seisol = bench / "scec_seisol"
    seisol.mkdir(parents=True)
    run = tmp_path / "run_trunc"
    run.mkdir()

    mfem_name, ref_label, _s, _d = vr.SCEC_STATIONS[0]
    ref_t = np.linspace(0.0, 10.0, 101)
    ref_v = 2.0 + np.sin(ref_t)
    _write_ref_txt(
        seisol / f"tpv31_seisol_{ref_label}.txt",
        ref_t, ref_v, np.full_like(ref_t, 1.0),
        np.full_like(ref_t, 35.0), np.full_like(ref_t, 60.0))

    mfem_t = np.linspace(0.0, 4.0, 41)             # only 40% coverage
    _write_mfem_dat(
        run / f"tpv31_station_{mfem_name}.dat",
        mfem_t, 2.0 + np.sin(mfem_t), np.full_like(mfem_t, 1.0),
        np.full_like(mfem_t, 35.0), np.full_like(mfem_t, 60.0))

    rc = _run_main(
        ["--mfem", str(run), "--seisol",
         "--benchmark-dir", str(bench), "--stations", "1",
         "--tol-peak", "0.10"])
    assert rc == 1


# ---------------------------------------------------------------------------
# Supporting checks: mu_eff is never gated; dip channels gated only on dip
# motion.
# ---------------------------------------------------------------------------

def test_mu_eff_never_gated():
    ref = {
        "V_dip": np.zeros(10),
    }
    chans = vr.gated_channels_for_station(ref)
    assert "mu_eff" not in chans
    assert set(chans) == set(vr.STRIKE_GATE_CHANNELS)


def test_dip_channels_added_only_with_motion():
    ref_no_dip = {"V_dip": np.full(10, 1.0e-9)}
    ref_dip = {"V_dip": np.full(10, 0.5)}
    assert "V_dip" not in vr.gated_channels_for_station(ref_no_dip)
    assert "V_dip" in vr.gated_channels_for_station(ref_dip)


# ---------------------------------------------------------------------------
# Helpers to invoke main() with a synthetic argv (no subprocess needed).
# ---------------------------------------------------------------------------

def _run_main(argv):
    """Call visualize_results.main() with sys.argv patched to `argv`.

    main() inspects sys.argv directly to preserve typed-source ordering, so
    we patch the real sys.argv around the call.
    """
    old = sys.argv
    sys.argv = ["visualize_results.py"] + list(argv)
    try:
        return vr.main()
    finally:
        sys.argv = old


def test_subprocess_help_lists_new_flags():
    """Smoke: the CLI advertises the new flags (also confirms headless
    import: the module is imported with matplotlib guarded)."""
    out = subprocess.run(
        [sys.executable, str(SCRIPT), "--help"],
        capture_output=True, text=True, check=True)
    assert "--tol-peak" in out.stdout
    assert "--tol-rms" in out.stdout


def test_f1_all_gated_channels_nan_fails_not_silent_pass():
    """F-1 (Rev I-10): a station that HAS data but whose every gated reference
    channel is all-NaN must FAIL (exit 1), not silently pass.  Symmetric to the
    BUG-18 mu_eff all-NaN class: a gate over zero comparable channels is not
    agreement.
    """
    t = np.linspace(0.0, 10.0, 50)
    nan = np.full_like(t, np.nan)
    mfem = {
        "time_s": t, "V_strike": 1.0 + np.sin(t),
        "slip_strike": np.full_like(t, 1.0), "tau_strike": np.full_like(t, 35.0),
        "slip_dip": np.zeros_like(t), "V_dip": np.zeros_like(t),
        "tau_dip": np.zeros_like(t), "sigma_n": np.full_like(t, 60.0),
        "mu_eff": np.full_like(t, 0.5),
    }
    # Reference present (valid time) but ALL gated channels are NaN.
    ref = {
        "time_s": t, "V_strike": nan.copy(), "slip_strike": nan.copy(),
        "tau_strike": nan.copy(), "slip_dip": np.zeros_like(t),
        "V_dip": np.zeros_like(t), "tau_dip": np.zeros_like(t),
        "sigma_n": nan.copy(), "mu_eff": nan.copy(),
    }
    rc = vr.run_tolerance_gate(
        [("st", mfem, ref)], tol_peak=0.10, tol_rms=0.10)
    assert rc == 1
