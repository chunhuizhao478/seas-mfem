"""
resolve_fault_stress.py
========================

Resolve a constant bulk stress tensor onto a fault surface of arbitrary
orientation, returning effective normal stress, shear stress, and apparent
friction coefficient.

Standard use case: SAF system dynamic rupture simulations where you set a
single bulk sigma^0 tensor (parameterized by SH_max, Sh_min, S_v and the
azimuth of SH_max) and let the geometry of each fault segment determine the
local tractions.

Conventions
-----------
- Coordinate frame: x = east, y = north, z = up.
- Compression is NEGATIVE for stress tensor components (continuum-mechanics
  convention). Pore pressure P_p is reported as a positive number; effective
  stress is sigma_eff = sigma_total + P_p * I (compression reduced).
- Principal stress magnitudes (sigma1, sigma3, sigma_v) are entered as
  POSITIVE numbers (magnitudes); the sign is applied internally.
- Geological azimuth: clockwise from north, 0 to 360 deg.
- Mathematical angle alpha: counterclockwise from east. Conversion is
  alpha = 90 - azimuth (with wrapping to [0, 360)).
- Fault strike: geological azimuth of the strike direction. For a right-
  lateral fault, the hanging-wall side is to the right when looking along
  strike; the unit normal is constructed accordingly.
- Vertical strike-slip faults (dip = 90 deg) are the default; other dips
  are supported.

References
----------
- Hickman & Zoback (2004) GRL 31, L15S12 -- SAFOD pilot hole stress
  orientation and magnitude; source for the reference parameter set below.
- Townend & Zoback (2004) GRL 31, L15S11 -- regional NNE-SSW SH_max.
- Aochi & Madariaga (2003) BSSA 93, 1249-1266 -- stress resolution onto
  non-planar faults, the methodology this code implements.

Author: written for Andy Zhao (USC / SCEC), May 2026.
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# Canonical default path for the H&Z σ⁰ dump file consumed by Phases 3/5/8
# of the on-fault stress preprocessing pipeline (see
# code_preprocess/data_projection_onfaultstress/PLAN_onfaultstress.md
# Phase 0).  Declared here so downstream code can import it without
# per-test configuration.
DEFAULT_DUMP_PATH = Path(__file__).resolve().parent / "hickman_zoback_sigma0.json"

# -----------------------------------------------------------------------------
# Angle utilities
# -----------------------------------------------------------------------------


def azimuth_to_math_angle(azimuth_deg: float) -> float:
    """Convert geological azimuth (clockwise from north) to mathematical
    angle alpha (counterclockwise from east). Result is in [0, 360)."""
    alpha = 90.0 - azimuth_deg
    return alpha % 360.0


def math_angle_to_azimuth(alpha_deg: float) -> float:
    """Inverse of azimuth_to_math_angle, result in [0, 360)."""
    az = 90.0 - alpha_deg
    return az % 360.0


# -----------------------------------------------------------------------------
# Build the bulk stress tensor
# -----------------------------------------------------------------------------


def build_bulk_stress_tensor(
    SHmax: float,
    Shmin: float,
    Sv: float,
    SHmax_azimuth_deg: float,
) -> np.ndarray:
    """
    Construct the bulk initial stress tensor sigma^0 (3x3, total stress,
    compression negative) in the global (East, North, Up) Cartesian frame.

    Parameters
    ----------
    SHmax, Shmin, Sv : float
        Magnitudes (POSITIVE values in MPa) of the maximum horizontal,
        minimum horizontal, and vertical principal stresses.
    SHmax_azimuth_deg : float
        Geological azimuth of the SH_max direction, clockwise from north, in
        degrees. For the SAFOD deep result (Psi = 69 deg from local SAF
        strike of N46W), this is approximately 023 (N23E).

    Returns
    -------
    sigma0 : (3, 3) ndarray
        Stress tensor in the global frame. Compression is negative.

    Notes
    -----
    The construction starts from a diagonal tensor in the principal frame
    (SH_max along x', Sh_min along y', S_v along z) and rotates about the
    vertical axis by alpha = 90 - azimuth so that the principal x'-axis
    points in the SH_max direction in the global (x = east, y = north) plane.
    """
    alpha = np.deg2rad(azimuth_to_math_angle(SHmax_azimuth_deg))
    c, s = np.cos(alpha), np.sin(alpha)

    # Principal-frame stress (compression negative)
    sigma_principal = np.diag([-SHmax, -Shmin, -Sv])

    # Rotation about vertical axis. Rz rotates the principal frame so that
    # its x'-axis ends up at angle alpha in the global frame.
    Rz = np.array(
        [
            [c, -s, 0.0],
            [s, c, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )

    sigma0 = Rz @ sigma_principal @ Rz.T
    return sigma0


# -----------------------------------------------------------------------------
# Build fault basis vectors
# -----------------------------------------------------------------------------


def fault_basis_vectors(
    strike_azimuth_deg: float,
    dip_deg: float = 90.0,
    rake_sense: str = "right-lateral",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Construct the orthonormal fault basis (strike-vector, dip-vector,
    normal-vector) in the global Cartesian frame.

    Parameters
    ----------
    strike_azimuth_deg : float
        Geological azimuth of the fault strike, clockwise from north.
        For the SAF near Parkfield: strike = N46W = azimuth 314.
    dip_deg : float, default 90 (vertical)
        Dip angle of the fault, from horizontal. 90 = vertical, 0 = flat.
    rake_sense : {"right-lateral", "left-lateral"}, default "right-lateral"
        Determines which side of the fault is the hanging wall, which sets
        the sign of the normal vector. For the SAF, this is "right-lateral".

    Returns
    -------
    s_hat, d_hat, n_hat : (3,) ndarray each
        Unit vectors along strike, up-dip, and fault-normal, respectively.
    """
    # Strike direction: horizontal unit vector at given azimuth
    az_rad = np.deg2rad(strike_azimuth_deg)
    s_hat = np.array([np.sin(az_rad), np.cos(az_rad), 0.0])

    # Fault normal: rotate strike 90 deg in the horizontal plane. We orient
    # the normal so that for a right-lateral fault under regional compression
    # at high angle to strike, the resolved shear along strike is POSITIVE
    # (i.e., consistent with right-lateral slip). For the SAF with SH_max at
    # high angle, this means the normal points to the LEFT of s_hat in map
    # view (rotate counterclockwise 90 deg).
    if rake_sense == "right-lateral":
        n_horizontal = np.array([-np.cos(az_rad), np.sin(az_rad), 0.0])
    elif rake_sense == "left-lateral":
        n_horizontal = np.array([np.cos(az_rad), -np.sin(az_rad), 0.0])
    else:
        raise ValueError(
            f"rake_sense must be right-lateral or left-lateral, got {rake_sense}"
        )

    # Apply dip: if dip < 90, the fault leans, so the normal tilts from
    # horizontal toward the vertical. Up-dip direction in the fault plane.
    dip_rad = np.deg2rad(dip_deg)
    n_hat = np.sin(dip_rad) * n_horizontal + np.cos(dip_rad) * np.array([0.0, 0.0, 1.0])
    n_hat = n_hat / np.linalg.norm(n_hat)

    # Up-dip vector: perpendicular to both s_hat and n_hat, pointing upward
    d_hat = np.cross(n_hat, s_hat)
    d_hat = d_hat / np.linalg.norm(d_hat)

    return s_hat, d_hat, n_hat


# -----------------------------------------------------------------------------
# Resolve the tractions
# -----------------------------------------------------------------------------


@dataclass
class ResolvedTraction:
    """Container for the resolved fault traction quantities."""

    sigma_n_total: float  # Total normal stress (compression negative), MPa
    sigma_n_eff: float  # Effective normal stress (after subtracting P_p), MPa
    tau_strike: float  # Shear along strike, MPa (+ for right-lateral on the fault)
    tau_dip: float  # Shear up-dip, MPa
    tau_magnitude: float  # Shear magnitude sqrt(tau_s^2 + tau_dip^2), MPa
    rake_deg: (
        float  # Rake angle in degrees (0 = pure strike-slip RL, 90 = pure reverse)
    )
    mu_apparent: float  # tau / |sigma_n_eff|

    def summary(self) -> str:
        return (
            f"  sigma_n_total = {self.sigma_n_total:+8.2f} MPa  (compression negative)\n"
            f"  sigma_n_eff   = {self.sigma_n_eff:+8.2f} MPa  (after pore pressure)\n"
            f"  tau_strike    = {self.tau_strike:+8.2f} MPa\n"
            f"  tau_dip       = {self.tau_dip:+8.2f} MPa\n"
            f"  |tau|         = {self.tau_magnitude:8.2f} MPa\n"
            f"  rake          = {self.rake_deg:+8.2f} deg\n"
            f"  mu_apparent   = {self.mu_apparent:8.3f}\n"
        )


def resolve_traction(
    sigma0: np.ndarray,
    s_hat: np.ndarray,
    d_hat: np.ndarray,
    n_hat: np.ndarray,
    P_p: float = 0.0,
) -> ResolvedTraction:
    """
    Resolve the bulk stress tensor onto a fault plane and return effective
    normal stress, shear components, and apparent friction.

    Parameters
    ----------
    sigma0 : (3, 3) ndarray
        Bulk stress tensor (compression negative), MPa.
    s_hat, d_hat, n_hat : (3,) ndarray
        Orthonormal fault basis (strike, up-dip, normal).
    P_p : float, default 0
        Pore pressure (POSITIVE value in MPa).

    Returns
    -------
    ResolvedTraction
        Resolved quantities with the conventions documented in the dataclass.
    """
    # Traction vector on the fault plane: t = sigma . n
    t_vec = sigma0 @ n_hat

    # Decompose
    sigma_n_total = float(n_hat @ t_vec)
    tau_s = float(s_hat @ t_vec)
    tau_d = float(d_hat @ t_vec)

    sigma_n_eff = sigma_n_total + P_p  # compression reduced

    tau_mag = float(np.sqrt(tau_s**2 + tau_d**2))
    rake = float(np.rad2deg(np.arctan2(tau_d, tau_s)))

    # Apparent friction: |shear| / |effective normal| (effective normal is
    # negative under compression, so we take its absolute value)
    if abs(sigma_n_eff) < 1e-9:
        mu_app = float("nan")
    else:
        mu_app = tau_mag / abs(sigma_n_eff)

    return ResolvedTraction(
        sigma_n_total=sigma_n_total,
        sigma_n_eff=sigma_n_eff,
        tau_strike=tau_s,
        tau_dip=tau_d,
        tau_magnitude=tau_mag,
        rake_deg=rake,
        mu_apparent=mu_app,
    )


# -----------------------------------------------------------------------------
# Top-level convenience wrapper
# -----------------------------------------------------------------------------


def compute_fault_stress(
    SHmax: float,
    Shmin: float,
    Sv: float,
    P_p: float,
    SHmax_azimuth_deg: float,
    fault_strike_azimuth_deg: float,
    fault_dip_deg: float = 90.0,
    rake_sense: str = "right-lateral",
) -> tuple[np.ndarray, ResolvedTraction]:
    """
    Single-call interface: takes geological inputs and returns the bulk
    stress tensor plus the resolved fault traction.

    Parameters
    ----------
    SHmax, Shmin, Sv : float
        Principal stress magnitudes (POSITIVE values, MPa). For SAF strike-slip,
        the ordering should be Shmin < Sv < SHmax. If your inputs give a
        different ordering, the faulting regime is reverse (Sv < Shmin < SHmax)
        or normal (SHmax < Sv with Shmin smallest).
    P_p : float
        Pore pressure (POSITIVE value, MPa). Use 0 if working in effective stress.
    SHmax_azimuth_deg : float
        Geological azimuth of SH_max (clockwise from north). For the SAFOD
        deep result, this is approximately 023 (N23E).
    fault_strike_azimuth_deg : float
        Geological azimuth of the fault strike. For the SAF near Parkfield,
        this is 314 (N46W).
    fault_dip_deg : float, default 90 (vertical)
    rake_sense : str, default "right-lateral"

    Returns
    -------
    sigma0 : (3, 3) ndarray
        Bulk stress tensor (MPa, compression negative).
    traction : ResolvedTraction
        Resolved normal stress, shear, and apparent friction.

    Examples
    --------
    SAFOD deep result (Hickman & Zoback 2004) at 1671 m depth:

    >>> sigma0, t = compute_fault_stress(
    ...     SHmax=113.0, Shmin=49.0, Sv=45.0, P_p=16.0,
    ...     SHmax_azimuth_deg=23.0,
    ...     fault_strike_azimuth_deg=314.0,
    ... )
    >>> print(t.summary())
    """
    # Sanity check on faulting regime (informational only; the SAF at SAFOD
    # depth is famously transitional between strike-slip and reverse, so this
    # often reports "reverse" for literal pilot-hole magnitudes -- that's a
    # known feature of the data, not a code error).
    ordering = sorted(
        [("SHmax", SHmax), ("Shmin", Shmin), ("Sv", Sv)], key=lambda x: -x[1]
    )
    sigma1, sigma2, sigma3 = ordering[0][0], ordering[1][0], ordering[2][0]
    if (sigma1, sigma2, sigma3) == ("SHmax", "Sv", "Shmin"):
        regime = "strike-slip"
    elif (sigma1, sigma2, sigma3) == ("SHmax", "Shmin", "Sv"):
        regime = "reverse (or transitional SS-RF if Sv ~ Shmin)"
    elif (sigma1, sigma2, sigma3) == ("Sv", "SHmax", "Shmin"):
        regime = "normal"
    else:
        regime = f"unusual: sigma1={sigma1}, sigma2={sigma2}, sigma3={sigma3}"

    sigma0 = build_bulk_stress_tensor(SHmax, Shmin, Sv, SHmax_azimuth_deg)
    s_hat, d_hat, n_hat = fault_basis_vectors(
        fault_strike_azimuth_deg, fault_dip_deg, rake_sense
    )
    traction = resolve_traction(sigma0, s_hat, d_hat, n_hat, P_p)

    # Attach the regime as metadata via a print on demand (not stored in dataclass)
    compute_fault_stress.last_regime = regime
    return sigma0, traction


# -----------------------------------------------------------------------------
# Demonstration: reproduce the Hickman-Zoback (2004) SAFOD result
# -----------------------------------------------------------------------------


def demo_safod():
    """Reproduce the SAFOD pilot hole result at 1671 m depth.

    Expected outputs (from Hickman & Zoback 2004, paragraphs [11] and [14]):
    - sigma_n_eff approx -88 to -89 MPa (effective normal stress on a vertical
      plane parallel to the SAF)
    - mu_apparent approx 0.24 (deep value; lower end of the depth range)

    Note: The 0.24 value in the paper is for the deepest interval (2050-2200 m,
    where Psi = 69 deg). At 1671 m the Psi is intermediate (~37-55 deg from
    the paper's depth table), so the literal 1671-m magnitudes combined with
    the deep Psi = 69 deg give the deep mu_app of about 0.24.
    """
    print("=" * 70)
    print(" SAFOD pilot hole reference case (Hickman & Zoback 2004)")
    print(" Magnitudes at 1671 m; orientation Psi = 69 deg (deep result)")
    print("=" * 70)

    sigma0, traction = compute_fault_stress(
        SHmax=113.0,  # MPa, max horizontal principal stress
        Shmin=49.0,  # MPa, min horizontal principal stress
        Sv=45.0,  # MPa, vertical stress at 1671 m (rho * g * z)
        P_p=16.0,  # MPa, hydrostatic pore pressure at 1671 m
        SHmax_azimuth_deg=23.0,  # SH_max points N23E
        fault_strike_azimuth_deg=314.0,  # SAF strikes N46W = azimuth 314
        fault_dip_deg=90.0,
        rake_sense="right-lateral",
    )

    print(f"\n Faulting regime: {compute_fault_stress.last_regime}")
    print("\n Bulk stress tensor sigma^0 (MPa, compression negative):")
    np.set_printoptions(precision=2, suppress=True, sign="+")
    print(sigma0)

    print("\n Resolved fault traction:")
    print(traction.summary())

    print(" Expected (Hickman & Zoback 2004, deep SAF):")
    print("   mu_apparent ~ 0.24  (paragraph [14] of the paper)")
    print()


# -----------------------------------------------------------------------------
# Bonus: sweep over fault strike to visualize geometry effects
# -----------------------------------------------------------------------------


def sweep_strike(
    SHmax: float = 113.0,
    Shmin: float = 49.0,
    Sv: float = 80.0,  # bumped above Shmin to enforce strike-slip
    P_p: float = 16.0,
    SHmax_azimuth_deg: float = 23.0,
    n_points: int = 37,
):
    """Sweep fault strike from 0 to 360 deg and print the resolved traction
    on a vertical right-lateral fault at each orientation.

    Useful for visualizing how a kinked or curved fault picks up varying
    sigma_n and tau from a constant bulk sigma^0.
    """
    print("=" * 70)
    print(" Sweep of fault strike under constant bulk stress")
    print(f" SH_max = {SHmax} MPa at azimuth {SHmax_azimuth_deg} deg")
    print(f" Sh_min = {Shmin} MPa, S_v = {Sv} MPa, P_p = {P_p} MPa")
    print("=" * 70)
    print(f"{'strike(az)':>10s}  {'sigma_n_eff':>12s}  {'tau':>10s}  {'mu_app':>8s}")

    for az in np.linspace(0, 360, n_points):
        _, t = compute_fault_stress(
            SHmax,
            Shmin,
            Sv,
            P_p,
            SHmax_azimuth_deg,
            float(az),
            fault_dip_deg=90.0,
            rake_sense="right-lateral",
        )
        print(
            f"{az:10.1f}  {t.sigma_n_eff:+12.2f}  {t.tau_magnitude:10.2f}  {t.mu_apparent:8.3f}"
        )


# -----------------------------------------------------------------------------
# Dump entry point (Phase 0 of the on-fault stress preprocessing pipeline)
# -----------------------------------------------------------------------------


def dump_safod_sigma0(
    out_path: Path,
    *,
    SHmax: float = 113.0,
    Shmin: float = 49.0,
    Sv: float = 45.0,
    P_p: float = 16.0,
    SHmax_azimuth_deg: float = 23.0,
    fault_strike_azimuth_deg: float = 314.0,
    fault_dip_deg: float = 90.0,
    rake_sense: str = "right-lateral",
) -> None:
    """Dump the H&Z σ⁰ tensor and input parameters to a JSON sidecar.

    The defaults reproduce ``demo_safod()`` exactly (SAFOD pilot hole at
    1671 m depth, Hickman & Zoback 2004).  The output schema is
    ``hickman_zoback_sigma0_v1`` and is consumed by Phase 3 of the
    on-fault stress preprocessing pipeline via ``--hz-dump-file``.

    Parameters
    ----------
    out_path : Path
        Destination JSON file.  Parent directory is created if missing.
    SHmax, Shmin, Sv : float
        Principal-stress magnitudes (positive values, MPa).
    P_p : float
        Pore pressure (positive value, MPa).
    SHmax_azimuth_deg : float
        Geological azimuth of SH_max (clockwise from north).
    fault_strike_azimuth_deg : float
        Geological azimuth of the fault strike (default SAF N46W = 314°).
    fault_dip_deg : float
        Fault dip (default vertical = 90°).
    rake_sense : str
        ``"right-lateral"`` or ``"left-lateral"``.

    Notes
    -----
    σ⁰ is written in the H&Z **continuum-mechanics** convention
    (compression NEGATIVE).  The split-site sign-flip contract used by
    the downstream pipeline (Phase 3 / 4 / 5 of the on-fault stress
    plan) flips this to compression-positive at well-defined
    boundaries; the dump file itself preserves the H&Z native sign so
    a reader of this script sees the same numbers it prints in
    ``demo_safod()``.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    sigma0, _traction = compute_fault_stress(
        SHmax=SHmax,
        Shmin=Shmin,
        Sv=Sv,
        P_p=P_p,
        SHmax_azimuth_deg=SHmax_azimuth_deg,
        fault_strike_azimuth_deg=fault_strike_azimuth_deg,
        fault_dip_deg=fault_dip_deg,
        rake_sense=rake_sense,
    )

    payload = {
        "schema": "hickman_zoback_sigma0_v1",
        "params": {
            "SHmax_MPa": float(SHmax),
            "Shmin_MPa": float(Shmin),
            "Sv_MPa": float(Sv),
            "P_p_MPa": float(P_p),
            "SHmax_azimuth_deg": float(SHmax_azimuth_deg),
            "fault_strike_azimuth_deg": float(fault_strike_azimuth_deg),
            "fault_dip_deg": float(fault_dip_deg),
            "rake_sense": str(rake_sense),
        },
        "sigma0_MPa": [[float(v) for v in row] for row in sigma0],
        "convention": (
            "compression negative; continuum-mechanics; MPa; "
            "frame (east, north, up)"
        ),
        "regime": str(getattr(compute_fault_stress, "last_regime", "unknown")),
    }

    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "--dump":
        if len(sys.argv) < 3:
            print(
                "error: --dump requires a path argument\n"
                "usage: hickman_and_zoback_regional_stress_projection.py "
                "--dump <out_path>",
                file=sys.stderr,
            )
            sys.exit(1)
        dump_safod_sigma0(Path(sys.argv[2]))
        print(f"Wrote {sys.argv[2]}")
    else:
        demo_safod()
        print()
        sweep_strike()
