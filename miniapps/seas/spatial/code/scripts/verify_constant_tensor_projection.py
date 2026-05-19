#!/usr/bin/env python3
"""verify_constant_tensor_projection.py — Phase 3b cross-check.

Standalone Python reference for the per-DOF projection of a constant
background Cauchy tensor onto a fault.  Mirrors the C++ logic in
mfem::seas::spatial::ConstantTensorStressSource +
FaultGeometry::ComputeParams<StressSource>.

Used by the Phase 3b acceptance criterion that compares the C++ output
at 64 randomly sampled fault DOFs against this Python reference (L∞
≤ 1e-8 relative).  That comparison runs as an offline step in Phase 4,
not in the standalone unit test; this script is the reference
implementation.

Usage (standalone smoke check):

    python3 verify_constant_tensor_projection.py

Run with --help for the full CLI.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from typing import Tuple

import numpy as np


@dataclass(frozen=True)
class CauchyTensor:
    """Symmetric Cauchy tensor (compression POSITIVE) in EAST-NORTH-UP."""

    sxx: float
    syy: float
    szz: float
    sxy: float
    syz: float
    sxz: float

    def matrix(self) -> np.ndarray:
        return np.array(
            [
                [self.sxx, self.sxy, self.sxz],
                [self.sxy, self.syy, self.syz],
                [self.sxz, self.syz, self.szz],
            ],
            dtype=np.float64,
        )


def resolve_traction(
    sigma: CauchyTensor,
    n: np.ndarray,
    t1: np.ndarray,
    t2: np.ndarray,
    P_p: float = 0.0,
) -> Tuple[float, float, float]:
    """Per-DOF projection: returns (sigma_n_eff, tau_dip, tau_strike).

    Sign convention (CLAUDE.md project-wide; matches the native
    TPV102/104/205 drivers): positive tau{1,2} represents driving
    stress in the +dip / +strike direction.  The raw Cauchy projection
    T = σ·n returns the traction the n-side block exerts on the
    opposite block — the Newton's-3rd-law mirror of the driving stress
    — so flip the sign on tau1 / tau2 to match the native convention.
    sigma_n_total = n·S·n is unchanged (sign-invariant under n → −n).
    See R-001 in tpv102_tpv104_review.md.
    """
    S = sigma.matrix()
    Sn = S @ n
    sigma_n_total = float(n @ Sn)
    tau1 = -float(t1 @ Sn)
    tau2 = -float(t2 @ Sn)
    sigma_n_eff = sigma_n_total - P_p
    return sigma_n_eff, tau1, tau2


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sxx", type=float, default=80e6)
    parser.add_argument("--syy", type=float, default=80e6)
    parser.add_argument("--szz", type=float, default=80e6)
    parser.add_argument("--sxy", type=float, default=20e6)
    parser.add_argument("--syz", type=float, default=0.0)
    parser.add_argument("--sxz", type=float, default=0.0)
    parser.add_argument("--P_p_pa", type=float, default=0.0)
    args = parser.parse_args()

    sigma = CauchyTensor(args.sxx, args.syy, args.szz,
                         args.sxy, args.syz, args.sxz)

    # Three test normals (planar-fault sanity).
    normals = [
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
        np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0),
    ]
    # Up = (0, 0, 1); choose dip = n × up and strike = up × dip per
    # FaultBasis (Tandem) convention, then normalize.
    up = np.array([0.0, 0.0, 1.0])
    for n in normals:
        # If n is parallel to up, fall back to a 2-D in-plane convention.
        if abs(np.dot(n, up)) > 0.9999:
            t2 = np.array([1.0, 0.0, 0.0])
        else:
            t2 = np.cross(up, n)
            t2 /= np.linalg.norm(t2)
        t1 = np.cross(n, t2)
        t1 /= np.linalg.norm(t1)

        sn, td, ts = resolve_traction(sigma, n, t1, t2, P_p=args.P_p_pa)
        print(
            f"n = {n.tolist()} -> "
            f"sigma_n_eff = {sn:.3e} Pa, "
            f"tau_dip = {td:.3e} Pa, "
            f"tau_strike = {ts:.3e} Pa"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
