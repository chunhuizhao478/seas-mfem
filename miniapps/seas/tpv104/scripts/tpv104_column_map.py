"""tpv104_column_map.py — coordinate / sign-convention translation between
MFEM's TPV104 probe output (§5.1 format) and the reference FVW runtime's
TPV104 benchmark trace.

The two convention differences that matter (plan §3.2 + §3.11):
    1. Tangent-component ordering:
        MFEM  (BP5 / Tandem)       : tangent1 = dip,    tangent2 = strike
        Reference FVW              : tangent1 = strike, tangent2 = dip
       So MFEM ``slip1`` (dip) maps to reference ``slip2``, and
       MFEM ``slip2`` (strike) maps to reference ``slip1``.  The
       probe-diff tool swaps these channels before comparing.
    2. Normal-stress sign:
        MFEM               : sigma_n > 0 = compression (geology)
        Reference FVW      : sigma_n < 0 = compression (elasticity)
       The diff tool flips the sign of the reference's ``normalStress``
       channel before comparing against MFEM's ``sigma_n_*``.

This module is a thin, stateless collection of pure-Python helpers so
it runs under the project's default ``pythonenv`` conda env without
additional dependencies (only ``numpy`` is required).
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Channel-name canonicalisation.
#
# The MFEM iterator emits five probe files (see
# dynamic/tpv104_substep_iterator.cpp):
#     probe=trial_traction     : sigma_n_trial, tau1_trial, tau2_trial
#     probe=state_evolution    : psi_in, V, L, dt, V_w, a, b, V0, f0, muW, psi_out
#     probe=friction_coeff     : V, psi, a, V0, mu
#     probe=slip_rate          : Theta, psi, sigma_n, eta_s, a, V_abs
#     probe=corrected_imposed  : sigma_n_corr, tau1_corr, tau2_corr,
#                                Q_imp_plus_{vn, vt1, vt2},
#                                Q_imp_minus_{vn, vt1, vt2}
#
# Each field key below marks:
#   - "dip"    : tangent1 in BP5 convention (= tangent2 in reference)
#   - "strike" : tangent2 in BP5 convention (= tangent1 in reference)
#   - "sign_flip_for_ref" : channel whose sign must be flipped when
#                           comparing against the reference runtime
# ---------------------------------------------------------------------------
CHANNEL_ROLE: Dict[str, Dict[str, str]] = {
    "trial_traction": {
        "sigma_n_trial": "sign_flip_for_ref",
        "tau1_trial":    "dip",
        "tau2_trial":    "strike",
    },
    "friction_coeff": {
        # No coordinate-dependent channels; mu, V, psi, a, V0 are scalars.
    },
    "slip_rate": {
        "sigma_n": "sign_flip_for_ref",
        # Theta, psi, eta_s, a, V_abs are frame-independent scalars.
    },
    "corrected_imposed": {
        "sigma_n_corr": "sign_flip_for_ref",
        "tau1_corr":    "dip",
        "tau2_corr":    "strike",
        # Q_imp velocities are in fault-local frame: VX = normal,
        # VY = tangent1 (dip), VZ = tangent2 (strike).
        "Q_imp_plus_vt1":  "dip",
        "Q_imp_plus_vt2":  "strike",
        "Q_imp_minus_vt1": "dip",
        "Q_imp_minus_vt2": "strike",
    },
    "state_evolution": {
        # Pure scalars — no coordinate-dependent channel.
    },
}


def apply_reference_to_mfem(probe_name: str,
                            fields: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Translate reference-runtime probe columns into the MFEM column
    space so ``fields[k]`` lines up with MFEM's same-named channel.

    Concretely:
        - swap dip ↔ strike channel names (tangent1 ↔ tangent2),
        - negate every ``sign_flip_for_ref`` channel in-place.

    The returned dict has the MFEM-side channel names.  Unknown channels
    pass through unchanged so the caller can compose maps.
    """
    if probe_name not in CHANNEL_ROLE:
        raise KeyError(f"unknown probe name: {probe_name!r}")

    roles = CHANNEL_ROLE[probe_name]
    out: Dict[str, np.ndarray] = {}

    # Build a name-rewrite map for tangent swap.  Only channel names that
    # both exist in the reference and are marked "dip" / "strike" are
    # swapped.  The rewrite matches on a suffix so "tau1_..." ↔ "tau2_..."
    # and "Q_imp_plus_vt1" ↔ "Q_imp_plus_vt2" work.
    def swap_name(name: str) -> str:
        # Explicit swap table — only the specific TPV104 channel names.
        swap_pairs = {
            "tau1_trial": "tau2_trial",
            "tau2_trial": "tau1_trial",
            "tau1_corr":  "tau2_corr",
            "tau2_corr":  "tau1_corr",
            "Q_imp_plus_vt1":  "Q_imp_plus_vt2",
            "Q_imp_plus_vt2":  "Q_imp_plus_vt1",
            "Q_imp_minus_vt1": "Q_imp_minus_vt2",
            "Q_imp_minus_vt2": "Q_imp_minus_vt1",
        }
        return swap_pairs.get(name, name)

    for name, arr in fields.items():
        role = roles.get(name, "")
        mfem_name = swap_name(name)
        col = np.array(arr, dtype=float, copy=True)
        if role == "sign_flip_for_ref":
            col = -col
        out[mfem_name] = col

    return out


def invariant_channels(probe_name: str,
                       all_channels: Iterable[str]) -> List[str]:
    """Return the channels that are frame-independent (no swap, no sign
    flip) — useful for direct numeric comparison without any map.
    """
    if probe_name not in CHANNEL_ROLE:
        return list(all_channels)
    roles = CHANNEL_ROLE[probe_name]
    return [c for c in all_channels if roles.get(c, "") == ""]
