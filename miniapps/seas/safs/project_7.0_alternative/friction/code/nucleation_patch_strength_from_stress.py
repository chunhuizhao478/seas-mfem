#!/usr/bin/env python3
"""
nucleation_patch_strength_from_stress.py — compute the rate-and-state SHEAR
STRENGTH inside the gradual_overstress nucleation patch ANALYTICALLY, the same
way the driver does, WITHOUT reading the (ZFP-compressed) fault.vtkhdf.

  sigma_n on the fault = n . sigma_const . n  (constant_tensor stress, the
      driver's geom.sigma_n_per_dof input), sigma_n_eff = sigma_n - P_p.
  steady-state RS friction at V_init:  f_ss = f0 + (a - b) * ln(V_init / V0)
      with (a-b) = param_a_minus_b(depth), a = param_a(depth), b = a - (a-b).
  shear strength:  tau_str = f_ss * sigma_n_eff.

Compares tau_str to the projected shear |tau_pre| and the gradual_overstress
delta_tau (strike) to answer "is the overstress strong enough to nucleate?".

Reuses the stress-projection mesh + traction-resolution helpers.

Usage:
    conda activate pythonenv     # (python: numpy + meshio)
    python nucleation_patch_strength_from_stress.py
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

_HERE = Path(__file__).resolve().parent
_STRESS = _HERE.parent.parent / "stress" / "code"
sys.path.insert(0, str(_STRESS))
from project_to_fault_stress import (  # noqa: E402
    load_fault_mesh, build_fault_basis, extract_cells_by_physical,
    triangle_geometry, resolve_traction_per_cell,
)
from project_friction_to_fault import read_profile_csv, eval_profile  # noqa: E402

# --- run / config constants (from the depthprofile RS config + job 7751862) ---
FAULT_VTU = (_HERE.parent.parent / "meshing" / "results" / "vtu"
             / "safs_fault_box_nwcut_500m_lcfar3000_z0embed_fault.vtu")
A_CSV = _HERE.parent / "rate-and-state" / "param_a_vwvs11km.csv"
AMB_CSV = _HERE.parent / "rate-and-state" / "param_a_minus_b_vwvs11km.csv"

# constant_tensor stress [Pa], compression-POSITIVE SEAS convention.
SIGMA = np.array([
    [5.056619147855508e+07, 9.888543819998317e+06, 0.0],
    [9.888543819998317e+06, 1.1143380852144492e+08, 0.0],
    [0.0, 0.0, 4.50e+07],
])
P_P = 16.0e6          # [Pa]
MIN_SIGMA_N = 1.0e6   # [Pa] effective-normal-stress floor

# nucleation patch (gradual_overstress)
CX, CY, CZ = 606971.0, 3707270.0, -4965.62
RADIUS = 4000.0
DELTA_TAU = 20.0e6    # [Pa] strike overstress amplitude

# rate-and-state scalars
F0, V0, V_INIT = 0.6, 1.0e-6, 1.0e-12
SIGMA_N_FLOOR = 10.0e6  # the cap (only bites if sigma_n_eff < this)


def _mpa(x):
    return f"{x/1e6:.4g}"


def main():
    a_d, a_v = read_profile_csv(A_CSV)
    amb_d, amb_v = read_profile_csv(AMB_CSV)

    mesh = load_fault_mesh(FAULT_VTU)
    tri = extract_cells_by_physical(mesh, "triangle", "fault")
    centroids, _n, _A = triangle_geometry(mesh.points, tri)
    geom = build_fault_basis(mesh)   # normals/strikes/dips harmonised

    d = np.sqrt((centroids[:, 0] - CX) ** 2 + (centroids[:, 1] - CY) ** 2
                + (centroids[:, 2] - CZ) ** 2)
    patch = d <= RADIUS
    npatch = int(patch.sum())

    res = resolve_traction_per_cell(
        SIGMA / 1e6,                              # MPa (resolver works in MPa)
        geom.strikes, geom.dips, geom.normals,
        P_p_per_cell=P_P / 1e6,
    )
    sigma_n_eff = res["sigma_n_eff"][patch] * 1e6     # back to Pa
    tau_strike = res["tau_strike"][patch] * 1e6
    tau_dip = res["tau_dip"][patch] * 1e6
    tau_mag = res["tau_magnitude"][patch] * 1e6
    depth_km = np.maximum(0.0, -centroids[patch, 2]) / 1000.0

    a_patch = eval_profile(a_d, a_v, depth_km)
    amb_patch = eval_profile(amb_d, amb_v, depth_km)        # (a-b)
    f_ss = F0 + amb_patch * np.log(V_INIT / V0)             # steady-state f
    tau_str = f_ss * sigma_n_eff
    # driving with the overstress added to the strike channel (the nucleation
    # channel); |tau_pre + delta_tau|.
    tau_drive = np.sqrt((tau_strike + DELTA_TAU) ** 2 + tau_dip ** 2)
    overshoot = tau_drive - tau_str

    def stat(name, arr, unit=" MPa", scale=1e6):
        a = np.asarray(arr)
        print(f"  {name:26s} min={a.min()/scale:.4g}  mean={a.mean()/scale:.4g}"
              f"  max={a.max()/scale:.4g}{unit}")

    print(f"fault VTU      : {FAULT_VTU.name}")
    print(f"patch          : centre depth ~{abs(CZ)/1000:.2f} km, radius "
          f"{RADIUS:.0f} m -> {npatch} fault triangles")
    print(f"depth range    : [{depth_km.min():.2f}, {depth_km.max():.2f}] km "
          f"(all VW: a-b<0, transition at 11 km)")
    print("--- stress / friction in the nucleation patch ---")
    stat("sigma_n_eff (=n.s.n - Pp)", sigma_n_eff)
    print(f"  param a                    min={a_patch.min():.4g}  "
          f"mean={a_patch.mean():.4g}  max={a_patch.max():.4g}")
    print(f"  (a-b)                      min={amb_patch.min():.4g}  "
          f"mean={amb_patch.mean():.4g}  max={amb_patch.max():.4g}")
    print(f"  f_ss(V_init)               min={f_ss.min():.4g}  "
          f"mean={f_ss.mean():.4g}  max={f_ss.max():.4g}")
    stat("SHEAR STRENGTH f_ss*sn_eff", tau_str)
    print(f"  reference f0*sn_eff (0.6)  min={ (F0*sigma_n_eff).min()/1e6:.4g}  "
          f"mean={(F0*sigma_n_eff).mean()/1e6:.4g}  "
          f"max={(F0*sigma_n_eff).max()/1e6:.4g} MPa")
    stat("|tau_pre| (projected shear)", tau_mag)
    print("--- overstress sufficiency ---")
    print(f"  delta_tau (overstress)     = {_mpa(DELTA_TAU)} MPa")
    print(f"  delta_tau / tau_str (mean) = {DELTA_TAU/tau_str.mean():.3f}")
    stat("|tau_pre+delta_tau| driving", tau_drive)
    stat("overshoot (drive - tau_str)", overshoot)
    print(f"  patch cells with overshoot > 0 (will accelerate): "
          f"{int((overshoot>0).sum())} / {npatch}")
    # e-folds of slip-rate drive: delta_tau / ((b-a)*sigma_n)
    bma_sn = (-amb_patch) * sigma_n_eff
    efold = DELTA_TAU / bma_sn
    print(f"  overstress e-folds delta_tau/((b-a)*sn_eff): "
          f"min={efold.min():.3g} mean={efold.mean():.3g} max={efold.max():.3g}")
    print(f"  sigma_n cap (10 MPa) active in patch? "
          f"{'YES' if (sigma_n_eff < SIGMA_N_FLOOR).any() else 'no (sigma_n_eff >> 10 MPa)'}")


if __name__ == "__main__":
    main()
