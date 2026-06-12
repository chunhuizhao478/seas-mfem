#!/usr/bin/env python3
"""Compute the SAFS QuakeWorx "hypocenter facts" sheet from the actual
run inputs of a SeisSol case folder.

At the nucleation location (hypocenter) this reports:
  - on-fault stress  : effective normal stress sigma_N, shear stress
    magnitude |tau|, apparent friction mu_app = |tau| / sigma_N
    (projection of the easi !ConstantMap initial stress onto the
    nearest fault facet of the PUML mesh — same math/conventions as
    ../on_fault_stress_projection/project_stress_to_vtu.py, which is
    imported as a library)
  - material         : rho, shear modulus mu, lambda from the CVM
    ASAGI netCDF (trilinear interpolation, matching easi
    `interpolation: linear`), and derived Vs = sqrt(mu/rho),
    Vp = sqrt((lambda+2 mu)/rho), nu = lambda / (2 (lambda+mu))

Spec-sheet reference values (safs_seisol_v2_0_0_RSSRW, 2026-06-11) are
built in and printed alongside for a quick drift check:

    Nucleation Location: (606971, 3707270, -4965.62) [Hypocenter]
    Effective normal stress sigma_N: 63.8 MPa
    Shear stress magnitude: 21.2 MPa
    Apparent friction (mu): 0.332
    Shear wave speed (Vs): 3032 m/s
    Density of the material (rho): 2551 kg/m^3
    Shear modulus (mu): 23.5 GPa
    Poisson's ratio (nu): 0.275

Usage:
    python3 hypocenter_facts.py                       # default case dir
    python3 hypocenter_facts.py --case-dir DIR [--hypocenter X Y Z]
        [--P-p-MPa 20] [--json OUT.json]

Requires numpy + h5py only (the .nc is netCDF4 = HDF5).
"""
import argparse
import json
import os
import sys

import h5py
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "on_fault_stress_projection"))
from project_stress_to_vtu import (  # noqa: E402
    read_constant_stress_yaml, stress_tensor_seas_MPa, load_puml,
    sanity_check_face_ordering, extract_fault, triangle_geometry,
    harmonise_normals, tandem_basis, resolve_tractions,
)

# quakeworx case-folder root (= the parent of toolbox/), resolved
# relative to this script so the folder stays relocatable
_ROOT = os.path.normpath(os.path.join(_HERE, "..", ".."))
DEFAULT_CASE_DIR = os.path.join(_ROOT, "safs_seisol_v2_0_0_RSSRW")
DEFAULT_HYPOCENTER = (606971.0, 3707270.0, -4965.62)


def portable_path(p):
    """Path relative to the quakeworx root if inside it (keeps uploaded
    artifacts free of machine-specific absolute paths), else as given."""
    rel = os.path.relpath(os.path.abspath(p), _ROOT)
    return p if rel.startswith("..") else rel

# Spec sheet as circulated 2026-06-11 (source: safs_fault.yaml header
# comments + QuakeWorx case description).
SPEC = {
    "sigma_n_eff_MPa": 63.8,
    "tau_magnitude_MPa": 21.2,
    "mu_apparent": 0.332,
    "Vs_m_s": 3032.0,
    "rho_kg_m3": 2551.0,
    "shear_modulus_GPa": 23.5,
    "poisson_ratio": 0.275,
}

NOTE = ("Note: Nucleation must occur near the south-east end of the fault; "
        "that is the only hard constraint on its location. You may vary the "
        "hypocenter freely within that region, but expect to adjust the "
        "friction parameters to achieve nucleation under the Section 5 "
        "constraints, since both the velocity and stress fields are "
        "heterogeneous along strike and with depth.")


def sample_cvm(nc_path, xyz, mode):
    """Sample (rho, mu, lambda) at xyz from the ASAGI netCDF (compound
    variable data[z, y, x]).

    mode='trilinear' mirrors easi's `interpolation: linear` (what the
    SeisSol run actually sees, clamp-to-boundary); mode='nearest' takes
    the closest grid node."""
    x0, y0, z0 = xyz
    with h5py.File(nc_path, "r") as f:
        xs, ys, zs = f["x"][:], f["y"][:], f["z"][:]
        if mode == "nearest":
            ix = int(np.argmin(np.abs(xs - x0)))
            iy = int(np.argmin(np.abs(ys - y0)))
            iz = int(np.argmin(np.abs(zs - z0)))
            node = f["data"][iz, iy, ix]
            out = {k: float(node[k]) for k in ("rho", "mu", "lambda")}
            out["_grid"] = {"ix": ix, "iy": iy, "iz": iz, "mode": mode,
                            "node_xyz": [float(xs[ix]), float(ys[iy]),
                                         float(zs[iz])]}
            return out

        def bracket(coords, v):
            i = int(np.clip(np.searchsorted(coords, v) - 1,
                            0, len(coords) - 2))
            t = (v - coords[i]) / (coords[i + 1] - coords[i])
            return i, float(np.clip(t, 0.0, 1.0))

        ix, tx = bracket(xs, x0)
        iy, ty = bracket(ys, y0)
        iz, tz = bracket(zs, z0)
        corner = f["data"][iz:iz + 2, iy:iy + 2, ix:ix + 2]  # (2,2,2) compound

    out = {}
    wz = np.array([1 - tz, tz])
    wy = np.array([1 - ty, ty])
    wx = np.array([1 - tx, tx])
    w = wz[:, None, None] * wy[None, :, None] * wx[None, None, :]
    for name in ("rho", "mu", "lambda"):
        out[name] = float(np.sum(w * corner[name].astype(np.float64)))
    out["_grid"] = {"ix": ix, "iy": iy, "iz": iz, "mode": mode,
                    "tx": tx, "ty": ty, "tz": tz}
    return out


def derive_material(mat):
    rho, mu, lam = mat["rho"], mat["mu"], mat["lambda"]
    return {
        "Vs_m_s": float(np.sqrt(mu / rho)),
        "Vp_m_s": float(np.sqrt((lam + 2.0 * mu) / rho)),
        "shear_modulus_GPa": mu / 1e9,
        "lambda_GPa": lam / 1e9,
        "rho_kg_m3": rho,
        "poisson_ratio": lam / (2.0 * (lam + mu)),
    }


def fault_stress_at(stress_yaml, mesh_path, xyz, strike_hint_az, P_p_MPa):
    """sigma_N_eff / |tau| / mu_app at the fault facet nearest xyz."""
    sigma_eff = stress_tensor_seas_MPa(read_constant_stress_yaml(stress_yaml))
    geom, conn, _group, bc = load_puml(mesh_path)
    sanity_check_face_ordering(geom, conn, bc)
    pts, tris = extract_fault(geom, conn, bc)
    centroids, normals_raw, _areas = triangle_geometry(pts, tris)
    normals, _ = harmonise_normals(normals_raw, strike_hint_az)
    strikes, dips, _degen = tandem_basis(normals)
    r = resolve_tractions(sigma_eff, strikes, dips, normals, P_p_MPa)
    i = int(np.nanargmin(np.linalg.norm(centroids - np.asarray(xyz), axis=1)))
    return {
        "facet_index": i,
        "facet_centroid": centroids[i].tolist(),
        "facet_distance_m": float(np.linalg.norm(centroids[i] - xyz)),
        "sigma_n_eff_MPa": float(r["sigma_n_eff"][i]),
        "sigma_n_total_MPa": float(r["sigma_n_total"][i]),
        "tau_magnitude_MPa": float(r["tau_magnitude"][i]),
        "tau_strike_MPa": float(r["tau_strike"][i]),
        "tau_dip_MPa": float(r["tau_dip"][i]),
        "rake_deg": float(r["rake_deg"][i]),
        "mu_apparent": float(r["mu_apparent"][i]),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--case-dir", default=DEFAULT_CASE_DIR,
                    help="SeisSol case folder (default: "
                    "../../safs_seisol_v2_0_0_RSSRW relative to this script)")
    ap.add_argument("--stress-yaml", default=None,
                    help="override <case-dir>/safs_initial_stress.yaml")
    ap.add_argument("--mesh", default=None,
                    help="override <case-dir>/safs_mesh.puml.h5")
    ap.add_argument("--material-nc", default=None,
                    help="override <case-dir>/safs_material_cvm.nc")
    ap.add_argument("--hypocenter", type=float, nargs=3,
                    default=list(DEFAULT_HYPOCENTER), metavar=("X", "Y", "Z"))
    ap.add_argument("--P-p-MPa", type=float, default=20.0, dest="P_p",
                    help="pore pressure baked into the effective stress "
                    "tensor (default 20, per safs_initial_stress.yaml)")
    ap.add_argument("--strike-hint-az", type=float, default=314.0)
    ap.add_argument("--sampling", choices=("trilinear", "nearest"),
                    default="trilinear",
                    help="CVM sampling for the headline material values; "
                    "trilinear = what easi `interpolation: linear` feeds "
                    "SeisSol (default); the other mode is always shown "
                    "in the provenance lines")
    ap.add_argument("--json", default=None, help="also write full-precision "
                    "facts + provenance to this JSON file")
    args = ap.parse_args()

    cd = args.case_dir
    stress_yaml = args.stress_yaml or os.path.join(cd, "safs_initial_stress.yaml")
    mesh = args.mesh or os.path.join(cd, "safs_mesh.puml.h5")
    material_nc = args.material_nc or os.path.join(cd, "safs_material_cvm.nc")
    for p in (stress_yaml, mesh, material_nc):
        if not os.path.isfile(p):
            sys.exit(f"ERROR: input not found: {p}")
    hyp = tuple(args.hypocenter)

    st = fault_stress_at(stress_yaml, mesh, hyp, args.strike_hint_az, args.P_p)
    other_mode = "nearest" if args.sampling == "trilinear" else "trilinear"
    mat = sample_cvm(material_nc, hyp, args.sampling)
    mat_other = sample_cvm(material_nc, hyp, other_mode)
    derived = derive_material(mat)
    derived_other = derive_material(mat_other)

    print(f"Nucleation Location: ({hyp[0]:.10g}, {hyp[1]:.10g}, "
          f"{hyp[2]:.10g}) [Hypocenter]")
    print()
    print(NOTE)
    print()
    rows = [
        ("Effective normal stress sigma_N", st["sigma_n_eff_MPa"], "MPa",
         "{:.1f}", "sigma_n_eff_MPa"),
        ("Shear stress magnitude", st["tau_magnitude_MPa"], "MPa",
         "{:.1f}", "tau_magnitude_MPa"),
        ("Apparent friction (mu)", st["mu_apparent"], "",
         "{:.3f}", "mu_apparent"),
        ("Shear wave speed (Vs)", derived["Vs_m_s"], "m/s",
         "{:.0f}", "Vs_m_s"),
        ("Density of the material (rho)", derived["rho_kg_m3"], "kg/m^3",
         "{:.0f}", "rho_kg_m3"),
        ("Shear modulus (mu)", derived["shear_modulus_GPa"], "GPa",
         "{:.1f}", "shear_modulus_GPa"),
        ("Poisson's ratio (nu)", derived["poisson_ratio"], "",
         "{:.3f}", "poisson_ratio"),
    ]
    for label, val, unit, fmt, key in rows:
        shown = fmt.format(val) + (f" {unit}" if unit else "")
        spec = SPEC.get(key)
        spec_s = fmt.format(spec) + (f" {unit}" if unit else "")
        flag = "" if fmt.format(val) == fmt.format(spec) \
            else f"   <-- spec sheet says {spec_s}"
        print(f"{label}: {shown}{flag}")
    print()
    print(f"[provenance] fault facet #{st['facet_index']} at "
          f"({st['facet_centroid'][0]:.1f}, {st['facet_centroid'][1]:.1f}, "
          f"{st['facet_centroid'][2]:.1f}), {st['facet_distance_m']:.1f} m "
          f"from the hypocenter; rake {st['rake_deg']:.1f} deg "
          f"(tau_strike {st['tau_strike_MPa']:.2f}, tau_dip "
          f"{st['tau_dip_MPa']:.2f} MPa); P_p = {args.P_p:g} MPa")
    print(f"[provenance] CVM {args.sampling} at grid cell "
          f"(iz,iy,ix)=({mat['_grid']['iz']},{mat['_grid']['iy']},"
          f"{mat['_grid']['ix']}), Vp = {derived['Vp_m_s']:.0f} m/s, "
          f"lambda = {derived['lambda_GPa']:.1f} GPa")
    print(f"[provenance] CVM {other_mode} (for comparison): "
          f"rho {derived_other['rho_kg_m3']:.0f} kg/m^3, "
          f"Vs {derived_other['Vs_m_s']:.0f} m/s, "
          f"mu {derived_other['shear_modulus_GPa']:.1f} GPa, "
          f"nu {derived_other['poisson_ratio']:.3f}")

    if args.json:
        out = {
            "hypocenter": list(hyp),
            "note": NOTE,
            "stress": st,
            "material_sampling": args.sampling,
            "material": {k: v for k, v in mat.items() if k != "_grid"},
            "material_grid": mat["_grid"],
            "derived": derived,
            "derived_" + other_mode: derived_other,
            "spec_sheet": SPEC,
            "inputs": {"stress_yaml": portable_path(stress_yaml),
                       "mesh": portable_path(mesh),
                       "material_nc": portable_path(material_nc),
                       "P_p_MPa": args.P_p,
                       "strike_hint_az_deg": args.strike_hint_az},
        }
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"[provenance] JSON written to {args.json}")


if __name__ == "__main__":
    main()
