#!/usr/bin/env python3
"""Project a SeisSol easi initial-stress tensor onto a PUML mesh fault
and write ParaView VTUs (fault tractions + bulk tensor) + summary JSON.

Mirrors the output of stress/code/project_to_fault_stress.py (same field
names, same Tandem fault basis, same conventions) but takes the SeisSol
run inputs directly:

    input stress : easi YAML !ConstantMap with [s_xx..s_xz]
                   (SeisSol convention: compression NEGATIVE, Pa,
                   EFFECTIVE tensor — pore pressure baked in)
    input mesh   : pumgen .puml.h5 (fault = boundary-code-3 faces)

Outputs (into --out-dir):
    <prefix>_fault_stress.vtu  per-facet + per-vertex resolved tractions
    <prefix>_bulk_stress.vtu   constant tensor on all tets (6 components)
    <prefix>_summary.json      params, global tensor, per-field stats,
                               hypocenter cross-check

Conventions (matching the reference pipeline / summary "convention" key):
    compression POSITIVE (SEAS internal) in all outputs, units MPa;
    sigma_n_eff = sigma_n_total - P_p;
    fault basis: Tandem  s = up x n,  d = s x n (down-dip);
    normals harmonised to the SW half-space via --strike-hint-az
    (right-lateral SAF) =>  tau_strike + = right-lateral,
                            tau_dip    + = reverse;
    frame (east, north, up), UTM Zone 11N, metres.

Sign conversion at the single source site:
    sigma_seas_eff[MPa] = -(easi tensor [Pa]) / 1e6
Shear tractions are P_p-invariant (s.(P_p I).n = 0), so projecting the
effective tensor gives sigma_n_eff directly and the same tau as the
total tensor;  sigma_n_total = sigma_n_eff + P_p (--P-p-MPa).

Usage:
    python3 project_stress_to_vtu.py <stress.yaml> <mesh.puml.h5>
        [--out-dir DIR] [--P-p-MPa 20] [--strike-hint-az 314]
        [--hypocenter X Y Z] [--no-bulk] [--prefix NAME]

Requires only numpy + h5py (easi YAML is parsed with a strict
ConstantMap-only reader; any LuaMap/FunctionMap/ASAGI stress aborts).
"""
import argparse
import json
import os
import re
import struct
import sys

import h5py
import numpy as np

EPS = 1.0e-9
UP = np.array([0.0, 0.0, 1.0])

# PUML Numbering<TETRAHEDRON> face->local-vertex map, verified against
# safs_mesh.puml.h5 in ../h5_to_vtu/puml_h5_to_vtu.py (free-surface
# flatness test); the relaxed check is re-run here on every conversion.
FACE_VERTS = ((1, 0, 2), (0, 1, 3), (1, 2, 3), (2, 0, 3))
BC_FAULT = 3

STRESS_KEYS = ("s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz")


# ------------------------------------------------------------ easi YAML
def read_constant_stress_yaml(path):
    """Parse the six s_* values from an easi !ConstantMap stress file.

    Strict by design: only the ConstantMap form used by this project is
    accepted.  A spatially-varying map (LuaMap/FunctionMap/ASAGI/
    AffineMap...) cannot be evaluated here and aborts with a clear
    message instead of silently mis-projecting.
    """
    text = open(path).read()
    body = re.sub(r"#.*", "", text)  # strip comments
    for tag in ("LuaMap", "FunctionMap", "ASAGI", "AffineMap",
                "EvalModel", "PolynomialMap"):
        if "!" + tag in body:
            sys.exit(f"ERROR: {os.path.basename(path)} uses !{tag}; this "
                     "tool only evaluates !ConstantMap stress tensors.")
    if "!ConstantMap" not in body:
        sys.exit(f"ERROR: no !ConstantMap found in {path}")
    vals = {}
    for key in STRESS_KEYS:
        m = re.search(rf"^\s*{key}\s*:\s*([-+0-9.eE]+)\s*$", body, re.M)
        if not m:
            sys.exit(f"ERROR: ConstantMap key {key!r} not found in {path}")
        vals[key] = float(m.group(1))
    return vals


def stress_tensor_seas_MPa(vals):
    """SeisSol compression-negative Pa -> SEAS compression-positive MPa.
    This is the single sign-flip site of the tool."""
    s = np.array([
        [vals["s_xx"], vals["s_xy"], vals["s_xz"]],
        [vals["s_xy"], vals["s_yy"], vals["s_yz"]],
        [vals["s_xz"], vals["s_yz"], vals["s_zz"]],
    ])
    return -s / 1.0e6


# ------------------------------------------------------------ PUML mesh
def load_puml(path):
    with h5py.File(path, "r") as f:
        geom = f["geometry"][:]
        conn = f["connect"][:].astype(np.int64)
        group = f["group"][:].astype(np.int32)
        fmt = f.attrs.get("boundary-format", "i32")
        if isinstance(fmt, bytes):
            fmt = fmt.decode()
        bnd = f["boundary"][:]
    bits = {"i32": 8, "i64": 16}.get(fmt)
    if np.asarray(bnd).ndim == 2:
        bc = np.asarray(bnd).astype(np.int64)
    elif bits is None:
        sys.exit(f"ERROR: unsupported boundary-format {fmt!r}")
    else:
        b = bnd.astype(np.int64)
        mask = (1 << bits) - 1
        bc = np.stack([(b >> (bits * k)) & mask for k in range(4)], axis=1)
    return geom, conn, group, bc


def sanity_check_face_ordering(geom, conn, bc):
    """BC-1 (free surface) faces must be flat at z=z_max; a wrong
    FACE_VERTS would pick wrong triangles of the boundary tets."""
    sel = bc == 1
    if not sel.any():
        return
    ztop = geom[:, 2].max()
    tol = 1e-6 * max(float(geom[:, 2].max() - geom[:, 2].min()), 1.0)
    bad = 0
    for k in range(4):
        tri = conn[sel[:, k]][:, FACE_VERTS[k]]
        bad += int(np.sum(~np.all(np.abs(geom[tri][:, :, 2] - ztop) < tol,
                                  axis=1)))
    if bad:
        sys.exit(f"ERROR: face-ordering sanity check failed ({bad} "
                 "free-surface faces off the top plane); PUML face "
                 "numbering mismatch — output would be wrong.")


def extract_fault(geom, conn, bc):
    """Deduplicated BC-3 triangles with a compacted point array."""
    tris = np.concatenate(
        [conn[bc[:, k] == BC_FAULT][:, FACE_VERTS[k]] for k in range(4)],
        axis=0)
    if len(tris) == 0:
        sys.exit("ERROR: mesh has no BC-3 (dynamic rupture) faces")
    _, keep = np.unique(np.sort(tris, axis=1), axis=0, return_index=True)
    tris = tris[np.sort(keep)]
    used, inv = np.unique(tris.ravel(), return_inverse=True)
    return geom[used], inv.reshape(tris.shape)


# ----------------------------------------------------- fault basis math
def triangle_geometry(points, tris):
    p = points[tris]
    centroids = p.mean(axis=1)
    cross = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    twice_area = np.linalg.norm(cross, axis=1)
    areas = 0.5 * twice_area
    degen = twice_area < EPS
    normals = cross / np.where(degen, 1.0, twice_area)[:, None]
    normals[degen] = np.nan
    return centroids, normals, areas


def harmonise_normals(normals, strike_hint_az_deg):
    """Flip facet normals into the half-space 'left of strike in map
    view' (SW for the NW-striking right-lateral SAF):
    n_global = up x s_hint, s_hint = (sin az, cos az, 0)."""
    az = np.radians(strike_hint_az_deg)
    n_global = np.array([-np.cos(az), np.sin(az), 0.0])
    out = normals.copy()
    finite = ~np.isnan(normals).any(axis=1)
    flip = finite & (out @ n_global < 0)
    out[flip] = -out[flip]
    return out, n_global


def tandem_basis(normals):
    """s = up x n (strike), d = s x n (down-dip); NaN where degenerate."""
    s_raw = np.cross(np.broadcast_to(UP, normals.shape), normals)
    s_norm = np.linalg.norm(s_raw, axis=1)
    degen = ~(s_norm > np.sqrt(2e-6))      # near-horizontal facets
    strikes = np.full_like(normals, np.nan)
    dips = np.full_like(normals, np.nan)
    good = ~degen
    strikes[good] = s_raw[good] / s_norm[good, None]
    dips[good] = np.cross(strikes[good], normals[good])
    return strikes, dips, degen


def resolve_tractions(sigma_eff, strikes, dips, normals, P_p_MPa):
    """Project the constant effective SEAS tensor onto each facet."""
    t = normals @ sigma_eff.T                      # traction vec (N,3), MPa
    sigma_n_eff = np.einsum("ki,ki->k", normals, t)
    sigma_n_total = sigma_n_eff + P_p_MPa
    tau_strike = np.einsum("ki,ki->k", strikes, t)
    tau_dip = np.einsum("ki,ki->k", dips, t)
    tau_mag = np.hypot(tau_strike, tau_dip)
    rake = np.degrees(np.arctan2(tau_dip, tau_strike))
    mu = np.full(len(t), np.nan)
    ok = sigma_n_eff > 0
    mu[ok] = tau_mag[ok] / sigma_n_eff[ok]
    return {
        "sigma_n_total": sigma_n_total, "sigma_n_eff": sigma_n_eff,
        "tau_strike": tau_strike, "tau_dip": tau_dip,
        "tau_magnitude": tau_mag, "rake_deg": rake, "mu_apparent": mu,
        "traction_vec": t,
    }


def cell_to_node_average(n_pts, tris, cell_values, areas):
    """Area-weighted cell->vertex average; NaN cells excluded."""
    cv = np.asarray(cell_values, dtype=np.float64)
    trailing = cv.shape[1:]
    finite = np.isfinite(cv.reshape(len(cv), -1)).all(axis=1) \
        & np.isfinite(areas)
    w_cell = np.where(finite, areas, 0.0)
    cv_safe = np.where(finite.reshape((-1,) + (1,) * len(trailing)), cv, 0.0)
    w = np.zeros(n_pts)
    np.add.at(w, tris.ravel(), np.repeat(w_cell, 3))
    acc = np.zeros((n_pts,) + trailing)
    contrib = (w_cell.reshape((-1,) + (1,) * len(trailing)) * cv_safe)
    np.add.at(acc, tris.ravel(),
              np.repeat(contrib[:, None], 3, axis=1
                        ).reshape((-1,) + trailing))
    w_b = w.reshape((-1,) + (1,) * len(trailing))
    return np.divide(acc, w_b, out=np.full_like(acc, np.nan), where=w_b > 0)


def basis_to_node(n_pts, tris, strikes, normals, areas):
    """Per-vertex orthonormal (s, d, n): average n, Gram-Schmidt s, d=sxn."""
    n_avg = cell_to_node_average(n_pts, tris, normals, areas)
    n_norm = np.linalg.norm(n_avg, axis=1)
    n_node = np.full_like(n_avg, np.nan)
    ok = n_norm > EPS
    n_node[ok] = n_avg[ok] / n_norm[ok, None]
    s_avg = cell_to_node_average(n_pts, tris, strikes, areas)
    dot = np.einsum("ij,ij->i", np.nan_to_num(s_avg), np.nan_to_num(n_node))
    s_proj = s_avg - dot[:, None] * n_node
    s_norm = np.linalg.norm(s_proj, axis=1)
    s_node = np.full_like(s_avg, np.nan)
    ok2 = ok & (s_norm > EPS)
    s_node[ok2] = s_proj[ok2] / s_norm[ok2, None]
    d_node = np.cross(s_node, n_node)
    return s_node, d_node, n_node


# ------------------------------------------------------------ VTU writer
_VTK_TYPE = {np.dtype("<f8"): "Float64", np.dtype("<f4"): "Float32",
             np.dtype("<i8"): "Int64", np.dtype("<i4"): "Int32",
             np.dtype("u1"): "UInt8"}


def write_vtu(path, points, cells, cell_type, cell_data=None,
              point_data=None):
    """XML VTU, appended raw binary, little endian, UInt64 headers."""
    points = np.ascontiguousarray(points, dtype="<f8")
    cells = np.ascontiguousarray(cells, dtype="<i8")
    n_pts, n_cells = len(points), len(cells)
    offsets = np.arange(1, n_cells + 1, dtype="<i8") * cells.shape[1]
    types = np.full(n_cells, cell_type, dtype="u1")

    def prep(arr):
        a = np.ascontiguousarray(arr)
        if a.dtype not in _VTK_TYPE:
            a = a.astype("<f8" if a.dtype.kind == "f" else "<i4")
        ncomp = 1 if a.ndim == 1 else a.shape[1]
        return a, _VTK_TYPE[a.dtype], ncomp

    # connectivity MUST be a flat single-component array — declaring it
    # with NumberOfComponents=k makes vtkXMLUnstructuredGridReader fail
    # ("could not be created with one component")
    blocks = [("Points", *prep(points)),
              ("connectivity", *prep(cells.reshape(-1))),
              ("offsets", *prep(offsets)), ("types", *prep(types))]
    n_fixed = len(blocks)
    cd = [(k, *prep(v)) for k, v in (cell_data or {}).items()]
    pd = [(k, *prep(v)) for k, v in (point_data or {}).items()]
    raw_all = blocks + cd + pd
    offs, pos = [], 0
    for _, a, _, _ in raw_all:
        offs.append(pos)
        pos += 8 + a.nbytes

    def da(name, vtype, ncomp, off):
        comp = f' NumberOfComponents="{ncomp}"' if ncomp > 1 else ""
        return (f'        <DataArray type="{vtype}" Name="{name}"{comp}'
                f' format="appended" offset="{off}"/>\n')

    xml = ['<?xml version="1.0"?>\n',
           '<VTKFile type="UnstructuredGrid" version="1.0" '
           'byte_order="LittleEndian" header_type="UInt64">\n',
           '  <UnstructuredGrid>\n',
           f'    <Piece NumberOfPoints="{n_pts}" NumberOfCells="{n_cells}">\n',
           '      <Points>\n', da("Points", "Float64", 3, offs[0]),
           '      </Points>\n', '      <Cells>\n',
           da("connectivity", "Int64", 1, offs[1]),
           da("offsets", "Int64", 1, offs[2]),
           da("types", "UInt8", 1, offs[3]), '      </Cells>\n']
    for tag, items, off0 in (("CellData", cd, n_fixed),
                             ("PointData", pd, n_fixed + len(cd))):
        if items:
            xml.append(f'      <{tag}>\n')
            for (name, _, vtype, ncomp), off in zip(
                    items, offs[off0:off0 + len(items)]):
                xml.append(da(name, vtype, ncomp, off))
            xml.append(f'      </{tag}>\n')
    xml += ['    </Piece>\n', '  </UnstructuredGrid>\n',
            '  <AppendedData encoding="raw">\n_']
    with open(path, "wb") as fh:
        fh.write("".join(xml).encode())
        for _, a, _, _ in raw_all:
            fh.write(struct.pack("<Q", a.nbytes))
            fh.write(a.tobytes())
        fh.write(b"\n  </AppendedData>\n</VTKFile>\n")


# ----------------------------------------------------------------- main
def field_stats(arr):
    finite = np.isfinite(arr)
    return {"min": float(np.min(arr[finite])),
            "median": float(np.median(arr[finite])),
            "max": float(np.max(arr[finite])),
            "nan_count": int((~finite).sum())}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("stress_yaml", help="easi !ConstantMap stress file "
                    "(SeisSol convention: compression-negative Pa, effective)")
    ap.add_argument("mesh", help="pumgen .puml.h5 mesh")
    ap.add_argument("--out-dir", default=".")
    ap.add_argument("--prefix", default=None,
                    help="output prefix (default <yaml>_on_<mesh>)")
    ap.add_argument("--P-p-MPa", type=float, default=0.0, dest="P_p",
                    help="pore pressure baked into the effective tensor; "
                    "used only to reconstruct sigma_n_total (default 0)")
    ap.add_argument("--strike-hint-az", type=float, default=314.0,
                    help="global strike azimuth (deg) for normal "
                    "harmonisation; 314 = SAF N46W right-lateral")
    ap.add_argument("--hypocenter", type=float, nargs=3, default=None,
                    metavar=("X", "Y", "Z"),
                    help="report resolved tractions at the facet nearest "
                    "this point (cross-check against the fault yaml)")
    ap.add_argument("--no-bulk", action="store_true",
                    help="skip the (constant) bulk tensor VTU")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    ystem = os.path.splitext(os.path.basename(args.stress_yaml))[0]
    mstem = os.path.basename(args.mesh)
    for suf in (".h5", ".puml"):
        if mstem.endswith(suf):
            mstem = mstem[:-len(suf)]
    prefix = args.prefix or f"{ystem}_on_{mstem}"

    vals = read_constant_stress_yaml(args.stress_yaml)
    sigma_eff = stress_tensor_seas_MPa(vals)        # SEAS comp-positive MPa
    sigma_total = sigma_eff + args.P_p * np.eye(3)
    print("sigma_eff (SEAS, MPa):\n", np.round(sigma_eff, 3))

    geom, conn, group, bc = load_puml(args.mesh)
    sanity_check_face_ordering(geom, conn, bc)
    pts, tris = extract_fault(geom, conn, bc)
    print(f"fault: {len(tris)} triangles, {len(pts)} vertices")

    centroids, normals_raw, areas = triangle_geometry(pts, tris)
    normals, n_global = harmonise_normals(normals_raw, args.strike_hint_az)
    strikes, dips, degen = tandem_basis(normals)
    n_degen = int(degen.sum())
    if n_degen:
        print(f"warning: {n_degen} degenerate/near-horizontal facets "
              "(NaN basis)", file=sys.stderr)

    r = resolve_tractions(sigma_eff, strikes, dips, normals, args.P_p)

    # nodal (continuous) fields: average components, recompute derived
    n_pts = len(pts)
    node = {k: cell_to_node_average(n_pts, tris, r[k], areas)
            for k in ("sigma_n_total", "sigma_n_eff",
                      "tau_strike", "tau_dip", "traction_vec")}
    tau_mag_node = np.hypot(node["tau_strike"], node["tau_dip"])
    rake_node = np.degrees(np.arctan2(node["tau_dip"], node["tau_strike"]))
    mu_node = np.full(n_pts, np.nan)
    okn = node["sigma_n_eff"] > 0
    mu_node[okn] = tau_mag_node[okn] / node["sigma_n_eff"][okn]
    s_node, d_node, n_node = basis_to_node(n_pts, tris, strikes,
                                           normals, areas)

    fault_path = os.path.join(args.out_dir, f"{prefix}_fault_stress.vtu")
    write_vtu(
        fault_path, pts, tris, 5,
        cell_data={
            "sigma_n_total_MPa_cell": r["sigma_n_total"],
            "sigma_n_eff_MPa_cell": r["sigma_n_eff"],
            "tau_strike_MPa_cell": r["tau_strike"],
            "tau_dip_MPa_cell": r["tau_dip"],
            "tau_magnitude_MPa_cell": r["tau_magnitude"],
            "rake_deg_cell": r["rake_deg"],
            "mu_apparent_cell": r["mu_apparent"],
            "strike_vec_cell": strikes, "dip_vec_cell": dips,
            "normal_vec_cell": normals,
        },
        point_data={
            "sigma_n_total_MPa": node["sigma_n_total"],
            "sigma_n_eff_MPa": node["sigma_n_eff"],
            "tau_strike_MPa": node["tau_strike"],
            "tau_dip_MPa": node["tau_dip"],
            "tau_magnitude_MPa": tau_mag_node,
            "rake_deg": rake_node,
            "mu_apparent": mu_node,
            "strike_vec": s_node, "dip_vec": d_node, "normal_vec": n_node,
            "traction_vec_MPa": node["traction_vec"],
        })
    print(f"wrote {fault_path}")

    bulk_n_cells = len(conn)
    if not args.no_bulk:
        comps = {"sigma_xx_MPa": sigma_total[0, 0],
                 "sigma_yy_MPa": sigma_total[1, 1],
                 "sigma_zz_MPa": sigma_total[2, 2],
                 "sigma_xy_MPa": sigma_total[0, 1],
                 "sigma_xz_MPa": sigma_total[0, 2],
                 "sigma_yz_MPa": sigma_total[1, 2]}
        bulk_cd = {k: np.full(bulk_n_cells, v, dtype="<f4")
                   for k, v in comps.items()}
        bulk_cd["group"] = group
        bulk_path = os.path.join(args.out_dir, f"{prefix}_bulk_stress.vtu")
        write_vtu(bulk_path, geom, conn, 10, cell_data=bulk_cd)
        print(f"wrote {bulk_path}  (constant total tensor, float32)")

    summary = {
        # paths recorded as given on the command line (no machine-
        # specific absolute paths in the artifact)
        "input_stress_yaml": args.stress_yaml,
        "input_mesh": args.mesh,
        "fault_n_cells": int(len(tris)),
        "fault_n_degenerate_basis": n_degen,
        "bulk_n_cells": int(bulk_n_cells),
        "params": {
            "P_p_MPa": args.P_p,
            "fault_strike_azimuth_hint_deg": args.strike_hint_az,
            "rake_sense": "right-lateral",
            "harmonisation_normal_global": n_global.tolist(),
            "source_convention": "SeisSol easi ConstantMap, "
                                 "compression-NEGATIVE Pa, EFFECTIVE "
                                 "(P_p baked in); converted here via "
                                 "sigma_seas_eff = -sigma/1e6",
        },
        "sigma_eff_global_MPa": sigma_eff.tolist(),
        "sigma_total_global_MPa": sigma_total.tolist(),
        "stats": {
            "sigma_n_total_MPa": field_stats(r["sigma_n_total"]),
            "sigma_n_eff_MPa": field_stats(r["sigma_n_eff"]),
            "tau_strike_MPa": field_stats(r["tau_strike"]),
            "tau_dip_MPa": field_stats(r["tau_dip"]),
            "tau_magnitude_MPa": field_stats(r["tau_magnitude"]),
            "rake_deg": field_stats(r["rake_deg"]),
            "mu_apparent": field_stats(r["mu_apparent"]),
        },
        "convention": "compression POSITIVE (SEAS internal); P_p positive; "
                      "sigma_n_eff = sigma_n_total - P_p; tau_strike + = "
                      "right-lateral (under Tandem fault basis s = up x n, "
                      "d = s x n); frame (east, north, up) UTM Zone 11N; "
                      "units MPa",
    }

    if args.hypocenter is not None:
        hc = np.asarray(args.hypocenter)
        i = int(np.nanargmin(np.linalg.norm(centroids - hc, axis=1)))
        summary["hypocenter_check"] = {
            "requested_xyz": hc.tolist(),
            "nearest_facet_centroid": centroids[i].tolist(),
            "distance_m": float(np.linalg.norm(centroids[i] - hc)),
            "sigma_n_eff_MPa": float(r["sigma_n_eff"][i]),
            "tau_magnitude_MPa": float(r["tau_magnitude"][i]),
            "tau_strike_MPa": float(r["tau_strike"][i]),
            "tau_dip_MPa": float(r["tau_dip"][i]),
            "mu_apparent": float(r["mu_apparent"][i]),
        }
        print("hypocenter facet:", json.dumps(
            summary["hypocenter_check"], indent=2))

    sum_path = os.path.join(args.out_dir, f"{prefix}_summary.json")
    with open(sum_path, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"wrote {sum_path}")
    for k, st in summary["stats"].items():
        print(f"  {k:22s} min {st['min']:+9.4f}  med {st['median']:+9.4f}  "
              f"max {st['max']:+9.4f}  nan {st['nan_count']}")


if __name__ == "__main__":
    main()
