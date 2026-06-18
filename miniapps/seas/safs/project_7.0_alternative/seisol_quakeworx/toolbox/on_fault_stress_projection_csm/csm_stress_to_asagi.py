#!/usr/bin/env python3
"""csm_stress_to_asagi.py — CSM C1 k=2.39 on-fault stress field -> SeisSol
easi !ASAGI NetCDF (volumetric global Cartesian effective-stress tensor).

V2 fed SeisSol a single CONSTANT stress tensor (ConstantMap) and let SeisSol
project it per fault facet.  V3 adopts the depth-dependent, CSM-orientation
C1 k=2.39 field, which is SPATIALLY VARYING, so the supported mechanism is a
gridded easi !ASAGI field (same workflow as the material safs_material_cvm.nc).

The CSM effective tensor is a genuine VOLUMETRIC field:
  - orientation (principal axes U) and shape ratio R: the CSM (YHSM-2013)
    tensor, linear-interpolated in (x,y), depth-invariant;
  - magnitude: sig2 = Sv_eff(z) = MUSCAL lithostat - hydrostatic P_p,
    sig3 = Sv_eff / ((1-R)*k + R), sig1 = k*sig3   (closure C1, k=2.39).
So sigma_eff(x,y,z) is evaluated on a rectilinear grid and SeisSol projects
it per facet exactly as it did the V2 constant tensor (BaseDRInitializer
rotateStressToFaultCS).  Output is comp-NEGATIVE Pa (SeisSol convention);
P_p is already removed (effective tensor) -> baked-in, exact under projection.

Self-check (always runs): trilinearly sample the written nc at the fault
facet centroids, project onto the (harmonised) facet normals, and compare
sigma_n_eff / tau to the C1 projection VTU (ground truth).  Non-zero exit if
the per-facet medians disagree beyond --selfcheck-tol.

Reuses project_csm_stress_to_vtu.py (same tensor builder, same conventions).
Requires numpy + h5py + scipy + pyproj + netCDF4 (conda env 'pythonenv').

Usage:
  python3 csm_stress_to_asagi.py \
    --csm-csv ../../raw_data/yang_and_hauksson_orientation/CSM_data_1781290479811.csv \
    --material-nc ../../safs_seisol_v2_0_0_RSSRW/safs_material_cvm.nc \
    --mesh ../../safs_seisol_v2_0_0_RSSRW/safs_mesh.puml.h5 \
    --verify-vtu csm_yhsm2013_stress_on_safs_mesh_fault_stress.vtu \
    --out ../../safs_seisol_v3_0_0_LSW/safs_stress_csm.nc \
    [--k-ratio 2.39] [--dx 1000] [--dz 250]
    [--xmin 350000 --xmax 630000 --ymin 3680000 --ymax 3850000
     --zmin -17000 --zmax 0]
"""
import argparse
import os
import re
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from project_csm_stress_to_vtu import (  # noqa: E402
    read_csm_csv, csm_tensors_tension, csm_axes_and_shape,
    interpolate_csm_field, geographic_to_utm11n, muscal_sv_total_profile,
    sv_total_at, pore_pressure, magnitudes_C1, build_tensor_from_axes,
    load_puml, sanity_check_face_ordering, extract_fault, triangle_geometry,
    harmonise_normals, tandem_basis, resolve_tractions, field_stats, EPS)
from fault_local_overpressure import (  # noqa: E402
    load_overpressure_field, delta_pp_mpa, apply_overpressure_to_tensor)

_VTK_RD = {"Float64": "<f8", "Float32": "<f4", "Int64": "<i8",
           "Int32": "<i4", "UInt8": "u1"}
STRESS_FIELDS = ["s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz"]


def read_vtu_cell(path, names):
    with open(path, "rb") as fh:
        raw = fh.read()
    marker = b'<AppendedData encoding="raw">\n_'
    mk = raw.find(marker)
    start = mk + len(marker)
    hdr = raw[:mk].decode("latin1")
    out = {}
    for typ, nm, nc, off in re.findall(
            r'<DataArray type="([^"]+)" Name="([^"]+)"'
            r'(?: NumberOfComponents="(\d+)")? format="appended" '
            r'offset="(\d+)"/>', hdr):
        if nm not in names:
            continue
        off = int(off)
        nc = int(nc) if nc else 1
        nbytes = struct.unpack_from("<Q", raw, start + off)[0]
        a = np.frombuffer(raw[start + off + 8:start + off + 8 + nbytes],
                          dtype=_VTK_RD[typ])
        out[nm] = a.reshape(-1, nc) if nc > 1 else a
    return out


def write_stress_asagi(path, x, y, z, comps, attrs):
    """Write comps={s_xx..s_xz} each shaped (nx,ny,nz) -> compound nc
    data{...}(z,y,x), float32, mirroring convert_cvm_to_asagi."""
    from netCDF4 import Dataset
    nx, ny, nz = len(x), len(y), len(z)
    st = np.dtype([(f, np.float32) for f in STRESS_FIELDS])
    with Dataset(str(path), "w", format="NETCDF4") as ds:
        ds.createDimension("x", nx)
        ds.createDimension("y", ny)
        ds.createDimension("z", nz)
        ds.createVariable("x", "f8", ("x",))[:] = x
        ds.createVariable("y", "f8", ("y",))[:] = y
        ds.createVariable("z", "f8", ("z",))[:] = z
        mtype = ds.createCompoundType(st, "stress")
        data = ds.createVariable("data", mtype, ("z", "y", "x"))
        buf = np.empty((nz, ny, nx), dtype=st)
        for f in STRESS_FIELDS:
            buf[f] = np.transpose(comps[f], (2, 1, 0)).astype(np.float32)
        data[:] = buf
        for k, v in attrs.items():
            ds.setncattr(k, v)


def trilinear_sample(ncpath, qx, qy, qz):
    """Trilinear sample of the compound stress nc at points (qx,qy,qz);
    returns dict field -> values (NaN-clamped to grid edges)."""
    from netCDF4 import Dataset
    with Dataset(ncpath, "r") as ds:
        xg = np.asarray(ds.variables["x"][:], float)
        yg = np.asarray(ds.variables["y"][:], float)
        zg = np.asarray(ds.variables["z"][:], float)
        raw = ds.variables["data"][:]
        flds = {f: np.asarray(raw[f], float) for f in STRESS_FIELDS}  # (z,y,x)

    def locate(g, q):
        i = np.clip(np.searchsorted(g, q) - 1, 0, len(g) - 2)
        w = (q - g[i]) / (g[i + 1] - g[i])
        return i, np.clip(w, 0.0, 1.0)
    ix, wx = locate(xg, qx)
    iy, wy = locate(yg, qy)
    iz, wz = locate(zg, qz)
    out = {}
    for f, arr in flds.items():
        acc = np.zeros(len(qx))
        for dz_ in (0, 1):
            for dy_ in (0, 1):
                for dx_ in (0, 1):
                    w = ((wx if dx_ else 1 - wx) * (wy if dy_ else 1 - wy)
                         * (wz if dz_ else 1 - wz))
                    acc += w * arr[iz + dz_, iy + dy_, ix + dx_]
        out[f] = acc
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--csm-csv", required=True)
    ap.add_argument("--material-nc", required=True,
                    help="MUSCAL nc for the Sv_eff lithostatic profile")
    ap.add_argument("--mesh", required=True, help="PUML mesh for the self-check")
    ap.add_argument("--verify-vtu", required=True,
                    help="C1 projection fault VTU (ground truth for self-check)")
    ap.add_argument("--out", required=True, help="output ASAGI nc")
    ap.add_argument("--k-ratio", type=float, default=2.39, dest="k_ratio")
    ap.add_argument("--pp-model", default="hydrostatic")
    ap.add_argument("--strike-hint-az", type=float, default=314.0)
    ap.add_argument("--dx", type=float, default=1000.0, help="x,y spacing (m)")
    ap.add_argument("--dz", type=float, default=250.0, help="z spacing (m)")
    ap.add_argument("--xmin", type=float, default=350000.0)
    ap.add_argument("--xmax", type=float, default=630000.0)
    ap.add_argument("--ymin", type=float, default=3680000.0)
    ap.add_argument("--ymax", type=float, default=3850000.0)
    ap.add_argument("--zmin", type=float, default=-17000.0)
    ap.add_argument("--zmax", type=float, default=0.0)
    ap.add_argument("--selfcheck-tol", type=float, default=0.03,
                    help="max allowed p95 rel err vs the projection VTU")
    ap.add_argument("--overpressure-field", default=None,
                    help="optional .npz fault-local overpressure field "
                    "(fault_local_overpressure): subtract isotropic DeltaPp(s,z)"
                    " from the effective tensor on the grid. OFF by default "
                    "(nc then byte-identical to the no-overpressure build).")
    args = ap.parse_args()
    if args.k_ratio < 1.0:
        sys.exit(f"ERROR: --k-ratio must be >= 1 (got {args.k_ratio})")
    op_field = (load_overpressure_field(args.overpressure_field)
                if args.overpressure_field else None)

    # ---- CSM orientation (depth-invariant) ----
    lon, lat, _dep, S, _sh, _V, _R, _ap, _v2 = read_csm_csv(args.csm_csv)
    cx, cy = geographic_to_utm11n(lon, lat)

    # ---- grid ----
    gx = np.arange(args.xmin, args.xmax + 0.5 * args.dx, args.dx)
    gy = np.arange(args.ymin, args.ymax + 0.5 * args.dx, args.dx)
    gz = np.arange(args.zmin, args.zmax + 0.5 * args.dz, args.dz)
    nx, ny, nz = len(gx), len(gy), len(gz)
    print(f"grid: nx={nx} ny={ny} nz={nz} ({nx*ny*nz:,} nodes); "
          f"dx={args.dx} dz={args.dz} m")

    # ---- interpolate the CSM tensor at the grid (x,y), build axes + R ----
    GX, GY = np.meshgrid(gx, gy, indexing="ij")          # (nx, ny)
    qx2 = GX.ravel()
    qy2 = GY.ravel()
    S_xy, n_fb = interpolate_csm_field(cx, cy, S, qx2, qy2)   # (nxy, 6)
    if n_fb:
        print(f"warning: {n_fb}/{len(qx2)} (x,y) grid nodes outside the CSM "
              "hull -> nearest-neighbour fallback (grid margin)")
    U, _gaps, R_eig = csm_axes_and_shape(csm_tensors_tension(S_xy))
    R_xy = np.clip(R_eig, 0.0, 1.0)                      # (nxy,)
    nxy = len(R_xy)

    # ---- Sv_eff(z) from MUSCAL + hydrostatic P_p ----
    dgrid, svgrid, _rho = muscal_sv_total_profile(args.material_nc)

    # ---- assemble sigma_eff(x,y,z), comp-NEGATIVE Pa ----
    comps = {f: np.empty((nx, ny, nz), dtype=np.float64) for f in STRESS_FIELDS}
    cidx = {"s_xx": (0, 0), "s_yy": (1, 1), "s_zz": (2, 2),
            "s_xy": (0, 1), "s_yz": (1, 2), "s_xz": (0, 2)}
    for iz, zz in enumerate(gz):
        depth = max(-zz, 0.0)
        sv_total = sv_total_at(np.array([depth]), dgrid, svgrid)[0]
        pp = pore_pressure(np.array([depth]), args.pp_model)[0]
        sv_eff = sv_total - pp
        if sv_eff <= 0:                       # surface node: isotropic ~0
            for f in STRESS_FIELDS:
                comps[f][:, :, iz] = 0.0
            continue
        sv_eff_arr = np.full(nxy, sv_eff)
        sig1, sig2, sig3 = magnitudes_C1(sv_eff_arr, R_xy, args.k_ratio)
        sigma = build_tensor_from_axes(U, sig1, sig2, sig3)   # (nxy,3,3) MPa, comp+
        if op_field is not None:
            # fault-local overpressure: subtract isotropic DeltaPp(s,z) from the
            # EFFECTIVE tensor (raises mu_app at the gate; tau unchanged).  In
            # comp-NEGATIVE Pa this adds DeltaPp*1e6 to the diagonals only.
            dpp = delta_pp_mpa(op_field, qx2, qy2, np.full(nxy, zz))
            sigma = apply_overpressure_to_tensor(sigma, dpp)
        for f, (a, b) in cidx.items():
            # comp-NEGATIVE Pa for SeisSol
            comps[f][:, :, iz] = (-sigma[:, a, b] * 1.0e6).reshape(nx, ny)

    attrs = {
        "title": "SAFS V3 initial effective stress (CSM C1 k=2.39, comp-negative)",
        "convention": "SeisSol compression-NEGATIVE; effective (hydrostatic "
                      "P_p removed); global Cartesian x=E,y=N,z=Up; Pa",
        "source": "csm_stress_to_asagi.py from CSM YHSM-2013 orientation + "
                  "MUSCAL lithostat + hydrostatic P_p + closure C1",
        "k_ratio": float(args.k_ratio),
        "pp_model": args.pp_model,
        "crs": "EPSG:32611",
        "z_positive": "elevation",
        "units": "SI: stress Pa, axes m",
        "dz": float(args.dz),
    }
    if op_field is not None:
        m = op_field["meta"]
        attrs["overpressure"] = (
            f"fault-local graded overpressure: mu_s={m.get('recipe_mu_s')}, "
            f"floor={m.get('floor'):.3f}, lambda_max={m.get('lambda_max')}, "
            f"DPp_max={m.get('DPp_max_MPa'):.1f} MPa, gates="
            f"{m.get('gate_bands_km')}; field {os.path.basename(args.overpressure_field)}")
    # P-002 (REVIEW_nc): refuse to write a non-finite nc (a degenerate-R grid
    # node or NaN overpressure would otherwise ship NaN stress to SeisSol).
    n_nan = sum(int(np.sum(~np.isfinite(comps[f]))) for f in STRESS_FIELDS)
    if n_nan:
        sys.exit(f"ERROR: {n_nan} non-finite stress samples across "
                 f"{STRESS_FIELDS}; refusing to write a NaN nc (degenerate CSM "
                 "tensor or bad overpressure). Raise --dx or filter the input.")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    write_stress_asagi(args.out, gx, gy, gz, comps, attrs)
    print(f"wrote {args.out} ({os.path.getsize(args.out)/1e6:.1f} MB)")

    # ---- self-check vs the projection VTU ----
    geom, conn, bc = load_puml(args.mesh)
    sanity_check_face_ordering(geom, conn, bc)
    pts, tris = extract_fault(geom, conn, bc)
    cent, normals_raw, _areas = triangle_geometry(pts, tris)
    normals, _ng = harmonise_normals(normals_raw, args.strike_hint_az)
    strikes, dips, _degen = tandem_basis(normals)
    # P-001 (REVIEW_nc): trilinear_sample clamps out-of-grid facets to an edge
    # slice (plausible-but-wrong); abort instead so an under-covering grid fails
    # loudly rather than silently passing the median test.
    from netCDF4 import Dataset as _DS
    with _DS(args.out) as _ds:
        _xg = np.asarray(_ds["x"][:]); _yg = np.asarray(_ds["y"][:])
        _zg = np.asarray(_ds["z"][:])
    oob = ((cent[:, 0] < _xg.min()) | (cent[:, 0] > _xg.max())
           | (cent[:, 1] < _yg.min()) | (cent[:, 1] > _yg.max())
           | (cent[:, 2] < _zg.min()) | (cent[:, 2] > _zg.max()))
    if oob.any():
        sys.exit(f"ERROR: {int(oob.sum())} fault facets outside the nc grid box "
                 f"(deepest facet {cent[:, 2].min():.0f} m vs z-floor "
                 f"{_zg.min():.0f} m); widen --zmin/--xmin/... so the grid covers "
                 "the whole fault, or they read the ConstantMap fallback.")
    smp = trilinear_sample(args.out, cent[:, 0], cent[:, 1], cent[:, 2])
    sigma_f = np.zeros((len(tris), 3, 3))
    for f, (a, b) in cidx.items():
        sigma_f[:, a, b] = -smp[f] / 1.0e6          # back to comp+ MPa
        sigma_f[:, b, a] = sigma_f[:, a, b]
    # P_p is baked in (effective); resolve with P_p=0 so sigma_n_eff = n.sigma.n
    r = resolve_tractions(sigma_f, strikes, dips, normals, 0.0)
    truth = read_vtu_cell(args.verify_vtu,
                          ["sigma_n_eff_MPa_cell", "tau_magnitude_MPa_cell"])
    fin = (np.isfinite(r["sigma_n_eff"]) & np.isfinite(r["tau_magnitude"])
           & (truth["sigma_n_eff_MPa_cell"] > 1.0))
    relsn = np.abs(r["sigma_n_eff"][fin] - truth["sigma_n_eff_MPa_cell"][fin]
                   ) / truth["sigma_n_eff_MPa_cell"][fin]
    reltau = np.abs(r["tau_magnitude"][fin] - truth["tau_magnitude_MPa_cell"][fin]
                    ) / np.maximum(truth["tau_magnitude_MPa_cell"][fin], 1.0)
    print("self-check vs projection VTU (per-facet, %d facets):" % int(fin.sum()))
    print("  sigma_n_eff rel err: median %.2e  p95 %.2e  max %.2e"
          % (np.median(relsn), np.percentile(relsn, 95), relsn.max()))
    print("  tau         rel err: median %.2e  p95 %.2e  max %.2e"
          % (np.median(reltau), np.percentile(reltau, 95), reltau.max()))
    # P-005 (REVIEW_nc): gate on the p95, not the median, so a localized wrong
    # region (a coherent patch a few % off) fails instead of being masked.
    ok = (np.percentile(relsn, 95) < args.selfcheck_tol
          and np.percentile(reltau, 95) < args.selfcheck_tol)
    print("SELF-CHECK:", "PASS" if ok else
          "FAIL (p95 rel err exceeds tol; raise grid resolution)")
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
