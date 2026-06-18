#!/usr/bin/env python3
"""Project a CSM (Community Stress Model) orientation field onto a PUML
mesh fault and write a ParaView VTU (fault tractions + mu) + summary JSON.

Sibling of on_fault_stress_projection/project_stress_to_vtu.py (same
Tandem fault basis, same normal harmonisation, same output field names)
but the stress tensor is SPATIALLY VARYING, built from the CSM
"YHSM-2013" model of Yang & Hauksson (CSV export by K. Luttrell):

    input stress : CSM csv with per-(lon,lat) tensor components
                   See,Sen,Seu,Snn,Snu,Suu in the (east,north,up) frame,
                   TENSION-POSITIVE, MPa.  The model is ORIENTATION
                   ONLY -- magnitudes are NOT meaningful (deviatoric,
                   non-uniformly normalised), and depth-invariant.
    input mesh   : pumgen .puml.h5 (fault = boundary-code-3 faces)

Reconstruction (the whole point of this tool): per CSM point, take the
principal AXES of the tension-positive tensor (eigendecomposition) and
prescribe effective principal magnitudes from the CLI:

    sigma_eff(x) = s1 * u_c u_c^T + s2 * u_i u_i^T + s3 * u_t u_t^T

with u_c / u_i / u_t the most-compressional / intermediate / most-
tensional CSM axes and s1 >= s2 >= s3 the compression-positive
effective principal stresses (MPa):

    --sigma1-eff-MPa  (default 80)   most compressive
    --sigma2-eff-MPa  (default 40)   intermediate
    --sigma3-eff-MPa  (default 35)   least compressive

Defaults reproduce the previous CONSTANT-tensor run
(on_fault_stress_projection: eff tensor [[46.107,14.387,0],
[14.387,73.893,0],[0,0,35]] => principals 80/40/35, P_p 20 MPa);
adjust them later as the magnitude model firms up.

Spatial sampling: CSM (lon,lat) -> UTM 11N via pyproj EPSG:4326 ->
EPSG:32611 always_xy (the velocity-pipeline convention, see
velocity/code/crs.py and generate_velocity_nc_from_raw); per-facet
nearest-neighbour in (x,y) -- the model is depth-invariant so z is
ignored.  UTM grid convergence (<~1 deg here) is neglected, i.e. the
mesh (x,y,z) axes are identified with (east,north,up), consistent with
the rest of the project.

Outputs (into --out-dir):
    <prefix>_fault_stress.vtu  per-facet + per-vertex resolved tractions,
                               mu, sampled eff tensor components, CSM
                               SHmax + sampling diagnostics
    <prefix>_summary.json      params, per-field stats, orientation
                               validation (eigvec vs csv V-columns,
                               SHmax reconstruction), hypocenter check

Conventions (matching the reference pipeline / summary "convention" key):
    compression POSITIVE (SEAS internal) in all outputs, units MPa;
    sigma_n_eff = sigma_n_total - P_p (magnitudes are EFFECTIVE);
    fault basis: Tandem  s = up x n,  d = s x n (down-dip);
    normals harmonised to the SW half-space via --strike-hint-az
    (right-lateral SAF) =>  tau_strike + = right-lateral,
                            tau_dip    + = reverse;
    frame (east, north, up), UTM Zone 11N, metres.

Usage:
    python3 project_csm_stress_to_vtu.py <csm.csv> <mesh.puml.h5>
        [--out-dir DIR] [--prefix NAME]
        [--sigma1-eff-MPa 80] [--sigma2-eff-MPa 40] [--sigma3-eff-MPa 35]
        [--P-p-MPa 20] [--strike-hint-az 314] [--hypocenter X Y Z]
        [--max-nn-dist-km 5]

Requires numpy + h5py + scipy (cKDTree) + pyproj.
"""
import argparse
import csv as csv_mod
import json
import os
import struct
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fault_local_overpressure import (  # noqa: E402
    load_overpressure_field, delta_pp_mpa, apply_overpressure_to_tensor)

EPS = 1.0e-9
UP = np.array([0.0, 0.0, 1.0])
G_ACCEL = 9.81                 # m/s^2 (matches the plan's Sv/P_p numbers)
RHO_W = 1000.0                 # kg/m^3, water (hydrostatic pore pressure)

# PUML Numbering<TETRAHEDRON> face->local-vertex map (see sibling tool;
# the free-surface flatness check below re-validates it on every run).
FACE_VERTS = ((1, 0, 2), (0, 1, 3), (1, 2, 3), (2, 0, 3))
BC_FAULT = 3

# CSM csv column indices (header line "# LON,LAT,DEP,See,...")
COL_LON, COL_LAT, COL_DEP = 0, 1, 2
COL_S = slice(3, 9)            # See, Sen, Seu, Snn, Snu, Suu
COL_SHMAX = 9
COL_PHI = 11                  # phi shape parameter (Angelier 1979), tension+
COL_R = 12                    # R stress ratio (Gephart & Forsyth 1984), tens+
COL_APHI = 13                 # Aphi (Simpson 1997)
COL_V = slice(20, 29)          # V1x..V1z, V2x..V2z, V3x..V3z
COL_V2_PLUNGE = 30            # V2 plunge of intermediate axis (deg below horiz)

# Luttrell differential-stress csv: column 16 (1-based) = "differential
# stress (S1-S3) in MPa".  Same lon/lat columns as the YHSM csv.
COL_DIFF = 15


# ------------------------------------------------------------- CSM csv
def read_csm_csv(path):
    """Parse the CSM csv: lon/lat/dep, tensor components, SHmax, the
    author principal-axis vectors (for validation), and the shape /
    regime scalars R (col 13), Aphi (col 14) and the intermediate-axis
    plunge (col 31) used by the depth-dependent magnitude model.
    '#' lines skipped, empty numeric fields -> NaN."""
    lon, lat, dep, S, shmax, V = [], [], [], [], [], []
    R_gf, aphi, v2pl = [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            row = next(csv_mod.reader([line]))
            vals = [float(v) if v.strip() not in ("", "NaN", "nan")
                    else np.nan for v in row]
            lon.append(vals[COL_LON])
            lat.append(vals[COL_LAT])
            dep.append(vals[COL_DEP])
            S.append(vals[COL_S])
            shmax.append(vals[COL_SHMAX])
            V.append(vals[COL_V])
            R_gf.append(vals[COL_R])
            aphi.append(vals[COL_APHI])
            v2pl.append(vals[COL_V2_PLUNGE])
    lon, lat, dep = map(np.asarray, (lon, lat, dep))
    S = np.asarray(S)
    if len(S) == 0:
        sys.exit(f"ERROR: no data rows in {path}")
    bad = np.isnan(S).any(axis=1)
    if bad.any():
        print(f"warning: dropping {int(bad.sum())} CSM rows with NaN "
              "tensor components", file=sys.stderr)
    keep = ~bad
    return (lon[keep], lat[keep], dep[keep], S[keep],
            np.asarray(shmax)[keep], np.asarray(V)[keep].reshape(-1, 3, 3),
            np.asarray(R_gf)[keep], np.asarray(aphi)[keep],
            np.asarray(v2pl)[keep])


def csm_tensors_tension(S):
    """(N,6) See,Sen,Seu,Snn,Snu,Suu -> (N,3,3) symmetric tensors in the
    (e,n,u) frame, TENSION-POSITIVE (as in the csv)."""
    T = np.empty((len(S), 3, 3))
    T[:, 0, 0] = S[:, 0]
    T[:, 0, 1] = T[:, 1, 0] = S[:, 1]
    T[:, 0, 2] = T[:, 2, 0] = S[:, 2]
    T[:, 1, 1] = S[:, 3]
    T[:, 1, 2] = T[:, 2, 1] = S[:, 4]
    T[:, 2, 2] = S[:, 5]
    return T


def reconstruct_eff_tensors(T_tension, s1, s2, s3):
    """Keep CSM principal AXES, prescribe effective magnitudes.

    eigh(ascending) of the tension-positive tensor orders the axes
    most-compressional (index 0, csv V3) -> most-tensional (index 2,
    csv V1).  Assign s1 (largest compression) to index 0, s2 to 1,
    s3 to 2 and rebuild a compression-positive SEAS tensor.

    Returns (sigma_eff (N,3,3), eigvecs (N,3,3) columns ascending,
    eig_gap_min (N,) smallest eigenvalue gap of the source tensor --
    near-zero gap means the corresponding axes are poorly constrained).
    """
    w, U = np.linalg.eigh(T_tension)          # ascending eigenvalues
    gaps = np.minimum(w[:, 1] - w[:, 0], w[:, 2] - w[:, 1])
    mags = np.array([s1, s2, s3])
    sigma = np.einsum("k,nik,njk->nij", mags, U, U)
    return sigma, U, gaps


def csm_axes_and_shape(T_tension):
    """Eigendecomposition of the tension-positive CSM tensor for the
    depth-dependent magnitude model.  Returns (U, gaps, R_eig):
      U      (N,3,3) eigenvectors, columns ASCENDING eigenvalue ->
             [most-compressional, intermediate, most-tensional];
      gaps   (N,) smallest eigenvalue gap (axis-constraint quality);
      R_eig  (N,) shape ratio (w1-w0)/(w2-w0) == csv column 13
             (verified to csv rounding); NaN where the gap collapses."""
    w, U = np.linalg.eigh(T_tension)
    gaps = np.minimum(w[:, 1] - w[:, 0], w[:, 2] - w[:, 1])
    denom = w[:, 2] - w[:, 0]
    R_eig = np.where(denom > EPS, (w[:, 1] - w[:, 0]) / np.where(
        denom > EPS, denom, 1.0), np.nan)
    return U, gaps, R_eig


# ---------------------------------------------- depth-dependent magnitudes
def muscal_sv_total_profile(nc_path, g=G_ACCEL):
    """Lithostatic total overburden from the MUSCAL material file.

    Reads the bulk density rho(z,y,x) (compound variable `data.rho`,
    kg/m^3), forms the LATERAL-MEAN profile rho_mean(z) (plan Section
    2.1: not per-column, to avoid non-equilibrium lateral Sv jumps),
    and integrates Sv_total(depth) = int_0^depth rho_mean g d(depth).

    Returns (depth_m ascending, Sv_total_MPa, rho_mean_kgm3) on the nc
    z-grid (depth = -elevation, surface first)."""
    try:
        from netCDF4 import Dataset
    except ImportError:
        sys.exit("ERROR: netCDF4 is required for --density-source muscal "
                 "(conda env 'pythonenv')")
    if not os.path.exists(nc_path):
        sys.exit(f"ERROR: MUSCAL material file not found: {nc_path}")
    ds = Dataset(nc_path, "r")
    if "z" not in ds.variables or "data" not in ds.variables:
        ds.close()
        sys.exit(f"ERROR: {nc_path} missing 'z'/'data' variables")
    z = np.asarray(ds.variables["z"][:], dtype=float)          # elevation, m
    rho = np.asarray(ds.variables["data"][:]["rho"], dtype=float)  # (z,y,x)
    ds.close()
    if rho.shape[0] != len(z):
        sys.exit("ERROR: MUSCAL rho first axis does not match z")
    rho_mean = np.nanmean(rho.reshape(len(z), -1), axis=1)
    order = np.argsort(-z)                  # z descending: 0, -250, ...
    zt, rt = z[order], rho_mean[order]
    depth = -zt                             # ascending 0, 250, ... (m)
    if depth[0] > 1.0e-6:
        sys.exit("ERROR: MUSCAL z grid does not reach the surface (z=0)")
    sv = np.concatenate(
        [[0.0], np.cumsum(0.5 * (rt[1:] + rt[:-1]) * g * np.diff(depth))])
    return depth, sv / 1.0e6, rt            # m, MPa, kg/m^3 (depth-ordered)


def sv_total_at(depth_m, depth_grid, sv_grid):
    """Linear interpolation of Sv_total(MPa) at facet depths (m)."""
    return np.interp(np.clip(depth_m, depth_grid[0], depth_grid[-1]),
                     depth_grid, sv_grid)


def pore_pressure(depth_m, model, g=G_ACCEL, rho_w=RHO_W):
    """Pore-pressure model P_p(depth) in MPa (depth >= 0, m)."""
    if model == "hydrostatic":
        return rho_w * g * np.maximum(depth_m, 0.0) / 1.0e6
    sys.exit(f"ERROR: unknown --pp-model {model!r} (only 'hydrostatic')")


def read_diff_csv(path):
    """Differential-stress csv (Luttrell-2017 CSM export): lon, lat,
    D = S1-S3 (column 16).  '#' lines skipped; rows with empty/NaN D or
    coordinates dropped."""
    lon, lat, dvals = [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            row = next(csv_mod.reader([line]))
            try:
                lo = float(row[COL_LON])
                la = float(row[COL_LAT])
                ds = row[COL_DIFF].strip()
                d = float(ds) if ds not in ("", "NaN", "nan") else np.nan
            except (ValueError, IndexError):
                continue
            lon.append(lo)
            lat.append(la)
            dvals.append(d)
    lon, lat, dvals = map(np.asarray, (lon, lat, dvals))
    good = np.isfinite(dvals) & np.isfinite(lon) & np.isfinite(lat)
    if not good.any():
        sys.exit(f"ERROR: no finite differential-stress rows in {path}")
    return lon[good], lat[good], dvals[good]


def magnitudes_C2(Sv_eff, R, D, cap_k_max):
    """Closure C2 (plan Section 4): sig2 = Sv_eff, sig1 = Sv_eff + R*D,
    sig3 = Sv_eff - (1-R)*D.

    FRICTIONAL cap: the brittle crust cannot sustain a differential
    stress larger than the Mohr-Coulomb / Byerlee frictional strength on
    optimally oriented faults, i.e. the stress ratio sig1/sig3 is bounded
    by k = (sqrt(mu^2+1)+mu)^2 (Byerlee 1978; Sibson 1974; observed in
    the SAF region by Zoback & Healy 1992 at Cajon Pass and generalised
    by Townend & Zoback 2000).  That ceiling, expressed as a maximum
    differential stress, is the plan's closure-C1 value at k = cap_k_max:

        D_max(z) = (k-1) * Sv_eff(z) / ((1-R)*k + R)

    Because D_max is PROPORTIONAL to Sv_eff, the sustainable differential
    stress -> 0 as the effective overburden -> 0 at the free surface
    (Brace & Kohlstedt 1980: Byerlee's law as a bound on in-situ stress
    under hydrostatic pore pressure).  D is capped to D_max so the stress
    ratio never exceeds cap_k_max and the near-surface stress is no
    longer super-critical.  Guarantees 0 < sig3 <= sig2 for k >= 1,
    R in [0, 1].  Returns (sig1, sig2, sig3, D_applied, capped)."""
    sig2 = Sv_eff
    D_max = (cap_k_max - 1.0) * Sv_eff / ((1.0 - R) * cap_k_max + R)
    capped = D > D_max
    D_applied = np.where(capped, D_max, D)
    sig1 = Sv_eff + R * D_applied
    sig3 = Sv_eff - (1.0 - R) * D_applied
    return sig1, sig2, sig3, D_applied, capped


def magnitudes_C1(Sv_eff, R, k):
    """Closure C1 (plan Section 4): sig1 = k*sig3, sig2 = Sv_eff,
    sig3 = Sv_eff/((1-R)*k + R).  Automatically ordered for k>=1,
    R in [0,1] (no cap needed)."""
    sig2 = Sv_eff
    sig3 = Sv_eff / ((1.0 - R) * k + R)
    sig1 = k * sig3
    return sig1, sig2, sig3


def build_tensor_from_axes(U, sig1, sig2, sig3):
    """Option A (--axes csm): keep the CSM principal AXES, assign the
    effective magnitudes.  U columns ascending = [most-compressional,
    intermediate, most-tensional]; assign sig1->col0, sig2->col1,
    sig3->col2 (compression-positive)."""
    mags = np.stack([sig1, sig2, sig3], axis=1)        # (N,3)
    return np.einsum("nk,nik,njk->nij", mags, U, U)


def build_tensor_andersonian(az_deg, sig1, sig2, sig3):
    """Option B (--axes andersonian): one principal axis exactly
    vertical.  e_H = (sin az, cos az, 0) carries sig1, e_h = (cos az,
    -sin az, 0) carries sig3, e_v = (0,0,1) carries sig2 (Anderson
    1951; az = SHmax deg E of N)."""
    az = np.radians(az_deg)
    z0 = np.zeros_like(az)
    eH = np.stack([np.sin(az), np.cos(az), z0], axis=1)
    eh = np.stack([np.cos(az), -np.sin(az), z0], axis=1)
    ev = np.broadcast_to(UP, eH.shape)
    return (sig1[:, None, None] * np.einsum("ni,nj->nij", eH, eH)
            + sig3[:, None, None] * np.einsum("ni,nj->nij", eh, eh)
            + sig2[:, None, None] * np.einsum("ni,nj->nij", ev, ev))


def apply_deep_taper(sigma, depth_m, zseis_km):
    """Optional deep taper (plan Section 2.3, default OFF): below the
    seismogenic depth taper the DEVIATORIC part of each facet tensor to
    zero while keeping the isotropic (mean) effective part:
        sigma_t = Omega(z)*(sigma - p I) + p I,  p = trace(sigma)/3.
    Omega = 1 for depth <= z_seis, cosine ramp to 0 over a 0.3*z_seis
    width below it.  (The plan leaves the isotropic continuation
    'Sv_iso_eff' unspecified in form; the built-tensor mean p is the
    natural choice and is what continues below z_seis.)"""
    z_km = depth_m / 1000.0
    width = 0.3 * zseis_km
    t = np.clip((z_km - zseis_km) / max(width, EPS), 0.0, 1.0)
    omega = 0.5 * (1.0 + np.cos(np.pi * t))            # 1 above, 0 below
    p = np.einsum("nii->n", sigma) / 3.0
    iso = p[:, None, None] * np.eye(3)[None, :, :]
    dev = sigma - iso
    return omega[:, None, None] * dev + iso


def shape_recovery_check(sigma):
    """R_GF recovered from the BUILT compression-positive tensor:
    eigvalsh ascending = [sig3, sig2, sig1]; R = (sig1-sig2)/(sig1-sig3).
    Should reproduce the input R to ~1e-3 (closed form -> failure = bug)."""
    w = np.linalg.eigvalsh(sigma)
    denom = w[:, 2] - w[:, 0]
    return np.where(denom > EPS, (w[:, 2] - w[:, 1]) / np.where(
        denom > EPS, denom, 1.0), np.nan)


def interpolate_csm_field(cx, cy, values, qx, qy):
    """Barycentric (Delaunay-linear) interpolation of per-CSM-point
    `values` (shape (N,) or (N, k)) at query points (qx, qy), used by
    `--sample linear` to smooth the otherwise piecewise-constant
    (nearest-neighbour) CSM orientation field.  LinearNDInterpolator
    returns NaN outside the convex hull of the CSM points; those query
    points fall back to nearest neighbour.  Returns (interp, n_fallback)
    with interp reshaped to (len(q),) + values.shape[1:]."""
    try:
        from scipy.interpolate import (LinearNDInterpolator,
                                       NearestNDInterpolator)
    except ImportError:
        sys.exit("ERROR: scipy.interpolate is required for --sample linear "
                 "(conda env 'pythonenv')")
    pts = np.stack([cx, cy], axis=1)
    q = np.stack([qx, qy], axis=1)
    vals = np.asarray(values, dtype=float)
    flat = vals.reshape(len(vals), -1)
    out = LinearNDInterpolator(pts, flat)(q)        # (Nq, k); NaN outside hull
    bad = ~np.isfinite(out).all(axis=1)
    n_fallback = int(bad.sum())
    if n_fallback:
        out[bad] = NearestNDInterpolator(pts, flat)(q[bad])
    return out.reshape((len(q),) + vals.shape[1:]), n_fallback


def shmax_from_sigma(sigma):
    """Azimuth (deg E of N, in [0,180)) of the most-compressive
    eigenvector of the horizontal 2x2 block of the compression-positive
    tensor (the standard CSM-style SHmax proxy)."""
    h = sigma[:, :2, :2]
    w, U = np.linalg.eigh(h)                  # ascending
    v = U[:, :, 1]                            # most compressive (e, n)
    az = np.degrees(np.arctan2(v[:, 0], v[:, 1]))   # E of N
    return np.mod(az, 180.0)


def angdiff180(a, b):
    """Minimal angular distance between bi-directional azimuths (deg)."""
    d = np.abs(np.mod(a - b, 180.0))
    return np.minimum(d, 180.0 - d)


# ----------------------------------------------------- lon/lat -> UTM
def geographic_to_utm11n(lon, lat):
    """EPSG:4326 -> EPSG:32611, always_xy — the velocity-pipeline
    convention (velocity/code/crs.py / generate_velocity_nc_from_raw)."""
    try:
        from pyproj import Transformer
    except ImportError:
        sys.exit("ERROR: pyproj is required (conda env 'pythonenv')")
    tr = Transformer.from_crs("EPSG:4326", "EPSG:32611", always_xy=True)
    x, y = tr.transform(np.asarray(lon).ravel(), np.asarray(lat).ravel())
    return (np.asarray(x).reshape(np.shape(lon)),
            np.asarray(y).reshape(np.shape(lat)))


# ------------------------------------------------------------ PUML mesh
def load_puml(path):
    with h5py.File(path, "r") as f:
        geom = f["geometry"][:]
        conn = f["connect"][:].astype(np.int64)
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
    return geom, conn, bc


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
    """Project the PER-FACET effective SEAS tensors onto each facet."""
    t = np.einsum("kij,kj->ki", sigma_eff, normals)    # traction (N,3) MPa
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
    if not finite.any():
        return {"min": None, "median": None, "max": None,
                "nan_count": int((~finite).sum())}
    return {"min": float(np.min(arr[finite])),
            "median": float(np.median(arr[finite])),
            "max": float(np.max(arr[finite])),
            "nan_count": int((~finite).sum())}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("csm_csv", help="CSM csv (YHSM-2013 export: "
                    "tension-positive (e,n,u) tensor, orientation only)")
    ap.add_argument("mesh", help="pumgen .puml.h5 mesh (UTM 11N, m)")
    ap.add_argument("--out-dir", default=".")
    ap.add_argument("--prefix", default=None,
                    help="output prefix (default <csv>_on_<mesh>)")
    ap.add_argument("--sigma1-eff-MPa", type=float, default=80.0,
                    dest="s1", help="most compressive effective principal "
                    "stress, compression-positive MPa (default %(default)s "
                    "= previous constant-tensor run)")
    ap.add_argument("--sigma2-eff-MPa", type=float, default=40.0,
                    dest="s2", help="intermediate effective principal "
                    "stress (default %(default)s)")
    ap.add_argument("--sigma3-eff-MPa", type=float, default=35.0,
                    dest="s3", help="least compressive effective principal "
                    "stress (default %(default)s)")
    ap.add_argument("--P-p-MPa", type=float, default=20.0, dest="P_p",
                    help="pore pressure; magnitudes are effective, "
                    "sigma_n_total = sigma_n_eff + P_p (default %(default)s)")
    ap.add_argument("--strike-hint-az", type=float, default=314.0,
                    help="global strike azimuth (deg) for normal "
                    "harmonisation; 314 = SAF N46W right-lateral")
    ap.add_argument("--hypocenter", type=float, nargs=3, default=None,
                    metavar=("X", "Y", "Z"),
                    help="report resolved tractions at the facet nearest "
                    "this point")
    ap.add_argument("--max-nn-dist-km", type=float, default=5.0,
                    help="warn if any facet is farther than this from the "
                    "nearest CSM sample (default %(default)s km)")

    # ---- depth-dependent magnitude model (plan rev 5) -----------------
    ap.add_argument("--magnitude-mode", choices=("constant", "depth-R"),
                    default="constant",
                    help="constant: prescribed sigma1/2/3 on the CSM axes "
                    "(default, reproduces the existing artifact). depth-R: "
                    "depth-dependent Sv_eff(z) + closure + per-point R.")
    ap.add_argument("--closure", choices=("diff", "ratio"), default="diff",
                    help="depth-R closure between S1 and S3: diff = "
                    "differential stress D=S1-S3 (DEFAULT, Luttrell C2-s); "
                    "ratio = effective ratio k=S1/S3 (C1).")
    ap.add_argument("--diff-csv", default=None,
                    help="differential-stress csv (Luttrell-2017 CSM "
                    "export; D=S1-S3 in column 16) for closure diff. If "
                    "omitted, the scalar --diff-MPa is used everywhere.")
    ap.add_argument("--diff-MPa", type=float, default=40.0, dest="diff_mpa",
                    help="scalar differential stress D when --diff-csv is "
                    "not given (default %(default)s)")
    ap.add_argument("--k-ratio", type=float, default=3.12, dest="k_ratio",
                    help="effective stress ratio k=S1/S3 for --closure "
                    "ratio (default %(default)s = mu 0.6 frictional)")
    ap.add_argument("--pp-model", choices=("hydrostatic",),
                    default="hydrostatic",
                    help="depth-R pore-pressure model (default hydrostatic)")
    ap.add_argument("--density-source", choices=("muscal", "constant"),
                    default="muscal",
                    help="depth-R overburden density: muscal lateral-mean "
                    "rho(z) (default) or constant --rho-const-kgm3")
    ap.add_argument("--material-nc",
                    default="safs_seisol_v2_0_0_RSSRW/safs_material_cvm.nc",
                    help="MUSCAL material .nc for --density-source muscal "
                    "(default %(default)s)")
    ap.add_argument("--rho-const-kgm3", type=float, default=2670.0,
                    dest="rho_const",
                    help="constant overburden density (default %(default)s)")
    ap.add_argument("--axes", choices=("csm", "andersonian"), default="csm",
                    help="depth-R tensor axes: csm = full CSM eigenvectors "
                    "(default, Option A); andersonian = SHmax azimuth + "
                    "vertical (Option B)")
    ap.add_argument("--sample", choices=("nearest", "linear"),
                    default="nearest",
                    help="depth-R CSM orientation sampling: nearest "
                    "(default; piecewise-constant in the ~2 km CSM grid "
                    "cells -> blocky) or linear (barycentric interpolation "
                    "of the CSM stress tensor onto facet centroids -> "
                    "smooth orientation field, NN fallback outside the data "
                    "hull). Affects depth-R only.")
    ap.add_argument("--cap-k-max", type=float, default=5.83, dest="cap_k_max",
                    help="closure diff frictional cap: limit D so the stress "
                    "ratio sig1/sig3 <= this Byerlee ceiling (default "
                    "%(default)s = mu 1.0; use 3.12 for mu 0.6). Bounds the "
                    "near-surface apparent friction and tapers the "
                    "differential to 0 as Sv_eff -> 0 at the free surface "
                    "(Byerlee 1978; Brace & Kohlstedt 1980; Zoback & Healy "
                    "1992; Townend & Zoback 2000).")
    ap.add_argument("--taper-zseis-km", type=float, default=None,
                    help="optional deep taper: taper the deviatoric part to "
                    "zero below this seismogenic depth (km); off by default")
    ap.add_argument("--overpressure-field", default=None,
                    help="optional .npz fault-local overpressure field "
                    "(fault_local_overpressure.build_overpressure_field): "
                    "subtract isotropic DeltaPp(s,z) from the effective tensor "
                    "to raise mu_app at a gate. OFF by default (output then "
                    "byte-identical to the no-overpressure run).")
    args = ap.parse_args()

    if args.magnitude_mode == "constant" and not (
            args.s1 >= args.s2 >= args.s3 > 0):
        sys.exit("ERROR: require sigma1 >= sigma2 >= sigma3 > 0 "
                 f"(got {args.s1}, {args.s2}, {args.s3}); the magnitudes "
                 "are assigned to the most->least compressional CSM axes "
                 "in that order.")
    if args.magnitude_mode == "depth-R":
        if args.closure == "ratio" and args.k_ratio < 1.0:
            sys.exit(f"ERROR: --k-ratio must be >= 1 (got {args.k_ratio})")
        if args.cap_k_max < 1.0:
            sys.exit(f"ERROR: --cap-k-max must be >= 1 (got {args.cap_k_max})")
        if args.taper_zseis_km is not None and args.taper_zseis_km <= 0:
            sys.exit("ERROR: --taper-zseis-km must be > 0")

    os.makedirs(args.out_dir, exist_ok=True)
    op_field = (load_overpressure_field(args.overpressure_field)
                if args.overpressure_field else None)
    cstem = os.path.splitext(os.path.basename(args.csm_csv))[0]
    mstem = os.path.basename(args.mesh)
    for suf in (".h5", ".puml"):
        if mstem.endswith(suf):
            mstem = mstem[:-len(suf)]
    prefix = args.prefix or f"{cstem}_on_{mstem}"

    # ---- CSM model: axes from data, magnitudes prescribed
    (lon, lat, dep, S, shmax_model, V_model,
     R_csv, aphi_model, v2pl_model) = read_csm_csv(args.csm_csv)
    print(f"CSM: {len(lon)} samples, lon [{lon.min():.3f}, {lon.max():.3f}]"
          f", lat [{lat.min():.3f}, {lat.max():.3f}], "
          f"depths {np.unique(dep)} km (model is depth-invariant)")
    T_tension = csm_tensors_tension(S)
    if args.magnitude_mode == "constant":
        sigma_all, U, eig_gaps = reconstruct_eff_tensors(
            T_tension, args.s1, args.s2, args.s3)
        R_eig = None
    else:
        sigma_all = None
        U, eig_gaps, R_eig = csm_axes_and_shape(T_tension)

    # orientation validation 1: our eigenvectors vs the csv V1/V3 columns
    # (V rows: V1 = most tensional = our column 2; V3 = most compr = col 0)
    val = {}
    vok = np.isfinite(V_model.reshape(len(V_model), -1)).all(axis=1)
    if vok.any():
        d3 = np.abs(np.einsum("ni,ni->n", V_model[vok, 2], U[vok, :, 0]))
        d1 = np.abs(np.einsum("ni,ni->n", V_model[vok, 0], U[vok, :, 2]))
        val["eigvec_dot_V3_compressional"] = field_stats(d3)
        val["eigvec_dot_V1_tensional"] = field_stats(d1)
        print(f"axis check vs csv eigenvectors: |dot| median "
              f"V3 {np.median(d3):.4f}, V1 {np.median(d1):.4f}")

    # ---- mesh + fault
    geom, conn, bc = load_puml(args.mesh)
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

    # ---- nearest-neighbour sampling in UTM (2D; model depth-invariant)
    from scipy.spatial import cKDTree
    cx, cy = geographic_to_utm11n(lon, lat)
    tree = cKDTree(np.stack([cx, cy], axis=1))
    nn_dist, nn_idx = tree.query(centroids[:, :2])
    nn_dist_km = nn_dist / 1000.0
    n_far = int((nn_dist_km > args.max_nn_dist_km).sum())
    if n_far:
        print(f"warning: {n_far} facets farther than "
              f"{args.max_nn_dist_km} km from the nearest CSM sample "
              f"(max {nn_dist_km.max():.2f} km)", file=sys.stderr)
    print(f"sampling: nn distance min {nn_dist_km.min():.3f} / median "
          f"{np.median(nn_dist_km):.3f} / max {nn_dist_km.max():.3f} km")

    shmax_facet = shmax_model[nn_idx]
    eiggap_facet = eig_gaps[nn_idx]

    # ---- magnitudes -> per-facet effective tensor + pore pressure
    depth_cell = {}          # extra cell fields (depth-R mode only)
    mag_params = {}          # extra summary params (depth-R mode only)
    if args.magnitude_mode == "constant":
        sigma_facet = sigma_all[nn_idx]            # (Nf, 3, 3) eff, MPa
        P_p_facet = args.P_p
        # SHmax validation on the per-CSM-point reconstructed tensor (as before)
        shmax_recon = shmax_from_sigma(sigma_all)
        sok = np.isfinite(shmax_model)
        dsh = angdiff180(shmax_recon[sok], shmax_model[sok])
    else:
        depth_facet = np.maximum(-centroids[:, 2], 0.0)        # m, >= 0
        if args.density_source == "muscal":
            dgrid, svgrid, _rho = muscal_sv_total_profile(args.material_nc)
            Sv_total = sv_total_at(depth_facet, dgrid, svgrid)
            dens_desc = f"muscal lateral-mean rho(z) from {args.material_nc}"
        else:
            Sv_total = args.rho_const * G_ACCEL * depth_facet / 1.0e6
            dens_desc = f"constant rho = {args.rho_const} kg/m^3"
        P_p_facet = pore_pressure(depth_facet, args.pp_model)
        Sv_eff = Sv_total - P_p_facet
        if np.any(Sv_eff <= 0):
            sys.exit(f"ERROR: Sv_eff <= 0 on {int(np.sum(Sv_eff <= 0))} "
                     "facets (pore pressure exceeds overburden)")
        # ---- sample the CSM orientation + shape ratio onto the facets ----
        # per-point R: csv col 13 where finite, else eigenvalue R; clip [0,1]
        R_point = np.clip(np.where(np.isfinite(R_csv), R_csv, R_eig), 0.0, 1.0)
        if args.sample == "linear":
            # interpolate the 6 tension-positive tensor components, then
            # rebuild axes + shape ratio per facet (smooth orientation)
            S_facet, n_interp_fb = interpolate_csm_field(
                cx, cy, S, centroids[:, 0], centroids[:, 1])
            T_facet = csm_tensors_tension(S_facet)
            w_facet, U_facet = np.linalg.eigh(T_facet)
            denom = w_facet[:, 2] - w_facet[:, 0]
            R_facet = (w_facet[:, 1] - w_facet[:, 0]) / np.where(
                denom > EPS, denom, np.nan)
            # SHmax + eig-gap from the interpolated tensor (override the NN
            # values set above so the cell diagnostics match the smoothed
            # field); SHmax wants compression-positive -> negate
            shmax_facet = shmax_from_sigma(-T_facet)
            eiggap_facet = np.minimum(w_facet[:, 1] - w_facet[:, 0],
                                      w_facet[:, 2] - w_facet[:, 1])
            aphi_facet, _ = interpolate_csm_field(
                cx, cy, aphi_model, centroids[:, 0], centroids[:, 1])
            v2pl_facet, _ = interpolate_csm_field(
                cx, cy, v2pl_model, centroids[:, 0], centroids[:, 1])
            if n_interp_fb:
                print(f"warning: {n_interp_fb} facets outside the CSM convex "
                      "hull -> nearest-neighbour fallback", file=sys.stderr)
            print(f"orientation: linear interpolation of the CSM tensor "
                  f"({n_interp_fb} NN-fallback facets)")
        else:                                  # nearest (default, blocky)
            U_facet = U[nn_idx]
            R_facet = R_point[nn_idx]
            aphi_facet = aphi_model[nn_idx]
            v2pl_facet = v2pl_model[nn_idx]
            n_interp_fb = 0
        R_facet = np.clip(R_facet, 0.0, 1.0)
        # closure between sig1 and sig3
        if args.closure == "diff":
            if args.diff_csv:
                dlon, dlat, dval = read_diff_csv(args.diff_csv)
                dx, dy = geographic_to_utm11n(dlon, dlat)
                dtree = cKDTree(np.stack([dx, dy], axis=1))
                dd, didx = dtree.query(centroids[:, :2])
                D_diff = dval[didx]
                d_nn_km = dd / 1000.0
                n_dfar = int((d_nn_km > args.max_nn_dist_km).sum())
                if n_dfar:
                    print(f"warning: {n_dfar} facets farther than "
                          f"{args.max_nn_dist_km} km from the nearest "
                          "Luttrell sample", file=sys.stderr)
                print(f"Luttrell D: nn dist median {np.median(d_nn_km):.3f} "
                      f"/ max {d_nn_km.max():.3f} km; D median "
                      f"{np.median(D_diff):.2f} MPa")
                d_src = (f"Luttrell csv {args.diff_csv} (col 16), nn median "
                         f"{np.median(d_nn_km):.3f} km")
            else:
                D_diff = np.full(len(tris), args.diff_mpa)
                d_nn_km = None
                d_src = f"scalar --diff-MPa {args.diff_mpa}"
            sig1, sig2, sig3, D_applied, capped = magnitudes_C2(
                Sv_eff, R_facet, D_diff, args.cap_k_max)
            n_capped = int(capped.sum())
            if n_capped:
                print(f"closure diff: capped D on {n_capped}/{len(tris)} "
                      f"facets to the frictional ceiling k<={args.cap_k_max}")
        else:  # ratio (C1)
            sig1, sig2, sig3 = magnitudes_C1(Sv_eff, R_facet, args.k_ratio)
            D_diff = sig1 - sig3
            D_applied = D_diff
            capped = np.zeros(len(tris), dtype=bool)
            n_capped = 0
            d_nn_km = None
            d_src = f"frictional ratio k = {args.k_ratio}"
        # build per-facet effective tensor from the chosen axes
        if args.axes == "csm":
            sigma_facet = build_tensor_from_axes(U_facet, sig1, sig2, sig3)
        else:
            sigma_facet = build_tensor_andersonian(
                shmax_facet, sig1, sig2, sig3)
        if args.taper_zseis_km is not None:
            sigma_facet = apply_deep_taper(
                sigma_facet, depth_facet, args.taper_zseis_km)
        # SHmax validation on the BUILT per-facet tensor
        shmax_recon = shmax_from_sigma(sigma_facet)
        sok = np.isfinite(shmax_facet)
        dsh = angdiff180(shmax_recon[sok], shmax_facet[sok])
        # acceptance checks: shape recovery + ordering/positivity
        R_recovered = shape_recovery_check(sigma_facet)
        ord_ok = bool(np.all(sig1 + 1e-6 >= sig2)
                      and np.all(sig2 + 1e-6 >= sig3)
                      and np.all(sig3 > 0.0))
        depth_cell = {
            "sigma1_eff_MPa_cell": sig1,
            "sigma2_eff_MPa_cell": sig2,
            "sigma3_eff_MPa_cell": sig3,
            "Sv_total_MPa_cell": Sv_total,
            "Sv_eff_MPa_cell": Sv_eff,
            "P_p_MPa_cell": P_p_facet,
            "depth_m_cell": depth_facet,
            "csm_R_cell": R_facet,
            "csm_Aphi_cell": aphi_facet,
            "csm_V2_plunge_deg_cell": v2pl_facet,
            "D_diff_MPa_cell": D_diff,
            "D_applied_MPa_cell": D_applied,
            "capped_cell": capped.astype(np.int32),
        }
        mag_params = {
            "closure": args.closure,
            "differential_source": d_src,
            "k_ratio": args.k_ratio if args.closure == "ratio" else None,
            "pp_model": args.pp_model,
            "density_source": dens_desc,
            "axes": args.axes,
            "orientation_sampling": args.sample,
            "orientation_nn_fallback_facets": n_interp_fb,
            "cap_k_max": args.cap_k_max,
            "cap_method": "frictional ceiling sig1/sig3 <= cap_k_max "
                "(Byerlee/Mohr-Coulomb; D_max = (k-1)*Sv_eff/((1-R)*k+R), "
                "-> 0 at the free surface)",
            "taper_zseis_km": args.taper_zseis_km,
            "n_capped_facets": n_capped,
            "ordering_positivity_ok_all_facets": ord_ok,
            "R_recovery_max_abs_err": float(np.nanmax(
                np.abs(R_recovered - R_facet))),
            "sigma_zz_eff_minus_Sv_eff_MPa": field_stats(
                sigma_facet[:, 2, 2] - Sv_eff),
            "g_accel_m_s2": G_ACCEL,
        }
        print(f"depth-R: closure={args.closure} axes={args.axes} "
              f"capped={n_capped} ordering_ok={ord_ok} "
              f"R_recovery_maxerr={mag_params['R_recovery_max_abs_err']:.2e}")

    val["shmax_recon_vs_model_deg"] = field_stats(dsh)
    print(f"SHmax check: |recon - model| median {np.median(dsh):.3f} deg, "
          f"max {dsh.max():.3f} deg")

    # ---- fault-local overpressure (approach 2): subtract isotropic DeltaPp
    # from the EFFECTIVE tensor at the gate (raises mu_app = tau/sigma_n_eff;
    # tau unchanged because pore pressure is isotropic).  Applied AFTER the
    # depth-R closure diagnostics (SHmax/R recovery are invariant to an
    # isotropic shift) and BEFORE traction resolution. ----
    dpp_op = None
    if op_field is not None:
        dpp_op = delta_pp_mpa(op_field, centroids[:, 0], centroids[:, 1],
                              centroids[:, 2])
        sigma_facet = apply_overpressure_to_tensor(sigma_facet, dpp_op)
        print(f"fault-local overpressure: DeltaPp max {dpp_op.max():.1f} MPa, "
              f"applied on {int((dpp_op > 0).sum())}/{len(tris)} facets "
              f"(field {os.path.basename(args.overpressure_field)})")

    r = resolve_tractions(sigma_facet, strikes, dips, normals, P_p_facet)

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

    sig6 = {  # sampled effective tensor components on the fault (cell)
        "sigma_xx_eff_MPa_cell": sigma_facet[:, 0, 0],
        "sigma_yy_eff_MPa_cell": sigma_facet[:, 1, 1],
        "sigma_zz_eff_MPa_cell": sigma_facet[:, 2, 2],
        "sigma_xy_eff_MPa_cell": sigma_facet[:, 0, 1],
        "sigma_xz_eff_MPa_cell": sigma_facet[:, 0, 2],
        "sigma_yz_eff_MPa_cell": sigma_facet[:, 1, 2],
    }

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
            **sig6,
            "csm_SHmax_deg_cell": shmax_facet,
            "csm_nn_dist_km_cell": nn_dist_km,
            "csm_eig_gap_cell": eiggap_facet,
            **depth_cell,
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

    if args.magnitude_mode == "constant":
        params = {
            "sigma1_eff_MPa": args.s1,
            "sigma2_eff_MPa": args.s2,
            "sigma3_eff_MPa": args.s3,
            "magnitude_provenance": "previous constant-tensor run "
                "(on_fault_stress_projection: eff principals 80/40/35, "
                "P_p 20); CSM magnitudes are orientation-only and NOT "
                "used — adjust via --sigma{1,2,3}-eff-MPa",
            "P_p_MPa": args.P_p,
            "fault_strike_azimuth_hint_deg": args.strike_hint_az,
            "rake_sense": "right-lateral",
            "harmonisation_normal_global": n_global.tolist(),
            "source_convention": "CSM YHSM-2013 csv (Yang & Hauksson), "
                "tension-positive MPa (e,n,u), ORIENTATION ONLY, "
                "depth-invariant; principal axes kept, effective "
                "magnitudes prescribed (compression-positive); "
                "lon/lat -> UTM via EPSG:4326->32611 always_xy; "
                "grid convergence neglected",
        }
    else:
        params = {
            "magnitude_mode": "depth-R",
            "magnitude_provenance": "depth-dependent (plan rev 5): "
                "sig2 = Sv_eff(z) = MUSCAL lithostatic overburden - "
                "hydrostatic P_p; sig1/sig3 closed-form from the closure "
                "and per-point CSM shape ratio R (col 13); built on the "
                "CSM principal axes (compression-positive)",
            "fault_strike_azimuth_hint_deg": args.strike_hint_az,
            "rake_sense": "right-lateral",
            "harmonisation_normal_global": n_global.tolist(),
            "source_convention": "orientations from CSM YHSM-2013 csv "
                "(tension-positive (e,n,u), depth-invariant); magnitudes "
                "depth-dependent (compression-positive); lon/lat -> UTM "
                "via EPSG:4326->32611 always_xy; grid convergence neglected",
            **mag_params,
        }
    summary = {
        # paths recorded as given on the command line (no machine-
        # specific absolute paths in the artifact)
        "input_csm_csv": args.csm_csv,
        "input_mesh": args.mesh,
        "fault_n_cells": int(len(tris)),
        "fault_n_degenerate_basis": n_degen,
        "csm_n_samples": int(len(lon)),
        "params": params,
        "sampling": {
            "method": "nearest neighbour in UTM (x, y); depth ignored "
                      "(model depth-invariant)",
            "nn_dist_km": field_stats(nn_dist_km),
            "n_facets_beyond_max_nn_dist": n_far,
            "csm_eig_gap_sampled": field_stats(eiggap_facet),
        },
        "orientation_validation": val,
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
    # depth-R only: record the differential-stress source path (kept out
    # of the constant-mode summary so that artifact stays byte-identical)
    if args.magnitude_mode == "depth-R" and args.closure == "diff":
        summary["input_diff_csv"] = args.diff_csv
    # overpressure provenance is recorded ONLY when active, so a no-overpressure
    # run's summary.json stays byte-identical to the pre-change tool.
    if op_field is not None:
        summary["fault_local_overpressure"] = {
            "field_npz": args.overpressure_field,
            "DeltaPp_MPa_stats": field_stats(dpp_op),
            "field_meta": op_field["meta"]}

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
            "csm_SHmax_deg": float(shmax_facet[i]),
            "csm_nn_dist_km": float(nn_dist_km[i]),
        }
        if args.magnitude_mode == "depth-R":
            summary["hypocenter_check"].update({
                "depth_m": float(depth_cell["depth_m_cell"][i]),
                "Sv_eff_MPa": float(depth_cell["Sv_eff_MPa_cell"][i]),
                "sigma1_eff_MPa": float(depth_cell["sigma1_eff_MPa_cell"][i]),
                "sigma2_eff_MPa": float(depth_cell["sigma2_eff_MPa_cell"][i]),
                "sigma3_eff_MPa": float(depth_cell["sigma3_eff_MPa_cell"][i]),
                "D_diff_MPa": float(depth_cell["D_diff_MPa_cell"][i]),
                "D_applied_MPa": float(depth_cell["D_applied_MPa_cell"][i]),
                "capped": bool(depth_cell["capped_cell"][i]),
            })
        print("hypocenter facet:", json.dumps(
            summary["hypocenter_check"], indent=2))

    sum_path = os.path.join(args.out_dir, f"{prefix}_summary.json")
    with open(sum_path, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"wrote {sum_path}")
    for k, st in summary["stats"].items():
        if st["min"] is None:
            print(f"  {k:22s} all-NaN")
            continue
        print(f"  {k:22s} min {st['min']:+9.4f}  med {st['median']:+9.4f}  "
              f"max {st['max']:+9.4f}  nan {st['nan_count']}")


if __name__ == "__main__":
    main()
