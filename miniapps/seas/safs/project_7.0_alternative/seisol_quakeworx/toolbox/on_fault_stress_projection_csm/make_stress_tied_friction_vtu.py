#!/usr/bin/env python3
"""make_stress_tied_friction_vtu.py — visualize the stress-tied LSW friction
fields (mu_s, mu_d) for safs_seisol_v3, derived from the C1 k=2.39 on-fault
stress projection.

Per fault facet the projection gives the apparent friction mu_app = tau/sigma_n
(the ratio the fault is already loaded to at t=0). The V3 design (PLAN.md
Section 4.2 / 5.2) ties the linear-slip-weakening friction to it:

    mu_s = mu_app + delta_s                       (static; > mu_app -> no t=0 pre-slip)
    mu_d = max(mu_app - delta_d, mu_d_floor)       (dynamic / residual)
    strength excess SE = (mu_s - mu_app)*sigma_n = delta_s*sigma_n
    stress drop     dt = (mu_app - mu_d)*sigma_n
    S-ratio         S  = SE/dt                     (controlled & low -> propagates)

This tool reads the C1 projection fault VTU (mu_apparent_cell,
sigma_n_eff_MPa_cell) + the mesh fault geometry, computes the friction fields
and writes a VTU with both per-cell and (area-averaged) per-vertex fields for
ParaView.  It does NOT apply the deep barrier (mu_s=1e6 for z in
(-20000,-15000) m) — that would blow out the colour scale — but flags those
facets with `deep_barrier_cell` so the locked band is visible.

Usage (conda env 'pythonenv'):
  python3 make_stress_tied_friction_vtu.py \
      --c1-vtu csm_yhsm2013_stress_on_safs_mesh_fault_stress.vtu \
      --mesh ../../safs_seisol_v2_0_0_RSSRW/safs_mesh.puml.h5 \
      --out ../../safs_seisol_v3_0_0_LSW/csm_friction_stress_tied_fault.vtu \
      [--delta-s 0.05] [--delta-d 0.15] [--mu-d-floor 0.0]
      [--barrier-ztop -15000] [--barrier-zbot -20000]
"""
import argparse
import os
import re
import struct
import sys

import numpy as np

# reuse the projection tool's mesh + VTU machinery (guarded main, safe import)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from project_csm_stress_to_vtu import (  # noqa: E402
    load_puml, sanity_check_face_ordering, extract_fault, triangle_geometry,
    cell_to_node_average, write_vtu, field_stats)

_VTK_RD = {"Float64": "<f8", "Float32": "<f4", "Int64": "<i8",
           "Int32": "<i4", "UInt8": "u1"}


def read_vtu_cell_arrays(path, names):
    """Read named appended-binary DataArrays from a VTU written by
    project_csm_stress_to_vtu.write_vtu (little-endian, UInt64 size prefix)."""
    with open(path, "rb") as fh:
        raw = fh.read()
    marker = b'<AppendedData encoding="raw">\n_'
    mk = raw.find(marker)
    if mk < 0:
        sys.exit(f"ERROR: {path} is not the expected appended-raw VTU")
    start = mk + len(marker)
    hdr = raw[:mk].decode("latin1")
    arrs = re.findall(
        r'<DataArray type="([^"]+)" Name="([^"]+)"'
        r'(?: NumberOfComponents="(\d+)")? format="appended" offset="(\d+)"/>',
        hdr)
    out = {}
    for typ, nm, nc, off in arrs:
        if nm not in names:
            continue
        off = int(off)
        nc = int(nc) if nc else 1
        nbytes = struct.unpack_from("<Q", raw, start + off)[0]
        buf = raw[start + off + 8:start + off + 8 + nbytes]
        a = np.frombuffer(buf, dtype=_VTK_RD[typ])
        out[nm] = a.reshape(-1, nc) if nc > 1 else a
    missing = [n for n in names if n not in out]
    if missing:
        sys.exit(f"ERROR: {path} missing cell arrays {missing}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--c1-vtu", required=True,
                    help="C1 k=2.39 projection fault VTU (mu_apparent_cell, "
                    "sigma_n_eff_MPa_cell)")
    ap.add_argument("--mesh", required=True, help="PUML mesh (.puml.h5)")
    ap.add_argument("--out", required=True, help="output friction VTU")
    ap.add_argument("--delta-s", type=float, default=0.05, dest="delta_s",
                    help="static margin: mu_s = mu_app + delta_s (default %(default)s)")
    ap.add_argument("--delta-d", type=float, default=0.15, dest="delta_d",
                    help="dynamic drop: mu_d = max(mu_app-delta_d, floor) (default %(default)s)")
    ap.add_argument("--mu-d-floor", type=float, default=0.0, dest="mu_d_floor",
                    help="lower bound on mu_d (default %(default)s)")
    ap.add_argument("--barrier-ztop", type=float, default=-15000.0,
                    help="deep-barrier band top z (m, default %(default)s)")
    ap.add_argument("--barrier-zbot", type=float, default=-20000.0,
                    help="deep-barrier band bottom z (m, default %(default)s)")
    args = ap.parse_args()

    if args.delta_s <= 0 or args.delta_d <= 0:
        sys.exit("ERROR: --delta-s and --delta-d must be > 0")

    # --- fault geometry (same extract_fault order the C1 VTU was written in) ---
    geom, conn, bc = load_puml(args.mesh)
    sanity_check_face_ordering(geom, conn, bc)
    pts, tris = extract_fault(geom, conn, bc)
    centroids, _normals, areas = triangle_geometry(pts, tris)
    n_facets = len(tris)

    # --- read mu_app, sigma_n_eff from the C1 projection (per-facet, same order) ---
    cd = read_vtu_cell_arrays(
        args.c1_vtu, ["mu_apparent_cell", "sigma_n_eff_MPa_cell"])
    mu_app = cd["mu_apparent_cell"]
    sn = cd["sigma_n_eff_MPa_cell"]            # MPa
    if len(mu_app) != n_facets:
        sys.exit(f"ERROR: C1 VTU has {len(mu_app)} facets, mesh fault has "
                 f"{n_facets} — mesh/VTU mismatch (regenerate the C1 VTU on "
                 "this mesh)")

    # --- stress-tied friction (PLAN.md Section 4.2 / 5.2) ---
    mu_s = mu_app + args.delta_s
    mu_d = np.maximum(mu_app - args.delta_d, args.mu_d_floor)
    se_mpa = (mu_s - mu_app) * sn              # strength excess, MPa
    drop_mpa = (mu_app - mu_d) * sn            # stress drop, MPa
    with np.errstate(divide="ignore", invalid="ignore"):
        s_ratio = np.where(drop_mpa > 0, se_mpa / drop_mpa, np.nan)
    z = centroids[:, 2]
    barrier = ((z > args.barrier_zbot) & (z < args.barrier_ztop)
               ).astype(np.int32)

    n_pre = int((mu_s <= mu_app).sum())        # must be 0 by construction
    print(f"facets: {n_facets}; pre-slip (mu_s<=mu_app): {n_pre} (should be 0)")
    for nm, a in (("mu_app", mu_app), ("mu_s", mu_s), ("mu_d", mu_d),
                  ("strength_excess_MPa", se_mpa), ("stress_drop_MPa", drop_mpa),
                  ("S_ratio", s_ratio)):
        st = field_stats(a)
        print(f"  {nm:20s} min {st['min']:+8.3f}  med {st['median']:+8.3f}  "
              f"max {st['max']:+8.3f}")
    seis = (z < -3000) & (z > -12000)
    print(f"  seismogenic-band (z=-12..-3km) S-ratio: median "
          f"{np.nanmedian(s_ratio[seis]):.3f}  p95 "
          f"{np.nanpercentile(s_ratio[seis], 95):.3f}  "
          f"max {np.nanmax(s_ratio[seis]):.3f}")
    print(f"  deep-barrier facets (z in {args.barrier_zbot:.0f}.."
          f"{args.barrier_ztop:.0f}): {int(barrier.sum())} "
          "(mu_s shown stress-tied here; SeisSol input overrides to 1e6)")

    # --- per-vertex (area-averaged) for smooth ParaView display ---
    n_pts = len(pts)
    mu_s_n = cell_to_node_average(n_pts, tris, mu_s, areas)
    mu_d_n = cell_to_node_average(n_pts, tris, mu_d, areas)
    mu_app_n = cell_to_node_average(n_pts, tris, mu_app, areas)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    write_vtu(
        args.out, pts, tris, 5,
        cell_data={
            "mu_s_cell": mu_s,
            "mu_d_cell": mu_d,
            "mu_app_cell": mu_app,
            "sigma_n_eff_MPa_cell": sn,
            "strength_excess_MPa_cell": se_mpa,
            "stress_drop_MPa_cell": drop_mpa,
            "S_ratio_cell": s_ratio,
            "deep_barrier_cell": barrier,
        },
        point_data={
            "mu_s": mu_s_n,
            "mu_d": mu_d_n,
            "mu_app": mu_app_n,
        })
    print(f"wrote {args.out}")
    print(f"params: delta_s={args.delta_s}, delta_d={args.delta_d}, "
          f"mu_d_floor={args.mu_d_floor} (mu_s=mu_app+delta_s, "
          "mu_d=max(mu_app-delta_d,floor))")


if __name__ == "__main__":
    main()
