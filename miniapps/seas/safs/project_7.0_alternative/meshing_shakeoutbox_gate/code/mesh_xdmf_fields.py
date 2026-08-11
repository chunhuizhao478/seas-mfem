#!/usr/bin/env python3
"""mesh_xdmf_fields.py -- ParaView views carrying the resolution fields per ELEMENT.

THE EQUATION (locked convention for this project)

    f_min = (p/4) * Vs / dx          [Hz]

      dx = element MAXIMUM edge length  [m]  (longest of the 6 tet edges)
      Vs = sqrt(mu/rho), NEAREST-GRID at the element BARYCENTRE  [m/s]
      p  = polynomial order

`p/4` is the sampling correction: a degree-p element carries p sub-intervals
across its span and a sine needs ~4 points per wavelength, so the shortest
resolvable wavelength is lambda_min = (4/p)*dx and f = Vs/lambda_min.
Because dx is the MAX edge, f_min is the element's WORST-direction (minimum)
resolved frequency -- which is what the gate is judged on.

    intermediate, p3 @ 0.5 Hz : f_min = 0.75 * Vs/dx , gate Vs/dx >= 0.6667
    heavy,        p3 @ 0.5 Hz : f_min = 0.75 * Vs/dx , gate Vs/dx >= 0.6667
    heavy,        p5 @ 1.0 Hz : f_min = 1.25 * Vs/dx , gate Vs/dx >= 0.8000

Cell arrays written (all per element):

    f_min_p3   resolved frequency at p3   [Hz]
    f_min_p5   resolved frequency at p5   [Hz]
    Vs         nearest-grid Vs at the barycentre  [m/s]
    edge_min   shortest edge of the element  [m]
    edge_max   longest edge = dx, the one in the formula  [m]
    vs_over_dx the raw gate quantity (compare to 0.6667 / 0.8)
    bc         boundary code, surface/fault modes only (1 free, 3 fault, 5 absorbing)

Modes:
    volume   sidecar .h5 of the per-tet arrays + an XDMF that takes
             geometry/connect straight from the .puml.h5 (not copied) and the
             attributes from the sidecar. XDMF happily mixes files.
    fault    de-duplicated dynamic-rupture triangles, each carrying its OWNING
             tet's fields.
    surface  fault + free surface + absorbing hull, same.

Everything is chunked: at 157.8M tets a single P[T] would be a ~15 GB temporary.
"""
import argparse
import os

import h5py
import numpy as np

LOCAL_FACES = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2))
PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
EI = np.array([p[0] for p in PAIRS])
EJ = np.array([p[1] for p in PAIRS])
CH = 1_000_000

HDR = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
 <Domain>
  <Grid Name="{name}" GridType="Uniform">
   <Topology TopologyType="{topo}" NumberOfElements="{ne}">
    <DataItem NumberType="Int" Precision="8" Format="HDF" Dimensions="{ne} {nv}">
     {conn_h5}:/{conn}
    </DataItem>
   </Topology>
   <Geometry name="geo" GeometryType="XYZ" NumberOfElements="{np_}">
    <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="{np_} 3">
     {geo_h5}:/{geom}
    </DataItem>
   </Geometry>
{attrs}  </Grid>
 </Domain>
</Xdmf>
"""
ATTR = """   <Attribute Name="{an}" AttributeType="Scalar" Center="Cell">
    <DataItem NumberType="{nt}" Precision="4" Format="HDF" Dimensions="{ne}">
     {h5}:/{path}
    </DataItem>
   </Attribute>
"""
FIELDS = ("f_min_p3", "f_min_p5", "Vs", "edge_min", "edge_max", "vs_over_dx")


def face_code(b, s):
    return ((np.ascontiguousarray(b, np.int32).view(np.uint32) >> np.uint32(8 * s))
            & np.uint32(0xFF)).astype(np.int32)


def tet_fields(G, T, vs):
    """Per-element (edge_min, edge_max, Vs, f_p3, f_p5) for one chunk of tets."""
    P = G[T]
    E = np.linalg.norm(P[:, EI] - P[:, EJ], axis=2)
    emin = E.min(1)
    emax = E.max(1)
    v = vs.at(P.mean(1))
    r = np.where(emax > 0, v / emax, 0.0)
    return (emin.astype(np.float32), emax.astype(np.float32), v.astype(np.float32),
            (0.75 * r).astype(np.float32), (1.25 * r).astype(np.float32),
            r.astype(np.float32))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--mode", choices=["volume", "fault", "surface"], default="fault")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cvm", default=None)
    a = ap.parse_args()

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]
                           / "meshing_shakeoutbox_intermediate" / "code"))
    from collar_lib import VsGrid
    from muscal_vs import MuscalVs
    CVM = a.cvm or ("/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_"
                    "RSSRW_ALT_THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_"
                    "k1p40_attenuation_deep40km/safs_material_cvm.nc")
    vs = MuscalVs(backup=VsGrid(CVM))

    base = a.out if a.out.endswith(".xdmf") else a.out + ".xdmf"
    d = os.path.dirname(os.path.abspath(base))
    h5out = os.path.splitext(base)[0] + ".h5"

    with h5py.File(a.mesh, "r") as f:
        nt = f["connect"].shape[0]
        nvert = f["geometry"].shape[0]

    if a.mode == "volume":
        # geometry + connect stay in the .puml.h5; only the fields are written.
        with h5py.File(a.mesh, "r") as f, h5py.File(h5out, "w") as o:
            conn = f["connect"]
            G = f["geometry"][:]
            ds = {k: o.create_dataset(k, (nt,), np.float32, chunks=(min(nt, 1 << 20),),
                                      compression="gzip", compression_opts=1)
                  for k in FIELDS}
            for s0 in range(0, nt, CH):
                s1 = min(s0 + CH, nt)
                T = conn[s0:s1].astype(np.int64)
                vals = tet_fields(G, T, vs)
                for k, v in zip(("edge_min", "edge_max", "Vs", "f_min_p3",
                                 "f_min_p5", "vs_over_dx"), vals):
                    ds[k][s0:s1] = v
                del T
                if (s0 // CH) % 25 == 0:
                    print(f"    ... {s1:,}/{nt:,}", flush=True)
        rel_mesh = os.path.relpath(os.path.abspath(a.mesh), d)
        rel_side = os.path.basename(h5out)
        attrs = "".join(ATTR.format(an=k, nt="Float", ne=nt, h5=rel_side, path=k)
                        for k in FIELDS)
        attrs += ATTR.format(an="group", nt="Int", ne=nt, h5=rel_mesh, path="group")
        xml = HDR.format(name="volume", topo="Tetrahedron", ne=nt, nv=4, np_=nvert,
                         conn_h5=rel_mesh, conn="connect",
                         geo_h5=rel_mesh, geom="geometry", attrs=attrs)
        open(base, "w").write(xml)
        print(f"[volume] {base}")
        print(f"         {nt:,} tets; geometry/connect referenced in place, "
              f"fields sidecar {os.path.getsize(h5out)/1e6:,.0f} MB")
        vs.report()
        return

    # ---- boundary modes: keep the OWNING tet index so fields can be attached --
    codes = (3,) if a.mode == "fault" else (1, 3, 5)
    tris, bcs, fvals = [], [], []
    with h5py.File(a.mesh, "r") as f:
        B = f["boundary"][:].astype(np.int32)
        conn = f["connect"]
        G = f["geometry"][:]
        for s0 in range(0, nt, CH):
            s1 = min(s0 + CH, nt)
            Bc = B[s0:s1]
            hit = np.zeros(s1 - s0, bool)
            for s in range(4):
                fc = face_code(Bc, s)
                for c in codes:
                    hit |= fc == c
            if not hit.any():
                continue
            idx = np.nonzero(hit)[0]
            Tc = conn[s0:s1][idx].astype(np.int64)
            Bs = Bc[idx]
            # Fields are computed HERE, on the chunk's own connectivity, and
            # carried along.  Collecting owner indices and re-reading them
            # afterwards looks tidier but is a trap: h5py fancy-indexing a few
            # million SCATTERED rows out of a 10^7-10^8-row dataset is orders of
            # magnitude slower than the streaming pass that produced them, and
            # it stalled this tool for >20 min on the intermediate surface.
            got = tet_fields(G, Tc, vs)          # (emin, emax, Vs, f3, f5, ratio)
            for s in range(4):
                fc = face_code(Bs, s)
                m = np.zeros(len(idx), bool)
                for c in codes:
                    m |= fc == c
                if not m.any():
                    continue
                tris.append(Tc[m][:, list(LOCAL_FACES[s])])
                bcs.append(fc[m])
                fvals.append(np.stack([g[m] for g in got], axis=1))
            del Tc, got
    tri = np.vstack(tris)
    bc = np.concatenate(bcs)
    FV = np.vstack(fvals)                        # (ntri, 6) float32
    del tris, bcs, fvals

    # de-duplicate the fault: a PUML stores each fault triangle twice, once per
    # side. Keep the first occurrence, so each surviving triangle carries the
    # fields of ONE of its two neighbouring tets.
    key = np.sort(tri, axis=1)
    order = np.lexsort((key[:, 2], key[:, 1], key[:, 0]))
    key = key[order]
    keep_sorted = np.ones(len(key), bool)
    keep_sorted[1:] = (key[1:] != key[:-1]).any(1)
    keep = order[keep_sorted]
    keep.sort()
    tri, bc, FV = tri[keep], bc[keep], FV[keep]

    order_map = {k: FV[:, i] for i, k in enumerate(
        ("edge_min", "edge_max", "Vs", "f_min_p3", "f_min_p5", "vs_over_dx"))}

    used = np.unique(tri)
    remap = np.full(nvert, -1, np.int64)
    remap[used] = np.arange(len(used))
    with h5py.File(h5out, "w") as o:
        o.create_dataset("geometry", data=G[used].astype(np.float64))
        o.create_dataset("connect", data=remap[tri].astype(np.int64))
        o.create_dataset("bc", data=bc.astype(np.int32))
        for k in FIELDS:
            o.create_dataset(k, data=order_map[k])
    rel = os.path.basename(h5out)
    attrs = "".join(ATTR.format(an=k, nt="Float", ne=len(tri), h5=rel, path=k)
                    for k in FIELDS)
    attrs += ATTR.format(an="bc", nt="Int", ne=len(tri), h5=rel, path="bc")
    xml = HDR.format(name=a.mode, topo="Triangle", ne=len(tri), nv=3, np_=len(used),
                     conn_h5=rel, conn="connect", geo_h5=rel, geom="geometry",
                     attrs=attrs)
    open(base, "w").write(xml)
    print(f"[{a.mode}] {base}")
    print(f"        {len(tri):,} triangles, {len(used):,} verts, "
          f"{os.path.getsize(h5out)/1e6:,.1f} MB")
    u, c = np.unique(bc, return_counts=True)
    print("        bc: " + ", ".join(f"{int(k)}={v:,}" for k, v in zip(u, c)))
    f3 = order_map["f_min_p3"]
    print(f"        f_min_p3 [Hz]: min {f3.min():.3f}  p1 {np.percentile(f3,1):.3f}  "
          f"med {np.median(f3):.3f}")
    vs.report()


if __name__ == "__main__":
    main()
