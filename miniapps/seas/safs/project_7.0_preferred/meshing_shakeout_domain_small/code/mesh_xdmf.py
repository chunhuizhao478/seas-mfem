#!/usr/bin/env python3
"""mesh_xdmf.py -- ParaView views of a SeisSol PUML mesh, without duplicating it.

Three modes, cheapest first:

  wrap      a ~2 KB XDMF that POINTS AT the existing .puml.h5 datasets.  Zero
            extra bytes on disk.  ParaView then loads the FULL mesh, so this is
            the right choice up to ~10^7 cells and the wrong one at 10^8.

  fault     just the dynamic-rupture triangles (BC 3), de-duplicated -- each
            fault triangle appears twice in a PUML, once per side.  This is
            what you actually want to look at, and it is ~2 % of the cells.

  surface   fault + free surface + absorbing hull, carrying a `bc` scalar so
            ParaView can colour by boundary type.  Bigger, but it shows the
            domain box and the lid together with the fault.

`fault` and `surface` write a small companion .h5; `wrap` writes nothing but
the XML.  XDMF stores connectivity as Int and points as Float, and the paths
inside <DataItem> are resolved RELATIVE TO THE XDMF FILE, so the .xdmf must sit
next to the .h5 it names.

Usage:
    python mesh_xdmf.py --mesh m.puml.h5 --mode fault [--out view.xdmf]
"""

import argparse
import os

import numpy as np
import h5py

LOCAL_FACES = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2))
CH = 2_000_000
XDMF = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
 <Domain>
  <Grid Name="{name}" GridType="Uniform">
   <Topology TopologyType="{topo}" NumberOfElements="{ne}">
    <DataItem NumberType="Int" Precision="8" Format="HDF" Dimensions="{ne} {nv}">
     {h5}:/{conn}
    </DataItem>
   </Topology>
   <Geometry name="geo" GeometryType="XYZ" NumberOfElements="{np_}">
    <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="{np_} 3">
     {h5}:/{geom}
    </DataItem>
   </Geometry>
{attrs}  </Grid>
 </Domain>
</Xdmf>
"""
ATTR = """   <Attribute Name="{an}" Center="Cell">
    <DataItem NumberType="Int" Precision="4" Format="HDF" Dimensions="{ne}">
     {h5}:/{path}
    </DataItem>
   </Attribute>
"""


def face_code(b, s):
    return ((np.ascontiguousarray(b, np.int32).view(np.uint32) >> np.uint32(8 * s))
            & np.uint32(0xFF)).astype(np.int32)


def collect(path, codes):
    """De-duplicated boundary triangles carrying their BC code."""
    tri, bc = [], []
    with h5py.File(path, "r") as f:
        B = f["boundary"][:].astype(np.int32)
        conn = f["connect"]
        nt = conn.shape[0]
        for s0 in range(0, nt, CH):
            s1 = min(s0 + CH, nt)
            T = conn[s0:s1].astype(np.int64)
            Bc = B[s0:s1]
            for s in range(4):
                fc = face_code(Bc, s)
                for c in codes:
                    m = fc == c
                    if m.any():
                        tri.append(T[m][:, list(LOCAL_FACES[s])])
                        bc.append(np.full(int(m.sum()), c, np.int32))
            del T, Bc
    tri = np.vstack(tri)
    bc = np.concatenate(bc)
    # a fault triangle is stored on BOTH sides; keep one
    k = np.sort(tri, axis=1)
    o = np.lexsort((k[:, 2], k[:, 1], k[:, 0]))
    k, tri, bc = k[o], tri[o], bc[o]
    keep = np.ones(len(k), bool)
    keep[1:] = np.any(k[1:] != k[:-1], axis=1)
    return tri[keep], bc[keep]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--mode", choices=["wrap", "fault", "surface"], default="fault")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    base = a.out or (os.path.splitext(a.mesh)[0].replace(".puml", "")
                     + f"_{a.mode}.xdmf")
    base = base if base.endswith(".xdmf") else base + ".xdmf"
    d = os.path.dirname(os.path.abspath(base))

    if a.mode == "wrap":
        with h5py.File(a.mesh, "r") as f:
            ne, nv = f["connect"].shape
            np_ = f["geometry"].shape[0]
        rel = os.path.relpath(os.path.abspath(a.mesh), d)
        xml = XDMF.format(name="mesh", topo="Tetrahedron", ne=ne, nv=nv, np_=np_,
                          h5=rel, conn="connect", geom="geometry",
                          attrs=ATTR.format(an="group", ne=ne, h5=rel, path="group"))
        open(base, "w").write(xml)
        print(f"[wrap] {base}  ({os.path.getsize(base):,} B, ZERO data copied)")
        print(f"       references {rel}  -> {ne:,} tets, {np_:,} verts")
        print("       ParaView will load the FULL mesh; use --mode fault above ~1e7 cells")
        return

    codes = (3,) if a.mode == "fault" else (1, 3, 5)
    tri, bc = collect(a.mesh, codes)
    with h5py.File(a.mesh, "r") as f:
        P = f["geometry"][:]
    used = np.unique(tri)
    remap = np.full(P.shape[0], -1, np.int64)
    remap[used] = np.arange(len(used))
    tri2 = remap[tri]
    h5out = os.path.splitext(base)[0] + ".h5"
    with h5py.File(h5out, "w") as f:
        f.create_dataset("geometry", data=P[used].astype(np.float64))
        f.create_dataset("connect", data=tri2.astype(np.int64))
        f.create_dataset("bc", data=bc.astype(np.int32))
    rel = os.path.basename(h5out)
    xml = XDMF.format(name=a.mode, topo="Triangle", ne=len(tri2), nv=3,
                      np_=len(used), h5=rel, conn="connect", geom="geometry",
                      attrs=ATTR.format(an="bc", ne=len(tri2), h5=rel, path="bc"))
    open(base, "w").write(xml)
    tot = os.path.getsize(base) + os.path.getsize(h5out)
    src = os.path.getsize(a.mesh)
    print(f"[{a.mode}] {base}")
    print(f"        {len(tri2):,} triangles, {len(used):,} verts")
    print(f"        {tot/1e6:,.1f} MB vs the mesh's {src/1e9:,.2f} GB "
          f"({100*tot/src:.2f} %)")
    u, c = np.unique(bc, return_counts=True)
    print(f"        bc: " + ", ".join(f"{int(k)}={v:,}" for k, v in zip(u, c)))


if __name__ == "__main__":
    main()
