#!/usr/bin/env python3
"""make_view_xdmf.py -- ParaView-ready XDMF for a PUML mesh.

Writes three views next to each other:

  <tag>_full.xdmf     the VOLUME, referencing results/*.puml.h5 in place (no data
                      copied).  Carries `group` and -- the point of it -- a
                      per-cell `resolved_Hz`, so the resolution gate is visible
                      rather than merely asserted.  A 189 M-cell volume will not
                      open whole on a workstation; use it with a Clip/Slice.
  <tag>_surface.h5    every BOUNDARY face (BC != 0) as triangles, with `bc` and
                      `resolved_Hz` of the cell behind each face.  ~10 M tris.
  <tag>_fault.h5      the BC 3 dynamic-rupture surface only.  ~2.8 M tris.

`resolved_Hz` = (p/4) * Vs / dx, dx = element MAX edge, Vs nearest-grid at the
element barycentre -- the SAME convention gate_census2.py scores, so what you see
is what the gate measured.  Vs comes from native MUSCAL (deck value only where
MUSCAL is NaN), matching the measurement of record.

Usage:
  make_view_xdmf.py --mesh results/x.puml.h5 --tag heavy --order 4 \
      --cvm <deck.nc> --muscal <MUSCAL.nc> [--outdir view] [--no-volume-field]
"""
import argparse, sys
from pathlib import Path
import numpy as np, h5py

sys.path.insert(0, str(Path(__file__).resolve().parent))
from material import Material

PAIRS = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
LOCAL_FACES = ((0,2,1),(0,1,3),(1,2,3),(0,3,2))
CH = 4_000_000
BOX = (132000.0, 724000.0, 3490000.0, 4055000.0)


def face_code(B, slot):
    return ((np.ascontiguousarray(B, np.int32).view(np.uint32) >> np.uint32(8*slot))
            & np.uint32(0xFF)).astype(np.int32)


def grid(name, topo, ntri, nvert, conn_src, geo_src, attrs, dim):
    a = "".join(f"""   <Attribute Name="{n}" Center="Cell">
    <DataItem NumberType="{t}" Precision="{p}" Format="HDF" Dimensions="{ntri}">
     {s}
    </DataItem>
   </Attribute>\n""" for n, t, p, s in attrs)
    return f"""<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
 <Domain>
  <Grid Name="{name}" GridType="Uniform">
   <Topology TopologyType="{topo}" NumberOfElements="{ntri}">
    <DataItem NumberType="Int" Precision="8" Format="HDF" Dimensions="{ntri} {dim}">
     {conn_src}
    </DataItem>
   </Topology>
   <Geometry name="geo" GeometryType="XYZ" NumberOfElements="{nvert}">
    <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="{nvert} 3">
     {geo_src}
    </DataItem>
   </Geometry>
{a}  </Grid>
 </Domain>
</Xdmf>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True); ap.add_argument("--tag", required=True)
    ap.add_argument("--order", type=int, required=True)
    ap.add_argument("--cvm", required=True); ap.add_argument("--muscal")
    ap.add_argument("--outdir", default="view")
    ap.add_argument("--no-volume-field", action="store_true")
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(exist_ok=True)
    mesh = Path(a.mesh)
    rel = Path("..") / mesh.parent.name / mesh.name        # view/ -> results/

    with h5py.File(mesh) as f:
        nt = f["connect"].shape[0]; nv = f["geometry"].shape[0]
    print(f"[mesh] {nt:,} tets  {nv:,} verts   ORDER {a.order} (p{a.order-1}) -> resolved = {(a.order-1)/4:.2f} x Vs/dx")
    # source="muscal" is a REQUEST, not a guarantee: Material returns early when
    # muscal_nc is None, leaving self.M unset, and at() then silently falls back
    # to the DECK cube -- whose 250 m z-binning manufactures shallow failures and
    # is ~6x more demanding than native MUSCAL. Omitting --muscal once here made
    # this mesh look like 942,375 cells below gate when the true MUSCAL count is
    # 35,353. Fail loudly instead of quietly scoring on the wrong cube.
    if not a.muscal:
        raise SystemExit("--muscal is REQUIRED: source='muscal' silently degrades "
                         "to the deck cube without it (measured: 27x more apparent "
                         "gate failures). Pass MUSCAL.nc explicitly.")
    mat = Material(a.cvm, a.muscal, box=BOX, source="muscal")
    # resolved f = (p/4) * Vs/dx with p = ORDER - 1 (p3 -> 0.75, p5 -> 1.25).
    # Using order/4 here would inflate the displayed field by 1.33x.
    scale = (a.order - 1) / 4.0

    # ---- per-cell resolved frequency, and the boundary faces, in one sweep ---
    fh5 = out / f"{a.tag}_cellfreq.h5"
    tri_v, tri_bc, tri_f = [], [], []
    with h5py.File(mesh) as f, h5py.File(fh5, "w") as fo:
        # resolved_Hz alone is not enough to act on: when a cell fails the gate you
        # need to know WHICH term did it -- a slow-material cell that is correctly
        # sized reads the same as an oversized cell in fast rock.  dx and Vs are
        # already computed here, so emit them rather than make them re-derivable.
        d = fo.create_dataset("resolved_Hz", (nt,), np.float32)
        d_dx = fo.create_dataset("edge_max_m", (nt,), np.float32)
        d_vs = fo.create_dataset("vs_ms", (nt,), np.float32)
        for s in range(0, nt, CH):
            T = f["connect"][s:s+CH].astype(np.int64)
            B = f["boundary"][s:s+CH].astype(np.int32)
            V = f["geometry"][:] if s == 0 else V
            p = V[T]
            dx = np.max(np.stack([np.linalg.norm(p[:,j]-p[:,i],axis=1) for i,j in PAIRS],1),1)
            vs = mat.at(p.mean(1))
            d[s:s+len(T)] = (scale*np.where(dx > 0, vs/dx, 0.0)).astype(np.float32)
            d_dx[s:s+len(T)] = dx.astype(np.float32)
            d_vs[s:s+len(T)] = vs.astype(np.float32)
            fr = (scale*np.where(dx > 0, vs/dx, 0.0)).astype(np.float32)
            for slot in range(4):
                m = face_code(B, slot) != 0
                if not m.any(): continue
                tri_v.append(T[m][:, list(LOCAL_FACES[slot])])
                tri_bc.append(face_code(B, slot)[m].astype(np.int32))
                tri_f.append(fr[m])
            del T, B, p, dx, vs, fr
    print(f"[freq] wrote {fh5.name}")

    tri = np.vstack(tri_v); bc = np.concatenate(tri_bc); ff = np.concatenate(tri_f)
    del tri_v, tri_bc, tri_f
    print(f"[surface] {len(tri):,} boundary triangles  (BC "
          + ", ".join(f"{c}:{int((bc==c).sum()):,}" for c in np.unique(bc)) + ")")

    def emit(name, keep):
        t = tri[keep]; b = bc[keep]; fq = ff[keep]
        uv, inv = np.unique(t.ravel(), return_inverse=True)
        h5p = out / f"{a.tag}_{name}.h5"
        with h5py.File(mesh) as f, h5py.File(h5p, "w") as fo:
            fo.create_dataset("geometry", data=f["geometry"][:][uv])
            fo.create_dataset("connect", data=inv.reshape(-1, 3).astype(np.int64))
            fo.create_dataset("bc", data=b)
            fo.create_dataset("resolved_Hz", data=fq)
        (out / f"{a.tag}_{name}.xdmf").write_text(grid(
            name, "Triangle", len(t), len(uv), f"{h5p.name}:/connect",
            f"{h5p.name}:/geometry",
            [("bc", "Int", 4, f"{h5p.name}:/bc"),
             ("resolved_Hz", "Float", 4, f"{h5p.name}:/resolved_Hz")], 3))
        print(f"[{name}] {len(t):,} tris, {len(uv):,} verts -> {h5p.name}")

    emit("surface", np.ones(len(tri), bool))
    emit("fault", bc == 3)

    (out / f"{a.tag}_full.xdmf").write_text(grid(
        "volume", "Tetrahedron", nt, nv, f"{rel.as_posix()}:/connect",
        f"{rel.as_posix()}:/geometry",
        [("group", "Int", 4, f"{rel.as_posix()}:/group"),
         ("resolved_Hz", "Float", 4, f"{fh5.name}:/resolved_Hz"),
         ("edge_max_m", "Float", 4, f"{fh5.name}:/edge_max_m"),
         ("vs_ms", "Float", 4, f"{fh5.name}:/vs_ms")], 4))
    print(f"[full] {out / (a.tag + '_full.xdmf')}")


if __name__ == "__main__":
    main()
