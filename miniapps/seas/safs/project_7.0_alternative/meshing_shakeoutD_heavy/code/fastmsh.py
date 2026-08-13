#!/usr/bin/env python3
"""fastmsh.py -- memory-lean binary gmsh 2.2 reader/writer, drop-in for the
subset of the meshio API this toolchain uses.

Why this exists: `meshio.read` on the 1.67 GB binary .msh of the 48.9 M-tet
deep40km mesh exhausts 36 GB of RAM and kills the process (it materialises
per-element Python structures).  The same call sits in red_refine_fault.py,
mmg_refine_sizemap.py, reorient_negative_tets.py and msh_to_puml.py, so one
lean implementation unblocks the whole ladder.

Provides `read`, `write`, `Mesh`, `CellBlock` with exactly the attributes the
tools touch: `.points`, `.cells` (list of blocks with `.type` / `.data`) and
`.cell_data["gmsh:physical"]` (list parallel to `.cells`).

gmsh 2.2 binary layout:
    $MeshFormat / "2.2 1 8" / int32(1) / $EndMeshFormat
    $Nodes    / N / N x (int32 tag, 3 x float64)          / $EndNodes
    $Elements / N / blocks of [int32 type, int32 count, int32 ntags]
                    then count x int32(id, tags..., nodes...) / $EndElements
"""

import numpy as np

GMSH_TYPE = {2: ("triangle", 3), 4: ("tetra", 4)}
TYPE_GMSH = {"triangle": (2, 3), "tetra": (4, 4)}


class CellBlock:
    def __init__(self, type, data):
        self.type = type
        self.data = data

    def __iter__(self):                     # some call sites unpack (type, data)
        return iter((self.type, self.data))

    def __len__(self):
        return len(self.data)


class Mesh:
    def __init__(self, points, cells, cell_data=None, **kw):
        self.points = points
        self.cells = [c if isinstance(c, CellBlock) else CellBlock(c[0], c[1])
                      for c in cells]
        self.cell_data = cell_data or {}


def _line(f):
    return f.readline().decode("ascii", "replace").strip()


def read(path):
    with open(path, "rb") as f:
        if _line(f) != "$MeshFormat":
            raise ValueError("not a gmsh file")
        ver, ftype, dsize = _line(f).split()
        if ftype != "1":
            raise ValueError("fastmsh reads BINARY gmsh22 only")
        if np.frombuffer(f.read(4), "<i4")[0] != 1:
            raise ValueError("unexpected endianness marker")
        pos = f.tell()
        if f.read(1) != b"\n":
            f.seek(pos)
        if _line(f) != "$EndMeshFormat":
            raise ValueError("bad $EndMeshFormat")

        if _line(f) != "$Nodes":
            raise ValueError("expected $Nodes")
        nn = int(_line(f))
        nd = np.dtype([("tag", "<i4"), ("xyz", "<f8", (3,))])
        arr = np.fromfile(f, dtype=nd, count=nn)
        if len(arr) != nn:
            raise ValueError("truncated $Nodes")
        points = np.ascontiguousarray(arr["xyz"])
        tags = arr["tag"]
        contiguous = tags[0] == 1 and tags[-1] == nn
        lookup = None
        if not contiguous:                  # rare; build an explicit map
            lookup = np.zeros(int(tags.max()) + 1, np.int64)
            lookup[tags] = np.arange(nn)
        del arr, tags
        pos = f.tell()
        if f.read(1) != b"\n":
            f.seek(pos)
        if _line(f) != "$EndNodes":
            raise ValueError("bad $EndNodes")

        if _line(f) != "$Elements":
            raise ValueError("expected $Elements")
        ne = int(_line(f))
        blocks, physical, got = [], [], 0
        while got < ne:
            hdr = np.fromfile(f, dtype="<i4", count=3)
            etype, count, ntags = int(hdr[0]), int(hdr[1]), int(hdr[2])
            if etype not in GMSH_TYPE:
                raise ValueError(f"unsupported gmsh element type {etype}")
            name, npts = GMSH_TYPE[etype]
            stride = 1 + ntags + npts
            raw = np.fromfile(f, dtype="<i4", count=count * stride)
            if len(raw) != count * stride:
                raise ValueError("truncated $Elements")
            raw = raw.reshape(count, stride)
            conn = raw[:, 1 + ntags:].astype(np.int64)
            conn = (lookup[conn] if lookup is not None else conn - 1)
            blocks.append(CellBlock(name, conn))
            physical.append(raw[:, 1].astype(np.int64) if ntags >= 1
                            else np.zeros(count, np.int64))
            del raw
            got += count
    return Mesh(points, blocks, {"gmsh:physical": physical})


def write(path, mesh, file_format="gmsh22", binary=True, **kw):
    if not binary:
        raise ValueError("fastmsh writes BINARY gmsh22 only")
    phys = mesh.cell_data.get("gmsh:physical")
    pts = np.ascontiguousarray(mesh.points, np.float64)
    with open(path, "wb") as f:
        f.write(b"$MeshFormat\n2.2 1 8\n")
        f.write(np.array([1], "<i4").tobytes())
        f.write(b"\n$EndMeshFormat\n")

        f.write(b"$Nodes\n%d\n" % len(pts))
        nd = np.dtype([("tag", "<i4"), ("xyz", "<f8", (3,))])
        CH = 4_000_000
        for s in range(0, len(pts), CH):
            e = min(s + CH, len(pts))
            buf = np.empty(e - s, dtype=nd)
            buf["tag"] = np.arange(s + 1, e + 1, dtype="<i4")
            buf["xyz"] = pts[s:e]
            f.write(buf.tobytes())
        f.write(b"\n$EndNodes\n")

        total = sum(len(cb.data) for cb in mesh.cells)
        f.write(b"$Elements\n%d\n" % total)
        eid = 1
        for bi, cb in enumerate(mesh.cells):
            etype, npts = TYPE_GMSH[cb.type]
            data = np.ascontiguousarray(cb.data)
            n = len(data)
            pv = (np.asarray(phys[bi], np.int64) if phys is not None
                  else np.zeros(n, np.int64))
            f.write(np.array([etype, n, 2], "<i4").tobytes())
            for s in range(0, n, CH):
                e = min(s + CH, n)
                blk = np.empty((e - s, 3 + npts), "<i4")
                blk[:, 0] = np.arange(eid + s, eid + e)
                blk[:, 1] = pv[s:e]
                blk[:, 2] = pv[s:e]
                blk[:, 3:] = data[s:e] + 1
                f.write(blk.tobytes())
            eid += n
        f.write(b"\n$EndElements\n")
