#!/usr/bin/env python3
"""puml_io.py -- read/write SeisSol PUML (.puml.h5) and decode/encode the
packed boundary word.

PUML layout (SeisSol v1.x, the format every SAFS deck consumes):
    /geometry  (nv,3) float64   vertex coordinates
    /connect   (nt,4) uint64    tet -> vertex
    /boundary  (nt,)  int32     4 x 8-bit face codes, slot s = local face s
    /group     (nt,)  int32     volume attribute (rock = 1)

Local face s of a tet is spanned by LOCAL_FACES[s].  Verified by round-trip on
safalt_0d5Hz_p3.puml.h5: fault faces appear exactly twice, boundary faces once.

BC codes: 0 interior, 1 free surface, 3 dynamic rupture (fault), 5 absorbing.
"""

import h5py
import numpy as np

LOCAL_FACES = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2))

BC_INTERIOR = 0
BC_FREE_SURFACE = 1
BC_DYNAMIC_RUPTURE = 3
BC_ABSORBING = 5


def read_puml(path):
    """Return (geometry, connect, boundary, group)."""
    with h5py.File(path, "r") as f:
        geometry = f["geometry"][:]
        connect = f["connect"][:].astype(np.int64)
        boundary = f["boundary"][:].astype(np.int32)
        group = f["group"][:].astype(np.int32) if "group" in f else np.ones(len(connect), np.int32)
    return geometry, connect, boundary, group


def write_puml(path, geometry, connect, boundary, group):
    """Write a PUML file with the dtypes SeisSol expects."""
    with h5py.File(path, "w") as f:
        f.create_dataset("geometry", data=np.ascontiguousarray(geometry, np.float64))
        f.create_dataset("connect", data=np.ascontiguousarray(connect, np.uint64))
        f.create_dataset("boundary", data=np.ascontiguousarray(boundary, np.int32))
        f.create_dataset("group", data=np.ascontiguousarray(group, np.int32))


def face_code(boundary, slot):
    """8-bit BC code of local face `slot` for every tet.

    Shift through a uint32 view: an int32 word with a code in slot 3 would
    sign-extend under an arithmetic shift.
    """
    return ((np.ascontiguousarray(boundary, np.int32).view(np.uint32) >> np.uint32(8 * slot))
            & np.uint32(0xFF)).astype(np.int32)


def clear_face_code(boundary, slot_mask):
    """Zero the code of face `slot` for the tets flagged in slot_mask[slot]."""
    out = boundary.copy()
    for slot in range(4):
        idx = slot_mask[slot]
        if len(idx) == 0:
            continue
        out[idx] &= ~np.int32(0xFF << (8 * slot))
    return out


def faces_with_code(connect, boundary, code):
    """All faces carrying `code` as (tri_vertices, tet_idx, slot)."""
    tris, tets, slots = [], [], []
    for slot in range(4):
        idx = np.nonzero(face_code(boundary, slot) == code)[0]
        if idx.size:
            tris.append(connect[idx][:, LOCAL_FACES[slot]])
            tets.append(idx)
            slots.append(np.full(idx.size, slot, np.int8))
    if not tris:
        return np.zeros((0, 3), np.int64), np.zeros(0, np.int64), np.zeros(0, np.int8)
    return np.vstack(tris), np.concatenate(tets), np.concatenate(slots)


def tet_signed_volume(points, tets):
    """Signed volume of each tet."""
    p = points[tets]
    return np.einsum(
        "ij,ij->i", p[:, 1] - p[:, 0], np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0])
    ) / 6.0


def tet_edge_lengths(points, tets):
    """(nt,6) edge lengths of every tet."""
    p = points[tets]
    pairs = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    return np.stack([np.linalg.norm(p[:, a] - p[:, b], axis=1) for a, b in pairs], axis=1)


def tet_eta(points, tets):
    """Normalised Joe-Liu quality eta in (0,1]; the project's sliver metric."""
    vol = np.abs(tet_signed_volume(points, tets))
    e2 = (tet_edge_lengths(points, tets) ** 2).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        eta = 12.0 * (3.0 * vol) ** (2.0 / 3.0) / e2
    return np.nan_to_num(eta, nan=0.0, posinf=0.0)
