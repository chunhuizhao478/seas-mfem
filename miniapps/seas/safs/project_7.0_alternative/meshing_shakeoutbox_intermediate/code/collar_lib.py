#!/usr/bin/env python3
"""collar_lib.py -- shared pieces for the lateral domain extension ("collar").

The parent mesh is FROZEN.  Its vertical absorbing side wall becomes the INNER
boundary of a new collar that fills the region between that wall and an
enlarged outer box.  The collar welds on conformally, so the fault, the free
surface, the gate compliance and the dt floor are carried over bit-for-bit
instead of being re-derived by a global remesh.

This is the lateral analogue of meshing_deep40km/code/build_deep_slab.py, and
deliberately reuses its vocabulary (lid -> wall, slab -> collar).

Geometry facts this relies on (measured, not assumed):
  * the ALT footprint is a rectangle ROTATED 30 deg to the San Andreas strike,
    so the "wall" is four vertical planes, not an axis-aligned box;
  * the ALT lid is EXACTLY flat at z = 0 (0 of 187,161 lid vertex-uses have
    |z| > 1e-6), so the collar lid is flat too and the free-surface definition
    is unchanged;
  * the wall's own vertical grading (fine at the top, coarse at depth) is
    reused as the collar's inner sizing, which is what makes the collar cheap.
"""

import os
import sys
from contextlib import contextmanager

import numpy as np

LOCAL_FACES = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2))
BC_INTERIOR, BC_FREE_SURFACE, BC_DYNAMIC_RUPTURE, BC_ABSORBING = 0, 1, 3, 5
PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))

# --- the target box: union of the ALT footprint and the ShakeOut v1 bbox -----
# ShakeOut v1 (UTM11N, from 04_comparison/HEAVY/shakeout_comparison_HEAVY.json):
#   E 74,850.144 .. 781,468.938 m ; N 3,542,641.466 .. 3,993,303.385 m
# The ALT mesh already reaches N 3,996,866.9, i.e. 3.56 km FURTHER north than
# ShakeOut, so the north edge is KEPT (trimming it would destroy verified mesh).
SHAKEOUT_E = (74850.14406758995, 781468.9381141155)
SHAKEOUT_N = (3542641.465638498, 3993303.3849738695)


@contextmanager
def suppress_stdout(logfile):
    """tetgen/gmsh write pages of statistics to fd 1; park them in a log."""
    sys.stdout.flush()
    saved = os.dup(1)
    with open(logfile, "a") as fh:
        os.dup2(fh.fileno(), 1)
        try:
            yield
        finally:
            sys.stdout.flush()
            os.dup2(saved, 1)
            os.close(saved)


def face_code(boundary, slot):
    """8-bit BC code of local face `slot`; uint32 view so slot 3 cannot sign-extend."""
    return ((np.ascontiguousarray(boundary, np.int32).view(np.uint32) >> np.uint32(8 * slot))
            & np.uint32(0xFF)).astype(np.int32)


def faces_with_code(connect, boundary, code):
    tris, tets, slots = [], [], []
    for slot in range(4):
        idx = np.nonzero(face_code(boundary, slot) == code)[0]
        if idx.size:
            tris.append(connect[idx][:, LOCAL_FACES[slot]])
            tets.append(idx)
            slots.append(np.full(idx.size, slot, np.int8))
    if not tris:
        return (np.zeros((0, 3), np.int64), np.zeros(0, np.int64), np.zeros(0, np.int8))
    return np.vstack(tris), np.concatenate(tets), np.concatenate(slots)


def extract_wall(geometry, connect, boundary, vert_tol=0.5):
    """The VERTICAL absorbing faces = the parent's side wall.

    Returns (P_wall, T_wall_local, wall_global_idx, owner_tet, owner_slot).
    Vertical is decided by the face normal, exactly as the bottom is excluded
    by the same test in the deep-slab builder.
    """
    tris, tets, slots = faces_with_code(connect, boundary, BC_ABSORBING)
    V = geometry[tris]
    n = np.cross(V[:, 1] - V[:, 0], V[:, 2] - V[:, 0])
    ln = np.linalg.norm(n, axis=1, keepdims=True)
    n = n / np.where(ln == 0.0, 1.0, ln)
    vert = np.abs(n[:, 2]) < vert_tol
    g = tris[vert]
    uniq, inv = np.unique(g, return_inverse=True)
    return geometry[uniq], inv.reshape(-1, 3), uniq, tets[vert], slots[vert]


def boundary_loops(tri):
    """All boundary cycles of a triangulation, as ordered vertex index arrays.

    The side wall is an open band, so it has exactly two: the top rim (z = 0)
    and the bottom rim.
    """
    e = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
    uq, cnt = np.unique(e, axis=0, return_counts=True)
    adj = {}
    for a, b in uq[cnt == 1]:
        adj.setdefault(int(a), []).append(int(b))
        adj.setdefault(int(b), []).append(int(a))
    bad = [v for v, nb in adj.items() if len(nb) != 2]
    if bad:
        raise RuntimeError(f"boundary is not a union of simple cycles: "
                           f"{len(bad)} vertices with valence != 2")
    loops, seen = [], set()
    for start in adj:
        if start in seen:
            continue
        loop, prev, cur = [start], None, start
        seen.add(start)
        while True:
            a, b = adj[cur]
            nxt = a if a != prev else b
            if nxt == start:
                break
            loop.append(nxt)
            seen.add(nxt)
            prev, cur = cur, nxt
        loops.append(np.array(loop, np.int64))
    return loops


def signed_area(xy):
    return 0.5 * float(np.sum(xy[:, 0] * np.roll(xy[:, 1], -1)
                              - np.roll(xy[:, 0], -1) * xy[:, 1]))


def as_ccw(loop, xy):
    return (loop, xy) if signed_area(xy) > 0 else (loop[::-1], xy[::-1])


def orient(tris, pts, want_up):
    """Flip triangles so the z-component of the normal has the wanted sign."""
    p = pts[tris]
    nz = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])[:, 2]
    flip = (nz < 0) if want_up else (nz > 0)
    out = tris.copy()
    out[flip] = out[flip][:, [0, 2, 1]]
    return out


def check_closed(F, P=None, piece_of=None):
    """A PLC bounding a solid must use every edge exactly twice.

    On failure, report WHERE the open edges are and which pieces they belong
    to -- a bare count is useless for debugging a 6-piece assembly.
    """
    e = np.sort(np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1)
    uq, cnt = np.unique(e, axis=0, return_counts=True)
    bad = int(np.sum(cnt != 2))
    if not bad:
        return
    open_e = uq[cnt != 2]
    msg = [f"PLC not closed: {bad} edges with != 2 facets "
           f"(counts seen: {sorted(set(cnt.tolist()))})"]
    if P is not None:
        A, B = P[open_e[:, 0]], P[open_e[:, 1]]
        mid = 0.5 * (A + B)
        msg.append(f"  open-edge midpoints: x {mid[:,0].min():,.1f}..{mid[:,0].max():,.1f}  "
                   f"y {mid[:,1].min():,.1f}..{mid[:,1].max():,.1f}  "
                   f"z {mid[:,2].min():,.1f}..{mid[:,2].max():,.1f}")
        for zv in np.unique(np.round(mid[:, 2], 3)):
            n = int((np.round(mid[:, 2], 3) == zv).sum())
            msg.append(f"    z = {zv:>12,.3f} : {n:,} open edges")
        msg.append(f"  edge length min {np.linalg.norm(B-A,axis=1).min():,.1f} "
                   f"max {np.linalg.norm(B-A,axis=1).max():,.1f} m")
    if piece_of is not None:
        for k in (0, 1):
            names, counts = np.unique(piece_of[open_e[:, k]], return_counts=True)
            msg.append(f"  endpoint{k} pieces: "
                       + ", ".join(f"{n}={c:,}" for n, c in zip(names, counts)))
    raise RuntimeError("\n".join(msg))


def tet_signed_volume(points, tets):
    p = points[tets]
    return np.einsum("ij,ij->i", p[:, 1] - p[:, 0],
                     np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0])) / 6.0


def tet_edge_lengths(points, tets):
    p = points[tets]
    return np.stack([np.linalg.norm(p[:, a] - p[:, b], axis=1) for a, b in PAIRS], axis=1)


def tet_eta(points, tets):
    vol = np.abs(tet_signed_volume(points, tets))
    e2 = (tet_edge_lengths(points, tets) ** 2).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        eta = 12.0 * (3.0 * vol) ** (2.0 / 3.0) / e2
    return np.nan_to_num(eta, nan=0.0, posinf=0.0)


class VsGrid:
    """Nearest-grid Vs = sqrt(mu/rho) from the deck's CVM nc.

    THE LOCKED CONVENTION for the gate (check_fsurf_freq.py): dx = element MAX
    edge, Vs sampled at the element BARYCENTER, nearest grid node.  Outside the
    CVM hull the nearest node is used, which is what ASAGI's edge-clamp does at
    runtime -- so the mesh is sized against the material SeisSol will actually
    see there.
    """

    def __init__(self, path):
        import h5py
        with h5py.File(path, "r") as f:
            self.x, self.y, self.z = f["x"][:], f["y"][:], f["z"][:]
            d = f["data"]
            mu = np.asarray(d["mu"]).astype(np.float64)
            rho = np.asarray(d["rho"]).astype(np.float64)
        self.vs = np.sqrt(np.maximum(mu, 0.0) / rho).astype(np.float32)

    @staticmethod
    def _near(axis, q):
        i = np.clip(np.searchsorted(axis, q), 0, len(axis) - 1)
        j = np.clip(i - 1, 0, len(axis) - 1)
        return np.where(np.abs(axis[j] - q) <= np.abs(axis[i] - q), j, i)

    def at(self, P):
        return self.vs[self._near(self.z, P[:, 2]),
                       self._near(self.y, P[:, 1]),
                       self._near(self.x, P[:, 0])].astype(np.float64)

    def column_min(self, xy, z_lo, z_hi):
        """min Vs over [z_lo, z_hi] at each (x,y) -- the refinement-stable target.

        Pooling in z (not 3-D) is the rule proven on the 1 Hz p5 build: the CVM
        is 1500 m laterally but 250 m vertically, so bin changes are essentially
        always vertical, and a full 3-D min over-refines badly.
        """
        kz = (self.z >= z_lo) & (self.z <= z_hi)
        if not kz.any():
            kz = np.zeros(len(self.z), bool)
            kz[int(self._near(self.z, np.array([0.5 * (z_lo + z_hi)]))[0])] = True
        pool = self.vs[kz].min(axis=0)
        return pool[self._near(self.y, xy[:, 1]),
                    self._near(self.x, xy[:, 0])].astype(np.float64)
