#!/usr/bin/env python3
"""smooth_fault_volume_ratio.py — reduce the Eq.(18) fault volume ratio r_v of a
non-planar embedded-fault mesh by near-fault volume-balancing smoothing.

PLAN: experimental_mesh_refinement/PLAN_rv_volume_smoothing.md

Eq.(18) (Zhang et al. 2023, JGR Solid Earth 10.1029/2022JB025817):
    r_v = max{V_A/V_B, V_B/V_A}
per fault triangle, A/B the two tets sharing it.  The two tets share the SAME
triangle, so V = (1/3)*area*h and  r_v = max(h_A/h_B, h_B/h_A)  where h is the
apex node's perpendicular distance to the (fixed) fault triangle plane.  We
reduce r_v by moving ONLY the near-fault interior apex nodes to equalize the two
apex heights, keeping every fault node and every domain-boundary node fixed
(so fault geometry, the z=0 trace, the domain box, and the mesh topology are all
unchanged — a fully unstructured Delaunay mesh, just nudged near the fault).

Usage:
    conda activate pythonenv
    python smooth_fault_volume_ratio.py --in IN.msh --out OUT.msh \
        [--omega 0.5] [--floor-frac 0.2] [--max-sweeps 60] [--tol 1e-3]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# Gmsh element type codes
GMSH_LINE = 1
GMSH_TRI = 2
GMSH_TET = 4

ROCK_PHYS = 1
FAULT_PHYS = 101
BOUNDARY_PHYS = (102, 103, 104)   # top, bottom, sides


# --------------------------------------------------------------------------- #
#  Gmsh v2.2 ASCII reader / writer (preserves node ids, element block verbatim)
# --------------------------------------------------------------------------- #
class GmshMesh:
    """Minimal Gmsh v2.2 ASCII container that preserves everything needed to
    re-emit the mesh byte-for-byte except for moved node coordinates."""

    def __init__(self):
        self.phys_lines: list[str] = []     # raw $PhysicalNames body (or [])
        self.node_ids: np.ndarray = None    # (N,) original gmsh node ids
        self.node_lines: list[str] = []     # raw node lines (file order)
        self.coords: np.ndarray = None      # (N,3) float
        self.elem_lines: list[str] = []     # raw element lines (verbatim)
        self.id2idx: dict[int, int] = {}
        # derived (0-based node-index) connectivity
        self.tets: np.ndarray = None        # (M,4) phys==1 tetra
        self.fault_tris: np.ndarray = None  # (F,3) phys==101 triangles
        self.fault_node_idx: np.ndarray = None
        self.boundary_node_idx: np.ndarray = None


def read_gmsh22(path: Path, fault_phys: int = FAULT_PHYS,
                rock_phys: int = ROCK_PHYS,
                boundary_phys=BOUNDARY_PHYS) -> GmshMesh:
    """Read a Gmsh v2.2 ASCII mesh into a GmshMesh, parsing rock tetra
    (rock_phys), fault triangles (fault_phys) and boundary triangles
    (boundary_phys) into 0-based index arrays, while keeping the raw
    node/element lines.

    Defaults are the SAFS tags (rock=1, fault=101, boundary=102/103/104).  For
    TPV102 pass fault_phys=3, rock_phys=1, boundary_phys=(1, 5) — note a
    triangle is classified by (etype, phys), so tag 1 used by BOTH the rock
    tetra and the top-surface triangles is handled correctly."""
    boundary_phys = tuple(boundary_phys)
    g = GmshMesh()
    with open(path) as f:
        lines = f.read().splitlines()
    i, n = 0, len(lines)
    fault_ids: set[int] = set()
    bnd_ids: set[int] = set()
    tets: list[list[int]] = []
    fault_tris: list[list[int]] = []
    while i < n:
        tok = lines[i].strip()
        if tok == "$PhysicalNames":
            cnt = int(lines[i + 1])
            g.phys_lines = lines[i + 2: i + 2 + cnt]
            i += 2 + cnt + 1            # skip $EndPhysicalNames
            continue
        if tok == "$Nodes":
            cnt = int(lines[i + 1])
            ids = np.empty(cnt, dtype=np.int64)
            coords = np.empty((cnt, 3), dtype=np.float64)
            raw = lines[i + 2: i + 2 + cnt]
            for k, ln in enumerate(raw):
                p = ln.split()
                ids[k] = int(p[0])
                coords[k, 0] = float(p[1])
                coords[k, 1] = float(p[2])
                coords[k, 2] = float(p[3])
            g.node_ids = ids
            g.node_lines = list(raw)
            g.coords = coords
            g.id2idx = {int(t): k for k, t in enumerate(ids)}
            i += 2 + cnt + 1            # skip $EndNodes
            continue
        if tok == "$Elements":
            cnt = int(lines[i + 1])
            raw = lines[i + 2: i + 2 + cnt]
            g.elem_lines = list(raw)
            id2idx = g.id2idx
            for ln in raw:
                p = ln.split()
                etype = int(p[1])
                ntags = int(p[2])
                phys = int(p[3])
                conn_ids = p[3 + ntags:]
                if etype == GMSH_TET and phys == rock_phys:
                    tets.append([id2idx[int(c)] for c in conn_ids])
                elif etype == GMSH_TRI and phys == fault_phys:
                    idx = [id2idx[int(c)] for c in conn_ids]
                    fault_tris.append(idx)
                    fault_ids.update(idx)
                elif etype == GMSH_TRI and phys in boundary_phys:
                    bnd_ids.update(id2idx[int(c)] for c in conn_ids)
            i += 2 + cnt + 1            # skip $EndElements
            continue
        i += 1

    if g.coords is None:
        raise ValueError(f"{path}: no $Nodes section")
    if not tets:
        raise ValueError(f"{path}: no rock tetrahedra (physical {rock_phys})")
    if not fault_tris:
        raise ValueError(f"{path}: no fault triangles (physical {fault_phys})")
    g.tets = np.asarray(tets, dtype=np.int64)
    g.fault_tris = np.asarray(fault_tris, dtype=np.int64)
    g.fault_node_idx = np.fromiter(sorted(fault_ids), dtype=np.int64)
    g.boundary_node_idx = np.fromiter(sorted(bnd_ids), dtype=np.int64)
    return g


def write_gmsh22(path: Path, g: GmshMesh, coords_new: np.ndarray,
                 moved_mask: np.ndarray) -> None:
    """Re-emit the mesh with moved node coordinates updated.  Unmoved node lines
    and the entire element block are written verbatim, so fault/boundary nodes
    and all connectivity/tags are byte-identical to the input."""
    out = []
    out.append("$MeshFormat")
    out.append("2.2 0 8")
    out.append("$EndMeshFormat")
    if g.phys_lines:
        out.append("$PhysicalNames")
        out.append(str(len(g.phys_lines)))
        out.extend(g.phys_lines)
        out.append("$EndPhysicalNames")
    out.append("$Nodes")
    out.append(str(len(g.node_ids)))
    for k in range(len(g.node_ids)):
        if moved_mask[k]:
            nid = int(g.node_ids[k])
            x, y, z = coords_new[k]
            out.append(f"{nid} {x:.12g} {y:.12g} {z:.12g}")
        else:
            out.append(g.node_lines[k])          # verbatim (byte-identical)
    out.append("$EndNodes")
    out.append("$Elements")
    out.append(str(len(g.elem_lines)))
    out.extend(g.elem_lines)                     # verbatim
    out.append("$EndElements")
    Path(path).write_text("\n".join(out) + "\n")


# --------------------------------------------------------------------------- #
#  Geometry helpers
# --------------------------------------------------------------------------- #
def signed_volumes(coords: np.ndarray, tets: np.ndarray) -> np.ndarray:
    """Signed volume of every tet (sign follows node ordering)."""
    p0 = coords[tets[:, 0]]
    d1 = coords[tets[:, 1]] - p0
    d2 = coords[tets[:, 2]] - p0
    d3 = coords[tets[:, 3]] - p0
    return np.einsum("ij,ij->i", d1, np.cross(d2, d3)) / 6.0


_EDGE_PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))


def _tet_edge_lengths(coords: np.ndarray, tets: np.ndarray) -> np.ndarray:
    """(M,6) edge lengths for every tet (same edge ordering as check_mesh)."""
    a = coords[tets[:, [p[0] for p in _EDGE_PAIRS]]]
    b = coords[tets[:, [p[1] for p in _EDGE_PAIRS]]]
    return np.linalg.norm(a - b, axis=2)


def tet_etas(coords: np.ndarray, tets: np.ndarray) -> np.ndarray:
    """Joe-Liu shape quality eta = 12*(3V)^(2/3)/sum(edge^2) in [0,1] per tet.

    Identical definition to meshing/code/check_mesh_quality._tet_metrics, so the
    smoothing guard speaks the same quality language as the project gate."""
    e = _tet_edge_lengths(coords, tets)
    sumL2 = (e ** 2).sum(axis=1)
    V = np.abs(signed_volumes(coords, tets))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(sumL2 > 0.0,
                        12.0 * np.cbrt((3.0 * V) ** 2) / sumL2, 0.0)


def tet_min_edges(coords: np.ndarray, tets: np.ndarray) -> np.ndarray:
    """Shortest edge of every tet (for the Q1 min-edge guard)."""
    return _tet_edge_lengths(coords, tets).min(axis=1)


def build_node_tets_csr(tets: np.ndarray, n_nodes: int):
    """CSR node->incident-tet adjacency: (offsets[N+1], csr_tets).

    node n's incident tet indices are csr_tets[offsets[n]:offsets[n+1]]."""
    M = tets.shape[0]
    flat = tets.reshape(-1)
    tet_of = np.repeat(np.arange(M, dtype=np.int64), 4)
    order = np.argsort(flat, kind="stable")
    csr_tets = tet_of[order]
    counts = np.bincount(flat, minlength=n_nodes)
    offsets = np.zeros(n_nodes + 1, dtype=np.int64)
    offsets[1:] = np.cumsum(counts)
    return offsets, csr_tets


def _incident_ok(coords, inc_tets, sgn0_inc, eta_ok_inc, edge_ok_inc) -> bool:
    """True iff every tet in inc_tets keeps its sign and clears the eta/edge
    guards at the current coords (cheap sign check first)."""
    V = signed_volumes(coords, inc_tets)
    if np.any(np.sign(V) != sgn0_inc):
        return False
    if np.any(tet_min_edges(coords, inc_tets) < edge_ok_inc):
        return False
    if np.any(tet_etas(coords, inc_tets) < eta_ok_inc):
        return False
    return True


class FaultTopo:
    """Per-fault-face apex pairing + fixed face planes (constant geometry)."""

    def __init__(self, apex_a, apex_b, n_f, c_f):
        self.apex_a = apex_a    # (F,) node index of +side apex
        self.apex_b = apex_b    # (F,) node index of -side apex
        self.n_f = n_f          # (F,3) unit face normal (constant)
        self.c_f = c_f          # (F,3) face centroid (constant)


def build_fault_apex_map(coords: np.ndarray, tets: np.ndarray,
                         fault_tris: np.ndarray) -> FaultTopo:
    """For each fault triangle, find the two owning tets and their apex nodes.

    Raises if any fault triangle is not shared by exactly two tets (the embedding
    invariant)."""
    fault_keys = {frozenset(int(v) for v in tri) for tri in fault_tris}
    face_to_apex: dict[frozenset, list[int]] = {}
    faces = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
    for t in tets:
        t0, t1, t2, t3 = int(t[0]), int(t[1]), int(t[2]), int(t[3])
        nodes = (t0, t1, t2, t3)
        for fa in faces:
            key = frozenset((nodes[fa[0]], nodes[fa[1]], nodes[fa[2]]))
            if key in fault_keys:
                apex = (set(nodes) - key).pop()
                face_to_apex.setdefault(key, []).append(apex)

    F = fault_tris.shape[0]
    apex_a = np.empty(F, dtype=np.int64)
    apex_b = np.empty(F, dtype=np.int64)
    n_f = np.empty((F, 3), dtype=np.float64)
    c_f = np.empty((F, 3), dtype=np.float64)
    bad = 0
    for j, tri in enumerate(fault_tris):
        key = frozenset(int(v) for v in tri)
        apexes = face_to_apex.get(key, [])
        if len(apexes) != 2:
            bad += 1
            # leave a degenerate self-pair; flagged below
            a = apexes[0] if apexes else int(tri[0])
            b = apexes[1] if len(apexes) > 1 else a
        else:
            a, b = apexes
        v0, v1, v2 = coords[tri[0]], coords[tri[1]], coords[tri[2]]
        nrm = np.cross(v1 - v0, v2 - v0)
        ln = np.linalg.norm(nrm)
        n_f[j] = nrm / ln if ln > 0 else np.array([0.0, 0.0, 1.0])
        c_f[j] = (v0 + v1 + v2) / 3.0
        apex_a[j], apex_b[j] = a, b
    if bad:
        raise ValueError(
            f"{bad} fault triangle(s) are not shared by exactly 2 tets; the "
            f"fault is not cleanly embedded — refusing to smooth.")
    # Orient so apex_a is on the +n_f side, apex_b on the -n_f side (cosmetic;
    # the algorithm is sign-agnostic but this keeps reporting consistent).
    ha = np.einsum("ij,ij->i", coords[apex_a] - c_f, n_f)
    flip = ha < 0
    apex_a[flip], apex_b[flip] = apex_b[flip], apex_a[flip]
    return FaultTopo(apex_a, apex_b, n_f, c_f)


def rv_values(coords: np.ndarray, topo: FaultTopo) -> np.ndarray:
    """r_v per fault face = max(|h_a|/|h_b|, |h_b|/|h_a|)."""
    ha = np.abs(np.einsum("ij,ij->i", coords[topo.apex_a] - topo.c_f, topo.n_f))
    hb = np.abs(np.einsum("ij,ij->i", coords[topo.apex_b] - topo.c_f, topo.n_f))
    eps = np.finfo(float).tiny
    ha = np.maximum(ha, eps)
    hb = np.maximum(hb, eps)
    return np.maximum(ha / hb, hb / ha)


def rv_stats(coords: np.ndarray, topo: FaultTopo) -> dict:
    r = np.sort(rv_values(coords, topo))
    n = r.size
    return {
        "n": int(n), "min": float(r[0]), "mean": float(r.mean()),
        "median": float(r[n // 2]), "p95": float(r[int(0.95 * (n - 1))]),
        "p99": float(r[int(0.99 * (n - 1))]), "max": float(r[-1]),
        "n_gt15": int((r > 1.5).sum()), "n_gt20": int((r > 2.0).sum()),
    }


def build_free_axes(g: "GmshMesh", slide_boundary: bool = False) -> np.ndarray:
    """Return a (N,3) bool array: which coordinate axes each node may move along.

    - Fault nodes: no axis free (fault geometry frozen).
    - Boundary nodes: with slide_boundary=False, no axis free (domain frozen,
      the plan default). With slide_boundary=True, free to slide *within* their
      axis-aligned box face(s): the axis normal to each face the node lies on is
      locked, so the node stays exactly on that planar face (box edges keep one
      free axis, corners none) — the domain box is preserved exactly while the
      in-plane position relaxes.
    - All other (interior) nodes: all three axes free.
    """
    N = len(g.node_ids)
    free = np.ones((N, 3), dtype=bool)
    free[g.fault_node_idx] = False
    if not slide_boundary:
        free[g.boundary_node_idx] = False
        return free
    # slide_boundary: lock only the axis normal to each box face the node is on
    c = g.coords
    lo = c.min(axis=0)
    hi = c.max(axis=0)
    tol = 1e-6 * np.maximum(1.0, hi - lo)
    bnd = g.boundary_node_idx
    for ax in range(3):
        on_lo = np.abs(c[bnd, ax] - lo[ax]) <= tol[ax]
        on_hi = np.abs(c[bnd, ax] - hi[ax]) <= tol[ax]
        locked = on_lo | on_hi
        free[bnd[locked], ax] = False
    # fault nodes (incl. the z=0 trace, which are also boundary) stay fully fixed
    free[g.fault_node_idx] = False
    return free


# --------------------------------------------------------------------------- #
#  Targeted Gauss-Seidel ceiling enforcement
# --------------------------------------------------------------------------- #
def targeted_ceiling_pass(X, tets, topo, movable, free_axes, sgn0, eta_ok,
                          edge_ok, offsets, csr_tets, rv_ceiling, *,
                          omega_gs=0.5, margin=0.98, passes=30, n_bisect=18,
                          verbose=True):
    """Sequentially (Gauss-Seidel) drive every over-ceiling face to the ceiling.

    For each over-ceiling face (worst first), move each movable apex the MINIMAL
    distance along the fault normal that brings the face's r_v down to the band
    edge (rv_ceiling*margin) — not all the way to r_v=1, which over-corrects and
    unbalances neighbours sharing the apex.  The move is damped by omega_gs and
    line-searched (bisection) so every incident tet stays valid (sign + eta_ok +
    edge_ok).  GS immediate-update + worst-first + minimal-move avoids the Jacobi
    averaging/cancellation that pins the global sweeps.

    Returns (X_best, n_over_best): the fewest-over-ceiling state visited (faces
    that still cannot move are genuinely inversion-blocked)."""
    apex_a, apex_b, n_f, c_f = topo.apex_a, topo.apex_b, topo.n_f, topo.c_f
    band = rv_ceiling * margin
    X_best = X.copy()
    over_best = int((rv_values(X, topo) > rv_ceiling).sum())
    stall = 0
    for it in range(passes):
        rv = rv_values(X, topo)
        over = np.where(rv > rv_ceiling)[0]
        if over.size == 0:
            break
        over = over[np.argsort(-rv[over])]          # worst first
        for f in over:
            a, b, nrm, cen = apex_a[f], apex_b[f], n_f[f], c_f[f]
            for apex, other in ((a, b), (b, a)):
                if not movable[apex]:
                    continue
                fa = free_axes[apex]
                if not fa.any():
                    continue
                h_self = float((X[apex] - cen) @ nrm)   # fresh (other may have moved)
                h_other = float((X[other] - cen) @ nrm)
                a_self, a_other = abs(h_self), abs(h_other)
                if a_self <= 0.0 or a_other <= 0.0:
                    continue
                if a_self >= a_other:                    # self is the tall side
                    if a_self / a_other <= rv_ceiling:
                        continue
                    tgt = a_other * band                 # shrink toward fault
                else:                                    # self is the short side
                    if a_other / a_self <= rv_ceiling:
                        continue
                    tgt = a_other / band                 # grow away from fault
                sgn = 1.0 if h_self >= 0 else -1.0
                d = omega_gs * (sgn * tgt - h_self) * nrm * fa
                if not np.any(d):
                    continue
                inc = csr_tets[offsets[apex]:offsets[apex + 1]]
                inct = tets[inc]
                s_inc, e_inc, g_inc = sgn0[inc], eta_ok[inc], edge_ok[inc]
                old = X[apex].copy()
                step, ok = 1.0, False
                for _ in range(n_bisect):
                    X[apex] = old + step * d
                    if _incident_ok(X, inct, s_inc, e_inc, g_inc):
                        ok = True
                        break
                    step *= 0.5
                if not ok:
                    X[apex] = old
        n_over = int((rv_values(X, topo) > rv_ceiling).sum())
        if n_over < over_best:
            X_best, over_best, stall = X.copy(), n_over, 0
        else:
            stall += 1
        if verbose:
            print(f"   targeted GS pass {it + 1}: over={n_over} (best {over_best})")
        if over_best == 0 or stall >= 4:
            break
    return X_best, over_best


# --------------------------------------------------------------------------- #
#  The smoother
# --------------------------------------------------------------------------- #
def smooth(coords: np.ndarray, tets: np.ndarray, topo: FaultTopo,
           free_axes: np.ndarray, *, omega: float = 0.5,
           eta_floor: float | None = None, edge_floor: float | None = None,
           rv_ceiling: float | None = None, eta_gate: float = 0.105,
           edge_gate: float = 100.0, relax_factor: float = 0.85,
           stall_patience: int = 6, max_sweeps: int = 60, tol: float = 1e-3,
           max_revert: int = 6, targeted_passes: int = 10,
           verbose: bool = True):
    """Vectorized Jacobi volume-balancing smoothing (see PLAN §Algorithm).

    Moves each node only along its free axes (free_axes[idx]).  A tentative move
    is rejected (the node reverts) if any incident tet would (a) invert, (b) drop
    below the Joe-Liu eta guard, or (c) drop below the min-edge guard.  Both
    guards start at the input mesh's global minimum (eta_min / shortest edge).

    Hard r_v ceiling (rv_ceiling, e.g. 1.5):
      When set, the smoother keeps iterating until NO fault face exceeds the
      ceiling.  If it stalls (the over-ceiling count stops shrinking) while
      faces remain above the ceiling, it RELAXES the quality floors by
      relax_factor — but never below the Q1/Q2 project gates (eta_gate,
      edge_gate) — trading the minimum element quality (down to, but not below,
      the gates) for the r_v ceiling.  The caller checks the returned over-count
      and refuses to write if it is non-zero.

    Returns (coords_best, history, n_over_best): the best (fewest over-ceiling,
    then lowest p99) coordinates found, the per-sweep p99 history, and the
    over-ceiling face count of that best state (0 ⇒ ceiling met; None when no
    ceiling was requested).
    """
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("coords must be (N,3)")
    if free_axes.shape != coords.shape:
        raise ValueError("free_axes must match coords shape")
    X = coords.copy()
    apex_a, apex_b = topo.apex_a, topo.apex_b
    n_f, c_f = topo.n_f, topo.c_f
    movable = free_axes.any(axis=1)          # node can move at all
    mov_a = movable[apex_a]
    mov_b = movable[apex_b]

    V0 = signed_volumes(coords, tets)
    sgn0 = np.sign(V0)
    if np.any(V0 == 0.0):
        raise ValueError("input mesh has zero-volume tet(s); cannot smooth")
    eta0 = tet_etas(coords, tets)
    if eta_floor is None:
        eta_floor = float(eta0.min())
    edge0 = tet_min_edges(coords, tets)
    if edge_floor is None:
        edge_floor = float(edge0.min())
    # a tet is acceptable iff it keeps its sign AND stays >= min(its own
    # baseline, the floor) for BOTH eta and shortest edge.
    eta_ok = np.minimum(eta0, eta_floor)
    edge_ok = np.minimum(edge0, edge_floor)

    N = coords.shape[0]
    history = []

    def metrics(Xc):
        r = np.sort(rv_values(Xc, topo))
        p99 = float(r[int(0.99 * (r.size - 1))])
        mx = float(r[-1])
        n_over = int((r > rv_ceiling).sum()) if rv_ceiling else 0
        return p99, mx, n_over

    p99, mx, n_over = metrics(X)
    history.append(p99)
    X_best, p99_best, n_over_best = X.copy(), p99, n_over
    if verbose:
        tail = f"  over({rv_ceiling})={n_over}" if rv_ceiling else ""
        print(f"  sweep  0: p99 r_v = {p99:.4f}  max = {mx:.4f}  "
              f"(eta_floor={eta_floor:.4f} edge_floor={edge_floor:.1f}){tail}")

    stagnant = 0
    prev_key = n_over if rv_ceiling else p99
    for it in range(1, max_sweeps + 1):
        ha = np.einsum("ij,ij->i", X[apex_a] - c_f, n_f)   # signed
        hb = np.einsum("ij,ij->i", X[apex_b] - c_f, n_f)
        aha, ahb = np.abs(ha), np.abs(hb)

        both = mov_a & mov_b
        only_a = mov_a & ~mov_b
        only_b = mov_b & ~mov_a
        target = np.where(both, 0.5 * (aha + ahb),
                          np.where(only_a, ahb,
                                   np.where(only_b, aha, 0.0)))

        da = (np.sign(ha) * (target - aha))[:, None] * n_f
        db = (np.sign(hb) * (target - ahb))[:, None] * n_f

        disp = np.zeros((N, 3), dtype=np.float64)
        cnt = np.zeros(N, dtype=np.float64)
        np.add.at(disp, apex_a[mov_a], da[mov_a])
        np.add.at(cnt, apex_a[mov_a], 1.0)
        np.add.at(disp, apex_b[mov_b], db[mov_b])
        np.add.at(cnt, apex_b[mov_b], 1.0)
        nz = cnt > 0
        disp[nz] = omega * disp[nz] / cnt[nz, None]
        disp *= free_axes                    # project onto each node's free axes

        X_try = X.copy()
        X_try[nz] += disp[nz]

        # validity guard: revert movable nodes that spoil any incident tet
        for _ in range(max_revert):
            Vn = signed_volumes(X_try, tets)
            etan = tet_etas(X_try, tets)
            edgen = tet_min_edges(X_try, tets)
            bad = (np.sign(Vn) != sgn0) | (etan < eta_ok) | (edgen < edge_ok)
            if not bad.any():
                break
            bad_nodes = np.unique(tets[bad].ravel())
            bad_nodes = bad_nodes[movable[bad_nodes]]
            X_try[bad_nodes] = X[bad_nodes]

        X = X_try
        p99, mx, n_over = metrics(X)
        history.append(p99)
        # keep the best state seen (fewest over-ceiling, then lowest p99)
        if (n_over < n_over_best) or (n_over == n_over_best and p99 < p99_best):
            X_best, p99_best, n_over_best = X.copy(), p99, n_over
        if verbose and (it <= 5 or it % 5 == 0):
            tail = f"  over={n_over}" if rv_ceiling else ""
            print(f"  sweep {it:2d}: p99 r_v = {p99:.4f}  max = {mx:.4f}{tail}")

        if rv_ceiling:
            if n_over == 0:
                if verbose:
                    print(f"  ceiling met at sweep {it}: 0 faces > {rv_ceiling}")
                break
            improved = n_over < prev_key
            stagnant = 0 if improved else stagnant + 1
            prev_key = min(prev_key, n_over)
            if stagnant >= stall_patience:
                can_relax = (eta_floor > eta_gate + 1e-12) or \
                            (edge_floor > edge_gate + 1e-12)
                if can_relax:
                    eta_floor = max(eta_gate, eta_floor * relax_factor)
                    edge_floor = max(edge_gate, edge_floor * relax_factor)
                    eta_ok = np.minimum(eta0, eta_floor)
                    edge_ok = np.minimum(edge0, edge_floor)
                    stagnant = 0
                    if verbose:
                        print(f"  stalled at over={n_over}; relaxing floors -> "
                              f"eta>={eta_floor:.4f} edge>={edge_floor:.1f}")
                else:
                    if verbose:
                        print(f"  stalled at over={n_over} with floors at the "
                              f"Q1/Q2 gates; cannot meet ceiling by smoothing")
                    break
        else:
            if prev_key - p99 < tol:
                stagnant += 1
                if stagnant >= 3:
                    if verbose:
                        print(f"  converged at sweep {it} "
                              f"(Δp99 r_v < {tol} for 3 sweeps)")
                    break
            else:
                stagnant = 0
            prev_key = p99

    # Targeted Gauss-Seidel ceiling enforcement on the residual over-ceiling
    # faces.  The global Jacobi sweeps can stall with faces still over due to
    # averaging/cancellation across an apex's many faces; the GS pass moves each
    # violating face's apex directly (worst-first, per-apex line search) at the
    # most permissive (Q1/Q2 gate) floors.  Faces that still cannot move are
    # genuinely inversion-blocked.
    if rv_ceiling and n_over_best > 0 and targeted_passes > 0:
        if verbose:
            print(f"  global sweeps left {n_over_best} face(s) > {rv_ceiling}; "
                  f"running targeted Gauss-Seidel at the Q1/Q2 gates ...")
        eta_ok_t = np.minimum(eta0, eta_gate)
        edge_ok_t = np.minimum(edge0, edge_gate)
        offsets, csr_tets = build_node_tets_csr(tets, N)
        X_t = X_best.copy()
        X_t, n_over_t = targeted_ceiling_pass(
            X_t, tets, topo, movable, free_axes, sgn0, eta_ok_t, edge_ok_t,
            offsets, csr_tets, rv_ceiling, passes=targeted_passes,
            verbose=verbose)
        if n_over_t < n_over_best:
            X_best, n_over_best = X_t, n_over_t

    return X_best, history, (n_over_best if rv_ceiling else None)


# --------------------------------------------------------------------------- #
#  CLI
# --------------------------------------------------------------------------- #
def _print_stats(label: str, s: dict) -> None:
    print(f"  {label}: n={s['n']:,}  min={s['min']:.4f}  mean={s['mean']:.4f}  "
          f"median={s['median']:.4f}  p95={s['p95']:.4f}  p99={s['p99']:.4f}  "
          f"max={s['max']:.4f}")
    print(f"           r_v>1.5: {s['n_gt15']:,} ({100.0*s['n_gt15']/s['n']:.2f}%)"
          f"   r_v>2.0: {s['n_gt20']:,} ({100.0*s['n_gt20']/s['n']:.2f}%)")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", type=Path, required=True)
    ap.add_argument("--out", dest="out", type=Path, required=True)
    ap.add_argument("--omega", type=float, default=0.5,
                    help="relaxation factor in (0,1] (default 0.5)")
    ap.add_argument("--eta-floor", type=float, default=None,
                    help="min allowed Joe-Liu tet eta (default: the input "
                         "mesh's global eta_min, so Q2 cannot regress)")
    ap.add_argument("--edge-floor", type=float, default=None,
                    help="min allowed tet edge length [m] (default: the input "
                         "mesh's global shortest edge, so Q1 cannot regress)")
    ap.add_argument("--slide-boundary", action="store_true",
                    help="also let domain-boundary nodes slide IN-PLANE on "
                         "their box face (frees boundary-touching tets; the "
                         "domain box is still preserved exactly). Default off: "
                         "boundary frozen (plan invariant).")
    ap.add_argument("--rv-ceiling", type=float, default=None,
                    help="HARD guard: require every fault face r_v <= this "
                         "(e.g. 1.5, the Zhang et al. SSO threshold). The "
                         "smoother adaptively relaxes the eta/edge floors (never "
                         "below the Q1/Q2 gates) until the ceiling is met; if it "
                         "cannot, the tool ERRORS and writes nothing.")
    ap.add_argument("--best-effort", action="store_true",
                    help="with --rv-ceiling: if the ceiling cannot be met, WRITE "
                         "the best (fewest-over) mesh anyway with a warning, "
                         "instead of erroring. The residual over-ceiling faces "
                         "are reported.")
    ap.add_argument("--max-sweeps", type=int, default=60)
    ap.add_argument("--tol", type=float, default=1e-3)
    ap.add_argument("--max-revert", type=int, default=6)
    ap.add_argument("--fault-phys", type=int, default=FAULT_PHYS,
                    help="fault triangle physical tag (SAFS=101, TPV102=3)")
    ap.add_argument("--rock-phys", type=int, default=ROCK_PHYS,
                    help="rock tetra physical tag (default 1)")
    ap.add_argument("--boundary-phys", type=str,
                    default=",".join(str(b) for b in BOUNDARY_PHYS),
                    help="comma-separated boundary triangle tags "
                         "(SAFS='102,103,104', TPV102='1,5')")
    args = ap.parse_args(argv)
    boundary_phys = tuple(int(x) for x in args.boundary_phys.split(","))

    if not args.inp.is_file():
        print(f"ERROR: input mesh not found: {args.inp}", file=sys.stderr)
        return 1
    if not (0.0 < args.omega <= 1.0):
        print("ERROR: --omega must be in (0,1]", file=sys.stderr)
        return 1

    print(f"reading {args.inp} ...")
    g = read_gmsh22(args.inp, fault_phys=args.fault_phys,
                    rock_phys=args.rock_phys, boundary_phys=boundary_phys)
    print(f"  nodes={len(g.node_ids):,}  tets={g.tets.shape[0]:,}  "
          f"fault_tris={g.fault_tris.shape[0]:,}  "
          f"fault_nodes={g.fault_node_idx.size:,}  "
          f"boundary_nodes={g.boundary_node_idx.size:,}")

    topo = build_fault_apex_map(g.coords, g.tets, g.fault_tris)

    free_axes = build_free_axes(g, slide_boundary=args.slide_boundary)
    movable = free_axes.any(axis=1)
    n_apex = np.unique(np.concatenate([topo.apex_a, topo.apex_b]))
    print(f"  movable nodes: {int(movable.sum()):,} "
          f"(near-fault apex movable: {int(movable[n_apex].sum()):,})  "
          f"slide_boundary={args.slide_boundary}")

    print("before:")
    _print_stats("r_v", rv_stats(g.coords, topo))

    ceil_msg = f", rv_ceiling={args.rv_ceiling}" if args.rv_ceiling else ""
    print(f"smoothing (omega={args.omega}, max_sweeps={args.max_sweeps}"
          f"{ceil_msg}) ...")
    coords_new, _, n_over = smooth(
        g.coords, g.tets, topo, free_axes, omega=args.omega,
        eta_floor=args.eta_floor, edge_floor=args.edge_floor,
        rv_ceiling=args.rv_ceiling, max_sweeps=args.max_sweeps,
        tol=args.tol, max_revert=args.max_revert)

    print("after:")
    s_after = rv_stats(coords_new, topo)
    _print_stats("r_v", s_after)

    # r_v ceiling guard.
    if args.rv_ceiling is not None and n_over > 0:
        over = np.where(rv_values(coords_new, topo) > args.rv_ceiling)[0]
        zc = coords_new[g.fault_tris[over]][:, :, 2].mean(axis=1)
        msg = (f"hard r_v ceiling {args.rv_ceiling} NOT met — {n_over} fault "
               f"face(s) remain above it (max={s_after['max']:.4f}, depths "
               f"z=[{zc.min():.0f}, {zc.max():.0f}] m). These are not reducible "
               f"by node smoothing at the Q1/Q2 quality gates (inversion-blocked "
               f"/ fault-spanning tets); local fault re-triangulation is required")
        if not args.best_effort:
            print(f"ERROR: {msg}. Nothing written.", file=sys.stderr)
            return 3
        print(f"WARNING: {msg}. Writing best-effort mesh anyway "
              f"(--best-effort).", file=sys.stderr)

    moved_mask = np.linalg.norm(coords_new - g.coords, axis=1) > 0.0
    n_moved = int(moved_mask.sum())
    # invariant guard: fault nodes must NEVER move; locked boundary axes (the
    # axis normal to each box face) must be byte-unchanged.
    if np.any(moved_mask[g.fault_node_idx]):
        print("ERROR: a fault node moved — aborting write", file=sys.stderr)
        return 2
    locked = ~free_axes
    if np.any((coords_new[locked] != g.coords[locked])):
        print("ERROR: a locked (frozen-axis) coordinate moved — aborting write",
              file=sys.stderr)
        return 2
    print(f"moved {n_moved:,} nodes "
          f"({100.0*n_moved/len(g.node_ids):.1f}% of nodes)")

    write_gmsh22(args.out, g, coords_new, moved_mask)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
