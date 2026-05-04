"""Cross-fault triangle-triangle intersection geometry.

Implements PLAN_multifault_intersections.md Phase 1.

Public API:
    Segment, Polyline       — output dataclasses
    tri_tri_intersect_3d    — Möller's interval-overlap method on every
                              triangle-pair (i, j) ∈ T_A × T_B
    chain_segments          — connect end-to-end-coincident segments into
                              polylines
    scan_cross_fault_crossings  — one-shot fault-pair scan from STL files

Pure numpy.  ``gmpy2`` is an optional dependency used only for exact-
arithmetic resolution of ambiguous orientation predicates; if it is
unavailable, ambiguous cases are logged and treated as non-crossing
(no false positives, at most a missed near-tangency — acceptable for
non-adversarial CFM input).

Coordinate frame: SAFS local-Cartesian metres (per ``safs_origin.py``).
"""
from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np

try:
    from gmpy2 import mpq                                # noqa: F401
    _HAS_MPQ = True
except ImportError:
    _HAS_MPQ = False

_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public dataclasses
# ---------------------------------------------------------------------------
@dataclass
class Segment:
    """One triangle-triangle intersection segment in 3-D.

    Only segments with length >= eps_min_seg_len_m are emitted; the
    plan's earlier ``nondegenerate`` field has been removed because
    every emitted segment satisfied it (degenerates are filtered out
    upstream rather than emitted with a flag).
    """
    p0: tuple[float, float, float]
    p1: tuple[float, float, float]
    tri_a: int
    tri_b: int
    fault_a: str
    fault_b: str


@dataclass
class Polyline:
    """A chained sequence of intersection segments lying on two faults."""
    points: list[tuple[float, float, float]]
    closed: bool
    fault_a: str
    fault_b: str
    pierces_a: list[int] = field(default_factory=list)
    pierces_b: list[int] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Plane primitives — use UNIT normals so signed distances have units of
# metres; the eps_orient × mean_edge threshold is then dimensionally
# consistent (length × length-fraction = length).
# ---------------------------------------------------------------------------
def _plane_unit(V_tri: np.ndarray) -> tuple[np.ndarray | None, float]:
    """Return (n_unit, d) for plane n_unit · x + d = 0; ``None`` if the
    triangle is degenerate (zero area)."""
    n = np.cross(V_tri[1] - V_tri[0], V_tri[2] - V_tri[0])
    n_norm = float(np.linalg.norm(n))
    if n_norm == 0.0:
        return None, 0.0
    n_unit = n / n_norm
    d = -float(n_unit @ V_tri[0])
    return n_unit, d


def _mean_edge(V_tri: np.ndarray) -> float:
    e0 = float(np.linalg.norm(V_tri[1] - V_tri[0]))
    e1 = float(np.linalg.norm(V_tri[2] - V_tri[1]))
    e2 = float(np.linalg.norm(V_tri[0] - V_tri[2]))
    return (e0 + e1 + e2) / 3.0


def _categorise(d: np.ndarray, thresh: float) -> np.ndarray:
    """Map a signed distance vector to {-1, 0, +1} with band ±thresh
    treated as zero."""
    s = np.zeros_like(d, dtype=np.int8)
    s[d > thresh] = 1
    s[d < -thresh] = -1
    return s


# ---------------------------------------------------------------------------
# AABB pre-filter (vectorised).  Reduces candidate pairs from n_A × n_B
# to those whose axis-aligned bounding boxes overlap.
# ---------------------------------------------------------------------------
def _aabb_overlap(V_A: np.ndarray, T_A: np.ndarray,
                  V_B: np.ndarray, T_B: np.ndarray) -> np.ndarray:
    tri_A = V_A[T_A]                            # (n_A, 3, 3)
    tri_B = V_B[T_B]                            # (n_B, 3, 3)
    A_min = tri_A.min(axis=1)
    A_max = tri_A.max(axis=1)
    B_min = tri_B.min(axis=1)
    B_max = tri_B.max(axis=1)
    return (
        (A_min[:, None, :] <= B_max[None, :, :]).all(axis=-1)
        & (B_min[None, :, :] <= A_max[:, None, :]).all(axis=-1)
    )


# ---------------------------------------------------------------------------
# Möller core
# ---------------------------------------------------------------------------
def _interval_on_L(
    V_tri: np.ndarray,
    d_signed: np.ndarray,
    D: np.ndarray,
    thresh_zero: float,
) -> tuple[float, float, np.ndarray, np.ndarray] | None:
    """Compute (t_lo, t_hi, p_lo, p_hi) — the triangle's interval on the
    intersection line L plus the two 3-D ``p_hit`` points.  Returns
    ``None`` if the geometry is degenerate (e.g., an edge lying entirely
    in the other plane — see Phase 1 §Edge cases).

    ``thresh_zero`` is the signed-distance band inside which a vertex is
    treated as "on the other plane".
    """
    edges = ((0, 1), (1, 2), (2, 0))
    p_hits: list[np.ndarray] = []
    ts: list[float] = []
    edge_on_plane = False
    for i, j in edges:
        di, dj = float(d_signed[i]), float(d_signed[j])
        # Both endpoints on the plane → 1-D-intersection edge case.
        if abs(di) <= thresh_zero and abs(dj) <= thresh_zero:
            edge_on_plane = True
            continue
        # Same-strict-sign → does not cross.
        if di > thresh_zero and dj > thresh_zero:
            continue
        if di < -thresh_zero and dj < -thresh_zero:
            continue
        denom = di - dj
        if denom == 0.0:
            return None        # numerical breakdown, abandon this pair
        s = di / denom
        s = 0.0 if s < 0.0 else (1.0 if s > 1.0 else s)
        p_hit = V_tri[i] + s * (V_tri[j] - V_tri[i])
        p_hits.append(p_hit)
        ts.append(float(D @ p_hit))

    if edge_on_plane:
        # Positive-length 1-D intersection (B's edge in A's plane or
        # vice versa).  See §Edge cases — proper handling requires
        # clipping the on-plane edge against the other triangle's
        # interior.  For Phase 1 we log and skip; the cascade phase
        # should not produce this in practice for distinct CFM faults.
        # NOTE: in the gmpy2-absent default configuration the caller
        # `tri_tri_intersect_3d` runs the ambiguity check before
        # invoking `_interval_on_L`; a "vertex within thresh of plane"
        # triggers ambiguous-skip there, so this branch is only
        # reachable when gmpy2 is installed (and resolves to "intersect"
        # for an exactly-on-plane edge).
        _log.warning(
            "edge-on-plane intersection encountered (1-D positive-length "
            "case); not yet implemented, treating as non-crossing"
        )
        return None

    # We expect exactly two crossing edges.  If a vertex sat exactly on
    # the plane (within thresh_zero), three edges may register; collapse
    # to the two with the most-distinct t-values.
    if len(p_hits) < 2:
        return None
    if len(p_hits) > 2:
        order = sorted(range(len(ts)), key=lambda k: ts[k])
        lo, hi = order[0], order[-1]
        p_hits = [p_hits[lo], p_hits[hi]]
        ts = [ts[lo], ts[hi]]

    if ts[0] <= ts[1]:
        return ts[0], ts[1], p_hits[0], p_hits[1]
    return ts[1], ts[0], p_hits[1], p_hits[0]


def _interp_on_chord(
    t: float, t_lo: float, t_hi: float,
    p_lo: np.ndarray, p_hi: np.ndarray,
) -> np.ndarray:
    """Linear interpolation along a triangle's chord on L (PLAN
    Phase 1 step (6))."""
    span = t_hi - t_lo
    if span == 0.0:
        return p_lo.copy()
    alpha = (t - t_lo) / span
    if alpha < 0.0:
        alpha = 0.0
    elif alpha > 1.0:
        alpha = 1.0
    return (1.0 - alpha) * p_lo + alpha * p_hi


def _free_surface_clamp(
    p0: np.ndarray, p1: np.ndarray, clearance_m: float,
) -> tuple[np.ndarray, np.ndarray] | None:
    """PLAN Phase 1 step (8): if either endpoint has z > -clearance_m,
    walk it down the segment to z = -clearance_m.  Drop the segment
    only if it lies wholly above the clearance plane."""
    if clearance_m <= 0.0:
        return p0.copy(), p1.copy()
    z_thresh = -clearance_m
    z0 = float(p0[2])
    z1 = float(p1[2])
    above_0 = z0 > z_thresh
    above_1 = z1 > z_thresh
    if above_0 and above_1:
        return None                                  # wholly above; drop
    if not above_0 and not above_1:
        return p0.copy(), p1.copy()
    out0 = p0.copy()
    out1 = p1.copy()
    if above_0:
        dz = z1 - z0
        if dz == 0.0:
            return None
        alpha = (z_thresh - z0) / dz
        alpha = max(0.0, min(1.0, alpha))
        out0 = p0 + alpha * (p1 - p0)
    if above_1:
        dz = z0 - z1
        if dz == 0.0:
            return None
        alpha = (z_thresh - z1) / dz
        alpha = max(0.0, min(1.0, alpha))
        out1 = p1 + alpha * (p0 - p1)
    return out0, out1


# ---------------------------------------------------------------------------
# Optional gmpy2 fallback for ambiguous predicates
# ---------------------------------------------------------------------------
def _resolve_with_mpq(TriA: np.ndarray, TriB: np.ndarray) -> int:
    """Exact-rational resolution of the same-side / mixed-side test.
    Returns -1 (proven no intersection), +1 (proven mixed), 0 (still
    ambiguous — only happens for genuinely-coplanar input)."""
    from gmpy2 import mpq                      # type: ignore[import-not-found]
    SCALE = 1_000_000

    def to_mpq(p):
        return tuple(mpq(int(round(float(c) * SCALE)), SCALE) for c in p)

    A0, A1, A2 = (to_mpq(p) for p in TriA)
    B0, B1, B2 = (to_mpq(p) for p in TriB)

    def sub(p, q):
        return tuple(pi - qi for pi, qi in zip(p, q))

    def cross(u, v):
        return (u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0])

    def dot(u, v):
        return sum(ui * vi for ui, vi in zip(u, v))

    n_B = cross(sub(B1, B0), sub(B2, B0))
    d_B_off = -dot(n_B, B0)
    d_A = [dot(n_B, A) + d_B_off for A in (A0, A1, A2)]

    n_A = cross(sub(A1, A0), sub(A2, A0))
    d_A_off = -dot(n_A, A0)
    d_B_local = [dot(n_A, B) + d_A_off for B in (B0, B1, B2)]

    def all_strict_pos(xs):
        return all(x > 0 for x in xs)

    def all_strict_neg(xs):
        return all(x < 0 for x in xs)

    if all_strict_pos(d_A) or all_strict_neg(d_A):
        return -1
    if all_strict_pos(d_B_local) or all_strict_neg(d_B_local):
        return -1
    return +1


# ---------------------------------------------------------------------------
# Public: tri_tri_intersect_3d
# ---------------------------------------------------------------------------
def tri_tri_intersect_3d(
    V_A: np.ndarray, T_A: np.ndarray,
    V_B: np.ndarray, T_B: np.ndarray,
    fault_a: str, fault_b: str,
    eps_orient: float = 1e-9,
    eps_min_seg_len_m: float = 5e-2,
    clearance_m: float = 0.0,
    *,
    coplanar_dot_thresh: float = 1.0 - 1e-6,
    coplanar_offset_thresh: float = 1.0,
) -> list[Segment]:
    """Möller's interval-overlap intersection over every (i, j) ∈
    T_A × T_B (filtered by AABB).  See PLAN Phase 1 for the algorithm.

    Parameters
    ----------
    V_A, T_A : per-fault-A vertex array (N, 3) and triangle index array
        (M, 3); int indices into V_A.  V coordinates are in metres
        (SAFS local frame).
    V_B, T_B : same for fault B.
    fault_a, fault_b : short-name strings stamped onto every emitted
        Segment for downstream attribution.
    eps_orient : dimensionless tolerance applied as
        ``eps_orient * mean_edge`` to signed distances (which are in
        metres because we use unit-normalised plane normals).  Vertices
        whose signed distance falls inside this band are treated as
        on-plane and trigger the gmpy2 fallback (or a non-crossing
        decision if gmpy2 is absent).
    eps_min_seg_len_m : segments shorter than this are dropped.
        Default 5 cm — must satisfy ``eps_min_seg_len_m > 2 * snap_m``
        so that ``chain_segments`` can match endpoints without
        producing self-loops.
    clearance_m : free-surface clamp distance in metres.  Endpoints
        with ``z > -clearance_m`` are linearly interpolated down to
        ``z = -clearance_m``.  Default 0 disables the clamp.
    coplanar_dot_thresh, coplanar_offset_thresh : guards for the
        coplanar-input case (unexpected on CFM data); raise loudly
        rather than silently mishandle.
    """
    V_A = np.asarray(V_A, dtype=np.float64)
    V_B = np.asarray(V_B, dtype=np.float64)
    T_A = np.asarray(T_A, dtype=np.int64)
    T_B = np.asarray(T_B, dtype=np.int64)
    if V_A.ndim != 2 or V_A.shape[1] != 3:
        raise ValueError(f"V_A shape {V_A.shape} expected (N, 3)")
    if V_B.ndim != 2 or V_B.shape[1] != 3:
        raise ValueError(f"V_B shape {V_B.shape} expected (N, 3)")
    if T_A.ndim != 2 or T_A.shape[1] != 3:
        raise ValueError(f"T_A shape {T_A.shape} expected (M, 3)")
    if T_B.ndim != 2 or T_B.shape[1] != 3:
        raise ValueError(f"T_B shape {T_B.shape} expected (M, 3)")
    if eps_orient <= 0.0:
        raise ValueError(f"eps_orient must be > 0; got {eps_orient}")
    if eps_min_seg_len_m <= 0.0:
        raise ValueError(
            f"eps_min_seg_len_m must be > 0; got {eps_min_seg_len_m}"
        )
    if clearance_m < 0.0:
        raise ValueError(f"clearance_m must be >= 0; got {clearance_m}")

    if T_A.shape[0] == 0 or T_B.shape[0] == 0:
        return []

    overlap = _aabb_overlap(V_A, T_A, V_B, T_B)
    candidate_pairs = np.argwhere(overlap)

    out: list[Segment] = []
    for i_arr, j_arr in candidate_pairs:
        i = int(i_arr)
        j = int(j_arr)
        TriA = V_A[T_A[i]]
        TriB = V_B[T_B[j]]

        n_A_unit, d_A_off = _plane_unit(TriA)
        n_B_unit, d_B_off = _plane_unit(TriB)
        if n_A_unit is None or n_B_unit is None:
            continue                                  # zero-area triangle

        cos_dh = abs(float(n_A_unit @ n_B_unit))
        if cos_dh > coplanar_dot_thresh:
            sign = math.copysign(1.0, float(n_A_unit @ n_B_unit))
            offset = abs(d_A_off - sign * d_B_off)
            if offset < coplanar_offset_thresh:
                # The plan claims this should not occur ("different
                # geological structures"); empirically it does, on the
                # SAFS CFM all-8 fixture, where adjacent fault segments
                # (e.g. Mojave-SAF and SBMT-SAF) share triangulation
                # boundaries and produce bit-exactly coplanar triangle
                # pairs.  Coplanar pairs share a 2-D region, not a 1-D
                # crossing curve, so emitting no segment is the correct
                # geometric answer.  Logged at WARNING level so the
                # caller knows it happened.  (See implementation
                # report — deviation from the plan's "raise" wording.)
                _log.warning(
                    "coplanar triangle pair %s.tri_a=%d × %s.tri_b=%d "
                    "(|n_A·n_B|=%.6g, offset=%.6g m); skipping — no "
                    "1-D intersection segment exists for coplanar pairs.",
                    fault_a, i, fault_b, j, cos_dh, offset,
                )
                continue
            # Parallel non-coplanar: same-side check below catches it.

        d_A_to_B = TriA @ n_B_unit + d_B_off          # length, signed
        d_B_to_A = TriB @ n_A_unit + d_A_off

        thresh_A = eps_orient * _mean_edge(TriA)
        thresh_B = eps_orient * _mean_edge(TriB)

        s_A = _categorise(d_A_to_B, thresh_A)
        s_B = _categorise(d_B_to_A, thresh_B)

        # Quick reject: A wholly on one strict side of plane B?
        if (s_A > 0).all() or (s_A < 0).all():
            continue
        if (s_B > 0).all() or (s_B < 0).all():
            continue

        # Ambiguity: at least one signed distance fell inside the
        # closed band [-thresh, +thresh].  Using `<=` here matches
        # `_categorise`'s strict `>`/`<` (a vertex with |d| == thresh
        # is in-band per `_categorise`, so it must also be flagged
        # ambiguous here for the two checks to agree).
        ambiguous = (
            (np.abs(d_A_to_B) <= thresh_A).any()
            or (np.abs(d_B_to_A) <= thresh_B).any()
        )
        if ambiguous:
            if _HAS_MPQ:
                if _resolve_with_mpq(TriA, TriB) <= 0:
                    continue
            else:
                _log.warning(
                    "ambiguous orientation predicate at "
                    "%s.tri_a=%d × %s.tri_b=%d; gmpy2 not available, "
                    "treating as non-crossing",
                    fault_a, i, fault_b, j,
                )
                continue

        D = np.cross(n_A_unit, n_B_unit)
        D_norm = float(np.linalg.norm(D))
        if D_norm == 0.0:
            continue
        D = D / D_norm

        intervalA = _interval_on_L(TriA, d_A_to_B, D, thresh_A)
        if intervalA is None:
            continue
        intervalB = _interval_on_L(TriB, d_B_to_A, D, thresh_B)
        if intervalB is None:
            continue

        tA_lo, tA_hi, pA_lo, pA_hi = intervalA
        tB_lo, tB_hi, _pB_lo, _pB_hi = intervalB

        t_start = max(tA_lo, tB_lo)
        t_end = min(tA_hi, tB_hi)
        if t_start > t_end:
            continue

        p_start = _interp_on_chord(t_start, tA_lo, tA_hi, pA_lo, pA_hi)
        p_end = _interp_on_chord(t_end,   tA_lo, tA_hi, pA_lo, pA_hi)

        # PLAN Phase 1 step (7) — length filter — runs BEFORE step (8)
        # free-surface clamp.  Drop pre-clamp degenerate segments here;
        # the clamp may further shorten what's left, so apply a second
        # length check after the clamp to drop post-clamp degenerates.
        seg_len = float(np.linalg.norm(p_end - p_start))
        if seg_len < eps_min_seg_len_m:
            continue

        if clearance_m > 0.0:
            clamped = _free_surface_clamp(p_start, p_end, clearance_m)
            if clamped is None:
                continue
            p_start, p_end = clamped
            # Re-check length after clamp: a long pre-clamp segment can
            # become a sub-eps post-clamp sliver when most of it sits
            # above the free surface.
            if float(np.linalg.norm(p_end - p_start)) < eps_min_seg_len_m:
                continue

        out.append(Segment(
            p0=tuple(float(x) for x in p_start),
            p1=tuple(float(x) for x in p_end),
            tri_a=i, tri_b=j,
            fault_a=fault_a, fault_b=fault_b,
        ))

    return out


# ---------------------------------------------------------------------------
# Public: chain_segments
# ---------------------------------------------------------------------------
def chain_segments(
    segments: list[Segment],
    snap_m: float = 1e-2,
) -> list[Polyline]:
    """Chain end-to-end-coincident segments into polylines.

    Pre-condition: every input segment satisfies ``|p1 - p0| > 2 *
    snap_m`` (asserted on entry; the default ``eps_min_seg_len_m =
    5e-2`` in :func:`tri_tri_intersect_3d` keeps this invariant).
    Otherwise the chain step would collapse a segment to a self-loop
    after endpoint snapping, which downstream consumers reject.
    """
    if snap_m <= 0.0:
        raise ValueError(f"snap_m must be > 0; got {snap_m}")
    if not segments:
        return []

    fault_a = segments[0].fault_a
    fault_b = segments[0].fault_b
    for s in segments:
        if s.fault_a != fault_a or s.fault_b != fault_b:
            raise ValueError(
                f"chain_segments: mixed fault pairs in input "
                f"({fault_a, fault_b} vs {s.fault_a, s.fault_b})"
            )
        L = math.sqrt(
            (s.p0[0] - s.p1[0]) ** 2
            + (s.p0[1] - s.p1[1]) ** 2
            + (s.p0[2] - s.p1[2]) ** 2
        )
        assert L > 2.0 * snap_m, (
            f"chain_segments precondition violated: segment length "
            f"{L:.6e} m is not > 2 * snap_m = {2.0 * snap_m:.6e} m. "
            f"Filter shorter segments upstream "
            f"(eps_min_seg_len_m in tri_tri_intersect_3d)."
        )

    def key(p):
        return (
            int(round(p[0] / snap_m)),
            int(round(p[1] / snap_m)),
            int(round(p[2] / snap_m)),
        )

    node_to_id: dict[tuple[int, int, int], int] = {}
    nodes: list[tuple[float, float, float]] = []
    seg_endpoints: list[tuple[int, int]] = []
    for s in segments:
        for p in (s.p0, s.p1):
            k = key(p)
            if k not in node_to_id:
                node_to_id[k] = len(nodes)
                nodes.append(tuple(float(c) for c in p))
        seg_endpoints.append((node_to_id[key(s.p0)], node_to_id[key(s.p1)]))

    adj: list[list[tuple[int, int]]] = [[] for _ in nodes]
    for sidx, (u, v) in enumerate(seg_endpoints):
        adj[u].append((sidx, v))
        adj[v].append((sidx, u))

    polylines: list[Polyline] = []
    visited_segs: set[int] = set()

    for sidx_root in range(len(segments)):
        if sidx_root in visited_segs:
            continue
        # BFS to discover the connected component (set of segment indices).
        component: set[int] = set()
        stack = [sidx_root]
        while stack:
            sidx = stack.pop()
            if sidx in component:
                continue
            component.add(sidx)
            u, v = seg_endpoints[sidx]
            for n in (u, v):
                for nbr_sidx, _ in adj[n]:
                    if nbr_sidx not in component:
                        stack.append(nbr_sidx)

        # A component may have branches (degree-3+ junction nodes).  Each
        # outer iteration of this while-loop walks ONE chain from a leaf
        # through the component (or one closed loop), emits it as a
        # polyline, and marks the walked segments as visited.  Subsequent
        # outer-iterations of the enclosing for-loop pick up unvisited
        # branches as additional polylines.  This is the F-001 fix:
        # `visited_segs |= component` (the previous version) silently
        # dropped branches; we now mark only what was actually walked.
        while True:
            unvisited_in_component = component - visited_segs
            if not unvisited_in_component:
                break

            # Endpoint search: prefer a leaf in the *unvisited remainder*
            # of the component.  Counting only unvisited adjacencies
            # ensures we restart at the next branch's free end after a
            # previous walk consumed an arm of a Y or T.
            endpoint_node: int | None = None
            for sidx in unvisited_in_component:
                for node_id in seg_endpoints[sidx]:
                    deg_unvisited = sum(
                        1 for (sx, _) in adj[node_id]
                        if sx in component and sx not in visited_segs
                    )
                    if deg_unvisited == 1:
                        endpoint_node = node_id
                        break
                if endpoint_node is not None:
                    break

            is_closed = endpoint_node is None
            if is_closed:
                sample_sidx = next(iter(unvisited_in_component))
                endpoint_node = seg_endpoints[sample_sidx][0]

            ordered_pts: list[tuple[float, float, float]] = [
                nodes[endpoint_node]
            ]
            ordered_pierces_a: list[int] = []
            ordered_pierces_b: list[int] = []
            used: set[int] = set()
            cur = endpoint_node
            while True:
                next_sidx: int | None = None
                next_node: int | None = None
                for sidx, other in adj[cur]:
                    if (sidx in component
                            and sidx not in used
                            and sidx not in visited_segs):
                        next_sidx = sidx
                        next_node = other
                        break
                if next_sidx is None:
                    break
                used.add(next_sidx)
                ordered_pts.append(nodes[next_node])
                ordered_pierces_a.append(segments[next_sidx].tri_a)
                ordered_pierces_b.append(segments[next_sidx].tri_b)
                cur = next_node
                if is_closed and cur == endpoint_node:
                    break

            visited_segs |= used
            if not used:
                # Defensive: a component with no walkable segments
                # (shouldn't happen given the unvisited_in_component
                # check above, but guard against an infinite loop in
                # case of a cycle of zero-edge graph anomalies).
                break
            polylines.append(Polyline(
                points=ordered_pts,
                closed=is_closed,
                fault_a=fault_a,
                fault_b=fault_b,
                # One entry per polyline segment (NOT deduplicated):
                # pierces_a[i] is the fault_a triangle that segment i
                # (between points[i] and points[i+1]) lies in.  Phase 2
                # relies on this 1:1 segment-to-triangle mapping; if a
                # polyline visits the same triangle twice (loops back),
                # the index appears twice — that's correct.
                pierces_a=ordered_pierces_a,
                pierces_b=ordered_pierces_b,
            ))

    return polylines


# ---------------------------------------------------------------------------
# Verification helper: minimum 3-D distance from a point to a triangle set.
# Used by the all-8 polylines-on-source-surfaces sanity check
# (test_F010_*).  Pure numpy, vectorised over the triangle set; called
# once per polyline vertex (~hundreds of vertices on the all-8 fixture).
# ---------------------------------------------------------------------------
def _distance_point_to_segment(p: np.ndarray, a: np.ndarray,
                                b: np.ndarray) -> np.ndarray:
    """Vectorised perpendicular distance from point ``p`` (3,) to each
    segment ``a_i — b_i`` in arrays a, b of shape (M, 3).  Returns
    array of shape (M,)."""
    ab = b - a
    ap = p - a
    ab_dot_ab = (ab * ab).sum(axis=-1)
    # Avoid division-by-zero on degenerate segments (a == b); for those
    # we fall back to distance to the endpoint.
    safe = np.where(ab_dot_ab > 0.0, ab_dot_ab, 1.0)
    t = (ap * ab).sum(axis=-1) / safe
    t = np.clip(t, 0.0, 1.0)
    closest = a + t[:, None] * ab
    return np.linalg.norm(p - closest, axis=-1)


def _distance_point_to_triangle_set(
    p: np.ndarray, V: np.ndarray, T: np.ndarray,
) -> float:
    """Minimum 3-D Euclidean distance from point ``p`` (3,) to any
    triangle in (V, T).

    Algorithm: project ``p`` onto each triangle's plane; if the
    projection's barycentric coordinates are all >= 0 (inside the
    triangle), use the perpendicular distance; otherwise fall back to
    the minimum edge-clamped distance (the projection is closest to a
    vertex or edge of the triangle).  Vectorised over triangles.

    Degenerate triangles (zero area) contribute their edge-distance
    only.  The function returns ``inf`` for an empty triangle set.
    """
    p = np.asarray(p, dtype=np.float64).reshape(3)
    V = np.asarray(V, dtype=np.float64)
    T = np.asarray(T, dtype=np.int64)
    if T.shape[0] == 0:
        return float("inf")

    A = V[T[:, 0]]                              # (M, 3)
    B = V[T[:, 1]]
    C = V[T[:, 2]]

    # 1. Foot of perpendicular onto each triangle's plane.
    n = np.cross(B - A, C - A)                  # (M, 3) unnormalised
    n_norm_sq = (n * n).sum(axis=-1)            # (M,)
    valid_plane = n_norm_sq > 0.0
    safe_n = np.where(valid_plane[:, None],
                       n / np.sqrt(np.where(valid_plane, n_norm_sq, 1.0)
                                    )[:, None],
                       0.0)
    d_signed = ((p - A) * safe_n).sum(axis=-1)  # (M,)
    p_proj = p - d_signed[:, None] * safe_n     # (M, 3)

    # 2. Barycentric coords of p_proj relative to (A, B, C).
    # Solve [e1 e2] [u v]^T = p_proj - A where e1 = B - A, e2 = C - A.
    e1 = B - A
    e2 = C - A
    d11 = (e1 * e1).sum(axis=-1)
    d12 = (e1 * e2).sum(axis=-1)
    d22 = (e2 * e2).sum(axis=-1)
    rhs = p_proj - A
    r1 = (e1 * rhs).sum(axis=-1)
    r2 = (e2 * rhs).sum(axis=-1)
    det = d11 * d22 - d12 * d12
    safe_det = np.where(det > 0.0, det, 1.0)
    u = (d22 * r1 - d12 * r2) / safe_det
    v = (d11 * r2 - d12 * r1) / safe_det
    w = 1.0 - u - v

    # Inside iff all bary coords >= 0 (within FP tolerance) AND the
    # plane is well-defined.
    inside = (u >= -1e-12) & (v >= -1e-12) & (w >= -1e-12) & valid_plane

    # 3. Edge-clamped distances (always valid, even for degenerate tris).
    d_AB = _distance_point_to_segment(p, A, B)
    d_BC = _distance_point_to_segment(p, B, C)
    d_CA = _distance_point_to_segment(p, C, A)
    d_edges = np.minimum(np.minimum(d_AB, d_BC), d_CA)

    # 4. Pick perpendicular distance when inside, edge distance otherwise.
    d_perp = np.abs(d_signed)
    d_per_tri = np.where(inside, d_perp, d_edges)

    return float(d_per_tri.min())


# ---------------------------------------------------------------------------
# Public: scan_cross_fault_crossings
# ---------------------------------------------------------------------------
def _read_stl_arrays(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read an ASCII STL into (V, T) arrays.  Each triangle keeps its own
    three vertex slots — Phase 1 doesn't need a deduplicated mesh."""
    verts: list[tuple[float, float, float]] = []
    tris: list[tuple[int, int, int]] = []
    cur: list[tuple[float, float, float]] = []
    with open(path, "r") as fh:
        for raw in fh:
            tok = raw.strip().split()
            if len(tok) >= 4 and tok[0] == "vertex":
                cur.append((float(tok[1]), float(tok[2]), float(tok[3])))
            elif tok and tok[0] == "endloop":
                if len(cur) != 3:
                    raise ValueError(
                        f"{path}: outer loop with {len(cur)} vertices"
                    )
                a = len(verts)
                verts.extend(cur)
                tris.append((a, a + 1, a + 2))
                cur = []
    if not tris:
        raise ValueError(f"{path}: no triangles found")
    return (np.asarray(verts, dtype=np.float64),
            np.asarray(tris, dtype=np.int64))


def simplify_polyline(
    polyline: "Polyline",
    *,
    min_pierce_separation_m: float,
) -> "Polyline":
    """Merge consecutive polyline vertices closer than
    ``min_pierce_separation_m``.

    Polyline pierces produced by Phase 1 can be only metres apart when
    the intersection curve kinks at every triangle-edge crossing on the
    other fault.  Inserting both as Steiner vertices in Phase 2 produces
    sliver triangles with aspect ratios up to ~1000.  Merging close
    pierces removes these slivers; the cost is up to
    ``min_pierce_separation_m / 2`` of geometric precision in the
    intersection curve location.

    The simplification is applied in 3-D (not in either fault's local
    frame), and the merged point is the midpoint of the two collapsed
    vertices.  Both faults that share the polyline see the same
    simplified curve, so cross-fault polyline-edge coincidence is
    preserved.

    Returns a new Polyline with ``points``, ``pierces_a``, and
    ``pierces_b`` all consistently shortened.  The original is
    unmodified.

    Edge cases:
        - A polyline with ≤ 2 vertices (one segment) is never
          simplified.
        - A closed polyline (first == last) keeps that closure
          property.
        - When merging polyline vertex i+1 into vertex i, the segment
          from i to i+1 disappears, so ``pierces_a[i]`` and
          ``pierces_a[i+1]`` collapse: the merged vertex is shared by
          the segments that were on either side, and the merged
          ``pierces_a`` retains both attribution entries (the polyline
          still crosses both triangles, just with the kink at the
          midpoint instead of at the original vertex).  This matches
          the post-merge geometry: the polyline goes
          ``points[i-1] → merged → points[i+2]`` with two segments,
          in triangles ``pierces_a[i-1]`` and ``pierces_a[i+1]``
          respectively (i.e., ``pierces_a[i]`` is the dropped one).
    """
    if min_pierce_separation_m <= 0.0:
        return polyline
    pts = polyline.points
    if len(pts) <= 2:
        return polyline

    new_pts: list[tuple[float, float, float]] = [pts[0]]
    keep_pierces_a: list[int] = []
    keep_pierces_b: list[int] = []
    n = len(pts)
    n_seg = len(polyline.pierces_a)
    if n_seg != n - 1:
        # Should not happen with the post-Phase-1-fix Polyline, but guard.
        return polyline

    last = np.asarray(pts[0], dtype=np.float64)
    pending_seg_a = polyline.pierces_a[0]
    pending_seg_b = polyline.pierces_b[0]
    for i in range(1, n):
        cur = np.asarray(pts[i], dtype=np.float64)
        d = float(np.linalg.norm(cur - last))
        # Only consider merging when both segments meeting at pt[i]
        # are in the SAME (tri_a, tri_b) pair (i.e., the kink is
        # within a single triangle on each fault).  Merging across
        # triangle boundaries would re-attribute the merged segment to
        # only one of the two crossed triangles, leaving the other
        # without its share of the polyline as a mesh edge → T-junction.
        same_pair = False
        if i < n_seg:
            same_pair = (
                polyline.pierces_a[i - 1] == polyline.pierces_a[i]
                and polyline.pierces_b[i - 1] == polyline.pierces_b[i]
            )
        if d < min_pierce_separation_m and same_pair:
            # Merge cur into the last kept point.  The segment ending at
            # cur is dropped; pending_seg_* (the segment STARTING at
            # cur in the original polyline) carries forward as the
            # attribution for the next-emitted segment.
            mid = 0.5 * (last + cur)
            new_pts[-1] = (float(mid[0]), float(mid[1]), float(mid[2]))
            last = mid
            pending_seg_a = polyline.pierces_a[i]
            pending_seg_b = polyline.pierces_b[i]
        else:
            new_pts.append(pts[i])
            keep_pierces_a.append(pending_seg_a)
            keep_pierces_b.append(pending_seg_b)
            last = cur
            if i < n_seg:
                pending_seg_a = polyline.pierces_a[i]
                pending_seg_b = polyline.pierces_b[i]

    # Re-check the closure case: a closed polyline has points[0] == points[-1].
    # The merge logic above operates on consecutive pairs and respects
    # that invariant naturally.
    if not keep_pierces_a:
        # Whole polyline collapsed to one point — degenerate; return original.
        return polyline

    return Polyline(
        points=new_pts,
        closed=polyline.closed,
        fault_a=polyline.fault_a,
        fault_b=polyline.fault_b,
        pierces_a=keep_pierces_a,
        pierces_b=keep_pierces_b,
    )


def refine_uniform_midpoint(
    V: np.ndarray, T: np.ndarray,
    target_edge_length_m: float,
    *,
    max_rounds: int = 8,
) -> tuple[np.ndarray, np.ndarray]:
    """Refine a triangulated surface to ~uniform edge length via the
    standard 1-to-4 midpoint subdivision.

    Each iteration: every triangle that has at least one edge longer
    than ``target_edge_length_m`` is split into 4 children by inserting
    midpoints on each of its edges.  Edge midpoints are shared between
    adjacent triangles by construction, so the refinement is
    consistent (no T-junctions) without any explicit propagation step.

    Triangles whose max edge is already ≤ target are left unchanged
    in that round.  This produces a non-uniform mesh in the
    intermediate state, which is then smoothed by the next round.
    The loop terminates when every triangle has max-edge ≤ target,
    or after ``max_rounds`` (safeguard against degenerate input).

    Returns (V_new, T_new).  Original V and T are not modified.

    Note: a triangle that has only ONE long edge gets all 3 edges
    midpoint-split anyway (the 1-to-4 scheme).  This over-refines
    modestly but preserves shape quality (all children are similar
    to the parent).  The alternative — adaptive red-green
    subdivision — needs T-junction handling and is overkill for the
    ~2× edge-length-reduction we typically need (1500 m → 750 m for
    a 1000 m target on 2000 m CFM input).
    """
    if target_edge_length_m <= 0.0:
        raise ValueError("target_edge_length_m must be > 0")
    V_out = np.asarray(V, dtype=np.float64).copy()
    T_out = np.asarray(T, dtype=np.int64).copy()
    if V_out.shape[0] == 0 or T_out.shape[0] == 0:
        return V_out, T_out

    # Determine how many lock-step refinement rounds we need based on
    # the current global max edge.  Refining EVERY triangle in lock-
    # step is what guarantees no T-junctions: each shared edge gets the
    # same midpoint on both sides because both adjacent triangles
    # request the same midpoint.  An alternative — refine only
    # triangles exceeding target — would T-junction every shared edge
    # whose neighbour didn't refine.
    a = V_out[T_out[:, 0]]
    b = V_out[T_out[:, 1]]
    c = V_out[T_out[:, 2]]
    max_edge = float(np.maximum.reduce([
        np.linalg.norm(b - a, axis=1),
        np.linalg.norm(c - b, axis=1),
        np.linalg.norm(a - c, axis=1),
    ]).max())
    if max_edge <= target_edge_length_m:
        return V_out, T_out
    # 1-to-4 splits halve each edge length.  Rounds needed:
    #   k = ceil(log2(max_edge / target))
    k_needed = int(math.ceil(math.log2(max_edge / target_edge_length_m)))
    if k_needed > max_rounds:
        _log.warning(
            "refine_uniform_midpoint: max_edge=%.1f m far exceeds target "
            "%.1f m; needs %d rounds but capping at %d (loop will exit "
            "with some triangles still over target)",
            max_edge, target_edge_length_m, k_needed, max_rounds,
        )
        k_needed = max_rounds

    for _round in range(k_needed):
        edge_to_mid: dict[tuple[int, int], int] = {}
        V_list: list[tuple[float, float, float]] = [
            (float(V_out[i, 0]), float(V_out[i, 1]), float(V_out[i, 2]))
            for i in range(V_out.shape[0])
        ]

        def get_midpoint(u: int, v: int) -> int:
            key = (min(u, v), max(u, v))
            if key not in edge_to_mid:
                p = 0.5 * (V_out[u] + V_out[v])
                edge_to_mid[key] = len(V_list)
                V_list.append((float(p[0]), float(p[1]), float(p[2])))
            return edge_to_mid[key]

        new_tris: list[tuple[int, int, int]] = []
        for ti in range(T_out.shape[0]):
            ai = int(T_out[ti, 0])
            bi = int(T_out[ti, 1])
            ci = int(T_out[ti, 2])
            m_ab = get_midpoint(ai, bi)
            m_bc = get_midpoint(bi, ci)
            m_ca = get_midpoint(ci, ai)
            new_tris.extend([
                (ai, m_ab, m_ca),
                (m_ab, bi, m_bc),
                (m_bc, ci, m_ca),
                (m_ab, m_bc, m_ca),
            ])
        V_out = np.asarray(V_list, dtype=np.float64)
        T_out = np.asarray(new_tris, dtype=np.int64)

    return V_out, T_out


def dedup_mesh(
    V: np.ndarray, T: np.ndarray, *, snap_m: float = 1e-3,
) -> tuple[np.ndarray, np.ndarray]:
    """Collapse vertices closer than ``snap_m`` to a single index and
    rebuild T accordingly.  Returns (V_unique, T_remapped).

    STL files store each triangle's three vertices in their own slots
    even when adjacent triangles share a vertex; ``_read_stl_arrays``
    preserves that 3-fold duplication.  Phase 2 (conformalize) needs
    proper connectivity (shared vertices) to run the manifold and
    edge-coincidence gates, so the caller dedups after reading.
    """
    V = np.asarray(V, dtype=np.float64)
    T = np.asarray(T, dtype=np.int64)
    if snap_m <= 0.0:
        raise ValueError(f"snap_m must be > 0; got {snap_m}")
    if V.shape[0] == 0:
        return V, T
    snapped = np.round(V / snap_m).astype(np.int64)
    _unique, inverse = np.unique(snapped, axis=0, return_inverse=True)
    n_unique = int(inverse.max()) + 1
    V_out = np.empty((n_unique, 3), dtype=np.float64)
    seen = np.zeros(n_unique, dtype=bool)
    for i in range(V.shape[0]):
        r = int(inverse[i])
        if not seen[r]:
            V_out[r] = V[i]
            seen[r] = True
    T_out = inverse[T]
    return V_out, T_out


def tri_tri_intersect_3d_interior_only(
    V_A: np.ndarray, T_A: np.ndarray,
    V_B: np.ndarray, T_B: np.ndarray,
    fault_a: str, fault_b: str,
    *,
    eps_orient: float = 1e-9,
    eps_min_seg_len_m: float = 5e-2,
    clearance_m: float = 0.0,
    vertex_match_tol_m: float = 2e-2,
) -> list[Segment]:
    """Phase 2 validation gate: same as :func:`tri_tri_intersect_3d` but
    filters out triangle pairs that **share one or more vertices**.

    After Phase 2 conformalization, post-split children of fault A and
    fault B share *exact edges* along the polyline; the bare
    `tri_tri_intersect_3d` would report those shared edges as full-edge-
    length crossings (correct, but not what we want to validate).  This
    helper drops shared-vertex pairs so the post-split set should report
    **zero** crossings on a correctly conformalised input.

    A "shared vertex" is detected by Euclidean distance below
    ``vertex_match_tol_m`` (default 2 cm = 2 × snap_m).
    """
    raw = tri_tri_intersect_3d(
        V_A, T_A, V_B, T_B, fault_a, fault_b,
        eps_orient=eps_orient,
        eps_min_seg_len_m=eps_min_seg_len_m,
        clearance_m=clearance_m,
    )
    out: list[Segment] = []
    V_A = np.asarray(V_A, dtype=np.float64)
    V_B = np.asarray(V_B, dtype=np.float64)
    T_A = np.asarray(T_A, dtype=np.int64)
    T_B = np.asarray(T_B, dtype=np.int64)
    tol_sq = vertex_match_tol_m ** 2
    for s in raw:
        verts_A = V_A[T_A[s.tri_a]]
        verts_B = V_B[T_B[s.tri_b]]
        # Any vertex of A within tol of any vertex of B → shared.
        diff = verts_A[:, None, :] - verts_B[None, :, :]
        d_sq = (diff * diff).sum(axis=-1)
        if (d_sq < tol_sq).any():
            continue
        out.append(s)
    return out


def scan_cross_fault_crossings(
    stl_dir: Path | str,
    included: list[str],
    clearance_m: float = 0.0,
    *,
    eps_orient: float = 1e-9,
    eps_min_seg_len_m: float = 5e-2,
    snap_m: float = 1e-2,
) -> dict[tuple[str, str], list[Polyline]]:
    """For every unordered pair (fault_a, fault_b) in ``included``,
    run :func:`tri_tri_intersect_3d` and :func:`chain_segments`.

    Returns a dict mapping (short_a, short_b) — *with short_a coming
    earlier in the input list than short_b* — to the list of polylines
    along which they intersect.  Pairs without any crossing are omitted
    from the dict; callers can use ``len(out) == 0`` as the "disjoint
    subsets" signal.

    Used by ``generate_safs_mesh.py`` as the gate that refuses
    non-conformalised multi-fault input.
    """
    stl_dir = Path(stl_dir)
    if not stl_dir.is_dir():
        raise FileNotFoundError(f"stl_dir does not exist: {stl_dir}")
    if not included:
        return {}

    cache: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for short in included:
        p = stl_dir / f"{short}.stl"
        if not p.exists():
            raise FileNotFoundError(
                f"missing STL for fault {short!r}: {p}"
            )
        cache[short] = _read_stl_arrays(p)

    out: dict[tuple[str, str], list[Polyline]] = {}
    for ai in range(len(included)):
        for bi in range(ai + 1, len(included)):
            sa = included[ai]
            sb = included[bi]
            V_A, T_A = cache[sa]
            V_B, T_B = cache[sb]
            segs = tri_tri_intersect_3d(
                V_A, T_A, V_B, T_B, sa, sb,
                eps_orient=eps_orient,
                eps_min_seg_len_m=eps_min_seg_len_m,
                clearance_m=clearance_m,
            )
            if not segs:
                continue
            polylines = chain_segments(segs, snap_m=snap_m)
            if polylines:
                out[(sa, sb)] = polylines
    return out
