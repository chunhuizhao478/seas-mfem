#!/usr/bin/env python3
"""build_collar.py -- extend a SAFS ALT mesh's footprint out to the ShakeOut box.

The parent mesh is FROZEN.  Its vertical absorbing side wall becomes the INNER
boundary of a new collar; we build an enlarged outer box around it and let
tetgen fill the region between with `-Y`, so the wall triangulation is
reproduced vertex-for-vertex and the collar welds on conformally.  Nothing
inside the old footprint moves -- fault, free surface, gate compliance and the
min-edge/dt floor are preserved by construction.

PLC (6 pieces, every edge must be used exactly twice):
    inner wall     FROZEN, straight from the parent
    top annulus    z = 0, between the parent's top rim and the new outer rim
    bottom annulus z = z_bot, likewise
    4 outer strips the new box's vertical sides, meshed in unrolled (s, z)

Every shared 1-D discretisation (outer top rim, outer bottom rim, the four
corner verticals) is generated ONCE and pinned into each 2-D problem with
transfinite-2 curves, so the pieces agree edge-for-edge by construction rather
than by a tolerance match afterwards.

Sizing honours the parent's own rule.  With --h-mode gate the collar keeps
f = Vs/dx >= gate (0.6667 = 0.5 Hz at p3), Vs pooled over each column's own
VERTICAL extent -- the refinement-stable target proven on the 1 Hz p5 build.
With --h-mode const it just matches the parent's far-field spacing (the small
mesh has no volume gate).

Usage:
    python build_collar.py --parent <p.puml.h5> --out <collar.npz> \
        --cvm <cvm.nc> --h-mode gate --h-max 5000
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from collar_lib import (SHAKEOUT_E, SHAKEOUT_N, TARGET_E, TARGET_N, VsGrid, as_ccw, boundary_loops,
                        check_closed, extract_wall, orient, signed_area,
                        suppress_stdout, tet_edge_lengths, tet_eta,
                        tet_signed_volume)

REG_TET_VOL = 1.0 / (6.0 * np.sqrt(2.0))


# ----------------------------------------------------------------- sizing ---
class SizeField:
    """h(x, y, z) = clip(safety * Vs_pooled / gate, h_min, h_max), or a constant.

    THE SAFETY FACTOR IS NOT SLACK -- it corrects a units mismatch.  The size
    field sets the length of TRIANGLE edges on the lid and the spacing of the
    interior point cloud, but the gate is judged on each TET's MAXIMUM edge,
    which is always longer than the base triangle's edge (the vertical and
    diagonal edges of a well-shaped tet run ~1.2-1.5x its base spacing).
    Sizing triangles at exactly Vs/gate therefore produces tets that miss it.
    Measured with safety = 1.0: 2,079,219 of 14,763,187 collar cells below gate,
    and 2,027,111 of them -- 97.5 % -- in the top 500 m, barycentres clustered
    at z = -30 to -100 m, i.e. precisely the surface layer whose triangles were
    sized at the limit.  This is the same reasoning as the campaign rule "judge
    mmg compliance at 1.25x the spec, never at 1.0x".
    """

    def __init__(self, mode, cvm, gate, h_min, h_max, const, safety=1.0):
        self.mode, self.gate, self.safety = mode, gate, safety
        self.h_min, self.h_max, self.const = h_min, h_max, const
        self.vs = VsGrid(cvm) if mode == "gate" and cvm else None

    def at_xy(self, xy, z_lo, z_hi):
        if self.mode == "const":
            return np.full(len(xy), self.const)
        v = self.vs.column_min(np.asarray(xy, float), z_lo, z_hi)
        return np.clip(self.safety * v / self.gate, self.h_min, self.h_max)


def walk_polyline(a, b, h_of):
    """Nodes along segment a->b spaced by the local size field (a kept, b not).

    Marching with the LOCAL size keeps the outer rim compliant with the same
    gate the volume is judged by, instead of a single global spacing.
    """
    L = float(np.linalg.norm(b - a))
    d = (b - a) / L
    out, s = [0.0], 0.0
    while True:
        p = a + d * s
        h = float(h_of(p[None, :])[0])
        s_next = s + max(h, 1.0)
        if s_next >= L - 0.5 * h:
            break
        out.append(s_next)
        s = s_next
    return a + d * np.array(out)[:, None]


def seed_interior(corners, rect_xy, z_top, depth, h_of, h_min, h_max):
    """Graded interior points for the shallow band, as extra tetgen input.

    WHY THIS EXISTS.  tetgen's only size controls are the boundary point
    density and ONE global volume cap; it has no 3-D size field.  Given a fine
    top annulus and a coarse cap it grades away from the surface far faster
    than the 0.5 Hz gate allows, and the measured result was 1.17M of 8.80M
    collar cells below gate.  Trying to repair that afterwards by longest-edge
    bisection does not work -- bisection fixes SIZE, not GRADING: it tripled the
    collar to 28.9M tets, left MORE gate failures (1.31M) than it started with,
    and drove the minimum edge from 65 m to 8.0 m.

    Interior points fix it at the source.  tetgen inserts input points that lie
    in no facet, so a graded cloud sets the local element size directly.  Only
    the shallow band needs it: below ~2 km the CVM's Vs already permits 2.5-4 km
    cells, which the global cap delivers on its own.

    Points are laid on a per-level lattice: level h keeps lattice sites whose
    required size falls in [h, 2h), so each level occupies its own depth band
    and the cloud is graded rather than uniform.  Sites are kept clear of the
    PLC (0.4 h from the lid, the wall's footprint and the outer box) so tetgen
    is never asked to insert a point on top of a frozen facet.
    """
    lo = np.array([corners[:, 0].min(), corners[:, 1].min()])
    hi = np.array([corners[:, 0].max(), corners[:, 1].max()])
    levels = []
    h = h_max
    while h > h_min * 1.01:
        levels.append(h)
        h *= 0.5
    levels.append(h_min)

    def outside_rect(P):
        """True where P is OUTSIDE the parent footprint (i.e. in the collar)."""
        c = rect_xy.mean(0)
        order = np.argsort(np.arctan2(rect_xy[:, 1] - c[1], rect_xy[:, 0] - c[0]))
        R = rect_xy[order]
        inside = np.ones(len(P), bool)
        for i in range(4):
            a, b = R[i], R[(i + 1) % 4]
            e = b - a
            inside &= (e[0] * (P[:, 1] - a[1]) - e[1] * (P[:, 0] - a[0])) >= 0.0
        return ~inside

    def rect_clear(P, m):
        """Distance from the parent footprint boundary, for a margin test."""
        c = rect_xy.mean(0)
        order = np.argsort(np.arctan2(rect_xy[:, 1] - c[1], rect_xy[:, 0] - c[0]))
        R = rect_xy[order]
        d = np.full(len(P), np.inf)
        for i in range(4):
            a, b = R[i], R[(i + 1) % 4]
            e = b - a
            e = e / np.linalg.norm(e)
            n = np.array([-e[1], e[0]])
            d = np.minimum(d, np.abs((P - a) @ n))
        return d > m

    out = []
    for h in levels:
        zs = np.arange(z_top - 0.5 * h, z_top - depth - h, -h)
        if not len(zs):
            continue
        gx = np.arange(lo[0] + 0.5 * h, hi[0], h)
        gy = np.arange(lo[1] + 0.5 * h, hi[1], h)
        if len(gx) * len(gy) > 40_000_000:
            continue
        X, Y = np.meshgrid(gx, gy, indexing="xy")
        XY = np.column_stack([X.ravel(), Y.ravel()])
        keep_xy = (outside_rect(XY) & rect_clear(XY, 0.4 * h)
                   & (XY[:, 0] > lo[0] + 0.4 * h) & (XY[:, 0] < hi[0] - 0.4 * h)
                   & (XY[:, 1] > lo[1] + 0.4 * h) & (XY[:, 1] < hi[1] - 0.4 * h))
        XY = XY[keep_xy]
        if not len(XY):
            continue
        for z in zs:
            P = np.column_stack([XY, np.full(len(XY), z)])
            want = h_of(P)
            sel = (want >= h) & (want < 2.0 * h)
            if sel.any():
                out.append(P[sel])
    if not out:
        return np.zeros((0, 3))
    # Enforce a minimum spacing BETWEEN levels.  Each level is a clean lattice,
    # but two adjacent levels straddling a size transition sit on different
    # lattices and can land within metres of one another; tetgen then has to
    # honour both and emits a sliver.  Measured without this filter: collar
    # min edge 15 m and eta_min 0.0211, against 65 m / 0.0674 unseeded.
    # Coarse levels are accepted first, so the surviving cloud stays graded.
    keep = [out[0]]
    tree = cKDTree(out[0])
    for blk in out[1:]:
        h_blk = None
        d, _ = tree.query(blk)
        # spacing floor tied to this level's own lattice pitch
        pitch = np.median(np.linalg.norm(blk[1:] - blk[:-1], axis=1)) if len(blk) > 1 else 0.0
        ok = d > max(0.45 * pitch, 1.0)
        if ok.any():
            keep.append(blk[ok])
            tree = cKDTree(np.vstack(keep))
    return np.vstack(keep)


def z_ladder(z_top, z_bot, h_start, h_cap, ratio=1.3):
    """Geometrically graded z levels, landing EXACTLY on z_top and z_bot.

    The two end values are snapped, not merely rescaled to.  `z_top -
    cumsum(steps)[-1]` lands within ~1e-11 m of z_bot, and the strips stamp
    their bottom nodes with zs[-1] while the bottom annulus stamps z_bot -- a
    difference that small still misses the exact-key weld and leaves the PLC
    open along the entire outer bottom rim (measured: 976 open edges, all at
    z = -22,100.939).
    """
    thick = z_top - z_bot
    steps, h = [], min(h_start, h_cap)
    while sum(steps) < thick:
        steps.append(h)
        h = min(h * ratio, h_cap)
    steps = np.array(steps) * (thick / sum(steps))
    zs = np.concatenate([[z_top], z_top - np.cumsum(steps)])
    zs[0], zs[-1] = z_top, z_bot
    return zs


# ------------------------------------------------------------- gmsh pieces ---
def _pin(geo, gmsh, pts_tags):
    """Chain of transfinite-2 lines through pts_tags (closed if first==last)."""
    lines = []
    for i in range(len(pts_tags) - 1):
        ln = geo.addLine(pts_tags[i], pts_tags[i + 1])
        lines.append(ln)
    return lines


def mesh_planar(outer_xy, inner_xy, z, h_of, log, tag):
    """Triangulate the annulus between outer_xy (CCW) and inner_xy (hole) at z.

    Every boundary segment is pinned transfinite-2, so the returned mesh
    contains the given boundary nodes and NO others on the boundary -- that is
    what lets this piece agree exactly with its neighbours.
    """
    import gmsh
    gmsh.initialize()
    try:
        with suppress_stdout(log):
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add(tag)
            geo = gmsh.model.geo
            hmean = float(np.mean(h_of(outer_xy)))
            oz = [geo.addPoint(float(x), float(y), float(z), hmean) for x, y in outer_xy]
            iz = [geo.addPoint(float(x), float(y), float(z), hmean) for x, y in inner_xy]
            ol = _pin(geo, gmsh, oz + [oz[0]])
            il = _pin(geo, gmsh, iz + [iz[0]])
            surf = geo.addPlaneSurface([geo.addCurveLoop(ol), geo.addCurveLoop(il)])
            geo.synchronize()
            for ln in ol + il:
                gmsh.model.mesh.setTransfiniteCurve(ln, 2)
            _bg_field(gmsh, h_of, outer_xy, z)
            gmsh.model.mesh.generate(2)
            P, T, nfix = _harvest(gmsh, surf, oz + iz)
        # Re-stamp the pinned nodes with the EXACT input coordinates.  gmsh
        # round-trips them through its own storage, and a sub-mm drift is
        # enough to break the exact-key weld between PLC pieces.
        exact = np.column_stack([np.vstack([outer_xy, inner_xy]),
                                 np.full(len(outer_xy) + len(inner_xy), z)])
        P[:nfix] = exact
        return P, T, nfix
    finally:
        gmsh.finalize()


def mesh_strip(corner_a, corner_b, top_nodes, bot_nodes, zs, h_of, log, tag):
    """One vertical side of the outer box, meshed in unrolled (s, z).

    top_nodes/bot_nodes are the already-fixed rim discretisations of this side
    (including both corners); zs is the shared corner z-ladder.  Meshing in
    (s, z) lets the strip coarsen with depth, which a structured quad ring
    cannot -- a ring forced to carry the fine surface spacing down to 40 km
    would pin ~1000 m x 5000 m sliver facets along the whole outer boundary.
    """
    import gmsh
    L = float(np.linalg.norm(corner_b - corner_a))
    d = (corner_b - corner_a) / L
    s_top = (top_nodes - corner_a) @ d
    s_bot = (bot_nodes - corner_a) @ d

    def to_xyz(s, zz):
        p = corner_a[None, :] + d[None, :] * np.asarray(s)[:, None]
        return np.column_stack([p, zz])

    def h_sz(SZ):
        return h_of(to_xyz(SZ[:, 0], SZ[:, 1]))

    gmsh.initialize()
    try:
        with suppress_stdout(log):
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add(tag)
            geo = gmsh.model.geo
            hm = float(np.mean(h_sz(np.column_stack([s_top, np.zeros(len(s_top))]))))
            pt_top = [geo.addPoint(float(s), float(zs[0]), 0.0, hm) for s in s_top]
            pt_bot = [geo.addPoint(float(s), float(zs[-1]), 0.0, hm) for s in s_bot]
            pt_l = [geo.addPoint(0.0, float(zz), 0.0, hm) for zz in zs[1:-1]]
            pt_r = [geo.addPoint(float(L), float(zz), 0.0, hm) for zz in zs[1:-1]]
            left = [pt_top[0]] + pt_l + [pt_bot[0]]
            right = [pt_top[-1]] + pt_r + [pt_bot[-1]]
            ltop = _pin(geo, gmsh, pt_top)
            lright = _pin(geo, gmsh, right)
            lbot = _pin(geo, gmsh, pt_bot[::-1])
            lleft = _pin(geo, gmsh, left[::-1])
            loop = geo.addCurveLoop(ltop + lright + lbot + lleft)
            surf = geo.addPlaneSurface([loop])
            geo.synchronize()
            for ln in ltop + lright + lbot + lleft:
                gmsh.model.mesh.setTransfiniteCurve(ln, 2)
            _bg_field_sz(gmsh, h_sz, s_top, zs)
            gmsh.model.mesh.generate(2)
            P2, T, nfix = _harvest(gmsh, surf, pt_top + pt_bot + pt_l + pt_r)
        P = to_xyz(P2[:, 0], P2[:, 1])
        # Same re-stamp as mesh_planar, and here it matters more: the strip is
        # meshed in unrolled (s, z) and mapped BACK, so every pinned node would
        # otherwise carry the round-trip error of a dot product on ~4e6 m UTM
        # coordinates.  These four groups are exactly the nodes shared with the
        # annuli and with the neighbouring strips.
        nt_, nb_ = len(top_nodes), len(bot_nodes)
        exact = np.vstack([
            np.column_stack([top_nodes, np.full(nt_, zs[0])]),
            np.column_stack([bot_nodes, np.full(nb_, zs[-1])]),
            np.column_stack([np.repeat(corner_a[None, :], len(zs) - 2, 0), zs[1:-1]]),
            np.column_stack([np.repeat(corner_b[None, :], len(zs) - 2, 0), zs[1:-1]]),
        ])
        if len(exact) != nfix:
            raise RuntimeError(f"strip pin count {nfix} != expected {len(exact)}")
        P[:nfix] = exact
        return P, T
    finally:
        gmsh.finalize()


def _lat(extent, step, lo=32, hi=760):
    """Lattice count resolving `step` over `extent` (the CVM's own spacing).

    A fixed coarse lattice is wrong here: the size field is driven by the CVM,
    whose lateral grid is 1500 m and vertical 250 m.  Sampling the 706 x 464 km
    top annulus on a 160^2 lattice would be 4.4 km spacing -- it would smooth
    away exactly the narrow low-Vs basins that set the gate.
    """
    return int(np.clip(np.ceil(extent / step) + 1, lo, hi))


def _bg_field(gmsh, h_of, outer_xy, z):
    """Background size view on a lattice resolving the CVM's 1500 m spacing."""
    lo, hi = outer_xy.min(0), outer_xy.max(0)
    nx = _lat(hi[0] - lo[0], 1500.0)
    ny = _lat(hi[1] - lo[1], 1500.0)
    gx = np.linspace(lo[0], hi[0], nx)
    gy = np.linspace(lo[1], hi[1], ny)
    X, Y = np.meshgrid(gx, gy, indexing="xy")
    H = h_of(np.column_stack([X.ravel(), Y.ravel()])).reshape(ny, nx)
    _push_view(gmsh, gx, gy, H, z)


def _bg_field_sz(gmsh, h_sz, s_top, zs):
    ns = _lat(float(s_top[-1]), 1500.0)
    nz = _lat(float(zs[0] - zs[-1]), 250.0)
    gs = np.linspace(0.0, float(s_top[-1]), ns)
    gz = np.linspace(float(zs[-1]), float(zs[0]), nz)
    S, Z = np.meshgrid(gs, gz, indexing="xy")
    H = h_sz(np.column_stack([S.ravel(), Z.ravel()])).reshape(nz, ns)
    _push_view(gmsh, gs, gz, H, 0.0)


def _push_view(gmsh, ga, gb, H, third):
    """Register H(ga,gb) as a scalar-triangle view and use it as the size field.

    Vectorised: the lattice now resolves the CVM (up to ~480 x 320), so the
    naive Python double loop would build ~3 x 10^5 lists per piece.
    "ST" layout is x1 x2 x3 y1 y2 y3 z1 z2 z3 v1 v2 v3.
    """
    A, Bb = np.meshgrid(ga, gb, indexing="xy")

    def corner(dj, di):
        return (A[dj:dj + len(gb) - 1, di:di + len(ga) - 1].ravel(),
                Bb[dj:dj + len(gb) - 1, di:di + len(ga) - 1].ravel(),
                H[dj:dj + len(gb) - 1, di:di + len(ga) - 1].ravel())

    c00, c01, c11, c10 = corner(0, 0), corner(0, 1), corner(1, 1), corner(1, 0)
    n = len(c00[0])
    z3 = np.full(n, third)
    t1 = np.column_stack([c00[0], c01[0], c11[0], c00[1], c01[1], c11[1],
                          z3, z3, z3, c00[2], c01[2], c11[2]])
    t2 = np.column_stack([c00[0], c11[0], c10[0], c00[1], c11[1], c10[1],
                          z3, z3, z3, c00[2], c11[2], c10[2]])
    tri = np.vstack([t1, t2])
    v = gmsh.view.add("h")
    gmsh.view.addListData(v, "ST", len(tri), tri.ravel())
    f = gmsh.model.mesh.field.add("PostView")
    gmsh.model.mesh.field.setNumber(f, "ViewTag", v)
    gmsh.model.mesh.field.setAsBackgroundMesh(f)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)


def _harvest(gmsh, surf, fixed_tags):
    """Nodes+triangles of a 2-D surface, with the pinned points listed FIRST."""
    ntags, ncoord, _ = gmsh.model.mesh.getNodes()
    ncoord = np.asarray(ncoord).reshape(-1, 3)
    order = {int(t): i for i, t in enumerate(ntags)}
    fixed = [order[int(gmsh.model.mesh.getNodes(0, t)[0][0])] for t in fixed_tags]
    etypes, _, enodes = gmsh.model.mesh.getElements(2, surf)
    tri = None
    for et, en in zip(etypes, enodes):
        if int(et) == 2:
            tri = np.array([order[int(v)] for v in en], np.int64).reshape(-1, 3)
    if tri is None:
        raise RuntimeError("gmsh produced no triangles")
    perm = np.empty(len(ncoord), np.int64)
    perm[fixed] = np.arange(len(fixed))
    rest = np.setdiff1d(np.arange(len(ncoord)), fixed)
    perm[rest] = np.arange(len(rest)) + len(fixed)
    P = np.empty_like(ncoord)
    P[perm] = ncoord
    return P, perm[tri], len(fixed)


# ------------------------------------------------------------------- main ---
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parent", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cvm", default=None)
    ap.add_argument("--h-mode", choices=["gate", "const"], default="gate")
    ap.add_argument("--h-const", type=float, default=3000.0)
    ap.add_argument("--gate", type=float, default=0.6667)
    ap.add_argument("--h-min", type=float, default=250.0)
    ap.add_argument("--h-max", type=float, default=5000.0)
    ap.add_argument("--h-safety", type=float, default=0.75,
                    help="triangle-edge to tet-max-edge correction; see SizeField")
    ap.add_argument("--h-safety-lid", type=float, default=None,
                    help="tighter safety for the TOP annulus only.  The lid is "
                         "where the gate binds (93-97 %% of residual failures sit "
                         "in the top 500 m) because its safety must absorb TWO "
                         "overshoots: gmsh returns edges ~1.0-1.3x its size field, "
                         "and a tet max edge exceeds its base triangle edge.  "
                         "Tightening only the lid avoids paying for both over the "
                         "whole 40 km column.  Defaults to --h-safety.")
    ap.add_argument("--ratio", type=float, default=1.3)
    ap.add_argument("--margin", type=float, default=10000.0,
                    help="minimum collar width where the parent already reaches "
                         "or exceeds the ShakeOut box (see the box comment)")
    ap.add_argument("--minratio", type=float, default=1.414)
    ap.add_argument("--mindihedral", type=float, default=18.0)
    ap.add_argument("--seed-top", type=float, default=0.0,
                    help="seed graded interior points in the top N metres; see "
                         "seed_interior() for why tetgen needs them")
    ap.add_argument("--plc-only", action="store_true")
    ap.add_argument("--log", default=None)
    a = ap.parse_args()

    t0 = time.time()
    log = a.log or str(Path(a.out).with_suffix("")) + "_build.log"
    open(log, "w").close()
    import h5py
    with h5py.File(a.parent, "r") as f:
        G = f["geometry"][:]
        C = f["connect"][:].astype(np.int64)
        B = f["boundary"][:].astype(np.int32)
    print(f"[parent] {len(C):,} tets, {len(G):,} verts", flush=True)

    PW, TW, wall_global, wall_tet, wall_slot = extract_wall(G, C, B)
    loops = boundary_loops(TW)
    if len(loops) != 2:
        raise RuntimeError(f"expected 2 wall rims, got {len(loops)}")
    loops.sort(key=lambda L: -PW[L][:, 2].mean())
    top_rim, bot_rim = loops
    # The annuli are built as PLANES, so each rim must be exactly level or the
    # weld to the wall would miss.  The ALT lid is exactly flat (measured: 0 of
    # 187,161 lid vertex-uses have |z| > 1e-6); assert rather than assume.
    for nm, L in (("top", top_rim), ("bottom", bot_rim)):
        zz = PW[L][:, 2]
        if float(np.ptp(zz)) > 1e-6:
            raise RuntimeError(f"{nm} rim is not level: z spread {np.ptp(zz):.3e} m "
                               f"({zz.min():.6f} .. {zz.max():.6f})")
    z_top = float(PW[top_rim][:, 2][0])
    z_bot = float(PW[bot_rim][:, 2][0])
    print(f"[wall]   {len(TW):,} tris, {len(PW):,} verts, z {z_bot:,.1f} .. {z_top:,.1f}")
    print(f"         top rim {len(top_rim):,} verts, bottom rim {len(bot_rim):,} verts")

    # inner rims oriented CW (they are holes in the annuli, outer is CCW)
    top_rim, top_xy = as_ccw(top_rim, PW[top_rim][:, :2])
    bot_rim, bot_xy = as_ccw(bot_rim, PW[bot_rim][:, :2])

    # ---- the enlarged box ---------------------------------------------------
    # box = TARGET_E/N (see collar_lib) UNION (parent bbox + margin).
    #
    # TARGET already contains the ShakeOut grid box with margin AND the
    # PREFERRED domain box; the union with parent+margin is the safety net that
    # guarantees requirement (2)+(3) for whichever parent is passed -- a
    # frozen-parent extension can only ADD, and a wall drawn exactly at a
    # parent corner pinches the collar to zero width, where gmsh reports
    # "2 intersections in the 1D mesh" and emits no elements at all.
    ex = (min(TARGET_E[0], G[:, 0].min() - a.margin),
          max(TARGET_E[1], G[:, 0].max() + a.margin))
    ny_ = (min(TARGET_N[0], G[:, 1].min() - a.margin),
           max(TARGET_N[1], G[:, 1].max() + a.margin))
    corners = np.array([[ex[0], ny_[0]], [ex[1], ny_[0]], [ex[1], ny_[1]], [ex[0], ny_[1]]])
    print(f"[box]    E {ex[0]:,.1f} .. {ex[1]:,.1f}  ({(ex[1]-ex[0])/1e3:.3f} km)")
    print(f"         N {ny_[0]:,.1f} .. {ny_[1]:,.1f}  ({(ny_[1]-ny_[0])/1e3:.3f} km)")
    a_new = (ex[1] - ex[0]) * (ny_[1] - ny_[0]) / 1e6
    a_old = abs(signed_area(top_xy)) / 1e6
    print(f"         footprint {a_old:,.0f} -> {a_new:,.0f} km2  (collar {a_new-a_old:,.0f} km2)")

    sf = SizeField(a.h_mode, a.cvm, a.gate, a.h_min, a.h_max, a.h_const, a.h_safety)
    sf_lid = SizeField(a.h_mode, None, a.gate, a.h_min, a.h_max, a.h_const,
                       a.h_safety_lid if a.h_safety_lid else a.h_safety)
    sf_lid.vs = sf.vs
    h_top = lambda xy: sf_lid.at_xy(np.asarray(xy)[:, :2], -300.0, 0.0)
    h_bot = lambda xy: sf.at_xy(np.asarray(xy)[:, :2], z_bot, z_bot + 500.0)
    h_vol = lambda xyz: np.clip(
        sf.at_xy(np.asarray(xyz)[:, :2], -300.0, 0.0) if sf.mode == "const" else
        np.clip(sf.safety * sf.vs.at(np.asarray(xyz)) / sf.gate, sf.h_min, sf.h_max),
        sf.h_min, sf.h_max)

    # ---- shared 1-D discretisations ----------------------------------------
    zs = z_ladder(z_top, z_bot, float(np.median(np.linalg.norm(
        PW[top_rim][:, :2] - np.roll(PW[top_rim][:, :2], -1, 0), axis=1))),
        a.h_max, a.ratio)
    print(f"[ladder] {len(zs)} levels, first step {zs[0]-zs[1]:,.0f} m, "
          f"last {zs[-2]-zs[-1]:,.0f} m")

    otop_sides, obot_sides = [], []
    for i in range(4):
        A, Bc = corners[i], corners[(i + 1) % 4]
        otop_sides.append(walk_polyline(A, Bc, h_top))
        obot_sides.append(walk_polyline(A, Bc, h_bot))
    otop = np.vstack(otop_sides)
    obot = np.vstack(obot_sides)
    print(f"[rims]   outer top {len(otop):,} nodes, outer bottom {len(obot):,} nodes")

    # ---- the six PLC pieces -------------------------------------------------
    V, F = [], []
    nv = 0

    piece = []

    def push(P, T, name="?"):
        nonlocal nv
        V.append(P)
        F.append(T + nv)
        piece.append(np.full(len(P), name, dtype=object))
        nv += len(P)

    # 1. inner wall, frozen -- listed FIRST so indices 0..len(PW)-1 are the
    #    parent's wall vertices, which is what the weld map keys on.
    push(np.asarray(PW, float), np.asarray(TW, np.int64), "wall")
    wall_n = len(PW)

    # 2. top annulus (z = z_top)
    Pt, Tt, _ = mesh_planar(otop, PW[top_rim][:, :2], z_top, h_top, log, "top")
    Tt = orient(Tt, Pt, want_up=True)
    push(Pt, Tt, "top")
    print(f"[top]    {len(Tt):,} tris, {len(Pt):,} verts", flush=True)

    # 3. bottom annulus
    Pb, Tb, _ = mesh_planar(obot, PW[bot_rim][:, :2], z_bot, h_bot, log, "bot")
    Tb = orient(Tb, Pb, want_up=False)
    push(Pb, Tb, "bottom")
    print(f"[bottom] {len(Tb):,} tris, {len(Pb):,} verts", flush=True)

    # 4-7. the four outer strips
    for i in range(4):
        A, Bc = corners[i], corners[(i + 1) % 4]
        tn = np.vstack([otop_sides[i], corners[(i + 1) % 4][None, :]])
        bn = np.vstack([obot_sides[i], corners[(i + 1) % 4][None, :]])
        Ps, Ts = mesh_strip(A, Bc, tn, bn, zs, h_vol, log, f"side{i}")
        push(Ps, Ts, f"side{i}")
        print(f"[side{i}]  {len(Ts):,} tris, {len(Ps):,} verts", flush=True)

    piece_all = np.concatenate(piece)
    Vall = np.vstack(V)
    Fall = np.vstack(F).astype(np.int64)

    # Weld duplicate nodes across pieces.  The rims are shared BY CONSTRUCTION
    # (identical float64 coordinates handed to each gmsh problem), so a 0.1 mm
    # rounding key is exact here and O(n log n) -- a radius join on ~10^6 PLC
    # points is not affordable.  Genuinely distinct nodes are >= 100 m apart.
    key = np.round(Vall * 1e4).astype(np.int64)
    _, uniq, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    Vu = Vall[uniq]
    Fu = inv[Fall]
    keep = (Fu[:, 0] != Fu[:, 1]) & (Fu[:, 1] != Fu[:, 2]) & (Fu[:, 0] != Fu[:, 2])
    Fu = Fu[keep]
    wall_map = inv[np.arange(wall_n)]
    print(f"[plc]    {len(Vu):,} verts, {len(Fu):,} facets "
          f"(welded {len(Vall)-len(Vu):,} duplicates, dropped {int((~keep).sum())} degenerate)")
    check_closed(Fu, Vu, piece_all[uniq])
    print("[plc]    CLOSED: every edge used exactly twice")

    if a.plc_only:
        np.savez_compressed(a.out, plc_points=Vu, plc_facets=Fu, wall_map=wall_map,
                            wall_global=wall_global, wall_tet=wall_tet, wall_slot=wall_slot,
                            z_top=z_top, z_bot=z_bot, box=corners)
        print(f"[write] PLC only -> {a.out}  ({time.time()-t0:.0f} s)")
        return

    # ---- graded interior seed for the shallow band -------------------------
    n_plc = len(Vu)
    if a.seed_top > 0:
        seed = seed_interior(corners, PW[top_rim][:, :2], z_top, a.seed_top,
                             h_vol, a.h_min, a.h_max)
        print(f"[seed]   {len(seed):,} interior points in the top "
              f"{a.seed_top/1e3:.1f} km", flush=True)
        Vu = np.vstack([Vu, seed])

    # ---- fill ---------------------------------------------------------------
    import tetgen
    shift = Vu.mean(axis=0)
    vol_cap = REG_TET_VOL * a.h_max ** 3
    sw = f"pq{a.minratio}/{a.mindihedral}YAa{vol_cap:.6e}"
    print(f"[tetgen] switches {sw}", flush=True)
    tg = tetgen.TetGen(np.ascontiguousarray(Vu - shift), np.ascontiguousarray(Fu))
    with suppress_stdout(log):
        tg.tetrahedralize(switches=sw)
    SP = np.asarray(tg.node) + shift
    ST = np.asarray(tg.elem).astype(np.int64)
    sv = tet_signed_volume(SP, ST)
    if np.any(sv < 0):
        ST[sv < 0] = ST[sv < 0][:, [0, 1, 3, 2]]
    print(f"[collar] {len(ST):,} tets, {len(SP):,} verts  ({time.time()-t0:.0f} s)")

    # the wall vertices must have survived -Y untouched
    d, near = cKDTree(SP).query(Vu[:n_plc][wall_map])
    if d.max() > 1e-6:
        raise RuntimeError(f"tetgen moved wall vertices (max {d.max():.3e} m)")
    print(f"[collar] wall vertices preserved to {d.max():.2e} m")

    eta = tet_eta(SP, ST)
    ed = tet_edge_lengths(SP, ST)
    print(f"  quality: eta min {eta.min():.4f} med {np.median(eta):.3f}, "
          f"eta<0.1 {int((eta<0.1).sum()):,}, eta<0.05 {int((eta<0.05).sum()):,}")
    print(f"  edges: min {ed.min():,.0f} max {ed.max():,.0f} m")

    np.savez_compressed(a.out, points=SP, tets=ST,
                        wall_collar_idx=near, wall_parent_idx=wall_global,
                        wall_tet_idx=wall_tet, wall_slot=wall_slot,
                        z_top=z_top, z_bot=z_bot, box=corners)
    print(f"[write] {a.out}  ({time.time()-t0:.0f} s total)")


if __name__ == "__main__":
    main()
