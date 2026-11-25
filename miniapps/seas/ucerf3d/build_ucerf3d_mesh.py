#!/usr/bin/env python3
"""UCERF 3-D fault-opening model: geometry and gmsh mesh (plan rev5, phase P0).

Two jobs, both with gmsh's built-in kernel (Vista's gmsh has no OpenCASCADE):

  ucerf  the San Andreas / San Jacinto model.  The five R 70 km chains of the
         2.5D mesh (ucerf2d/mesh/ucerf2d_saf_sjf_500m_tilted_D415km_R70.msh,
         the 10 km depth slice of UCERF3 FM3.1) are put back at their depth:
           - saf1, saf3, sjf1, sjf2 are vertical: the chain extruded straight
             down;
           - saf2 (San Gorgonio Pass - Garnet Hill) dips 58 deg toward azimuth
             19.8 deg: at depth d the chain is moved along the dip direction
             by (d - 10 km) / tan 58 deg.
         Each chain gives a friction surface from 0 to --fault_depth and an
         option A creep surface from there to --creep_depth (E4: 25-50 km).
         The ends of saf2 are cut back until it is at least --gap from saf1
         and saf3 at every depth (D6: three separate surfaces, test M3).
  test   the small test meshes of plan section 10.1: T-plane, T-dip, T-bend,
         T-pml.

Boundary attributes (the fd_cycle convention, cycle_config.hpp):
   3  friction fault faces (interior)
   6  option A creep faces below the friction surface (interior)
   1  top, z = 0 (natural: free surface)
   2  bottom (natural or Dirichlet, per deck)
   5  side walls (Dirichlet: plate loading)
 Volume attribute 1.  The PML band is not tagged: the cycle code builds it
 from the wall faces and [cycle.pml] thickness_km.

Mesh size (plan 9.4): h = h_f within --d0 of a friction surface, then
h_f + --grade * (d - d0), capped at --hmax.  The creep surfaces and the volume
take their size from that field.

Units: metres; x east, y north, z up (z = 0 is the free surface).

Output: gmsh 2.2 ASCII (what MFEM reads) and a JSON record of the geometry,
the cut-back, the gaps and the counts.  gmsh runs on one thread, so a rebuild
with the same gmsh version gives the same file.
"""
import argparse
import collections
import hashlib
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SEAS = os.path.dirname(HERE)
SRC_2D = os.path.join(SEAS, "ucerf2d", "mesh", "ucerf2d_saf_sjf_500m_tilted_D415km_R70.msh")
SRC_2D_JSON = os.path.join(SEAS, "ucerf2d", "mesh", "ucerf2d_saf_sjf_500m_tilted_D415km_R70.json")
SRC_2D_SHA = "d3df54b3379f5042d79897a593a91af419647c6d5526c000f57cb67b2e4a7dbd"

# UCERF3 FM3.1, San Andreas (San Gorgonio Pass-Garnet Hill): aveDip 58,
# dipDirection 19.8 (plan rev5 section 4).
SAF2_DIP_DEG = 58.0
SAF2_DIPDIR_DEG = 19.8
SLICE_DEPTH = 10.0e3          # the 2.5D chains are the 10 km slice

# chain names by length (km) of the R 70 km chains, plan rev5 section 4
CHAIN_LENGTHS = {"saf1": 420.7, "saf2": 55.9, "saf3": 69.5, "sjf1": 173.1, "sjf2": 103.5}

ATTR_TOP, ATTR_BOTTOM, ATTR_FAULT, ATTR_WALL, ATTR_CREEP = 1, 2, 3, 5, 6


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ----------------------------------------------------------------------------
# chains
# ----------------------------------------------------------------------------
def read_chains_2d(path):
    """Fault chains (attribute 3 edges) of a gmsh 2.2 2-D mesh, in metres."""
    nodes, edges = {}, []
    with open(path) as f:
        line = f.readline()
        while line:
            if line.startswith("$Nodes"):
                for _ in range(int(f.readline())):
                    p = f.readline().split()
                    nodes[int(p[0])] = (float(p[1]), float(p[2]))
            elif line.startswith("$Elements"):
                for _ in range(int(f.readline())):
                    p = list(map(int, f.readline().split()))
                    if p[1] == 1 and p[3] == ATTR_FAULT:
                        edges.append((p[-2], p[-1]))
            line = f.readline()
    adj = collections.defaultdict(list)
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    seen, chains = set(), []
    for s in sorted(adj):
        if len(adj[s]) == 1 and s not in seen:
            c, cur = [s], s
            seen.add(s)
            while True:
                nx = [q for q in adj[cur] if q not in seen]
                if not nx:
                    break
                cur = nx[0]
                seen.add(cur)
                c.append(cur)
            chains.append(np.array([nodes[i] for i in c]))
    return chains


def arclen(c):
    return np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(c, axis=0).T))])


def name_chains(chains):
    named = {}
    for c in chains:
        L = arclen(c)[-1] / 1e3
        n = min(CHAIN_LENGTHS, key=lambda k: abs(CHAIN_LENGTHS[k] - L))
        if abs(CHAIN_LENGTHS[n] - L) > 1.0 or n in named:
            sys.exit(f"chain of {L:.1f} km does not match {CHAIN_LENGTHS}")
        named[n] = c
    if len(named) != 5:
        sys.exit(f"expected 5 chains, found {sorted(named)}")
    return named


def trim(c, s0, s1):
    """The part of polyline c between arc lengths s0 and L - s1."""
    s = arclen(c)
    a, b = s0, s[-1] - s1
    keep = (s > a) & (s < b)
    pa = np.array([np.interp(a, s, c[:, 0]), np.interp(a, s, c[:, 1])])
    pb = np.array([np.interp(b, s, c[:, 0]), np.interp(b, s, c[:, 1])])
    out = [pa] + list(c[keep]) + [pb]
    # drop a kept vertex closer than 1/4 of the median spacing to a new end
    out = np.array(out)
    ds = np.median(np.diff(s))
    if len(out) > 3 and np.hypot(*(out[1] - out[0])) < 0.25 * ds:
        out = np.delete(out, 1, axis=0)
    if len(out) > 3 and np.hypot(*(out[-2] - out[-1])) < 0.25 * ds:
        out = np.delete(out, -2, axis=0)
    return out


class Surface:
    """A ruled fault surface: chain c (the 10 km slice, or the trace for a
    vertical fault) moved by shift(d) at depth d >= 0 (metres)."""

    def __init__(self, name, chain, dip_deg=90.0, dipdir_deg=0.0, ref_depth=SLICE_DEPTH):
        self.name, self.chain = name, np.asarray(chain, float)
        self.dip, self.dipdir, self.ref = dip_deg, dipdir_deg, ref_depth
        if dip_deg < 90.0:
            t = 1.0 / math.tan(math.radians(dip_deg))
            u = np.array([math.sin(math.radians(dipdir_deg)), math.cos(math.radians(dipdir_deg))])
            self.rate = t * u          # horizontal shift per metre of depth
        else:
            self.rate = np.zeros(2)

    def at(self, d):
        """Horizontal polyline at depth d."""
        return self.chain + self.rate * (d - self.ref)

    def shift(self, d0, d1):
        """(dx, dy, dz) from depth d0 to depth d1 (z up)."""
        r = self.rate * (d1 - d0)
        return float(r[0]), float(r[1]), -(d1 - d0)


def seg_dist(P, Q):
    """Distance from each point of P (n,2) to the polyline Q (m,2)."""
    a, b = Q[:-1], Q[1:]
    ab = b - a
    L2 = np.maximum((ab * ab).sum(1), 1e-30)
    d = np.full(len(P), np.inf)
    for i in range(0, len(P), 256):
        p = P[i:i + 256, None, :]
        t = np.clip(((p - a) * ab).sum(2) / L2, 0.0, 1.0)
        r = p - a - t[..., None] * ab
        d[i:i + 256] = np.sqrt((r * r).sum(2)).min(1)
    return d


def densify(c, ds=50.0):
    s = arclen(c)
    t = np.linspace(0.0, s[-1], max(2, int(math.ceil(s[-1] / ds)) + 1))
    return np.column_stack([np.interp(t, s, c[:, 0]), np.interp(t, s, c[:, 1])])


def min_gap(A, B, d_lo, d_hi, step=500.0):
    """Minimum 3-D distance between ruled surfaces A and B over depths
    [d_lo, d_hi], sampled every step.  One of A, B must be vertical, so the
    distance at depth d is the horizontal distance of the slices (the vertical
    surface has the same slice at every depth)."""
    assert not A.rate.any() or not B.rate.any()
    out = []
    for d in np.arange(d_lo, d_hi + 0.5 * step, step):
        g = seg_dist(densify(A.at(d)), B.at(d)).min()
        out.append((float(d), float(g)))
    return out


def point_gaps(s, others, d_hi, step=500.0):
    """For each point of the densified chain of surface s: the smallest
    distance, over depths 0..d_hi sampled every step, to any of the vertical
    surfaces in others.  Returns (arc length along s, gap)."""
    P = densify(s.chain)
    g = np.full(len(P), np.inf)
    reach = np.abs(s.rate).max() * d_hi + 50e3
    lo, hi = P.min(0) - reach, P.max(0) + reach
    near = []
    for o in others:
        # the part of the other chain inside the reach box (and its neighbours)
        m = np.all((o.chain >= lo) & (o.chain <= hi), axis=1)
        m = m | np.roll(m, 1) | np.roll(m, -1)
        if m.any():
            idx = np.flatnonzero(m)
            near.append(Surface(o.name, o.chain[idx.min():idx.max() + 1], o.dip, o.dipdir, o.ref))
    for d in np.arange(0.0, d_hi + 0.5 * step, step):
        Pd = P + s.rate * (d - s.ref)
        for o in near:
            g = np.minimum(g, seg_dist(Pd, o.at(d)))
    return arclen(P), g


def cut_saf2(saf2, saf1, saf3, gap, d_hi, step=50.0):
    """Shortest cuts (start, end) of the saf2 chain, by arc length on the
    chain and to a multiple of step, so that saf2 stays >= gap from saf1 and
    saf3 at every depth in [0, d_hi].  The gap near one end does not depend
    on the cut at the other end (the chain is 56 km long), so each end is
    found by bisection with the other end at a full 40 % cut."""
    L = arclen(saf2.chain)[-1]
    big = 0.4 * L

    def ok(cuts, end):
        s = Surface("saf2", trim(saf2.chain, *cuts), saf2.dip, saf2.dipdir)
        a, g = point_gaps(s, (saf1, saf3), d_hi)
        half = 0.5 * a[-1]
        sel = a < half if end == 0 else a >= half
        return bool((g[sel] >= gap).all())

    cuts = []
    for end in (0, 1):
        other = [big, big]
        mk = (lambda c: [c, big]) if end == 0 else (lambda c: [big, c])
        lo, hi = 0, int(round(big / step))
        if ok(mk(0.0), end):
            cuts.append(0.0)
            continue
        if not ok(mk(hi * step), end):
            sys.exit("saf2 cut-back exceeds 40 % of its length")
        while hi - lo > 1:
            m = (lo + hi) // 2
            if ok(mk(m * step), end):
                hi = m
            else:
                lo = m
        cuts.append(hi * step)
    # the final pair must pass as a whole
    s = Surface("saf2", trim(saf2.chain, *cuts), saf2.dip, saf2.dipdir)
    if not (point_gaps(s, (saf1, saf3), d_hi)[1] >= gap).all():
        sys.exit(f"saf2 cuts {cuts} fail the gap {gap} m together")
    return cuts


# ----------------------------------------------------------------------------
# gmsh
# ----------------------------------------------------------------------------
def gmsh_box(gmsh, corners, depth):
    """Tilted box: corners (4, 2) counter-clockwise, z from 0 down to -depth.
    Returns (volume tag, top surface, bottom surface, side surfaces)."""
    g = gmsh.model.geo
    top = [g.addPoint(x, y, 0.0) for x, y in corners]
    bot = [g.addPoint(x, y, -depth) for x, y in corners]
    lt = [g.addLine(top[i], top[(i + 1) % 4]) for i in range(4)]
    lb = [g.addLine(bot[i], bot[(i + 1) % 4]) for i in range(4)]
    lv = [g.addLine(top[i], bot[i]) for i in range(4)]
    s_top = g.addPlaneSurface([g.addCurveLoop(lt)])
    s_bot = g.addPlaneSurface([g.addCurveLoop(lb)])
    sides = []
    for i in range(4):
        j = (i + 1) % 4
        sides.append(g.addPlaneSurface([g.addCurveLoop([lt[i], lv[j], -lb[i], -lv[i]])]))
    vol = g.addVolume([g.addSurfaceLoop([s_top, s_bot] + sides)])
    return vol, s_top, s_bot, sides


def gmsh_fault(gmsh, surf, d_fault, d_creep):
    """Friction surface 0..d_fault and creep surface d_fault..d_creep of a
    ruled surface, by translation extrusion of a spline through the top
    slice.  Returns (top curve, friction surfaces, creep surfaces)."""
    g = gmsh.model.geo
    top = surf.at(0.0)
    pts = [g.addPoint(x, y, 0.0) for x, y in top]
    c0 = g.addSpline(pts)
    out = g.extrude([(1, c0)], *surf.shift(0.0, d_fault))
    c1 = out[0][1]
    fric = [out[1][1]]
    creep = []
    if d_creep > d_fault:
        out2 = g.extrude([(1, c1)], *surf.shift(d_fault, d_creep))
        creep = [out2[1][1]]
    return c0, fric, creep


def mesh_and_write(gmsh, args, vol, s_top, s_bot, sides, tops, fric, creep, out, pairs=()):
    """pairs: (friction+creep surfaces of A, of B, smallest gap) for each pair
    of fault surfaces closer than 2 h_f somewhere; there the size is capped
    at 0.4 times the local distance between A and B (at least 0.4 gap), so the
    gap holds >= 2 elements."""
    g = gmsh.model.geo
    g.synchronize()
    gmsh.model.mesh.embed(1, tops, 2, s_top)
    gmsh.model.mesh.embed(2, fric + creep, 3, vol)
    gmsh.model.addPhysicalGroup(3, [vol], 1)
    gmsh.model.addPhysicalGroup(2, [s_top], ATTR_TOP)
    gmsh.model.addPhysicalGroup(2, [s_bot], ATTR_BOTTOM)
    gmsh.model.addPhysicalGroup(2, fric, ATTR_FAULT)
    if creep:
        gmsh.model.addPhysicalGroup(2, creep, ATTR_CREEP)
    gmsh.model.addPhysicalGroup(2, sides, ATTR_WALL)

    f = gmsh.model.mesh.field
    fd = f.add("Distance")
    f.setNumbers(fd, "SurfacesList", fric)
    f.setNumber(fd, "Sampling", args.sampling)
    expr = f"Min({args.hf} + {args.grade} * Max(F{fd} - {args.d0}, 0), {args.hmax})"
    for sa, sb, g in pairs:
        fa, fb = f.add("Distance"), f.add("Distance")
        f.setNumbers(fa, "SurfacesList", sa)
        f.setNumbers(fb, "SurfacesList", sb)
        f.setNumber(fa, "Sampling", args.sampling)
        f.setNumber(fb, "Sampling", args.sampling)
        expr = f"Min({expr}, Max({0.4 * g}, 0.4 * (F{fa} + F{fb})))"
    fm = f.add("MathEval")
    f.setString(fm, "F", expr)
    f.setAsBackgroundMesh(fm)
    o = gmsh.option
    o.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    o.setNumber("Mesh.MeshSizeFromPoints", 0)
    o.setNumber("Mesh.MeshSizeFromCurvature", 0)
    o.setNumber("Mesh.Algorithm", 6)        # Frontal-Delaunay 2-D
    o.setNumber("Mesh.Algorithm3D", 10)     # HXT (one thread); Delaunay (1) left a 9.1 deg tet in the 1 km mesh (M2)
    o.setNumber("General.NumThreads", 1)
    o.setNumber("Mesh.Optimize", 1)
    o.setNumber("Mesh.OptimizeThreshold", 0.45)   # default 0.3 left one 8.7 deg sliver at a fault (M2)
    o.setNumber("Mesh.MeshSizeMax", args.hmax)
    o.setNumber("Mesh.MshFileVersion", 2.2)
    o.setNumber("Mesh.Binary", 0)
    o.setNumber("Mesh.SaveAll", 0)
    o.setNumber("General.Verbosity", 3)
    t0 = time.time()
    gmsh.model.mesh.generate(3)
    t_mesh = time.time() - t0
    gmsh.write(out)
    return t_mesh


def record_counts(gmsh):
    """Element counts per physical group."""
    counts = {}
    for dim, tag in gmsh.model.getPhysicalGroups():
        n = 0
        for e in gmsh.model.getEntitiesForPhysicalGroup(dim, tag):
            types, tags, _ = gmsh.model.mesh.getElements(dim, e)
            n += sum(len(t) for t in tags)
        counts[f"{dim}:{tag}"] = n
    return counts


# ----------------------------------------------------------------------------
# jobs
# ----------------------------------------------------------------------------
def job_ucerf(args):
    import gmsh
    if sha256(SRC_2D) != SRC_2D_SHA:
        sys.exit(f"{SRC_2D}: sha256 is not the committed {SRC_2D_SHA}")
    ch = name_chains(read_chains_2d(SRC_2D))
    box = json.load(open(SRC_2D_JSON))["tilt"]["corners"]
    S = {n: Surface(n, c) for n, c in ch.items()}
    S["saf2"] = Surface("saf2", ch["saf2"], SAF2_DIP_DEG, SAF2_DIPDIR_DEG)
    gap = args.gap
    d_all = max(args.fault_depth, args.creep_depth)
    cuts = cut_saf2(S["saf2"], S["saf1"], S["saf3"], gap, d_all)
    S["saf2"] = Surface("saf2", trim(ch["saf2"], *cuts), SAF2_DIP_DEG, SAF2_DIPDIR_DEG)

    # gaps between every pair (one of each pair is vertical)
    gaps = {}
    names = list(S)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            A, B = (S[a], S[b]) if not S[b].rate.any() else (S[b], S[a])
            prof = min_gap(A, B, 0.0, d_all)
            gaps[f"{a}-{b}"] = {"min_km": min(g for _, g in prof) / 1e3,
                                "at_depth_km": min(prof, key=lambda x: x[1])[0] / 1e3,
                                "profile_km": [[d / 1e3, g / 1e3] for d, g in prof]}

    gmsh.initialize(["", "-nt", "1"])
    gmsh.model.add("ucerf3d")
    vol, s_top, s_bot, sides = gmsh_box(gmsh, box, args.depth)
    tops, fric, creep, tagmap = [], [], [], {}
    for n in ("saf1", "saf2", "saf3", "sjf1", "sjf2"):
        c0, f_, c_ = gmsh_fault(gmsh, S[n], args.fault_depth, args.creep_depth)
        tops.append(c0)
        fric += f_
        creep += c_
        tagmap[n] = {"friction": f_, "creep": c_}
    pairs, refined = [], []
    for k, v in gaps.items():
        if v["min_km"] * 1e3 < 2.0 * args.hf:
            a, b = k.split("-")
            pairs.append((tagmap[a]["friction"] + tagmap[a]["creep"],
                          tagmap[b]["friction"] + tagmap[b]["creep"], v["min_km"] * 1e3))
            refined.append(k)
    t_mesh = mesh_and_write(gmsh, args, vol, s_top, s_bot, sides, tops, fric, creep, args.out, pairs)
    counts = record_counts(gmsh)
    gmsh.finalize()

    rec = {
        "what": "UCERF 3-D SAF/SJF mesh, plan rev5 P0",
        "command": " ".join(sys.argv),
        "gmsh": gmsh.__version__,
        "source_2d_mesh": {"path": os.path.relpath(SRC_2D, SEAS), "sha256": SRC_2D_SHA},
        "chains": {n: {"vertices": int(len(S[n].chain)),
                       "length_km": float(arclen(S[n].chain)[-1] / 1e3),
                       "dip_deg": S[n].dip,
                       "dipdir_deg": S[n].dipdir if S[n].dip < 90 else None,
                       "slice_depth_km": SLICE_DEPTH / 1e3,
                       "top_trace_m": S[n].at(0.0).tolist()} for n in S},
        "saf2_cut_km": {"saf3_end": cuts[0] / 1e3, "saf1_end": cuts[1] / 1e3,
                        "length_before_km": float(arclen(ch["saf2"])[-1] / 1e3)},
        "gap_target_km": gap / 1e3,
        "refined_pairs": refined,
        "gaps": gaps,
        "box": {"corners_m": box, "depth_km": args.depth / 1e3},
        "depths_km": {"friction": args.fault_depth / 1e3, "creep_bottom": args.creep_depth / 1e3},
        "size": {"hf_m": args.hf, "d0_m": args.d0, "grade": args.grade, "hmax_m": args.hmax},
        "attributes": {"top": ATTR_TOP, "bottom": ATTR_BOTTOM, "friction": ATTR_FAULT,
                       "walls": ATTR_WALL, "creep": ATTR_CREEP, "volume": 1},
        "gmsh_surface_tags": tagmap,
        "counts": counts,
        "mesh_seconds": t_mesh,
        "msh": os.path.basename(args.out),
        "msh_sha256": sha256(args.out),
    }
    json.dump(rec, open(os.path.splitext(args.out)[0] + ".json", "w"), indent=1)
    print(json.dumps({k: rec[k] for k in ("saf2_cut_km", "gap_target_km", "counts", "mesh_seconds", "msh_sha256")}, indent=1))


TESTS = {
    # name: (fault polylines [(chain, dip)], box half width, band)
    "T-plane": dict(dip=90.0, bend=0.0, box=120e3, band=0.0),
    "T-dip": dict(dip=58.0, bend=0.0, box=120e3, band=0.0),
    "T-bend": dict(dip=90.0, bend=20.0, box=120e3, band=0.0),
    "T-pml": dict(dip=90.0, bend=0.0, box=120e3, band=30e3),
}


def job_test(args):
    import gmsh
    T = TESTS[args.name]
    L, W = 40e3, 20e3
    if T["bend"]:
        a = math.radians(T["bend"])
        chain = np.array([[-L / 2, 0.0], [0.0, 0.0], [L / 2 * math.cos(a), L / 2 * math.sin(a)]])
        # resample so the spline is the polyline (a corner at x = 0)
        n = int(L / 2 / 1e3)
        chain = np.vstack([np.column_stack([np.linspace(-L / 2, 0, n + 1), np.zeros(n + 1)]),
                           np.column_stack([np.linspace(0, L / 2 * math.cos(a), n + 1)[1:],
                                            np.linspace(0, L / 2 * math.sin(a), n + 1)[1:]])])
    else:
        chain = np.array([[-L / 2, 0.0], [L / 2, 0.0]])
    # dipping fault: dip direction +y, top trace at y = 0
    surf = Surface(args.name, chain, T["dip"], 0.0, ref_depth=0.0)
    half = T["box"] / 2 + T["band"]
    depth = T["box"] / 2 + T["band"]
    corners = [[-half, -half], [half, -half], [half, half], [-half, half]]
    gmsh.initialize(["", "-nt", "1"])
    gmsh.model.add(args.name)
    vol, s_top, s_bot, sides = gmsh_box(gmsh, corners, depth)
    g = gmsh.model.geo
    pts = [g.addPoint(x, y, 0.0) for x, y in surf.at(0.0)]
    tops, fric, creep = [], [], []
    if T["bend"]:
        # two straight segments: a polyline of lines keeps the corner sharp
        ip = len(pts) // 2
        c0 = [g.addLine(pts[0], pts[ip]), g.addLine(pts[ip], pts[-1])]
    else:
        c0 = [g.addLine(pts[0], pts[1])]
    for c in c0:
        out = g.extrude([(1, c)], *surf.shift(0.0, W))
        fric.append(out[1][1])
        if args.creep_depth > W:
            out2 = g.extrude([(1, out[0][1])], *surf.shift(W, args.creep_depth))
            creep.append(out2[1][1])
        tops.append(c)
    t_mesh = mesh_and_write(gmsh, args, vol, s_top, s_bot, sides, tops, fric, creep, args.out)
    counts = record_counts(gmsh)
    gmsh.finalize()
    rec = {"what": f"{args.name} test mesh, plan rev5 section 10.1", "command": " ".join(sys.argv),
           "gmsh": gmsh.__version__, "fault": {"length_km": L / 1e3, "depth_km": W / 1e3,
           "dip_deg": T["dip"], "dipdir": "+y" if T["dip"] < 90 else None, "bend_deg": T["bend"],
           "creep_bottom_km": args.creep_depth / 1e3, "top_trace_m": surf.at(0.0).tolist()},
           "box": {"interior_km": T["box"] / 1e3, "band_km": T["band"] / 1e3,
                   "half_width_km": half / 1e3, "depth_km": depth / 1e3},
           "size": {"hf_m": args.hf, "d0_m": args.d0, "grade": args.grade, "hmax_m": args.hmax},
           "attributes": {"top": ATTR_TOP, "bottom": ATTR_BOTTOM, "friction": ATTR_FAULT,
                          "walls": ATTR_WALL, "creep": ATTR_CREEP, "volume": 1},
           "counts": counts, "mesh_seconds": t_mesh,
           "msh": os.path.basename(args.out), "msh_sha256": sha256(args.out)}
    json.dump(rec, open(os.path.splitext(args.out)[0] + ".json", "w"), indent=1)
    print(json.dumps({k: rec[k] for k in ("counts", "mesh_seconds", "msh_sha256")}, indent=1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="job", required=True)
    for name in ("ucerf", "test"):
        p = sub.add_parser(name)
        p.add_argument("--out", required=True)
        p.add_argument("--hf", type=float, default=1000.0, help="fault element size h_f (m)")
        p.add_argument("--d0", type=float, default=5000.0, help="h = h_f within d0 of a fault (m)")
        p.add_argument("--grade", type=float, default=0.2, help="growth of h per metre beyond d0")
        p.add_argument("--hmax", type=float, default=20000.0, help="cap on h (m)")
        p.add_argument("--sampling", type=int, default=200, help="gmsh Distance field sampling")
        if name == "ucerf":
            p.add_argument("--fault_depth", type=float, default=25e3)
            p.add_argument("--creep_depth", type=float, default=50e3, help="option A creep bottom (E4)")
            p.add_argument("--depth", type=float, default=415e3, help="box depth (m)")
            p.add_argument("--gap", type=float, default=2000.0,
                           help="minimum surface gap; where two surfaces are closer than 2 h_f the mesh "
                                "is refined to 0.4 times their distance, so there are >= 2 elements across")
        else:
            p.add_argument("--name", required=True, choices=sorted(TESTS))
            p.add_argument("--creep_depth", type=float, default=30e3,
                           help="creep bottom (m); the friction fault is 0-20 km")
    args = ap.parse_args()
    (job_ucerf if args.job == "ucerf" else job_test)(args)


if __name__ == "__main__":
    main()
