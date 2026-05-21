# DEBUG: z0cut at finer near-fault resolution (500 m, 250 m) — 2026-05-20

Goal: build z0embed meshes like `safs_fault_box_nwcut_1000m_lcfar3000_z0embed.msh`
but with **near-fault resolution 500 m and 250 m**, keeping the same far-field
coarsening (lc_far=3000, dist_inner=3000, dist_outer=40000) to avoid an element
explosion, and compare quality across 1000/500/250 m.

Builds on `DEBUG_z0cut_option_a_failure.md` (the geometric-trace-embed + weld fix
that made z0cut work at 1000 m).

## Key constraint: z0cut does NOT re-mesh the fault

The fault is a *discrete* (merged-STL) surface; `run_z0cut_meshing.py` keeps its
triangulation verbatim (embed + weld; 10603 tris preserved at 1000 m). So the
**fault-triangle size = STL resolution**, and the near-fault *bulk* size = the
size field's `lc_near`.  To get a true 500/250 m near-fault mesh BOTH must be set:

| near-fault res | fault STL (fault-tri size) | lc_near (bulk) |
|---|---|---|
| 1000 m | `…_1000m_clean_clip_nwcut.stl` (10,603 tris) | 1500 (the lcfar3000 recipe) |
| 500 m  | `…_500m_clean_clip_nwcut.stl` (42,358 tris)  | 500 |
| 250 m  | **subdivided** (see below; 169,432 tris)     | 250 |

### No 250 m STL exists → geometry-preserving 1→4 subdivision of the 500 m STL
`raw/` has only 500/1000/2000 m `.ts`/`.stl`, and the STL-prep pipeline is
archived.  A **midpoint (1→4) subdivision** of the 500 m STL gives ~250 m tris:
midpoints lie on the parent facets, so the **surface and the z=0 trace are exact**
(no geometric change), and each sub-triangle is *similar* to its parent, so the
**shape quality (min_angle) is identical** — only edge length halves.  Verified:
500 m STL min_angle min/p10/med = 9.1/53.5/58.9° ; 250 m subdiv = **identical**.
Output: `stl_nwcut/…_250m_subdiv_clean_clip_nwcut.stl` (1199 nodes on z=0).
(Consequence: the 250 m fault-tri *shape* quality equals the 500 m by construction;
only the size and the bulk differ.)

## The wall: exact-z=0 PLC degeneracy at fine resolution

At 1000 m, z0cut meshes cleanly (3D in 18 s, 963k tets).  At 500 m it FAILED:
`PLC Error: A segment and a facet intersect at point` during TetGen boundary
recovery (and HXT: `a vertex lies in a segment` at a z=0 node on a trace segment).

**Feasibility matrix on the 500 m STL** (optimize off) isolated the cause:

| field | algorithm | result |
|---|---|---|
| uniform 3000 | Delaunay(1) | 780k tets **OK** |
| uniform 3000 | HXT(10)     | 664k tets **OK** |
| lc_near=500  | Delaunay(1) | **FAIL** (segment/facet intersect) |
| lc_near=500  | HXT(10)     | **FAIL** (vertex-in-segment) |

So it is **not the algorithm and not the bulk** — the *coarse* field works on the
same STL; the *fine* field fails.  Root cause: the SAFS fault dips shallowly right
at the surface trace, so some fault facets are **near-horizontal** (|n_z| up to
0.97 ≈ 14° dip, with 2 vertices on z=0).  These are nearly **coplanar with the
z=0 top face**.  A *fine* top-face mesh there places surface nodes on/near the
trace segments shared with those facets → degenerate PLC.  The 1000 m mesh
escapes it (coarser, steeper: |n_z|≤0.94, and lc_near=1500 keeps the top face
coarse).  Near-horizontal trace facets: 7 (1000 m) → 14 (500 m) → 28 (250 m).

## Fix: near-surface size buffer (Box field) — keeps EXACT z=0

`run_z0cut_meshing.py` gained `--surface-buffer-size S [--surface-buffer-depth D]`
(default off).  It Max's a **Box field** (VIn=S for z>−D, VOut=0 below, linear
transition of width D) into the background size field, so the **top ~D layer is
coarsened to S** while the seismogenic bulk below stays at lc_near.  The fault
still reaches **exactly z=0** (surface rupture preserved); only the near-surface
*bulk/top-face* is coarser — which keeps the top face away from the degeneracy
AND trims element count (the user's "avoid too many elements").  Used S=1500,
D=1000 for 500 m and 250 m.

A Box field is essential: an `exp()` MathEval buffer also clears the PLC but is
~5× slower (transcendental evaluated per refinement query: 32 min vs the Box's
few min on 500 m).

### Gotcha: launch the long runs as a DIRECT command
A `nohup bash -c "time python … &"` wrapper terminated the 500 m run at ~1m51s
mid-refinement (no mesh).  Re-launched as a **direct** background python
invocation (as the 1000 m and matrix runs were), it completed normally.

## Results

| near-fault | n_tet | fault tris (all ==2 tets) | tet_emin | eta_min | Q1 | Q2 | near-tip min_angle / tri_q |
|---|---:|---:|---:|---:|:-:|:-:|---|
| 1000 m | 963,322    | 10,603  | 150.2 m | 0.1247 | ✓  | ✓ | 13.3° / 0.385 |
| 500 m  | 3,693,060  | 42,358  | 112.5 m | 0.1838 | ✓  | ✓ | 12.7° / 0.371 |
| 250 m  | 15,367,950 | 169,432 | 56.2 m  | 0.1210 | ✗* | ✓ | 12.7° / 0.371 |

All three: **mesh_zmax = 0 exact**, **100 % fault embedded + welded** (every
fault tri borders exactly 2 tets), **0 bulk slivers** (η<0.1), η_median ≈ 0.855.

*Q1 (min tet edge ≥ 100 m) is a gate calibrated for ~1 km meshes.  At 250 m the
subdivided fault tris carry ~56 m min edges (= half the 500 m's 112 m), so
sub-100 m edges are EXPECTED; the **shape** stays healthy (η_min 0.121 > 0.1,
Q2 ✓), so this is a gate-threshold mismatch, not a sliver problem.  Lower the
Q1 floor (e.g. `--min-edge-floor-m 50`) when checking sub-km meshes.

Element growth ~3.8–4.2× per halving (16× total 1000→250 m), held down by the
coarse far field (lc_far=3000) + the surface buffer.  Fault-tri *shape* is
~constant across resolutions (STL-inherited; 250 m = subdivision of 500 m), so
the refinement buys element *size* / cohesive-zone (L_nuc/h) resolution, NOT
better fault-surface shape — the near-tip worst sliver stays ~12.7–13.3° at all
three.  (Improving that needs fault re-triangulation, per
`DEBUG_z0cut_option_a_failure.md`.)

## Reproduce
```bash
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
cd meshing/code
# 500 m
$PY run_z0cut_meshing.py --stl ../results/stl_nwcut/…_500m_clean_clip_nwcut.stl \
    --out ../results/msh/…_500m_lcfar3000_z0embed.msh \
    --lc-near 500 --surface-buffer-size 1500 --surface-buffer-depth 1000
# 250 m (subdivided STL)
$PY run_z0cut_meshing.py --stl ../results/stl_nwcut/…_250m_subdiv_clean_clip_nwcut.stl \
    --out ../results/msh/…_250m_lcfar3000_z0embed.msh \
    --lc-near 250 --surface-buffer-size 1500 --surface-buffer-depth 1000
```
