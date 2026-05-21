# DEBUG: z0cut "Option A" (gmsh trace-embed) produces an empty volume (2026-05-20)

Companion to `spatial_dynamic_rupture_postevent_debug_2026-05-20.md` (the
north-tip slip-rate blow-up) and `PLAN_mesh_quality.md` (the buried-fault
headroom fix that left residual fault-trace slivers). This doc records why the
in-progress surface-rupturing mesher `code/run_z0cut_meshing.py` fails, so the
next attempt does not re-tread it.

## Goal recap

Eliminate the fault∩free-surface **sliver triangles** (min_angle 6.4–9.9°,
aspect 6–8, tri_q 0.19–0.29) at the shallow north trace (z ≈ −141…−150 m) that
seed the blow-up. User requirement: a **surface-rupturing** fault — cut EXACTLY
at z = 0, nothing above, z = 0 IS the domain top boundary (removes the ~100 m
"headroom wedge" that the buried-fault production mesher creates and which is the
documented sliver source).

## Symptom (reproduced 2026-05-20)

`run_z0cut_meshing.py … --uniform-lc 3000` (fast feasibility pass) on
`stl_nwcut/…ALT6_1000m_clean_clip_nwcut.stl` (z ∈ [−16607, 0], 300 nodes on
z = 0):

```
Meshing 1D … done    Meshing 2D … done (6 surfaces)
Meshing 3D …
  Reconstructing mesh / Recovering boundary … done
  Warning : No tetrahedra in region 1
  Warning : No elements in volume 1
RESULT: volume tets = 0; mesh_zmax = 0.000000
ERROR: volume 1 is EMPTY (rc=2)
```

1D and 2D succeed; the z = 0 trace is found (299 segments / 300 nodes); 3D
boundary recovery runs but encloses **no volume** → 0 tets.

## Root cause (CONFIRMED by instrumentation — /tmp/diag_z0.py)

Reproduced the exact `main()` setup, stopped after `generate(2)`, and measured
trace conformity between the fault and the OCC box top face:

| check | result | meaning |
|---|---|---|
| (a) trace nodes present in **top-face** node set | **17 / 300** | top face barely touches the trace |
| (b) trace segments that are **top-face triangle edges** | **0 / 299** | top face does **not** conform to the trace at all |
| (c) coincident-but-distinct nodes near trace | 0 | not a duplicate-node "crack" |
| (d) trace tags still in **fault-surface** node set after `generate(2)` | **23 / 300** | discrete fault was **re-tagged/re-meshed** during meshing |

**Mechanism.** The script builds the trace as a *discrete* curve
(`addDiscreteEntity(1)` + `addElementsByType(LINE, …)`) referencing fault node
tags captured **before** meshing, then `embed(1, trace_curve, 2, top)`. Two
things break it:

1. **Node tags are not stable across `generate()`.** `gmsh.model.mesh.generate(2)`
   regenerates the discrete (STL-merged) fault surface's nodes — only 23/300 of
   the captured z = 0 tags survive (d). The trace curve's line elements then
   reference mostly **stale/dangling** tags.
2. **A pre-meshed discrete curve embedded into an OCC (CAD) surface does not
   constrain the OCC 2-D mesher.** The top face is re-meshed from its OCC
   parametrization and ignores the discrete trace (0/299 edges, (b)).

Net: the fault "curtain" reaches z = 0 but its top edge is **not welded** to the
top-face triangulation. The volume boundary has an open seam along the trace →
not watertight → Delaunay-3D recovers no enclosed region → 0 tets. This is a
geometry/topology conformity failure, independent of size field or resolution.

## Why a small patch is unlikely to rescue Option A as written

- Marking the fault "do not re-mesh" stabilizes tags (fixes (d)) but **(b) still
  stands**: embedding a discrete curve into an OCC face doesn't force the OCC
  mesher to conform. The seam stays open.
- The conformity has to be established **geometrically/topologically before
  meshing**, not by post-hoc embedding of a pre-meshed curve. That means a
  different mechanism (see candidates), i.e. an approach change — flagged for the
  user per project debugging rules rather than silently swapped in.

## Candidate paths (path 1 CHOSEN by user 2026-05-20 — see FIX below)

1. **Geometric-trace embed + weld (smallest change to z0cut).** Replace the
   discrete trace with *geometric* Points+Lines at the 300 z = 0 coords, embed
   those into the top OCC face (CAD embed IS honored → top face conforms), keep
   the fault discrete and **don't** re-mesh it, then `removeDuplicateNodes()` to
   weld the coincident top-face/fault z = 0 nodes before `generate(3)`. Plausible
   but still fragile (relies on the discrete fault's z = 0 nodes landing exactly
   on the geometric trace so the weld closes the seam).

2. **CGAL corefinement (already scoped: `PLAN_corefine_pipeline.md`,
   `EXPLORE_cgal_corefine.md`, `cgal-61` env).** Boolean-corefine the fault into
   the box so the trace becomes a **real shared edge by construction** →
   guaranteed conforming → then tetrahedralize. Highest robustness for "fault
   reaches the surface"; larger upfront effort.

3. **Stay buried, raise headroom to kill slivers (contradicts the z = 0
   requirement).** Raise `--fault-top-clamp` toward the 250 m physics ceiling
   (documented as ≪ seismogenic length, so physically acceptable) and/or get the
   Netgen optimizer working (currently SIGBUSes) or a TMOP polish
   (`EXPLORE_tmop_mesh_optimization.md`). Lowest effort to remove the blow-up
   seed, but gives up co-seismic surface rupture.

## Reproduce

```bash
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python   # NOTE: env is
   # "pythonenv" (gmsh 4.15, meshio 5.3.5, numpy 2.3.3); the docs' generic
   # "conda activate pythonenv" still applies — it exists, just not in the
   # truncated `conda env list | head`.
cd meshing/code
$PY run_z0cut_meshing.py \
    --stl ../results/stl_nwcut/SAFS-…-ALT6_1000m_clean_clip_nwcut.stl \
    --out /tmp/z0embed_feasibility.msh --uniform-lc 3000     # → rc=2, empty
# conformity instrumentation: /tmp/diag_z0.py  (stops at generate(2))
```

## FIX APPLIED + VERIFIED (path 1, 2026-05-20)

`run_z0cut_meshing.py` updated:
- Trace is now **geometric** OCC points+lines at the z = 0 trace coords,
  embedded in the top face (`embed(1, trace_lines, 2, top)`); the fault stays
  discrete (NOT re-meshed — 10603 tris preserved).
- Mesh sequence is `generate(2)` → `removeDuplicateNodes()` → `generate(3)`;
  the weld merges exactly the 300 coincident top-face/fault z = 0 nodes so the
  boundary is watertight.  (gmsh internally re-meshes 2-D surface *interiors*
  inside `generate(3)`, but the 1-D welded trace nodes are preserved, so the
  seam holds — confirmed by the 100 % embedding check below.)

Production mesh `results/msh/safs_fault_box_nwcut_1000m_lcfar3000_z0embed.msh`
(963,322 tets, 51.7 MB):

| metric | production (zclamp) | **z0embed** | gate / target |
|---|---|---|---|
| mesh_zmax | 0 | **0** | = 0 exact ✓ |
| fault embedded (tris bordering 2 tets) | n/a | **10603/10603 (100 %)** | watertight ✓ |
| bulk tet_emin | (varies) | **150.2 m** | ≥ 100 m (Q1) ✓ |
| bulk eta_min | residual η 0.03–0.09 | **0.1247** | > 0.1 (Q2) ✓ |
| bulk slivers (η<0.1) | present | **0** | 0 |
| north-tip fault min_angle (worst) | 8.67° | **13.27°** | ≳ 25–30° (not met) |
| north-tip fault aspect (worst) | 6.34 | **4.27** | ≲ 3 (not met) |
| north-tip fault tri_q (worst) | 0.257 | **0.385** | — |
| worst trace sliver mesh-wide | 6.40° (z=−141) | **11.50°** (z=−74) | — |

**Net:** the headroom-wedge **bulk** slivers and the zclamp-induced
**fault-trace** slivers (the documented blow-up seed) are gone / much improved;
the surface-rupturing z = 0 cut and all hard invariants hold.  **Residual**
fault-trace slivers remain (min_angle 11.5–15°, tri_q 0.34–0.43) plus one deep
z = −1903 m STL sliver (9.43°, tri_q 0.25) — these are **STL-inherent** (the
fault triangulation is not re-meshed) and are still short of the conservative
25–30° / aspect ≤ 3 target.  Closing that gap needs **fault re-triangulation**
(`PLAN_retriangulate.md`), a separate effort, OR a re-run to confirm the
improved (tri_q 0.39 vs 0.26) tip is already enough to suppress the blow-up.

Reproduce the verified mesh:
```bash
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
cd meshing/code
$PY run_z0cut_meshing.py \
    --stl ../results/stl_nwcut/SAFS-…-ALT6_1000m_clean_clip_nwcut.stl \
    --out ../results/msh/safs_fault_box_nwcut_1000m_lcfar3000_z0embed.msh
$PY check_mesh_quality.py ../results/msh/…_z0embed.msh        # Q1 ✓ Q2 ✓
python3 locate_fault_slivers.py ../results/msh/…_z0embed.msh \
    --target 448338 3798440 -249.7                            # tip tri quality
```

## Status

- CONFIRMED: production slivers; z0cut Option-A empty-volume failure + root
  cause (instrumented); **path-1 fix works** — watertight, fully embedded,
  Q1/Q2 pass, blow-up-seed slivers removed (data above).
- OPEN: (a) whether the residual STL-inherent fault slivers need
  re-triangulation, or the improved tip already suppresses the blow-up (decide
  via a short re-run); (b) regenerate 500 m / 2000 m z0embed variants if needed
  (only 1000m_lcfar3000 done — the variant tied to the blow-up).
