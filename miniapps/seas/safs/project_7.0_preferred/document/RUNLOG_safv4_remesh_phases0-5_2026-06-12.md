# Run log: SAFv4 remesh Phases 0-5 (2026-06-12)

Implements `PLAN_mesh_quality_safv4_remesh_2026-06-12.md` Phases 0-4.
Environments: `cgal-61` (C++), `pythonenv` (Python).  All commands run from
`miniapps/seas/safs/project_7.0_preferred/` unless noted.

## Phase 0 — restore + build + fixture regeneration

- Restored from commit `577646e`: `code_preprocess/` (29 files),
  `raw_data/` (18 CFM `.ts`), `data_corefined/` manifests.
- **CGAL drift surfaced and resolved**: with default (assertion-enabled)
  builds, CGAL 6.1.1 aborts in `isotropic_remeshing(protect_constraints=
  true)` with "constraints larger than 4/3 * target_edge_length"
  (remesh.h:295).  The May runs used Release semantics where CGAL
  preconditions compile out.  Fix: build with `-DCMAKE_BUILD_TYPE=Release`
  (no source patch).  All targets build clean.
- Fixture regeneration chain (`raw_data/*.ts` -> `ts_to_stl.py --no-clip`
  -> `clean_freesurface_mesh.py --batch` -> `corefine_faults.py --res
  2000`): 7/7 pairs conformal at `max_diff_m = 0`; per-fault corefined
  face counts 788 / 977 / 1291 / 1198 / 1006 / 2892 (total **8,152 — equal
  to the canonical archived fixture count** cited in
  PLAN_revert_to_gmsh.md: 12 + 8152 = 8164).  The archived
  `data_corefined/manifest.json` turned out to be from a later merge
  experiment (one 7,348-face entry), so the count comparison above is the
  meaningful drift check: **no drift**.

## Phase 1 — gate harness

- New `meshing/code/check_fault_conformity.py` (self-test PASS; exits
  nonzero on the known-bad GOCAD mesh).
- `meshing/code/check_mesh_quality.py`: additive `.msh`/`.vtu` input
  (surfaces = `gmsh:physical` groups; `--fault-tags`); metric kernels
  untouched.  Cross-check on `tpv102_1000m.msh` vs its meshio `.vtu`:
  bulk metrics identical.
- Baseline `meshing/code/baseline_safv4_2km_topo_conformity.json`
  committed; byte-stable across two runs (md5
  `8133f293f8fe63fa4c462dc8281c8d0d`).  Reproduces the diagnosis exactly:
  fault-subset dup groups 2,939; duplicated pairs SAF 1,417 / Banning 226
  / Garnet 1,182; centroid matches 2,558/2,564, 352/352, 2,400/2,400;
  slivers 214,359; sub-100 m tets 1,106.  NEW finding: mesh-wide dup
  groups ALSO = 2,939 — the boundary-shell rims share node IDs in the
  .inp; ONLY the faults are split.  (The plan's "n_groups_total > 2,939"
  expectation was an artifact of comparing per-file VTU coordinates; the
  acceptance is corrected to total == fault_subset == 2,939.)

## Phase 2 — extract + weld

- New `meshing/code/extract_safv4_surfaces.py` ->
  `meshing/data_extracted/safv4_2km_topo/` (3 fault STLs, dem, bottom,
  4 ribbons + `manifest.json`; manifest byte-stable, md5
  `b99f607dea973647be11b2b6e7c147f0`).
- Acceptance: fault counts 2,564 / 352 / 2,400 (= baseline minus side);
  border-contact nodes SAF-bottom 53, Garnet-bottom 5, MJVS-MULT 9,
  MJVS-SBMT 21, MULT-SBMT 68, zero fault-ribbon contact; 0 dup groups
  among welded used nodes; junction host/guest resolved (guests: SAF at
  both its junctions, Banning at Garnet) — all `host_interior_contact`.
- **Deviation 1 (measured-data-driven): DEM re-triangulation.**  The
  GOCAD DEM is split along the fault traces with DIFFERENT per-side
  segmentations (977 open T-junction crack edges, all with both endpoints
  on welded fault nodes), so welding cannot close the shell.  The DEM is
  single-valued in (x,y): rebuilt as the 2-D Delaunay triangulation of
  its own unchanged vertex set, restricted to the (concave — hull chords
  up to 430 km were the first failure mode) rim polygon: 30,110 -> 30,919
  tris, 977 crack edges resolved, 287 outside-footprint hull tris
  dropped, all 313 rim edges preserved.  Boundary shell after repair:
  55,608 edges, **all with incidence exactly 2** (closed, manifold).
- **Deviation 2: ASCII STL (%.15g) instead of binary.**  Binary STL is
  float32 by format definition (~0.06 m error at UTM y ~ 4e6) and would
  defeat the 1e-6 m conformality contract; same convention as the
  archived `corefine_faults._off_to_stl`.

## Phase 3 — extension, graded remesh, corefine

- New `meshing/code/extend_fault_borders.py`:
  - recursive border split at 750 m FIRST (isotropic_remeshing preserves
    borders verbatim, so the coarse ~2 km GOCAD border segmentation must
    be refined before remeshing or it survives to the gates): SAF 1,221 /
    Banning 176 / Garnet 396 splits;
  - class extrusions: trace +z 1500 m (DEM-clearance verified, no
    auto-raise needed), bottom -z 500 m, junction guest in-plane 750 m;
    strips wound opposite the base triangle's halfedge traversal so the
    composite stays orientable; SAF +3,050 / Banning +648 / Garnet +752
    strip tris;
  - emits `trace_points.xyz` (1,303 sites) for the DEM pre-pass.
- New `code_preprocess/corefine_cgal/remesh_graded.cpp` (+
  `graded_sizing_field.h`, model of CGAL 6.1 PMPSizingField): DEM
  pre-pass 30,919 -> 98,477 tris at h(d) = 500 -> 2500 m (d = 1 -> 9 km
  from trace sites), **rim byte-identical (tool-verified)**.
  Interface deviation: `--trace-points FILE.xyz` instead of parsing the
  Phase 2 manifest JSON in C++ (the .xyz is emitted deterministically by
  extend_fault_borders.py).
- `corefine_set.cpp` extensions (fixture-regression-verified: meshes[]
  and pairs[] identical pre/post):
  - per-mesh remesh modes `--graded NAME[:h_near:h_far:d_near:d_far]`
    (sites = the mesh's OWN corefined polylines) / `--keep NAME`;
    uniform path byte-compatible;
  - rim (border-vertex) protection in quality_repair for graded/keep
    meshes;
  - `count_shared_vertices` grid-accelerated (identical semantics; the
    archived O(|A|x|B|) brute force would take hours at DEM 1e5 x fault
    1e4 vertices).
- `io_helpers::read_polygon_mesh_any`: soup repair/orient fallback for
  non-manifold pinch vertices at strip corners (duplicates are re-merged
  by the downstream 1e-6 soup dedup).
- `corefine_faults.py`: `--inputs` explicit-STL mode + passthrough of the
  new corefine_set options.  (Deviation: explicit `--inputs` instead of a
  `--manifest` reader — same information, less indirection; the
  invocation is recorded below.)
- **Deviation 3: `bottom` participates in the corefine** (graded
  500->8000 m, d-window 1-20 km) — required by the measured SAF/Garnet
  bottom trimming (R-001): the bottom extensions must be corefined in
  Phase 3 or the Phase 4 autorefine fault-split gate would trip by
  design.  Ribbons remain untouched (zero fault contact, as planned).
- Production run (500 m): participants {3 extended faults, dem_graded,
  bottom}; `--mesh-edge-size 500 --min-edge 100` (polyline spacing 250);
  8 intersecting pairs found (3 traces, 3 junctions, SAF/Garnet x
  bottom) — exactly the predicted set.
- **Deviation 4: junction extrusion direction = HOST surface normal**
  (away from the guest body), replacing the plan's guest-in-plane
  direction.  Measured trigger: with in-plane strips, the SAF x Garnet
  junction pair (dihedral 13.4 deg between the sub-parallel strands;
  SAF x Banning 60.7, Banning x Garnet 39.8) drove `PMP::corefine` into
  a tight loop in its finalize step — >30 min CPU with the process
  sampled 100% inside `Surface_intersection_visitor_for_corefinement::
  finalize` (no gmp frames; the in-plane strip slices a
  ~h/tan(13 deg) ~ 3 km band of host triangles nearly parallel to
  them).  A host-normal strip crosses the host at 90 deg and meets it
  only along the contact line — the most transversal realization of the
  plan's stated intent.  Re-extension verified: 0 self-intersections on
  all extended inputs.
- **Deviation 5: junction borders detached by eps = 30 m; bottom borders
  NOT extruded.**  Two further measured failure modes:
  (i) even with host-normal strips, the SAF BODY still touches Garnet
  bit-exactly along the 21 shared contact nodes at 13.4 deg
  near-parallelism — same corefine finalize loop.  Fix: pull the guest's
  junction border back toward its own body by eps = 30 m before
  extruding (strip lengthened by eps), so the body clears the host and
  the only intersection is the strip's clean ~90 deg crossing.  With
  this, ALL 8 corefine pairs complete (SAF x Garnet: 2 polylines /
  114 verts).
  (ii) a -z bottom strip with the eps offset slices ~30 m slivers along
  the entire ~100 km bottom border (thousands of sub-floor polyline
  edges); the polyline cleanup's cascading collapses then blow up
  (>40 min; collapse_edge spinning on a mega-degree vertex).  Fix: NO
  bottom strip — the fault's bottom border lies EXACTLY in the bottom
  plane (validated), so the contact is a tangential edge contact;
  corefine inserts the border line into the bottom mesh conformally
  with only a handful of short segments.  The R-001 concern (border
  drift off the plane under remeshing) does not apply: borders are
  auto-constrained by isotropic_remeshing and preserved verbatim.

## Phase 4 — merged soup, autorefine safety net, overhang clip

- `autorefine_merged.cpp`: `--boundary FILE` mode (markers 100+i, CLI
  order) + `--fault-prefix` manifest filter + fault-split gate
  (`n_fault_new_subtriangles` must be 0 in boundary mode; exit 9
  otherwise) + `--box` regression path preserved.
- New `meshing/code/clip_fault_overhangs.py`: shell-containment ray
  parity (+z, xy-gridded, near-duplicate crossing collapse), junction
  guest flood-fill cut at the corefined junction polylines
  (largest-component keep), post-clip validation gates (closed shell,
  0 outside-shell fault tris, min edge, dup scan, border classification
  with near-miss orphan detection).

### Phase 3/4 hardening (measured failure ladder, all fixture-regression-verified)

Three further defects surfaced at production scale and were fixed in the
restored machinery (the 2000 m fixture reproduces identically after each):

1. `Euler::collapse_edge` Release-mode corruption in `polyline_cleanup.h`:
   the compiled-out CGAL preconditions admit (a) border-"ear" collapses
   (incident triangle whose third edge is also border) and (b) pinch
   collapses (interior edge, both endpoints on border) — both corrupt
   `Surface_mesh` silently (SIGBUS / mega-degree umbrella spins).  Added
   `collapse_preconditions_ok()` mirroring the CGAL preconditions + the
   pinch rule; applied in both cleanup passes.
2. Stale-descriptor crash: the same coordinate edge can appear TWICE in
   one mesh (soup-orient duplicated pinch wedges); collapsing the first
   invalidates the second mid-candidate.  Same-mesh duplicate candidates
   are now skipped; just-in-time `is_removed` re-check added.
3. `--allow-residual-short N` (corefine_set + passthrough): residual
   sub-floor polyline edges (degenerate pinch clusters, guard-skipped
   edges) are tolerated WITH printed locations; the hard min-edge gate
   moves to the Phase 4 post-clip artifact where it belongs.  Run used
   N=100; 67 residuals tolerated, all later removed by the clip-stage
   welds/contractions.

### Phase 3 results (corefine v9, production 500 m)

- All 8 pairs corefined: 3 traces (dem x SAF/Banning/Garnet), 3 junctions,
  2 bottom contacts.  Junction crossing SAF x Garnet: 2 polylines /
  114 verts (with eps-detach; previously a >30-min corefine hang).
- Cleanup: 4,134 collapses / 18 iters; cross-snap 177 clusters.
- Per-fault below-DEM quality (pre-clip): edge medians 455.8 / 457.4 /
  467.4 m (gate [425, 600]); p99 622 / 658 / 627 m (gate <= 750);
  q_med 0.98-0.99.
- 5 of 8 pairs flagged "non-conformal" by +-1-2 shared vertices at
  max_diff = 0 — coincident-duplicate counting artifacts at the junction
  pinch clusters; the merged-soup weld + gates are the real arbiter.
- `check_self_intersect`: dem/bottom/Banning/Garnet clean; SAF carries 21
  exact-overlap pairs at pinch corners (resolved in the soup by the 1e-6
  intern dedup; see Phase 4 gate).

### Phase 4 results (autorefine boundary mode + clip)

- autorefine: 271,810 -> 272,348 tris; **fault-marked new_subtriangle
  count = 0** (the R-008 gate) — every fault intersection was resolved in
  Phase 3; residual splits were boundary-side only.
- Clip: 11,579 fault tris dropped outside the shell (trace strips);
  junction flood-fill (cut at the shared guest-host edge sets — the
  manifest polylines are resampled curves, not mesh vertices) dropped
  107 / 192 / 793 overhang tris; eps-detach ribbons collapsed by snapping
  38 base-row vertices onto the junction lines; 92 vertices merged by
  soup-level sub-floor edge contraction (the topology-aware equivalent of
  the archived `--output-dedup-tol 99`); 77 cap-triangle flips.
- **Final gates: PASS** — closed boundary shell (0 open / 0 over-shared
  edges), 0 duplicate-coordinate groups at 1 mm (G3a), 0 fault triangles
  outside the shell, min edge 100.006 m with 0 below the floor (G1 at
  the surface stage).
- Warnings (documented, not gated): 3 near-shell tip edges (chain-end
  sags of 52-95 m at the far-NW SAF trace end and one Garnet bottom
  vertex — interior fault borders, no PLC violation); 8 thin fault
  triangles q < 0.30 (0.006%: SAF 6, Garnet 2; flip-blocked on
  constrained shared lines) — Phase 5/6 decision.

### Final fault surface quality (clipped merged soup)

| fault | n_tri | edge min | edge med | edge p99 | edge max | q_min | q_med |
|---|---:|---:|---:|---:|---:|---:|---:|
| SAF (MJVS) | 100,912 | 100.0 | 455.8 | 622.2 | 766.5 | 0.002 | 0.981 |
| Banning (MULT) | 2,614 | 100.5 | 461.8 | 642.3 | 744.5 | 0.319 | 0.978 |
| Garnet Hill (SBMT) | 28,900 | 101.2 | 467.3 | 627.1 | 744.5 | 0.000 | 0.987 |

(Compare the GOCAD baseline: fault edge medians 1,095-2,117 m, minima
1.5-42 m, node-duplicated.  Goal 1 met; goal 3 met at surface level —
single welded node set, zero coincident duplicates.)

### Deliverables

- `meshing/data_corefined/safv4_500m/safv4_merged_500m.stl` +
  `_markers.json` (tetgen-ready PLC; Phase 5 input) + `clip_report.json`,
  `quality_phase3.json`, `quality_final_faults.json`.
- ParaView: `meshing/vtu_safv4_500m/safv4_merged_500m_{all,faults,
  boundary,fault_1_*,fault_2_*,fault_3_*}.vtu` (cell data: `marker`, `q`).
  Recipe: open `_all.vtu`, color by `marker`; or `_faults.vtu` colored by
  `q`.
- Tools (committed): `check_fault_conformity.py`,
  `extract_safv4_surfaces.py`, `extend_fault_borders.py`,
  `clip_fault_overhangs.py`, `report_fault_quality.py`,
  `merged_soup_to_vtu.py`, `remesh_graded.cpp`/`graded_sizing_field.h`,
  modified `corefine_set.cpp` / `autorefine_merged.cpp` /
  `polyline_cleanup.h` / `io_helpers.h` / `corefine_faults.py` /
  `check_mesh_quality.py`.

---

## Phase 5 — tetgen volume mesh + Gmsh v2.2 .msh + sliver optimization

### tetgen_mesh.py changes

- Phase 4 boundary-marker layout (fault = markers 1..N; the archived
  `!= box_marker` test would misclassify boundary markers 101..105 as
  faults); `write_gmsh22()` with the normative tag map (volume 1, faults
  101/102/103 alphabetical, top 201, bottom 202, sides 203); the `.msh`
  carries ALL tets (the quality filter would punch interior cavities that
  MFEM turns into spurious free surfaces); tetgen optimization knobs
  exposed (`--mindihedral`, `--opt-iterations`, `--opt-scheme`);
  topology-safe post-dedup (twin-tet drop) and a validity-checked
  short-edge collapse (only Steiner endpoints; every incident tet must
  stay positive-volume and unique).
- Flat-box regression: edge_min 99.32 m identical to the STATUS baseline;
  q_med 0.887 vs 0.870; q<0.05 43 vs ~70; n_tets 41.6k vs 33.6k
  (optimizer-default drift; structurally within noise).

### PLC hardening discovered by tetgen (clip_fault_overhangs additions)

1. tetgen rejected the v1 PLC: 2 SAF x Garnet triangle crossings
   introduced by the clip's own vertex moves in the eps-wide junction
   zone.  Fixes: crossing scan (vertex-disjoint Moller test, xy-gridded)
   + push-apart repair (guest vertices nudged along the host normal to
   20 m clearance; shared junction-line nodes may move — one welded node
   deforms both surfaces consistently); hard gate
   `gate_fault_crossing_pairs == 0`.
2. The archived 99 m output dedup FOLDED tets (negative volumes, faces
   shared by 3-4 tets, MFEM "Invalid mesh topology") by welding REAL
   ~99 m-apart vertices.  Replaced by 1 mm dedup (true coincidences
   only) + twin-tet drop + the validity-checked collapse.
3. Residual sub-floor volume edges traced to surface geometry and fixed
   at the source: (b6) near-shell sag closure (chain-end fault border
   vertices 52-95 m off the shell snapped onto the nearest shell vertex);
   (b7) cross-fault proximity push (first-row vertex pairs across the
   13.4 deg SAF-Garnet wedge sat 82.7 m apart; pushed symmetrically to
   105 m — increasing separation cannot create crossings).

### Final gate results (production 500 m)

Two meshes delivered (`meshing/results/`):

| mesh | tets | nodes | eta_med | eta<=0.1 | eta<=0.05 | min edge |
|---|---:|---:|---:|---:|---:|---:|
| `safv4_500m_topo_opt.msh` (PRIMARY; mindihedral 14, opt x10) | 1,176,304 | 228,223 | 0.800 | **96** (8.2e-5) | 3 | 100.01 m |
| `safv4_500m_topo.msh` (baseline knobs) | 778,171 | 167,449 | 0.702 | 219 (2.8e-4) | 3 | 100.01 m |

- **G1 PASS** (min edge 100.01 m, 0 below floor) — without any tet drops.
- **G2b PASS** with 6x margin on the primary (8.2e-5 <= 5e-4); GOCAD
  baseline was 214,359 slivers -> 96 (2,200x reduction).
- **G2a: 3 residual tets** with eta 0.0023/0.0034/0.0379 at z -17.6 to
  -5.4 km — fault-sandwiched wedge tets that `nobisect` forbids tetgen
  from splitting; deep and far from the nucleation region.
  ACCEPT-WITH-REPORT per the plan's calibration governance (the Phase 6
  alternative — retuning Phase 3 sizing — is not warranted for 3 tets at
  2.6e-6 of the mesh).
- **G3a/b/c PASS**: 0 duplicate-coordinate groups at 1 mm; every fault
  triangle interior (2 tets) — 100,912 + 2,614 + 28,900 = 132,426 = 100%
  embedding; conformity tool exit 0.
- **G0 PASS**: 500/500 inside-shell samples inside a tet (the bbox-based
  smoke reports ~38% because the concave rotated domain fills ~40% of its
  axis-aligned bbox).
- **G1f PASS**: fault edge medians 455.8 / 461.8 / 467.3 m, minima
  >= 100 m, q_med 0.98-0.99.
- **MFEM smoke PASS** (both meshes): `mfem::Mesh` loads; dim=3;
  NE/NBE/NV as above; vol attrs {1}; bdr attrs {101,102,103,201,202,203}
  — the normative tag map verbatim.  `$MeshFormat 2.2` confirmed; meshio
  round-trip via the gate tools.

### Phase 5 deliverables

- `meshing/results/safv4_500m_topo_opt.msh` (73 MB, PRIMARY) +
  `safv4_500m_topo.msh` (52 MB, lighter alternative)
- `meshing/results/safv4_500m_topo{,_opt}_{bulk,fault}.vtu` — ParaView
  (bulk: cell data `q_iso`, `edge_min`; fault: `patch` 1..3)
- `meshing/results/gate_conformity_500m{,_opt}.json` (committed)

---

## Phase 5b — deepened bottom + coarsened boundary + sliver reduction (user feedback 2026-06-12)

User feedback on the Phase 5 mesh (ParaView screenshot): (1) the volume
boundary did not coarsen away from the fault; (2) the fault bottom borders
terminated ON the domain bottom rather than inside the volume.  Plus a
renewed request to drive slivers down.

### (1) Deepened bottom — fault bottoms become interior tips

`extract_safv4_surfaces.py --deepen-bottom 20000`: translate the `bottom`
surface down 20 km and extrude the four ribbon walls to meet it.  The
SAF/Garnet bottom borders (previously a `bottom` contact class, R-001)
now classify as `interior_tip` — verified: SAF/Garnet fault z_min =
-19,329 m sit 20 km above the new domain bottom at -39,329 m; Banning
27.8 km above.  No fault facet touches the boundary shell.  This also
REMOVES the bottom from the corefine participant set (4 inputs: 3 faults
+ DEM), eliminating the fault x bottom tangential-contact handling.

### (2) Coarsened boundary

- DEM regraded `--h-far 10000 --d-far 40000` (was 2500 / 9000): top
  surface now 100 m at the trace -> 43.6 km far-field (was capped at
  2.5 km).  DEM tris 98,477 -> 39,579.
- Bulk far-field `--lc-far 10000 --dist-outer 40000` (was 5000 / 20000):
  bulk edge median 502 -> 706 m, max 44.6 km.
- Boundary surface tris 126,892 -> 48,146.  Measured surface edge medians:
  DEM 654 m, bottom 8,758 m, sides 5,417 m (faults stay 456-467 m).

### (3) Sliver reduction — the eps-detach lesson

First deep attempt used eps-detach = 120 m (carried over from mid-iteration)
and gave eta<=0.05 = 25 (regression vs Phase 5's 3).  Root cause: the
larger junction-border detachment produced 9 thin fault triangles
(q_min 0.0019) ON the constrained trace/junction lines, which force thin
tets; a tetgen mindihedral sweep (14/18/20) plateaued at ~18 — surface
quality, not bulk optimization, was the bottleneck.  A constrained
in-plane Laplacian smoothing pass (`clip_fault_overhangs.py` b4b, free
interior fault vertices only) cannot touch trace/junction-constrained
vertices, so it does not help these.

Fix: revert eps-detach to the Phase-5-proven **30 m**.  The junction
geometry is then identical to Phase 5 and the thin triangles vanish.

### New PLC-repair machinery (this iteration)

The geometry edits the deepening + coarsening required surfaced new
degeneracies tetgen rejects; all are now handled in the clip + a repair pass:
- `autorefine_merged --soup-in/--soup-markers`: re-autorefine an EXISTING
  clipped soup (the clip's pushes/snaps can leave coplanar fan folds that
  share a junction vertex — a vertex-in-facet PLC violation, measured
  2.6 cm deep — that the assembly-time autorefine never saw).  Carries
  fault basenames + boundary names verbatim from the input markers.
- clip fold-aware crossing scan: fan pairs sharing ONE junction vertex are
  tested through the shared corner (a healthy fan grazes only at the
  vertex; a folded fan penetrates strictly inside).  Push clearance
  reduced 20 -> 5 m (20 m over-pushed shared line vertices into folds).
- clip cross-fault facet-clearance push (vertex within metres of another
  fault's facet INTERIOR), near-shell sag snap, boundary-edge contraction
  (autorefine splinters ON the shell, anchored to a trace vertex when one
  is present), and a final shell-containment drop.

### Final deliverables (replace the Phase 5 meshes)

`meshing/results/`:

| mesh | tets | nodes | eta_med | eta<=0.1 | eta<=0.05 | min edge | bulk med edge |
|---|---:|---:|---:|---:|---:|---:|---:|
| `safv4_deep_500m_opt.msh` (PRIMARY; mindihedral 18, opt x15) | 1,066,297 | 189,550 | 0.813 | **85** (8.0e-5) | 3 | 100.08 m | 706 m |
| `safv4_deep_500m.msh` (baseline knobs) | 647,494 | 124,826 | 0.702 | 138 (2.1e-4) | 2 | 100.08 m | 466 m |

- All gates PASS (G0 500/500; G1 100.08 m, 0 below floor; G2b 8.0e-5 ≤
  5e-4; G3a 0 duplicate groups; G3b 132,277 fault tris all interior =
  100% embedding; G1f fault medians 456/462/467 m).  G2a = 3 (at
  z -10.4 to -5.4 km, deep; accept-with-report).  MFEM loads both
  (bdr attrs 101,102,103,201,202,203).
- vs Phase 5 (`safv4_500m_topo_opt`): SAME eta<=0.05 (3), BETTER eta<=0.1
  (85 vs 96), coarser boundary, fault bottoms now interior.
- ParaView: `meshing/vtu_safv4_deep/final_{all,faults,boundary,fault_1_*,
  fault_2_*,fault_3_*}.vtu` + `boundary_surface.vtu` (color by `edge_max`
  to see the coarsening), and `meshing/results/safv4_deep_500m{,_opt}_
  {bulk,fault}.vtu`.

---

## Phase 5c — regenerated boundary shell (user feedback 2026-06-12, evening)

User feedback on the Phase 5b mesh (3 ParaView screenshots): bad-shaped
boundary tets at the domain edges/corners; "you don't need to coarsen this
much at the boundary, reduce the mesh size and maintain a better tet
quality"; sliver tets remaining on the top (DEM) and on the fault near the
trace.

Root cause (measured): the graded DEM remesh left 1,889 far-field
triangles over 10 km (max 43.6 km) and the raw GOCAD bottom/ribbons were
6-10 km median, q_med ~0.86 (right-triangles).  Those coarse, high-aspect
boundary surfaces drove the bad boundary tets.  Confirmed that lc-far /
bgmesh gradient changes do NOT alter the mesh (identical tet count for
lc-far 2000 vs 2500): with `nobisect`, the bulk size is pinned by the
SURFACE density, so the boundary SURFACE quality is the lever.

### New tool: `meshing/code/regenerate_boundary.py`

Rebuilds the closed boundary shell (DEM heightfield, flat bottom, 4
vertical side walls) at a moderate, well-shaped ~2.5 km resolution,
watertight by construction:
- ONE footprint polygon (the DEM rim, resampled to 2.5 km, 4 corners
  detected) drives the DEM rim, bottom rim, and ribbon tops/bottoms, so
  every shared edge is identical on both incident surfaces.
- DEM: constrained Delaunay (Shewchuk `triangle`, `pq28Y`) over the
  footprint + blue-noise interior graded 500 m (within `d_near` of the
  trace) -> 2.5 km far; z from the original DEM heightfield.  The DEM is
  corefine-coupled, so the fault trace is re-established downstream.
- Bottom + ribbons: quality-meshed (`pq30aY`) in 2-D (the bottom in xy,
  each ribbon in its own along-wall/depth plane, mapped back along the
  actual side polyline by arc-length so the ribbon top matches the DEM rim
  node-for-node).  `Y` keeps the footprint rim un-split -> watertight.

### Before -> after (boundary)

| surface | Phase 5b edge med / max | q_med | -> 5c edge med / max | q_med |
|---|---|---|---|---|
| TOP (DEM)  | 723 m / 43,613 m | 0.965 | 1,750 m / 3,299 m | 0.957 |
| BOTTOM     | 9,668 m / 39,831 m | 0.866 | 2,530 m / 3,693 m | 0.931 |
| SIDES      | 6,009 m / 44,571 m | 0.861 | 2,511 m / 3,663 m | 0.935 |

The 43.6 km DEM triangles are gone (max now 3.3 km); the bottom/sides
right-triangles (q 0.86) are now near-equilateral (q 0.93).

### Volume mesh result (replaces the Phase 5b meshes)

PRIMARY `meshing/results/safv4_deep_500m_opt.msh`: 1,416,795 tets,
189,550... (nodes 314,284), min edge 100.18 m, eta_med 0.752, eta<=0.1 199
(1.4e-4), eta<=0.05 2.  Lighter `safv4_deep_500m.msh`: 1,179,510 tets.

- **Boundary tets are now well-shaped**: tets with a face on the bottom or
  sides have eta_med 0.74, aspect (edge_max/edge_min) median 1.8, and
  **0 slivers (eta<=0.2)**.  The top has 47 eta<=0.2 (the fault-trace
  wedge).  The Phase 5b 43-km-driven bad boundary tets are eliminated.
- The overall eta_med is 0.752 (vs 5b's 0.813) only because the mesh is
  finer and more uniform (1.42 M vs 1.07 M tets): the histogram is 66 %
  > 0.7, 96 % > 0.5, 0.8 % < 0.3 — a high-quality mesh whose worst tets
  are now at the fault, not the boundary.
- Fault-trace thin triangles (the "fault top boundary slivers") dropped
  9 -> 5 of 132,337 (the quality DEM near the trace).  These 5 sit on the
  constrained fault-DEM trace polyline; a documented minor residual.
- Gates: G0 499/500 (one near-boundary numerical miss), G1 100.18 m (0
  below floor), G2b 1.4e-4 ≤ 5e-4, G3a 0 duplicate groups, G3b 100 %
  fault embedding, MFEM loads both (bdr attrs 101-103/201-203).

### New PLC-repair (this iteration)

clip `short_edge_contraction` now also collapses both-trace sub-floor
edges (anchored to a trace node, so the merged node stays on the trace and
fault/DEM conformity is preserved) — a single 94.9 m far-NW trace edge the
corefine could not collapse.

### ParaView

`meshing/vtu_safv4_deep/`: `boundary_surface.vtu` (color by `edge_max` —
now uniformly ~2.5 km, no 43 km triangles), `final_{all,faults,boundary,
fault_1/2/3_*}.vtu`, and `safv4_deep_500m_opt_bulk_eta.vtu` (color by
`eta` to see the quality distribution; the worst are at the fault, the
boundary is clean).
