# Implementation Plan: SAFv4 mesh-quality rebuild — uniform 500 m fault, sliver control, single-node fault trace

Anchor: `/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_preferred/`.
All relative paths below are relative to this anchor unless absolute.

Date: 2026-06-12.
Input under review: `meshing/model_SAFv4_mesh_2km_topo.inp` (SKUA-GOCAD export, 5.16 M nodes / 31.3 M tets, z in [-19329, +3116] m).
Quality log: `meshing/code/run_vtu_quality_model_SAFv4_2km_topo.log`.

## Overview

The SKUA-GOCAD-exported SAFv4 volume mesh fails every project quality gate and is structurally
unusable for the MFEM SEAS drivers because the volume mesh is **split (node-duplicated) along all
three fault surfaces and along the fault trace on the DEM**. This plan keeps SKUA-GOCAD as the
*geometry* source (the sealed SAFv4 model: 3 CFM6 fault patches trimmed to topography, DEM top,
4 ribbon side walls, flat bottom) but moves *all meshing* back in-repo: extract and weld the model
surfaces from the `.inp`, remesh the faults at a uniform 500 m target, and rebuild the volume mesh
with the proven CGAL-corefine -> autorefine -> tetgen pipeline (archived in git commit `577646e`),
extended from a flat box to a topography-bounded domain. Output is a Gmsh v2.2 `.msh` with physical
tags, plus VTUs, passing hard gates for all three goals.

### Diagnosis baseline (measured 2026-06-12, this session)

These numbers ground the plan and become the Phase 1 regression fixture. Measured directly on
`model_SAFv4_mesh_2km_topo.inp` with a coordinate-coincidence tolerance of 1e-3 m:

```
Goal-3 evidence (fault sides are node-duplicated, i.e. split mesh):
  SAF (MJVS)    : minus 2,564 tris / 1,527 nodes; plus 3,037 tris / 1,872 nodes
                  shared node IDs: 135; coincident coords: 1,527 (1,417 different-ID pairs)
                  minus tri centroids matching plus: 2,558 / 2,564
  Banning (MULT): minus 352 / 263; plus 449 / 360; shared IDs 51; 226 different-ID pairs
                  centroid match 352/352
  Garnet (SBMT) : minus 2,400 / 1,446; plus 2,643 / 1,608; shared IDs 206; 1,182 diff-ID pairs
                  centroid match 2,400/2,400
  Trace vs DEM  : 46..317 duplicated (coincident, different-ID) pairs per fault side
  Global        : 2,939 coincident groups among 6,490 fault-surface nodes (max multiplicity 4
                  -> at least one fault-fault or fault-fault-DEM junction exists)

Goal-2 evidence (slivers are a GOCAD volume-mesher artifact, NOT a fault problem):
  214,359 tets with Joe-Liu eta <= 0.1; eta_min = 1e-4
  z-histogram of sliver centroids: 58% in [-20 km, -15 km] (bottom layer);
  31% above z = -1 km (topo layer; 24.9% of all slivers touch a DEM node);
  only 0.4% within 1 km of a fault.

Goal-1 evidence (fault sizing):
  fault tri edge medians 1.10-2.12 km (target 500 m); minima 1.5-42 m (trace junction);
  all 1,106 sub-100 m-edge tets touch fault nodes (1,104) and DEM nodes (1,088)
  -> they sit at the fault-trace/DEM junction.

Fault border topology (measured 2026-06-12, coordinate coincidence at 1 mm on meshing/vtu/*.vtu):
  SAF (MJVS) x bottom : 53 coincident nodes (fault z_min = -19,329.4 m = the bottom plane)
  Garnet (SBMT) x bottom : 5 coincident nodes
  fault-fault junctions : MJVS x MULT 9, MJVS x SBMT 21, MULT x SBMT 68 coincident nodes
  no fault touches a ribbon; ribbons planar to < 1 m max plane-fit deviation;
  bottom exactly planar (0.000 m); boundary shell rims coincide at 1 mm
  (DEM x ribbons 47-142 nodes, ribbons x bottom 24-69, ribbon-ribbon corners 3-4)
  -> the shell welds watertight at tol 1e-3; faults are trimmed at the BOTTOM and at
  EACH OTHER (sealed model), not only at the DEM. Border handling must cover all classes.

Systematic, not one-off: the earlier GOCAD export SAFtopo_test1km.inp (meshing/code/run_quality.log,
2026-06-04, 47 M tets, 4 faults) shows the identical signature (224,163 slivers, min edge 3.62 m,
plus/minus tri-count mismatch). Fixing inside GOCAD is therefore not pursued: the split fault is
structural to SKUA's sealed-model export, and the bottom/topo sliver layers are its volume mesher's
behavior with no exposed quality control. The `.inp` would also need new conversion code anyway,
because MFEM consumes Gmsh v2.2 `.msh` only (see miniapps/seas/CLAUDE.md, "Known limitation").

### Why the tetgen pipeline and not Gmsh

- Gmsh multi-fault `Surface{} In Volume{}` fails with a tetgen PLC error; HXT silently drops
  surfaces ("Unknown surface 9-12") — documented in `document/PLAN_revert_to_gmsh.md` (Overview).
- The combined-STL Gmsh variant also failed PLC validation — `document/STATUS_tetgen_pipeline.md`
  ("Pre-existing Documents": "Failed at tetgen PLC validation; superseded").
- The direct CGAL-autorefine -> tetgen pipeline is the only proven multi-fault volume mesher in
  this project: 100% fault embedding, 0 cavities, 6-fault flat-box fixture
  (`document/STATUS_tetgen_pipeline.md`, "Current Output Status").
- The alternative project's Gmsh recipe (`../project_7.0_alternative/meshing/code/run_z0cut_meshing.py`)
  is proven only for ONE fault embedded in an OCC box with a FLAT top; SAFv4 needs 3 faults under a
  discrete DEM. Not a template here.

### Pipeline (target state)

```
meshing/model_SAFv4_mesh_2km_topo.inp        (SKUA-GOCAD; geometry source, frozen)
        |
        v
[P2] meshing/code/extract_safv4_surfaces.py   weld duplicate nodes (1e-3 m), emit canonical STLs:
        |                                     3 faults (minus side), dem, bottom, 4 ribbons + manifest
        v
[P3] meshing/code/extend_fault_borders.py     extrude fault borders past their contacted entities
        |                                     (trace +z above DEM; bottom -z; junctions in-plane)
[P3] code_preprocess/corefine_faults.py       remesh faults @ 500 m; remesh DEM (graded, rim
        |  + corefine_cgal binaries           protected); corefine fault-fault and fault-DEM;
        |                                     polyline cleanup; min-edge enforcement
        v
[P4] code_preprocess/corefine_cgal/autorefine_merged
        |                                     merge dem+bottom+ribbons+faults soup; safety-net
        |                                     autorefine; clip fault overhangs (outside shell /
        |                                     past host faults); markers JSON
        v
[P5] code_preprocess/tetgen_mesh.py           tetgen (plc, nobisect, quality, graded bgmesh);
        |                                     output dedup; quality gates
        v
     meshing/results/safv4_<R>m.msh (Gmsh v2.2, physical tags) + bulk/fault VTUs
        |
        v
[P1] meshing/code/check_fault_conformity.py + meshing/code/check_mesh_quality.py   (gates)
```

## Constraints

### Interface constraints
- `meshing/model_SAFv4_mesh_2km_topo.inp` is FROZEN — geometry source only, never edited.
  If the user updates the SAFv4 model in SKUA-GOCAD, a re-export of this file re-feeds the whole
  pipeline; everything downstream must be regenerable by script.
- `meshing/code/check_mesh_quality.py` metric definitions (Joe-Liu eta, triangle q, histogram
  bins, Q1/Q2 gate semantics) must NOT change — quality numbers must stay comparable with the
  existing logs. Only additive input-format support is allowed.
- `meshing/code/gocad_inp_to_vtu.py` untouched.
- Final mesh format: Gmsh v2.2 ASCII `.msh` (MFEM's reader is v2.2-only; a v4 file aborts at
  mesh/mesh_readers.cpp:1628 — miniapps/seas/CLAUDE.md "Known limitation").
- Restored pipeline files go back to their original archived paths (`code_preprocess/...`),
  so `git log --follow` keeps their history. Do NOT place new artifacts in any directory named
  `archive/` (globally gitignored by `miniapps/seas/.gitignore`).

### Dependency constraints
- CGAL 6.1.1 via conda env `cgal-61` (C++ builds: `cmake` against
  `code_preprocess/corefine_cgal/CMakeLists.txt`).
- Python steps in conda env `pythonenv`: meshio, numpy, scipy (cKDTree), tetgen 0.8.4 (verified
  installed), gmsh 4.15.0 (only if needed for `.msh` inspection — the pipeline itself does not
  call gmsh).
- The archived pipeline is recovered from commit `577646e` ("safs: archive project_7.0_preferred
  (CGAL-corefine multi-fault pipeline)") — restore via
  ```
  git checkout 577646e -- miniapps/seas/safs/project_7.0_preferred/code_preprocess \
      miniapps/seas/safs/project_7.0_preferred/raw_data \
      miniapps/seas/safs/project_7.0_preferred/data_corefined
  ```
  (`raw_data/` holds the CFM `.ts` inputs the flat-box regression fixture is regenerated from;
  the commit message excludes the corefined STLs as "regenerable", so the fixture MUST be
  regenerated — `data_corefined/` restores only the manifests.)
- Large artifacts (`*.stl`, `*.msh`, `*.vtu`, `*.inp`) are NOT committed (consistent with the
  existing `.gitignore` handling of 500 m meshes); manifests, scripts, and logs are committed.

### Convention constraints
- Per repo CLAUDE.md: if a phase fails, stop and report — do not silently fall back (this matters
  for Phase 5's PLC risk and Phase 6's sliver bar).
- Plain-ASCII math in all docs/comments (no LaTeX).
- Keep the alphabetical-basename ordering convention for stable per-fault tag assignment
  (established in `document/PLAN_cgal_corefine_multifault.md` Phase 3).

### Numerical constraints (the acceptance contract)
Metric definitions (identical to `check_mesh_quality.py`):

```
tet:      Joe-Liu eta = 12 * (3V)^(2/3) / sum_{i=1..6}(L_i^2),  regular tet -> 1
triangle: q = 4*sqrt(3)*A / sum_{i=1..3}(L_i^2),                equilateral -> 1
```

Hard gates on the final volume mesh (Phase 5):
- G1 (min edge):       min tet edge >= 100 m. 0 violations.
- G2a (hard slivers):  0 tets with eta <= 0.05.
- G2b (soft slivers):  count(eta <= 0.1) / n_tets <= 5e-4 (goal: 0; Phase 6 escalation if missed).
- G3a (no duplicates): 0 coincident-coordinate node pairs with distinct IDs at 1e-3 m tolerance,
                       mesh-wide (faults, trace, everywhere).
- G3b (embedding):     every fault triangle is the shared face of exactly 2 tets.
- G3c (trace):         every fault-trace node IS a DEM-surface node (same node ID).
                       Operational definition on .msh (no trace marker exists there): compute
                       fault border edges (exactly 1 incident fault-tagged triangle, per fault
                       tag); classify each — shared with a top(201) triangle edge -> trace; with
                       bottom(202)/sides(203) -> termination; shared with another fault tag ->
                       junction; else tip. G3c passes iff every border edge is classified (no
                       orphans) AND trace-edge endpoints appear in top-surface triangles with
                       the SAME node IDs. Report per-class counts per fault.
- G1f (fault sizing):  per fault patch: edge median in [425, 600] m; p99 <= 750 m; min >= 100 m.
- G0 (coverage):       500/500 random domain-interior points inside a tet (no cavities).

Surface-stage gates (Phase 3, before volume meshing):
- every surface edge >= 100 m after cleanup;
- fault tri q_min >= 0.30; DEM tri q_min >= 0.20;
- corefined pairs conformal at max_diff_m = 0 (manifest semantics identical to the archived
  pipeline's `data_corefined/manifest.json`).

Scale constraint: production target on-fault h = 500 m. Estimated fault triangle count
~70-90 k (scaling current counts by (median/500)^2: SAF 2,564 x (2117/500)^2 ~ 46 k, Garnet
2,400 x (1277/500)^2 ~ 15.7 k, Banning ~ 1.7 k, x2 safety for trace grading), volume mesh
~5-20 M tets depending on lc_far — same order as the current 31 M-tet GOCAD mesh, known to be
processable locally (this session) and within Frontera job sizes. 250 m on-fault is known to OOM
the driver — do not target below 500 m without explicit user approval.

---

## Phase 0: restore the archived pipeline and rebuild its binaries

### Goal
The CGAL-corefine/autorefine/tetgen toolchain that produced the proven 2000 m flat-box fixture
exists again on disk and compiles, unmodified.

### Files to Create
- None new; restored from git:
  ```
  git checkout 577646e -- miniapps/seas/safs/project_7.0_preferred/code_preprocess \
      miniapps/seas/safs/project_7.0_preferred/raw_data \
      miniapps/seas/safs/project_7.0_preferred/data_corefined
  ```
  (`raw_data/` = the CFM `.ts` regression inputs; `data_corefined/` = archived manifests only —
  the corefined STLs were excluded from the archive and are regenerated below.)
  The code_preprocess pathspec restores (29 files; document/ excluded):
  - `code_preprocess/corefine_cgal/{CMakeLists.txt, corefine_pair.cpp, corefine_set.cpp,
     autorefine_merged.cpp, check_self_intersect.cpp, mesh_volume.cpp, distance_sizing_field.h,
     intersection_graph.h, io_helpers.h, polyline_cleanup.h, polyline_resample.h, quality_repair.h}`
  - `code_preprocess/{corefine_faults.py, tetgen_mesh.py, ts_to_stl.py, clean_freesurface_mesh.py,
     generate_multifault_geo.py, merge_corefined_faults.py, msh_to_vtu.py, medit_to_vtu.py,
     stl_to_vtu.py, stitch_combined_stl.py, check_msh_quality.py, locate_slivers.py,
     sensitivity_sweep.py}`

### Files to Modify
- None in this phase (build only).

### Detailed Requirements
1. Restore exactly the paths above; verify `git status` shows them as new files at their original
   locations.
2. Build C++ tools:
   ```
   conda activate cgal-61
   cd code_preprocess/corefine_cgal
   cmake -S . -B build && cmake --build build -j \
       --target corefine_pair corefine_set autorefine_merged check_self_intersect
   ```
   (`mesh_volume` may be built too but is fallback-only; do not wire it into the new pipeline.)
3. Verify the Python entry points run: `python code_preprocess/corefine_faults.py --help`,
   `python code_preprocess/tetgen_mesh.py --help` (env `pythonenv`).
3b. Regenerate the flat-box regression fixture:
   `python code_preprocess/corefine_faults.py --res 2000`; record its output manifest as THE
   regression reference for Phases 4-5, and diff it against the restored archived manifest —
   any drift (CGAL env changes since 2026-05) is reported here, before any new code exists.
4. Restored-tool CLI reference (verified against the archived source in `577646e`):
   - `corefine_faults.py --in-dir D --out-dir D --res N --mesh-edge-size M --min-edge M
      --polyline-spacing M --features-angle-bound DEG --cgal-bin P --workdir D
      [--keep-intermediate] [--verbose]`
   - `tetgen_mesh.py --res N --merged-stl P --markers-json P --out-base P --lc-near M --lc-far M
      --dist-inner M --dist-outer M --bg-spacing M --minratio X --bulk-min-edge M --bulk-q-min Q
      --vertex-merge-tol M --output-dedup-tol M [--mmg-cleanup] [--verbose]`
   - `autorefine_merged MANIFEST OUT_STL OUT_MARKERS [--pad-xy M] [--pad-top M] [--pad-bottom M]
      [--box-edge-size M] [--box-marker N] [--verbose]`

### Edge Cases to Handle
- CGAL version drift since 2026-05: if the cgal-61 env was upgraded and the build breaks, STOP and
  report the compiler errors (do not patch CGAL calls without approval).

### Acceptance Criteria
- [ ] All 4 C++ targets build with zero errors in `cgal-61`.
- [ ] `corefine_faults.py --help` and `tetgen_mesh.py --help` exit 0 in `pythonenv`.
- [ ] `python -c "import tetgen, meshio, scipy"` exits 0 in `pythonenv`.
- [ ] `ls raw_data/*.ts` non-empty (regression inputs restored).
- [ ] 2000 m fixture regenerated; manifest diff vs the archived manifest committed to the run
      log (drift reported, not silently absorbed).

### Dependencies
- Depends on: nothing.
- Required by: Phases 3, 4, 5.

---

## Phase 1: conformity + sliver diagnostics tool (the gate harness)

### Goal
A committed, reusable checker proves/refutes goals 1-3 on any mesh (`.inp`, `.msh`, `.vtu`), locks
today's diagnosis as a regression fixture, and later serves as the acceptance gate for the rebuilt
mesh.

### Files to Create
- `meshing/code/check_fault_conformity.py`

### Files to Modify
- `meshing/code/check_mesh_quality.py` — additive only: accept `.msh`/`.vtu` input. For `.msh`,
  surface triangles come from `gmsh:physical` cell data (one "surface" per physical tag) instead
  of `*SURFACE` ELSETs; bulk tets and all metric kernels unchanged. New optional CLI flag
  `--fault-tags 101,102,103` to label which physical tags are faults in the report.

### Detailed Requirements
1. `check_fault_conformity.py` module API (all NumPy-vectorized; chunk tet loops at 2e6 to bound
   memory, same pattern as `check_mesh_quality.py`):
   ```python
   def load_surfaces(path: Path, fault_tags: list[int] | None = None) -> SurfaceModel
       # SurfaceModel: dataclass with points (N,3) float64, tets (M,4) int64,
       #   surfaces: dict[name -> (K,3) int64 triangle node-IDs]
       # .inp: reuse the FACE_NODES/_cell_set_global_tet_ids/_extract_surface_triangles helpers
       #   imported from check_mesh_quality (do not re-implement);
       # .msh/.vtu: triangles grouped by gmsh:physical tag.

   def duplicate_node_groups(pts: np.ndarray, ids: np.ndarray | None = None,
                             tol: float = 1e-3) -> list[list[int]]
       # quantized-key grouping: key = (round(x/tol), round(y/tol), round(z/tol));
       # returns groups with >= 2 distinct node IDs. ids=None means scan all points.

   def pair_conformity(model: SurfaceModel, name_a: str, name_b: str,
                       tol: float = 1e-3) -> dict
       # returns {n_tri_a, n_tri_b, n_nodes_a, n_nodes_b, shared_node_ids,
       #          coincident_coord_pairs, duplicated_pairs, centroid_matches}
       # (exactly the quantities in the Diagnosis baseline table above)

   def fault_embedding(model: SurfaceModel, fault_name: str) -> dict
       # builds the key set of FAULT triangles only (sorted node triples, ~1e4-1e5 entries),
       # then streams tets in 2e6 chunks: form each chunk's 4 face triples (vectorized np.sort
       # on the (4*chunk, 3) index array), pack each triple into a single int64 key, and count
       # matches against the fault key set (np.isin); accumulate per-fault-triangle incidence.
       # Memory O(n_fault_tri + chunk), never O(n_tet) — a global face->tet map over all
       # 4*31.3e6 faces would exceed 12 GB and is forbidden.
       # returns {n_tri, n_interior (==2 incident tets), n_boundary (==1), n_orphan (==0)}

   def sliver_report(model: SurfaceModel, eta_floor: float = 0.1) -> dict
       # eta per tet (chunked); sliver centroid z-histogram with the bin edges used in the
       # baseline ([-20000,-15000,-10000,-5000,-2500,-1000,0,+4000] clipped to mesh z-range),
       # counts touching fault / DEM node sets (>=1 and >=3 nodes), and centroid
       # distance-to-fault-node percentiles (p5,p25,p50,p75,p95) via scipy cKDTree.
   ```
2. CLI:
   ```
   python check_fault_conformity.py MESH [--fault-tags 101,102,103] [--dem-name NAME|--dem-tag N]
       [--pairs auto|none] [--tol 1e-3] [--json out.json]
   ```
   `--pairs auto` (default for `.inp`): pair every `*_minus*` surface with its `*_plus*` sibling.
   Exit code 0 iff: zero duplicate groups AND (if tets present) every fault tri is interior
   (n_interior == n_tri). Print a markdown report mirroring the baseline table.
3. Regression fixture: run the tool on `model_SAFv4_mesh_2km_topo.inp` and commit the JSON output
   as `meshing/code/baseline_safv4_2km_topo_conformity.json`. The numbers MUST reproduce the
   Diagnosis baseline above (SAF duplicated pairs 1,417; Banning 226; Garnet 1,182; global groups
   2,939; slivers 214,359; sub-100 m tets 1,106).
4. Unit fixture: a tiny self-contained test (`python check_fault_conformity.py --self-test`) that
   builds in-memory (a) two tets sharing a fault triangle with shared nodes -> PASS, (b) the same
   two tets with the shared triangle's 3 nodes duplicated -> detects 3 duplicate groups and an
   orphan fault face -> FAIL exit code. No external files.

### Interfaces
- Consumed by Phases 2-5 acceptance runs; `--json` output is the machine-readable gate record.

### Edge Cases to Handle
- Multiplicity > 2 coincident groups (junction nodes, observed max 4): report group size; any
  group with >= 2 distinct IDs counts as a violation.
- `.msh` with no tets (surface-only intermediate): skip embedding/sliver sections, still run
  duplicate scan.
- Surfaces whose plus set has extra triangles (observed: SAF plus 3,037 vs minus 2,564): report
  `extra_a`, `extra_b` counts; conformity is judged on the matched subset.

### Acceptance Criteria
- [ ] `--self-test` passes both sub-cases.
- [ ] Baseline JSON committed and bit-stable across two runs.
- [ ] Tool exits nonzero on `model_SAFv4_mesh_2km_topo.inp` (it is a known-bad mesh).
- [ ] `check_mesh_quality.py` on a `.msh` produces the same bulk metrics as on the equivalent
      `.vtu` (cross-check on an existing v2.2 fixture already in the repo, e.g. a tpv205/tpv102
      `.msh` under miniapps/seas/, converted to `.vtu` via meshio — NOT a Phase 5 artifact;
      Phase 1 must be acceptable before Phase 5 exists).

### Dependencies
- Depends on: nothing (parallel with Phase 0).
- Required by: Phases 2, 3, 4, 5 acceptance.

---

## Phase 2: extract + weld the SAFv4 surfaces from the GOCAD .inp

### Goal
Canonical, deduplicated surface meshes for the SAFv4 geometry exist as STLs + a manifest, with the
split fault collapsed to a single triangulation per fault and the trace welded to the DEM.

### Files to Create
- `meshing/code/extract_safv4_surfaces.py`
- Output artifacts (gitignored): `meshing/data_extracted/safv4_2km_topo/{fault_<basename>.stl x3,
  dem.stl, bottom.stl, ribbon_NE30.stl, ribbon_NE120.stl, ribbon_NW60.stl, ribbon_NW150.stl,
  manifest.json}`

### Files to Modify
- None.

### Detailed Requirements
1. CLI:
   ```
   python extract_safv4_surfaces.py INP --out-dir meshing/data_extracted/<stem>
       [--tol 1e-3] [--fault-side minus]
   ```
2. Read the `.inp` once (meshio); extract every `*SURFACE` triangle set with the same
   FACE_NODES convention as `check_mesh_quality.py` (import, don't copy).
3. **Canonical fault triangulation = the minus side** (`--fault-side minus`). Rationale (measured):
   minus is a subset of plus up to <= 6 triangles (SAF 2,558/2,564 matched; Banning 352/352;
   Garnet 2,400/2,400); the plus extras (479/97/243) are tip/edge wrap-around faces, not fault
   area. Requirement: compute the matched set via centroid keys at `tol`; report
   (a) minus triangles unmatched in plus, (b) plus extras. If unmatched-minus > 1% of the patch,
   ABORT with a report (the subset assumption broke; user decides).
4. **Global weld**: quantized-key dedup of ALL mesh nodes at `tol = 1e-3 m` (same key scheme as
   Phase 1). Representative = lowest node ID in each group. Re-index all extracted surface
   triangles through the weld map. Drop triangles that become degenerate (two equal indices) —
   count and report; expected ~0 because plus/minus duplicates never appear in the same triangle.
5. Per-surface STL output, binary STL via meshio, points pruned to used set per file.
6. `manifest.json` (committed): per surface — basename, n_tri, n_vertices, bbox, edge-length
   (min/med/max), and for each fault a **border classification**: every border edge assigned to
   exactly one of {dem (trace), bottom, ribbon_<name>, fault_<basename> (junction),
   interior_tip}, by post-weld shared-node test against each candidate surface (both endpoints
   shared -> contact). Emit ordered chains per class. For junction chains, record host/guest:
   the fault whose border lies ON the other's interior is the guest. Also record the weld stats,
   split by class: `n_groups_total` (informational; expected > 2,939 because boundary rims are
   also coincident — measured: DEM x ribbons 47-142 nodes, ribbons x bottom 24-69), and
   `n_groups_fault_subset` (groups containing >= 1 fault-surface node), which MUST equal the
   baseline 2,939 at tol = 1e-3 m. Plus the minus/plus match counts.
7. Determinism: byte-stable manifest across re-runs (sorted keys, fixed float formatting,
   no timestamps inside; meshio binary STL is deterministic for identical input).

### Interfaces
- `manifest.json` schema feeds Phases 3 and 4 (mirrors the archived `data_corefined/manifest.json`
  pattern: `meshes[]`, later `pairs[]` added by corefine).

### Edge Cases to Handle
- Junction nodes with multiplicity 4 (fault-fault-DEM): weld all to one representative.
- Trace chain endpoints that touch a ribbon (fault reaching the domain side): the chain ends
  there; record the endpoint as `on_boundary: true`.
- Ribbon/bottom/DEM rim: after the weld, every rim edge must be shared by exactly 2 boundary
  surfaces (DEM-ribbon, ribbon-ribbon, ribbon-bottom). Verify and report; if a gap exists at
  `tol`, ABORT (GOCAD export inconsistency — user decides).

### Acceptance Criteria
- [ ] Phase 1 tool on the welded surface set: 0 duplicate groups among all used nodes.
- [ ] Fault STL triangle counts equal the baseline minus counts (2,564 / 352 / 2,400).
- [ ] Border classification reproduces the measured contacts (1 mm tol): SAF-bottom 53 nodes,
      Garnet-bottom 5, MJVS-MULT 9, MJVS-SBMT 21, MULT-SBMT 68, zero fault-ribbon contact.
- [ ] Weld stats: `n_groups_fault_subset` == 2,939; `n_groups_total` > 2,939.
- [ ] Every fault trace node ID is in the DEM node-ID set (G3c at surface level).
- [ ] Boundary watertightness: every non-fault surface edge shared by exactly 2 boundary
      triangles (closed shell of dem+ribbons+bottom).
- [ ] Manifest byte-stable across two runs.

### Dependencies
- Depends on: Phase 1 (tool reuse).
- Required by: Phases 3, 4.

---

## Phase 3: fault remesh at 500 m, DEM graded remesh, corefine, cleanup

### Goal
Fault patches are uniform-500 m triangulations extended transversally past EVERY contacted
entity (DEM above, bottom below, host faults at junctions); the DEM is remeshed graded (fine near
traces, coarse far); all fault-fault and fault-DEM intersections are corefined conformal with
cleaned polylines; every surface edge >= 100 m.

### Files to Create
- `meshing/code/extend_fault_borders.py`
- `code_preprocess/corefine_cgal/remesh_graded.cpp` — graded isotropic remeshing of one surface:
  ```
  remesh_graded IN.stl OUT.stl --trace-polylines manifest.json --h-near 500 --h-far 2500
      --d-near 1000 --d-far 9000 [--protect-border] [--verbose]
  ```
  Implements `h(d) = h_near + (h_far - h_near) * clamp((d - d_near)/(d_far - d_near), 0, 1)` as a
  model of PMPSizingField for the `PMP::isotropic_remeshing` sizing-field overload (CGAL >= 5.6;
  present in 6.1), with distance d to the trace polyline vertices via a CGAL KD/AABB query.
  Border edges constrained (`edge_is_constrained_map`, `protect_constraints = true`) and border
  vertices fixed (`vertex_is_constrained_map`) — the rim vertex set must be byte-identical
  pre/post. Add the target to `corefine_cgal/CMakeLists.txt`.
- Working artifacts (gitignored): `meshing/data_corefined/safv4_<R>m/...` + updated manifest.

### Files to Modify
- `code_preprocess/corefine_faults.py` — generalize input discovery: accept an explicit
  `--manifest` (Phase 2 output) listing fault STLs AND the DEM as a corefine participant, instead
  of the original `--in-dir` raw-CFM `.ts` scan. The corefine/remesh core (corefine_set
  invocation, manifest pairs/polylines emission) is reused as-is.
- `code_preprocess/corefine_cgal/corefine_set.cpp` — only if needed: the DEM's graded remesh is
  owned by the NEW `remesh_graded` tool (runs before corefine), so corefine_set must be able to
  take the DEM as a corefine participant WITHOUT re-remeshing it (and faults at 500 m uniform).
  If the archived CLI already supports per-mesh edge sizes / a no-remesh flag, no change;
  otherwise add `--mesh-edge-size-per-input a.stl=500,dem.stl=keep`. Investigate before coding;
  report which path was taken.

### Detailed Requirements
1. **Extension** (`extend_fault_borders.py`), runs BEFORE remeshing. Input: fault STL + Phase 2
   manifest border classification + dem STL. For each border chain, extrude by class (quad strip
   between consecutive border vertices, each quad split into 2 triangles, strip base welded to
   the fault border):
   - `dem` (trace): extrude `+z` by `--h-ext` (default 1500 m). Verification (hard): sample the
     strip top polyline at 100 m spacing; every sample must be STRICTLY above the DEM surface
     (point-location KDTree on DEM triangles + barycentric z). If violated, raise `--h-ext`
     automatically to `max_violation + 500` and re-emit, then report.
   - `bottom`: extrude `-z` by `--h-ext-bot` (default 500 m) below `z_bottom = -19329.4 m`
     (the bottom is exactly planar, max plane-fit deviation 0.000 m, measured).
   - `fault junction` (guest side only, from the manifest host/guest assignment): extrude along
     the local fault-plane outward direction (average of border-edge in-plane normals) by
     `--h-ext-junc` (default 750 m), so the extended patch crosses the host fault transversally.
     The HOST fault is NOT extended at the shared junction.
   - `interior_tip`: no extension.
   - `ribbon`: none observed in SAFv4 — assert; abort with a report if the Phase 2
     classification produced one.
   - Output: `fault_<basename>_ext.stl`.
2. **Remeshing + corefine** via the restored pipeline on the set
   {3 extended faults, dem} (+ ribbons/bottom NOT remeshed and NOT corefined here — they stay
   exactly as extracted; their quality is acceptable: q_min 0.30-0.75, min edge 829-1448 m):
   - Faults: CGAL isotropic_remeshing, `mesh_edge_size = 500`, `min_edge = 100`,
     `polyline_spacing = 250` (resampling of intersection polylines; half the target edge),
     `features_angle_bound = 60` (archived defaults otherwise).
   - DEM: graded remesh via the NEW `remesh_graded` tool (see Files to Create) with
     `--h-near 500 --h-far 2500 --d-near 1000 --d-far 9000`, distance measured to the
     (pre-corefine, Phase 2 manifest) trace polyline vertices; rim border edges and vertices
     constrained so the welded ribbon contact is untouched (rim byte-identical pre/post).
     DEM rationale: current DEM median 2.47 km next to a 500 m fault forces the trace-junction
     slivers; grading bounds the neighbor-size ratio near the trace to ~1.
   - Corefine pairs: all pairs with overlapping bboxes among {faults} x {faults} and
     {faults} x {dem}. The three junction pairs (MJVS-MULT, MJVS-SBMT, MULT-SBMT — all measured
     in contact today) are transversal crossings AFTER the junction extension; the corefined
     intersection polyline replaces the old trimmed junction, exactly as the trace replaces the
     welded fault-DEM contact.
   - Polyline cleanup + cross-polyline snap + min-edge enforcement: archived
     `polyline_cleanup.h` machinery, unchanged semantics (this is what achieved
     `max_diff_m = 0` and >= 96 m edges on the 6-fault fixture).
3. **Trace replaces the weld**: after corefine, the fault/DEM intersection polyline is the new
   trace (it will be within ~remesh-displacement of the old welded trace). The Phase 2 welded
   trace is only used for the DEM sizing field and the extension step.
4. Order of operations (normative): extract(P2) -> extend(3.1) -> remesh_graded on the DEM ->
   fault remesh + corefine + cleanup(3.2) -> [Phase 4: merge/autorefine/clip].
5. Outputs: corefined STLs + `manifest.json` with `pairs[]` (conformality + polylines), same
   schema as the archived `data_corefined/manifest.json`.

### Interfaces
- Output manifest consumed by Phase 4's `autorefine_merged`.

### Edge Cases to Handle
- **Near-tangent fault-DEM intersection** (gently-dipping fault under gentle topo): produces
  splinter triangles along the trace. The cleanup pass must collapse trace segments < 100 m
  (polyline_cleanup) — acceptance gate below enforces it. If cleanup cannot reach 100 m without
  breaking conformality, STOP and report the offending trace window (coordinates + lengths).
- **Fault tip inside the domain** (tipline not reaching DEM or boundary): allowed; the tip border
  is a free border — no extension, no constraint.
- **Trace chain ending on a ribbon**: extension strip would poke through the ribbon. Clamp strip
  x,y at the ribbon plane (project the last strip column onto the ribbon), flagged from the
  Phase 2 `on_boundary` marker.
- **Remesh displacement off the DEM**: isotropic remeshing moves fault vertices slightly off the
  old trace; this is fine — conformity is re-established by corefine in the same phase, not by
  the weld.

### Acceptance Criteria
- [ ] Per-fault (post-corefine, pre-clip, measured on below-DEM triangles only): edge median in
      [425, 600] m, p99 <= 750 m, min >= 100 m, tri q_min >= 0.30.
- [ ] DEM: min edge >= 100 m, q_min >= 0.20, edge median within [400, 700] m inside 1 km of the
      trace and >= 1500 m beyond 9 km.
- [ ] All corefined pairs: `max_diff_m = 0` in the manifest (bit-conformal shared polylines).
- [ ] `check_self_intersect` clean on every output STL.
- [ ] Extension strips strictly past their contacted entity (above the DEM / below the bottom
      plane / across the host fault at junctions; verification step log committed).
- [ ] DEM rim vertex set byte-identical pre/post remesh_graded.

### Dependencies
- Depends on: Phases 0, 2.
- Required by: Phase 4.

---

## Phase 4: merged soup, autorefine, overhang clip, closed PLC

### Goal
One non-manifold conforming triangle soup (boundary shell dem+ribbons+bottom, faults embedded,
fault overhangs — above the DEM, below the bottom, past host faults — removed), with
per-triangle provenance markers, ready for tetgen.

### Files to Create
- `meshing/code/clip_fault_overhangs.py`

### Files to Modify
- `code_preprocess/corefine_cgal/autorefine_merged.cpp` — replace the generated flat box with
  user-supplied boundary surfaces. New CLI (keep the old box path behind `--box` for the flat-box
  fixture regression):
  ```
  autorefine_merged MANIFEST OUT_STL OUT_MARKERS
      --boundary dem.stl --boundary bottom.stl --boundary ribbon_NE30.stl ... (repeatable)
      [--box | (default: boundary mode)] [--verbose]
  ```
  Each `--boundary` surface gets its own marker (recorded in OUT_MARKERS, named); fault markers
  remain per-fault from the manifest, exactly as the archived marker-propagating visitor does.

### Detailed Requirements
1. Soup assembly order: boundary surfaces first (markers 100+i in CLI order), then faults in
   alphabetical-basename order (markers 1..N) — matches the archived convention.
2. `PMP::autorefine_triangle_soup` with marker propagation is a SAFETY NET: all fault x fault
   and fault x DEM (incl. extension strips) intersections are already corefined + cleaned in
   Phase 3. Expected residual work here: strip x ribbon clamps and exact-duplicate removal only
   (the mechanism proven to fix "1 duplicate triangle + 8 segment-facet intersections -> 0" on
   the fixture). Gate: autorefine's resolved-intersection count involving any FAULT marker must
   be 0; if > 0, STOP and report which Phase 3 pair leaked (a Phase 3 bug — not something to
   absorb here).
3. **Clip** (`clip_fault_overhangs.py`), applied per fault patch AFTER autorefine:
   - (a) shell containment: drop every fault-marked triangle whose centroid is OUTSIDE the
     closed boundary shell. Test = ray parity (+z ray from centroid, count crossings against
     boundary-marked triangles; outside if even). The DEM-heightfield shortcut
     (centroid z > z_DEM(x, y)) is valid for the top only and misses below-bottom overhangs —
     do NOT use it as the sole predicate.
   - (b) junction overhang: for each junction pair from the Phase 2 classification, flood-fill
     the extended (guest) fault's triangles from its extension strip across shared edges,
     stopping at the corefined intersection polyline; drop the filled (overhang) component.
   - (c) ambiguity guard at the polylines: 1e-6 margin; triangles with any vertex within 1e-6 of
     the bounding surface and centroid on the drop side are dropped; report counts per fault and
     per class. The strips are engineering scaffolding — every overhang must go.
4. Post-clip weld + validation (in the same script): re-run the quantized dedup at 1e-3 m
   (autorefine may emit exact-duplicate points), then:
   - boundary closure: every boundary-marked edge incident to exactly 2 boundary triangles;
   - fault borders: every fault-marked border edge (1 incident fault tri) must lie on the DEM
     (trace), on the bottom plane, on a corefined junction polyline, on a ribbon clamp, or be a
     tip border — classify and count each; orphans = hard failure;
   - all edges >= 100 m (after Phase 3 cleanup and the item-2 autorefine gate, nothing new
     should violate).
5. Output: `meshing/data_corefined/safv4_<R>m/safv4_merged_<R>m.stl` + `_markers.json`.

### Interfaces
- Output pair feeds `tetgen_mesh.py --merged-stl ... --markers-json ...` (existing interface).

### Edge Cases to Handle
- Fault triangles split along the trace/junction polylines (done by Phase 3 corefine; autorefine
  must find nothing new there per the item-2 gate): the clip operates on post-autorefine
  triangles, so pieces are correctly kept/dropped by predicates (a)/(b).
- A clipped fault leaving a sliver triangle just below the DEM (height < 100 m between trace and
  first interior vertex): detected by the edge >= 100 m check; if found, collapse via
  quality_repair.h pass or report (same stop-and-report rule).
- Exactly-coplanar strip/DEM patches (degenerate tangency): autorefine handles coplanar overlap,
  but the clip's above-test is ambiguous at 0; the 1e-6 margin plus "any-vertex-on-DEM and
  centroid-above -> drop" rule resolves it deterministically.

### Acceptance Criteria
- [ ] `--box` regression (structural — the archived corefined STLs were excluded from the
      archive and the fixture is REGENERATED in Phase 0, so exact-count equality is not
      verifiable): merged soup self-intersection-free; 6 fault markers present; per-marker
      triangle counts within 2% of the regenerated fixture's corefine output.
- [ ] Autorefine resolved-intersection count involving any fault marker == 0.
- [ ] Boundary shell closed (0 open boundary edges).
- [ ] 0 fault triangles outside the closed shell (ray-parity re-scan after clip) AND 0 above
      the DEM.
- [ ] Every junction border edge lies on its corefined junction polyline.
- [ ] `check_self_intersect` clean on the merged soup.
- [ ] All edges >= 100 m.

### Dependencies
- Depends on: Phase 3.
- Required by: Phase 5.

---

## Phase 5: tetgen volume mesh, output contract, full gates

### Goal
The final SAFv4 volume mesh exists as Gmsh v2.2 `.msh` + VTUs and passes G0-G3 with the fault
uniformly ~500 m, slivers within the bar, and zero duplicated nodes.

### Files to Create
- Output artifacts (gitignored): `meshing/results/safv4_500m_topo.msh`,
  `meshing/results/safv4_500m_topo_{bulk,fault}.vtu`, gate JSONs + run log (committed).

### Files to Modify
- `code_preprocess/tetgen_mesh.py`:
  1. Accept the boundary-marker layout from Phase 4 (named boundary markers instead of the single
     box marker; the fault-aware sliver filter must treat ONLY fault markers as protected).
  2. Add `.msh` v2.2 output: `write_gmsh22(out_path, points, tets, fault_tris_by_patch,
     boundary_tris_by_name, tag_map)` using `meshio.write(..., file_format="gmsh22")`
     with `gmsh:physical`/`gmsh:geometrical` cell-data. Tag map: default below, configurable via
     `--tag-map JSON`.
     Normative default tag map (faults ordered alphabetically by FULL CFM basename:
     `SAFS-SAFZ-MJVS-...` < `SAFS-SAFZ-MULT-...` < `SAFS-SAFZ-SBMT-...`):
     ```
     volume rock = 1
     Physical Surface 101 = SAFS-SAFZ-MJVS-San_Andreas_fault-CFM6
     Physical Surface 102 = SAFS-SAFZ-MULT-Banning_fault-CFM6
     Physical Surface 103 = SAFS-SAFZ-SBMT-Garnet_Hill_fault-CFM6
     Physical Surface 201 = top (dem)
     Physical Surface 202 = bottom
     Physical Surface 203 = sides (all 4 ribbons under one tag)
     ```
     (Open question 1 below: confirm against the spatial-driver config before first production
     use; the alternative project uses fault=101/top=102/bottom=103/sides=104 for ONE fault —
     that exact map cannot extend to 3 faults, hence the new 2xx block for boundaries.)
  3. Fault triangles in the `.msh` must reference the SAME node IDs as the adjacent tets (single
     shared node set — meshio does this naturally from one points array; the gate is G3a/G3b).
3. tetgen invocation (archived defaults unless stated): `plc=True, nobisect=True, quality=True,
   minratio=2.0`, graded background sizing `--lc-near 500 --lc-far 5000 --dist-inner 2000
   --dist-outer 20000 --bg-spacing 1000` (bgmesh machinery already in the archived script),
   `--bulk-min-edge 100 --bulk-q-min 0.05 --vertex-merge-tol 1e-3 --output-dedup-tol 99`.
   NOTE `--vertex-merge-tol` changes from the archived 10 m to 1e-3 m: Phase 2-4 already welded
   the soup exactly; a 10 m input merge could corrupt 100 m-scale trace features.
4. Run ladder (each step gated before the next):
   a. **Flat-box regression**: 2000 m fixture through the updated script -> metrics match
      STATUS_tetgen_pipeline.md within noise (n_tets ~33.6 k, embedding 100%, coverage 100/100).
   b. **1-fault smoke**: Banning only (smallest, 352 tris) + dem/ribbons/bottom, target 2000 m.
   c. **3-fault coarse**: full SAFv4 at 2000 m fault target.
      Rung (c) is also the G2 calibration point: record count(eta<=0.05) and count(eta<=0.1)
      with the Phase-1 localization report. If count(eta<=0.05) > 0, STOP and present the report
      with a recommendation (retune Phase-3 sizing vs accept-localized) BEFORE rung (d); G2a
      applies to rung (d) as agreed at that decision point. G2a remains the default target.
   d. **Production**: 500 m fault target.
   All resolution-dependent knobs are functions of the rung target R (so coarse rungs exercise
   the same code paths they de-risk): fault mesh_edge_size = R, polyline_spacing = R/2,
   DEM h_near = R / h_far = 5R (d-window unchanged), tetgen lc-near = R / lc-far = 10R,
   G1f bands scaled by R/500. Record the resolved values in each rung's gate JSON.
5. Gates on (c) and (d): run `check_mesh_quality.py` (G1, G2a/b, fault sizing table) and
   `check_fault_conformity.py` (G3a/b/c, G0 coverage via 500-point sampling — add the point-in-tet
   sample to the conformity tool if simpler than a separate script; STATUS doc already had this
   test). Commit both JSON outputs + the run log.
6. MFEM-readability smoke: load the coarse (c) `.msh` via the existing seas test pattern (any
   driver `--dry-run` config or a 20-line `mfem::Mesh m(path)` snippet compiled against the seas
   Makefile) to prove the v2.2 file parses. Do NOT run a simulation (production runs are
   Frontera-only per project policy).

### Interfaces
- `meshing/results/safv4_500m_topo.msh` is the deliverable consumed by future SAFv4 spatial-driver
  configs (tag map documented in the manifest and this plan).

### Edge Cases to Handle
- **tetgen PLC rejection** (the central risk): if tetgen reports a PLC error on (b) or (c), STOP.
  Capture tetgen's offending-facet diagnostics, map facet IDs back through the markers JSON to a
  surface + location, and report. Do not silently switch meshers (CLAUDE.md rule).
- **nobisect starving quality near the trace**: with `nobisect=True` tetgen cannot insert Steiner
  points on the input surfaces; if G2 fails concentrated at the trace, that indicates Phase 3
  sizing (DEM grading / polyline spacing) needs retuning — iterate Phase 3 knobs, not tetgen.
- **Output dedup dropping fault triangles** (fixture lost 8/8,202): re-run the embedding gate
  after dedup; > 0.2% loss on any single fault patch -> stop and report.

### Acceptance Criteria
- [ ] G0, G1, G2a, G3a, G3b, G3c all PASS on the production 500 m mesh (hard).
- [ ] G2b: eta <= 0.1 fraction <= 5e-4; report the exact count, locations (z-histogram,
      boundary attribution) regardless.
- [ ] G1f fault sizing table per patch in [425, 600] m median.
- [ ] Flat-box regression (4a) within noise of the STATUS baseline.
- [ ] `.msh` parses in MFEM (smoke) and in `meshio` round-trip; file declares `$MeshFormat 2.2`.
- [ ] Gate JSONs + log committed; large artifacts gitignored.

### Dependencies
- Depends on: Phases 1, 4.
- Required by: Phase 6 (only if G2b misses).

---

## Phase 6 (CONDITIONAL): sliver elimination escalation

Enter ONLY if Phase 5 G2b fails. Do not enter silently — present the Phase 5 sliver report first.

Ranked levers (evidence from STATUS_tetgen_pipeline.md "Mitigation attempted"):
1. Retune Phase 3 sizing (DEM grading window, polyline_spacing 250 -> 150, fault min_edge) —
   attacks the root cause (surface density mismatch), highest expected payoff.
2. tetgen `optim` + `mindihedral` sweep (documented marginal: q_med 0.871 -> 0.892, slivers
   70 -> 60 on the fixture) — cheap, bounded gain.
3. Targeted local repair via `quality_repair.h` (collapse shortest interior edge of each sliver,
   fault-protected) — exists, needs a driver loop; medium effort.
4. MMG3D is BANNED (drops fault triangles on non-manifold input — documented failure).
Accept-with-report is the final fallback: slivers localized in the far field (bottom layer) may be
acceptable for simulation; that is the user's call with the localization report in hand.

---

## Testing Strategy

- **Phase 1 is the test harness**: self-test (synthetic duplicate/embedding fixtures), committed
  baseline JSON on the known-bad GOCAD mesh (regression-locks the diagnosis), and gate JSONs on
  every produced mesh.
- **Golden regression**: the archived 2000 m flat-box fixture is re-run after every change to
  `autorefine_merged.cpp` / `tetgen_mesh.py` (Phases 4-5 acceptance) — it is the only proven
  end-to-end reference (STATUS_tetgen_pipeline.md numbers).
- **Run ladder** (Phase 5): flat-box regression -> 1-fault smoke -> 3-fault coarse (2000 m) ->
  production (500 m). Every rung gated; failures stop the ladder with a report.
- **Determinism**: manifests and gate JSONs byte-stable across re-runs (md5 in acceptance);
  no timestamps inside generated artifacts.
- **Visual checks**: per-fault VTUs + bulk VTU at each rung (ParaView recipe as in
  gocad_inp_to_vtu.py docstring); slice through each fault to confirm embedding and the trace.

## Risk Assessment

1. **tetgen PLC failure on the topo-bounded domain (highest risk).** Proven only on a flat box.
   Mitigated by the run ladder (1-fault smoke isolates), autorefine's intersection resolution
   (the fixture's proven fix for PLC inputs), and the closed-shell/self-intersection gates in
   Phase 4 BEFORE tetgen ever runs. If it still fails: stop-and-report with facet diagnostics.
2. **Near-tangent fault/DEM intersections** producing trace splinters that cleanup cannot collapse
   without breaking conformality. Detector: Phase 3 acceptance (all edges >= 100 m + max_diff 0).
   Fallback decision for the user: locally relax the 100 m floor along the trace vs. locally
   simplify the DEM.
3. **Minus-side-as-canonical assumption.** Measured to hold today (subset up to 6 tris); Phase 2
   aborts loudly if a re-exported model breaks it.
4. **isotropic_remeshing self-intersections on the extension strip** at sharp trace bends.
   Gate: `check_self_intersect` per STL (Phase 3); fix by raising `--h-ext` smoothing or local
   polyline resampling.
5. **meshio gmsh22 writer fidelity** (physical tags, element ordering). Gate: MFEM read smoke +
   meshio round-trip; known project gotcha is v4-vs-v2.2, addressed by writing `gmsh22` explicitly.
6. **Scale/memory.** 31 M-tet meshio reads cost ~8 GB / 2.5 min locally (measured); the rebuilt
   mesh is the same order. Python steps stay chunked (2e6 tets/chunk). No full simulation runs
   locally (Frontera only).
7. **CGAL env drift** since the May archive — Phase 0 surfaces it immediately, before any new code.

## Open Questions for the User

1. **Physical tag map** (Phase 5): proposed faults 101/102/103 (alphabetical) + top 201 /
   bottom 202 / sides 203. The single-fault alternative-project convention (101/102/103/104)
   cannot represent 3 faults. Confirm before the first production mesh feeds a driver config.
2. **Sliver bar** (G2b): recommended start = hard zero for eta <= 0.05 plus <= 5e-4 fraction in
   (0.05, 0.1]. Driving (0.05, 0.1] to literal zero may cost Phase 6 iterations. Acceptable?
3. **Production resolution**: 500 m on-fault confirmed (estimate ~5-20 M tets)? 250 m is known to
   OOM the driver and is out of scope without explicit approval.
4. **Ribbon/bottom reuse**: plan keeps GOCAD's ribbon/bottom triangulations (quality acceptable,
   far-field). If you prefer coarser/cleaner regenerated side walls, that is a Phase 3 extension —
   say so before implementation starts.
