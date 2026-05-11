# Implementation Plan: revert volume meshing from CGAL Mesh_3 to gmsh, using a pre-stitched non-manifold conform surface

Anchor: `/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/safs/project_7.0_preferred/`. All paths below are relative unless absolute.

## Overview

CGAL `Polyhedral_complex_mesh_domain_3` cannot reliably mesh our 6 open fault patches without leaving 5-15% of the rock domain as cavities (regardless of subdomain pair `(1,1)` or `(1,2)`); the graded sizing field then turns those cavities into visible holes (Screenshot 2026-05-10 09:05). Revert to gmsh for the 3-D step but **avoid the failed `Surface{} In Volume{}` multi-fault path** (tetgen PLC error). Instead, pre-stitch the 6 corefined STLs into ONE non-manifold combined surface mesh that includes the bounding box, then feed it to gmsh as a single `Merge`'d discrete entity. tetgen handles non-manifold edges *inside* a single surface differently from across multiple `Merge`d surfaces — the alternative project's single-fault `freesurface_clip.geo` (`code_meshing/safs_fault_box_freesurface_clip.geo:91-96`) is proof of concept for the "single discrete surface, embedded via `Surface In Volume`" pattern.

This plan **preserves all CGAL artifacts** (`code_preprocess/corefine_cgal/mesh_volume.cpp`, `distance_sizing_field.h`, `*.mesh` outputs) untouched. The CGAL path remains as a fallback.

## Constraints

### Interface constraints
- **`data_corefined/` is FROZEN.** The 6 corefined STL fault meshes and `manifest.json` (Phase 1 + Phase 2 corefine output, with `polyline_cleanup` + `cross_polyline_snap`) are the canonical surface input. They have `min_edge ≥ 100 m` per fault and `max_diff_m = 0` for every pair. Do NOT regenerate; do NOT modify.
- **No CGAL Mesh_3 artifacts get deleted.** `mesh_volume.cpp`, `distance_sizing_field.h`, `c3t3_min_edge_cleanup.h` (if added later), and the `*.mesh` outputs remain in place. They are now a fallback path, not the canonical one.
- **Output contract**: gmsh `.msh` v2 or v4 ASCII at `code_meshing/safs_multifault_box_<R>m.msh` with Physical groups: `Volume("rock", 1)`, `Surface("fault_<NAME>", 101..(100+N))` per fault, `Surface("top", 200)`, `Surface("bottom", 201)`, `Surface("sides", 202)`. The downstream `msh_to_vtu.py` (port from `../project_7.0_alternative/code_preprocess/msh_to_vtu.py`) splits to `<base>_bulk.vtu` + `<base>_fault.vtu`.

### Dependency constraints
- gmsh **4.15.0** in `pythonenv` (`/Users/chunhuizhao/miniforge/envs/pythonenv/bin/gmsh`). Python API also available (`import gmsh` works).
- `meshio` for STL/MSH/VTU I/O.
- No tetgen on PATH — option E (direct tetgen .poly) requires `conda install -c conda-forge tetgen` first.

### Convention constraints
- Per `CLAUDE.md` (project root): do not silently fall back to simpler approaches without asking. If a phase fails, stop and report.
- Preserve the alphabetical fault ordering for Physical Surface tag stability (101..(100+N) tied to alphabetical basenames; established in `PLAN_cgal_corefine_multifault.md` §Phase 3 lines 449-456).
- Build any C++ helper via the existing `code_preprocess/corefine_cgal/CMakeLists.txt`.

### Numerical constraints
- **Hard floor**: every tet edge ≥ 100 m. (Currently met at 96.37 m by the CGAL Phase 1 path; user previously accepted that as effectively zero — same bar applies here.)
- `min q_tet ≥ 0.05`, `≥ 99.5 % q_tet ≥ 0.3`, `q_med ≥ 0.85`.
- **Coverage**: 0 cavities — every box-interior point must be inside a tet labeled subdomain "rock".
- **Conformality** preserved: every input STL polyline vertex (from `manifest.pairs[k].polylines`) must match a vertex in the output `.msh` within 1e-3 m (gmsh `Geometry.Tolerance`).
- **Graded sizing** still required: `lc_near = 1500 m` near faults → `lc_far = 15000 m` far field, ramp 3 km → 40 km. The original gmsh plan §Phase 3.5a Distance + Threshold field is still the recipe.

## Survey of options (evaluated up front so the user can pick)

| Option | Idea | Pre-condition | Risk | Code change |
|---|---|---|---|---|
| **A (Recommended)** | Pre-stitch the 6 STLs + box-as-STL into ONE non-manifold combined STL; `Merge` once; `Surface{combined_tag} In Volume{1}` | corefine output (already conformal, `max_diff = 0`) | Low — same recipe as alternative project's single-fault, just with non-manifold seams along polylines | New stitcher script; modify `.geo` generator to use one Merge |
| B | `ClassifySurfaces(angle, …, true)` after `Merge` of the 6 STLs, then `BooleanFragments` to split the box volume by the faults | Needs OCC kernel + classified surfaces to be water-tight enough to fragment | High — banned for single-fault by the original plan; multi-fault failure mode unknown; OCC + discrete is fragile | Significant rewrite of the .geo generator |
| C | Define each fault as `Compound{Surface{...}}` and embed | Compound surfaces are an OCC concept; fault is discrete (STL) — Compound on discrete may not be supported | Unknown | Moderate |
| D | Switch `Mesh.Algorithm3D = 10` (HXT) or `7` (MMG3D) instead of tetgen `1` | Already tried HXT with `Surface In Volume` — silently dropped surfaces ("Unknown surface 9-12") | High — likely same failure | Low |
| E | Generate tetgen `.poly` directly + call tetgen | Need `tetgen` binary installed | Low risk per se, but adds dep + new code path | Significant new tooling |

**Choose Option A.** Reasoning: (i) it's a small, principled extension of the proven single-fault recipe; (ii) the conformality contract from corefine (`max_diff = 0`) makes vertex stitching a coordinate-key match — no geometric ambiguity; (iii) it leaves the 3-D mesher (tetgen, `Algorithm3D = 1`) on the path it works on; (iv) options B/C/D have higher unknowns and option E adds a dependency. Phases 2-4 below address followups.

---

## Phase 1: stitch the 6 corefined STLs + bounding box into one combined STL

### Goal
Produce one STL file (`data_corefined/safs_combined_<R>m.stl`) containing the union of: 6 corefined fault triangles + 12 bounding-box triangles, with vertex deduplication so coincident polyline vertices and box-fault contact vertices share a single index. Output is one closed-volume polyhedral surface with non-manifold internal edges along every fault-fault polyline. Built deterministically from `data_corefined/manifest.json` so the file is byte-stable across re-runs.

### Files to Create
- `code_preprocess/stitch_combined_stl.py` — Python script using `meshio` and `numpy`. Reads the 6 corefined STLs and `manifest.json`; emits one combined STL.

### Files to Modify
- None.

### Detailed Requirements
1. **Input parsing**:
   - CLI: `python stitch_combined_stl.py [--manifest path] [--out path] [--res 2000] [--pad-xy M] [--pad-top M] [--pad-bottom M]`. Defaults derived from `corefine_faults.py` pattern (`code_preprocess/corefine_faults.py:107-112`).
   - Default `--manifest = <project>/data_corefined/manifest.json`.
   - Default `--out = <project>/data_corefined/safs_combined_<R>m.stl`.
   - Defaults: `pad_xy = 50000`, `pad_top = 100`, `pad_bottom = 25000` (carried over from `generate_multifault_geo.py:30-32`).
2. **Vertex dedup**:
   - Build a global point list. For each input STL triangle, look up each vertex by an integer-quantized coordinate key `(round(x/tol), round(y/tol), round(z/tol))` with `tol = 1e-6 m`. If key exists, reuse the index; else append.
   - Do the same for the 12 bounding-box triangles.
   - Justification: corefined STLs share polyline vertex coords bit-identically (post Phase 2 cleanup). At 1e-6 m tolerance, no false matches and no missed matches.
3. **Bounding box construction**:
   - Compute `xmin..zmax` from `manifest.meshes[].bbox` ± paddings (same recipe as `generate_multifault_geo.py:67-75`).
   - 8 corners + 12 outward-oriented triangles (mirror of `mesh_volume.cpp:218-231`).
   - Box vertices joined to the global point list via the same coord-key dedup; 4 box corners that touch a fault perimeter (if any) get reused, else they are new vertices.
4. **Triangle ordering**:
   - Box triangles emitted FIRST (12 entries), so they get the lowest STL solid offset.
   - Then the 6 fault STLs in alphabetical basename order (matches `manifest.meshes[]` traversal in `generate_multifault_geo.py:90`).
5. **STL writing**:
   - Use `meshio.write(out_path, mesh, file_format="stl-ascii")` with `points` and a single `("triangle", connectivity)` cell block.
   - The output STL is a single `solid` with all triangles in one block. Per gmsh's STL reader, this becomes one discrete surface with one tag after `Merge`.
6. **Per-fault provenance file** (needed by Phase 2):
   - Also emit `data_corefined/safs_combined_<R>m_provenance.json`:
     ```json
     {
       "n_box_triangles": 12,
       "faults": [
         {"basename": "...", "first_tri": 12, "n_tri": 788},
         ...
       ]
     }
     ```
   - Ranges are `[first_tri, first_tri + n_tri)` in the combined STL's triangle order.
7. **Idempotency**:
   - Re-running with the same `manifest.json` produces a byte-identical STL + provenance. Verifiable via `md5`.

### Interfaces
- New Python module `stitch_combined_stl.py` with a `main(): int` entry point.
- New artifact `data_corefined/safs_combined_<R>m.stl` and `data_corefined/safs_combined_<R>m_provenance.json` — consumed by Phase 2's .geo generator.

### Edge Cases to Handle
- **A fault perimeter vertex coincides with a box corner**: the dedup at 1e-6 m matches them and they share an index. Behavior: the box face triangulation has a "pin" vertex — fine for tetgen as long as the box face triangle still has positive area. Verify in Phase 1 acceptance criterion 4.
- **A fault triangle has a degenerate (sub-1mm) edge**: per Phase 2 corefine cleanup, all input STL edges are ≥ 96.37 m. No degenerate triangles to worry about; assert `min_edge ≥ 50 m` as a safety check, fail loudly if violated.
- **Two faults share a polyline edge**: by corefine conformality (`max_diff = 0`), the two endpoints dedup to the same global indices. Two fault triangles thus share an edge in the global combined STL — that edge is **non-manifold** (4 incident triangles: 2 from each fault). This is the desired output. gmsh handles it.
- **Triple-junction vertex (3 faults meeting at one point)**: the dedup walks correctly, all three faults' triangles share that vertex. Resulting non-manifold structure is what Phase 2 needs.

### Acceptance Criteria
- [ ] `python stitch_combined_stl.py --res 2000` exits 0 and writes both files.
- [ ] Re-running gives byte-identical output (MD5 stable).
- [ ] Combined STL triangle count = 12 + sum of `manifest.meshes[k].n_faces` (for canonical 2000 m: 12 + 8152 = 8164).
- [ ] Combined STL vertex count = (8 box corners) + (sum of fault vertex counts) − (number of dedup matches). Verify number of dedups = number of box-fault coincidences + number of fault-fault polyline-vertex pairs (each shared once globally). For the canonical 2000 m fixture: ≈ 8 + (sum n_verts_in_corefined) − N_polyline_shared. Sanity-check against `manifest.pairs[k].shared_a`.
- [ ] `meshio.read(combined.stl)` then count edges with both endpoints in different fault provenance ranges = polyline-edge count from `manifest`. Cross-check pairwise; if mismatch, the dedup is broken.
- [ ] Provenance JSON parses with `python -m json.tool` and ranges sum correctly.

### Dependencies
- Depends on: corefine output (frozen).
- Required by: Phase 2.

### Risk
- **Low**. Vertex dedup at 1e-6 m on corefine-conformal input is mechanical. The only risk is the dedup tolerance — verify with the cross-check in acceptance criterion 4.

---

## Phase 2: regenerate `.geo` to use the single combined STL

### Goal
Modify `code_preprocess/generate_multifault_geo.py` (or add a sibling) so that the produced `.geo` does ONE `Merge "safs_combined_<R>m.stl"` and ONE `Surface{tag} In Volume{1}`, replacing the failed loop over 6 individual `Merge`s + 6 `Surface In Volume`s. Per-fault Physical Surface tags are still emitted, but they reference *triangle id ranges* of the single discrete surface (or, if gmsh doesn't support that directly, fall back to using `Physical Surface("fault_all", 100)` and let downstream split by patch_id).

### Files to Create
- None (modify existing).

### Files to Modify
- `code_preprocess/generate_multifault_geo.py` — replace the per-fault Merge loop (currently around `generate_multifault_geo.py:154-160`, the `# 4. Merge each corefined fault.` block) with a single `Merge "../data_corefined/safs_combined_<R>m.stl";` line. Drop the For loop emitting 6 `Surface{} In Volume{}`; emit one `Surface{7} In Volume{1};`. Keep everything else (Box, size field, Physical groups for box top/bottom/sides, Mesh.Algorithm = 6 / Algorithm3D = 1).
- `code_meshing/safs_multifault_box_2000m.geo` — regenerated artifact.

### Detailed Requirements
1. **STL filename**: replace `f"../data_corefined/{mesh['basename']}_corefined.stl"` with `f"../data_corefined/safs_combined_{res}m.stl"`. Drop the loop.
2. **Discrete surface tag**: after `Merge`, gmsh assigns the single STL-derived discrete surface tag `7` (one past box faces 1..6). Hard-code that.
3. **Embed**: replace the `For k In {0:N_FAULTS-1} Surface{fault_surfs[k]} In Volume{1}; EndFor` block with `Surface{7} In Volume{1};`.
4. **Physical groups (per-fault)** — there are TWO sub-options. The implementer picks one and DOCUMENTS the choice in the .geo header:
   - **2.1 — Single-tag fault group** (simplest): `Physical Surface("fault_all", 100) = {7};`. Downstream consumers read the provenance JSON to split by fault.
   - **2.2 — Per-fault tag via gmsh `setMassive` / triangle-id list**: requires gmsh's geometry-by-id syntax. Investigate if gmsh's `.geo` supports `Physical Surface("fault_X", 101) = elements{...};` over a triangle-id list. If yes, use it; if not, fall back to 2.1.
   - **Pick 2.1** as the default. Phase 4's Python post-processor uses the provenance JSON for per-fault splitting.
5. **Graded size field** (porting from the original PLAN_cgal_corefine_multifault.md §Phase 3.5a): emit a Distance/Threshold size field block. Distance is computed against the discrete surface tag 7 minus the box face contributions — need to think about how:
   - Plan A: `Field[1] = Distance; Field[1].SurfacesList = {7};` — but 7 includes the box faces, so the size near the box surface is also `lc_near` — not what we want.
   - Plan B: emit auxiliary `Point(...)` and `Line(...)` entities for the polyline vertices (from `manifest.pairs[k].polylines`) and use `Field[1].PointsList = {polyline_pt_ids[]}` instead. This is the original plan §Phase 3.5a' fallback (PLAN_cgal_corefine_multifault.md §Sub-phase 3.5a' lines 600-601).
   - **Pick Plan B** because it explicitly targets fault interfaces only and matches the CGAL Mesh_3 distance-to-fault-surface semantic we just removed.
6. **Mesh.Optimize block** (PLAN_cgal_corefine_multifault.md §Phase 3.5b lines 619-632):
   ```gmsh
   Mesh.Optimize = 1;
   Mesh.OptimizeNetgen = 1;
   Mesh.OptimizeThreshold = 0.3;
   Mesh.HighOrderOptimize = 0;
   ```
   Append at the end of the .geo. Mandatory per the original plan.
7. **Idempotency**: keep the existing `generate_multifault_geo.py` byte-stable property (no embedded timestamp).

### Interfaces
- `generate_multifault_geo.py` API unchanged (same CLI). The output `.geo` content changes.
- Phase 2 output is consumed by Phase 3 (`gmsh -3` invocation).

### Edge Cases to Handle
- **gmsh assigns a different tag than 7 to the merged STL** (e.g. if Mesh.STLOneSolidPerSurface is set). Workaround: query `gmsh.model.getEntities(2)` after Merge in a Python wrapper to find the discrete surface tag dynamically. This adds complexity; defer to Phase 4 unless Phase 3 trips on it.
- **Polyline-vertex Point() entities collide with OCC's Box(1) point IDs**: use base offset `100000` for polyline points and `200000` for polyline lines (PLAN_cgal_corefine_multifault.md lines 547-549).

### Acceptance Criteria
- [ ] `python generate_multifault_geo.py --res 2000` produces `code_meshing/safs_multifault_box_2000m.geo` whose MD5 is stable across re-runs.
- [ ] The .geo contains exactly one `Merge` line referencing `safs_combined_2000m.stl`.
- [ ] Diff against the previous (gmsh-failing) .geo shows: per-fault Merges replaced by one Merge; per-fault Physical Surfaces replaced by `fault_all`; Distance field on PointsList; Mesh.Optimize block added.

### Dependencies
- Depends on: Phase 1.
- Required by: Phase 3.

### Risk
- **Low** for the .geo regeneration itself. The unknown is whether gmsh's discrete-surface tagging cooperates (covered by Phase 3 acceptance).

---

## Phase 3: run gmsh -3 and verify mesh quality

### Goal
Produce `code_meshing/safs_multifault_box_2000m.msh` via `gmsh -3 safs_multifault_box_2000m.geo` exiting 0, with the volume mesh meeting the hard quality contract.

### Files to Create
- None (run gmsh).

### Files to Modify
- None.

### Detailed Requirements
1. **Invocation**: `cd code_meshing && gmsh -3 safs_multifault_box_2000m.geo -o safs_multifault_box_2000m.msh -v 3 > gmsh.log 2>&1`. Run in background per project convention.
2. **Diagnostics on failure** (PLAN_cgal_corefine_multifault.md §Phase 3.5c):
   - Inspect `gmsh.log` for "PLC Error", "No closed volume", "Invalid boundary mesh".
   - Run a 2-D sanity check first via the gmsh Python API: `gmsh.model.mesh.generate(2)` and assert the merged surface triangle count matches `12 + sum(manifest.meshes[].n_faces)`.
3. **Quality verification** via `code_preprocess/check_msh_quality.py` (already exists for MEDIT; extend or write a parallel tool for `.msh`):
   - shortest edge ≥ 100 m (hard floor)
   - `min q_tet ≥ 0.05`
   - `≥ 99.5 %` of tets have `q_tet ≥ 0.3`
   - `q_med ≥ 0.85`
   - 0 cavities (every random box-interior point inside a tet — same point-in-tet sampling as the CGAL debug)

### Interfaces
- gmsh CLI invocation; no new code unless 2-D diagnostic fails.

### Edge Cases to Handle
- **gmsh -3 fails with PLC error like the previous Surface-In-Volume run**: that would mean Option A also fails. Stop, report, and ask the user whether to try Option B (BooleanFragments), C (Compound), D (different mesher), or E (direct tetgen).
- **gmsh -3 succeeds but the .msh has no fault Physical Surface output**: the Physical Surface assignments use the wrong tag. Fix in Phase 2 (re-emit .geo) and re-run.
- **gmsh -3 produces a mesh with cavities** (same symptom as CGAL): `point-in-tet` test fails. Indicates the non-manifold internal edges weren't honored — probably a sign that `Surface In Volume` semantically can't handle non-manifold. Stop and escalate.

### Acceptance Criteria
- [ ] `gmsh -3` exits 0 in ≤ 30 min on the 2000 m fixture.
- [ ] gmsh.log contains `0 ill-shaped tets` after `Mesh.Optimize`.
- [ ] `meshio.read('safs_multifault_box_2000m.msh')` succeeds; mesh has `tetra` cells in Physical group `rock`.
- [ ] All quality thresholds above are met.
- [ ] Cavity test (500 random box-interior points): 500/500 inside a tet.
- [ ] Physical Surface "fault_all" exists with at least `sum(manifest.meshes[].n_faces) = 8152` triangles (or whatever count post-gmsh-2D-mesher).

### Dependencies
- Depends on: Phases 1, 2.
- Required by: Phase 4.

### Risk
- **Medium**. We don't know empirically whether single-merged-STL bypasses the multi-fault PLC issue. The hypothesis (single-discrete-surface non-manifold edges are handled differently from cross-Merge non-manifold edges) is plausible but unverified. **Smoke-test on a 2-fault subset before the full 6-fault run** — same approach we used for the CGAL debug (PLAN §Phase 3 line 503).

---

## Phase 4: per-fault VTU split via the provenance JSON

### Goal
Port `../project_7.0_alternative/code_preprocess/msh_to_vtu.py` to this project so the gmsh `.msh` is split into `_bulk.vtu` (all rock tets) + per-fault `_fault_<NAME>.vtu` files. Use the provenance JSON from Phase 1 to map gmsh triangles back to fault basenames.

### Files to Create
- `code_preprocess/msh_to_vtu.py` — port from alternative project (`../project_7.0_alternative/code_preprocess/msh_to_vtu.py`), adapted for the single-`fault_all` Physical group + provenance-JSON splitting.

### Files to Modify
- None.

### Detailed Requirements
1. **CLI**: `python msh_to_vtu.py INPUT.msh OUTPUT_BASE [--bulk-name rock] [--fault-name fault_all] [--provenance path.json]`.
2. **Bulk VTU**:
   - Read `.msh` via `meshio`. Extract tets in Physical Volume `rock`. Write to `<base>_bulk.vtu`.
   - Cell data: `quality` (tet quality, computed by the existing `tet_quality()` function in the alternative project's tool).
3. **Per-fault VTU split**:
   - Triangles in Physical Surface `fault_all`: each triangle has a position in gmsh's mesh. Match each triangle's CENTROID to the `safs_combined_<R>m.stl` triangle list using a kd-tree of pre-stitch centroids.
   - For each matched STL triangle index, look up the fault range in `safs_combined_<R>m_provenance.json`.
   - Group triangles by fault basename → write one VTU per fault: `<base>_<fault_basename>.vtu`.
4. **Box surface VTU** (optional, for visualization): emit `<base>_box.vtu` with the Physical Surface "top"/"bottom"/"sides" triangles. Optional flag `--with-box`.

### Interfaces
- CLI matches the alternative project's tool with the additions noted.
- Output VTUs consumable by ParaView; carry per-cell `quality` and `gmsh:physical` data fields.

### Edge Cases to Handle
- **A gmsh-output triangle's centroid does not match any pre-stitch STL triangle** (because gmsh's surface mesher may have refined the discrete surface): this happens when gmsh inserts Steiner points on the merged STL. Strategy: walk gmsh's triangle vertex list; if all 3 vertices match pre-stitch points, use STL-triangle-index lookup; if not, use NEAREST pre-stitch triangle by centroid distance and tag with that fault.
- **Number of gmsh fault triangles doesn't match the input STL count**: if gmsh refined or coarsened, the count differs. Report as a warning; the matching algorithm above still works as long as the underlying surface is preserved.

### Acceptance Criteria
- [ ] One `_bulk.vtu` + N `_<fault_basename>.vtu` files written for the 6-fault dataset.
- [ ] Sum of per-fault triangle counts equals total `fault_all` triangle count from the .msh.
- [ ] ParaView can open all files; coloring by `quality` shows distribution; per-fault VTU shows just one fault.

### Dependencies
- Depends on: Phase 3.
- Required by: nothing.

### Risk
- **Low**. Mostly a port + a kd-tree match.

---

## Phase 5 (CONDITIONAL): fall back to Option B/C/D/E if Phase 3 fails

If `gmsh -3` (Phase 3) fails with the same PLC error as the original Surface-In-Volume multi-fault attempt, halt and report. Do NOT autonomously try Options B-E — they have higher risk and the user should choose.

The ranking if escalation is needed:
1. **Option E** (direct tetgen via `.poly`): adds a tetgen install but gives full control over PLC marking. The clearest path to a working multi-fault tet mesh, with the ability to mark non-manifold edges explicitly.
2. **Option D** (different gmsh mesher): cheap to try (single line change). Likely to fail with similar symptoms but may surprise.
3. **Option B** (BooleanFragments): higher risk per the original plan; needs OCC-able classified surfaces.
4. **Option C** (Compound surfaces): unclear if discrete-mesh Compound is supported.

---

## Testing Strategy

### Per-phase smoke
- After Phase 1: assert MD5 stability + provenance JSON sums match expected fault triangle counts.
- After Phase 2: visually diff the new `.geo` against the old; confirm `Merge` count = 1.
- After Phase 3: 2-fault subset smoke test FIRST (subset = SBMT-MC + SBMT-SAF, which previously triggered the PLC error). If 2-fault works, run 6-fault.
- After Phase 4: load every output VTU in ParaView; check fault coloring.

### Regression baseline
The CGAL Phase 1 mesh (4.25M tets, q ≥ 0.3 fraction = 99.95%, 0 cavities-by-quality-test, 16 q < 0.05 tets, edge_min = 96.37 m) is the regression baseline from the user's previous "actually good" verdict. Any new mesh worse than that on **any** of (cavity coverage, q ≥ 0.3 fraction, edge_min, q_med) is rejected.

### Reference / cross-verify
The alternative project's single-fault `freesurface_clip.geo` produces a working mesh of `[179..622] km × [3.69..3.99] km × [-16..0] km` box with one fault. Use it as a recipe template for our Phase 2 .geo emitter.

## Risk Assessment

### Specific risks per phase
- **Phase 1**: dedup tolerance miss → false matches or missed matches → broken non-manifold structure. Mitigated by acceptance criterion 4 (cross-check polyline-edge count against manifest).
- **Phase 2**: graded size field via Points-only Distance gives a coarser approximation than Lines (worst-case error = polyline_spacing/2 ≈ 375 m at default — still under lc_near = 1500 m). Acceptable per PLAN_cgal_corefine_multifault.md §Variant 3.5a' lines 600-601.
- **Phase 3**: tetgen PLC error reappears (the central bet of this whole plan). Mitigation: 2-fault smoke first; if it fails, the Option-B/C/D/E escalation path is documented.
- **Phase 4**: gmsh's surface mesher may insert Steiner points on the discrete combined STL, breaking the index-based provenance lookup. Centroid-NN fallback handles this. If fault refinement is severe (e.g. > 2× input triangle count), the centroid match becomes ambiguous near polylines — warn but accept.

### Cross-cutting risks
- **`Surface{} In Volume{}` may silently fail with HXT/Frontal but tetgen rejects multi-fault**: we already saw this. Phase 3 sticks with tetgen (Algorithm3D=1) and bets that the single-discrete-surface formulation removes the PLC error. If it doesn't, the bet is lost.
- **gmsh discrete-surface tag drift**: gmsh may assign the merged STL a tag other than 7 (e.g. if `Mesh.STLOneSolidPerSurface = 1`). Phase 2 hard-codes 7; if gmsh disagrees, Phase 3 errors out clearly.

## What is explicitly OUT OF SCOPE

- Modifying corefine output (`data_corefined/`) — frozen.
- Modifying or deleting CGAL artifacts (`mesh_volume.cpp`, `distance_sizing_field.h`, `medit_to_vtu.py`, `*.mesh`) — kept as fallback.
- Changing the geological model (fault perimeters, fault-fault intersections).
- Implementing Phase 5 (Option B/C/D/E) without explicit user approval.

## Open Questions for the User

1. **Fault Physical Surface labeling style** (Phase 2 sub-option 2.1 vs 2.2): single `fault_all` tag with provenance-JSON-driven Python split (simple, recommended), or per-fault tags emitted directly in the .geo if gmsh syntax allows? Default = 2.1.
2. **Smoke-test on 2-fault subset before 6-fault** (Phase 3): yes (recommended, ~3 min cost) or no (just run 6-fault directly, ~30 min cost on failure)? Default = yes.
3. **Time budget for Phase 5 escalation**: if Phase 3 fails, are we authorized for Option E (install tetgen + write a `.poly` emitter, ~1 day)? Default = ask first.
