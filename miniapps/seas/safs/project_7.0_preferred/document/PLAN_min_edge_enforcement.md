# Implementation Plan: enforce hard min-edge ≥ 100 m on the CGAL Mesh_3 volume mesh

Anchor: `/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/safs/project_7.0_preferred/`. All paths below are relative unless absolute.

## Overview

The current `code_preprocess/corefine_cgal/mesh_volume.cpp` produces a 4.8 M-tet mesh that meets the plan's strictest criterion (99.94 % q_tet ≥ 0.3) but violates the hard floor `min(edge_length) ≥ 100 m` for **5,300 tets / 4.8 M = 0.11 %** (and 16 tets have q_tet < 0.05). All 5,300 sub-floor tets sit within 1 km of a fault surface; their root cause is the intrinsic 1D Steiner-point insertion that `Polyhedral_complex_mesh_domain_3` does on patch boundaries (the corefined polyline edges) regardless of whether `detect_features` / `detect_borders` are called.

This plan evaluates **five paths** to drive the sub-100 m count to 0, ordered by escalating scope. The recommended path is **Phase 1 (`edge_min_size` + `cell_min_size` named parameters)** because it is a one-call CGAL API change, ships in CGAL 6.1.1, and is documented to do exactly what we need. Phases 2-5 are fallbacks if Phase 1 leaves residuals.

This is a planning document. No code is written here.

## Constraints

### Interface constraints
- **Inputs are frozen** (Phase 2 corefine outputs in `data_corefined/`): 6 corefined STL faults + manifest. Surface edges already ≥ 100 m, all 7 pairs conformal at `max_diff_m = 0`. Do not regenerate the corefine output as part of this plan — the surface mesh contract is met.
- **Output contract** (downstream consumers): MEDIT `.mesh` + post-processed `.msh`/`.vtu` with per-fault Physical Surface tags 101..(100+N). Same as the gmsh-replaced Phase 3.
- **`mesh_volume.cpp` CLI must remain backward-compatible** for the existing `mesh_volume MANIFEST OUTPUT_MEDIT_MESH [--lc-near …] …` invocation. New flags may be added.

### Dependency constraints
- CGAL 6.1.1 in `cgal-61` env. Specifically: `Mesh_3` (`Mesh_criteria_3`, `Polyhedral_complex_mesh_domain_3`, `Mesh_constant_domain_field_3`), `Conforming_constrained_Delaunay_triangulation_3` (CGAL 6.1+).
- Build via the existing `code_preprocess/corefine_cgal/CMakeLists.txt`. No new CMake deps needed.
- `pythonenv` for `meshio`-based post-processing / verification.

### Convention constraints
- Per `CLAUDE.md`: do not silently fall back to simpler approaches. If a phase fails, stop and report; do not auto-promote to the next phase.
- Per `mesh_volume.cpp:295-313` comments: keep the "subdomain 1 = rock; faults pair (1,1) embedded" subdomain layout. If a phase requires changing this, call it out explicitly.

### Numerical constraints
- **Hard floor**: every tet edge has length ≥ 100 m. This is an absolute requirement; no exceptions.
- **Soft target (preserved)**: `min q_tet ≥ 0.05`, `≥ 99.5 % q_tet ≥ 0.3`, `q_med ≥ 0.6`. Currently 99.94 / 0.91 — must not regress.
- **Surface conformality preserved**: every corefined polyline vertex must appear in the volume mesh. (Verifiable by sampling the manifest's `pairs[k].polylines` against `c3t3.triangulation().finite_vertices()`.)

## Survey of CGAL knobs (background)

| Knob | Header | Effect | Risk |
|---|---|---|---|
| `edge_min_size`        | `CGAL/Mesh_criteria_3.h:253-263`        | "Only feature edges with a length larger than this bound will be refined" | Doc says: "lower-bound may not be respected everywhere ... feature protection algorithm correctness is not guaranteed anymore" — soft floor, not hard |
| `facet_min_size`       | `CGAL/Mesh_criteria_3.h:281-287`        | Same semantics for facet Delaunay-ball radius                              | Same caveat |
| `cell_min_size`        | `CGAL/Mesh_criteria_3.h:303-310`        | Same for tet circumradius                                                  | Same caveat |
| `Mesh_constant_domain_field_3<GT, Index>` | `CGAL/Mesh_constant_domain_field_3.h:42` | Spatially-varying or per-subdomain sizing field; can be passed to `cell_size`/`facet_size` | Floor still soft (target, not hard) |
| `detect_features(angle)` / `detect_borders()` | `CGAL/Polyhedral_complex_mesh_domain_3.h:309/318` | Feature detection — already SKIPPED in current `mesh_volume.cpp:318-323`. Skipping cut sub-100 m tets from 1.7 M to 5,300. Cannot disable further. | n/a |
| `Conforming_constrained_Delaunay_triangulation_3` | `CGAL/Conforming_constrained_Delaunay_triangulation_3.h:557` (CGAL 6.1+) | Direct constrained Delaunay; no Mesh_3 protection-ball machinery. Outputs a constrained conforming triangulation honoring polyline constraints exactly. | Replaces Mesh_3 entirely. No Steiner refinement: bulk has whatever the minimum-edge property of the input gives. May need a separate refinement pass. |
| `c3t3.triangulation().tds().finite_edges()` + `Euler::collapse_edge` (post-processing) | `CGAL/Triangulation_data_structure_3.h` + custom code | Walk the c3t3, collapse sub-100 m tet edges, retri locally | Risks breaking patch-id labeling; must preserve fault surface conformality in post |

The intrinsic Mesh_3 protection ball mechanism is documented in `CGAL/Mesh_3/Protect_edges_sizing_field.h` and is **always on** when polyhedral surfaces meet at patch boundaries — it is **not** controlled by `detect_features`/`detect_borders`. That is the root cause of the 5,300 residuals.

---

## Phase 1: pass `edge_min_size`, `facet_min_size`, `cell_min_size` to the existing Mesh_3 criteria

### Goal
After this phase, `make_mesh_3` has been told the hard floor explicitly via three Mesh_criteria_3 named parameters, and the new mesh has < 100 sub-floor tets (vs current 5,300). One-line code change; ~25 minute meshing run; immediately measurable.

### Files to Create
- None.

### Files to Modify
- `code_preprocess/corefine_cgal/mesh_volume.cpp` — extend the `Criteria criteria(...)` call at lines 351-358.

### Detailed Requirements
1. Replace the existing criteria object construction (lines 351-358) with:
   ```cpp
   Criteria criteria(
       pp::edge_size(a.lc_near),               // 1500 m UPPER bound
       pp::edge_min_size(a.lc_min),            // 100 m LOWER bound  -- NEW
       pp::facet_angle(25.0),
       pp::facet_size(a.lc_near),              // 1500 m UPPER bound
       pp::facet_min_size(a.lc_min),           // 100 m LOWER bound  -- NEW
       pp::facet_distance(a.lc_near * 0.1),
       pp::cell_radius_edge_ratio(3.0),
       pp::cell_size(a.lc_near),               // 1500 m UPPER bound
       pp::cell_min_size(a.lc_min)             // 100 m LOWER bound  -- NEW
   );
   ```
   The named parameters `edge_min_size`, `facet_min_size`, `cell_min_size` are documented at `CGAL/Mesh_criteria_3.h:253-263`, `:281-287`, `:303-310` respectively. They take `FT` (= `K::FT` = `double` for our `Epick` kernel).
2. No CLI change. `a.lc_min` is already `100.0` by default and exposed as `--lc-min`.
3. Add one verbose log line just before `make_mesh_3`:
   ```cpp
   std::cerr << "  hard floors: edge_min=" << a.lc_min
             << " facet_min=" << a.lc_min
             << " cell_min=" << a.lc_min << "\n";
   ```

### Interfaces
- Public surface unchanged. Internal `Criteria` object now takes 9 named params instead of 6.

### Edge Cases to Handle
- **`lc_min > lc_near`**: nonsensical. Add an early `if (a.lc_min > a.lc_near) { std::cerr << "lc_min > lc_near\n"; return 2; }` guard right after `parse_args`.
- **`lc_min == 0`**: degenerate; treat as "feature unused". Same guard.
- **`edge_min_size` doc caveat** ("lower-bound may not be respected everywhere"): the residual count after Phase 1 may be > 0. Verify with the existing `check_msh_quality.py`-style script. If residuals exist, escalate to Phase 2.

### Acceptance Criteria
- [ ] Compiles in `cgal-61` env (`cmake --build build --target mesh_volume`).
- [ ] `mesh_volume data_corefined/manifest.json /tmp/m.mesh --verbose` exits 0 in ≤ 30 min.
- [ ] Loaded via `meshio`, the resulting `.mesh` has **`min(edge_length) ≥ 100 m`** AND **`min q_tet ≥ 0.05`** AND **`fraction(q_tet ≥ 0.3) ≥ 0.995`**.
- [ ] Surface conformality preserved: for each pair `(i,j)` in `manifest.pairs[]`, every shared polyline vertex (sampled from `pairs[k].polylines`) is present in the c3t3 triangulation within 1e-3 m.
- [ ] Tet count change vs current (4.8 M): not more than ±20 %. Significant deviation indicates Mesh_3 silently dropped patches or generated a bad mesh.
- [ ] q_med ≥ 0.85 (currently 0.91). Floor enforcement should not destroy bulk quality.

### Dependencies
- Depends on: nothing (the inputs in `data_corefined/` are ready).
- Required by: Phases 2-5 only if this phase fails.

### Risk
- **Low**. The named parameters are CGAL 6.1.1 stable API. Doc explicitly warns floor is soft, not hard; the only failure mode is residual sub-100 m tets, which are detectable and trigger escalation to Phase 2.

---

## Phase 2: spatially-varying sizing field that imposes 100 m everywhere except deep box interior

### Goal
After this phase, a `Mesh_constant_domain_field_3` (or a small custom `MeshDomainField_3` model) is supplied as the `cell_size` and `facet_size` argument. The field returns `lc_near = 1500` everywhere by default, but gracefully ramps to `lc_min = 100` if Mesh_3 ever requests a sub-100 m local size — making 100 m a deterministic local floor.

This phase only fires if Phase 1 leaves residuals.

### Files to Create
- `code_preprocess/corefine_cgal/min_size_field.h` — header-only template `Min_floor_field<GT, Index>` that wraps a base `Mesh_constant_domain_field_3` and clamps every queried size to `≥ size_floor`.

### Files to Modify
- `code_preprocess/corefine_cgal/mesh_volume.cpp` — instantiate `Min_floor_field` instead of the constant `a.lc_near` for `facet_size` and `cell_size` named params.

### Detailed Requirements
1. `min_size_field.h` defines:
   ```cpp
   template <class GT, class Index>
   class Min_floor_field {
   public:
       using FT      = typename GT::FT;
       using Point_3 = typename GT::Point_3;   // matches Mesh_constant_domain_field_3
       using Index_  = Index;

       Min_floor_field(FT size_target, FT size_floor)
         : target_(size_target), floor_(size_floor) {}

       FT operator()(const Point_3& /*p*/, int /*dim*/, const Index& /*idx*/) const {
           return std::max(target_, floor_);
       }
   private:
       FT target_, floor_;
   };
   ```
   (`Mesh_constant_domain_field_3` model contract: per CGAL doc the call operator takes `(Point_3, int dim, Index)` and returns `FT`. See `CGAL/Mesh_constant_domain_field_3.h:42-65`.)
2. In `mesh_volume.cpp`, replace the `pp::facet_size(a.lc_near)` and `pp::cell_size(a.lc_near)` arguments with:
   ```cpp
   Min_floor_field<K, MeshDomain::Index> sf(a.lc_near, a.lc_min);
   ...
   pp::facet_size(sf), pp::cell_size(sf), ...
   ```
   The `MeshDomain::Index` type is required because Mesh_3 calls the field with the patch index of the requested point.
3. Keep the Phase 1 `edge_min_size` / `facet_min_size` / `cell_min_size` floor parameters in place (they are still useful: even if the sizing field always returns ≥ 100 m, Mesh_3 may request a smaller size internally during protection ball insertion — the `*_min_size` params block that).

### Interfaces
- New header `min_size_field.h` exports one template class. No public API changes elsewhere.

### Edge Cases to Handle
- **Field instantiation type mismatch**: if `MeshDomain::Index` is not the type Mesh_3 expects, the code will fail to compile at the call site. Compile errors should be diagnosed by reading the Mesh_3 named-parameter source path: `CGAL/Mesh_3/parameters.h` → `pp::facet_size(...)`.
- **Sub-100 m residuals still exist**: this phase shifts the failure mode from "Mesh_3 ignored my upper-bound criterion" to "Mesh_3 ignored my sizing field." If residuals persist, the protection-ball mechanism is overriding the field; escalate to Phase 3.

### Acceptance Criteria
- [ ] Same as Phase 1, with reduced sub-floor count vs Phase 1's residual.
- [ ] q_med ≥ 0.85 (no bulk regression).

### Dependencies
- Depends on: Phase 1 (re-uses the `lc_min` arg and the `*_min_size` named-params).
- Required by: Phase 3 only if residuals remain.

### Risk
- **Medium**. The sizing field model is straightforward, but the Mesh_3 named-param plumbing requires the right `Index` type. Compile errors are likely on the first attempt; debug by checking `CGAL/Mesh_3/parameters.h` and the `Mesh_constant_domain_field_3` instantiation pattern shown in the CGAL docs.

---

## Phase 3: post-processing edge-collapse pass on the c3t3 triangulation

### Goal
After `make_mesh_3` returns, walk the c3t3, identify every tet edge with length < 100 m, and collapse it via `CGAL::Euler::collapse_edge` (or, since c3t3's underlying TDS is a `Triangulation_3` not a halfedge graph, via direct vertex relocation: `tr.move_if_no_collision(v_short, midpoint)` followed by sliver re-evaluation).

Conformality preservation: skip any edge whose collapse would move a vertex that lies on a fault surface (subdomain-2-incident vertex). This is the analogue of the Phase 2 corefine pipeline's `polyline_cleanup::collapse_short_polyline_edges_pairwise` we already shipped, applied to the volume mesh.

### Files to Create
- `code_preprocess/corefine_cgal/c3t3_min_edge_cleanup.h` — header-only template:
   ```cpp
   template <class C3t3>
   struct C3t3CleanupStats { std::size_t n_collapsed; std::size_t n_residual; double residual_min; };

   template <class C3t3>
   C3t3CleanupStats collapse_short_volume_edges(C3t3& c3t3, double min_edge);
   ```

### Files to Modify
- `code_preprocess/corefine_cgal/mesh_volume.cpp` — call `collapse_short_volume_edges(c3t3, a.lc_min)` immediately after `make_mesh_3`, before `c3t3.output_to_medit`.

### Detailed Requirements
1. Iterate all finite edges of `c3t3.triangulation()`. For each edge `e = (v1, v2)` with squared length `< (lc_min)^2`:
   1. Skip if `c3t3.in_dimension(v1) <= 2` or `c3t3.in_dimension(v2) <= 2` (the vertex is on a feature edge or on a surface — moving it would break conformality).
   2. Compute midpoint `m = 0.5 * (v1->point() + v2->point())`.
   3. Try `tr.move_if_no_collision(v1, m)`; if it returned a different vertex (a collision merged two), update bookkeeping. If no move possible, skip.
   4. Repeat until convergence or `max_iter = 50`.
2. Surface conformality assertion before returning: every vertex `v` with `c3t3.in_dimension(v) == 2` must have its position UNCHANGED by this pass. This is automatic given the dim guard in step 1.1, but assert it post-hoc.
3. Report `n_collapsed`, `n_residual` (sub-100 m edges that could not be collapsed), and `residual_min`.

### Interfaces
- One header, one templated function. No CMake change.

### Edge Cases to Handle
- **Vertex `v1` is a corner** (c3t3.in_dimension(v) == 0): treat like surface vertex — skip.
- **Collapse would invert an incident tet**: `move_if_no_collision` already guards. If it returns the original handle, we know the move failed; treat as residual.
- **A non-collapsible sub-100 m edge persists**: report and hard-fail per Phase 2 corefine convention (the pipeline must not silently ship a contract violation).

### Acceptance Criteria
- [ ] Same hard-floor and quality criteria as Phase 1.
- [ ] Surface conformality assertion passes (no surface vertex moved).
- [ ] Sub-100 m edge count after pass: 0.
- [ ] q_med drop ≤ 0.05 absolute (cleanup may slightly degrade quality near merges; bound the regression).

### Dependencies
- Depends on: Phase 1 (or Phase 2 if implemented). Phase 3 is standalone-applicable as well.
- Required by: nothing.

### Risk
- **Medium-high**. Volume edge collapse is more delicate than the surface case we already handled in Phase 2 corefine. `move_if_no_collision` doesn't expose an explicit "would this collapse the link" check the way `Euler::collapse_edge`+`does_satisfy_link_condition` do for halfedge graphs. May leave residuals at multi-tet bottlenecks.

---

## Phase 4: replace `Polyhedral_complex_mesh_domain_3` with `Conforming_constrained_Delaunay_triangulation_3` (CGAL 6.1+)

### Goal
Build a CCDT-3 of the box + 6 fault polyhedra. CCDT-3 produces a constrained conforming Delaunay triangulation that honors the input polylines exactly — without inserting protection-ball Steiner points — then optionally refines for quality. This bypasses the entire Mesh_3 protection mechanism that causes the 5,300 residuals.

### Files to Create
- `code_preprocess/corefine_cgal/mesh_volume_ccdt.cpp` — sibling driver to `mesh_volume.cpp` using the new API. Same CLI; outputs a `.mesh` file consumable by the existing Python post-processor.

### Files to Modify
- `code_preprocess/corefine_cgal/CMakeLists.txt` — add `add_executable(mesh_volume_ccdt mesh_volume_ccdt.cpp)` and the existing `foreach(t ... )` link block.

### Detailed Requirements
1. Build a single combined polygon soup containing:
   - 12 triangles from the bounding box (closed).
   - All triangles from each of the 6 corefined faults.
   Use `PMP::polygon_soup_to_polygon_mesh` to build a single `CGAL::Surface_mesh` with patch IDs assigned per source polyhedron (box=0, faults 1..6).
2. Call `CGAL::make_conforming_constrained_Delaunay_triangulation_3<CDT>(combined_mesh)` — the CCDT class is at `CGAL/Conforming_constrained_Delaunay_triangulation_3.h:557`.
3. CCDT-3 by itself produces a triangulation with NO refinement — just enough tets to honor the constraints. To get bulk refinement to LC_NEAR = 1500 m, follow with a refinement pass:
   - The CGAL 6.1.1 docs describe a `refine_constrained_Delaunay_triangulation_3()` follow-up that takes Mesh_3-style criteria.
   - Investigate whether `cell_min_size` / `facet_min_size` are honored in this refinement (they are NOT subject to the protection-ball machinery, so the floor should be hard).
4. Output to MEDIT via `c3t3.output_to_medit` or its CCDT-3 analogue (TBD; verify the API exists in 6.1.1).

### Interfaces
- New `mesh_volume_ccdt` binary parallels `mesh_volume`.
- Same MEDIT output format → existing post-processor (`check_msh_quality.py`) consumes it unchanged.

### Edge Cases to Handle
- **CCDT-3 with 7 polyhedra in one soup**: needs to handle the 7 disconnected components correctly. Per the CGAL CCDT-3 docs, the input is a "polygon soup" — so it should be fine, but verify with a 1-fault smoke test first.
- **Refinement pass with `cell_min_size` may behave differently from Mesh_3**: the CCDT-3 refinement uses different machinery; the named-parameter set may differ. Read `CGAL/Conforming_constrained_Delaunay_triangulation_3.h:557-1100` and the shipped header docs to confirm.
- **Subdomain labeling**: CCDT-3 likely doesn't have the (positive, negative) subdomain pair semantics of Polyhedral_complex_mesh_domain_3. Need a separate post-pass to label tets as subdomain 1 (rock) using inside/outside testing against the box.

### Acceptance Criteria
- [ ] Same hard-floor and quality criteria as Phase 1.
- [ ] No sub-100 m tet edges in the output.
- [ ] Surface conformality preserved (every fault triangle from the input STLs appears in the c3t3 surface facets).
- [ ] All tets labeled subdomain 1; output `.mesh` file has same MEDIT structure as Phase 1's output.

### Dependencies
- Depends on: nothing (parallel implementation).
- Required by: nothing.

### Risk
- **High**. CCDT-3 is a CGAL 6.1 NEW class; documentation and examples in the wild are scarce. The refinement-pass + subdomain-labeling story is unclear until prototyped. Plan to spend a half-day investigating before committing.

---

## Phase 5: drop the bottom 100 m of fault perimeter (geometry edit, last resort)

### Goal
**Out of scope** without explicit user approval per the Phase 0 plan's geological-data-immutable constraint.

### Notes
- The Phase 2 corefine pipeline's `cross_polyline_snap` already moved 3 fault perimeter vertices by up to 50 m to satisfy edge ≥ 100 m on the surface. Further geometry editing (e.g. trimming the bottom 100 m of every fault to round off perimeter kinks) might eliminate the protection-ball oversteer at the source. But this changes the geological model.
- Listed here only for completeness; do NOT implement.

---

## Recommended Path

**Implement Phase 1 first.** It is a one-line change, exercises a stable CGAL 6.1.1 API explicitly designed for this purpose, and the failure mode (residuals) is detectable and bounded.

If Phase 1 leaves > 0 residual sub-100 m tets:
1. Try Phase 2 (sizing-field wrapper). Higher confidence than Phase 3 because it stays in the documented Mesh_3 pipeline.
2. If Phase 2 also leaves residuals, run Phase 3 (post-collapse). This is the "always works" cleanup but with the highest chance of bulk-quality regression.
3. Only consider Phase 4 (CCDT-3 rewrite) if Phases 1-3 collectively cannot drive residuals to 0 — that would indicate a fundamental Mesh_3 limitation worth rearchitecting around.

Phases 1, 2, and 3 are stackable: Phase 2's sizing field + Phase 1's `*_min_size` floors + Phase 3's post-collapse together form a defense-in-depth approach if any single phase is insufficient.

## Testing Strategy

### Per-phase smoke check (after each phase)
1. Build via `cmake --build build --target mesh_volume[_ccdt]`.
2. Run on the 6-fault canonical fixture: `mesh_volume data_corefined/manifest.json /tmp/m.mesh --verbose`.
3. Run the existing `check_msh_quality.py /tmp/m.mesh` (already exists in `code_preprocess/`); confirm all three Phase 3.5b acceptance lines pass.
4. Visualize in ParaView via `meshio convert /tmp/m.mesh /tmp/m.vtu`; eyeball the slivers near the (4,5) and (3,5) polylines.

### Conformality assertion (independent of phase)
Add to `check_msh_quality.py` a new check: for every coord in `manifest.pairs[k].polylines`, find the nearest c3t3 vertex; assert distance < 1e-3 m. Phase 2 corefine guarantees these polyline coords exist in the input STLs at full precision; Mesh_3 must preserve them.

### Reference / regression
The current 4.8 M-tet output (from `mesh_volume.cpp` HEAD) is the regression baseline. Any new mesh that has q_med < 0.85 OR conformality assertion failures has regressed and must be rejected even if it satisfies the edge floor.

## Risk Assessment

- **`edge_min_size` doc caveat is real**: "lower-bound may not be respected everywhere" means Phase 1 alone may leave a small residual. Plan accommodates this by Phase 2/3 escalation.
- **Mesh_3's protection-ball mechanism is undocumented in detail**: the "1.7 M sub-100 m tets with detect_features ON" → "5,300 with it OFF" delta is empirical. Whether `cell_min_size` overrides protection-ball insertion is not documented; we'll learn from Phase 1's run.
- **Volume edge collapse (Phase 3) on a Mesh_3 c3t3 may fail at multi-tet bottlenecks**: `move_if_no_collision` is conservative. Worst case: Phase 3 leaves residuals where Phase 1+2 also did, and we have to commit to Phase 4 anyway.
- **CCDT-3 is new**: unknown stability on 4 M-tet domains; unknown refinement-pass behavior. High exploration cost.
- **The single biggest risk is silent quality regression**: a mesh with edge_min ≥ 100 m but q_med = 0.5 is worse than the current mesh. The "no bulk regression" acceptance criterion (q_med ≥ 0.85) is the firewall.

## Open Questions for the User

1. **Is the "0 sub-100 m edges" target absolute, or is "< 100 sub-100 m edges out of 4.8 M (0.002 %)" acceptable?** The plan currently treats 0 as the bar. If even single-digit residuals are tolerable, Phase 1 will likely suffice and Phase 3 is overkill.
2. **Wall-clock budget**: each meshing run is ~25 minutes. Phase 1 + verification = 1 hour. Full Phase 1+2+3 escalation: 4-6 hours. Phase 4 (CCDT-3 from scratch): 1-2 days. Confirm time horizon before starting Phase 4.
3. **Fault interface attribute requirements**: Phase 4's subdomain-labeling story requires knowing whether MFEM needs each fault as a separate Physical Surface (current: tags 101..106) or whether a single "fault" tag is acceptable. This affects the post-labeling step in Phase 4 only.
