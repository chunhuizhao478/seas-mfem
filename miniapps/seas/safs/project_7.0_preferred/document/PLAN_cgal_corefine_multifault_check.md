# Plan Review: PLAN_cgal_corefine_multifault.md (2026-05-08)

## Review Scope
- Plan: `/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/safs/project_7.0_preferred/document/PLAN_cgal_corefine_multifault.md`
- Reference docs consulted:
  - `/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/safs/project_7.0_alternative/code_preprocess/EXPLORE_cgal_corefine.md`
  - `/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/safs/project_7.0_alternative/code_preprocess/EXPLORE_pymeshlab_intersections.md`
  - `/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/safs/project_7.0_alternative/code_preprocess/PLAN_corefine_pipeline.md`
  - `/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/safs/project_7.0_preferred/code_preprocess/clean_freesurface_mesh.py`
  - `/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/safs/project_7.0_alternative/code_meshing/safs_fault_box_freesurface_clip.geo`
- CGAL 6.1.1 headers verified in `/Users/chunhuizhao/miniforge/envs/cgal-61/include/CGAL/Polygon_mesh_processing/`:
  - `tangential_relaxation.h:87,96,106` — confirms `edge_is_constrained_map` and `vertex_is_constrained_map` named parameters exist (and are `false` by default)
  - `surface_Delaunay_remeshing.h:94,101,180,219` — confirms `polyline_constraints` / `protect_constraints` named parameters; default `protect_constraints=false`
  - `intersection.h:1625,1664,1674` — confirms `intersecting_meshes(range, out, np, ...)` overloads

## Findings

### [R-001] [CRITICAL] [PLAN] [§Phase 1, `quality_repair::run`] — `tangential_relaxation` final pass will silently move polyline vertices and break pairwise conformality

**Category:** BUG (scope creep, missing constraint propagation)

**Description:**
Phase 1's algorithm step 9 specifies that `quality_repair::run` "iterates `PMP::remove_almost_degenerate_faces` until convergence or `max-iterations`, then runs `PMP::tangential_relaxation` once for final smoothing." The function signature in `quality_repair.h` takes only `min_edge_size`, `q_floor`, `cap_threshold_cos`, `needle_threshold_ratio`, `max_iterations` (inherited from the alternative plan PLAN_corefine_pipeline.md:L218-L227) and a single `ECM` (edge constraint map). It does **not** take a vertex constraint map and does not specify how `tangential_relaxation` is called.

By default `PMP::tangential_relaxation`'s `vertex_is_constrained_map` is the all-`false` map (verified at `tangential_relaxation.h:96,106`). When called on `A_out` after `surface_Delaunay_remeshing`, **every polyline vertex (including the shared intersection vertices in `B_out`) becomes free to drift along its tangent plane**. Polyline vertices in `A_out` and `B_out` will then move independently and pairwise conformality is destroyed — this is precisely the failure mode the whole pipeline exists to prevent.

This is also scope creep relative to the source plan. The alternative plan's `quality_repair::run` (PLAN_corefine_pipeline.md:L208-L227) ends after `remove_almost_degenerate_faces`; it does NOT call `tangential_relaxation`. The new plan added the final smoothing pass without specifying the constraint propagation needed to keep it safe.

**Trigger:**
Any pair `(A, B)` with a non-trivial intersection polyline, after the standard pipeline runs Phase 1 step 9 with the default parameters.

**Actual behavior:**
`tangential_relaxation` moves every polyline vertex onto the local 1-ring centroid (projected onto the tangent plane). Since `A_out` and `B_out` have different 1-rings around shared polyline vertices, the relaxation moves their copies to different positions. Phase 1 step 10 ("conformality check") then fails with exit 3, OR — worse — the writer's tolerance allows the drift through and gmsh dedup silently merges nearby-but-not-coincident vertices and produces ill-shaped tets along the intersection.

**Expected behavior:**
Either (a) drop the `tangential_relaxation` pass entirely (the alternative plan does without it), or (b) require `quality_repair::run` to construct a `vertex_is_constrained_map` that marks every polyline vertex constrained, and pass it to `tangential_relaxation`'s named parameter pack along with the existing `edge_is_constrained_map`.

**Suggested fix:**
Apply (a). It is the simpler, safer path and matches the source plan. Replace the second sentence of Phase 1 algorithm step 9 with:

```diff
-9. **Quality repair** on each remeshed mesh, with the resampled
-   polyline marked as constrained (re-derive the constraint map by
-   matching vertex coordinates against `resampled`; tolerance 1e-9):
-   ```cpp
-   auto ecmA2 = mark_polyline_as_constrained(A_remeshed, resampled);
-   auto stA = quality_repair::run(A_remeshed, ecmA2, args.min_edge);
-   // same for B
-   ```
+9. **Quality repair** on each remeshed mesh, with the resampled
+   polyline marked as constrained (re-derive the constraint map by
+   matching vertex coordinates against `resampled`; tolerance 1e-9):
+   ```cpp
+   auto ecmA2 = mark_polyline_as_constrained(A_remeshed, resampled);
+   auto stA = quality_repair::run(A_remeshed, ecmA2, args.min_edge);
+   // same for B
+   ```
+   `quality_repair::run` performs `PMP::remove_almost_degenerate_faces`
+   only (matching the alternative plan PLAN_corefine_pipeline.md:L208-L227);
+   it does NOT call `PMP::tangential_relaxation`.  Tangential relaxation
+   would move polyline vertices unless an explicit
+   `vertex_is_constrained_map` is provided, and the marginal smoothing
+   gain does not justify the conformality risk.
```

Also delete the parenthetical in Phase 1 §Detailed Requirements that introduces the smoothing pass:

```diff
-9. `auto stA = quality_repair::run(A_out, ecm_out_A, args.min_edge, /*q_floor=*/0.3, /*cap_threshold_cos=*/cos(160°), /*needle_threshold=*/4.0, args.max_iterations);` same for B. The function iterates `PMP::remove_almost_degenerate_faces` until convergence or `max-iterations`, then runs `PMP::tangential_relaxation` once for final smoothing.
+9. `auto stA = quality_repair::run(A_out, ecm_out_A, args.min_edge, /*q_floor=*/0.3, /*cap_threshold_cos=*/cos(160°), /*needle_threshold=*/4.0, args.max_iterations);` same for B. The function iterates `PMP::remove_almost_degenerate_faces` until convergence or `max-iterations`. It does NOT call `PMP::tangential_relaxation`: that function would move polyline vertices unless an explicit `vertex_is_constrained_map` is provided, which would break pairwise conformality. (See finding R-001.)
```

**Test case:**
```python
def test_R001_polyline_vertices_unchanged_by_quality_repair():
    """After quality_repair::run, polyline vertices in A_out and B_out
    must remain at the same coordinates as in `resampled`."""
    a_out, b_out, resampled = run_corefine_pair_smoke(
        in_a="synth_rect_xy.off", in_b="synth_rect_xz.off",
        polyline_spacing=2.0, min_edge=0.5, mesh_edge_size=2.0)
    for v_resampled in flatten(resampled):
        assert vertex_exists(a_out, v_resampled, tol=1e-9), \
            f"polyline vertex {v_resampled} drifted in A_out"
        assert vertex_exists(b_out, v_resampled, tol=1e-9), \
            f"polyline vertex {v_resampled} drifted in B_out"
```

---

### [R-002] [CRITICAL] [PLAN] [§Phase 2, `corefine_faults.py` step 5] — `corefine_set` re-entrancy bug: re-running picks up its own previous outputs as inputs

**Category:** BUG (logic error: same dir for input and output)

**Description:**
Phase 2 step 5 invokes `corefine_set` with **the same path** for input and output:
```python
subprocess.run([str(cgal_bin), str(workdir), str(workdir), "--ext", ".off",
                ...])
```
Phase 1 specifies that `corefine_set` "scans `IN_DIR` for files matching `*<ext>` and treats each as an input mesh" — the glob is `workdir/*.off`. After the first run, `workdir/` contains both `<basename>_2000m.off` (the input) **and** `<basename>_2000m_corefined.off` (the output). On the next run, the glob picks up *both* sets, and `corefine_set` will (a) try to corefine 12 meshes (6 inputs × 2 generations) instead of 6, (b) produce nonsensical pairwise conformality between an input and its own corefined version, and (c) write `<basename>_corefined_corefined.off`.

The plan does say "Re-running with the same args is idempotent (same number of vertices/faces in each output to within 1 %)" as an acceptance criterion in Phase 2, but the spec as written cannot satisfy that criterion.

**Trigger:**
Running `python code_preprocess/corefine_faults.py --res 2000` twice without `--keep-intermediate` removed in between, OR running once with `--keep-intermediate` and re-running. Also triggered if Phase 5's `test_corefine_set_smoke_2000m` runs more than once in the same workdir.

**Actual behavior:**
Second invocation reads 12 OFF files (6 inputs + 6 previous-run outputs), corefines all C(12,2)=66 candidate pairs, writes 12 `_corefined.off` files. Manifest reports nonsense.

**Expected behavior:**
Re-runs are idempotent: input glob excludes previous outputs, OR input and output are in separate directories.

**Suggested fix:**
Use separate input and output directories under `workdir`. Apply both edits:

```diff
   [`--workdir DIR`]              (default: ../work/corefine_2000m)
+  Two subdirectories are created beneath `--workdir`:
+    in/   — OFF intermediates produced by the .ts/.stl → .off converter (step 2)
+    out/  — `*_corefined.off` produced by `corefine_set` (step 5)
```

```diff
-4. Verify `cgal-bin` exists; if not, print:
+4. Verify `cgal-bin` exists; if not, print:
    ...
-5. `subprocess.run([str(cgal_bin), str(workdir), str(workdir), "--ext", ".off",
+5. `subprocess.run([str(cgal_bin), str(workdir / "in"), str(workdir / "out"), "--ext", ".off",
    "--mesh-edge-size", str(mesh_edge_size), "--min-edge", str(min_edge),
    "--polyline-spacing", str(polyline_spacing),
-   "--manifest", str(workdir / "manifest.json")] + (["--verbose"] if verbose else []),
+   "--manifest", str(workdir / "out" / "manifest.json")] + (["--verbose"] if verbose else []),
    check=True, capture_output=False)`.
```

And update step 3 ("Convert each `.stl` to `.off` in `workdir`") to write into `workdir / "in"`, and step 6 ("Convert each `*_corefined.off` in `workdir`") to read from `workdir / "out"`.

**Test case:**
```python
def test_R002_idempotent_under_repeated_runs(tmp_path):
    """Two consecutive invocations on the same fixture must yield the
    same number of output STLs and the same vertex/face counts ±1%."""
    out_dir = tmp_path / "data_corefined"
    work_dir = tmp_path / "work"
    fixture = "data_cleanfreesurf"  # 6 STLs at 2000m
    run_corefine_faults(fixture, out_dir, work_dir, res=2000)
    counts1 = read_manifest_counts(out_dir / "manifest.json")
    run_corefine_faults(fixture, out_dir, work_dir, res=2000)
    counts2 = read_manifest_counts(out_dir / "manifest.json")
    assert len(counts1) == 6 and len(counts2) == 6, \
        "Manifest must contain exactly 6 meshes on both runs"
    for n1, n2 in zip(counts1, counts2):
        assert abs(n1 - n2) / max(n1, 1) < 0.01, \
            f"Re-run vertex/face count drifted by >1%: {n1} → {n2}"
```

---

### [R-003] [CRITICAL] [PLAN] [§Phase 2 CLI defaults] — `--mesh-edge-size` default contradicts Phase 1 default and the `.geo` `LC_NEAR`, producing an end-to-end size mismatch

**Category:** DEVIATION (internal inconsistency in defaults)

**Description:**
- Phase 1 (C++ binaries) default `--mesh-edge-size = 1500.0`.
- Phase 2 (Python driver) default `--mesh-edge-size = --res (e.g. 2000)`.
- Phase 3 generated `.geo` template has `LC_NEAR = 1500.0` (hardcoded).
- Phase 4 sensitivity sweep `baseline` config: `mesh-edge-size = 1500`.

Running the documented end-to-end pipeline with default flags (`python code_preprocess/corefine_faults.py --res 2000`) hands `--mesh-edge-size 2000` to `corefine_set`, which produces fault triangles up to ~2000 m. The generated `.geo` then requests `LC_NEAR = 1500` near the fault. Because gmsh embeds the discrete STL surfaces *as-is*, the actual local mesh size is dictated by the embedded triangles (~2000 m), not by `LC_NEAR` (1500 m). The user expects 1500 m near-fault tets and gets 2000 m+ tets.

This also makes Phase 4's `baseline` config differ from "what running the pipeline with defaults produces" — the sensitivity sweep is supposed to *defend* the defaults but the defaults are not aligned with the baseline config it claims to defend.

**Trigger:**
`python code_preprocess/corefine_faults.py --res 2000` with no other flags. Also triggers when comparing Phase 4 `baseline` results against Phase 5 `test_per_fault_quality` results (the latter uses the Python driver default).

**Actual behavior:**
Three different effective `mesh-edge-size` defaults across three phases (1500 / 2000 / 1500). End-to-end run produces 2000 m output that the `.geo` then mis-describes as 1500 m near-fault.

**Expected behavior:**
A single canonical default for `mesh-edge-size` across Phases 1, 2, 3, and 4 baseline. Recommend `1500.0`.

**Suggested fix:**
Align Phase 2's default to `1500.0` (matching Phase 1 and Phase 3):

```diff
-| `--mesh-edge-size` | float, m | `--res` (e.g. 2000) | passed through to `corefine_set` |
+| `--mesh-edge-size` | float, m | 1500.0 | passed through to `corefine_set` (matches Phase 1 default and Phase 3 `LC_NEAR`) |
```

And clarify in Phase 2 prose that `--res` selects only the **input fixture resolution** (which `_<R>m_clean_clip.stl` files to read), not the **output remesh target** — these are independent knobs.

```diff
-| `--res` | int | 2000 | resolution suffix to scan; only 2000 is validated by this plan |
+| `--res` | int | 2000 | input fixture resolution to scan (`*_<R>m_clean_clip.stl`); only 2000 is validated by this plan. Independent of `--mesh-edge-size`. |
```

**Test case:**
```python
def test_R003_default_mesh_edge_size_consistent():
    """The default --mesh-edge-size from the Python driver must equal
    the default from corefine_pair, and must equal LC_NEAR in the .geo."""
    py_default = get_default_arg("corefine_faults.py", "--mesh-edge-size")
    cpp_default = get_default_arg("corefine_pair", "--mesh-edge-size")
    geo_lc_near = parse_geo_constant(
        "code_meshing/safs_multifault_box_2000m.geo", "LC_NEAR")
    assert py_default == cpp_default == geo_lc_near, \
        f"Inconsistent defaults: py={py_default}, cpp={cpp_default}, geo={geo_lc_near}"
```

---

### [R-004] [CRITICAL] [PLAN] [§Phase 1 `io_helpers.h`] — OFF writer precision unspecified; CGAL default writer (`<<` with default ostream precision) loses sub-mm conformality at UTM scale

**Category:** EDGE_CASE / ASSUMPTION (precision contract incomplete)

**Description:**
Phase 1 lists three I/O helpers:
- `read_polygon_mesh_any(path, mesh)` — reads any format
- `write_off(path, mesh)` — undefined precision
- `write_stl_high_precision(path, mesh, "%.10g")` — explicit `%.10g`

OFF is the **inter-stage** format used between `corefine_set`'s in-place writes and the Python driver's STL conversion. CGAL's default OFF writer uses `operator<<` on a `std::ostream` with the default precision (6 significant digits — verified by inspection of the standard `std::ios::precision()` default). For a UTM Y coordinate of `3937744.123456789`, six-digit output is `3.93774e+06`, precision ~1 m. Polyline conformality, which requires sub-mm match, is **silently destroyed at the OFF write step** before the high-precision STL writer ever runs.

The plan's "STL ASCII round-trip preserves polyline conformality" acceptance criterion (Phase 1, criterion 4) tests the STL writer but not the OFF writer; the test would pass on a fresh in-memory pipeline yet still fail in the actual end-to-end flow because intermediate OFF I/O has already discarded precision.

**Trigger:**
Any end-to-end run that writes corefined meshes to OFF and reads them back (i.e., the documented Phase 2 algorithm: corefine_set OFF → meshio reads OFF → meshio writes STL). Most apparent in pairwise conformality check after STL conversion.

**Actual behavior:**
Default-precision OFF writer rounds to 6 significant digits → ~1 m precision at UTM Y. Polyline vertex coordinates differ between A and B by O(1 m) after the OFF round-trip. Phase 4's `Geometry.Tolerance = 1e-3` (1 mm) is now far smaller than the actual coordinate drift; gmsh sees A and B's polylines as different vertex sets and either fails to dedup (non-conformal merge → ill-shaped tets) or the size-field `LC_NEAR=1500` triangles wrap each near-but-distinct vertex pair separately.

**Expected behavior:**
`write_off` matches or exceeds the precision contract of `write_stl_high_precision`. Both should write at least 12 significant digits to keep conformality safely below the 1 mm `Geometry.Tolerance`.

**Suggested fix:**
Specify `write_off` explicitly in Phase 1's "Files to Create" list:

```diff
-- `code_preprocess/corefine_cgal/io_helpers.h` — wrappers for
-  `CGAL::IO::read_polygon_mesh` / `write_polygon_mesh` with explicit
-  `OFF`/`STL` format detection by extension and a high-precision ASCII
-  STL writer (`%.10g`).
+- `code_preprocess/corefine_cgal/io_helpers.h` — wrappers for
+  `CGAL::IO::read_polygon_mesh` / `write_polygon_mesh` with explicit
+  `OFF`/`STL` format detection by extension, **high-precision ASCII
+  OFF writer (`%.15g`), and high-precision ASCII STL writer (`%.15g`)**.
+  Both writers must call `out.precision(15)` (or the `printf` equivalent)
+  before emitting any vertex coordinate.  Default `operator<<` precision
+  (6 digits) is **insufficient** for UTM coordinates ~10⁶ m and silently
+  destroys polyline conformality at the inter-stage OFF round-trip.
```

Add a matching acceptance criterion in Phase 1:

```diff
 - [ ] STL ASCII round-trip preserves polyline conformality: writing
       output as ASCII STL with `%.10g`, reading back via `meshio`,
       every polyline vertex coordinate diff between A and B < 1e-6.
+- [ ] **OFF round-trip preserves polyline conformality**: writing
+      output as ASCII OFF with `%.15g`, reading back via `meshio`,
+      every polyline vertex coordinate diff between A and B < 1e-9 m
+      (absolute), independent of bbox magnitude.
```

**Test case:**
```python
def test_R004_off_writer_precision_at_utm_scale():
    """A vertex at UTM coords (5.4e5, 3.9e6, -1.2e4) round-tripped
    through write_off / read_off must round-trip to within 1e-9 m."""
    p_in = (541234.567890123, 3937744.123456789, -12345.678901234)
    p_out = roundtrip_off("/tmp/r004.off", p_in)
    for c_in, c_out in zip(p_in, p_out):
        assert abs(c_in - c_out) < 1e-9, \
            f"OFF writer lost precision: {c_in} → {c_out}"
```

---

### [R-005] [MODERATE] [PLAN] [§Phase 1 STL writer / §Phase 3 `Geometry.Tolerance`] — `%.10g` writer precision is at-the-edge for UTM Y-coordinates; tolerance/precision margin is 0× to 10×

**Category:** ASSUMPTION (numerical margin too tight)

**Description:**
The plan specifies the STL writer at `%.10g` (Phase 1 `io_helpers.h`) and `Geometry.Tolerance = 1e-3` (1 mm) in Phase 3's `.geo` template. `%.10g` formats with 10 significant digits. For a UTM Y coordinate `3937744.123456789`, `%.10g` produces `3937744.123` — last printed digit at the 10⁻³ m position, i.e. **1 mm precision**. Tolerance is **also** 1 mm. The two are equal, so a polyline vertex coordinate that rounds in one direction in mesh A and in the opposite direction in mesh B can drift by 1 mm — exactly at the `Geometry.Tolerance` threshold. Whether gmsh dedups or splits is then non-deterministic.

The plan's own §Defence section claims "tolerance > writer precision" but the arithmetic above shows they are equal at the worst-case Y coordinate. The alternative-plan EXPLORE_cgal_corefine.md:L426-L431 explicitly recommends "≥10 significant decimal digits" — `%.10g` meets that letter but not the spirit when STL coordinates approach 10⁷.

**Trigger:**
Any polyline vertex with Y ≥ 10⁶ m (which is every vertex in the SAFS UTM data; Y ranges 3.69 × 10⁶ to 3.99 × 10⁶ per EXPLORE_pymeshlab_intersections.md:L62-L69) where the unrounded coordinate's 11th significant digit is ≥ 5 in mesh A and < 5 in mesh B (or vice versa).

**Actual behavior:**
A's and B's polyline coordinate copies can differ by 1 mm in the LSB, causing inconsistent dedup behavior in gmsh.

**Expected behavior:**
Writer precision strictly less than `Geometry.Tolerance`. Recommend `%.15g` (~10⁻¹⁰ m at UTM Y) and tolerance `1e-3` so the margin is 10⁷×.

**Suggested fix:**
Tighten the writer precision in Phase 1:

```diff
-- `code_preprocess/corefine_cgal/io_helpers.h` — wrappers for
-  `CGAL::IO::read_polygon_mesh` / `write_polygon_mesh` with explicit
-  `OFF`/`STL` format detection by extension and a high-precision ASCII
-  STL writer (`%.10g`).
+- `code_preprocess/corefine_cgal/io_helpers.h` — wrappers for
+  `CGAL::IO::read_polygon_mesh` / `write_polygon_mesh` with explicit
+  `OFF`/`STL` format detection by extension and a high-precision ASCII
+  STL writer (`%.15g` — chosen so that the writer's LSB at UTM Y
+  ~ 4×10⁶ m is ~10⁻⁹ m, four orders of magnitude tighter than the
+  gmsh `Geometry.Tolerance = 1e-3` of Phase 3).
```

(R-004 also touches this section; the two fixes compose: both writers should be `%.15g`.)

Update §Convention constraints and §Numerical constraints accordingly:

```diff
-- UTM coordinates ~10⁶ m. ASCII I/O must use ≥ 10 significant decimal digits (`%.10g` or 15-digit `printf`); STL conformality is lost if the writer rounds to 6 digits.
+- UTM coordinates ~10⁶ m (Y up to 4×10⁶). ASCII I/O **must** use ≥ 15 significant decimal digits (`%.15g`); 10 digits leaves the LSB at the same magnitude as `Geometry.Tolerance` and the dedup behavior becomes non-deterministic.
```

**Test case:**
```python
def test_R005_writer_lsb_below_geometry_tolerance():
    """For a UTM-scale vertex, two writes that differ in the unrounded
    11th decimal (one rounds up, one rounds down) must produce printed
    coords whose absolute difference is < 0.1 × Geometry.Tolerance."""
    geometry_tolerance_m = 1e-3
    # Vertex magnitudes that bracket the LSB at the writer's precision.
    p_a = (541234.5678901234567, 3937744.1234567894, -12345.6789)
    p_b = (541234.5678901234999, 3937744.1234567899, -12345.6789)
    s_a = format_stl_vertex(p_a)
    s_b = format_stl_vertex(p_b)
    diff = max_coord_diff(parse_vertex(s_a), parse_vertex(s_b))
    assert diff < 0.1 * geometry_tolerance_m, \
        f"Writer precision {diff} m is too coarse for tolerance {geometry_tolerance_m} m"
```

---

### [R-006] [MODERATE] [PLAN] [§Defence section] — `protect_constraints=true` justification mis-cites a precondition that belongs to a different CGAL function (`isotropic_remeshing`, not `surface_Delaunay_remeshing`)

**Category:** QUALITY (incorrect technical rationale; risks future regression)

**Description:**
The §Defence section ("Why `protect_constraints=true` in `surface_Delaunay_remeshing` rather than `false`?") states:

> Because the polyline is *already resampled* to satisfy the `≤ 4/3·mesh_edge_size` precondition; `true` then preserves the polyline exactly...

The "≤ 4/3·target_edge_length" precondition documented at EXPLORE_cgal_corefine.md:L142-L145 belongs to `PMP::isotropic_remeshing`, not `surface_Delaunay_remeshing`. CGAL 6.1.1's `surface_Delaunay_remeshing.h:101-129` documents `protect_constraints` as "if true, enable protection of features / constrained edges / polyline constraints" — there is no `4/3` precondition; the remesher's Delaunay refinement engine resamples constrained polylines internally regardless.

The actual reason `protect_constraints=true` is correct here is that **both meshes A and B must see the SAME polyline vertex sequence after remeshing for pairwise conformality** (as the alternative plan PLAN_corefine_pipeline.md:L262-L268 explicitly states). Pre-resampling guarantees both Delaunay-remeshing calls receive bit-identical inputs; `protect_constraints=true` then preserves those vertices through the remesh.

This is QUALITY rather than CRITICAL because the chosen flag value is correct; only the rationale is wrong. But a future reviewer who reads the Defence section may incorrectly conclude that `protect_constraints=false` would also be valid (it would not — it would let the remesher collapse polyline sub-edges and break conformality between A_out and B_out).

**Trigger:**
A future plan revision swaps in `isotropic_remeshing` based on the (incorrect) belief that the same precondition reasoning applies, leading to silent conformality loss.

**Actual behavior:**
Plan cites a precondition that does not apply to the function being invoked.

**Expected behavior:**
Cite the actual reason: shared polyline vertex sequence between A_out and B_out is the hard requirement; `protect_constraints=true` is what enforces it.

**Suggested fix:**
Rewrite the bullet in §Defence:

```diff
-- **Why `protect_constraints=true` in `surface_Delaunay_remeshing` rather than `false`?** Because the polyline is *already resampled* to satisfy the `≤ 4/3·mesh_edge_size` precondition; `true` then preserves the polyline exactly (highest quality), whereas `false` would re-collapse the resampled polyline (sliver-recovery, but at the cost of polyline drift between A and B). The trap (PLAN_corefine_pipeline.md:L342-L351) is "raw corefine output as the constraint with `protect_constraints=true`", which we avoid by resampling first.
+- **Why `protect_constraints=true` in `surface_Delaunay_remeshing` rather than `false`?** Pairwise conformality requires that A_out and B_out share a bit-identical polyline-vertex sequence. Pre-resampling (Phase 1 step 5) ensures both `surface_Delaunay_remeshing` calls receive the same input polyline; `protect_constraints=true` then preserves that vertex sequence through remeshing. With `protect_constraints=false`, the remesher would re-collapse sub-edges and the two meshes would drift apart (each call's collapse pattern depends on the local 1-ring, which differs between A and B). The "≤ 4/3·target" precondition cited in EXPLORE_cgal_corefine.md:L142-L145 belongs to `isotropic_remeshing`, not `surface_Delaunay_remeshing` — the latter has no analogous precondition (verified at `cgal-61/include/CGAL/Polygon_mesh_processing/surface_Delaunay_remeshing.h:87-129`). The true risk we avoid by pre-resampling is *cross-mesh polyline drift*, not a documented precondition violation.
```

**Test case:** N/A (rationale-only change; no behavioral test).

---

### [R-007] [MODERATE] [PLAN] [§Phase 1 `intersection_graph::sort_by_overlap_volume`] — Tie-breaker rule unspecified; "deterministic ordering on tied volumes" cannot be enforced as written

**Category:** ASSUMPTION (specification gap)

**Description:**
Phase 5's testing strategy includes a unit test "`intersection_graph::sort_by_overlap_volume`: synthetic 4-mesh case; assert deterministic order on tied volumes." But Phase 1's spec for the function is:
> Sort pairs by descending bbox-overlap volume (a heuristic; see requirement 3 above).

It does not specify how to break ties. `std::sort` is not stable; the behavior on equal keys depends on the underlying algorithm and input order. Two meshes whose pairwise bbox-overlap volumes are equal (or differ only at the floating-point LSB) will sort in implementation-defined order.

This matters because the plan's §Risk Assessment explicitly says "Multi-pair corefine ordering" affects which polyline gets corefined first, and "exact predicates may produce a slightly different polyline than expected" depending on order. Determinism of the order is therefore a correctness requirement, not just a test-of-test concern.

**Trigger:**
Two or more pairs with bbox-overlap volume equal to within the floating-point representation, or any input where `std::sort` happens to reorder the tied pairs differently across compilers/STL versions.

**Actual behavior:**
Implementation-defined order; bit-different (but valid) outputs across rebuilds.

**Expected behavior:**
Total order. Recommended tie-breaker: lexicographic on `(i, j)` indices ascending.

**Suggested fix:**
Add tie-breaker spec to Phase 1's `intersection_graph.h` interface:

```diff
 // Sort pairs by descending bbox-overlap volume (a heuristic; see
 // requirement 3 above).
+// Tie-breaker: when two pairs have bbox-overlap volume equal within
+// 1e-12 × max-volume, break ties by lexicographic ascending order on
+// (i, j).  Implementation: use std::sort with a comparator that
+// returns (vol_a > vol_b) when |vol_a - vol_b| > tol; otherwise
+// (i_a, j_a) < (i_b, j_b).
 void sort_by_overlap_volume(std::vector<std::pair<std::size_t, std::size_t>>& pairs,
                              const std::vector<Mesh>& meshes);
```

And add a matching acceptance criterion in Phase 1:

```diff
+- [ ] `intersection_graph::sort_by_overlap_volume` is deterministic across rebuilds: synthetic 4-mesh fixture with two pairs of equal bbox-overlap volume sorts identically across two independent invocations.
```

**Test case:**
```python
def test_R007_sort_by_overlap_volume_deterministic():
    """Synthetic 4-mesh case where pairs (0,1) and (2,3) have identical
    bbox-overlap volumes must sort in lexicographic order on indices."""
    pairs = compute_overlap_pairs(meshes=[m0, m1, m2, m3])
    sort_by_overlap_volume(pairs, meshes)
    pairs2 = compute_overlap_pairs(meshes=[m0, m1, m2, m3])
    sort_by_overlap_volume(pairs2, meshes)
    assert pairs == pairs2, "Sort is non-deterministic across invocations"
    # And lexicographic on the tied subset:
    tied = [p for p in pairs if overlap_volume(p, meshes) == max_volume]
    assert tied == sorted(tied), "Tie-breaker is not lexicographic on (i,j)"
```

---

### [R-008] [MODERATE] [PLAN] [§Phase 3 generator] — `LC_NEAR = 1500.0` hardcoded in `.geo` template; generator does not propagate from manifest, so non-default `mesh-edge-size` runs produce a `.geo` whose size field is wrong

**Category:** DEVIATION (Phase 4 sweep cannot regenerate matched `.geo` without manual edits)

**Description:**
Phase 3's `.geo` template hardcodes:
```gmsh
LC_MIN     =  100.0;
LC_NEAR    = 1500.0;
LC_FAR     = 10000.0;
```
The generator's CLI flags are listed as `--manifest`, `--res`, `--out-geo`, `--pad-xy`, `--pad-top`, `--pad-bottom`. There is **no flag to override `LC_NEAR`**, and the template substitutes only the bbox values (`<XMIN>`, ..., `<ZMIN>`) and the fault filenames.

Phase 4's sensitivity sweep tries four configs that vary `mesh-edge-size` (1500, 2000, 1000, 1500) but reuses the same `.geo` template. The `coarser` (2000) and `finer` (1000) configs would produce STLs whose triangle edges don't match `LC_NEAR = 1500` in the `.geo`. The sweep cannot honestly compare bulk-tet quality across configs because the `.geo`'s size field is constant while the embedded surface size varies.

**Trigger:**
Any Phase 4 sweep run with `mesh-edge-size != 1500`. Also any future `--res 1000` or `--res 500` run.

**Actual behavior:**
Generated `.geo` has `LC_NEAR=1500` regardless of the input STL triangle size. Size-field-vs-embedded-surface mismatch makes the bulk tet quality stats config-dependent in ways the sweep doesn't isolate.

**Expected behavior:**
The generator reads `mesh_edge_size` from `manifest.json` and writes that value as `LC_NEAR` (or accepts `--lc-near` CLI flag).

**Suggested fix:**
Add `mesh_edge_size` propagation to the generator. Edit Phase 3 §Detailed Requirements:

```diff
 - Reads `data_corefined/manifest.json`.
 - Computes the union bbox (`xmin_fault`, `xmax_fault`, ...) by walking `meshes[].bbox`.
+- Reads `mesh_edge_size` from `manifest.json` (top-level field, written by
+  Phase 1 `corefine_set` per its manifest schema).  Substitutes this value
+  for `LC_NEAR` in the template.  This makes the `.geo`'s near-fault size
+  field consistent with the surface STL triangulation produced by corefine.
 - Fills the `<XMIN>`, `<XMAX>`, etc. placeholders.
 - Fills the six `<fault_k>` placeholders in the merged-files block, in alphabetical order of basename.
 - Writes the `.geo` to `code_meshing/safs_multifault_box_<R>m.geo`.
-- CLI: `--manifest`, `--res`, `--out-geo`, `--pad-xy`, `--pad-top`, `--pad-bottom`. Defaults match the template.
+- CLI: `--manifest`, `--res`, `--out-geo`, `--pad-xy`, `--pad-top`, `--pad-bottom`, `--lc-near` (default: read from manifest), `--lc-min` (default: read manifest's `min_edge`), `--lc-far` (default: 10·`lc-near`). Defaults match the template only when the manifest's `mesh_edge_size = 1500`.
```

And update the `.geo` template's size-controls block to be templated:

```diff
 // 2. Mesh size controls (HARD floor: LC_MIN must equal Phase 1 --min-edge).
-LC_MIN     =  100.0;
-LC_NEAR    = 1500.0;
-LC_FAR     = 10000.0;
+LC_MIN     = <LC_MIN>;     // = manifest.min_edge
+LC_NEAR    = <LC_NEAR>;    // = manifest.mesh_edge_size
+LC_FAR     = <LC_FAR>;     // = 10·LC_NEAR by default; override via --lc-far
```

**Test case:**
```python
def test_R008_generator_propagates_mesh_edge_size():
    """Running the generator on a manifest with mesh_edge_size=2000 must
    produce a .geo whose LC_NEAR equals 2000."""
    manifest = make_manifest(mesh_edge_size=2000.0, min_edge=200.0)
    geo_path = tmp_path / "out.geo"
    generate_multifault_geo(manifest, geo_path, res=2000)
    assert parse_geo_constant(geo_path, "LC_NEAR") == 2000.0
    assert parse_geo_constant(geo_path, "LC_MIN")  == 200.0
```

---

### [R-009] [LOW] [PLAN] [§Phase 5 `test_pairwise_conformality`] — Tolerance `1e-6 · bbox_diag` is 0.5 m for the SAFS bbox; far looser than the actual conformality goal

**Category:** EDGE_CASE (test threshold misaligned with the contract)

**Description:**
Phase 5 step 5 (`test_pairwise_conformality`) asserts polyline vertex match within `1e-6 · bbox_diag`. The SAFS bbox is approximately `(443 km, 296 km, 18 km)` per `safs_fault_box_freesurface_clip.geo:L36-L40`; bbox_diag ≈ √(443² + 296² + 18²) km ≈ 533 km = 5.33×10⁵ m. `1e-6 × bbox_diag ≈ 0.53 m`.

The Phase 4 `.geo` sets `Geometry.Tolerance = 1e-3` (1 mm) — gmsh dedups vertices within 1 mm. The plan's whole conformality story depends on polyline coordinates matching to **better than 1 mm**. A test that allows 0.53 m of drift will pass even when conformality is broken at the gmsh-tolerance level.

This is LOW because the OFF round-trip (post R-004 fix) and `%.15g` STL writer (post R-005 fix) will give sub-mm precision in practice, so the test will still pass when conformality is healthy. But a test threshold larger than the contract makes the test useless for catching the failure modes it's supposed to catch.

**Trigger:**
Any polyline-conformality regression that drifts vertices by > 1 mm but < 0.5 m. Test reports green; gmsh fails to dedup.

**Actual behavior:**
Test threshold (0.5 m) is ~500× looser than the gmsh tolerance (1 mm).

**Expected behavior:**
Test threshold ≤ `Geometry.Tolerance` so that a passing test guarantees gmsh will dedup correctly.

**Suggested fix:**
Tighten the threshold in Phase 5 step 5 from relative to absolute:

```diff
-5. `test_pairwise_conformality`: for every pair in `manifest.pairs`, load both meshes; for every polyline vertex, assert it appears in both meshes within `1e-6 · bbox_diag`.
+5. `test_pairwise_conformality`: for every pair in `manifest.pairs`, load both meshes; for every polyline vertex, assert it appears in both meshes within `1e-6 m` absolute (this is 1000× tighter than the Phase 3 `Geometry.Tolerance = 1e-3`, so a passing test guarantees gmsh will dedup correctly).
```

Also tighten Phase 1 acceptance criterion 4:

```diff
-- [ ] STL ASCII round-trip preserves polyline conformality: writing
-      output as ASCII STL with `%.10g`, reading back via `meshio`,
-      every polyline vertex coordinate diff between A and B < 1e-6.
+- [ ] STL ASCII round-trip preserves polyline conformality: writing
+      output as ASCII STL with `%.15g`, reading back via `meshio`,
+      every polyline vertex coordinate diff between A and B < 1e-6 m
+      (absolute; not bbox-relative).
```

**Test case:**
```python
def test_R009_conformality_tolerance_below_geometry_tolerance():
    """The conformality test threshold must be at least 100x tighter
    than the gmsh Geometry.Tolerance to guarantee dedup."""
    geometry_tolerance_m = 1e-3
    test_threshold_m = 1e-6  # post-R009 fix
    assert test_threshold_m < geometry_tolerance_m / 100, \
        f"Test threshold {test_threshold_m} not tight enough vs gmsh tol {geometry_tolerance_m}"
```

---

### [R-010] [LOW] [PLAN] [§Phase 3 `.geo` template] — Hardcoded Physical Surface ID-to-fault mapping (101..106) is fragile to fault renaming/reordering

**Category:** QUALITY (downstream MFEM BC mapping silently breaks if fault list changes)

**Description:**
The `.geo` template hardcodes:
```gmsh
Physical Surface("fault_COAV_Mission_Creek",        101) = {fault_surfs[0]};
Physical Surface("fault_MJVS_San_Andreas",          102) = {fault_surfs[1]};
...
Physical Surface("fault_SBMT_San_Andreas",          106) = {fault_surfs[5]};
```
The plan's prose says "Fills the six `<fault_k>` placeholders ... in alphabetical order of basename" and "the per-Physical-Surface IDs (101..106) are stable across runs so MFEM's BC-attribute mappings remain valid". But the generator is described as substituting `<fault_k>` placeholders only; it does not regenerate the Physical Surface lines from the actual fault basenames. If the manifest's mesh order is ever non-alphabetical (e.g. after intentional reordering, or if Phase 2's input glob returns files in filesystem order on a different OS), the IDs and names misalign silently — `Physical Surface("fault_COAV_Mission_Creek", 101)` could end up pointing at the SBMT-Garnet-Hill geometry.

The plan says "If a fault is removed or added in a future run, the IDs above must be updated by hand". This is OK for the locked 2000 m fixture but fragile.

**Trigger:**
Manifest with a non-alphabetical mesh order (e.g., a future generator improvement that orders by bbox volume). Or any change to the SAFS fault list.

**Actual behavior:**
Hardcoded physical-surface labels can mismatch the underlying STL.

**Expected behavior:**
The generator emits the `Physical Surface(name, id) = {fault_surfs[k]};` lines using the actual basenames from the manifest, in alphabetical order, with `id = 100 + k + 1`.

**Suggested fix:**
Edit Phase 3 §Detailed Requirements to make the Physical-Surface lines generator-emitted, not template-static:

```diff
 - Fills the six `<fault_k>` placeholders in the merged-files block, in alphabetical order of basename.
+- **Generator-emits** (not template-static) the `Physical Surface(name, 100+k+1) = {fault_surfs[k]};` block, with `name` derived from the alphabetical basename (strip `SAFS-SAFZ-` prefix and `_<R>m_corefined` suffix; replace `-` and `_` with `_` to keep gmsh-legal identifiers).  This guarantees the name-to-tag mapping always reflects the actual STL list, so a future fault add/remove updates the `.geo` automatically.
 - Writes the `.geo` to `code_meshing/safs_multifault_box_<R>m.geo`.
```

And the corresponding template block becomes:

```diff
 // 7. Physical groups.
 Physical Volume("rock", 1) = {1};
-Physical Surface("fault_COAV_Mission_Creek",        101) = {fault_surfs[0]};
-Physical Surface("fault_MJVS_San_Andreas",          102) = {fault_surfs[1]};
-Physical Surface("fault_MULT_S_SAF_Banning",        103) = {fault_surfs[2]};
-Physical Surface("fault_SBMT_Garnet_Hill",          104) = {fault_surfs[3]};
-Physical Surface("fault_SBMT_Mission_Creek",        105) = {fault_surfs[4]};
-Physical Surface("fault_SBMT_San_Andreas",          106) = {fault_surfs[5]};
+<PHYSICAL_FAULT_SURFACES>      // emitted by the generator: one line per
+                               // alphabetically-sorted fault basename,
+                               // tag = 100 + k + 1
 Physical Surface("top",    200) = {6};
```

**Test case:**
```python
def test_R010_generator_emits_physical_surfaces_from_manifest():
    """Adding a 7th fault to the manifest must yield 7 Physical Surface
    lines in the .geo with consecutive tags 101..107."""
    manifest = make_manifest_with_n_faults(7)
    geo = tmp_path / "out.geo"
    generate_multifault_geo(manifest, geo, res=2000)
    surf_lines = grep_geo(geo, r"^Physical Surface\(\"fault_")
    assert len(surf_lines) == 7
    tags = [parse_physical_surface_tag(s) for s in surf_lines]
    assert tags == [101, 102, 103, 104, 105, 106, 107]
```

---

## Summary
- Critical issues: 4 (R-001, R-002, R-003, R-004)
- Moderate issues: 4 (R-005, R-006, R-007, R-008)
- Low issues: 2 (R-009, R-010)
- Plan compliance: PARTIAL (the plan covers the alternative-plan structure but introduces 4 critical errors of its own — `tangential_relaxation` scope creep, IO re-entrancy, default mismatch, OFF precision gap)
- Verdict: **FAIL — must fix before proceeding**. R-001 alone destroys the conformality contract that motivates the entire pipeline; R-002 makes Phase 5 acceptance criterion 5 ("Re-running with the same args is idempotent") impossible to satisfy as written; R-003 makes the documented end-to-end pipeline produce a different mesh than the sensitivity-sweep baseline claims to defend; R-004 silently destroys conformality before the high-precision STL writer is reached.

## Unreviewed Areas
- The CGAL 6.1.1 multi-mesh corefine question (plan §Open questions): not investigated. The plan flagged this for the implementer to verify; that flag is appropriate.
- Phase 4 sensitivity sweep's expected quality numbers: the plan does not commit to specific quality bounds per config — the report is descriptive. Cannot review numbers that don't exist.
- gmsh ≥ 30 minute time budget: untestable without running the pipeline. Not a correctness issue.
