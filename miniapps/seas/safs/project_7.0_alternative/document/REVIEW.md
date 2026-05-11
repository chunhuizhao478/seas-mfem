# Code Review: NW Hard-Cut Plan (PLAN_nw_hard_cut.md)

**Date:** 2026-05-08
**Reviewer:** /code-review (adversarial)

## Review Scope

- Plan: `project_7.0_alternative/document/PLAN_nw_hard_cut.md`
- Implementation: **none yet** — only the plan exists. This review treats the
  plan as a specification that `/code-implement` will follow mechanically;
  bugs in the spec become bugs in the code.
- Files reviewed: the plan itself, plus the upstream code the plan
  references (`code_preprocess/clean_freesurface_mesh.py`,
  `code_preprocess/ts_to_stl.py`) for consistency-of-conventions checks.
- Domain context: `seas-mfem-safs/CLAUDE.md`,
  `miniapps/seas/CLAUDE.md`, the user's auto-memory entries on
  BP5/dynamic-folder editability.
- Real STL data was sampled to validate concrete numerical claims in the
  plan (Phase 1 acceptance criteria).

## Findings

### [R-001] [CRITICAL] Phase 1 / nw_anchor — `@` operand order is invalid for `(2,) @ (N, 2)`

**Category:** BUG

**Description:**
Plan §Phase 1, step 3 specifies the implementation as:

> compute `s = NW_DIRECTION_XY @ verts[:, :2]`

`NW_DIRECTION_XY` has shape `(2,)`. `verts[:, :2]` has shape `(N, 2)`.
NumPy's `@` rule for a 1-D LHS promotes it to `(1, 2)`; for the matmul
`(1, 2) @ (N, 2)` to be valid, the inner dim of the LHS (= 2) must equal
the second-to-last dim of the RHS (= N). For any real fault mesh
(N = 1479, 13722, 54168, 214938) this raises
`ValueError: matmul: Input operand 1 has a mismatch in its core dimension 0`.

I verified this empirically:

```text
NW.shape = (2,);  verts[:, :2].shape = (3, 2)
NW @ verts[:, :2]   →  ValueError (size 3 is different from 2)
verts[:, :2] @ NW   →  array([0.707…, 0.707…, 0.707…])   ✓
```

**Trigger:**
First call to `nw_anchor(verts)` on the preferred mesh (or any real STL).
Phase 1 acceptance criterion #2 is impossible to satisfy as written.

**Actual behavior (per the spec):**
`ValueError` at the dot-product line; the entire driver crashes before
reaching Phase 2 or Phase 3.

**Expected behavior:**
Compute the row-wise dot product of every vertex's `(x, y)` with
`NW_DIRECTION_XY`, returning a 1-D array of shape `(N,)`.

**Suggested fix:**
Reverse the operand order so the matmul is `(N, 2) @ (2,) → (N,)`.

```diff
-   Implementation: compute `s = NW_DIRECTION_XY @ verts[:, :2]`; let
-   `s_max = s.max()`; `mask = s >= s_max - 1e-9`; among the masked rows,
-   pick `argmin(z)`. Return that row and `s_max`.
+   Implementation: compute `s = verts[:, :2] @ NW_DIRECTION_XY` (shape
+   `(N, 2) @ (2,) → (N,)`); let `s_max = float(s.max())`;
+   `mask = s >= s_max - EPS`; among the masked rows, pick `argmin(z)`.
+   Return that row (as a `(3,)` float64 array) and `s_max`.
```

The same fix should be applied wherever else the plan reuses the `n_h ·
verts[:, :2]` shorthand. Phase 1 step 3 is the only place this exact
expression appears in the plan, but Phase 4 §5 ("post-cut s_max ≤
s_max_pref + 1e-6") will re-compute the projection over a freshly loaded
mesh; the implementer must use the correct operand order there too.

**Test case:**

```python
def test_R001_nw_anchor_does_not_crash_on_real_n():
    import numpy as np
    from nw_cut_strip import nw_anchor, NW_DIRECTION_XY
    rng = np.random.default_rng(0)
    verts = rng.uniform(-1, 1, size=(1000, 3))   # N=1000 != 2
    p, s = nw_anchor(verts)
    assert p.shape == (3,)
    assert isinstance(s, float)
    # Brute-force the same answer to verify the row-wise dot is correct.
    s_brute = (-verts[:, 0] + verts[:, 1]) / np.sqrt(2.0)
    j = int(np.argmax(s_brute))
    np.testing.assert_allclose(p, verts[j], atol=0)
    np.testing.assert_allclose(s, s_brute[j], atol=0)
```

---

### [R-002] [MODERATE] Phase 3 / CLI — `--print-anchor`, `--batch`, and positional `INPUT OUTPUT` lack mutual-exclusion semantics

**Category:** ASSUMPTION

**Description:**
Plan §Phase 3 step 2 documents the CLI as:

```
nw_cut_strip.py [--reference PATH] [--alt-dir DIR] [--glob PAT]
                [--suffix SFX] [--batch | INPUT OUTPUT]
                [--print-anchor PATH] [--keep-orphan-verts]
                [--verbose | -q]
```

The pipe in `[--batch | INPUT OUTPUT]` suggests these are mutually
exclusive, but `--print-anchor PATH` is shown as a separate optional flag
with no relationship documented. In Phase 1 step 5 it says
`--print-anchor` should "exit 0 immediately after printing", which
implies the positional `INPUT OUTPUT` and `--batch` are **not** required
when `--print-anchor` is given. As written, an `argparse` parser that
follows the synopsis literally will either:
- require `INPUT OUTPUT` always (when no `--batch` is set), making
  `--print-anchor` unusable on its own, or
- silently accept `--print-anchor` together with `--batch` and run both
  paths, with undefined ordering.

**Trigger:**
User runs `python nw_cut_strip.py --print-anchor <preferred-stl>` with no
positional args (the example given in the plan).

**Actual behavior (per the spec):**
Either `argparse` errors with "the following arguments are required:
input, output", contradicting Phase 1 acceptance criterion #4, or the
parser silently allows ambiguous combinations.

**Expected behavior:**
`--print-anchor PATH` is a "do this and exit" flag, mutually exclusive
with both `--batch` and the positional `INPUT OUTPUT`. The parser should
reject any combination that supplies more than one mode at once.

**Suggested fix:**
Make the modes explicit and mutually exclusive in the spec, e.g. with an
`argparse.add_mutually_exclusive_group()`. Add this paragraph after the
CLI synopsis in Phase 3 step 2:

```diff
+   The three top-level modes — `--print-anchor PATH`, `--batch`, and the
+   positional `INPUT OUTPUT` form — are mutually exclusive. The
+   implementer must use `argparse.add_mutually_exclusive_group()` (or
+   equivalent post-parse validation) and exit 2 with a clear message if
+   the user combines them. `INPUT` and `OUTPUT` together count as a
+   single mode (both required if either is given).
```

**Test case:**

```python
def test_R002_print_anchor_alone_works():
    import subprocess, sys
    result = subprocess.run(
        [sys.executable, "nw_cut_strip.py", "--print-anchor",
         str(REFERENCE_STL)],
        capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.startswith("anchor:")

def test_R002_print_anchor_with_batch_rejected():
    import subprocess, sys
    result = subprocess.run(
        [sys.executable, "nw_cut_strip.py", "--print-anchor",
         str(REFERENCE_STL), "--batch"],
        capture_output=True, text=True)
    assert result.returncode != 0
    assert "mutually exclusive" in result.stderr.lower() \
        or "not allowed with" in result.stderr.lower()
```

---

### [R-003] [MODERATE] Phase 3 / loading — meshio may return `float32` vertices, silently degrading clipping precision

**Category:** ASSUMPTION

**Description:**
Plan §Phase 3 step 3.4 says:

> Load the alternative file as `(verts_alt, faces_alt)` using `meshio` …
> Use the *vertex-matrix-and-face-matrix* form, not the dict form,
> because Phase 2 expects ndarrays.

Phase 2's `clip_mesh_at_plane` signature in step 4 says
`verts: (N, 3) float64`. But `meshio.read` for an ASCII STL can yield
`m.points` of dtype `float32` or `float64` depending on the meshio
version and the file's content (binary STLs are float32 in the spec;
ASCII STLs default to float64 in modern meshio but some versions still
return float32). The plan does not specify a dtype-coercion step.

If `verts_alt.dtype == np.float32`, `signed_dist(p, n, c) = n·p − c`
loses ~5 decimal digits (rel-eps ≈ 1.2e-7) on coordinates of magnitude
~4e6. The cutting plane is at `c ≈ 2.456e6`; any vertex within ~0.5 m of
the plane could end up on the wrong side after the float32 round-off,
producing a jagged edge instead of a clean cut.

**Trigger:**
A meshio version that decodes the alternative STL as float32; Phase 3
step 3.4 says nothing about dtype.

**Actual behavior:**
Silent precision loss; the clipped boundary may be ~1 m off and
non-deterministic across meshio versions.

**Expected behavior:**
Always cast input vertices to float64 before clipping. Match the
`np.array(..., dtype=np.float64)` style used in
`clean_freesurface_mesh.clip_to_arrays`.

**Suggested fix:**
Insert into Phase 3 step 3.4:

```diff
   4. Load the alternative file as `(verts_alt, faces_alt)` using `meshio`
      (or pymeshlab MeshSet → `current_mesh().vertex_matrix()` /
      `face_matrix()`; either is fine). Use the *vertex-matrix-and-face-
      matrix* form, not the dict form, because Phase 2 expects ndarrays.
+     Immediately after loading, coerce dtypes:
+     `verts_alt = np.ascontiguousarray(verts_alt, dtype=np.float64)`,
+     `faces_alt = np.ascontiguousarray(faces_alt, dtype=np.int32)`.
+     This is required because `meshio.read` of an ASCII STL can return
+     `float32` on some versions, and `clip_mesh_at_plane`'s contract
+     specifies `float64` (precision loss on a ~4e6-magnitude UTM coord
+     base would otherwise cause non-deterministic mis-classification of
+     vertices within ~1 m of the cutting plane).
```

The same `dtype=np.float64` requirement should be added explicitly to
the Phase 1 `load_vertices` return-value docstring (it already says
"float64" in the docstring, but the spec doesn't enforce it for the
non-`.ts` branch).

**Test case:**

```python
def test_R003_clip_uses_float64_even_if_meshio_gives_float32(tmp_path):
    import numpy as np
    from nw_cut_strip import nw_cut_file, cutting_plane
    # Build a tiny synthetic mesh in float32 and verify the clip plane
    # is hit to ~1e-9 precision in float64, not ~1e-7 (float32).
    V32 = np.array([[0, 0, 0], [3e6, 0, 0], [0, 3e6, 0]], dtype=np.float32)
    F = np.array([[0, 1, 2]], dtype=np.int32)
    # Force-load via meshio if needed; here we exercise the function
    # directly. The relevant assertion: signed distances of new clip
    # vertices < 1e-6.
    n, c = cutting_plane(np.array([1.0e6, 1.0e6, 0.0]))
    # ... call clip_mesh_at_plane with V32; verify dtype after coercion
```

---

### [R-004] [MODERATE] Phase 3 step 3.6 — empty-output handling specified for CLI but not for the importable entry point

**Category:** EDGE_CASE

**Description:**
Plan §Phase 3 step 3.6:

> If `len(F_out) == 0`: print a clear error … and exit non-zero.

But Phase 3 §Interfaces also exposes a programmatic entry point:

```python
def nw_cut_file(input_path, output_path, plane, cleanup=True, verbose=True) -> dict:
```

The spec doesn't say what `nw_cut_file` does when every triangle is
clipped away. Options the implementer might pick (and any of which is a
silent contract violation if guessed wrong):

- raise an exception (preferred — matches "exit non-zero" intent)
- write an empty `solid` STL and return `stats` with `n_out=0`
- skip the file and return `stats` with `n_out=0`

The downstream meshing pipeline (`safs_fault_box*.geo`) would silently
load the empty STL and produce an unmeshable input — exactly the kind
of silent failure that CLAUDE.md tells us to avoid.

**Trigger:**
User passes `--reference` pointing at the wrong file (e.g., the
alternative mesh itself, whose anchor sits NW of every other vertex →
nothing kept). This is the same scenario Phase 3 step 3 lists under
"Edge Cases", but only for the CLI path.

**Suggested fix:**

```diff
   6. If `len(F_out) == 0`: print a clear error
-     `"refusing to write empty STL: every triangle was clipped away"` and
-     exit non-zero. (This shouldn't happen on the real data … but a typo
-     in `--reference` could trigger it.)
+     `"refusing to write empty STL: every triangle was clipped away"` and
+     exit non-zero in the CLI path. In the `nw_cut_file()` programmatic
+     path, raise `ValueError(<same message>)`. Do **not** write the empty
+     STL in either case — the downstream gmsh `.geo` files would silently
+     load it and emit an unmeshable input.
```

**Test case:**

```python
def test_R004_nw_cut_file_raises_on_empty_clip(tmp_path):
    import numpy as np, pytest
    from nw_cut_strip import nw_cut_file, cutting_plane
    # Anchor "before" all of the alternative mesh → everything dropped.
    bad_anchor = np.array([0.0, 1e9, 0.0])  # huge NW projection
    plane = cutting_plane(bad_anchor)
    with pytest.raises(ValueError, match="empty STL"):
        nw_cut_file(ALT_2000M_STL, tmp_path / "out.stl", plane)
    assert not (tmp_path / "out.stl").exists()
```

---

### [R-005] [MODERATE] Phase 3 step 3.9 — ambiguity on whether `s_max_after` is measured before or after pymeshlab cleanup

**Category:** ASSUMPTION

**Description:**
Phase 3 step 3.9:

> Print the post-cut bbox, the number of triangles in / out, and the
> max-NW-projection of the surviving mesh, which by construction must
> satisfy `s_max_after ≤ s_max + EPS`.

Step 3.7 (immediately before) runs the pymeshlab cleanup chain
(`meshing_remove_duplicate_vertices`, `meshing_remove_null_faces`,
`meshing_repair_non_manifold_edges(method="Remove Faces")`, …). The
spec doesn't say whether `s_max_after` is computed from `V_out`
(pre-cleanup) or from the cleaned `MeshSet.current_mesh().vertex_matrix()`
(post-cleanup, what's actually written to disk).

These can differ subtly: `repair_non_manifold_edges(method="Remove
Faces")` deletes faces but keeps the underlying vertex array; if
`meshing_remove_unreferenced_vertices` runs after that, then orphaned
vertices are stripped. So the "max NW projection of the surviving mesh"
is always less-than-or-equal to the pre-cleanup projection.

By itself this is a small ambiguity, but it interacts with the Phase 3
acceptance criterion ("post-cut max NW projection ≤ 2456425.1 + 1e-6"):
if the implementer computes from `V_out` they may report a value that
does not match what's on disk.

**Suggested fix:**

```diff
   9. Print the post-cut bbox, the number of triangles in / out, and the
-     max-NW-projection of the surviving mesh, which by construction must
-     satisfy `s_max_after ≤ s_max + EPS`.
+     max-NW-projection of the surviving mesh **after the pymeshlab
+     cleanup chain has run** (i.e. read it from
+     `ms.current_mesh().vertex_matrix()`, not from the pre-cleanup
+     `V_out`). By construction this must satisfy
+     `s_max_after ≤ s_max + EPS`. The acceptance criteria in this Phase
+     refer to the same post-cleanup quantity.
```

**Test case:**

```python
def test_R005_smax_after_is_measured_post_cleanup(tmp_path):
    import numpy as np, meshio
    from nw_cut_strip import nw_cut_file, cutting_plane, NW_DIRECTION_XY
    plane = cutting_plane(np.array([365072.3, 3838982.0, 0.0]))
    out = tmp_path / "alt_2000m_nwcut.stl"
    stats = nw_cut_file(ALT_2000M_STL, out, plane)
    # Reload the saved file and recompute s_max independently.
    m = meshio.read(str(out))
    s = m.points[:, :2] @ NW_DIRECTION_XY
    assert s.max() <= plane[1] + 1e-6
    assert abs(stats["s_max_out"] - s.max()) <= 1e-6
```

---

### [R-006] [MODERATE] Phase 2 step 4 — leftover "let me restate" prose in the `clip_mesh_at_plane` docstring

**Category:** QUALITY (but actively confusing → could mask a bug)

**Description:**
Phase 2 step 4 contains a self-correction block left in the
specification:

```text
The semantics of '2to1' and '1to2' refer to keep-side counts: '2to1' = 2
keep + 1 drop → 1 output triangle is wrong; let me restate: as in the
existing code, the case names track *number-above-the-plane*, where
'above' = drop side. So:
  - n_above == 0 (keep all 3): kept_whole
  - n_above == 1: 2 below + 1 above → emit 2 triangles ('1to2')
  - n_above == 2: 1 below + 2 above → emit 1 triangle ('2to1')
  - n_above == 3: dropped.
```

A `/code-implement` agent reading this verbatim is likely to:
1. include the contradictory "is wrong; let me restate" prose in the
   docstring of `clip_mesh_at_plane`, polluting the user-facing
   documentation; or
2. mis-name the keys in `stats` (e.g. swap `n_clipped_2to1` and
   `n_clipped_1to2`), since the conflicting "keep-side" framing in the
   first sentence is plausible enough that a hurried implementer might
   pick it.

The reference implementation `clean_freesurface_mesh.clip_to_arrays`
defines:
- `n_clipped_2to1` ↔ `n_above == 2` (1 keep + 2 drop → 1 triangle)
- `n_clipped_1to2` ↔ `n_above == 1` (2 keep + 1 drop → 2 triangles)

So the second framing in the plan is the correct one and must be the
only framing that survives.

**Suggested fix:**

```diff
       - n_clipped_2to1, n_clipped_1to2
     (same five counts that `clean_freesurface_mesh.clip_to_arrays`
-    reports). The semantics of '2to1' and '1to2' refer to keep-side
-    counts: '2to1' = 2 keep + 1 drop → 1 output triangle is wrong; let
-    me restate: as in the existing code, the case names track
-    *number-above-the-plane*, where 'above' = drop side. So:
+    reports). The case names track *number-above-the-plane*
+    (= number on the drop side), exactly as in
+    `clean_freesurface_mesh.clip_to_arrays`:
       - n_above == 0 (keep all 3): kept_whole
       - n_above == 1: 2 below + 1 above → emit 2 triangles ('1to2')
       - n_above == 2: 1 below + 2 above → emit 1 triangle ('2to1')
       - n_above == 3: dropped.
```

**Test case:**

```python
def test_R006_stats_keys_match_reference_implementation():
    import numpy as np
    from nw_cut_strip import clip_mesh_at_plane, cutting_plane
    # Build a 4-triangle mesh with 1 of each case (kept whole, dropped,
    # 2to1, 1to2) and assert the four counts.
    n, c = cutting_plane(np.array([0.0, 0.0, 0.0]))   # plane: -x+y = 0
    verts = np.array([
        [-1, -2, 0], [-1, -1, 0], [0, -1, 0],   # all keep (s<0): kept_whole
        [ 2,  1, 0], [ 1,  2, 0], [ 1,  1, 0],   # all drop (s>0): dropped
        [-1, -1, 0], [ 1, -3, 0], [ 1,  3, 0],   # 1 keep + 2 drop: 2to1
        [-1, -3, 0], [ 1,  3, 0], [-1,  1, 0],   # 2 keep + 1 drop: 1to2
    ], dtype=np.float64)
    faces = np.array([[0,1,2],[3,4,5],[6,7,8],[9,10,11]], dtype=np.int32)
    V_out, F_out, stats = clip_mesh_at_plane(verts, faces, n, c)
    assert stats["n_kept_whole"]   == 1
    assert stats["n_dropped"]      == 1
    assert stats["n_clipped_2to1"] == 1   # n_above==2 → 1 triangle
    assert stats["n_clipped_1to2"] == 1   # n_above==1 → 2 triangles
    # Total output triangles = 1 + 0 + 1 + 2 = 4
    assert len(F_out) == 4
```

---

### [R-007] [MODERATE] Phase 1 step 2 — `load_vertices` for `.ts` skips the orphan-vertex filter, breaking the contract

**Category:** DEVIATION

**Description:**
Phase 1 step 2 docstring promises an `(N, 3)` array of "vertex
coordinates from any mesh file readable by `meshio.read`. ASCII STL,
binary STL, PLY, OBJ, and `.ts` … must all work". The body of the spec
then says:

> - `.ts` branch: import `parse_ts` from `ts_to_stl` … and stack the
>   values of the returned `verts` dict.
> - All other suffixes: `meshio.read(str(path))` and stack `m.points`.
> - Vertices that are not referenced by any triangle are ignored …
>   (gather the unique vertex indices used by `m.cells` of type
>   `"triangle"` and slice `m.points` with that index array).

The orphan filter is described only in the third bullet, which by
flow-of-prose applies to "all other suffixes" (meshio path) but not to
`.ts`. Yet `parse_ts` returns `(verts: dict, tris: list)` and tris uses
1-indexed IDs that may not cover all of `verts.keys()` (the upstream
`.ts` files are known to carry orphan ATOM definitions that no triangle
references). The user's Phase 1 acceptance criterion #2 (returns a
specific anchor for the preferred 2000 m STL) is fine as written — but
if anyone later calls `load_vertices` on a `.ts` raw fault file, the
anchor will be polluted by orphan ATOMs.

This is a silent contract violation: the docstring promises an `(N, 3)`
array of triangle-referenced vertices, but the `.ts` branch returns all
declared vertices.

**Suggested fix:**

```diff
   - `.ts` branch: import `parse_ts` from `ts_to_stl` (same `sys.path.insert`
     trick already used at the top of `clean_freesurface_mesh.py`) and
-    stack the values of the returned `verts` dict.
+    stack the values of the returned `verts` dict, then apply the same
+    orphan filter as the meshio branch: collect the union of vertex IDs
+    referenced by any triangle in `tris` and keep only those rows.
   - All other suffixes: `meshio.read(str(path))` and stack `m.points`.
   - Vertices that are not referenced by any triangle are ignored …
```

**Test case:**

```python
def test_R007_load_vertices_drops_orphans_from_ts(tmp_path):
    # Synthesize a tiny .ts file with one orphan ATOM far NW of any
    # triangle, then assert nw_anchor() does NOT pick the orphan.
    ts = tmp_path / "tiny.ts"
    ts.write_text(
        "GOCAD TSurf 1\n"
        "VRTX 1 100 100 0\n"
        "VRTX 2 200 100 0\n"
        "VRTX 3 100 200 0\n"
        "VRTX 4 -1e9 1e9 0\n"   # orphan ATOM, very NW
        "TRGL 1 2 3\n"
        "END\n"
    )
    from nw_cut_strip import load_vertices, nw_anchor
    V = load_vertices(ts)
    assert V.shape == (3, 3)            # orphan stripped
    p, s = nw_anchor(V)
    # If orphan filter failed, anchor would be (-1e9, 1e9, 0).
    assert p[0] >= 100 - 1.0
```

---

### [R-008] [LOW] Plan mixes "below/above" (z=0 plane terminology) with "keep/drop" (generic plane terminology)

**Category:** QUALITY

**Description:**
Phase 4 §2 unit tests:

> - `test_one_keep_two_drop`: 1 below, 2 above → 1 output triangle …
> - `test_two_keep_one_drop`: 2 below, 1 above → 2 output triangles
>   whose union (as a polygon) equals the **SE-side quadrilateral**

The test name uses keep/drop, the body uses below/above (a z=0 framing
that does not generalise to the rotated cut planes the test is supposed
to exercise per Phase 4 §2 first sentence: "run all four with one
randomly-rotated cutting plane"). Worse, the test asserts the union
"equals the **SE-side** quadrilateral", but for a *random* plane the
keep half-space is not "SE" at all — it's the "negative-d half-space".
A literal-minded implementer may write a test that fails on rotated
planes.

This is purely terminological but it could mask the bug in R-001:
if the test for R-001's correctness uses "SE side" literally, and
R-001 is then "fixed" to swap the operand order, the test still fails
because the keep side is no longer geographically SE.

**Suggested fix:**

```diff
   2. Unit tests for `clip_triangle_at_plane` (run all four with one
      randomly-rotated cutting plane to break any axis-aligned coincidence):
-     - `test_all_keep`: triangle entirely inside `H_keep` → 1 output
-       triangle = input.
-     - `test_all_drop`: triangle entirely in NW half → 0 output triangles.
-     - `test_one_keep_two_drop`: 1 below, 2 above → 1 output triangle whose
-       two new vertices satisfy
-       `|signed_dist(v, n, c)| < 1e-6 * (|c| + 1)`.
-     - `test_two_keep_one_drop`: 2 below, 1 above → 2 output triangles
-       whose union (as a polygon) equals the SE-side quadrilateral; verify
-       by area sum.
+     - `test_all_keep`: triangle with all 3 vertices on the keep side
+       (`d ≤ 0`) → 1 output triangle = input.
+     - `test_all_drop`: triangle with all 3 vertices on the drop side
+       (`d > 0`) → 0 output triangles.
+     - `test_one_keep_two_drop`: 1 keep + 2 drop → 1 output triangle
+       whose two new vertices satisfy
+       `|signed_dist(v, n, c)| < 1e-6 * (|c| + 1)`.
+     - `test_two_keep_one_drop`: 2 keep + 1 drop → 2 output triangles
+       whose union (as a polygon) equals the **keep-side quadrilateral
+       in the cutting plane's local frame** (do not refer to compass
+       directions here — the test plane is randomly oriented). Verify by
+       polygon area sum: it must equal the area of the input triangle
+       intersected with `H_keep`, computed independently from the
+       parametric clip.
```

**Test case:** subsumed by the corrected wording above.

---

### [R-009] [LOW] Phase 4 §1 fixture spec uses ellipsis (`…`) for parameters

**Category:** QUALITY

**Description:**
Phase 4 §1: `make_unit_triangle(d_keep, d_drop, …)`. The trailing `…`
omits parameters the implementer must invent — at minimum the third
vertex's signed distance, and how the triangle is embedded in 3-D. A
mechanically-following implementer will pick something arbitrary and
the resulting tests may not exercise the same cases the plan intends.

**Suggested fix:**

```diff
   1. Synthetic test fixtures (no I/O against the real data):
-    - `make_unit_triangle(d_keep, d_drop, …)`: build a triangle in 3-D with
-      prescribed signed distances from a chosen plane. Use this to drive
-      the four cases below.
+    - `make_triangle_with_signed_dists(d0, d1, d2, n, c, *, rng=None)`:
+      build a triangle with three signed distances `d0, d1, d2` from the
+      plane `n·p = c`. Implementation: pick three random in-plane points
+      `q0, q1, q2` (project `rng.uniform(-10, 10, (3, 3))` onto the
+      plane) and translate each `qi` along `n` by `di` to land at the
+      requested signed distance. Returns `(p0, p1, p2)` as 3-tuples.
```

**Test case:** N/A (the fixture is itself a test helper).

---

### [R-010] [LOW] [POSSIBLE] Phase 2 step 2 — `lerp_to_plane` snap branch may fire for a pathological non-degenerate input

**Category:** EDGE_CASE [POSSIBLE]

**Description:**
Phase 2 step 2:

> If `abs(d_drop − d_keep) < EPS` (degenerate: both ~on plane), return
> `p_keep` snapped onto the plane (subtract `d_keep * n` from `p_keep`).

`EPS = 1e-9`. The condition `|d_drop − d_keep| < 1e-9` *could* fire for
two non-degenerate keep/drop vertices if both signed distances happen
to be near-equal but on opposite sides of the plane — e.g.
`d_keep = -5e-10`, `d_drop = +4e-10`. Their difference is `9e-10 < EPS`,
so the snap branch runs and the function returns `p_keep` snapped to
the plane. That's actually the right answer (the segment really is at
the plane to numerical noise), but the returned point is `p_keep -
d_keep * n` regardless of `d_drop`'s magnitude — so a `p_drop` that
sits 0.4 nm on the drop side gets effectively ignored. In practice the
result is fine to ~1 nm. But the spec's wording ("degenerate: both ~on
plane") could mislead a reviewer or a future implementer into thinking
the branch only fires when both are *on the keep side*; that's not what
the code will do.

**Suggested fix:**

```diff
   Implementation: `d_keep = signed_dist(p_keep, n, c)`,
   `d_drop = signed_dist(p_drop, n, c)`. If
-  `abs(d_drop − d_keep) < EPS` (degenerate: both ~on plane), return
+  `abs(d_drop − d_keep) < EPS` (segment is parallel to or in the plane
+  to within float64 noise — both endpoints have signed distances within
+  `EPS` of each other regardless of sign), return
   `p_keep` snapped onto the plane (subtract `d_keep * n` from `p_keep`).
```

**Test case:**

```python
def test_R010_lerp_handles_near_parallel_segment():
    import numpy as np
    from nw_cut_strip import lerp_to_plane, signed_dist, EPS
    n = np.array([-1.0, 1.0, 0.0]) / np.sqrt(2.0)
    c = 0.0
    # Two points whose signed distances differ by < EPS but straddle 0.
    p_keep = (1.0, 1.0 - 5e-10, 0.0)        # d ≈ -3.5e-10
    p_drop = (1.0, 1.0 + 4e-10, 0.0)        # d ≈ +2.8e-10
    q = lerp_to_plane(p_keep, p_drop, n, c)
    # Result should land exactly on the plane.
    assert abs(signed_dist(q, n, c)) <= 1e-9
```

---

## Summary

- Critical issues: **1** (R-001)
- Moderate issues: **6** (R-002, R-003, R-004, R-005, R-006, R-007)
- Low issues: **3** (R-008, R-009, R-010)
- Plan compliance: **N/A** (no implementation yet — this is a plan
  review, not a plan-vs-code review)
- Verdict: **FAIL — must fix before proceeding to /code-implement**

R-001 alone is a hard blocker: the very first numerical operation in the
driver (`s = NW_DIRECTION_XY @ verts[:, :2]`) raises `ValueError` on any
real input. R-002 through R-007 are spec ambiguities that will produce
silent-but-wrong implementations if the plan is followed literally;
fixing them is cheap (text edits in the plan). The Low items are
clarifications that don't affect correctness but will reduce review
friction once the implementation lands.

## Unreviewed Areas

- **No source code exists yet**, so Pass 2 (bugs in the implementation)
  and Pass 3 (code quality) only inspected the plan's pseudocode and
  example snippets. A second review pass is required after
  `/code-implement` produces `nw_cut_strip.py` and
  `test_nw_cut_strip.py`.
- The pymeshlab API calls
  (`meshing_remove_duplicate_vertices`, `meshing_repair_non_manifold_edges`)
  were not exercised against the actual installed pymeshlab version in
  `pythonenv`; the plan inherits these from
  `clean_freesurface_mesh.clean_freesurface_clip` which is known to work,
  so we did not re-validate.
- The interaction with the downstream meshing pipeline
  (`code_meshing/safs_fault_box*.geo`) was not reviewed; the plan
  declares this out of scope (outputs use a new `_nwcut.stl` suffix so
  existing `.geo` consumers are unaffected).
