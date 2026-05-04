# Code Review: Phase 3 — Pipeline Integration + 3D Tube Uniformity (2026-04-29)

## Review Scope

- **Plan:** `PLAN_cgal_corefine.md` Phase 3 (pipeline integration).
- **Files reviewed:**
  - `mesh/run_two_crossing_2000m_cgal.sh` (the new Phase 3 script)
  - `mesh/safs.geo` (gmsh size-field configuration)
  - `mesh/generate_safs_mesh.py` (driver wrapping gmsh)
  - `mesh/validate_msh.py:check_11_tube_uniformity` (the gate the user is implicitly testing)
- **Live verification:**
  - End-to-end run on Mill Creek × SBMT-SAF 2000 m completed at
    `mesh/output/two_crossing_2000m_cgal/` with **11/11 validator
    checks passing**.
  - 166,891 tets, 7,932 fault triangles, 28,596 vertices.
  - `algo3d=10` (HXT), `tube_radius=5000`, `res_f=1000`,
    `res_ff=20000`, `ramp_dist=30000`, `target_edge_m=1000`.
- **User concern:** *"are we sure the mesh size near the fault
  traces (must be 3D like a tube) are all within uniform mesh?"*

The headline answer is **NO** — direct measurement of the tet
edge distribution by distance from the fault shows the tube is
**not** uniform at `res_f`. The validator says PASS because its
bound (±50% of `res_f`) is too loose to detect non-uniformity.

I assumed at least 3 bugs and found 6 (2 CRITICAL, 2 MODERATE,
2 LOW). The CRITICAL findings directly answer the user's
question: the tube tet edges average **1383 m** vs target
**1000 m** (+38%), and only **6 % of in-tube tets are within
±15 % of `res_f`**. SEAS friction-zone resolution depends on
this number; current Phase-3 output ships with a 40 %-coarser
near-fault mesh than the plan specifies.

---

## Findings

### [R-401] [CRITICAL] [safs.geo + check_11_tube_uniformity] — In-tube tets are NOT uniform at `res_f`; mean = 1383 m vs target 1000 m, only 6 % within ±15 %

**Category:** BUG (silent acceptance — validator gate passes
because its tolerance is too loose)

**Description:**
Direct measurement of the post-Phase-3 mesh
(`output/two_crossing_2000m_cgal/output/safs_two_crossing_2000m.msh`)
by binning tet centroid distance to the nearest fault triangle:

```
distance from fault | n tets  | mean | p5   | p50  | p95  | max  | aspect_max
[    0,   500) m   |  17520 | 1140 |  788 | 1173 | 1377 | 1634 |  13.7
[  500,  1000) m   |  17586 | 1389 | 1150 | 1394 | 1619 | 2023 |  12.5
[ 1000,  2000) m   |  24527 | 1413 | 1193 | 1408 | 1645 | 2079 |   3.7
[ 2000,  3000) m   |  26316 | 1409 | 1197 | 1405 | 1631 | 1949 |   2.7
[ 3000,  4000) m   |  27336 | 1410 | 1202 | 1406 | 1629 | 1918 |   2.5
[ 4000,  5000) m   |  25914 | 1458 | 1224 | 1449 | 1720 | 2134 |   2.8
```

`tube_radius` is **5000 m** and `res_f` is **1000 m**.  The
plan and `safs.geo` specify *"distance ∈ [0, tube_radius] →
size = res_f (uniform tube)"*.  The actual mesh size in the
tube is **1383 m mean** (38 % coarser than `res_f`), and the
distribution is essentially CONSTANT at ~1400 m across all
distance bands inside the tube — there is no "uniform res_f
shell" at all.

Detailed in-tube quality:

```
In-tube (139199 tets):   mean=1383   sd=174   cv=0.13
   fraction within ±15% of res_f:  6.0%
   fraction within ±30% of res_f: 28.6%
   fraction within ±50% of res_f: 75.6%
```

The validator's `check_11_tube_uniformity` (validate_msh.py:594)
bounds are:
```python
if not (0.5 * res_f <= mean_e <= 1.5 * res_f):  # 500..1500 m
    return FAIL
if p95_e > 2.5 * res_f:                          # 2500 m
    return FAIL
```
1383 m is inside [500, 1500] and p95 1893 m is below 2500 m, so
the check PASSES.  But the check tolerates a 50 % deviation —
that is not a uniformity check, it is a "tube exists at all"
check.

**Why this matters for SEAS:**
- The cohesive zone scale `L_nuc = G * Dc / (sigma_n * (b - a))`
  is typically a few hundred metres; SEAS BP5 expects
  resolution ≤ `L_nuc / 2`.  At `res_f = 1000 m` the plan
  already has only marginal cohesive-zone resolution; at the
  actual mesh size of **1383 m** the resolution is materially
  worse than the plan promised.
- Time-step CFL is bounded by `min(edge) / c_s` for explicit
  schemes.  Mean edge 1383 m is ~38 % larger than 1000 m, but
  `min` 343 m (inner band) is half target.  CFL is constrained
  by the smallest edge → 343 m / `c_s` ≈ small dt; throughput
  dominated by inner-band slivers.  This is the worst-of-both:
  coarse resolution AND small time step.

**Trigger:** Always — every Phase-3 run with default flags.

**Actual behavior:** Validator reports 11/11 PASS, mesh ships
with non-uniform tube.  User-visible: high mean tet edge in
the tube zone.

**Expected behavior:** Tube tet edge mean within ±20 % of
`res_f`.  Either (a) HXT must be tuned to actually honour the
size field, or (b) the size-field input must be biased to
compensate for HXT's ~1.4× coarsening factor.

**Suggested fix (option (b) — size-field bias, simplest):**

HXT empirically produces tets ~1.4× the requested size.  Bias
the size field by setting `MeshSizeFactor = 0.7` and tightening
the Threshold field's `SizeMin`:

```diff
@@ mesh/safs.geo  — Field[2] block
 Field[2] = Threshold;
 Field[2].InField  = 1;
-Field[2].SizeMin  = res_f;
+Field[2].SizeMin  = res_f * 0.7;   // R-401 bias: HXT produces
+                                    // ~1.4× the requested size on
+                                    // open-surface input.  Bias
+                                    // the field down so actual
+                                    // post-mesh size is ~ res_f.
 Field[2].SizeMax  = res_ff;
@@ — at the bottom, before the Mesh.Algorithm3D line
+// R-401: global mesh-size factor bias.  HXT produces tets
+// ~1.4× the requested size near surface constraints; this
+// factor scales every field's output by 0.7 so the actual
+// near-fault tets land at ~ res_f.
+Mesh.MeshSizeFactor = 0.7;
```

A more principled fix is to use a **distance-aware sizing
field** (e.g., MathEval over the distance):
```
Field[2] = MathEval;
Field[2].F = Sprintf("%.3f * %.3f^min(F1/%.0f, 1)", res_f, res_ff/res_f, ramp_dist);
```
but this is more invasive.  The MeshSizeFactor bias is the
1-line stop-gap; cite this finding in the followup ticket for
proper distance-aware sizing.

**Tighten the validator** so this gap can't recur:

```diff
@@ mesh/validate_msh.py  check_11_tube_uniformity
-    if not (0.5 * res_f <= mean_e <= 1.5 * res_f):
+    # R-401: tightened from [0.5, 1.5] to [0.8, 1.25] so the
+    # check actually verifies "uniform-at-res_f", not just
+    # "tube exists".
+    if not (0.8 * res_f <= mean_e <= 1.25 * res_f):
         return CheckResult(
             "11_tube_uniformity", False,
-            f"mean tube tet edge {mean_e:.0f} m outside "
-            f"[{0.5*res_f:.0f}, {1.5*res_f:.0f}] (res_f={res_f}); "
+            f"mean tube tet edge {mean_e:.0f} m outside "
+            f"[{0.8*res_f:.0f}, {1.25*res_f:.0f}] (res_f={res_f}); "
             f"sizing field may not be reaching res_f inside the tube",
             metrics,
         )
+    # R-401: also assert at least 70 % of tube tets are within
+    # ±25 % of res_f (per-tet uniformity, not just the mean).
+    n_uniform = int(np.sum((edges >= 0.75 * res_f) & (edges <= 1.25 * res_f)))
+    frac_uniform = n_uniform / max(edges.size, 1)
+    metrics["frac_in_tube_within_25pct"] = frac_uniform
+    if frac_uniform < 0.70:
+        return CheckResult(
+            "11_tube_uniformity", False,
+            f"only {frac_uniform:.0%} of in-tube tets are within "
+            f"±25% of res_f (need >= 70%); the tube has the right "
+            f"mean but is internally non-uniform",
+            metrics,
+        )
```

**Test case:**
```python
def test_R401_tube_uniformity_actually_checks():
    """Synthetic mesh: build a fault + tet mesh where mean tube
    tet edge = 1.4 × res_f.  Old (loose) bound: PASS.  New
    (tight) bound: FAIL with informative message."""
    # Build a synthetic mesh with known tube_uniformity stats.
    points, tets, fault_tris = make_uniform_tube_mesh(
        res_f=1000, mean_tube_edge=1400, tube_radius=5000)
    # ... call check_11 ...
    result = check_11_tube_uniformity(points, tets, ...,
                                       res_f=1000.0,
                                       tube_radius=5000.0)
    assert not result.passed
    assert "outside [800, 1250]" in result.message
```

---

### [R-402] [CRITICAL] [safs.geo + check_10_tet_quality] — Inner-band (d ≤ 500 m) tets have aspect ratio up to 13.7; check_10 PASS-with-warning hides this

**Category:** BUG (silent acceptance of bad-quality tets near
the fault — exactly where SEAS friction is most sensitive)

**Description:**
The inner band (d ≤ 500 m from any fault triangle) has 17 520
tets with:
- mean edge 1140 m,
- **aspect_max = 13.7** (max-edge / min-edge per tet).

These are needle-shaped tets sandwiched between adjacent fault
triangles whose edges are tightly constrained.  At aspect 13.7,
basis functions are ill-conditioned and SEAS DG fluxes will
produce numerical noise.

The validator's `check_10_tet_quality` reports:
```
[PASS] 10_tet_quality: min gamma=0.0378, mean gamma=0.796
                       — WARNING: 1 sliver tet(s) below 0.05
        gamma_min: 0.0378
        n_gamma_below_0p05: 1
        frac_gamma_below_0p05: 5.99e-06
```

`gamma = 0.0378` is **catastrophically low** (the SEAS
acceptance threshold per debug history is 0.05).  The check
"passes with warning" because the FRACTION is small (1 / 166 891).
But **a single gamma=0.038 tet near a fault is enough to
inject ~10× more noise than a gamma=0.5 tet** in the local
flux assembly.  The pass-with-warning is inappropriate.

Looking at the code:
```python
# validate_msh.py — check_10_tet_quality (inferred from output)
if frac_gamma_below_0p05 < some_threshold:
    return PASS_WITH_WARNING
```
Without seeing the exact threshold, the tolerance is too loose
for SEAS's needs.

**Trigger:** Any Phase-3 run.  Always.

**Actual behavior:** Pass with warning at `gamma=0.038`.

**Expected behavior:** FAIL if any tet within 1 km of the
fault has gamma < 0.05.  Sliver tets in the far field (where
they don't affect the fault) are tolerable; sliver tets at
the fault are not.

**Suggested fix:** Promote sliver tets within `tube_radius` of
the fault to fatal, while tolerating far-field slivers:

```diff
@@ mesh/validate_msh.py  check_10_tet_quality
+    # R-402: split the gamma check by distance from fault.
+    # Sliver tets in the far field are tolerable; sliver tets
+    # within the tube_radius are fatal because they live in the
+    # SEAS friction-active zone.
+    from scipy.spatial import cKDTree
+    fault_centroids = points[fault_tris].mean(axis=1)
+    tree = cKDTree(fault_centroids)
+    tet_centroids = points[tets].mean(axis=1)
+    d_to_fault, _ = tree.query(tet_centroids, k=1)
+    near_fault_slivers = int(np.sum(
+        (d_to_fault <= tube_radius) & (gamma < 0.05)))
+    metrics["near_fault_slivers"] = near_fault_slivers
+    if near_fault_slivers > 0:
+        return CheckResult(
+            "10_tet_quality", False,
+            f"{near_fault_slivers} sliver tet(s) within "
+            f"tube_radius={tube_radius} m of the fault have "
+            f"gamma < 0.05; SEAS DG flux assembly will produce "
+            f"numerical noise on these.",
+            metrics,
+        )
```

(`fault_tris` and `tube_radius` need to be threaded through to
`check_10`; they're already available in the driver.)

**Test case:**
```python
def test_R402_near_fault_sliver_is_fatal():
    """A single gamma=0.038 sliver within 500 m of the fault
    must FAIL the check, not pass-with-warning."""
    # Build a mesh with one sliver at d=400 m, one at d=20 km.
    # Old: pass-with-warning.  New: FAIL because of the d=400 sliver.
    ...
    assert not result.passed
    assert "within tube_radius" in result.message
```

---

### [R-403] [MODERATE] [validate_msh.py:check_10_tet_quality + safs.geo] — `Mesh.OptimizeNetgen` is left at default; HXT has a follow-up optimizer pass that improves quality but isn't aggressively configured

**Category:** ASSUMPTION (configuration left at default)

**Description:**
`safs.geo` sets:
```
Mesh.Optimize           = 1;
Mesh.OptimizeNetgen     = 1;
Mesh.OptimizeThreshold  = 0.3;
```

`OptimizeThreshold = 0.3` means the optimizer only attempts to
improve tets below 0.3 quality.  The `gamma=0.038` outlier IS
below 0.3, so the optimizer should have tried to fix it — but
clearly didn't.

This is because `Mesh.OptimizeNetgen` operates on the LOCAL
patch and can't move CONSTRAINED vertices (the fault
triangulation).  A 13.7-aspect tet that touches a fault
triangle has one face whose three vertices are pinned; the
optimizer can only move the FOURTH vertex.  Often this is
insufficient.

The remediation isn't in the optimizer — it's in the size
field producing better-shaped tets up front.  R-401's
`MeshSizeFactor = 0.7` bias should help indirectly by giving
the optimizer more vertex-space to work in.

**Suggested fix:** Pass `Mesh.QualityType = 2` (gamma — already
set) AND `Mesh.OptimizeThreshold = 0.5` (more aggressive):

```diff
@@ mesh/safs.geo
 Mesh.Optimize           = 1;
 Mesh.OptimizeNetgen     = 1;
-Mesh.OptimizeThreshold  = 0.3;
+Mesh.OptimizeThreshold  = 0.5;   // R-403: cast a wider net for
+                                  // sliver removal; HXT's
+                                  // post-mesh optimizer can fix
+                                  // tets up to gamma=0.5 — the
+                                  // 0.3 default leaves slivers
+                                  // at gamma=0.3..0.5 untouched.
 Mesh.QualityType        = 2;
```

---

### [R-404] [MODERATE] [safs.geo:106 — Field[1].Sampling = 100] — Distance-field sampling is per-surface; with 7 932 fault triangles the sample density is fine, but the per-surface=100 cap is fragile

**Category:** ASSUMPTION

**Description:**
```
Field[1] = Distance;
Field[1].SurfacesList = {fault_surfs[]};
Field[1].Sampling     = 100;
```

`Field[1].Sampling = 100` means the Distance field samples each
fault SURFACE 100 times to compute the distance.  With 7 932
fault triangles and Sampling=100, the effective sample density
is high enough.  But:
- If a future run has FEWER fault triangles (e.g. a coarser CFM
  resolution), Sampling=100 might miss fault detail and produce
  patchy distance fields → patchy tube refinement.
- The sample count is a per-SURFACE constant — it doesn't scale
  with the surface area or triangle count.

**Trigger:** Coarse fault input (e.g., 500 m CFM resolution
producing < 100 triangles per fault).

**Suggested fix:** Make sampling proportional to triangle count
in the surface, with a floor:

```diff
@@ mesh/safs.geo
 Field[1] = Distance;
 Field[1].SurfacesList = {fault_surfs[]};
-Field[1].Sampling     = 100;
+// R-404: 100-per-surface is fine when each fault has 100s of
+// triangles, but degrades on coarse input.  Increase to 200
+// for safety.  Cost is negligible — distance field is computed
+// once.
+Field[1].Sampling     = 200;
```

A more principled fix would compute Sampling = max(100,
total_fault_triangles / N_surfaces) at driver time and inject
via `-setnumber`.  Defer to a follow-up.

---

### [R-405] [LOW] [run_two_crossing_2000m_cgal.sh:106 — exit propagation] — Validator FAIL doesn't surface to the script's exit code

**Category:** QUALITY

**Description:**
```bash
echo "==> 6/7 validate_msh"
python validate_msh.py \
    --msh ... \
    --report ...
```

`validate_msh.py` writes the report and may exit non-zero on
FAIL.  `set -euo pipefail` propagates the exit code.  But the
script doesn't print the report to stderr — the user has to
manually `cat` the report file to see which check failed.  In
CI / automation, this is fine; for interactive use the path-only
output is opaque.

**Suggested fix:**

```diff
@@ mesh/run_two_crossing_2000m_cgal.sh — step 6
 echo "==> 6/7 validate_msh"
 python validate_msh.py \
     --msh "$OUTDIR/output/safs_two_crossing_${RES}m.msh" \
     ... \
-    --report "$OUTDIR/output/validation_report.txt"
+    --report "$OUTDIR/output/validation_report.txt" || {
+        echo "" >&2
+        echo "validate_msh.py FAILED — full report follows:" >&2
+        cat "$OUTDIR/output/validation_report.txt" >&2
+        exit 1
+    }
+# Even on success, print the summary tail so the user sees
+# pass/fail counts and any per-check warnings.
+tail -3 "$OUTDIR/output/validation_report.txt"
```

---

### [R-406] [LOW] [generate_safs_mesh.py + safs.geo — `Mesh.MeshSizeMin` global floor] — Global `MeshSizeMin = res_f` may FORCE the size field below res_f near the fault but without strictly enforcing res_f as a floor

**Category:** ASSUMPTION

**Description:**
```
Mesh.MeshSizeMin                = res_f;
Mesh.MeshSizeMax                = res_ff;
```

`Mesh.MeshSizeMin` is a GLOBAL minimum — gmsh refuses to
produce edges shorter than `res_f`.  But the actual mesh has
edges as short as **343 m** within 1 km of the fault (per the
Phase-3 measurement at the top of this review).  So either
gmsh is ignoring `MeshSizeMin`, or HXT's constraint-honoring
code creates short edges to match the fault triangulation
without consulting the size field.

The latter is the documented HXT behaviour: surface
constraints take precedence over size fields.  So
`MeshSizeMin = res_f` is informational only, not enforced.

This isn't a bug per se, but the user expecting "all edges in
the tube are ~ res_f" may be confused by edges down to 343 m.
Document this:

```diff
@@ mesh/safs.geo  — comment near MeshSizeMin
+// R-406: MeshSizeMin is a soft floor, not a hard constraint.
+// HXT honours fault-triangle constraints (which on Phase-2
+// output have edges as low as ~ res_f / 3 = 333 m near
+// polyline kinks) before the global size-field floor, so
+// in-tube tets adjacent to those short fault edges may be
+// shorter than MeshSizeMin.  This is expected; see R-406.
 Mesh.MeshSizeMin                = res_f;
 Mesh.MeshSizeMax                = res_ff;
```

---

## Summary

- Critical issues: **2** (R-401 tube non-uniform at +38 % vs target,
  R-402 sliver tet near fault accepted)
- Moderate issues: **2** (R-403 optimizer threshold, R-404
  distance-field sampling)
- Low issues: **2** (R-405 validator exit propagation, R-406
  MeshSizeMin documentation)
- Plan compliance: **PARTIAL**.
  - Phase 3 acceptance "9/11 must-pass": ✅ all 11 PASS.
  - Phase 3 acceptance "stretch — uniform near-fault shell at
    res_f": ❌ tube mean is 1383 m, not 1000 m.
  - The plan's strict reading would call this a **Phase-3 stretch
    failure**, but the validator's loose tolerance hides it.
- Verdict: **PASS WITH FIXES**. Phase 3 produces a SEAS-runnable
  mesh (R-401 doesn't break the simulation, just makes it less
  accurate than the plan promised).  The user is correct to push
  back — fix R-401 + R-402 before considering the tube uniform.

## Direct answer to the user's question

> *"are we sure the mesh size near the fault traces (must be 3D
> like a tube) are all within uniform mesh, this is very
> important to ensure the code can run correctly?"*

**Currently NO.** Direct measurement of the post-Phase-3 mesh:

| Tube band | n tets | mean tet edge | vs `res_f`=1000 |
|---|---:|---:|---:|
| 0–500 m from fault    | 17 520 | **1140** | +14 % |
| 500–1 000 m            | 17 586 | **1389** | +39 % |
| 1 000–2 000 m          | 24 527 | **1413** | +41 % |
| 2 000–3 000 m          | 26 316 | **1409** | +41 % |
| 3 000–4 000 m          | 27 336 | **1410** | +41 % |
| 4 000–5 000 m          | 25 914 | **1458** | +46 % |
| **In-tube total**     |139 199 | **1383** | **+38 %** |

Only **6 %** of in-tube tets are within ±15 % of `res_f`.  The
tube is "approximately uniform" at 1.4 × `res_f`, not at `res_f`.

The validator says PASS because its bound (`[0.5, 1.5] × res_f`)
is too loose — 1383 m is inside [500, 1500].  R-401's fix
applies a `Mesh.MeshSizeFactor = 0.7` bias to drop the actual
tube mean to ~1000 m, and tightens `check_11`'s bound to
[0.8, 1.25] × `res_f` so future runs can't silently regress
beyond ±25 %.

After R-401 + R-402 land, expected post-Phase-3 numbers:
- in-tube mean tet edge: 1000 m ± 100 m,
- 70 %+ of in-tube tets within ±25 % of `res_f`,
- 0 sliver tets (gamma < 0.05) within `tube_radius` of any fault.

For SEAS this matters because:
- the cohesive zone is `O(L_nuc) ≈ 200–500 m`; current
  tube resolution leaves `< 1 element per cohesive zone`,
- aspect-ratio-13.7 tets at the fault inject DG flux noise
  into the friction-rate equation,
- CFL is bounded by min-edge (343 m) so the simulation pays
  for the slivers but doesn't get the benefit of uniform
  near-fault resolution.

## Suggested fix order

1. **R-401** (one-line `Mesh.MeshSizeFactor = 0.7` in safs.geo +
   tighter `check_11` bound) — directly answers the user's
   question.  Land first.
2. **R-402** (per-band sliver gate in `check_10`) — catches
   the gamma=0.038 tet that R-401 may not eliminate.
3. **R-403** (raise `OptimizeThreshold` from 0.3 to 0.5) —
   trivial; lands together with R-401.
4. **R-404, R-405, R-406** — quality cleanups.

Total effort: ~30 LOC across `safs.geo`, `validate_msh.py`,
and `run_two_crossing_2000m_cgal.sh`.

## Unreviewed Areas

- **The validator's exact `check_10_tet_quality` threshold for
  pass-with-warning vs FAIL.** I read the output line but did
  not trace the exact bound.  R-402's fix bypasses this —
  add a near-fault sliver gate that's stricter than the global
  one, regardless of the global threshold.
- **`Mesh.MeshSizeFactor` interaction with Distance field.**
  R-401 assumes the factor scales the Distance field's output
  uniformly.  It does (per gmsh docs), but verify with one
  small build before committing the bias ratio of 0.7.
- **HXT's `Mesh.AlgorithmSwitchOnFailure`.**  Not set.  If HXT
  rejects the input on a future run (corefine output that
  R-101 onward isn't quite catching), the script would error
  rather than fall back.  Acceptable for now; document.
- **Cross-fault polyline-conformity inside the bulk mesh.** I
  didn't verify that fault triangles from A and B share the
  *same vertex* on the polyline in the bulk mesh's `points`
  array.  The combine step's `n_dropped_duplicate: 0` is
  consistent with welded vertices, but a direct vertex-id
  check would close this loop.

---

## Addendum (2026-04-30): Phase 4 + Phase 5 verification

The following observations follow the implementer's report
on Phases 4–5.  Per the user's request, they are appended here
rather than written to a separate review file.

### Phase 4 — verified COMPLETE

All 6 acceptance items the implementer reported were checked
directly on disk and at runtime:

| Claim | Verification |
|---|---|
| `mesh/run_two_crossing_2000m.sh` is an 8-line wrapper | `wc -l` = 8 lines (1 shebang + 5 comment + 1 `export` + 1 `exec`).  Body is `exec "$(dirname "$0")/run_two_crossing_2000m_cgal.sh" "$@"`. ✓ |
| `pytest mesh/tests/` skips legacy by default | Live run: `42 deselected, 1 warning in 0.11s`. ✓ |
| `pytest -m legacy` runs them | Live run: `42 passed, 1 warning in 36.82s`. ✓ |
| `import conformalize_faults` emits exactly one DeprecationWarning | Live run with `warnings.simplefilter('always')` + `catch_warnings(record=True)` → DeprecationWarning count = 1, message starts with "conformalize_faults.py is the legacy …". ✓ |
| `mesh/tests/pytest.ini` registers the marker + sets `addopts` | File present; `markers = legacy: …` and `addopts = -m "not legacy"`. ✓ |
| `mesh/tests/README.md` documents the 90-day retention clock | File present; clock-start, reset condition, removal eligibility all documented. ✓ |
| `PLAN_multifault_intersections.md` SUPERSEDED banner | First non-blank line: `> **Status: SUPERSEDED by PLAN_cgal_corefine.md (Phase 4, 2026-04-29).**`. ✓ |
| `mesh/conformalize_faults.py` deprecation hook | `warnings.warn(..., DeprecationWarning)` at line 41–45. ✓ |

**Phase 4 verdict: PASS.** No new findings against the Phase-4
changes.

### Phase 5 — wrapper landed; canonical script exits 2 on real
all-8 input as designed; the user-experience around the
"destructive change" needs review

The implementer correctly mirrored the Phase-4 P-009 wrapper
pattern for Phase 5: replaced `run_all8_2000m.sh` with an
8-line `exec` wrapper that delegates to a new canonical
`run_all8_2000m_cgal.sh`.  On Mill Creek × SBMT-SAF the
canonical CGAL pipeline succeeds in 11/11 (Phase 3); on the
all-8 input it exits 2 at step 3 with
`FAIL_SELF_INTERSECT` on `safs_mult_banning`.

This matches the plan's documented at-risk gate (`Pinto
Mountain` and `Banning standalone` are CFM-source-defective
per `REVIEW_phase1.md` and the all-8 README).  The
implementer's "Option 2" (keep the wrapper, treat exit-2 as
the documented blocker) is **the correct interpretation of
the plan**.

But two issues with the way it landed are worth flagging.

#### [R-401a] [MODERATE] [run_all8_2000m_cgal.sh:83 — `rm -rf "$OUTDIR"` is destructive and inconsistent with the Phase-3 script]

**Category:** EDGE_CASE (destructive operation; no equivalent
in the Phase-3 sibling script)

**Description:**
`run_all8_2000m_cgal.sh:83` does:

```bash
OUTDIR="output/all8_${RES}m${OUTPUT_SUFFIX}"
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR/stl_raw" "$OUTDIR/stl_conformal" "$OUTDIR/output"
```

The two-fault sibling `run_two_crossing_2000m_cgal.sh:71` does
NOT have the `rm -rf` step — it just calls `mkdir -p` and
overwrites files in place.  The all-8 script's `rm -rf`
unconditionally wipes the output directory, which:

- **Loses the prior all-8 build's diagnostics** if the user
  was relying on them.  Concretely, the user previously had
  `output/all8_2000m/output/safs_all8_2000m.msh` (1.3 MB,
  dated Apr 25) — that mesh is preserved because the wrapper
  defaults `ALL8_OUTPUT_SUFFIX=""` and the `rm -rf` path
  resolves to `output/all8_2000m/`, but only after the
  wrapper invocation runs to step 3 and exits 2.  The
  pre-step-3 file removal happens before any failure can
  abort the run, so the prior outputs are wiped on every
  invocation.
- **Is silent.**  No warning, no `--force` flag, no prompt.

The two scripts should behave the same way on this point.

**Trigger:** Any invocation of the wrapper with
`ALL8_OUTPUT_SUFFIX` resolving to a path that contains
already-good artefacts (the prior `output/all8_2000m/`
directory).

**Verification on this machine:** I did NOT re-run the
wrapper (deferred to user direction per the implementer's
report).  But the script lines are inspected directly; the
`rm -rf` is unambiguous.

**Suggested fix:**

```diff
@@ run_all8_2000m_cgal.sh:82–84
 OUTDIR="output/all8_${RES}m${OUTPUT_SUFFIX}"
-rm -rf "$OUTDIR"
-mkdir -p "$OUTDIR/stl_raw" "$OUTDIR/stl_conformal" "$OUTDIR/output"
+# R-401a: do NOT `rm -rf` the output directory unconditionally.
+# The Phase-3 sibling (`run_two_crossing_2000m_cgal.sh`) just
+# `mkdir -p`s and overwrites; we do the same here for consistency.
+# Users who want a clean slate should `rm -rf` manually before
+# running.
+mkdir -p "$OUTDIR/stl_raw" "$OUTDIR/stl_conformal" "$OUTDIR/output"
```

**Test case:**
```bash
# Confirm `rm -rf` is gone.  After the fix, this should produce
# an empty diff between the two scripts on the directory-creation
# logic.
diff <(grep -E 'rm -rf|mkdir -p' run_all8_2000m_cgal.sh) \
     <(grep -E 'rm -rf|mkdir -p' run_two_crossing_2000m_cgal.sh)
# Expected: only path differences; no `rm -rf` mismatch.
```

#### [R-402a] [MODERATE] [run_all8_2000m_cgal.sh — no preflight on at-risk faults]

**Category:** QUALITY (UX — diagnostic friction)

**Description:**
The plan documents that `safs_mult_banning` and
`safs_pmfz_pinto` are CFM-source-defective and will fail the
self-intersection gate inside `corefine_faults`.  The script
includes them in `FAULTS=(...)` unconditionally, runs
ts_to_stl on them (~30 s), then fails at corefine_faults step
3 with a CGAL exception.

The user-visible failure is:
```
[corefine_faults] manifold gate FAILED on fault safs_mult_banning
after corefine with safs_coav_missioncreek: FAIL_SELF_INTERSECT
```

This is correct, but takes ~30+ seconds to surface.  The
information about WHICH faults are at-risk is in the plan and
in `REVIEW_phase1.md` but NOT in the script's failure path.

**Suggested fix:** Add a preflight that lists the at-risk
faults included in this run and warns the user to expect
exit-2 if the upstream CFM source isn't repaired:

```diff
@@ run_all8_2000m_cgal.sh — between FAULTS array and step 1
+# Preflight: warn if at-risk faults are included.  Per
+# REVIEW_phase1.md and the all-8 README, these CFM faults have
+# self-intersections in the upstream .ts source; the corefine
+# pipeline will exit 2 at step 3 unless the source is repaired.
+AT_RISK=(safs_mult_banning safs_pmfz_pinto)
+included_at_risk=()
+for ar in "${AT_RISK[@]}"; do
+    for f in "${FAULTS[@]}"; do
+        if [ "$f" = "$ar" ]; then included_at_risk+=( "$ar" ); fi
+    done
+done
+if [ ${#included_at_risk[@]} -gt 0 ]; then
+    echo "WARN: at-risk fault(s) included: ${included_at_risk[*]}" >&2
+    echo "      these have self-intersections in the upstream CFM" >&2
+    echo "      .ts source; corefine_faults is expected to exit 2" >&2
+    echo "      at step 3 unless the source is repaired upstream." >&2
+    echo "      To skip them, edit the FAULTS array near the top" >&2
+    echo "      of this script.  Continuing anyway in 3 s..." >&2
+    sleep 3
+fi
```

This makes the failure mode legible and actionable from the
script's output, not just from `REVIEW_phase1.md`.

#### Re. the user's open question (Option 1 vs Option 2)

The implementer is asking whether to keep the wrapper or
revert to the legacy Frontal-Delaunay path.  Two facts the
user should weigh:

1. **The legacy `output/all8_2000m/output/safs_all8_2000m.msh`
   is still on disk** (1.3 MB, dated Apr 25 — not yet
   overwritten by the new wrapper because the new wrapper
   exits 2 at step 3, before generating a new bulk mesh
   artefact would happen at step 4).  The wrapper change has
   removed the *script* that produces such a mesh, but the
   prior artefact is preserved unless the user invokes the
   wrapper without R-401a's `rm -rf` removal (the
   `safs_all8_2000m.msh` is INSIDE `$OUTDIR`, so the `rm -rf`
   in the canonical script WOULD wipe it on first invocation).
2. **The legacy mesh failed `validate_msh.py:check_5` 100 %**
   per `REVIEW.md` R-001 (every fault triangle was detached
   from the tet mesh).  The previous "working" path was
   producing a SEAS-unusable mesh — the script's silent
   pass-through of bad output is what the Phase-1 review
   identified as the original bug worth replacing CGAL for.

   So the choice is:
   - **Option 1 (revert the wrapper)** preserves the legacy
     SCRIPT but the legacy script produced a SEAS-unusable
     mesh (n_offending = 4917, every fault triangle wrong).
   - **Option 2 (keep the wrapper)** matches the plan's P-009
     pattern AND is honest about the all-8 input's CFM-source
     defects.  The user gains a fail-fast diagnostic; they
     lose the silent-pass-through of garbage output.

   Recommendation: **Option 2**, with R-401a applied so the
   prior `safs_all8_2000m.msh` artefact survives the
   exit-2 path, and R-402a applied so the failure mode is
   legible.  The user can keep the prior `.msh` for visual
   reference in ParaView (with the caveat that it's
   SEAS-unusable).  When the upstream CFM defects are
   repaired, the canonical CGAL pipeline produces a
   correct all-8 mesh without script changes.

   If the user disagrees and wants Option 1 (legacy script
   restored), that's a one-command revert of
   `mesh/run_all8_2000m.sh` from git history; the canonical
   `_cgal` script can stay in the tree for the user-when-ready
   path.

**Phase 5 verdict (corrective): PASS WITH FIXES.**  The plan
compliance is correct; R-401a + R-402a make the wrapper
behaviour match the user's stated need ("retain a working
8-fault build path or know exactly why it failed") without
re-introducing the silent-pass-through bug Option 1 would
require.
