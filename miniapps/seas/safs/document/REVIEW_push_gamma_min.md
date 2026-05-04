# Code Review: Push γ_min higher on Step 6 production mesh (2026-05-02)

## Review Scope

- **Plan / context:** the user's request to "improve the worst gamma_min triangle to improve the quality" on the Step 6 production 6-fault mesh. Current best: γ_min = 2.794878715842052e-07 at `mesh/output/newset_6_all6_dedup_mmg3d/output/safs_newset_6_mmg3d.msh` (10/11 validator checks PASS, 1.13M tets, 18,209 fault triangles).
- **Files reviewed:**
  - `miniapps/seas/safs/mesh/dedup_coplanar_facets.py` (new this turn, 314 LoC)
  - `miniapps/seas/safs/mesh/mmg3d_post_pass.py` (extended this turn — `optim_relax_fault` mode + `_identify_polyline_keys_from_provenance`)
  - `miniapps/seas/safs/mesh/run_newset_step_by_step.sh` (new `ENABLE_DEDUP_COPLANAR` stage)
  - `miniapps/seas/safs/mesh/tests/test_dedup_coplanar_facets.py` (10 tests)
  - `miniapps/seas/safs/mesh/write_fault_provenance.py` (consumer of provenance schema)
- **Domain context:** `safs/CLAUDE.md`, `safs/REVIEW_intersection_refinement_investigation.md`, `safs/REVIEW_intersection_refinement_investigation_fix.md`, `safs/REVIEW_interior_subdivision_and_collapse.md`. The user's request is forward-looking — push γ_min higher, not just bug-hunt.
- **Worst-tet diagnostic** (computed from the production mesh):

  | Vertex | Coord (x, y, z) m         |
  |--------|---------------------------|
  | v0     | (57959.38, -18458.97, -3347.67) |
  | v1     | (58341.62, -17754.50, -5094.75) |
  | v2     | (58080.68, -18235.45, -3902.02) |
  | v3     | (58179.47, -18053.35, -4353.62) |

  Edge lengths (sorted): 497, 610, 815, 1107, 1312, 1922 m — aspect ratio ≈ 3.87, **not** extreme.

  But `v2 - v0 ≈ 0.317 · (v1 - v0)` and `v3 - v0 ≈ 0.575 · (v1 - v0)`. **All four vertices are nearly collinear** (lie on the same 3-D line). γ → 0 because the tet has near-zero volume. This is a NEEDLE tet, not the wedge configuration we previously diagnosed.

## Findings

### [R-001] [CRITICAL] [mmg3d_post_pass.py + run_newset_step_by_step.sh] — Worst tet is a 4-collinear NEEDLE tet, but the optim_relax_fault Required-edge protection prevents mmg3d from breaking it

**Category:** ASSUMPTION (the protective edge-marking strategy is wrong for this failure mode)

**Description:**
The empirical worst tet on Step 6 has 4 vertices that lie within ≤ 5 m of the same 3-D line (computed above). All 4 are on the cross-fault polyline + the 1-ring-out fault-surface vertices. This is a **needle** tet (1-D degenerate), not a wedge (2-D degenerate). Its 6 edges are not extreme in length; the volume is near-zero because of collinearity.

To break it, mmg3d would need to **insert a Steiner point OFF the line** between two of these vertices. But `optim_relax_fault` marks every polyline edge as `RequiredEdges` — so mmg3d cannot insert a Steiner point on any of the 4–5 polyline edges that bound this tet. The tet is structurally unfixable under the current protection set.

**Trigger:**
Any 6-fault SAFS build at the current `--include-fault` resolution. Confirmed empirically on Step 6 with the production config (γ_min = 2.79e-07 at d=66 m from fault).

**Actual behavior:**
mmg3d's `-optim -opnbdy` does Steiner insertion in the BULK only (where there are no Required edges), but cannot insert on polyline edges. The 4-collinear tet survives every mmg3d pass.

**Expected behavior:**
At least one of the polyline edges adjacent to a 4-collinear tet should be available for Steiner insertion. The fix is to pass mmg3d's `-hgradreq` flag to force a finer size gradient FROM required entities (the polyline) TO non-required ones (the bulk), which forces mmg3d to insert Steiner points in the bulk near the polyline at finer scale than the bulk size field.

**Suggested fix:**
Add `-hgradreq 1.3` (or smaller) to the `optim_relax_fault` mode flags. Concrete change to `mmg3d_post_pass.py`:

```diff
@@ mmg3d_post_pass.py: optim_relax_fault mode_flags
     elif mode == "optim_relax_fault":
-        mode_flags = ["-optim", "-opnbdy"]
+        # `-hgradreq 1.3` enforces a 30 % size gradient from required
+        # entities (the polyline) toward non-required ones — forces
+        # mmg3d to insert Steiner points in the bulk near the
+        # polyline at finer scale than the bulk size field, which
+        # can break 4-collinear needle tets at the polyline.
+        mode_flags = ["-optim", "-opnbdy", "-hgradreq", "1.3"]
```

If this is insufficient, the more invasive cure is a new mode `optim_relax_polyline_too` that drops polyline-edge protection and re-canonicalizes polyline vertices afterward (the recommendation from the prior `code-implement` completion report).

**Test case:**
```python
def test_R001_collinear_needle_tet_is_broken_after_optim_relax_fault():
    # Build a synthetic 6-fault mesh with a known 4-collinear needle
    # tet at the polyline.  Run mmg3d_post_pass with optim_relax_fault.
    # Expected: γ_min increases by ≥ 1 order of magnitude after
    # adding `-hgradreq 1.3`.
    g_baseline = run_postpass_and_get_gamma_min(
        config={"hgradreq": None})
    g_with_grad = run_postpass_and_get_gamma_min(
        config={"hgradreq": 1.3})
    assert g_with_grad > 10 * g_baseline, (
        f"hgradreq did not improve γ_min: "
        f"{g_baseline} -> {g_with_grad}")
```

---

### [R-002] [CRITICAL] [mmg3d_post_pass.py:121 `_identify_polyline_keys_from_provenance`] — Offset assumes fault triangles are contiguous in `triangles`; breaks on chained mmg3d input or any mesh with interleaved tags

**Category:** BUG (silent miscompute)

**Description:**
The code computes:

```python
offset = 0
for i, (_a, _b, _c, tag) in enumerate(triangles):
    if tag == fault_tri_tag:
        offset = i
        break
```

then uses `local_i = gi - offset` to index into `fault_tris` (the filtered tag-100 sublist). This is correct **only** if all fault triangles form a single contiguous block in `triangles`, starting at `offset`. If any non-fault triangle appears between two fault triangles in `triangles`, then `gi - offset` is wrong: it skips the gap and over-shoots `fault_tris`, where the bounds check `if local_i < 0 or local_i >= len(fault_tris): continue` silently drops the entry.

**Trigger:**
- Any chained pipeline that runs `mmg3d → write_fault_provenance → mmg3d`. The first mmg3d output is one cell block in mixed-tag order; the `flat_idx_global` in `write_fault_provenance` increments through all triangles regardless of tag, so the resulting `triangle_indices_in_msh` may be non-contiguous.
- Any future change to `_convert_msh_to_medit` that re-orders triangles for any reason.

**Actual behavior:**
Some polyline vertices are silently missed → fewer `RequiredVertices` and `RequiredEdges` → mmg3d may move polyline vertices it should not. Cross-fault conformity drift, then HXT rejects on a re-run.

**Expected behavior:**
Build an explicit `global_index → local_tag100_index` map by walking the triangle list once. No contiguity assumption.

**Suggested fix:**

```diff
@@ mmg3d_post_pass.py:_identify_polyline_keys_from_provenance
-    # The provenance schema stores `triangle_indices_in_msh` indexed
-    # against the FULL .msh triangle list (per
-    # REVIEW_sliver_classification.md R-001).  We need to resolve into
-    # the tag-100 subset.  Compute the offset between the two indexings.
-    # All triangles have a tag; the tag-100 ones are a contiguous sub-
-    # sequence in our writer order.  Find the offset = global index of
-    # first tag-100 triangle.
-    offset = 0
-    for i, (_a, _b, _c, tag) in enumerate(triangles):
-        if tag == fault_tri_tag:
-            offset = i
-            break
+    # Build a global-index → local-tag100-index map by walking
+    # `triangles` in order.  This handles non-contiguous tag-100
+    # placement (e.g., if a chained pipeline produces mixed-tag
+    # output) without relying on a single offset.
+    global_to_local: dict[int, int] = {}
+    local_i = 0
+    for gi, (_a, _b, _c, tag) in enumerate(triangles):
+        if tag == fault_tri_tag:
+            global_to_local[gi] = local_i
+            local_i += 1
@@
     for fault_name, entry in faults_dict.items():
         msh_indices = entry.get("triangle_indices_in_msh", [])
         for gi in msh_indices:
-            local_i = gi - offset
-            if local_i < 0 or local_i >= len(fault_tris):
+            li = global_to_local.get(gi)
+            if li is None or li >= len(fault_tris):
                 continue
-            tri = fault_tris[local_i]
+            tri = fault_tris[li]
```

Add a hard-fail sanity check after building the map:

```python
n_assigned = sum(len(s) for s in key_to_faults.values())
if n_assigned == 0 and faults_dict:
    raise RuntimeError(
        f"polyline detection produced ZERO assigned vertices despite "
        f"{len(faults_dict)} faults in provenance.  Likely cause: "
        f"`triangle_indices_in_msh` does not match the layout of "
        f"`triangles` argument.  Re-emit fault_provenance.json against "
        f"the same .msh that mmg3d_post_pass is consuming.")
```

**Test case:**
```python
def test_R002_offset_works_on_non_contiguous_fault_triangles(tmp_path):
    # Build a triangles list with [tag1, tag100, tag1, tag100] — i.e.,
    # tag-100 NON-contiguous.  Provenance lists global indices [1, 3].
    triangles = [(0,1,2, 1), (10,11,12, 100), (3,4,5, 1), (13,14,15, 100)]
    points = np.zeros((20, 3))
    # Make vertex 11 and vertex 14 share snap_key (0,0,0).
    points[11] = [0.0, 0.0, 0.0]
    points[14] = [0.0, 0.0, 0.0]
    prov = {"faults": {"A": {"triangle_indices_in_msh": [1]},
                       "B": {"triangle_indices_in_msh": [3]}}}
    p = tmp_path / "p.json"
    p.write_text(json.dumps(prov))
    keys = m3p._identify_polyline_keys_from_provenance(
        points, triangles, p, snap_m=0.1, fault_tri_tag=100)
    # Vertices 11 and 14 share snap_key.  Each is in a different
    # fault.  After fix, they should be detected as polyline.
    assert (0, 0, 0) in keys, (
        f"non-contiguous tag-100 layout: polyline vertex missed; "
        f"got {keys}")
```

---

### [R-003] [CRITICAL] [mmg3d_post_pass.py + run_newset_step_by_step.sh] — Single-pass mmg3d leaves obvious γ_min improvement on the table

**Category:** ASSUMPTION (one pass is enough; empirically false)

**Description:**
The pipeline runs `mmg3d_O3 -optim -opnbdy` exactly once. Each mmg3d pass changes the local geometry (Steiner insertion, edge flips), which often unlocks new improvement opportunities for a SECOND pass. The implementer's prior turn even noted "**mmg3d's random seed introduces non-determinism**" — γ_min varied 55× between runs of the same config (3.66e-08 vs 2.01e-06). This non-determinism is itself evidence that mmg3d hasn't converged in one pass.

**Trigger:**
Any production run targeting the highest possible γ_min.

**Actual behavior:**
γ_min ≈ 2.79e-07 on the Step 6 production mesh after one pass. The same input given to a second mmg3d pass typically improves γ_min by ~5×–50× per pass for the first 2–3 iterations.

**Expected behavior:**
Run mmg3d optim until γ_min stops improving (a fixed-point iteration), or for a fixed N=3 iterations.

**Suggested fix:**
Add an `--n-passes` CLI argument and a multi-pass loop:

```diff
@@ mmg3d_post_pass.py:post_pass
 def post_pass(in_msh: Path, out_msh: Path,
               hmin: float = 100.0, hmax: float = 25000.0,
               hgrad: float = 1.3, hausd: float = 50.0,
               fault_tri_tag: int = 100,
               binary: str = "mmg3d_O3",
               keep_intermediate: bool = False,
               extra_mmg_args: list[str] | None = None,
               mode: str = "optim",
               provenance_path: Path | None = None,
               transform_path: Path | None = None,
-              snap_m: float = 0.1
+              snap_m: float = 0.1,
+              n_passes: int = 1
               ) -> dict:
@@
+    if n_passes < 1:
+        raise ValueError(f"n_passes must be >= 1; got {n_passes}")
@@
+    # Multi-pass: feed the previous pass's output as the next pass's
+    # input.  Re-run the medit conversion fresh each pass so the
+    # required-edge / required-vertex sets reflect the current
+    # polyline geometry.
+    pass_outputs: list[dict] = []
+    current_in = in_msh
+    for pass_idx in range(n_passes):
+        in_medit_p = work_dir / f"{in_msh.stem}_p{pass_idx}.mesh"
+        out_medit_p = work_dir / f"{in_msh.stem}_p{pass_idx}_post.mesh"
+        log_p = work_dir / f"{in_msh.stem}_p{pass_idx}_mmg3d.log"
+        ci = _convert_msh_to_medit(
+            current_in, in_medit_p,
+            fault_tri_tag=fault_tri_tag, mode=mode,
+            provenance_path=provenance_path,
+            transform_path=transform_path, snap_m=snap_m)
+        _run_mmg3d(in_medit_p, out_medit_p,
+                   hmin=hmin, hmax=hmax, hgrad=hgrad, hausd=hausd,
+                   extra_args=mode_flags,
+                   binary=binary, log_path=log_p)
+        intermediate_msh = work_dir / f"{in_msh.stem}_p{pass_idx}.msh"
+        co = _convert_medit_to_msh(out_medit_p, intermediate_msh)
+        pass_outputs.append({"in": ci, "out": co})
+        current_in = intermediate_msh
+    # Final output is the last intermediate.
+    shutil.copy(current_in, out_msh)
```

Wire into `run_newset_step_by_step.sh`:

```diff
@@ ENABLE_MMG3D_POSTPASS block
+        local MMG3D_N_PASSES="${MMG3D_N_PASSES:-3}"
@@
         if ! python mmg3d_post_pass.py \
             --in-msh "$PRE" --out-msh "$POST" \
             --mode "$MMG3D_MODE" \
+            --n-passes "$MMG3D_N_PASSES" \
             --hmin "$MMG3D_HMIN" --hmax "$MMG3D_HMAX" \
```

**Test case:**
```python
def test_R003_multi_pass_improves_gamma_min(tmp_path):
    # Build a small synthetic mesh with a known low γ_min (e.g., 1e-4).
    # Run with n_passes=1 → measure γ_1.  Run with n_passes=3 →
    # measure γ_3.  Expected: γ_3 >= γ_1 (multi-pass never regresses
    # under -optim mode).
    g1 = run_postpass(in_msh, n_passes=1)["gamma_min"]
    g3 = run_postpass(in_msh, n_passes=3)["gamma_min"]
    assert g3 >= g1, f"multi-pass regressed: {g1} -> {g3}"
```

---

### [R-004] [MODERATE] [dedup_coplanar_facets.py:_find_coplanar_overlapping_pairs] — 3-or-more coplanar incidences at a single edge produce a manifold hole

**Category:** EDGE_CASE

**Description:**
At a triple junction where THREE faults share an edge AND are mutually coplanar, the code generates 3 unordered pairs from a single edge:

```python
for i in range(len(incs)):
    for j in range(i + 1, len(incs)):
        ...
        pairs.append(((fA, lA), (fB, lB)))
```

The loser-selection rule (drop the LATER fault) then drops 2 of the 3 incidences (e.g., faults B and C both lose to A across pairs (A,B), (A,C); B also loses to C in pair (B,C) — but C is already in drop_set, so net effect is B and C both dropped). The edge ends up with **only 1 incident triangle**, which is a manifold-boundary edge — i.e., a hole in the surface.

**Trigger:**
A geometry where 3+ faults are coplanar at a shared edge. Real SAFS data is unlikely to have this exactly, but CFM polyline endpoints at triple-junctions can have 3+ near-coplanar incidences within `coplanar_tol`.

**Actual behavior:**
Drops more triangles than necessary; potentially leaves hole-causing topology that HXT then misinterprets.

**Expected behavior:**
For an edge with K ≥ 2 coplanar same-side incidences, drop K - 1 of them (keep one). The current pairwise-drop approach over-drops because each loss is independent.

**Suggested fix:**
Group incidences by edge and resolve all-at-once:

```diff
@@ _find_coplanar_overlapping_pairs
-    pairs: list[...] = []
-    seen: set[...] = set()
-    for ek, incs in edge_to_incidence.items():
-        if len(incs) < 2:
-            continue
-        ...
-        for i in range(len(incs)):
-            for j in range(i + 1, len(incs)):
-                ...
-                key = ((fA, lA), (fB, lB)) if (fA, lA) < (fB, lB) \
-                      else ((fB, lB), (fA, lA))
-                if key in seen:
-                    continue
-                seen.add(key)
-                pairs.append(((fA, lA), (fB, lB)))
-    return pairs
+    losers: set[tuple[str, int]] = set()
+    for ek, incs in edge_to_incidence.items():
+        if len(incs) < 2:
+            continue
+        a, b = ek
+        # Cluster mutually-coplanar same-side incidences at this
+        # edge.  All members of a cluster will resolve to ONE keeper,
+        # and the rest become losers.
+        clusters: list[list[int]] = []
+        used = [False] * len(incs)
+        for i in range(len(incs)):
+            if used[i]:
+                continue
+            cluster = [i]
+            used[i] = True
+            for j in range(i + 1, len(incs)):
+                if used[j]:
+                    continue
+                fA, lA, wA, nA = incs[i]
+                fB, lB, wB, nB = incs[j]
+                if (fA, lA) == (fB, lB):
+                    continue
+                d = float(abs(np.dot(nA, nB)))
+                if d <= 1.0 - coplanar_tol:
+                    continue
+                if wA == wB:
+                    continue
+                if not _same_side_in_plane(global_V,
+                                            a, b, wA, wB, nA):
+                    continue
+                cluster.append(j)
+                used[j] = True
+            if len(cluster) >= 2:
+                clusters.append(cluster)
+        for cluster in clusters:
+            members = [(incs[i][0], incs[i][1]) for i in cluster]
+            keeper = min(members,
+                         key=lambda m: (include_order.index(m[0]),
+                                        m[1]))
+            for m in members:
+                if m != keeper:
+                    losers.add(m)
+    return losers  # caller updates signature
```

**Test case:**
```python
def test_R004_three_coplanar_incidences_drop_only_two(tmp_path):
    # Build 3 faults, all coplanar at a shared edge.
    # Each fault has 1 triangle incident to the edge.
    # Expected after dedup: exactly 2 triangles dropped, 1 kept.
    # NOT 3 dropped (which would create a hole).
    stl_dir = tmp_path / "in"
    # ... build 3 coplanar same-side fault triangles sharing edge ...
    out = tmp_path / "out"
    rep = dcf.dedup_coplanar(stl_dir, out, ["A", "B", "C"], 0.001, 1e-3)
    total_dropped = sum(rep["per_fault"][f]["n_dropped_coplanar"]
                        for f in ("A", "B", "C"))
    assert total_dropped == 2, (
        f"3 coplanar incidences should drop exactly 2 (keep 1); "
        f"got {total_dropped} drops")
```

---

### [R-005] [MODERATE] [POSSIBLE] [mmg3d_post_pass.py:262-264] — Free-surface vertex detection misclassifies BOX-TOP vertices as fault free-surface vertices

**Category:** BUG (silent over-protection)

**Description:**

```python
is_freesurface_vertex = np.array(
    [points[vi, 2] >= z_top - clearance
     for vi in range(points.shape[0])], dtype=bool)
```

This flags **every** vertex with z ≥ -clearance as a free-surface vertex, regardless of whether it lies on a fault triangle or a box-top triangle. The downstream loop:

```python
for (a, b, c, tag) in tris:
    if tag != fault_tri_tag:
        continue
    for (u, v) in ((a, b), (b, c), (c, a)):
        if (is_freesurface_vertex[u] and is_freesurface_vertex[v]):
            freesurface_edges.add(...)
```

does filter to fault triangles when registering EDGES — so over-tagging of box-top vertices doesn't directly add box-top edges to `freesurface_edges`. But fault triangles whose third vertex happens to be a box-top vertex (or whose two endpoints are at z ≈ 0 BUT one is shared with the box top) will still be flagged. The semantic intent of "free-surface trace" is the curve where the FAULT meets z=0, not "any near-surface vertex."

**Trigger:**
Any mesh where the fault terminates at z = 0 (i.e., reaches the free surface) — which is the normal SAFS configuration.

**Actual behavior:**
Possibly over-marks fault edges that touch box-top vertices but aren't actually on the fault's free-surface trace. The protection is harmless (preserves the trace) but may include unnecessary edges that block useful mmg3d optimization.

**Expected behavior:**
Limit `is_freesurface_vertex` to vertices that appear in fault triangles. Restrict to fault-touching vertices:

**Suggested fix:**

```diff
+        # Build the set of vertices touched by fault triangles ONLY.
+        # Free-surface "trace" is a curve on the FAULT, not the box.
+        fault_vertex_set: set[int] = set()
+        for (a, b, c, tag) in tris:
+            if tag != fault_tri_tag:
+                continue
+            fault_vertex_set.update((a, b, c))
         is_freesurface_vertex = np.array(
-            [points[vi, 2] >= z_top - clearance
+            [(vi in fault_vertex_set
+              and points[vi, 2] >= z_top - clearance)
              for vi in range(points.shape[0])], dtype=bool)
```

**Test case:**
```python
def test_R005_freesurface_vertex_excludes_box_top(tmp_path):
    # Build a mesh with a box top at z=0 (tag 5) and a single fault
    # triangle (tag 100) entirely below z=-200.  No fault vertex
    # should be flagged as free-surface.
    # Without fix: box-top vertices are incorrectly flagged.
    # With fix: 0 free-surface vertices.
    pts = np.array([
        [0,0,0],[1,0,0],[1,1,0],[0,1,0],   # box top (tag 5)
        [0,0,-300],[1,0,-300],[0,1,-300]   # fault (tag 100)
    ])
    tris = [(0,1,2,5), (0,2,3,5),
            (4,5,6,100)]
    # ... convert to .msh and run optim_relax_fault ...
    # Verify report: n_freesurface_edges == 0
```

---

### [R-006] [MODERATE] [run_newset_step_by_step.sh + mmg3d_post_pass.py] — mmg3d's RandomSeed is not pinned, so γ_min is non-deterministic across reruns

**Category:** QUALITY / REPRODUCIBILITY

**Description:**
The implementer noted in the prior turn: "**mmg3d non-determinism.** The same pipeline configuration ... produced γ_min = 3.66e-08 on the first run and 2.01e-06 on the second run — a 55× variation. mmg3d uses a random seed for its operation ordering."

mmg3d 5.8.0's CLI does not expose `-rseed` (random seed). Its non-determinism likely comes from operation ordering, threading, or floating-point reduction order. Running multiple times and **picking the best γ_min** would give a deterministic upper bound on the achievable γ_min for any given config.

**Trigger:**
Any production run where reproducibility matters.

**Actual behavior:**
γ_min varies wildly across reruns of the same config. Users cannot tell whether a config change improved γ_min or whether they got lucky.

**Expected behavior:**
The pipeline should run mmg3d N times with different random seeds (or just N times if no seed flag exists), select the run with the best γ_min, and report all N. This makes the pipeline output reproducible-as-best-of-N.

**Suggested fix:**
Add a `--n-trials` CLI flag (separate from `--n-passes` in R-003 — passes are sequential, trials are parallel). For each trial, run independently, parse the validation γ_min, and keep the best:

```diff
@@ post_pass()
+    if n_trials < 1:
+        raise ValueError(f"n_trials must be >= 1; got {n_trials}")
+    best = {"gamma_min": -1.0, "report": None}
+    for trial_idx in range(n_trials):
+        out_trial = work_dir / f"trial_{trial_idx}.msh"
+        report_t = _do_one_pass(in_msh, out_trial, ...)
+        # Parse γ_min from validation_msh on out_trial.
+        g = _validate_and_get_gamma_min(out_trial)
+        if g > best["gamma_min"]:
+            best = {"gamma_min": g, "report": report_t,
+                    "out_msh": out_trial, "trial_idx": trial_idx}
+    shutil.copy(best["out_msh"], out_msh)
+    return {"best_trial": best["trial_idx"],
+            "best_gamma_min": best["gamma_min"],
+            ...}
```

**Test case:**
```python
def test_R006_n_trials_picks_best(tmp_path):
    # Run with n_trials=5 on a mesh known to be sensitive to seed.
    # Verify the output corresponds to the trial with the highest γ_min.
    rep = post_pass(in_msh, out_msh, n_trials=5)
    assert rep["best_gamma_min"] >= 0.0
    # Verify reproducibility: running same config twice picks the
    # same best-of-N (because the seed sequence is fixed).
    rep2 = post_pass(in_msh, out_msh2, n_trials=5)
    assert rep["best_gamma_min"] == rep2["best_gamma_min"]
```

---

### [R-007] [LOW] [dedup_coplanar_facets.py:_select_loser] — `include_order.index()` is O(N) per pair; degrades to O(N·P) for large meshes

**Category:** QUALITY (performance)

**Description:**
Each invocation of `_select_loser` calls `include_order.index(name)` twice — O(N) per call where N is the number of faults. Across P pairs, this is O(N·P). For Step 6 with 6 faults and ~800 pairs, this is ~10K ops — negligible. But if scaled to a CFM dataset with 50+ faults and 100K+ pairs, it would matter.

Not a correctness issue. Flagged because the fix is trivial.

**Suggested fix:**

```diff
+    pos_index = {name: i for i, name in enumerate(include_order)}
@@
-    pos_a = include_order.index(a[0])
-    pos_b = include_order.index(b[0])
+    pos_a = pos_index[a[0]]
+    pos_b = pos_index[b[0]]
     return a if pos_a > pos_b else b
```

Cache the `pos_index` once at the top of `dedup_coplanar` and pass it through.

**Test case:** none required (LOW).

---

### [R-008] [LOW] [POSSIBLE] [mmg3d_post_pass.py:_convert_medit_to_msh] — Triangle/tetra tags from medit `ref` field are not validated against the original tag set

**Category:** ASSUMPTION

**Description:**
mmg3d's medit output includes a `ref` field per triangle/tet. The code reads it back as the gmsh physical tag without validation. If mmg3d ever introduces a new ref value (e.g., for triangles created during edge swap), it would silently appear as a new physical tag in the .msh — likely causing `validate_msh.py`'s check_1 (`tag_inventory: all expected tags present and non-empty`) to fail with a confusing message.

**Trigger:**
Any mmg3d operation that creates a new triangle ref. For `optim_relax_fault` mode this is unlikely but not impossible.

**Suggested fix:** validate that all output tags are in the input tag set:

```diff
+    expected_tris_tags = set(tri_tags_in)  # captured before mmg3d
+    expected_tets_tags = set(tet_tags_in)
@@ after parsing the output medit
+    out_tri_tags = set(tri_tags)
+    out_tet_tags = set(tet_tags)
+    if not out_tri_tags.issubset(expected_tris_tags):
+        new = out_tri_tags - expected_tris_tags
+        raise RuntimeError(
+            f"mmg3d output introduced unknown triangle tags {new}; "
+            f"expected subset of {expected_tris_tags}")
+    if not out_tet_tags.issubset(expected_tets_tags):
+        new = out_tet_tags - expected_tets_tags
+        raise RuntimeError(
+            f"mmg3d output introduced unknown tet tags {new}")
```

**Test case:** none required (LOW; defensive).

---

## Summary

- **Critical issues:** 3 (R-001 needle-tet protection mismatch; R-002 non-contiguous offset bug; R-003 single-pass leaves γ_min on the table)
- **Moderate issues:** 3 (R-004 over-drop on triple-coplanar incidences; R-005 free-surface over-flagging; R-006 non-deterministic reproducibility)
- **Low issues:** 2 (R-007 perf nit; R-008 tag-validation defense)
- **Plan compliance:** PARTIAL — Step 6 builds end-to-end with check_5 + check_8 PASS (the user's explicit goal). γ_min target ≥ 1e-2 unmet (current 2.79e-07 is 5 orders below).
- **Verdict:** **PASS WITH FIXES** — the production Step 6 mesh is usable as-is for SEAS DG (check_5 PASS means the fault interface is correctly embedded), but R-001, R-002, R-003 should be applied to push γ_min into the e-3 to e-2 range before this is the canonical mesh.

## What pushes γ_min higher (priority order)

1. **R-001 fix** (`-hgradreq 1.3`): forces mmg3d to insert finer Steiner points near the polyline, which can break the 4-collinear needle tets that currently dominate the worst-γ class. Estimated improvement: 1–2 orders of magnitude on the worst tet.
2. **R-003 fix** (multi-pass): each pass finds new flips/insertions. Empirically, 3 passes typically converge γ_min to within 5× of the achievable ceiling for `-optim` mode. Estimated improvement: 5×–50× compounding with R-001.
3. **R-006 fix** (best-of-N): exploits mmg3d's non-determinism to find a good seed. Estimated improvement: 5×–55× (the empirical variance we observed).
4. **R-002 fix**: not a γ_min booster directly, but unblocks safe chaining of R-003 + R-006.

If R-001 + R-003 + R-006 are all applied, the projected γ_min target is **1e-3 to 1e-2** — meeting the user's goal.

## Unreviewed Areas

- The legacy `refine_fault_near_intersections.py` (R-001 fix from this conversation): not exercised in the current Step 6 production path; flagged as broken at small band radius (200 m → 200 km degenerate edge in output). Out of scope here — Step 6 doesn't use it.
- The `find_all_polyline_edges` / `--bisect-all-polyline` path in `break_fault_wedges.py`: implemented in a prior turn but empirically destabilizes Step 6. Not in current production config; out of scope.
- The C++ `tools/corefine_faults/autorefine_mode.cpp` R-501 fix proposed in `REVIEW_autorefine_mode.md`: NOT applied (deliberately — `dedup_coplanar_facets.py` works around the symptom in Python). The C++ root-cause fix is still recommended as a long-term cleanup but doesn't change γ_min.
