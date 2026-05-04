# Code Review: New-set sbmt_saf addition (2026-04-30)

## Review Scope
- Plan: implicit (conversation-driven incremental build of newset 6-fault dataset)
- Files reviewed:
  - `mesh/run_newset_step_by_step.sh` (the step-by-step builder)
  - `/tmp/pair_test.sh`, `/tmp/pair_autorefine_test.sh` (one-off pair-intersection diagnostics)
  - The conclusion chain "sbmt_saf cannot be added with the current pipeline" derived from 5 specific attempts
- Domain context: `CLAUDE.md` + prior debug-session record from this conversation

## Premise of this review

The implementer concluded that **sbmt_saf cannot be added to the new-set 5-fault build with the current CGAL/HXT pipeline** based on these specific attempts:
1. Cascade default (fails the post-remesh polyline-coincidence gate at main.cpp:799)
2. Cascade + `--target-edge-m 0` (62 slivers, γ_min=7.85e-4 in the prior 3-fault test — but **never re-tested on this exact 5-fault input**)
3. Autorefine + Fix P at target=1000m (γ_min=1.14e-9, HXT works)
4. Autorefine + Fix P at target=2000m (γ_min=1.80e-10, HXT works for 5-fault but garnethill-first only)
5. Frontal-Delaunay (algo3d=4) — bulk γ_min=0.490 but fault-tet conformity broken

The reason this is a *premature conclusion*: several algorithmically distinct paths were not tested. The findings below identify those paths and a couple of code issues that contributed to the premature stopping.

## Findings

### [R-001] CRITICAL [conversation/algorithm choice] — sbmt_saf added jointly with mjvs_saf and missioncreek; never tested standalone

**Category:** ASSUMPTION

**Description:**
Every Step-5-with-sbmt_saf attempt added sbmt_saf as part of a 3-fault component **`{mjvs_saf, sbmt_saf, missioncreek}`**. The corner-touch between mjvs_saf and sbmt_saf at `(-72013, 43382, -2875)` was independently confirmed by HXT but missed by CGAL 6.1. **What was never tried**: process `{missioncreek, sbmt_saf}` as a clean 2-fault cascade pair (the working 2-fault recipe), and process `mjvs_saf` as a separate operation.

The corner-touch between `mjvs_saf` and `sbmt_saf` is geometrically real but is a single point (or a sub-meter polyline). A 2-fault cascade on `{missioncreek, sbmt_saf}` does not need to know about mjvs_saf at all.

**Trigger:**
Any Step-5 attempt that includes sbmt_saf as part of an autorefine group of 3+ faults including mjvs_saf.

**Actual behavior:**
Autorefine processes 3 faults jointly; HXT then either fails outright (target=1000) or produces a usable but γ_min-degenerate mesh.

**Expected behavior:**
The `missioncreek × sbmt_saf` pair was previously shown to be a valid 2-fault cascade target. Process it that way; treat `mjvs_saf` as a disjoint singleton; rely on HXT to handle the (sub-meter or genuinely-zero) corner touch as a coincident-vertex situation rather than a refinement target.

**Suggested fix:**
Run the following experiment before concluding sbmt_saf is unusable:

```bash
# Stage A: cascade on (missioncreek, sbmt_saf) — same recipe as 2-fault build
$BIN --in-stl-dir RAW --out-stl-dir CONFORMAL \
    --include-fault safs_sbmt_missioncreek --include-fault safs_sbmt_saf \
    --target-edge-m 1000 --remesh-iters 3
# Stage B: cascade on (coav, banning) — also a 2-fault working pair
$BIN --in-stl-dir RAW --out-stl-dir CONFORMAL \
    --include-fault safs_coav_missioncreek --include-fault safs_mult_ssaf_banning \
    --target-edge-m 1000 --remesh-iters 3
# Stage C: singleton remesh for mjvs_saf
$BIN --in-stl-dir RAW --out-stl-dir CONFORMAL \
    --include-fault safs_mjvs_saf --target-edge-m 1000 --remesh-iters 3
# Stage D: generate_safs_mesh on all 5 conformal STLs
```

Outcome to record: does the prior `cascade FAIL_SELF_INTERSECT` on `missioncreek × sbmt_saf` (observed at clean6 dataset) reproduce on the *new-set* sbmt_saf STL? The .ts source files have identical MD5 — but `ts_to_stl` post-processing can introduce tiny differences across runs.

**Test case:**
```bash
def test_R001_pair_isolated_cascade():
    # Run cascade on JUST (missioncreek, sbmt_saf) at target=1000m.
    # PASS criterion: no manifold gate failure; both fault outputs
    # have post-corefine triangle counts > inputs (real intersection).
    # If this passes, then the 5-fault build via stages A-D may
    # produce a working mesh that we never tested.
```

---

### [R-002] CRITICAL [conversation/algorithm choice] — Untested: per-fault remesh ONLY around the polyline neighbourhood

**Category:** ASSUMPTION

**Description:**
The "Fix P" pre-autorefine remesh applies isotropic_remeshing to every triangle of every fault, even those far from any intersection. The whole-fault remesh is what changes garnethill from 950 → 680 triangles (a 28% reduction) — and that aggressive smoothing is what shifts triangulation off the original CFM surface and produces the bit-zero edges.

**Never tested**: a *band-limited* remesh that only smooths triangles within K × `target_edge_m` of any other fault's bbox. The far-from-intersection majority of each fault is left at the CFM triangulation. The polyline-adjacent zone is what needs the regularisation; the rest doesn't need to be touched.

**Trigger:**
The whole-fault Fix P call in `autorefine_mode.cpp:fixP_pre_remesh_soup` runs `PMP::isotropic_remeshing(m_pre.faces(), target_edge_m, m_pre, ...)` — passing **all faces** as the patch.

**Actual behavior:**
Every triangle is potentially smoothed, including triangles 50+ km from any other fault. This is wasteful and contributes to bit-zero edges far from intersection geometry.

**Expected behavior:**
Pass only the *bbox-overlap-neighbourhood* faces as the patch. Faces outside that region are auto-treated as a constrained boundary by `isotropic_remeshing` — they stay exactly as-input, and only the band gets regularised.

**Suggested fix:**
In `tools/corefine_faults/autorefine_mode.cpp:fixP_pre_remesh_soup`, accept an optional list of `Vector_3` bbox centres + radii from the OTHER faults, and pass only the in-band faces:

```cpp
void fixP_pre_remesh_soup(
        std::vector<Kernel::Point_3>& pts,
        std::vector<std::vector<std::size_t>>& polys,
        const std::string& short_name,
        double target_edge_m,
        int    n_iters,
        const std::vector<std::pair<Kernel::Point_3, double>>& other_bboxes_inflated = {}) {
    // ... existing setup ...
    std::vector<Mesh::Face_index> faces_in_band;
    if (other_bboxes_inflated.empty()) {
        // Full-mesh fallback (current behaviour).
        for (auto f : m_pre.faces()) faces_in_band.push_back(f);
    } else {
        // Only smooth faces whose centroid is within K * target of any other fault.
        for (auto f : m_pre.faces()) {
            const auto& cent = face_centroid(m_pre, f);
            for (const auto& [c, r] : other_bboxes_inflated) {
                if (CGAL::squared_distance(cent, c) <= r * r) {
                    faces_in_band.push_back(f);
                    break;
                }
            }
        }
    }
    PMP::isotropic_remeshing(faces_in_band, target_edge_m, m_pre, ...);
    // ...
}
```

**Test case:**
```python
def test_R002_band_limited_remesh():
    # Run band-limited Fix P on garnethill against bboxes of (coav, banning,
    # missioncreek, sbmt_saf).  Compare:
    #   (a) full-mesh Fix P: garnethill 950 -> 680 faces
    #   (b) band-limited:    garnethill 950 -> ~860-900 faces (only band shrinks)
    # Then run autorefine + HXT.
    # PASS: gamma_min > 1e-3 (vs prior 1.8e-10 with full-mesh Fix P).
```

---

### [R-003] CRITICAL [conversation/algorithm choice] — Untested: explicit corner-touch handling via vertex insertion before autorefine

**Category:** ASSUMPTION

**Description:**
The `mjvs_saf × sbmt_saf` corner touch at `(-72013, 43382, -2875)` is a single near-tangent point that:
- CGAL 6.1 autorefine reports as DISJOINT (no surface intersection).
- HXT's PLC-recovery detects as a segment-facet intersection and refuses.

**Never tried**: explicitly insert a shared vertex at that corner in BOTH mjvs_saf's and sbmt_saf's STLs *before* corefine. With a shared vertex at the touch point, the corner becomes a topologically-resolved 0-cell rather than a near-tangency that confuses both algorithms.

A bbox-corner pre-detection step would identify near-tangencies (sub-meter approach distance between any vertex of fault A and any triangle of fault B), and add explicit vertices at those points to both faults.

**Trigger:**
Any pair where bbox overlap is YES but CGAL surface-intersection is NO yet HXT detects a PLC-recovery error.

**Actual behavior:**
The pair is treated as DISJOINT (autorefine path skips it) and HXT later trips on it.

**Expected behavior:**
Pre-detect near-tangencies; insert the touch point as a shared vertex in both STLs. Then either autorefine sees it as a degenerate intersection (with snap rounding it becomes a real shared vertex), or HXT can recover the PLC because both surfaces explicitly agree on the vertex.

**Suggested fix:**
Add a Python helper `mesh/inject_corner_touches.py` that:
1. Takes a list of fault STLs.
2. For each pair, computes min vertex-to-triangle distance across all (vertex_A, triangle_B) pairs.
3. If min distance < ε (e.g., 100 m), inserts a vertex at the closest point on triangle B equal to vertex A, and a vertex on A's nearest triangle equal to the same coordinate.
4. Outputs modified STLs ready for corefine_faults.

**Test case:**
```python
def test_R003_corner_touch_insertion():
    # Take mjvs_saf and sbmt_saf, run inject_corner_touches.
    # Verify: at least one new vertex is added at coord (-72013, 43382, -2875).
    # Then run cascade or autorefine on the modified STLs.
    # PASS: HXT generate_safs_mesh succeeds on the resulting mesh.
```

---

### [R-004] MODERATE [run_newset_step_by_step.sh] — Disjoint singleton-remesh runs *every* time, even when the conformal output is already correct

**Category:** BUG (efficiency / determinism)

**Description:**
After a cascade group succeeds, the `processed[]` array tracks which faults are in groups. Disjoint faults are then individually remeshed via `corefine_faults --include-fault X`. **But** if a fault that's a "disjoint singleton" in *this* step was a "disjoint singleton" in a prior step, its conformal STL is NOT idempotent — re-running the singleton remesh loads from `RAW_STL_DIR` again and re-overwrites the conformal STL with a freshly-remeshed copy. Because `isotropic_remeshing` is deterministic but operates on a randomly-seeded refinement schedule via `m_pre.faces()` iteration order, the resulting STL **can differ between runs**.

This means the same fault appears with slightly different triangulations across step-by-step runs of the script, depending on whether you run from `START_FROM=1` or `START_FROM=5`.

**Trigger:**
Re-running with different `START_FROM` values, or restarting after a partial run.

**Actual behavior:**
Per-step output dirs may have different fault conformal STLs for the same fault, depending on script invocation history.

**Expected behavior:**
Either (a) cache the per-fault conformal STL globally and reuse it across steps, or (b) document that per-step output dirs are independent and shouldn't be compared.

**Suggested fix:**
At the top of `run_step`:
```diff
+    # Cache: if a fault's conformal STL exists in a global cache and the
+    # step doesn't put it in any group, copy from cache rather than re-remesh.
+    local CACHE_DIR="output/newset_2000m_cgal/stl_conformal_cache"
+    mkdir -p "$CACHE_DIR"
```

Then in the singleton loop:
```diff
     for f in "${ALL_FAULTS[@]}"; do
         local already=0
         for p in "${processed[@]:-}"; do
             if [ "$p" = "$f" ]; then already=1; fi
         done
         if [ "$already" -eq 1 ]; then continue; fi
+        if [ -f "$CACHE_DIR/$f.stl" ]; then
+            echo "  -- singleton: $f (CACHED — copy from $CACHE_DIR)"
+            cp "$CACHE_DIR/$f.stl" "$OUTDIR/stl_conformal/$f.stl"
+            continue
+        fi
         echo "  -- singleton: $f (remesh to $TARGET_EDGE m)"
         "$COREFINE_BIN_56" \
             --in-stl-dir "$RAW_STL_DIR" \
             --out-stl-dir "$OUTDIR/stl_conformal" \
             --include-fault "$f" \
             --target-edge-m "$TARGET_EDGE" --remesh-iters 3 \
             > "$OUTDIR/singleton_$f.log" 2>&1
+        cp "$OUTDIR/stl_conformal/$f.stl" "$CACHE_DIR/$f.stl"
     done
```

**Test case:**
```bash
def test_R004_singleton_idempotent():
    # Run START_FROM=1; record md5 of mjvs_saf conformal STL after step 3.
    # Run START_FROM=3 (skip 1+2); record md5 of mjvs_saf after step 3.
    # PASS: md5 hashes are identical.
```

---

### [R-005] MODERATE [run_newset_step_by_step.sh:run_step] — `processed[@]:-` empty-array guard hides genuine empty-cascade-group bug

**Category:** EDGE_CASE

**Description:**
The `processed=()` array is referenced as `${processed[@]:-}`. With `set -uo pipefail` removed, that's syntactically fine. But when the cascade-group loop is empty (`pair_groups_csv` is empty for steps 1, 3, 4), `processed` stays as `()` and the singleton loop iterates over all faults. **This is the intended path for steps with no groups.** However, the empty-default `:-` operator silently hides the case where someone passes a malformed `pair_groups_csv` that produces an empty `GROUPS` array (e.g., trailing comma, all whitespace). The loop body never runs, the user thinks all faults are singletons.

**Trigger:**
Pass `pair_groups_csv=","` or `pair_groups_csv=" "` (trailing comma or whitespace).

**Actual behavior:**
Silent fallback to all-singletons mode; no error.

**Expected behavior:**
Reject malformed `pair_groups_csv` with an explicit error.

**Suggested fix:**
```diff
     if [ -n "$pair_groups_csv" ]; then
         local groups_lines
         groups_lines=$(echo "$pair_groups_csv" | tr ',' '\n')
+        # Reject groups with fewer than 2 faults (singleton "groups" make no sense)
+        while IFS= read -r g; do
+            [ -z "$g" ] && continue
+            local n_in_group=$(echo "$g" | tr '/' '\n' | grep -c .)
+            if [ "$n_in_group" -lt 2 ]; then
+                echo "  ERROR: group '$g' has < 2 faults; ill-formed pair_groups_csv"
+                return 4
+            fi
+        done <<< "$groups_lines"
         while IFS= read -r g; do
             [ -z "$g" ] && continue
```

**Test case:**
```bash
def test_R005_malformed_pair_groups_rejected():
    # Call run_step with pair_groups_csv="safs_a/,,safs_b/safs_c"
    # PASS: returns rc=4 with error message about ill-formed group.
```

---

### [R-006] MODERATE [conversation/algorithm choice] — Untested: cascade `--allow-polyline-asymmetry` on `missioncreek × sbmt_saf` for the new-set

**Category:** ASSUMPTION

**Description:**
The conclusion that cascade fails on `missioncreek × sbmt_saf` was based on the `manifold gate FAILED on fault safs_sbmt_missioncreek after corefine with safs_sbmt_saf: FAIL_SELF_INTERSECT` error from `/tmp/pair_test.sh`. The `--allow-polyline-asymmetry` flag exists to demote some gate failures to warnings but **was not tested on this specific pair in the new dataset**.

The `manifold gate FAILED ... FAIL_SELF_INTERSECT` is a different gate from the `polyline-edge-coincidence` gate — but the implementer declared the pair "cascade FAILS" without testing whether the asymmetry flag affects this manifold gate too.

**Trigger:**
Running cascade on `missioncreek × sbmt_saf` from the new-set's raw STLs.

**Actual behavior:**
The pair test script bails out with FAIL on the manifold gate; conclusion drawn that cascade is unusable here.

**Expected behavior:**
Re-test with `--allow-polyline-asymmetry` to see if there's a path that produces *some* output (even if degraded) — that output may be HXT-meshable when combined with downstream fixes.

**Suggested fix:**
Add the flag to the pair re-test:
```diff
 ../tools/build/corefine_faults \
     --include-fault safs_sbmt_missioncreek --include-fault safs_sbmt_saf \
-    --target-edge-m 0 --remesh-iters 0
+    --target-edge-m 0 --remesh-iters 0 --allow-polyline-asymmetry
```

If that still fails with the manifold gate (which it likely will since `--allow-polyline-asymmetry` only relaxes the polyline-coincidence gate, not the manifold gate), at least we have evidence that the `manifold gate` is the real blocker and we can reason about its specific cause (a self-intersection introduced by corefine_pair on these specific surfaces).

**Test case:**
```bash
def test_R006_asymmetry_flag_on_missioncreek_saf_pair():
    # Run cascade with --allow-polyline-asymmetry on missioncreek + sbmt_saf.
    # Record: rc, output STL face counts, gate-failure message text.
    # PASS criterion: any non-empty per-fault STL output (even if HXT later
    # fails) gives us a starting point that we never had.
```

---

### [R-007] LOW [/tmp/pair_test.sh] — Diagnostic script lives in /tmp; not reproducible across sessions

**Category:** QUALITY

**Description:**
The pair-intersection diagnostic script `/tmp/pair_test.sh` was used for a key piece of analysis (the 8-pair survey that drove component-decomposition decisions). It will be deleted on next reboot. The implementer's case-by-case findings are based on this script's output but the script itself is gone.

**Suggested fix:**
Move to `mesh/diag_pair_intersect.sh` (or convert to Python and place in `mesh/`). Same content, just a stable path.

**Test case:** none required (LOW).

---

### [R-008] LOW [conversation/conclusion] — `n_internal_with_2_tets: 0` from FD run not investigated

**Category:** QUALITY

**Description:**
The FD run reported 6670 fault triangles in the .msh BUT 0 of them have 2 incident bulk tets. The implementer attributed this to "FD doesn't enforce internal-face conformity" without verifying with the validator's actual logic.

It's worth confirming whether the validator sees 0 because:
(a) FD really left every fault triangle dangling (in which case the .msh file should have those triangles in a separate disconnected region — possibly visible in the .vtu)
(b) The FD output put the fault triangles on the BULK MESH BOUNDARY (meaning each has 1 incident tet, not 0)
(c) FD did include them but the validator's adjacency lookup is buggy

Quick check: open the `.vtu` and visually count fault triangles that don't have visible bulk tets on both sides. The conclusion "FD is unusable" is correct *if* (a) holds. If (b) or (c), the FD path may yet be partially salvageable.

**Suggested fix:**
Read the .vtu in ParaView; visually verify whether fault triangles are connected to bulk tets on both sides, one side, or neither. Document the actual mode in the README rather than relying on the validator's number.

**Test case:** none required (LOW).

---

## Summary

- Critical issues: 3 (R-001 missioncreek-sbmt_saf 2-fault path untested; R-002 band-limited Fix P untested; R-003 explicit corner-touch insertion untested)
- Moderate issues: 3 (R-004 cache idempotency; R-005 malformed-input guard; R-006 asymmetry flag untested on this pair)
- Low issues: 2 (R-007 /tmp script; R-008 FD failure mode unverified)
- Plan compliance: PARTIAL (the conclusion "sbmt_saf can't be added" was reached without exhausting the algorithmic search space)
- Verdict: **PASS WITH FIXES** — the recent code is functionally correct, but the *algorithmic conclusion* is premature. Three classes of approach (R-001, R-002, R-003) were never tested. Each is a small, contained experiment (≤ 100 LOC change for R-002; ~150 LOC for R-003; just shell-glue for R-001), and any one of them succeeding would change the verdict from "Step 4 only" to a usable 5- or 6-fault build.

## Recommended next experiments (in order of expected payoff and lowest cost)

1. **R-001 first** — pure 2-fault cascade on `missioncreek × sbmt_saf` from the new-set's raw STLs. Cost: 2 binary invocations + 1 generate_safs_mesh. ~5 min. Either produces a working 5-fault build (Stage A+B+C+D recipe), or confirms the pair fails at the manifold-gate level on the new dataset specifically.
2. **R-006 next** — re-run R-001 with `--allow-polyline-asymmetry`. Cost: tiny.
3. **R-003 third** — implement `mesh/inject_corner_touches.py` and force a shared vertex at `(-72013, 43382, -2875)`. Cost: ~150 LOC. May resolve the corner-touch class of failures.
4. **R-002 last** — band-limited Fix P. Cost: ~80 LOC in `autorefine_mode.cpp`. Most invasive but most general.

## Unreviewed Areas

- The choice of `combine-snap-m=0.1` in the disjoint-cascade path. This was tuned for the prior 6-fault clean6 build; whether it's still correct for the new-set's fault interaction patterns hasn't been verified.
- The geometric specifics of garnethill's "triple intersection" failure mode — we know CGAL 5.6 reports it; we never inspected which 3 input triangles trigger it. Knowing the exact triple may inform a targeted CFM-source repair.
- The interaction of `combine-snap-m` with the disjoint banning fault (added raw, no remesh). HXT's STL Merge tolerance is independent of this; potential silent vertex-merging issues across components were not audited.
