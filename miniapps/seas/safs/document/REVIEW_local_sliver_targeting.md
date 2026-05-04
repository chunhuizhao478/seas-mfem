# Code Review: Local sliver targeting on Step 6 production mesh (2026-05-02)

## Review Scope

- **User question:** "is γ_min = 2.23e-06 only on a few tets or many tets, maybe we can do something just locally"
- **Mesh inspected:** `mesh/output/newset_6_all6_dedup_mmg3d_p3t3/output/safs_newset_6_mmg3d.msh` (1,099,142 tets, 18,209 fault triangles, check_5 + check_8 PASS).
- **Files reviewed:** `mmg3d_post_pass.py` (post-pass orchestration), `dedup_coplanar_facets.py` (NOT a sliver source), `run_newset_step_by_step.sh` (pipeline wiring).
- **Domain context:** `safs/REVIEW_push_gamma_min.md` (prior turn, R-001..R-008).

## Empirical answer to the user's question

**It's a few tets, very localized.** The full γ distribution:

| Threshold | N tets | Fraction |
|-----------|--------|----------|
| γ < 1e-6  | 0      | 0.0000 % |
| γ < 1e-5  | 7      | 0.0006 % |
| γ < 1e-4  | 16     | 0.0015 % |
| γ < 1e-3  | 54     | 0.0049 % |
| γ < 1e-2  | **94** | **0.0086 %** |
| γ < 5e-2  | 296    | 0.0269 % |
| **γ_mean** | **0.849** | (median 0.863) |

The mesh is overall **excellent** — mean γ = 0.849, median 0.863. 99.99 % of tets have γ ≥ 1e-3. Only **94 tets** (out of 1.1M) have γ < 0.01, and only **7** are below 1e-5.

**Spatial clustering of the 94 worst tets (2 km buckets):**

| Centroid bucket (km) | # slivers | Cumulative % |
|----------------------|-----------|--------------|
| (58, -18, -10)       | 20        | 21.3 %       |
| (58, -16, -12)       | 13        | 35.1 %       |
| (58, -16, -10)       | 9         | 44.7 %       |
| (58, -18, -6)        | 7         | 52.1 %       |
| (58, -18, -8)        | 6         | 58.5 %       |
| ... (17 more buckets, ≤ 5 each) | 39 | 100 % |

**58.5 % of the worst slivers are in 5 spatial buckets, all clustered around `x ≈ 58 km, y ≈ -17 km, z ∈ [-12, -6] km`** — a single ~6 km × 4 km × 6 km region. This is the cross-fault polyline at coav_missioncreek × sbmt_garnethill (depth 6-12 km), the same triple-junction region we identified in the prior turn as the source of the worst tet.

A second smaller cluster (~14 % of slivers) is at `(-72, +43, z)` — a different polyline intersection on the opposite side of the model.

The remaining 27 % are spread across 17 buckets, ≤ 5 slivers each.

**Conclusion: this is fundamentally a local problem**, not a mesh-wide quality issue. A targeted local-remesh on ~5 spatial regions should fix > 50 % of the slivers; tightening on ~22 regions catches all of them.

## Findings

### [R-001] [CRITICAL] [improvement] [mmg3d_post_pass.py + new file] — Add a "local sliver patch" mode that operates ONLY on the spatially-clustered worst tets

**Category:** ASSUMPTION — global multi-pass mmg3d is not the right tool for a problem that is 0.01 % of the mesh, geographically clustered.

**Description:**
The current pipeline runs `mmg3d_O3 -optim -opnbdy -hgradreq 1.3` on **the entire mesh** (1.1M tets), `n_passes × n_trials = 9` times in the recommended 3×3 config. mmg3d optimizes globally — its quality functional weights every tet equally, and per-pass it spends most of its work touching well-conditioned tets that don't need help. The 296 slivers are 0.027 % of the mesh; mmg3d is doing 99.97 % wasted work per pass.

A local-only mode would:
1. Identify the spatial region(s) containing the worst slivers (the diagnostic above already does this).
2. Extract a **subdomain mesh** containing those tets + a "halo" of tets within a buffer distance.
3. Run a **dedicated mmg3d invocation on the subdomain** with much tighter tolerances (`-hausd` 5–10 m instead of 50, `-hgradreq` 1.1 instead of 1.3, `-hsiz` 200 m instead of res_f).
4. Stitch the subdomain output back into the full mesh, gluing along the halo boundary.

The same 9-mmg3d-invocation cost would then be focused entirely on the bad region; the per-tet effort there is ~10 000× higher.

**Trigger:**
Any production mesh where slivers are spatially clustered (which is the SAFS case — confirmed empirically above).

**Actual behavior:**
mmg3d's global pass touches every tet, but its quality functional is dominated by the 99.97 % well-conditioned mass. The few bad tets see a tiny share of the optimization effort.

**Expected behavior:**
A `mode="optim_local_patch"` flag that:
- accepts a list of bounding boxes (or a γ threshold) defining the sliver region,
- extracts the corresponding subdomain via `meshio` / numpy,
- runs mmg3d on the subdomain with tighter parameters,
- merges the result back.

**Suggested fix:** new file `safs/mesh/mmg3d_local_patch.py` (~250 LoC). Sketch:

```python
def local_patch(in_msh: Path, out_msh: Path,
                gamma_thresh: float = 0.01,
                halo_radius_m: float = 2000.0,
                hmin: float = 50.0, hmax: float = 1000.0,
                hgrad: float = 1.1, hausd: float = 10.0,
                fault_tri_tag: int = 100,
                provenance_path: Path | None = None,
                transform_path: Path | None = None,
                ...) -> dict:
    """Find tets with γ < gamma_thresh; extract a subdomain
    containing those tets plus all tets within `halo_radius_m`
    of any sliver centroid; run mmg3d in the subdomain with
    tight parameters; merge back."""
    # 1. Identify bad tet indices.
    bad_tet_ids = _find_bad_tets(in_msh, gamma_thresh)
    # 2. Build subdomain by tet-by-tet flood-fill from each bad
    #    tet outward within halo_radius_m.  Result: subdomain_tet_ids.
    # 3. Write subdomain .msh — preserve tags (incl. tag-100 fault
    #    tris that fall within), mark subdomain BOUNDARY triangles
    #    with a special tag (200) so mmg3d treats them as
    #    RequiredTriangles (locks the halo boundary).
    # 4. Run mmg3d on subdomain.msh with same `optim_relax_fault`
    #    protections + the new tag-200 boundary lock.
    # 5. Stitch: replace the subdomain tets in the input mesh with
    #    the mmg3d output.  Halo boundary vertices match by snap_key.
    ...
```

The halo-boundary lock (tag 200 RequiredTriangles) is what makes the stitch geometrically correct: every boundary vertex of the subdomain stays bit-fixed during mmg3d, so the merge re-attaches without geometric drift.

**Test case:**
```python
def test_R001_local_patch_only_modifies_subdomain(tmp_path):
    # Build a bipyramid-like mesh with one known sliver tet and
    # 100 well-conditioned tets surrounding it.  Run local_patch
    # with halo_radius covering only the sliver region (~5 tets).
    # Expected:
    #   - tets outside the halo have IDENTICAL coords + connectivity
    #     before/after.
    #   - tets inside the halo MAY change (and the sliver should
    #     have improved γ).
    #   - check_5 (fault adjacency) PASS.
    bad_count_before = count_tets_below_gamma(in_msh, 0.01)
    bad_count_after = count_tets_below_gamma(out_msh, 0.01)
    assert bad_count_after < bad_count_before
    # Coordinate-stability check on outside-halo tets.
    assert outside_halo_tets_unchanged(in_msh, out_msh, halo_set)
```

**Estimated impact:** 5–50× improvement in γ_min for the worst tet, with `n_tets` change < 1 % of total mesh size. Compute time per invocation ≈ 1/200 of a global pass.

---

### [R-002] [MODERATE] [improvement] [mmg3d_post_pass.py:_run_one_mmg3d_pass] — `_run_one_mmg3d_pass` does not pass `-hsiz` (uniform target size); the cluster region's bulk size is determined by the pre-existing field, but it could benefit from a tighter local target

**Category:** OPPORTUNITY (smaller-scope alternative to R-001)

**Description:**
mmg3d's `-hsiz <s>` flag forces a uniform target edge size. For the global pass we deliberately don't use it (it would re-mesh the whole bulk to size `s`, defeating the bulk size field). But for a LOCAL invocation around the sliver cluster, `-hsiz 250` (smaller than the local 1000 m bulk) gives mmg3d explicit guidance to refine in that region.

This is achievable WITHOUT R-001's full subdomain extraction: pass an mmg3d `-met` (metric) file with per-vertex isotropic size that's tight (200 m) at sliver-cluster centroids and loose (1000 m) elsewhere. mmg3d honours the metric.

**Trigger:**
Any production mesh where local size adjustment would help concentrated slivers.

**Suggested fix:** add an optional metric-file generator:

```python
def _build_local_size_metric(in_msh: Path, met_path: Path,
                               sliver_centroids: np.ndarray,
                               size_at_sliver_m: float = 200.0,
                               size_far_m: float = 1000.0,
                               falloff_radius_m: float = 2000.0
                               ) -> None:
    """Write a mmg3d isotropic metric .sol file.  At each sliver
    centroid the requested size is `size_at_sliver_m`; far away
    it's `size_far_m`; in between, smooth Gaussian falloff."""
    import meshio
    m = meshio.read(in_msh)
    P = m.points
    sizes = np.full(P.shape[0], size_far_m, dtype=np.float64)
    for c in sliver_centroids:
        d = np.linalg.norm(P - c, axis=1)
        contribution = (size_at_sliver_m
                        + (size_far_m - size_at_sliver_m)
                          * (1.0 - np.exp(-(d / falloff_radius_m) ** 2)))
        sizes = np.minimum(sizes, contribution)
    # Write medit-format `.sol` (one scalar per vertex).
    with met_path.open("w") as f:
        f.write("MeshVersionFormatted 2\n\nDimension 3\n\n")
        f.write(f"SolAtVertices\n{P.shape[0]}\n1 1\n")
        for s in sizes:
            f.write(f"{s:+.17e}\n")
        f.write("End\n")
```

Wire into `_run_one_mmg3d_pass`:

```diff
+    if met_path is not None:
+        cmd_extra = ["-met", str(met_path)]
+    else:
+        cmd_extra = []
     _run_mmg3d(in_medit, out_medit,
                hmin=hmin, hmax=hmax, hgrad=hgrad, hausd=hausd,
-               extra_args=mode_flags,
+               extra_args=mode_flags + cmd_extra,
                binary=binary, log_path=log_path)
```

**Test case:**
```python
def test_R002_local_metric_concentrates_size_at_centroid():
    centroids = np.array([[0.0, 0.0, 0.0]])
    met = tmp_path / "test.sol"
    _build_local_size_metric(input_msh, met, centroids,
                             size_at_sliver_m=100,
                             size_far_m=1000,
                             falloff_radius_m=500)
    # Read back the .sol; verify size at vertex (0,0,0) ≈ 100,
    # at vertex (10000,0,0) ≈ 1000.
    sizes = parse_medit_sol(met)
    assert abs(sizes[0] - 100) < 1.0
    assert abs(sizes[far_vertex_idx] - 1000) < 5.0
```

---

### [R-003] [MODERATE] [diagnostic helper] — Add a γ-distribution diagnostic to `validate_msh.py` so the user gets the histogram automatically

**Category:** QUALITY (diagnostic gap)

**Description:**
The user had to ask for the sliver-distribution analysis manually. The pipeline already runs `validate_msh.py`, but it only reports `gamma_min, gamma_mean, n_gamma_below_0p05` — it doesn't surface the full distribution, the spatial clustering, OR the worst-N coordinates that would drive a local-patch decision.

Adding these to the existing report would make the local-vs-global decision automatic.

**Suggested fix:** extend `validate_msh.py` check_10 with an optional `--gamma-histogram` flag that writes:

```json
{
  "gamma_histogram": {
    "thresholds": [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.2, 0.3],
    "n_below": [0, 7, 16, 54, 94, 296, 909, 3523, 5867]
  },
  "worst_tets_top_20": [
    {"gamma": 2.23e-6, "centroid_m": [58140, -18125, -4174]},
    ...
  ],
  "sliver_clusters_2km": [
    {"bucket_km": [58, -18, -10], "n_slivers": 20},
    ...
  ]
}
```

Out of scope for this review since it's a `validate_msh.py` change, but the agent should know the diagnostic is reusable.

**Test case:** none required (LOW; reporting addition).

---

### [R-004] [LOW] [POSSIBLE] [mmg3d_post_pass.py:_gamma_min_of_msh] — Sequential γ recomputation in `_gamma_min_of_msh` is O(N) per trial; for n_trials=10 on 1M-tet mesh, that's 10s of CPU

**Category:** QUALITY (perf)

**Description:**
`_gamma_min_of_msh` iterates every tet in pure Python, computing γ via numpy ops per-tet. On a 1.1M-tet mesh, this is ~3 s per call. For `n_trials=10`, that's 30 s of pure-Python computation just to pick the best.

Vectorize via numpy broadcasting:

```diff
-    for i in range(tets.shape[0]):
-        p0, p1, p2, p3 = (P[tets[i, 0]], ...)
-        e_lens_sq = (np.sum((p1-p0)**2) + ...)
-        vol = abs(np.dot(p1-p0, np.cross(p2-p0, p3-p0))) / 6.0
-        ...
+    p0 = P[tets[:, 0]]; p1 = P[tets[:, 1]]
+    p2 = P[tets[:, 2]]; p3 = P[tets[:, 3]]
+    e_lens_sq = (np.sum((p1-p0)**2, axis=1)
+                  + np.sum((p2-p0)**2, axis=1)
+                  + ...)  # 6 edges
+    vol = np.abs(np.einsum('ij,ij->i',
+                             p1-p0,
+                             np.cross(p2-p0, p3-p0))) / 6.0
+    g = 12.0 * (3.0 * vol)**(2.0/3.0) / e_lens_sq
+    return float(g[g > 0].min()) if (g > 0).any() else 0.0
```

**Test case:** existing `test_gamma_min_helper_on_regular_tet` passes for both implementations.

---

## Summary

- **Critical issues:** 1 (R-001 local-patch mode is the natural next architecture, not "another global pass")
- **Moderate issues:** 2 (R-002 local-metric file, R-003 diagnostic)
- **Low issues:** 1 (R-004 perf)
- **Plan compliance:** N/A — this review is forward-looking
- **Verdict:** **PASS WITH IMPROVEMENT OPPORTUNITY** — the current mesh is production-quality; pushing γ_min from 2.23e-06 → 1e-3 to 1e-2 requires shifting strategy from global mmg3d passes to local-patch refinement.

## What pushes γ_min higher (priority order)

1. **R-001 local patch (CRITICAL improvement)** — the slivers are 0.027 % of the mesh and 58.5 % cluster in 5 spatial buckets. Targeted local mmg3d with tight tolerances (`-hausd 5`, `-hgradreq 1.1`, `-hsiz 250`) on a small subdomain is the right tool. Estimated γ_min improvement: 5×–50× (the same factor as a global pass, but applied where it counts).

2. **R-002 local-metric file** — smaller cousin of R-001. Don't extract a subdomain; just give mmg3d a per-vertex size hint. Cheaper to implement, smaller improvement (likely 2×–10×).

3. **Stop global passes at 3×3.** Per the diagnostic, mmg3d's `-optim` cannot push γ_min below the polyline-locked structural floor on this mesh. Adding more trials (e.g., 10) would only narrow the variance, not break through the floor.

## Unreviewed Areas

- The actual implementation of R-001 local patch (subdomain extraction + halo locking + stitch-back) is non-trivial — ~250 LoC + tests + boundary-conformity verification. Out of scope here; recommended as a follow-up implementation pass.
- The C++ `tools/corefine_faults` autorefine R-501 fix: still not applied; not a γ_min limit on this dataset.
- Whether the sliver clusters at `(58, -17, [-12, -6])` reflect a CFM source-data quirk (e.g., near-coincident polyline geometry between coav_missioncreek and sbmt_garnethill) that would benefit from upstream fault-data cleaning rather than mesh-side fixes.
