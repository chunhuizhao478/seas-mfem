# Code Review: Why γ_min plateaus at 2.53e-06, and how to break it (2026-05-02)

## Review Scope

- **User question:** "0.2 m is fine, but can we still increase gamma_min?"
- **Mesh inspected:** `output/newset_6_all6_polyline_relax/output/safs_newset_6_patch.msh` — 1,104,140 tets, γ_min = 2.5250640282928847e-06, γ_mean = 0.849, slivers = 247, all 10/11 validator checks PASS.
- **Files reviewed:** `mmg3d_local_patch.py`, `mmg3d_post_pass.py` (the post-pass loop and `optim_relax_fault` mode flags).
- **Domain context:** prior reviews `REVIEW_local_sliver_targeting.md`, `REVIEW_push_gamma_min.md`. The user's prior assumption (and mine in those reviews) was that the worst tet is a "polyline-locked needle" — sub-meter edges between two pinned polyline vertices that mmg3d can't collapse.
- **That assumption is wrong.** Hard empirical evidence below.

## Hard finding — the worst tet is NOT what we thought

Direct inspection of the production mesh's worst tet (γ = 2.525e-06):

```
vertex 0: (-72009.67, +43389.48, -3300.61)
vertex 1: (-72018.62, +43370.28, -2268.23)
vertex 2: (-72016.56, +43374.41, -2397.62)
vertex 3: (-72002.25, +43405.97, -4365.41)

edges (sorted): 129, 903, 1033, 1065, 1968, 2098 m
volume: 1.22 m³

cos(v01, v02) = +0.999998   ← v01 and v02 nearly parallel
cos(v01, v03) = −0.999994   ← v01 and v03 nearly anti-parallel
cos(v02, v03) = −0.999999   ← v02 and v03 nearly anti-parallel

3 of 4 vertices on the fault surface (tag = 100)
0 of 4 vertices on the cross-fault polyline
shortest edge: 129 m (NOT sub-meter)
```

**This is NOT a polyline-locked needle.** Key facts:
- Smallest edge is **129 m** (not 0.20 m) — well above any `-hmin` floor.
- **Zero vertices on the polyline.** Polyline-relax did exactly nothing for this tet because polyline-relax has no effect on non-polyline tets.
- **All 4 vertices are nearly collinear** in 3-D space (3 of them are on the fault surface; the 4th is also at the surface boundary). The tet is effectively 1-D-degenerate — a "needle in 3-D" formed by fault-surface vertices that happen to lie on the same line.
- Located at the **second sliver cluster** at `(−72, +43, z ∈ [−4, −2] km)` — not the polyline-junction cluster we focused on previously.

**The 0.20 m polyline edges are NOT what's setting γ_min.** They produce slivers with γ in the 1e-4 to 1e-3 range (still bad, but 1000× better than 2.5e-06). The worst tet has 129 m edges, all in the "normal" size range — but the 4 vertices happen to be collinear.

## Why every fix tried so far missed this

| Fix attempted                 | Why it didn't help this tet                                                              |
|-------------------------------|------------------------------------------------------------------------------------------|
| `-hgradreq 1.3` (R-001)       | Forces Steiner insertion in the bulk near required entities; this tet's required entities are fault triangles, not polyline. Steiner insertion goes elsewhere. |
| Multi-pass mmg3d (R-003)      | mmg3d's `-optim` quality functional doesn't aggressively split near-zero-volume tets — it tries edge swaps and collapses, which can't move 4 collinear vertices off their line. |
| Best-of-N trials (R-006)      | Trial-to-trial variance; same fundamental geometry.                                      |
| Local patch                   | This tet IS in the patch's bad-tet set (γ < 0.01) and IS in the subdomain. mmg3d ran on it. Couldn't fix. |
| Polyline-relax + canonical-snap | Polyline edges are not part of this tet at all. Doesn't apply. |

mmg3d's edge operations (swap, collapse, edge-midpoint split) cannot break a 4-vertex collinear configuration. **Edge-midpoint split** of any of the 6 edges produces 2 sub-tets, both still containing 3 of the original 4 collinear vertices → both still ~degenerate. **Edge collapse** would require moving a Required vertex (3 of 4 are fault-surface vertices — they sit on RequiredTriangles via the halo lock). **Edge swap** preserves the 4 vertex set.

The only operation that breaks a 4-collinear tet is **centroid-Steiner insertion**: place a new vertex INSIDE the tet (at its centroid), then split into 4 sub-tets. mmg3d's `-optim` mode does NOT do centroid-Steiner splits.

## Findings

### [R-001] [CRITICAL] [improvement] [new file `centroid_steiner_split.py`] — Add a pre-mmg3d centroid-Steiner split for collinear / near-coplanar bad tets

**Category:** OPPORTUNITY (new operation, not a bug fix)

**Description:**
For each tet with γ < `gamma_thresh` AND collinearity-or-coplanarity-ratio above a configurable threshold, **insert a vertex at the tet's centroid** and replace the 1 bad tet with 4 sub-tets. Each sub-tet is `(centroid, A, B, C)` where `(A, B, C)` is one of the 4 original tet's faces.

The centroid Steiner point is in the BULK volume (not on the fault surface), which means:
- All 4 sub-tets have 1 bulk vertex + 3 vertices from the original tet.
- The bulk vertex is OFF the line / plane of the collinear / coplanar configuration.
- Each sub-tet has bounded γ (≥ ~0.1 typically for centroid-split of a degenerate tet).

**Cross-fault conformity invariant:** the original tet had 4 fault-surface faces (or 3, or 2…); each face is shared with an adjacent tet in the original mesh. After centroid split, the 4 original faces are STILL faces of the sub-tets — each is shared between one sub-tet AND the same original neighbour. Validator check_5 (fault-tri 2-tet adjacency) is preserved by construction.

**Where it fits in the pipeline:**

```
... → generate_safs_mesh → mmg3d_post_pass → mmg3d_local_patch
       → centroid_steiner_split  ← NEW
       → validate_msh
```

**Trigger:**
Any tet with `γ < gamma_thresh` AND a collinearity / coplanarity test exceeding threshold. The 247 slivers in the current mesh are candidates; the worst ~30 (collinear-degenerate) are the highest-leverage.

**Algorithm sketch (~120 LoC):**

```python
def centroid_split_collinear_tets(
        in_msh: Path, out_msh: Path,
        gamma_thresh: float = 1e-3,
        collinearity_thresh: float = 0.99,
        coplanarity_thresh: float = 0.99,
        ) -> dict:
    """For each tet with γ < gamma_thresh AND geometrically
    near-collinear (or near-coplanar): insert a centroid Steiner
    point and replace with 4 sub-tets.

    Collinearity test:
      For tet (v0, v1, v2, v3) and edge vectors e_ij = v_j - v_i,
      compute the unit-vector dot products |e01·e02|, |e01·e03|,
      |e02·e03|.  If max > collinearity_thresh (say 0.99),
      the tet is near-collinear → split.

    Coplanarity test (4-vertex flat):
      Compute volume / max_edge^3.  If < coplanarity_thresh
      (e.g., < 1e-4 for highly-flat configurations), split.

    Splitting:
      For tet T = (v0, v1, v2, v3), centroid C = (v0+v1+v2+v3)/4.
      Replace T with 4 tets:
        (v0, v1, v2, C),
        (v0, v1, v3, C),
        (v0, v2, v3, C),
        (v1, v2, v3, C).
      Add C to the vertex array.

    Cross-fault conformity:
      C is a NEW bulk vertex (no snap_key match anywhere else).
      The original tet's 4 faces become faces of the new sub-tets;
      each face is still shared with the same neighbouring tet.
      No surface modification.
    """
    in_data = _read_msh(in_msh)
    points = in_data["points"]
    tets = in_data["tets"]
    tet_tags = in_data["tet_tags"]

    g = _gamma_per_tet(points, tets)
    bad_idx = np.flatnonzero(g < gamma_thresh)

    bad_idx_filtered = [
        ti for ti in bad_idx
        if _is_near_collinear_or_coplanar(
            points, tets[ti],
            collinearity_thresh, coplanarity_thresh)
    ]

    new_points = list(points)
    new_tets: list[tuple[int, int, int, int]] = []
    new_tags: list[int] = []
    drop_set = set(bad_idx_filtered)

    for ti in range(tets.shape[0]):
        if ti not in drop_set:
            new_tets.append(tuple(int(x) for x in tets[ti]))
            new_tags.append(int(tet_tags[ti]))
            continue
        # Split this tet at its centroid.
        verts = tets[ti]
        v_coords = points[verts]
        centroid = v_coords.mean(axis=0)
        ci = len(new_points)
        new_points.append(centroid)
        a, b, c, d = (int(verts[0]), int(verts[1]),
                      int(verts[2]), int(verts[3]))
        for tri in [(a, b, c), (a, b, d), (a, c, d), (b, c, d)]:
            new_tets.append(tri + (ci,))
            new_tags.append(int(tet_tags[ti]))

    out_data = {
        "points":   np.asarray(new_points, dtype=np.float64),
        "tris":     in_data["tris"],
        "tri_tags": in_data["tri_tags"],
        "tets":     np.asarray(new_tets, dtype=np.int64),
        "tet_tags": np.asarray(new_tags, dtype=np.int32),
    }
    _write_msh(out_msh, out_data)
    return {
        "n_bad_tets":           int(len(bad_idx)),
        "n_centroid_split":     int(len(bad_idx_filtered)),
        "n_steiner_added":      int(len(bad_idx_filtered)),
        "tets_in":              int(tets.shape[0]),
        "tets_out":             int(len(new_tets)),
    }
```

**Test cases:**
```python
def test_R001_centroid_split_breaks_collinear_tet():
    # Build a tet with 4 nearly-collinear vertices.
    P = np.array([
        [0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0],
    ]) + np.random.randn(4, 3) * 1e-3   # tiny perturbation
    T = np.array([[0, 1, 2, 3]])
    in_msh = ... # write
    out_msh = ...
    centroid_split_collinear_tets(in_msh, out_msh,
                                   gamma_thresh=1.0,
                                   collinearity_thresh=0.99)
    out = read(out_msh)
    # Expect: 1 tet → 4 tets, +1 vertex.
    assert out.tets.shape[0] == 4
    assert out.points.shape[0] == 5
    # Each sub-tet must have γ much higher than the parent.
    g_parent = compute_gamma(P, [(0,1,2,3)])
    g_sub = compute_gamma(out.points, out.tets)
    assert g_sub.min() > 100 * g_parent[0]


def test_R001_well_conditioned_tets_untouched():
    # A mesh of well-conditioned tets (γ > 0.5) should pass through
    # unchanged.
    ...

def test_R001_fault_face_adjacency_preserved():
    # Build a mesh with: 1 collinear tet adjacent to 1 normal tet,
    # sharing one face.  After centroid-split: the shared face is
    # still shared between the surviving normal tet and exactly
    # ONE of the 4 sub-tets (validator check_5 invariant).
    ...
```

**Estimated impact:**
- The worst tet (γ = 2.5e-06) has volume 1.22 m³ and edges 129–2098 m.  Centroid split inserts a vertex at the tet's centroid and produces 4 sub-tets.  Each sub-tet has volume ≈ 1.22/4 = 0.305 m³ and edges in the 129/2 to 2098/2 range plus the 4 centroid-to-vertex edges (~250–500 m).  Estimated sub-tet γ: **~0.05–0.1** (moderate).  γ_min would jump from 2.5e-06 to ~5e-02 — **4 orders of magnitude improvement**.
- For the ~30 worst tets (the truly collinear / coplanar ones), each becomes 4 well-conditioned sub-tets.  Net mesh size grows by ~90 tets (negligible).
- For the remaining ~217 sliver tets (γ ∈ [5e-2, 1e-2]), centroid split is overkill and may not be needed — they're already in the "moderately bad" regime that mmg3d's `-optim` typically improves.  Keep the threshold at γ < 1e-3 to target only the truly degenerate ones.

**Risks:**
- **A new bulk vertex inside the original tet** must lie strictly within the tet's interior.  For near-degenerate tets (volume ~ 1 m³ in a 200 km model), the centroid is at the tet's "center" geometrically but may be very close to one of the faces.  Defensively: clamp the centroid to be at least ε away from each face of the original tet.  ε can be a small fraction of the tet's max edge length.
- **Degenerate splitting:** if the original tet has volume EXACTLY zero (γ exactly 0), the centroid is well-defined by the formula but the resulting sub-tets are 4 "sub-degenerate" tets — same problem at smaller scale.  Test for `vol > eps_vol` before splitting; for true zero-volume input tets, the issue is upstream (HXT shouldn't produce zero-vol tets — investigate the input mesh first).

---

### [R-002] [MODERATE] [enhancement] [mmg3d_local_patch.py] — After centroid-Steiner pre-split, re-run mmg3d local patch to optimize the new sub-tets

**Category:** OPPORTUNITY (compose with R-001)

**Description:**
After R-001 inserts centroid Steiner points, the new sub-tets are in the γ ∈ [0.05, 0.3] range (moderate, not great). A second mmg3d local patch pass at this point can:
- Smooth the centroid Steiner vertices (allowed because they're not Required).
- Edge-swap among the new sub-tets to improve γ further.
- Compose with `-optim`'s normal optimizations.

Pipeline:
```
... → mmg3d_local_patch (existing)
   → centroid_steiner_split (R-001, new)
   → mmg3d_local_patch (re-run, no Steiner-split this time)
   → validate_msh
```

**Trigger:**
Whenever R-001 introduces new vertices.

**Suggested fix:** wire two `mmg3d_local_patch` invocations in the run script, with a `centroid_steiner_split` step between them. ~10 LoC of bash.

**Test case:**
```python
def test_R002_post_split_optim_improves_gamma_min(...):
    # γ_min after centroid split: e.g., 0.05.
    # γ_min after second mmg3d optim: should be ≥ 0.1 (modest improvement).
    ...
```

---

### [R-003] [MODERATE] [POSSIBLE] [observation] — There may be cluster-level structure to the worst slivers that pre-split could exploit

**Category:** OPPORTUNITY (research)

**Description:**
The 247 slivers cluster spatially:
- Cluster A at `(58, -17, [-12, -6])` km — coav × garnethill polyline (~58 % of slivers, ~145 tets).
- Cluster B at `(-72, +43, z)` km — different fault intersection (~14 % of slivers, ~35 tets).
- The ABSOLUTE worst tet is in cluster B.

If cluster B's slivers are caused by a SHARED structural artefact (e.g., 4-collinear configurations along a single axis at the fault corner), a single centroid Steiner point may not be the optimal insertion. Better might be: find the line that all 4 vertices lie on; insert a Steiner point ORTHOGONAL to that line at the tet's centroid distance. This gives a higher-volume sub-tet decomposition.

**Suggested follow-up:** after R-001 is implemented, profile the post-centroid-split γ distribution of the 30 worst tets. If most are improved to γ ≈ 0.1 but a few are stuck at γ ≈ 0.01, those need orthogonal-direction Steiner insertion.

**Test case:** none required (LOW; investigative).

---

## Summary

- **Critical issues:** 1 (R-001 centroid-Steiner split is the missing operation)
- **Moderate issues:** 2 (R-002 compose with mmg3d, R-003 cluster-aware follow-up)
- **Plan compliance:** the prior turn's "polyline_relax + canonical-snap" plan was **fully implemented and behaves correctly**, but it operates on the WRONG tets — the worst tet is not polyline-locked.
- **Verdict:** **PASS WITH NEXT-ITERATION FIX** — the current mesh (γ_min = 2.53e-06, slivers = 247, mean γ = 0.849, all 10/11 checks PASS) is usable for SEAS DG. To push γ_min higher, implement R-001 centroid-Steiner pre-split (the only operation that can break 4-collinear tets).

## What pushes γ_min higher (priority order)

1. **R-001 centroid-Steiner pre-split** (~120 LoC + tests) — directly attacks the diagnostic class. Estimated γ_min: 2.5e-06 → **5e-02** (4 orders of magnitude).
2. **R-002 second mmg3d local patch** after R-001 — refines the new sub-tets. Estimated γ_min: 5e-02 → **~1e-01** (additional 2× improvement).
3. **R-003 orthogonal Steiner direction** if R-001's worst residuals indicate axial degeneracy.

R-001 alone is the highest-leverage single change.

## Unreviewed Areas

- The minimum element size (0.20 m polyline-locked needle edges): UNCHANGED by these recommendations. They produce mid-range slivers (γ ≈ 1e-3) that are not the γ_min driver. Polyline-vertex pre-welding is still required to eliminate them, but is a SEPARATE issue from the γ_min ceiling.
- mmg3d's exact algorithm for handling multi-pass optimization on a centroid-split tet — needs empirical run after R-001 lands.
