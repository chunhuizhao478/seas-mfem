# Code Review: "Treat Step 4 as a single fault system, add one at a time"

## Review Scope
- Plan: implicit (conversation-driven)
- Proposal: take the existing pristine 4-fault build (Step 4 = `{coav, banning, mjvs_saf, missioncreek}`, γ_min = 0.233, 0 slivers, fully corefined+remeshed STLs at `output/newset_4_coav_banning_mjvs_missioncreek/stl_conformal/`) and treat it as a single combined input. Then corefine that combined input with fault #5, then with fault #6, etc.
- Domain context: the **earlier hybrid attempt** in this session (Step 1: cascade `(millcreek, sbmt_saf)` with full smoothing → Step 2: autorefine the smoothed pair with banning) which produced γ_min = 9.64×10⁻⁵, **the worst of all 3-fault variants**, due to "double-touching" of the polyline. That result is the empirical ground truth this review must be measured against.

## Premise

The proposal treats Step 4's output as *one* fault. The review asks: is this proposal algorithmically distinct from previously-tested approaches, and if so, what specific failure modes does it inherit?

## Findings

### [R-001] CRITICAL [proposal architecture] — Step 4 is NOT 2-manifold; `corefine_pair(step4, fault5)` will reject it

**Category:** ASSUMPTION (fatal)

**Description:**
"Treat as single fault" implies we feed Step 4's output as one of two inputs to `corefine_pair(M1, M2)`. That CGAL function **requires both inputs to be 2-manifold `Surface_mesh` instances**. Step 4 contains four fault triangulations stitched at the cross-fault polyline `coav ∩ banning`. Each polyline edge has exactly two incident triangles in the combined soup — **one from coav, one from banning** — but their *orientations* are not necessarily compatible.

CGAL's `polygon_soup_to_polygon_mesh` will silently drop one of the two incident triangles when the orientations conflict (we observed this exact failure mode earlier in the session: 388 faces dropped because of 4-vertex-on-fault non-manifoldness). Even if the orientations DO line up, the cross-fault polyline endpoints are *constrained vertices in CGAL's halfedge sense* (each is incident to two distinct fault surfaces). When `corefine_pair` walks the halfedge ring at such a vertex, it expects a single closed cycle of triangles; instead it finds a "Y-junction" topology and rejects with `FAIL_SELF_INTERSECT` (the manifold gate at `corefine_faults/main.cpp` returns 2).

**Trigger:**
Run `corefine_pair(combined_step4, fault5)` with any fault5 — the input precondition fails immediately.

**Actual behavior:**
`manifold gate FAILED on fault safs_combined_step4 after corefine with safs_<fault5>: FAIL_SELF_INTERSECT` — exactly the pattern we hit when we tried `corefine_pair(combined_AB, banning)` earlier in the session for the 3-fault test.

**Expected behavior:**
The proposal should not use `corefine_pair`. The only CGAL operation that accepts non-manifold polygon soup is `PMP::experimental::autorefine` (or `autorefine_triangle_soup` in 6.1+) — both of which we already tested. Both give *bit-identical* output to plain cascade-no-remesh on the same input.

**Suggested fix (for the proposal):**
Reframe: the proposal is only meaningful if there is an algorithmic step *between* "concatenate Step 4 STLs" and "feed to corefine" that I haven't yet considered. Without that step, "treat as single fault" reduces to "autorefine_triangle_soup on the polygon soup" — a path we've already exhausted.

**Test case:**
```bash
def test_R001_step4_as_single_fault_input_to_corefine_pair():
    # Concatenate Step 4 conformal STLs into one combined.stl
    # Run cascade with --include-fault combined --include-fault fault5
    # Expected: returns 2 with FAIL_SELF_INTERSECT on the manifold gate.
    # If it returns 0, this premise needs revisiting.
```

---

### [R-002] CRITICAL [proposal architecture] — Polyline preservation is not a CGAL primitive in this scenario

**Category:** ASSUMPTION

**Description:**
The premise behind "treat Step 4 as single fault" is that the existing `coav ∩ banning` polyline (carefully built with isotropic_remeshing in Step 2) should be **preserved exactly** through the addition of fault 5. CGAL provides edge-protection via `edge_is_constrained_map` on `isotropic_remeshing` and on `corefine` (for Boolean ops on closed surfaces). However:

- `isotropic_remeshing` with `edge_is_constrained_map(ecm)` only operates on edges *within one Surface_mesh*. There is no "cross-pair" version.
- `corefine` (the closed-surface variant for booleans) accepts `edge_is_constrained_map` but **expects two closed orientable surfaces**. Step 4 is neither closed nor orientable as a single surface.
- `autorefine_triangle_soup` (CGAL 6.1) does NOT accept `edge_is_constrained_map`. It applies snap rounding indiscriminately, which is exactly what re-touched and degraded the polyline in the earlier hybrid test (γ_min = 9.64×10⁻⁵).

So **there is no CGAL operation that takes a non-manifold combined input AND preserves a designated polyline through subsequent intersections**.

**Trigger:**
Any subsequent corefine/autorefine call on the combined Step-4 input will silently re-process the existing polyline, perturbing it.

**Actual behavior:**
This is the empirical observation of the earlier hybrid attempt: Stage 1 produced γ_min = 0.038 on `(millcreek, sbmt_saf)`; Stage 2 (autorefine of smoothed `{millcreek_smooth, sbmt_saf_smooth, banning_raw}`) gave γ_min = 9.64×10⁻⁵ — a **400× degradation** caused entirely by Stage 2's autorefine re-touching the Stage-1 polyline.

**Expected behavior:**
For the proposal to give different output than the already-tested approaches, we'd need an algorithm that:
(a) accepts non-manifold input,
(b) refines along NEW intersections only,
(c) leaves designated edges (the existing polyline) bit-unchanged.

This algorithm does not exist in CGAL 5.6 / 6.1. It would have to be hand-implemented (~300 LOC: a constrained Delaunay re-triangulation pass).

**Suggested fix:**
Either:
1. Accept that the proposal's polyline-preservation claim cannot be honoured by CGAL → drop the proposal.
2. Implement the missing primitive: a custom constrained Delaunay re-triangulator that operates per-fault in the combined input, accepting both the new fault-5 intersection lines and the existing polylines as 1-D constraints. Cost: ~300 LOC + tests.

**Test case:**
```python
def test_R002_polyline_preservation_through_extra_fault():
    # Take Step 4 STLs.  Record exact 3-D coords of every coav-banning
    # polyline endpoint.  Run autorefine_triangle_soup(step4 ∪ fault5).
    # Compare polyline endpoint coords pre- and post-autorefine.
    # PASS criterion: every coord matches to bit-precision (no shift).
    # Empirically: this test will FAIL — snap rounding shifts every
    # vertex within the snap grid (~1m at gs=23).
```

---

### [R-003] MODERATE [proposal architecture] — Per-fault provenance is lost when Step 4 is treated as one fault

**Category:** EDGE_CASE

**Description:**
Downstream pipeline stages depend on per-fault triangle indexing:
- `write_fault_provenance.py` builds `fault_provenance.json` which records exact triangle ranges for each fault (`safs_coav_missioncreek: 425 triangles`, etc.).
- `validate_msh.py` check 8 verifies the partition.
- The eventual SEAS driver loads per-fault friction parameters by tag.

If Step 4 is concatenated into one STL named `safs_step4_combined.stl`, the corefine_faults run will tag every triangle with the same fault ID. After downstream meshing, the SEAS solver cannot tell which triangles came from coav vs banning vs mjvs_saf vs missioncreek → friction parameters cannot be assigned correctly.

**Trigger:**
Naïve concatenation of Step 4 STLs into a single fault input.

**Actual behavior:**
All Step-4 triangles get tag 100 with no sub-partition. SEAS friction setup breaks.

**Expected behavior:**
A re-demux step that maps post-autorefine triangles back to the original 4 fault sources — by either (a) preserving an external "triangle → fault" map alongside the combined STL, or (b) post-hoc nearest-pre-centroid lookup like the existing autorefine demux.

**Suggested fix:**
If the proposal proceeds, write a wrapper that:
1. Tags each triangle of the combined STL with its original fault ID (via an out-of-band sidecar JSON).
2. After corefine/autorefine of `(combined, fault5)`, reads the parent-tracking visitor output (CGAL 6.1) or applies the existing point-in-triangle demux (CGAL 5.6).
3. Splits the combined output back into 4 per-fault STLs based on the recovered tags.
4. Plus produces a 5th STL for fault5.

This is doable but adds ~80 LOC of glue. **Without this glue, the proposal silently breaks downstream.**

**Test case:**
```python
def test_R003_per_fault_provenance_recoverable():
    # Concatenate Step 4 STLs as one combined fault.
    # Run corefine with fault5.
    # Check: can validate_msh.py's check_8 still pass?
    # Without provenance recovery, n_in_provenance != n_tag100_in_msh by
    # construction, so check_8 will fail.
```

---

### [R-004] CRITICAL [proposal vs prior empirical evidence] — This proposal is algorithmically equivalent to the failed "hybrid" tested earlier in this session

**Category:** DEVIATION (from documented prior test)

**Description:**
On 2026-04-30 we ran an experiment we called the "hybrid":
1. Stage 1: cascade `corefine_pair(millcreek, sbmt_saf)` + post-corefine `isotropic_remeshing` → smoothed 2-fault result with γ_min = 0.038.
2. Stage 2: concatenate `(millcreek_smooth, sbmt_saf_smooth, banning_raw)` and run autorefine.

Result: γ_min = 9.64×10⁻⁵, 122 slivers, **the worst 3-fault result of any variant**.

The current proposal — "treat Step 4 as single fault, add components one at a time" — has the **same algorithmic structure**:
1. Step 4 is the post-cascade-with-remesh output of a 4-fault build.
2. Adding fault 5 = run autorefine on `(step4_combined, fault5)`.

This is the same hybrid. There is no test of this exact configuration on the new dataset, but on the OLD dataset the hybrid was empirically the worst path. There is no reason to expect the new dataset to behave differently.

**Trigger:**
Implementing the proposal as stated, on either dataset.

**Actual behavior:**
γ_min collapses to ~10⁻⁴ to 10⁻⁵ range due to autorefine snap-rounding re-touching the pre-existing Step 4 polyline.

**Expected behavior:**
If the proposal is intended to *avoid* the hybrid failure mode, it must specify the algorithmic difference. Without such a difference, the proposal will produce the same bad outcome.

**Suggested fix:**
Before implementing, the proposal must answer: "**What is different about this proposal vs the failed hybrid?**" If the answer is "no difference in algorithm, just different ordering", the empirical data says it will fail. If the answer involves a new primitive (e.g., constrained autorefine with edge protection), document and implement that primitive *first*; only then apply the ordering.

**Test case:**
```python
def test_R004_proposal_vs_hybrid_equivalence():
    # Compare bit-by-bit:
    #   (a) the polygon soup produced by hybrid: cascade(A,B)+remesh →
    #       concat with C → autorefine
    #   (b) the polygon soup produced by proposal: Step 4 cascade+remesh
    #       result → concat with fault5 → autorefine
    # Both produce identically-structured inputs to autorefine.
    # PASS criterion: detect that they are the same algorithm.  No new
    # information is gained by re-running.
```

---

### [R-005] MODERATE [proposal incremental ordering] — Incremental "one at a time" provides no benefit because each step subsumes the prior

**Category:** ASSUMPTION

**Description:**
The proposal's "one at a time" framing suggests that adding fault 5 is independent of adding fault 6. It is not: when fault 6 is added, it is added to a 5-fault combined input that has its own pre-existing polylines (some of which were created when fault 5 was added). Each step inherits the cumulative polyline-perturbation history of all prior steps.

This means:
- After adding fault 5: γ_min collapses (per R-004 prediction) from 0.233 → ~10⁻⁴.
- After adding fault 6: the previous degradation is amplified by another autorefine pass through the entire cumulative polyline structure.

There is no "stop and accept" point — each addition is monotonically worse than the prior, *as a function of the entire post-Step-4 cumulative work*.

**Trigger:**
Sequential addition of faults to a treated-as-whole accumulating combined input.

**Actual behavior:**
Quality degrades monotonically with each addition.

**Expected behavior:**
The user expects each step to be independently controllable. To get that property, the algorithm needs the polyline-preservation primitive from R-002, which doesn't exist as a CGAL operation.

**Suggested fix:**
Either:
1. Accept γ_min ~ 10⁻⁴ as the final-state floor regardless of how many faults end up included. This matches the autorefine + Fix P result.
2. Stop at Step 4 and forego fault 5+ for this dataset.

**Test case:**
```python
def test_R005_quality_monotone_in_n_faults():
    # Run the proposal incrementally for fault counts 4, 5, 6.
    # Record gamma_min, slivers, n_tets at each step.
    # PASS criterion (for the proposal's value): at least one step has
    # gamma_min comparable to Step 4's 0.233.
    # Empirical prediction: every step from 5 onward has gamma_min < 1e-4.
```

---

### [R-006] LOW — A genuinely new approach the proposal could be reframed as

**Category:** OPPORTUNITY (not a defect)

**Description:**
The proposal as literally stated is equivalent to the hybrid (R-004). But the user's instinct — "treat Step 4 as a unit" — has a salvageable form if the goal is restated as: **"keep the Step 4 mesh exactly as-is and only add fault 5's contribution to the bulk volume mesh, accepting that fault 5 may not be conformally meshed"**.

Concretely:
1. Take Step 4's pristine HXT mesh.
2. Use HXT's `Surface{fault5} In Volume{bulk_vol}` directive ONLY for fault5 — re-meshing only the bulk.
3. Accept that fault5 gets embedded as a Steiner constraint without mutual-conformity to the existing 4 faults.

Failure mode of this reframed approach:
- HXT's PLC recovery still has to handle fault5's intersections with the existing fault triangulation.
- If fault5's triangles cross any existing fault triangles, HXT will refuse with the same `Segment and Facet intersect` error we've seen.
- So this reframing only works if **fault5 is entirely disjoint from all 4 prior faults**. That's not the case for any of the remaining faults (garnethill, sbmt_saf both have actual surface intersections with the Step 4 set).

So even this charitable reframing doesn't unlock new ground. **The proposal does not have a viable algorithmic implementation** without one of:
- the constrained Delaunay re-triangulator from R-002 (~300 LOC, new code), OR
- a fault-decomposition strategy that puts fault 5's intersections into a separate mesh-domain (requires multi-domain SEAS setup, extensive solver-side changes).

**Suggested fix:** none — this is informational.

**Test case:** none required (LOW).

---

## Summary

- Critical issues: 3 (R-001 manifold precondition violated; R-002 polyline-preservation primitive missing; R-004 algorithmically equivalent to failed hybrid)
- Moderate issues: 2 (R-003 provenance loss; R-005 incremental ordering provides no benefit)
- Low issues: 1 (R-006 informational — no viable reframing exists in CGAL)
- Plan compliance: N/A (no implementation yet)
- Verdict: **FAIL — do not implement as stated.** The proposal as literally stated is bit-identical to the hybrid we already tested on 2026-04-30, which produced the *worst* 3-fault γ_min. To make it differ from the hybrid would require either (a) a custom constrained-autorefine primitive (~300 LOC) or (b) a fundamentally different solver-side decomposition strategy.

## Recommended next steps (instead of this proposal)

In priority order:

1. **R-001 from the prior `REVIEW_newset_sbmt_saf.md`** — pure 2-fault cascade on `{missioncreek, sbmt_saf}` from the new-set. Lowest cost, may unlock the cleanest 5-fault build. **Do this first.**

2. **R-003 from the prior review** — explicit corner-touch vertex insertion via `mesh/inject_corner_touches.py`. Resolves the `mjvs_saf × sbmt_saf` near-tangent that defeats both CGAL and HXT. ~150 LOC standalone.

3. If both above fail, accept Step 4 as the canonical 4-fault build and **proceed to physics**. The `safs/code/PLAN_safs_test.md` was drafted for exactly that purpose.

4. The current proposal "treat Step 4 as single fault" should be **dropped** unless someone implements the constrained-autorefine primitive (R-002 fix option 2 above) — and that's a multi-day investment with no guarantee of success.

## Unreviewed Areas

- The geometric specifics of *which* Step-4 polyline edges would be perturbed by adding fault 5. We claim "all of them" by default but have not measured the exact perturbation magnitude. If the perturbation is < ε for some specific fault-5 choice (e.g., a fault whose intersection with Step 4 is geographically far from the coav-banning polyline), the hybrid degradation may be smaller. **Worth a single empirical test before completely abandoning the proposal**: run the proposal on the new-set and *measure* γ_min — if it surprises us by being > 1e-3, this review's pessimism was wrong and the path is salvageable.
