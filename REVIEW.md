# Code Review: 2026-04-27 — Mixed-Flux Fault Identification (Attribute vs Coordinates)

## Question
> "I wonder if the current mixed flux formulation identifies the fault from
> attribute or coordinates."

## Direct Answer
**Production code identifies the fault EXCLUSIVELY by boundary-element
attribute.** Cross-rank, fault face IDENTITY is exchanged via global vertex
indices (`HYPRE_BigInt`s from `pmesh.GetGlobalVertexIndices(...)`), NOT by
coordinate proximity. The dispatch path is entirely index-based set lookup.

Coordinate-based fault tagging exists ONLY in the test fixture mesh
builders, where it is used as a one-shot setup step to insert a `BdrTriangle`
with `attr=3`; after that step the attribute is the sole source of truth.

## Review Scope
- Plan: `miniapps/seas/debug_document/tpv104_debug_document/MIXED_FLUX_PLAN.md`
- Files audited (focus: every place fault-identity is consulted):
  - `dynamic/wave_operator.inl` ctor (L286–311 fault list build,
    L179–241 cross-rank Allgatherv merge), `BuildCentralFluxFaceSet_`
    (L1303–1556), 4 dispatch sites (L2547 / L3113 / L4061 / L4573).
  - `drivers/tpv104_driver.cpp` — fault attribute parameter (`bc.fault_attr`).
  - 9 test fixtures: `BuildSmallFaultTetMesh`, `BuildOrthogonalFaultMesh`,
    `BuildSeamFaultMesh`, `BuildLinearFaultMesh`, `BuildTwoTetFaultMesh`
    in the 6 test files plus the dispatch_adjacent_mpi fixture.
- Domain context: `CLAUDE.md`, prior reviews R-1100..R-1605.
- ID range: **R-1700..R-1702**.

---

## Identification Chain — Step-by-step

For traceability — every step of how a fault face becomes "known" to the
mixed-flux dispatch:

| Step | Code site | Mechanism | Coordinate-based? |
|------|-----------|-----------|-------------------|
| (1) Mesh construction (Gmsh production OR test fixture) | mesh file `Physical Surface` tag, OR `mesh.AddBdrTriangle(v0, v1, v2, attr=3)` after a coordinate sieve | Attribute (3) | Production: **no**, set in Gmsh `.msh`. Test fixtures: yes — used **once** to decide which triangle gets `attr=3`. |
| (2) WaveOperator ctor — local fault face list | `wave_operator.inl:286–298`: `for (int b=0; b<mesh_.GetNBE(); ++b) if (mesh_.GetBdrAttribute(b) == bc_.fault_attr) fault_interior_faces_.Append(...)` | **Attribute** (`mesh_.GetBdrAttribute(b) == bc_.fault_attr`) | No |
| (3) Cross-rank merge (R-002) | `wave_operator.inl:179–241`: `pmesh.GetGlobalVertexIndices(gvi)`, `make_global_key(verts)`, `MPI_Allgatherv` of `std::array<HYPRE_BigInt, 4>` | **Global vertex indices** (sorted) | No (GVI is a permutation of mesh-vertex IDs, not coordinates) |
| (4) `BuildCentralFluxFaceSet_::is_fault_face` lambda | `wave_operator.inl:1371–1396` | **Mesh face index ∈ `fault_mesh_face_idx_set`** (built from steps 2+3) | No |
| (5) Adjacent walk (per Zhang 2023 Fig 5) | `wave_operator.inl:1411–1467` | **`E_fault_adj` = local elements that have a face in `fault_interior_faces_ ∪ fault_shared_faces_`** | No |
| (6) Dispatch sites (4 of them) | `wave_operator.inl:2547 / 3113 / 4061 / 4573` | **`central_flux_face_set_.count(f) > 0`** (set lookup) | No |

**Conclusion**: production code is end-to-end attribute/identity-driven.
No coordinate test, distance threshold, or `<` against a numeric position
appears anywhere in the dispatch chain. The only coordinate-based check
in the entire mixed-flux code is in the **test fixture mesh builders**
(see R-1701 below for fixture-only consistency notes).

---

## Findings

### [R-1700] [INFO] [dynamic/wave_operator.inl:1371–1396 + ctor:286–311 + cross-rank:179–241] — No issue: fault identification is purely attribute-based; coordinates are NOT consulted in production

**Category:** N/A (clarification of existing behavior)

**Description:**
The mixed-flux formulation identifies fault faces via the boundary-element
attribute `bc_.fault_attr` (default `3` for TPV104) at three layers:

1. The `WaveOperator` constructor scans `mesh_.GetNBE()` boundary elements
   and appends each one whose `GetBdrAttribute(b) == bc_.fault_attr` to
   `fault_interior_faces_` (interior 2-sided fault) or `fault_shared_faces_`
   (rank-seam fault) (`wave_operator.inl:286–298, 299–311`).
2. For ParMesh, the BE-owner asymmetry on shared fault faces is corrected
   by an `MPI_Allgatherv` of sorted global-vertex-index keys
   (`wave_operator.inl:179–241`); the merged `global_fault_keys` is then
   used to populate `shared_face_bdr_attr_[sf]` symmetrically across ranks.
3. `BuildCentralFluxFaceSet_::is_fault_face(f)` (`wave_operator.inl:1371`)
   reads ONLY from the canonical rank-symmetric `fault_mesh_face_idx_set`
   (built from `fault_interior_faces_` and `pmesh.GetSharedFace(sf)` for
   each `sf ∈ fault_shared_faces_`). The dispatch sites at
   `wave_operator.inl:2547` (RK4 local), `:3113` (RK4 shared), `:4061`
   (ADER local), and `:4573` (ADER shared) consult only
   `central_flux_face_set_.count(f) > 0` — pure set lookup.

No coordinate comparison, distance threshold, or numerical proximity test
appears in this chain. **No bug.**

This is the correct, robust design: it works for any fault geometry the
mesh tags (planar, branching, non-planar — the dispatch is geometry-agnostic
once the attribute is set). The implementation also tolerates
non-orthogonal fault planes, faults of any orientation, and meshes with
multiple disconnected fault patches sharing the same attribute.

**Trigger:** N/A — informational.

**Actual behavior:** As described.

**Expected behavior:** As described.

**Suggested fix:** None.

**Test case:** Existing `test_mixed_flux_face_set` already verifies this
(it asserts `central_flux_face_set_ ∩ fault_interior_faces_ == ∅` and
the AllContinuous variant on the 2-tet fixture, both keyed on attribute
`3`).

---

### [R-1701] [LOW] [tests/unit/test_mixed_flux_dispatch_none.cpp:74, test_mixed_flux_face_set.cpp:155,195, test_mixed_flux_dispatch_adjacent.cpp:123, test_mixed_flux_dispatch_all_continuous.cpp:124, tests/parallel/*:91/103/92] — Inconsistent absolute epsilon for coordinate-based fault tagging across test fixtures (1e-6 vs 1e-8), with no relative-scale guard

**Category:** QUALITY (fixture consistency)

**Description:**
Test fixtures use a coordinate sieve to decide which interior triangle gets
the fault `BdrTriangle` attribute. The epsilon for "is this triangle's
y-centroid on the fault plane" varies across fixtures:

```
test_mixed_flux_dispatch_none.cpp:74            std::abs(cy)             < 1e-6     // y=0 plane
test_mixed_flux_face_set.cpp:155                std::abs(cy - 0.5*L)     < eps      // L=1000, eps=1e-6
test_mixed_flux_face_set.cpp:195                std::abs(cy)             < 1e-6
test_mixed_flux_dispatch_adjacent.cpp:123       std::abs(cy - 0.5*L)     < eps      // eps=1e-6
test_mixed_flux_dispatch_all_continuous.cpp:124 std::abs(cy - 0.5*L)     < eps      // eps=1e-6
parallel/test_mixed_flux_shared_fault_face_excluded_mpi.cpp:91   < 1e-8
parallel/test_mixed_flux_adjacent_mpi.cpp:103                    < 1e-8
parallel/test_mixed_flux_dispatch_adjacent_mpi.cpp:92            < 1e-8
```

The unit-test fixtures use `1e-6` absolute; the parallel-test fixtures use
`1e-8` absolute. All three fixtures use `kL = 1000.0`, so:
- `1e-6 / 1000 = 1e-9` relative.
- `1e-8 / 1000 = 1e-11` relative.

The 3-orders-of-magnitude difference in absolute tolerance is unjustified
by any difference in mesh resolution (all fixtures use `MakeCartesian3D`
or hand-built meshes with vertices at exact coordinates `0`, `L/2`, `L`,
etc., so the centroid `cy` is exactly `L/2` to FP precision — both `1e-6`
and `1e-8` are far from triggering false negatives).

**This is a fixture-only concern** (no production code path uses coordinates
to identify the fault), and on the current 1000m fixtures both tolerances
are correct. But:
- A future fixture with `L = 1.0` (unit cube — common for diagnostic
  tests) would give `eps_relative = 1e-6 / 1 = 1e-6`. Hand-built meshes
  with vertices at `0.5` exactly would still pass, but a mesh perturbed
  by `1e-7` in y would silently be mis-tagged.
- A future fixture with `L = 1e7` (subduction-zone scale) would give
  `eps_relative = 1e-15`, which is below FP precision of vertex
  coordinates — a strict equality check that may FAIL on legitimate
  fault triangles.

The eps should either be (a) made consistent (`1e-8` matches the parallel
tests' precedent), or (b) computed relative to `L` (e.g.,
`std::abs(cy - 0.5*L) / std::max(L, 1.0) < 1e-9`).

**Trigger:** Future fixture authored at non-1000m scale, or a mesh
perturbed by FP noise (e.g., a hand-built mesh that uses computed
coordinates rather than literal `0.5*L`).

**Actual behavior:** Hard-coded absolute eps; values diverge across files.

**Expected behavior:** Uniform relative eps, OR a documented absolute eps
that's robust at the fixture's `L` scale.

**Suggested fix (apply to all 8 sites):**

```diff
-         if (std::abs(cy - 0.5 * L) < 1e-6)   // or 1e-8 in parallel files
+         // R-1701: relative epsilon handles fixtures of any scale.
+         const real_t eps_rel = 1e-9;
+         if (std::abs(cy - 0.5 * L) < eps_rel * std::max<real_t>(L, 1.0))
```

Alternatively, since these are hand-built fixtures with vertices at
exact-integer coordinates, the epsilon could be tightened to `1e-12`
relative — far stricter than any mesh-perturbation noise but still safe
for FP-equality on integer-derived vertex positions.

**Test case:** N/A (fixture consistency cleanup; behavior on the current
fixtures is unchanged).

---

### [R-1702] [POSSIBLE] [LOW] [dynamic/wave_operator.inl:286–311] — `WaveOperator` ctor walks `mesh_.GetNBE()` ONCE and assumes the BE table is FROZEN for the lifetime of the wave operator; if a future driver adds boundary elements post-construction (e.g., dynamic crack propagation marking new fault triangles), the cached `fault_interior_faces_` becomes stale silently — same lifecycle hazard as R-1408 but at the ctor level

**Category:** EDGE_CASE / ASSUMPTION

**Description:**
The fault-face-identification chain caches `fault_interior_faces_` and
`fault_shared_faces_` at the WaveOperator constructor. Any subsequent
modification of the underlying mesh's BE table (which is not currently
exposed by any driver, but is theoretically possible) would not be
reflected in the cached lists, and the mixed-flux dispatch would
silently mis-classify newly-added fault faces.

The R-1408 doc-comment in `wave_operator.hpp:130–145` documents this for
`SetMixedFluxMode`-time caching of `central_flux_face_set_`. The same
invariant applies one layer deeper: `fault_interior_faces_` itself is
ctor-frozen.

**Trigger:** Future SEAS extension that mutates the mesh (e.g., dynamic
fault tracking via `mesh_.AddBdrTriangle(...)` after WaveOperator
construction).

**Actual behavior:** New fault faces silently treated as non-fault →
mixed-flux dispatch incorrectly fires Central on them.

**Expected behavior:** Documented invariant + (optionally) a debug-only
assertion that `mesh_.GetNBE()` hasn't changed since construction.

**Suggested fix:**

```diff
@@ wave_operator.hpp ~L80-83
    /// @brief Construct a WaveOperator on the given mesh.
+   ///
+   /// IMPORTANT (R-1702 lifecycle invariant): the constructor scans
+   /// `mesh.GetNBE()` ONCE to populate `fault_interior_faces_` and
+   /// `fault_shared_faces_`.  If the mesh's BE table changes after this
+   /// call (e.g., dynamic crack propagation, AMR refinement, post-ctor
+   /// `AddBdrTriangle`), the cached lists become stale and the
+   /// downstream `is_fault_face` predicate, FaultBasis frame, and
+   /// mixed-flux dispatch will mis-classify the newly-added faces.
+   /// The wave operator currently assumes the mesh is constant for its
+   /// lifetime — this is true in production TPV104 but should be
+   /// verified before using on a refining/cracking mesh.
    WaveOperator(MeshType &mesh, int order, ...)
```

And, optionally, in debug builds:

```diff
@@ wave_operator.inl Mult / AdvanceADER (top of each)
+#ifndef NDEBUG
+   MFEM_ASSERT(mesh_.GetNBE() == cached_nbe_at_ctor_,
+               "R-1702: mesh BE table mutated since WaveOperator ctor "
+               "(was=" << cached_nbe_at_ctor_ << ", now="
+               << mesh_.GetNBE() << "); fault face lists are stale.  "
+               "If dynamic mesh modification is intended, the wave "
+               "operator must be reconstructed.");
+#endif
```

(`cached_nbe_at_ctor_` would be a private member set in the ctor.)

**Test case:** N/A in production (no driver mutates the mesh post-ctor);
add a defensive test only if dynamic-mesh support is planned.

---

## Summary
- Critical issues: **0**
- Moderate issues: **0**
- Low issues: **2** (R-1701 fixture-eps consistency, R-1702 ctor-time
  freeze documentation)
- Informational: **1** (R-1700 — direct answer to the user's question)
- Plan compliance: **FULL**
- Verdict: **PASS**

The mixed-flux formulation identifies the fault **by boundary-element
attribute** (production) and **by global vertex indices** (cross-rank
merge of the BE-owner asymmetry). No coordinate comparison appears in
the production dispatch chain. The only coordinate-based code is in the
test fixture mesh builders, where it is a one-shot setup step that
imposes the `attr=3` BdrTriangle; after that step, every downstream
consumer reads attribute, not position.

The two LOW findings (R-1701, R-1702) are quality / lifecycle
documentation hardening; neither indicates an actual bug in the
current implementation.

## Unreviewed Areas
- The Gmsh production mesh `tpv104/mesh/tpv104_repro.msh` and its variants
  — verified via `bc.fault_attr = 3` driver-side wiring; trust that
  Gmsh's `Physical Surface` tagging is correct (the mesh file itself was
  not re-verified this round).
- `FaultBasis` frame computation (`fault_basis.hpp`) — uses vertex
  COORDINATES to compute tangent vectors at fault QPs. This is
  geometry-derivative, not fault-IDENTIFICATION; out of scope for the
  user's question.
- Future SEAS dynamic-mesh / AMR support — see R-1702.
