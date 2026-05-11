# Code Review: 2026-05-10 — round 3, post-fix audit

## Review Scope
- Plan: `miniapps/seas/safs/project_7.0_alternative/document/heterogeneous_material_plan.md` (after round-2 fixes)
- Files reviewed (re-read fresh, all three passes):
  - `miniapps/seas/safs/project_7.0_alternative/document/heterogeneous_material_plan.md`
  - `miniapps/seas/io/field_coefficient.cpp` (R-009 round-2 fix)
- Domain context consulted:
  - `miniapps/seas/dynamic/godunov_flux.hpp` (confirmed 6 internal `DenseMatrix` members per `GodunovFlux`: `Ax_`, `Ax_plus_`, `Ax_minus_`, plus `std::array<DenseMatrix, 3> ref_star_`)
  - `miniapps/seas/dynamic/wave_operator.hpp` / `wave_operator.inl` (`template <typename MeshType = Mesh>` design)
  - Previous `REVIEW.md` (round 2, overwritten — but found via conversation history)

This pass treats the post-fix state as fresh and hunts for NEW issues introduced or left latent.

## Findings

### [R-001] CRITICAL `heterogeneous_material_plan.md` Phase 3 §Goal / §Acceptance Criteria — memory budget arithmetic understates `GodunovFlux` size by ~6x; the acceptance assertion uses `sizeof(GodunovFlux)` which IGNORES the heap-allocated DenseMatrix data

**Category:** BUG

**Description:**
The plan justifies "Option (a) — caching one `GodunovFlux` per element" with this estimate (line 296):
> *"per-element matrices are small (9×9 doubles = 648 bytes; 1 M elements = 648 MB which is acceptable for an explicit DG run)"*

This treats `GodunovFlux` as if it holds ONE 9×9 matrix. But the actual class (verified in `godunov_flux.hpp:235-243`) holds **six** 9×9 matrices:
```cpp
DenseMatrix Ax_;        // 9x9
DenseMatrix Ax_plus_;   // 9x9
DenseMatrix Ax_minus_;  // 9x9
std::array<DenseMatrix, 3> ref_star_;  // three more 9x9
```
Plus scalars (`lambda_, mu_, rho_, cp_, cs_, Zp_, Zs_`).

Per-instance heap cost: 6 × 81 × 8 bytes ≈ **3.9 KB per GodunovFlux** (just the matrix payload).  Each `mfem::DenseMatrix` also carries a small object header (capacity + height/width + data pointer ≈ 32 bytes per matrix).  Realistic per-instance total ≈ **4 KB**.

For the SAFS z-graded 3.2 M-tet fixture, even with aggressive dedup down to 1 M unique fluxes (an optimistic factor of 3), memory = **4 GB**, NOT 1 GB.  For zero dedup (worst case in basin layers where every element has a unique sidecar sample), the cost is **12+ GB**.  The plan's "648 MB / 1 M elements" justification for picking Option (a) over Option (b) is off by 6x.

Worse, the acceptance test at line 439 reads:
> *"`pool.NumUniqueFluxes() * sizeof(GodunovFlux) <= 1e9` is the concrete assertion."*

`sizeof(GodunovFlux)` returns only the in-struct layout (`DenseMatrix` object headers + scalar members), NOT the heap-allocated `data` arrays inside each `DenseMatrix`.  In practice `sizeof(GodunovFlux)` is on the order of 400 bytes, NOT 4000.  The assertion as specified ALLOWS `NumUniqueFluxes()` up to ~2.5 M (1e9 / 400) — corresponding to **10 GB of actual memory** — and silently passes.  An implementer running the test will believe the budget is satisfied while real memory consumption is an order of magnitude over.

**Trigger:** Any Mode::Coefficient / Mode::GridFunction run on a real SAFS-scale mesh.  The acceptance test passes (`sizeof`-based check) and the real memory blows past 1 GB and possibly past the host's RAM.

**Actual behavior:** Test passes vacuously; real memory budget overshoots by 5-10x.

**Expected behavior:** Either (a) re-derive a realistic budget that accounts for the 6 matrices × heap payload, and update both the justification text and the acceptance criterion to use a per-instance constant like `kGodunovFluxBytes = 4096`, or (b) drop the cached `ref_star_[3]` triple in the pool (these are computed for ADER and can be re-derived from `Ax_`/`Ay_`/`Az_` on demand) — that halves the per-instance cost.

**Suggested fix:**
```diff
- This phase is structurally larger than Phase 2 because `GodunovFlux` was designed around a SINGLE constant `(λ, μ, ρ)` triple — its `Ax`, `Ax_plus_`, `Ax_minus_` matrices are precomputed once at construction. Going heterogeneous means either (a) caching one `GodunovFlux` per element, or (b) computing the 9×9 Jacobians per element on the fly. **Option (a)** is chosen because the per-element matrices are small (9×9 doubles = 648 bytes; 1 M elements = 648 MB which is acceptable for an explicit DG run) and the lookup is then a flat array index, preserving the hot-path access pattern.
+ This phase is structurally larger than Phase 2 because `GodunovFlux` was designed around a SINGLE constant `(λ, μ, ρ)` triple — its `Ax_`, `Ax_plus_`, `Ax_minus_` matrices (and the 3 cached `ref_star_` matrices used by the ADER recursion) are precomputed once at construction.  Each `GodunovFlux` instance therefore owns SIX 9×9 `DenseMatrix` objects, each carrying an 81-double heap payload — roughly **4 KB per instance** (R-001 round-3 — earlier "648 bytes per matrix; 1 M elements = 648 MB" estimate was off by 6× because only ONE matrix was counted).  Going heterogeneous means either (a) caching one `GodunovFlux` per element, or (b) computing the 9×9 Jacobians per element on the fly.  **Option (a)** is chosen because the lookup is a flat array index, preserving the hot-path access pattern, AND because the uniqueness map (Detailed Req. 2) typically dedups 5–10× across smoothly-varying layers, bringing realistic memory down to **400 MB – 2 GB** for the SAFS 3.2 M-tet fixture.

- - [ ] **Memory budget (R-011):** on the SAFS z-graded 3.2 M-tet fixture with `Mode::Coefficient` and the 39-slice sidecar, `GodunovFluxPool` allocates AT MOST 1.0 GB total for `flux_storage_`.  If exceeded, the constructor's `dedup_sig_figs` is loosened (e.g., from 6 to 4) and re-run.  `pool.NumUniqueFluxes() * sizeof(GodunovFlux) <= 1e9` is the concrete assertion.
+ - [ ] **Memory budget (R-011 / R-001 round-3):** on the SAFS z-graded 3.2 M-tet fixture with `Mode::Coefficient` and the 39-slice sidecar, `GodunovFluxPool` allocates AT MOST 2 GB total of HEAP-INCLUSIVE memory (NOT `sizeof(GodunovFlux)`, which ignores the heap-allocated `DenseMatrix::data` arrays inside each instance).  The concrete assertion is:
+   ```cpp
+   // Realistic per-instance bytes including the 6 × (9×9 × sizeof(double)) heap payload.
+   constexpr size_t kGodunovFluxBytes = 6 * 81 * sizeof(double) + /*headers + scalars*/ 256;
+   assert(pool.NumUniqueFluxes() * kGodunovFluxBytes <= 2'000'000'000);
+   ```
+   If exceeded, either (a) loosen `dedup_sig_figs` (from 6 to 4, etc.), or (b) re-engineer the pool to drop the cached `ref_star_[3]` matrices on heterogeneous-mode paths (they are re-derivable from `Ax_`/`Ay_`/`Az_` per call at modest CPU cost) and revisit Option (b) — JIT Jacobian build.
```

**Test case:**
```python
def test_R001_memory_budget_uses_heap_inclusive_bytes():
    pool = GodunovFluxPool[ParMesh](mat_coef, safs_3M_pmesh, dedup_sig_figs=6)
    n_unique = pool.NumUniqueFluxes()
    # NOT sizeof(GodunovFlux) — that ignores the 6 × 648 bytes of heap.
    bytes_per_instance = 6 * 81 * 8 + 256
    assert n_unique * bytes_per_instance <= 2_000_000_000
```

---

### [R-002] MODERATE `heterogeneous_material_plan.md` Phase 3 §Files to Create — `GodunovFluxPool` ctor takes `const MeshType& mesh` but Detailed Req. 7a (option (a)) builds a `ParFiniteElementSpace` that requires non-const `ParMesh*`

**Category:** BUG (re-emerged after round-2 fix)

**Description:**
The pool ctor signature (lines 310-312):
```cpp
GodunovFluxPool(const MaterialField& material,
                const MeshType& mesh,
                int dedup_sig_figs = 6);
```
is `const MeshType&`. But Detailed Req. 7a option (a) — the recommended option — instructs:
> *"Build a **new** `mfem::L2_FECollection` of order 0 (DG0) and a **new** `mfem::ParFiniteElementSpace` over the supplied `ParMesh`."*

`mfem::ParFiniteElementSpace` ctor takes `ParMesh*` (NON-const, see mfem `fespace.hpp`).  Building one from a `const ParMesh&` requires `const_cast<ParMesh*>(&mesh)` — exactly the same const-correctness footgun that round-1 R-013 removed from `ComputeMeshBBoxParallel`.  An implementer will either:
1. Add the `const_cast` (compiles, but inherits the same maintenance hazard the prior review explicitly called out).
2. Or change `const MeshType&` → `MeshType&`, and inherit cascading non-const correctness changes through the call site in `WaveOperator`.

The plan doesn't pick.

**Trigger:** Phase 3 implementation; the first attempt to build the temporary DG0 ParFES inside the pool ctor.

**Actual behavior:** Either compile error (if the implementer doesn't use `const_cast`) or a duplicated R-013 footgun.

**Expected behavior:** Plan picks a consistent const policy. The cleanest is to drop the `const` on the ctor's mesh argument — `MeshType& mesh` — because building the DG0 ParFES is a legitimate, expected use of the mesh and the rest of `WaveOperator` already holds the mesh as a `MeshType&` reference.

**Suggested fix:**
```diff
   template <typename MeshType = mfem::Mesh>
   class GodunovFluxPool
   {
   public:
       GodunovFluxPool(const MaterialField& material,
-                      const MeshType& mesh,
+                      MeshType& mesh,
                       int dedup_sig_figs = 6);
       ...
   };
```
Update the WaveOperator ctor's call site accordingly — `flux_pool_ = std::make_unique<GodunovFluxPool<MeshType>>(material, mesh_, ...)` — where `mesh_` is the operator's existing non-const `MeshType&` member.

**Test case:**
```python
def test_R002_pool_ctor_compiles_for_parmesh_without_const_cast():
    # Pool ctor signature must accept a ParMesh& without needing const_cast.
    pmesh = ParMesh(...)
    pool = GodunovFluxPool[ParMesh](material, pmesh)   # no const_cast
    # Bug pre-fix: code review shows const_cast required.  Post-fix: clean.
    assert pool.NumUniqueFluxes() > 0
```

---

### [R-003] MODERATE `heterogeneous_material_plan.md` Phase 3 §Edge Cases vs §Detailed Req. 5 — internal contradiction about MPI `dt` reduction

**Category:** DEVIATION (internal inconsistency)

**Description:**
§Detailed Req. 5 (lines 372-391) explicitly REPLACES the legacy MPI `h_min_` reduction with a new per-call MPI_Allreduce of `dt_local` inside `ComputeMaxDt`:
> *"a per-element walk **followed by an MPI_Allreduce(MIN) of `dt_local`** (the existing code did NOT reduce inside `ComputeMaxDt`; it only reduced `h_min_` once at construction)"*

But §Edge Cases (line 430) still says:
> *"`MaxCp()` does an `MPI_Allreduce(MAX)` over local maxima; `MIN dt` reduction over `dt_local` continues to use the existing MPI flow."*

The phrase "existing MPI flow" is stale relative to the new Det. Req. 5 — the existing flow at lines 111-114 of `wave_operator.inl` reduces `h_min_`, not `dt_local`.  An implementer reading §Edge Cases in isolation may think no new MPI_Allreduce is needed in `ComputeMaxDt`, contradicting Det. Req. 5.

**Trigger:** Phase 3 implementer reads §Edge Cases for MPI semantics and incorrectly believes the existing MPI flow already covers the new code.

**Actual behavior:** Implementer either does nothing in `ComputeMaxDt` (breaking parallel CFL), or implements the Det. Req. 5 reduction but is confused by the contradiction.

**Expected behavior:** The two sections agree.

**Suggested fix:**
```diff
- - **MPI element ownership.** The pool indexes by LOCAL element ID. Each MPI rank constructs its own pool from its local elements. `MaxCp()` does an `MPI_Allreduce(MAX)` over local maxima; `MIN dt` reduction over `dt_local` continues to use the existing MPI flow.
+ - **MPI element ownership.** The pool indexes by LOCAL element ID. Each MPI rank constructs its own pool from its local elements. `MaxCp()` does an `MPI_Allreduce(MAX)` over local maxima (only in the parallel specialisation; serial `GodunovFluxPool<Mesh>` returns the local maximum directly).  The `MIN dt` reduction over `dt_local` lives in `WaveOperator::ComputeMaxDt` per Detailed Req. 5 — a new explicit `MPI_Allreduce(MIN)` call inside that function (the legacy reduction of `h_min_` at construction is REMOVED).
```

**Test case:**
```python
def test_R003_edge_cases_section_matches_detailed_req_5():
    # Static review of the plan: both sections must reference the
    # SAME post-Phase-3 MPI flow (per-call MPI_Allreduce inside
    # ComputeMaxDt, NOT a one-shot ctor-time reduce of h_min_).
    plan_text = open("...heterogeneous_material_plan.md").read()
    # Det. Req. 5 wording
    assert "MPI_Allreduce(MIN) of `dt_local`" in plan_text
    # Edge Cases wording must match (not "existing MPI flow")
    assert "existing MPI flow" not in plan_text  # stale wording removed
```

---

### [R-004] MODERATE `heterogeneous_material_plan.md` Phase 3 §Files to Create — `MaxCp()` MPI reduction is unconditional but the serial `GodunovFluxPool<Mesh>` specialisation has no comm

**Category:** BUG

**Description:**
The Phase 3 pool interface promises:
> *"`real_t MaxCp() const;                                    // global maximum, MPI-reduced"*

And §Detailed Req. 4:
> *"`MaxCp()` is computed and MPI-reduced once at construction so callers can use it without per-call MPI."*

But the pool is templated `<typename MeshType = mfem::Mesh>`.  For `GodunovFluxPool<mfem::Mesh>` (serial specialisation), there is no `comm` and no `MPI_Allreduce`.  The plan doesn't guard the MPI call with an `if constexpr (IsParallelMesh<MeshType>::value)` block — an implementer following the spec literally writes uncompilable code in the serial specialisation (`mesh.GetComm()` not a member of `Mesh`).

**Trigger:** Phase 3 implementation; first build of the serial pool specialisation.

**Actual behavior:** Build error in `GodunovFluxPool<Mesh>::MaxCp()`.

**Expected behavior:** Plan explicitly says MaxCp's MPI reduction is gated on the parallel specialisation.

**Suggested fix:** Update §Files to Create note and §Detailed Req. 4:
```diff
   real_t MaxCp() const;                                    // global maximum, MPI-reduced
+                                                           // (parallel specialisation only;
+                                                           //  serial returns local max).
```
And Det. Req. 4:
```diff
 4. **CFL helper**:
    ```cpp
    real_t GodunovFluxPool::CpForElement(int elem) const
    { return At(elem).GetCp(); }
    real_t GodunovFluxPool::MaxCp() const
-   { /* return max over flux_storage_; MPI_Allreduce(MAX) at construction time */ }
+   { /* return cached `max_cp_` field.  Computed at ctor as max over
+      * flux_storage_; if MeshType==ParMesh (constexpr-if), the value
+      * is MPI_Allreduce(MAX) reduced across ranks at construction.
+      * For MeshType==Mesh the local maximum is returned directly. */ }
    ```
```

**Test case:**
```python
def test_R004_serial_pool_compiles_and_max_cp_returns_local_max():
    # Serial path must compile and produce the local pool's max c_p.
    pool = GodunovFluxPool[Mesh](mat_const, serial_mesh)
    cp_max = pool.MaxCp()
    assert cp_max > 0
    # In serial, no MPI calls.  Verify by linking against an MPI-free
    # MFEM build path (or by checking MaxCp() does not appear in
    # the compiled binary's MPI symbol references on the serial test).
```

---

### [R-005] LOW `field_coefficient.cpp::ComputeMeshBBoxParallel` — R-009 abort message prints only `xmin`/`xmax` even when the failing axis is y or z

**Category:** QUALITY

**Description:**
The R-009 round-2 guard reads:
```cpp
if (xmin > xmax || ymin > ymax || zmin > zmax)
{
   MFEM_ABORT(
      "FieldProjector::ComputeMeshBBoxParallel: global mesh bbox "
      "is degenerate (xmin=" << xmin << " > xmax=" << xmax << ", "
      "or analogous in y/z).  This typically indicates the mesh "
      "has zero elements on every MPI rank — check that the mesh "
      "file is non-empty and that the ParMesh partitioning did "
      "not silently drop all elements.");
}
```

If only the y-axis triggers (e.g., a contrived rank-distribution that produces `xmin < xmax` but `ymin = +inf, ymax = -inf`), the abort message says `"xmin=X > xmax=Y"` where `X <= Y` — confusing.  The phrase "or analogous in y/z" hints at the issue but doesn't print actual y/z values.

**Trigger:** Hypothetical degenerate mesh where only y- or z-axis fails.

**Suggested fix:** Print all six bounds in the abort message and let the user identify which axis failed:
```diff
   if (xmin > xmax || ymin > ymax || zmin > zmax)
   {
      MFEM_ABORT(
         "FieldProjector::ComputeMeshBBoxParallel: global mesh bbox "
-        "is degenerate (xmin=" << xmin << " > xmax=" << xmax << ", "
-        "or analogous in y/z).  This typically indicates the mesh "
+        "is degenerate.  Global bounds: "
+        "x=[" << xmin << ", " << xmax << "], "
+        "y=[" << ymin << ", " << ymax << "], "
+        "z=[" << zmin << ", " << zmax << "].  "
+        "(A degenerate axis has min > max.)  This typically indicates the mesh "
         "has zero elements on every MPI rank — check that the mesh "
         "file is non-empty and that the ParMesh partitioning did "
         "not silently drop all elements.");
   }
```

**Test case:** N/A (diagnostic message quality).

---

### [R-006] LOW `heterogeneous_material_plan.md` Phase 1 §Files to Modify — `MakeGridFunction` tightening (round-2 Det. Req. 7) is not listed in §Files to Modify

**Category:** QUALITY (visibility)

**Description:**
The round-2 fix added a new Phase 1 Detailed Requirement 7 (lines ~140-160) that tightens the existing `MakeGridFunction` factory by adding an `MFEM_VERIFY` for null shared_ptrs.  But §Files to Modify (lines 52-64) only mentions:
- `Mode::Coefficient = 2` enum addition
- three new `Coefficient*` members
- `MakeCoefficient` factory
- `EvalAt` accessor
- `MaxCpInElement` extension
- delete-inline-body directive

The `MakeGridFunction` tightening is buried in Detailed Requirements and not surfaced in the file-level summary.  An implementer skimming §Files to Modify and jumping to §Acceptance Criteria may miss it entirely.

**Trigger:** Implementer reads Phase 1 §Files to Modify, doesn't read all of §Detailed Requirements.

**Suggested fix:**
```diff
   - Extend `MaxCpInElement(int elem)` to accept an optional `ElementTransformation*` ...
   - **The existing header-only inline body of `MaxCpInElement` MUST BE DELETED** (R-004). ...
+  - Tighten the existing `MakeGridFunction` factory with an `MFEM_VERIFY` that all three shared_ptr arguments are non-null (R-004 round-2; details in Detailed Req. 7).  The body remains header-inline.
```

**Test case:** N/A (visibility / specification organisation).

---

### [R-007] LOW `heterogeneous_material_plan.md` Phase 4 §Detailed Req. 6 — audit step refers callers to `MaterialField::EvalAt(elem, T, ip)` but does not say how to obtain `(elem, T, ip)` for a FAULT DOF

**Category:** ASSUMPTION

**Description:**
The audit step (added in round 2 to fix R-003) tells the implementer:
> *"Any match that uses the returned value in a numerical formula (friction impedance, pre-stress amplitude, CFL helpers, etc.) MUST be rewritten to evaluate the operator's `MaterialField` at the relevant location (fault QP, element centroid, etc.) when the active mode is `Mode::Coefficient` or `Mode::GridFunction`."*

But `MaterialField::EvalAt(int elem, ElementTransformation& T, const IntegrationPoint& ip, ...)` requires a 3-tuple.  For a **fault DOF** (the typical friction-init context), the DOF is a 2D point on the fault surface; its corresponding 3D `(elem, T, ip)` triple is non-obvious:
1. Find the bulk element that owns the fault QP (could be on either side of the fault).
2. Get that element's `ElementTransformation` (`mesh.GetElementTransformation(elem_idx)`).
3. Compute the local `IntegrationPoint` for the fault QP (requires inverse-mapping the fault QP's physical coords back to reference coords of the bulk element).

This is feasible but non-trivial.  The plan should either provide a helper (`MaterialField::EvalAtFaultQP(fault_qp_index, ...)`) or point the implementer to existing SEAS code that does the lookup (`FaultBasis::GetFaultQP3DCoords()` or similar).

**Trigger:** Phase 4 implementer writes the friction-impedance audit fix and stalls on step 3 above.

**Suggested fix:** Add a sub-bullet to Phase 4 §Detailed Req. 6:
```diff
   The expected audit set today includes:
   - friction state init in `solver/seas_operator.hpp`
   - pre-stress / nucleation in `drivers/tpv102_driver.cpp`, `tpv104_driver.cpp`, `tpv205_driver.cpp`
   - any utility in `friction/` that hard-codes a Constant-coefficient model
+
+  Practical recipe for fault-DOF material lookup (the most common
+  consumer): given a fault QP's 3D physical coords `(x, y, z)` and
+  the local bulk element index `elem`, do
+  ```cpp
+  ElementTransformation* T = mesh.GetElementTransformation(elem);
+  IntegrationPoint ip;
+  T->TransformBack(Vector{x, y, z}, ip);   // inverse-map physical -> reference
+  material.EvalAt(elem, *T, ip, lambda, mu, rho);
+  real_t eta_p = std::sqrt((lambda + 2*mu) * rho);
+  ```
+  The fault basis already stores per-QP `(elem, x, y, z)` in
+  `FaultBasis::GetFaultQPCoords()`; reuse that for the inverse-map
+  inputs.
```

**Test case:** N/A (specification clarification).

---

### [R-008] POSSIBLE LOW `heterogeneous_material_plan.md` Phase 3 §Edge Cases — "Empty rank" note says `MaxCp()` initialises to `0.0`, but `0.0` will collide with finite c_p values on `MPI_Allreduce(MAX)`

**Category:** EDGE_CASE

**Description:**
§Edge Cases line 431:
> *"`MaxCp()` initialises to `0.0` and the `MPI_Allreduce(MAX)` recovers the correct global value."*

If a rank is empty, its local `max_cp = 0.0` is sent to `MPI_Allreduce(MAX)`.  The reduce picks `max(0.0, finite_c_p_on_other_ranks)` = the finite value. ✓

But the same comment also could mislead a reader: 0.0 is a VALID c_p value (a degenerate void medium), and 0.0 ≤ any positive c_p, so the MAX reduction always discards 0.0 in favour of any finite value.  Could `MaxCp()` legitimately return 0.0?  Only if EVERY rank has empty pools — same all-empty-mesh case that R-009 in `field_coefficient.cpp` already aborts on.  So the case is unreachable in practice IF the bbox check fires first.  But the plan doesn't cross-reference, and an implementer auditing the empty case may worry.

**Trigger:** Hypothetical all-empty mesh that somehow escapes the bbox check.

**Suggested fix:** Cross-reference the bbox guard, or use `-infinity` as the sentinel:
```diff
- **Empty rank.** If a rank has zero local elements (unusual but possible at high MPI counts), `MaxCp()` initialises to `0.0` and the `MPI_Allreduce(MAX)` recovers the correct global value.
+ **Empty rank.** If a rank has zero local elements (unusual but possible at high MPI counts), `MaxCp()` initialises to `-std::numeric_limits<real_t>::infinity()` and the `MPI_Allreduce(MAX)` recovers the correct global value from non-empty ranks.  Using `-inf` (instead of `0.0`) ensures any legitimate `c_p = 0.0` value would still survive the reduce; in practice an all-empty global mesh is impossible because `ComputeMeshBBoxParallel` aborts upstream (`field_coefficient.cpp` R-009 guard).
```

**Test case:** N/A (defensive sentinel choice).

---

### [R-009] LOW `heterogeneous_material_plan.md` Phase 4 §LoadSidecarMaterialBundle — abort message must name the failing field, but the plan doesn't reference the existing helper that already does so

**Category:** QUALITY

**Description:**
§Implementation step 3 (post-fix) reads:
> *"Run `DataField3D::ContainsBBox` on EACH of the three loaded fields and abort if any of the three fails the containment check (R-006 round-2 ...).  The abort message must name which field's bbox failed."*

The existing `FieldProjector::AbortContainmentFailure(field, ...)` in `field_coefficient.cpp` already takes a `DataField3D&` and prints `field.FieldName()` in the message.  The plan should direct the implementer to reuse this helper rather than write a new abort path.

**Trigger:** Phase 4 implementer rolls their own abort message and loses the field-name detail.

**Suggested fix:**
```diff
- 3. Compute the mesh bbox via `FieldProjector::ComputeMeshBBoxParallel`.  Run `DataField3D::ContainsBBox` on EACH of the three loaded fields and abort if any of the three fails the containment check ... The abort message must name which field's bbox failed.
+ 3. Compute the mesh bbox via `FieldProjector::ComputeMeshBBoxParallel`.  For each of the three loaded fields, if `DataField3D::ContainsBBox(...)` returns false, call the existing `FieldProjector::AbortContainmentFailure(field, mxmin, ..., mzmax)` helper (in `field_coefficient.cpp`).  That helper prints both the mesh bbox table and the field's bbox table with the field name (via `field.FieldName()`), giving the user actionable diagnostic info without rolling a new abort path.
```

**Test case:** N/A (specification clarification).

---

### [R-010] LOW `heterogeneous_material_plan.md` Phase 3 §Detailed Req. 2 — `dedup_sig_figs` rounding algorithm still unspecified

**Category:** ASSUMPTION (remains from round-2 R-005 LOW)

**Description:**
The round-2 fix introduced `dedup_sig_figs = 6` as a tunable knob but does not state how to round a `real_t` value to N significant figures.  Two common approaches:
1. `std::round(value / std::pow(10.0, std::floor(std::log10(std::abs(value))) + 1 - sig_figs)) * std::pow(...)`.
2. `snprintf` to a temp buffer with `%.6g` format, parse back.

The two are NOT bit-equivalent at boundary values (e.g., `3.999999` rounds to `4.00000` in `%g` but to `3.99999` in numeric form for 6 sig figs).  Different implementations dedup differently → `NumUniqueFluxes()` is implementation-defined → test failures across compilers.

**Trigger:** Different compilers / platforms produce different `NumUniqueFluxes()` for the same input.

**Suggested fix:**
```diff
 2. **Memory optimisation — uniqueness map** (R-011): ...  Rounding precision is controlled by the constructor's `dedup_sig_figs` parameter (default 6).  ...
+
+    Rounding algorithm:
+    ```cpp
+    inline real_t round_to_sig_figs(real_t v, int sig_figs)
+    {
+       if (v == 0.0 || std::isnan(v)) return v;
+       const real_t magnitude = std::pow(10.0,
+          std::floor(std::log10(std::abs(v))) + 1 - sig_figs);
+       return std::round(v / magnitude) * magnitude;
+    }
+    ```
+    Hash key is the std::tuple<real_t, real_t, real_t> of the rounded
+    triple, using `std::hash<real_t>` element-wise XOR-mixed (or
+    `boost::hash_combine` if MFEM uses it elsewhere).
```

**Test case:** N/A (specification clarification).

---

## Summary
- Critical issues: **1** (R-001 — 6x memory underestimate + sizeof-based assertion that ignores heap; the budget acceptance test passes vacuously)
- Moderate issues: **3** (R-002 const-correctness on pool ctor; R-003 internal contradiction between Edge Cases and Det. Req. 5; R-004 MaxCp MPI in serial specialisation)
- Low issues: **6** (R-005 abort-message diagnostics; R-006 MakeGridFunction not listed in §Files to Modify; R-007 fault-DOF EvalAt recipe; R-008 MaxCp sentinel value; R-009 reuse existing AbortContainmentFailure; R-010 dedup rounding algorithm)
- Plan compliance: **PARTIAL** — round-2 CRITICAL/MODERATE findings (R-001..R-004 from the previous round) are correctly addressed; the post-fix state introduces 1 new CRITICAL (memory math) and 3 new MODERATE issues (mostly spec-internal contradictions and unhandled serial-specialisation gaps).
- Verdict: **PASS WITH FIXES** — R-001 must be resolved (the memory acceptance test is misleading and could let implementers ship a 10x-budget overrun); R-002..R-004 are smaller spec gaps that the implementer can plausibly catch but should be tightened before Phase 3 starts.

### Why PASS WITH FIXES (not FAIL)
- The single CRITICAL (R-001 round-3) is a memory accounting error confined to the plan document; no code is yet implemented; it surfaces during the first SAFS-scale acceptance test, NOT silently in production.
- R-002 (const_cast) is the same hazard the prior round-1 review explicitly flagged; the fix-round successfully fixed it once in `field_coefficient.cpp` but it re-emerged in the Phase 3 spec — the implementer can recognise the pattern and fix it again.
- R-003 / R-004 are textual contradictions in the plan; an implementer who reads top-to-bottom in order will follow the Detailed Requirements (latest authoritative spec) and skip the stale Edge Cases note.

### Why this is the strictest standard so far
With three review-fix rounds in, the remaining findings are concentrated in the plan document and have no compiled code impact yet. The plan's net trajectory is improving (round 1: 3 CRITICAL → round 2: 1 CRITICAL → round 3: 1 CRITICAL), but the round-3 CRITICAL is a NEW class of issue (a memory budget that is wrong by 6x, with a `sizeof`-based assertion that hides the discrepancy). Catching this before Phase 3 implementation begins saves a debug cycle later.

## Unreviewed Areas
- `MaterialField::EvalAt` body in `Mode::Coefficient` (not yet implemented; covered by Phase 1 spec).
- Phase 4 driver wiring scaffolding (not yet implemented; covered by Phase 4 spec).
- The actual `data_field_3d.cpp` trilinear evaluator (already verified bit-perfect earlier).
- Plot scripts under `code_preprocess/data_projection/` (research artifacts, not in scope).
- `FaultFaceFlux::InitializeImpedancesFromMaterial` (explicitly deferred by the plan).
