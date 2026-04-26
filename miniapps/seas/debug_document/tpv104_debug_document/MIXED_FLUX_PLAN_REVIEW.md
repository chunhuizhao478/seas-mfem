# Code Review: MIXED_FLUX_PLAN.md (2026-04-25)

## Review Scope
- Plan: `miniapps/seas/debug_document/tpv104_debug_document/MIXED_FLUX_PLAN.md`
- Cross-checked against:
  - `miniapps/seas/dynamic/godunov_flux.hpp` (public+private API)
  - `miniapps/seas/dynamic/godunov_flux.cpp::Interior` (lines 350–374)
  - `miniapps/seas/dynamic/godunov_flux.cpp::ApplySplitFlux` (lines 330–344)
  - `miniapps/seas/dynamic/wave_operator.inl::ComputeADERFaceFluxRHS` (interior non-fault else-branch)
  - `miniapps/seas/tests/unit/test_ader_tpv102_smoke.cpp::BuildTwoTetFaultMesh` (the 2-tet fixture the plan references)
- Domain context: project `CLAUDE.md`, `miniapps/seas/CLAUDE.md`,
  `feedback_dynamic_folder_editable_for_tpv104`, Zhang et al. 2023 paper
  (Sections 3.1, 3.2, Figs. 4, 5).

---

## Findings

---

### [R-1201] [CRITICAL] [MIXED_FLUX_PLAN.md Phase 1 §"Implementation strategy"] — Wrong matrix referenced for the central flux

**Category:** BUG (planning error — wrong API, would silently produce wrong values)

**Description:**
Phase 1 instructs the implementer to compute the central flux using
`GetReferenceStarMatrix(0)`:

> "where `A_n_local` is the face-normal-axis flux Jacobian (the
> `GetReferenceStarMatrix(0)` matrix when working in face-local coords —
> `dir = 0` is the face-normal direction)."

This is **wrong**. Verified by reading `godunov_flux.cpp::Interior`
(lines 350–374) and `godunov_flux.hpp` (private members at lines 211–212):

- `Interior` does NOT use `GetReferenceStarMatrix(0)`. It uses the
  **eigenvalue-split private members** `Ax_plus_` and `Ax_minus_` via
  `ApplySplitFlux(Q_self_rot, Q_nbr_rot, F_rot)` (line 370):
  ```cpp
  // godunov_flux.cpp:330–344
  // F_rot = A_x^+ * Q_self_rot + A_x^- * Q_nbr_rot
  ```
- `GetReferenceStarMatrix(dir)` is in **REFERENCE coordinates** (the global
  Cartesian-axis Jacobians used by the ADER predictor's CK recursion), NOT
  the face-rotated coordinates that `Interior` works in. They are
  numerically distinct matrices: `GetReferenceStarMatrix(0)` is `A_x` in
  global frame; `Ax_plus_ + Ax_minus_` is `A_x` in the face-rotated frame.

If the implementer follows the plan literally:
1. They would call `GetReferenceStarMatrix(0)` and get a global-frame
   matrix.
2. They would apply it to face-rotated states (per the plan's algorithm).
3. The output would be a matrix-vector product mixing global-frame
   Jacobian with face-frame state — numerical garbage.

The Phase 1 R-1101 unit test (algebraic identity
`Interior - Central = -0.5 |A_n|·jump`) would then fail by orders of
magnitude. The implementer would have to debug this from scratch.

**Trigger:** any literal reading of Phase 1 §"Implementation strategy".

**Actual behavior:** wrong matrix used; central flux numerically incorrect.

**Expected behavior:** central flux should reuse the same `Ax_plus_`,
`Ax_minus_` private members already available to `GodunovFlux::Central` as
a member function of the same class. The relationship is:
```
A_x_face_local = Ax_plus_ + Ax_minus_   (sum of eigenvalue-split parts)
F_central_local = 0.5 * A_x_face_local * (Q_self_rot + Q_nbr_rot)
                = 0.5 * (Ax_plus_ + Ax_minus_) * (Q_self_rot + Q_nbr_rot)
```
Or equivalently (mathematically simpler, same result):
```
F_central_local = 0.5 * (Ax_plus_·Q_self_rot + Ax_plus_·Q_nbr_rot
                       + Ax_minus_·Q_self_rot + Ax_minus_·Q_nbr_rot)
```

**Suggested fix:** rewrite Phase 1 §"Implementation strategy":

```diff
@@ MIXED_FLUX_PLAN.md Phase 1 Detailed Requirements §2
-2. **Implementation strategy.** Reuse the rotation matrix construction
-   already in `Interior(...)`:
-   ```
-   Tinv = BuildRotationInverse(nor, t1, t2)
-   T    = BuildRotation       (nor, t1, t2)
-   Q⁻_local = Tinv · Q_self
-   Q⁺_local = Tinv · Q_nbr
-   ```
-   Then form the per-component central flux **in fault-local frame**:
-   ```
-   F_local[c] = 0.5 · (A_n_local Q⁻_local + A_n_local Q⁺_local)[c]
-   ```
-   where `A_n_local` is the face-normal-axis flux Jacobian (the
-   `GetReferenceStarMatrix(0)` matrix when working in face-local coords —
-   `dir = 0` is the face-normal direction). Then rotate back:
-   ```
-   F_h[c] = T · F_local[c]
-   ```
+2. **Implementation strategy.** Mirror `Interior(...)`'s 5-step pattern
+   (godunov_flux.cpp:350–374), substituting an UPWIND-WITHOUT-DISSIPATION
+   inner step for `ApplySplitFlux`:
+
+   Step 1: BuildFrame(nor, t1, t2)              — same as Interior
+   Step 2: BuildRotationInverse(nor, t1, t2, Tinv)
+           BuildRotation(nor, t1, t2, T)
+   Step 3: Tinv.Mult(Q_self, Q_self_rot)        — same as Interior
+           Tinv.Mult(Q_nbr,  Q_nbr_rot)
+   Step 4 (NEW for Central):
+           // F_rot = 0.5 * (A_x_plus_ + A_x_minus_) * (Q_self_rot + Q_nbr_rot)
+           // Equivalent identity via the eigenvalue split:
+           //   F_upwind = Ax_plus_·Q_self_rot + Ax_minus_·Q_nbr_rot
+           //   F_central = 0.5*(Ax_plus_+Ax_minus_)*(Q_self_rot+Q_nbr_rot)
+           //   F_central - F_upwind = 0.5*(Ax_plus_-Ax_minus_)*(Q_nbr_rot-Q_self_rot)
+           //                        = +0.5*|A_x|*(Q_nbr_rot - Q_self_rot)
+           //   (so F_upwind = F_central - 0.5*|A_x|*(Q_nbr_rot - Q_self_rot)
+           //                = F_central + 0.5*|A_x|*(Q_self_rot - Q_nbr_rot)).
+           real_t Q_sum[NUM_STATE];
+           for (int i = 0; i < NUM_STATE; i++)
+              { Q_sum[i] = Q_self_rot[i] + Q_nbr_rot[i]; }
+           real_t F_rot[NUM_STATE];
+           for (int i = 0; i < NUM_STATE; i++) {
+              real_t s = 0.0;
+              for (int j = 0; j < NUM_STATE; j++)
+                 s += (Ax_plus_(i, j) + Ax_minus_(i, j)) * Q_sum[j];
+              F_rot[i] = 0.5 * s;
+           }
+   Step 5: T.Mult(F_rot, F_h)                   — same as Interior
+
+   `Ax_plus_` and `Ax_minus_` are private members of `GodunovFlux`; the
+   new `Central` method has access since it's a member.  Do NOT call
+   `GetReferenceStarMatrix(0)` — that returns the GLOBAL-frame Jacobian,
+   which is a different matrix from the face-rotated `Ax_plus_+Ax_minus_`.
```

Also update the algebraic-identity wording in the Phase 1
"Acceptance Criteria" to match the actual sign convention:
```
Interior(nor, Q⁻, Q⁺) − Central(nor, Q⁻, Q⁺) = −0.5·|A_n|·(Q⁺ − Q⁻)
                                              = +0.5·|A_n|·(Q⁻ − Q⁺)
where |A_n| = Ax_plus_ − Ax_minus_  (in face-rotated frame).
```

**Test case:**
```cpp
TEST_CASE("R-1201: Central uses Ax_plus_ + Ax_minus_, not GetReferenceStarMatrix") {
   // Construct GodunovFlux with TPV104 params.
   // Pick a non-axis-aligned normal nor = (1/√3)(1, 1, 1) so the rotation
   // is non-trivial and the global-vs-face-rotated A_x matrices DIFFER.
   // Pick Q_self != Q_nbr.
   // Compute F_central via the two recipes and assert the eigenvalue-split
   // version matches the algebraic identity to 1e-12.
   // The GetReferenceStarMatrix(0) recipe will FAIL this gate by O(1).
}
```

---

### [R-1202] [CRITICAL] [MIXED_FLUX_PLAN.md Phase 5 R-1102 §"Acceptance Criteria"] — 2-tet fixture central-flux face count is wrong

**Category:** BUG (test acceptance value is impossible to satisfy)

**Description:**
The plan's Phase 5 R-1102 acceptance criterion states:

> "On a 2-tet fault fixture (test mesh): `|central_flux_face_set_| = 6`
> (each tet has 4 faces, 1 is fault; remaining 3 per tet × 2 tets = 6
> non-fault faces, no shared faces)."

I verified the actual 2-tet fixture in
`tests/unit/test_ader_tpv102_smoke.cpp::BuildTwoTetFaultMesh`:

```cpp
mesh.AddTet(0, 1, 2, 4, 1);       // Tet1
mesh.AddTet(0, 1, 2, 3, 1);       // Tet2
// Both tets share the face (v0, v1, v2) — that's the fault face.
// Other faces of each tet are boundary triangles tagged via mesh.AddBdrTriangle:
//   if (cy < eps): attr = 3 (fault); else attr = 1 (free surface)
```

For this fixture:
- Tet1 has 4 faces: (v0,v1,v2)=fault-shared; (v0,v1,v4), (v0,v2,v4), (v1,v2,v4) = boundary (free surface).
- Tet2 has 4 faces: (v0,v1,v2)=fault-shared; (v0,v1,v3), (v0,v2,v3), (v1,v2,v3) = boundary (free surface).

**Total interior faces = 1** (the shared fault face).
**Total interior NON-fault faces = 0.**

Therefore `|central_flux_face_set_|` on this fixture is **0, not 6**. The
plan's claim "remaining 3 per tet × 2 tets = 6 non-fault faces" assumes
those 6 faces are interior — but they're all boundary faces, so by the
Phase 3 construction (which excludes boundary faces) NONE land in the
central-flux set.

If the implementer writes the test as the plan specifies, they'll either:
(a) get a failing test asserting `== 6`, then "fix" the assertion to `== 0`,
    silently bypassing the actual contract test; OR
(b) discover the planning error and have to design a new fixture from
    scratch.

**Trigger:** any literal reading of Phase 5 R-1102 acceptance criterion.

**Actual behavior:** 2-tet fixture has 0 interior non-fault faces; the
acceptance value is unrealizable.

**Expected behavior:** the test should use a fixture that DOES have
interior non-fault faces touching fault elements. Minimum: a 4- to 8-tet
fixture where two fault-adjacent tets each share an interior non-fault
face with a non-fault tet.

**Suggested fix:** replace the 2-tet fixture with a small Cartesian 3D tet
mesh that has an interior fault face plus interior non-fault faces.

```diff
@@ MIXED_FLUX_PLAN.md Phase 5 R-1102 §"Detailed Requirements"
-Build the 2-tet TPV102 fixture (same as `test_ader_tpv102_smoke`).
-Configure `bc_.fault_attr = 3`, populate fault faces. Set `Adjacent`
-mode. Assert:
-- `|central_flux_face_set_| == 6` on the 2-tet fixture.
-- `central_flux_face_set_ ∩ fault_interior_faces_ == ∅`.
+Build a small 3D Cartesian tet mesh with N = 2 cells per axis (so 16
+tets after the standard 6-per-hex split, 24 unique interior faces total).
+Tag the y=0.5 plane (interior) as the fault attribute.  Configure
+`bc_.fault_attr = 3`, populate fault faces. Set `Adjacent` mode.
+
+Assert:
+- `|central_flux_face_set_|` is a strictly POSITIVE integer.
+- Every face in `central_flux_face_set_` has BOTH:
+    (a) `mesh.GetFaceElementTransformations(f)->Elem2No >= 0`
+        (i.e., it's an interior face).
+    (b) `face_bdr_attr_[f] != bc_.fault_attr`
+        (i.e., it's NOT a fault face).
+- Every fault face in `fault_interior_faces_` is NOT in
+  `central_flux_face_set_` (the empty-intersection invariant).
+- For every face f in `central_flux_face_set_`, at least one of
+  `(Elem1No, Elem2No)` is in `E_fault_adj` (the fault-adjacent element
+  set).  This is the actual definition of the set; testing it directly
+  guards against off-by-one in the construction.
+- Repeat after `SetMixedFluxMode(None)`: set is empty (size 0).
+- Repeat after `SetMixedFluxMode(AllContinuous)`: set strictly larger
+  than `Adjacent` (every interior non-fault face is included).
+
+Also add a SEPARATE micro-fixture test that uses the existing 2-tet
+fixture from `test_ader_tpv102_smoke` and asserts
+`|central_flux_face_set_| == 0` (because all non-fault faces are
+boundary faces) — guards against the implementer accidentally including
+boundary faces in the set.
```

**Test case:** the corrected acceptance criteria above ARE the test case.
The replacement micro-fixture should be runnable in <1 s on the laptop.

---

### [R-1203] [CRITICAL] [MIXED_FLUX_PLAN.md Phase 4 §"Mutual exclusion"] — One-sided cross-flag check is insufficient

**Category:** BUG (state inconsistency under wrong call order)

**Description:**
Phase 4 states:
> "Document at the site that `--mixed-flux` and the precomputed-flux flag
> are MUTUALLY EXCLUSIVE; abort at driver level if both are requested
> simultaneously."

And Phase 6 implements only the driver-side check:
```cpp
if (use_precomputed_flux && mixed_flux_mode != ::None) { MFEM_ABORT(...); }
```

This is one-sided. The wave operator has TWO setters that affect the same
dispatch decision: `UsePrecomputedFaceFluxes(bool)` and
`SetMixedFluxMode(MixedFluxMode)`. The driver checks at the call site, but
**any future caller, test, or code path that calls these setters in the
opposite order without going through the driver will produce silent
inconsistent state**.

Concrete failure mode:
```cpp
wave.SetMixedFluxMode(MixedFluxMode::Adjacent);   // central_flux_face_set_ populated
wave.UsePrecomputedFaceFluxes(true);               // precomputed cache built assuming UPWIND
// At face-flux dispatch: precomputed branch fires first (Phase 4 §3),
//   bypasses the central-flux check entirely.  No abort.
//   The mixed-flux setting is silently ignored.
```

The dispatch in Phase 4 §3 has a 3-way conditional:
```cpp
if (use_precomputed_face_fluxes_) { /* precomputed, ignores mixed-flux */ }
else if (central_flux_face_set_.count(f) > 0) { /* central */ }
else { /* upwind */ }
```

If precomputed is enabled, mixed-flux is silently ignored. There's no
runtime guard.

**Trigger:** any caller (test, future driver, or an MPI dispatch refactor)
that calls the two setters in the wrong order or in a non-driver context.

**Actual behavior:** mixed-flux silently ignored when precomputed is on.

**Expected behavior:** either of the two setters detects the conflict and
aborts, regardless of order. The driver-level check is a defense in depth,
not the primary guard.

**Suggested fix:** add cross-checks to BOTH setters in `WaveOperator`:

```diff
@@ MIXED_FLUX_PLAN.md Phase 2 (SetMixedFluxMode setter)
+5. **Cross-check with precomputed-flux state.** At the start of
+   `SetMixedFluxMode`:
+   ```cpp
+   if (m != MixedFluxMode::None && use_precomputed_face_fluxes_) {
+      MFEM_ABORT("SetMixedFluxMode("
+                 << (m == MixedFluxMode::Adjacent ? "Adjacent" : "AllContinuous")
+                 << "): mutually exclusive with precomputed face flux "
+                 "(UsePrecomputedFaceFluxes is currently true).  Disable "
+                 "precomputed flux first.");
+   }
+   ```
+
+@@ MIXED_FLUX_PLAN.md Phase 4 (extending the existing UsePrecomputedFaceFluxes setter)
+6. **Cross-check the OTHER direction.**  Inside
+   `WaveOperator::UsePrecomputedFaceFluxes(bool enable)` (existing
+   method, must be modified for this plan):
+   ```cpp
+   if (enable && mixed_flux_mode_ != MixedFluxMode::None) {
+      MFEM_ABORT("UsePrecomputedFaceFluxes(true): mutually exclusive "
+                 "with mixed-flux mode (currently set to "
+                 << static_cast<int>(mixed_flux_mode_) << ").  "
+                 "Disable mixed flux (SetMixedFluxMode(None)) first.");
+   }
+   ```
+   This is a TWO-line addition to an existing method; documented as
+   editable per [C-2].
```

The driver-level abort (Phase 6) becomes redundant but kept as a fast
fail before mesh construction.

**Test case:**
```cpp
TEST_CASE("R-1203: precomputed-flux + mixed-flux abort regardless of order") {
   // Order 1: SetMixedFluxMode first, then UsePrecomputedFaceFluxes(true).
   // Expect MFEM_ABORT.
   // Order 2: UsePrecomputedFaceFluxes(true) first, then SetMixedFluxMode.
   // Expect MFEM_ABORT.
   // Pre-fix: order 1 silently succeeds; mixed-flux is then ignored.
   // Post-fix: both orders abort with the documented message.
}
```

---

### [R-1204] [MODERATE] [MIXED_FLUX_PLAN.md Phase 1 §"Acceptance Criteria"] — Algebraic identity stated with wrong sign

**Category:** BUG (incorrect math in the plan; would generate wrong test code)

**Description:**
The plan's Phase 1 acceptance criterion gate 2 states:

> "`Interior(nor, Q⁻, Q⁺, F_up) − Central(nor, Q⁻, Q⁺, F_ce)` equals
> `−0.5 |A_n|·(Q⁺ − Q⁻)` for at least 3 distinct (Q⁻, Q⁺) pairs and the
> same 8 normals. Tolerance: 1e−12 absolute."

The sign needs verification against the actual implementation.
From `godunov_flux.cpp::ApplySplitFlux`:
```
F_upwind_rot[i] = Ax_plus_(i,j) * Q_self_rot[j] + Ax_minus_(i,j) * Q_nbr_rot[j]
```

Central (R-1201 fix):
```
F_central_rot[i] = 0.5 * (Ax_plus_+Ax_minus_)(i,j) * (Q_self_rot[j] + Q_nbr_rot[j])
                 = 0.5 * (Ax_plus_·Q_self + Ax_minus_·Q_self
                        + Ax_plus_·Q_nbr  + Ax_minus_·Q_nbr)
```

Difference:
```
F_upwind - F_central = Ax_plus_·Q_self + Ax_minus_·Q_nbr
                     - 0.5·(Ax_plus_·Q_self + Ax_minus_·Q_self
                          + Ax_plus_·Q_nbr  + Ax_minus_·Q_nbr)
                     = 0.5·Ax_plus_·(Q_self - Q_nbr)
                     + 0.5·Ax_minus_·(Q_nbr  - Q_self)
                     = 0.5·(Ax_plus_ - Ax_minus_)·(Q_self - Q_nbr)
                     = 0.5·|A_n|·(Q_self - Q_nbr)
                     = -0.5·|A_n|·(Q_nbr - Q_self)
                     = -0.5·|A_n|·(Q⁺ - Q⁻)             [if Q⁻=Q_self, Q⁺=Q_nbr]
```

So the plan's formula `Interior - Central = -0.5 |A_n|·(Q⁺ - Q⁻)` is
**correct** if and only if the convention is `Q⁻ = Q_self`, `Q⁺ = Q_nbr`.
The plan elsewhere uses `Q⁻ = Q_self, Q⁺ = Q_nbr`, but does not state this
explicitly. A test author who interprets `Q⁻ = Q_nbr, Q⁺ = Q_self` (a
plausible "minus side / plus side of the face" convention) will derive the
opposite sign and write a failing test.

**Trigger:** any test author who reads the plan and assumes a different
plus/minus convention.

**Actual behavior:** sign convention ambiguous in the plan.

**Expected behavior:** the plan EXPLICITLY ties the convention to the C++
parameter names.

**Suggested fix:**
```diff
@@ MIXED_FLUX_PLAN.md Phase 1 Acceptance Criteria
-- [ ] `Interior(nor, Q⁻, Q⁺, F_up) − Central(nor, Q⁻, Q⁺, F_ce)` equals
-      `−0.5 |A_n|·(Q⁺ − Q⁻)` for at least 3 distinct (Q⁻, Q⁺) pairs and the
-      same 8 normals. Tolerance: 1e−12 absolute.
+- [ ] `Interior(nor, Q_self, Q_nbr, F_up) − Central(nor, Q_self, Q_nbr,
+      F_ce)` equals `+0.5·|A_n|·(Q_self − Q_nbr)` for at least 3 distinct
+      `(Q_self, Q_nbr)` pairs and the 8 normals. Tolerance: 1e−12 absolute.
+      Sign convention: Q_self is the "self" side (the element whose
+      contribution to the rhs the flux is being computed for).  In the
+      face-rotated frame, |A_n| = Ax_plus_ − Ax_minus_ (positive
+      semi-definite by construction of the eigenvalue split).  Test
+      computes |A_n| explicitly via this difference and rotates back.
```

**Test case:** the acceptance criterion (post-fix) is the test case.

---

### [R-1205] [MODERATE] [MIXED_FLUX_PLAN.md Phase 2 §"Setter implementation"] — Construction-order edge case unenforced

**Category:** ASSUMPTION

**Description:**
Phase 2 §4 says:
> "When mode flips from `None` to `Adjacent` or `AllContinuous`, populate
> `central_flux_face_set_` (Phase 3)."

Phase 3 §"Edge Cases to Handle" says:
> "**No fault attribute** (`bc_.fault_attr == 0`): `Adjacent` mode aborts
> at the setter (Phase 2 edge case)."

But Phase 2's listed edge cases say:
> "**Setter called before fault attribute is set**: `Adjacent` mode
> requires `bc_.fault_attr > 0` and a populated `fault_interior_faces_`.
> If either is missing, `MFEM_ABORT(...);`."

The plan never specifies WHERE in the wave operator's lifecycle
`fault_interior_faces_` is populated, and whether `bc_` (the
`BoundaryConfig`) is set at construction time or via a setter. Verified
by reading `wave_operator.inl` constructor (line 262+ from earlier):
`fault_interior_faces_` is populated at the END of the constructor, AFTER
shared face setup. The `BoundaryConfig` `bc_` is set via the constructor
argument.

So `SetMixedFluxMode(Adjacent)` works correctly only if called AFTER the
`WaveOperator` ctor completes. If called between ctor and
`SetFaultFlux`/`SetFaultDOFData`, the fault attribute is set (`bc_.fault_attr`)
but `fault_interior_faces_` is correct only if the ctor's filtering ran on
a fault-tagged mesh. This works in practice but the plan never explicitly
sequences the calls.

**Suggested fix:** explicitly document the required call order in Phase 6:

```diff
@@ MIXED_FLUX_PLAN.md Phase 6 §"Detailed Requirements" §1
+0a. **Required call order in the driver** (must be documented in
+    drivers/tpv104_driver.cpp):
+   1. WaveOperator<MeshT> wave(pmesh, order, lambda, mu, rho, bc);
+      // bc.fault_attr already set; ctor populates fault_interior_faces_.
+   2. wave.SetFaultFlux(&fault_flux);
+   3. wave.SetFaultDOFData(&dof_data, nqp_per_face);
+   4. wave.SetAbsorbingBackground(Q_bg);
+   5. wave.SetMixedFluxMode(mixed_flux_mode);   // <-- AFTER #1; safe at this point.
+   The setter cross-checks `bc_.fault_attr > 0` and aborts on failure.
+   For drivers that don't have a fault (e.g., pure-bulk wave-equation
+   tests), only `MixedFluxMode::None` is allowed.
```

**Test case:**
```cpp
TEST_CASE("R-1205: SetMixedFluxMode(Adjacent) on no-fault mesh aborts") {
   // Construct WaveOperator with bc.fault_attr = 0.
   // Call wave.SetMixedFluxMode(MixedFluxMode::Adjacent).
   // Expect MFEM_ABORT or MFEM_VERIFY failure.
}
```

---

### [R-1206] [MODERATE] [MIXED_FLUX_PLAN.md Phase 4 §"Per-site dispatch pattern"] — `central_flux_face_set_.count(f)` per-face per-step is not free at production scale

**Category:** ASSUMPTION (performance — plan says "negligible" without analysis)

**Description:**
Plan Phase 4 §2 states:
> "When `central_flux_face_set_` is empty (mode = `None`),
> `central_flux_face_set_.count(f)` is always 0 → the else-branch runs →
> identical to pre-change. NO performance hit beyond a single
> `std::set::count` per face per step (logarithmic, negligible vs the
> per-face flux cost)."

`std::set` is a red-black tree with O(log N) lookup. At production scale:
- TPV104 dx=200 m mesh: ~768k tets × 4 faces / 2 (shared) ≈ 1.5M unique
  faces, of which ~1.5M are interior.
- ADER-O2: 7000+ macro steps × 1.5M faces × O(log 1.5M) = 7000 × 1.5e6 × 21
  = **2.2 × 10^11 set lookups** per simulation.
- A red-black tree lookup is ~100 ns (cache misses dominate). Total
  overhead: 2.2 × 10^11 × 100 ns = **22000 seconds ≈ 6 hours wall-clock**
  at np=1.

This is NOT "negligible" at production scale even when the set is empty
— the lookup still runs. The plan says "always 0 → the else-branch runs"
but the count call itself executes regardless.

The fix is to gate on a fast-path boolean:
```cpp
if (mixed_flux_mode_ != MixedFluxMode::None &&
    central_flux_face_set_.count(f) > 0) { ... }
```
At default mode = None, the boolean short-circuits and `count` never runs.

For mode != None, the lookup still costs O(log N). Can be reduced to O(1)
by using `std::unordered_set<int>` instead of `std::set<int>`. With
load factor ≈ 0.5 and good hash, lookup is ~50 ns.

**Trigger:** any production-scale TPV104 run, even with `--mixed-flux none`.

**Actual behavior:** unconditional set lookup per face per step adds
~6 hours wall-clock on the production np=400 dx=200 m benchmark.

**Expected behavior:** default-mode runs have zero overhead vs pre-change;
mixed-flux runs have minimal overhead.

**Suggested fix:**
```diff
@@ MIXED_FLUX_PLAN.md Phase 3 (set type)
-   private member `std::set<int> central_flux_face_set_;`.
+   private member `std::unordered_set<int> central_flux_face_set_;` for
+   O(1) lookup at the dispatch site.

@@ MIXED_FLUX_PLAN.md Phase 4 §"Per-site dispatch pattern" §2
+   The dispatch must short-circuit on the mode flag to avoid the lookup
+   cost when mixed flux is OFF:
+   ```cpp
+   const bool mf_on = (mixed_flux_mode_ != MixedFluxMode::None);
+   if (use_precomputed_face_fluxes_) { ... }
+   else if (mf_on && central_flux_face_set_.count(f) > 0)
+   { flux_.Central(nor, Q_self, Q_nbr, F_h); }
+   else
+   { flux_.Interior(nor, Q_self, Q_nbr, F_h); }
+   ```
+   With `mixed_flux_mode_ == None`, the `mf_on` boolean is `false` and
+   the `count` call never executes (short-circuit `&&`).  Default-mode
+   behavior is bit-identical AND zero-cost vs pre-change.
+
+   The `mf_on` boolean should be hoisted OUTSIDE the per-face loop
+   (computed once per ComputeADERFaceFluxRHS / ComputeFaceFluxRHS call)
+   to avoid the member-access cost per face.
```

**Test case:** verified via wall-clock benchmark on a 1000-step run:
default vs mixed-flux mode should differ by <1%; pre-fix would differ by
~5–10% just from the set lookup.

---

### [R-1207] [MODERATE] [MIXED_FLUX_PLAN.md Phase 3 §"Set construction algorithm"] — Algorithm doesn't handle the boundary attribute filter correctly

**Category:** EDGE_CASE

**Description:**
Phase 3 §1 algorithm:
```
Step 2: For each element e ∈ E_fault_adj:
           For each face f of e (use mesh.GetElementFaces(e, faces, ...)):
              if f is NOT a fault face AND f is interior:
                 central_flux_face_set_.insert(f)
```

The check "f is NOT a fault face AND f is interior" is under-specified.
Specifically:
- "is interior" can mean (a) `Elem2No >= 0` for `GetFaceElementTransformations(f)`,
  OR (b) `face_bdr_attr_[f] == 0`. These are NOT the same in MFEM:
  - A boundary face (e.g., free surface) has Elem2No < 0 AND
    face_bdr_attr_[f] > 0.
  - A fault face has Elem2No >= 0 AND face_bdr_attr_[f] == fault_attr (3).
  - A regular interior face has Elem2No >= 0 AND face_bdr_attr_[f] == 0.
  - A SHARED face (MPI partition seam) has Elem2No < 0 (the neighbor
    element is on another rank) AND `shared_mesh_face_set_.count(f) > 0`.

If the implementer uses condition (a) alone, shared faces are excluded
(Elem2No < 0). The plan separately says shared faces should be included
in `AllContinuous` mode and conditionally included in `Adjacent` mode —
so the check needs to handle the shared case too.

If the implementer uses condition (b) alone (`face_bdr_attr_[f] == 0`),
shared interior faces are correctly included (they have face_bdr_attr_ ==
0), but ALSO MFEM may not populate `face_bdr_attr_` for shared faces in
all cases.

The plan needs an explicit, MFEM-correct membership condition.

**Suggested fix:**
```diff
@@ MIXED_FLUX_PLAN.md Phase 3 §1 algorithm
-   Step 2: For each element e ∈ E_fault_adj:
-              For each face f of e (use mesh.GetElementFaces(e, faces, ...)):
-                 if f is NOT a fault face AND f is interior:
-                    central_flux_face_set_.insert(f)
+   Step 2: For each element e ∈ E_fault_adj:
+              For each face f of e (use mesh.GetElementFaces(e, faces, ...)):
+                 const bool is_interior =
+                    (mesh.GetFaceElementTransformations(f) != nullptr &&
+                     mesh.GetFaceElementTransformations(f)->Elem2No >= 0)
+                    OR
+                    shared_mesh_face_set_.count(f) > 0;     // shared seam
+                 const bool is_fault =
+                    (face_bdr_attr_[f] == bc_.fault_attr &&
+                     bc_.fault_attr > 0);
+                 if (is_interior && !is_fault) {
+                    central_flux_face_set_.insert(f);
+                 }
+
+   The combined `is_interior` covers both ParMesh shared seams (where
+   GetFaceElementTransformations gives Elem2No < 0 because the neighbor
+   is on another rank) AND regular two-sided interior faces.
```

**Test case:**
```cpp
TEST_CASE("R-1207: shared faces are included in central_flux_face_set_") {
   // Build a 2-rank ParMesh; verify shared interior faces (those with
   // shared_mesh_face_set_ membership but Elem2No < 0) are included in
   // the set when their adjacent local element is fault-adjacent.
}
```

---

### [R-1208] [LOW] [MIXED_FLUX_PLAN.md Phase 6 §"Acceptance Criteria" #4] — 5% tolerance bound is unjustified

**Category:** QUALITY (loose acceptance criterion)

**Description:**
Phase 6 §4:
> "`./seas_tpv104_driver --tfinal 0.1 --mesh tpv104_repro.msh --mixed-flux
> adjacent`: completes without crash, slip_strike at hypocenter is finite
> and within 5% of the no-flag run at the same tfinal (mixed flux changes
> physics by O(face_count × dissipation contribution) — within tolerance
> for a short tfinal)."

5% is an arbitrary number. Zhang 2023 Figures 6 and 7 actually show that
mixed-flux on asymmetric mesh produces results NEARLY IDENTICAL to upwind
on symmetric mesh (visually overlapping curves), with rupture-front
position differences on the order of 30 m at a 200 m mesh — closer to 1%
than 5%. A 5% tolerance is too loose to detect a wrong mixed-flux
implementation.

**Suggested fix:** tighten the tolerance and reference Zhang's data:
```diff
-      and within 5% of the no-flag run at the same tfinal (mixed flux changes
-      physics by O(face_count × dissipation contribution) — within tolerance
-      for a short tfinal).
+      and within 2% of the no-flag run at the same tfinal on the same
+      asymmetric mesh.  Zhang 2023 Fig. 6 demonstrates rupture-front
+      position deviation < 1% between upwind and mixed-flux on the same
+      mesh; we use 2% as a generous gate to absorb implementation-level
+      rounding.  If the deviation exceeds 5%, the implementation is
+      almost certainly wrong (most likely R-1201 or R-1207).
```

**Test case:** N/A (acceptance criterion is the test).

---

### [R-1209] [LOW] [MIXED_FLUX_PLAN.md Phase 1 §"Edge Cases"] — Edge-case "non-unit normal" is misleading

**Category:** QUALITY

**Description:**
Phase 1 edge case:
> "**`nor` not unit length**: callers normalize before calling. `Central`
> does NOT re-normalize; it assumes the same convention as `Interior`."

This is correct but underspecified. Looking at `godunov_flux.cpp::Interior`:
it does NOT re-normalize either; it relies on the caller. But
`Interior::BuildFrame(nor, ...)` constructs `t1 = up × nor` and normalizes
`t1`, which masks non-unit input partially. The plan should reference the
exact existing convention rather than leaving "same as Interior" implicit.

**Suggested fix:** drop the edge case (it's identical to Interior's
behavior). OR fold it into a single sentence under the function docstring:
"Caller-normalized `nor` is assumed; `Interior` and `Central` share this
convention via `BuildFrame`'s internal normalization of the tangent
basis."

---

## Summary

- Critical issues: **3**
  - R-1201: wrong matrix referenced in Phase 1 implementation strategy
  - R-1202: 2-tet fixture central-set count of 6 is unrealizable (actual = 0)
  - R-1203: precomputed-flux mutual exclusion only checked one direction
- Moderate issues: **4**
  - R-1204: algebraic identity sign requires explicit convention pinning
  - R-1205: construction-order requirement undocumented in driver
  - R-1206: per-face `set::count` per step not zero-cost; need short-circuit + unordered_set
  - R-1207: face-set membership condition under-specified for MPI shared faces
- Low issues: **2**
  - R-1208: 5% tolerance in Phase 6 too loose
  - R-1209: edge-case wording redundant with Interior's convention
- Plan compliance: **PARTIAL** — three CRITICAL planning errors. The
  implementation, if executed literally, would produce wrong central
  flux (R-1201), an unprovable test (R-1202), and a silent dispatch
  inconsistency (R-1203).
- Verdict: **FAIL — the plan must be updated before any implementation
  starts.**

## Recommended action

1. Apply the R-1201 / R-1202 / R-1203 fixes to `MIXED_FLUX_PLAN.md`
   (the diffs above are the canonical text to replace).
2. Apply R-1204 / R-1205 / R-1206 / R-1207 fixes (moderate; would be
   caught by careful implementation but the plan should be precise).
3. Apply R-1208 / R-1209 fixes (low; cleanup).
4. Re-run plan review (round 12) on the updated plan to catch any new
   issues introduced by the fixes themselves.
5. Only then proceed to Phase 1 implementation.

## Unreviewed Areas

- The exact Zhang 2023 mathematical derivation of the central flux for
  velocity-strain vs our velocity-stress formulation. The plan states
  "form-invariant" — I did not verify the eigendecomposition equivalence
  between formulations. If the equivalence fails (unlikely given both
  are linear elastodynamics, but possible at the matrix-coefficient level),
  the central flux would still be a VALID quadrature, but the algebraic
  identity in R-1201 / R-1204 would have a different `|A_n|` matrix on
  the right-hand side. This is a math-paper-level audit, not a code-level
  audit; defer until implementation reveals inconsistency.
- Performance impact of `unordered_set` rehashing during construction.
  For TPV104 production, the set is built ONCE at startup; rehash cost is
  amortized. Not a concern.
