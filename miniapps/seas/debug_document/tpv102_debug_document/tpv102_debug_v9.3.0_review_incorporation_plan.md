# Implementation Plan: v9.3.0 Plan Review-Incorporation (R-001 … R-014)

> **Author:** planning agent, 2026-04-21.
> **Target of this plan:** `tpv102_debug_v9.3.0_debug_plan.md`
> (the planning document under review).
> **Input:** `REVIEW.md` (1061 lines) — 14 findings: 4 CRITICAL, 7 MODERATE, 3 LOW.
> **Deliverable:** an edited v9.3.0 plan document that resolves every
> finding, plus a new `##` **"Review-incorporation log (rev 2)"** section
> at the end of the plan recording every change with finding IDs.
> No source-code changes land in this plan — only edits to the plan
> document.

## Overview

This plan edits `tpv102_debug_v9.3.0_debug_plan.md` in four phases
ordered by severity (CRITICALs first so blocking issues are resolved
before MODERATE polish).  Every edit is a textual diff against the
existing plan; every edit cites the finding ID it addresses; every
moderate/critical edit is accompanied by either a plan-level
acceptance gate ("after editing, the plan now says X") or a new
unit-test sketch inserted inline in the plan for a future
implementation agent to execute.

The plan explicitly does NOT write any new C++ code.  Unit-test
C++ snippets added here are plan-level descriptions (they go into the
Phase acceptance criteria of the v9.3.0 plan) and will be implemented
later by the `/code-implement` agent when v9.3.0 implementation begins.

As of 2026-04-22, this review-incorporation plan follows a superseding
product decision for I-06: TPV102 migrates to total-stress Q as the
only supported runtime representation.  Any earlier review text that
would preserve a dual `fluctuation|total` runtime mode or add a
`--state-representation` CLI must be rewritten into a one-way migration
plan instead.

## Constraints

### Document-integrity constraints
- **Additive-preferring:** where possible, replace or insert — do NOT
  delete entire blocks of the plan that carry history-value context.
  The existing "Overview", "Constraints", "Risk Assessment",
  "Deferred / explicitly NOT in scope", and "Quick-reference for
  implementation agent" sections remain intact.
- **No renumbering:** Phases 1–7 in the existing plan keep their
  numbers.  New Phase acceptance criteria are APPENDED; existing
  acceptance criteria are AMENDED only where a finding requires it.
- **Review-incorporation log is mandatory:** every CRITICAL and
  MODERATE finding gets a 1-line entry in the new §"Review-
  incorporation log (rev 2)" section at the end of the document.
- **CLI coordination with the existing driver:** the final v9.3.0 plan
  must NOT add a `--state-representation` or `--no-nucleation` CLI for
  TPV102.  If Phase 4 needs a locked-fault harness before Phase 5 lands,
  that suppression belongs in the test fixture / setup path, not in a
  shipped runtime flag.

### Convention constraints
- Use the review's exact textual-patch syntax (```diff ... ```) where
  a verbatim replacement is prescribed — do not paraphrase.
- Where the review says "Suggested fix:" and provides a diff, use
  that diff.
- Where the review says "Test case:" and provides a C++ snippet,
  append that snippet to the relevant phase's "Acceptance Criteria"
  subsection under a new "**Proposed unit tests**" bullet.

### Scope constraints
- R-011 changes Phase 7's scope (deferring `EvaluateADERTotal` and the
  ADER-vs-RK4 total-stress acceptance test to v9.3.1).  The
  plan must reflect this scope reduction in the §"Scope and timeline"
  table (reduce Phase 7's LOC estimate and wall-clock) AND in the
  §"Deferred / explicitly NOT in scope" section (add
  `EvaluateADERTotal` explicitly as a v9.3.1 prerequisite).
- R-014 is superseded by the total-only migration decision.  The edited
  plan should remove the need for a TPV102 `StateRepresentation` enum
  altogether rather than relocate it.

## Phase 1: Address CRITICAL findings R-001..R-004

### Goal
After this phase, the v9.3.0 plan no longer contains any
specification-level bug that blocks implementation.  All 4
CRITICAL findings are resolved by direct textual patches.

### Files to Modify
- `miniapps/seas/debug_document/tpv102_debug_document/tpv102_debug_v9.3.0_debug_plan.md`
  — apply 4 textual patches (one per CRITICAL finding) using the
  `Edit` tool.  No other file is touched in this phase.

### Detailed Requirements

**1. R-001 (self-contradictory pre-stress sign).**  Replace the
Phase 4 "Files to Create" bullet that currently lists
`sigma_yy = -sigma_n0` with the corrected sign-consistent version.
Exact patch to apply:

```diff
  `dynamic/tpv102_setup_total.hpp` — companion to `tpv102_setup.hpp`
  exposing `InitializeStateTotal(Q, ndof_total, pre_stress_at_dof)`
  which fills every DOF's Q[SXX/SYY/SZZ/SXY/…] with the rotated
  pre-stress tensor in GLOBAL coordinates.  For TPV102 (homogeneous
- half-space under a vertical fault) the pre-stress is:
- ```
- sigma_xx = 0,   sigma_yy = -sigma_n0,   sigma_zz = 0,
- sigma_xy = -tau_ini,   sigma_yz = 0,   sigma_xz = 0
- ```
- — i.e. compressive normal stress of 120 MPa on the y=const plane
- plus shear of 75 MPa on the x-y plane.
+ half-space under a vertical fault), under MFEM's compression-
+ positive convention (miniapps/seas/CLAUDE.md §"Normal stress"), the
+ bulk pre-stress written to Q is:
+ ```
+ Q[SXX] = 0,          Q[SYY] = +sigma_n0,   Q[SZZ] = 0,
+ Q[SXY] = -tau_ini,   Q[SYZ] = 0,           Q[SXZ] = 0
+ ```
+ Rotation check: with BP5 canonical `n=(0,-1,0)`, `t1_dip=(0,0,-1)`,
+ `t2_strike=(+1,0,0)`, the outer-product rotation gives
+ `sigma_yy_global = +sigma_nn_local = +sigma_n0`  (compressive, as
+ expected) and `sigma_xy_global = -sigma_nt2_local = -tau_ini`
+ (strike shear, sign flip from canonical-local).
```

Then ADD a "**Proposed unit tests**" bullet to Phase 4's Acceptance
Criteria containing the `test_R001_initstate_total_sign.cpp` snippet
verbatim from REVIEW.md lines 104-127.

**2. R-002 (nucleation sign flip).**  Rewrite Phase 5's "Projection
step" comment block to explicitly include the sign-flip when mapping
fault-local `tau_nt2` to global `sigma_xy`.  Exact patch:

```diff
 // In the canonical fault-local frame, nucleation adds to tau_nt2
-// (strike shear).  In global frame this maps to sigma_xy (see I-03
-// note).
+// (strike shear).  In global frame this maps to -sigma_xy (I-03
+// note: sigma_xy_global = -tau_nt2_local under BP5 canonical
+// t2=+x, n=-y).  THE SIGN MATTERS — adding dtau to tau_nt2_local
+// requires *subtracting* dtau from Q[SXY, QP].
 //
 // Because Q is stored per ELEMENT (not per QP), and the QP lies
 // strictly INSIDE the element, we use the local L2 Lagrange basis:
-//   sigma_xy(ip) = Σ_j shape_j(ip) * Q[SXY, elem_dof(j)]
-// Adding dtau at the QP requires adding dtau * M_ref^-1 * shape to
-// the element DOFs, where M_ref is the element face-to-element
-// coupling mass.
+//   sigma_xy(ip) = Σ_j shape_j(ip) * Q[SXY, elem_dof(j)]
+// Injecting tau_nt2_local += dtau at the QP therefore requires
+//   Q[SXY, elem_dof(j)] -= dtau   (for interpolatory GaussLobatto
+//                                   nodal basis with j = node at QP)
+// or equivalently adding -dtau * M_ref^-1 * shape for the L2
+// projection path.
```

Append REVIEW.md's `test_R002_nucleation_sign_match.cpp` (lines
196-210) to Phase 5's Acceptance Criteria under "**Proposed unit
tests**".

**3. R-003 (missing psi-pristineness + DIAG guard).**  Insert the
`#ifndef NDEBUG` psi-invariant guard and the `#ifdef
SEAS_DIAG_FAULT_FLUX` C-1 diagnostic block into the Phase 3
`EvaluateTotal` pseudocode.  Exact patches:

(a) At the top of the `EvaluateTotal` body (just after the opening
brace, before the `homog_ok` lambda):
```diff
 void FaultFaceFlux::EvaluateTotal(DOFData &data, ...) const
 {
+#ifndef NDEBUG
+   const real_t psi_at_entry = data.psi;   // R-V92-H07 invariant
+#endif
    auto homog_ok = [](real_t a, real_t b) { ... };
    MFEM_VERIFY(homog_ok(data.Zp_plus, data.Zp_minus) && ...);
```

(b) Right after the `ComputeTrialTraction` call:
```diff
    real_t sigma_n_trial, tau1_trial, tau2_trial;
    ComputeTrialTraction(data, Q_plus, Q_minus,
                         sigma_n_trial, tau1_trial, tau2_trial);
+
+#ifdef SEAS_DIAG_FAULT_FLUX
+   if (data.diag_print)
+   {
+      std::fprintf(stderr,
+         "[C-1 EVAL-TOTAL] rank=%d  tau1_trial=%+.3e Pa  tau2_trial=%+.3e Pa  "
+         "psi=%.3e\n",
+         g_seas_my_rank, tau1_trial, tau2_trial, data.psi);
+   }
+#endif
```

(c) Just before the closing brace (after Step 5's DOFData updates):
```diff
    data.sigma_n_corr = sigma_n_corr;
+
+#ifndef NDEBUG
+   MFEM_ASSERT(data.psi == psi_at_entry,
+               "FaultFaceFlux::EvaluateTotal mutated data.psi "
+               "(before = " << psi_at_entry << ", after = " << data.psi
+               << ").  R-V92-H07: the driver's coupled-RK4-on-psi "
+               "integrator assumes this function is psi-pure.");
+#endif
 }
```

Append REVIEW.md's `test_R003_evaluatetotal_psi_invariant.cpp`
(lines 293-301) to Phase 3's Acceptance Criteria under "**Proposed
unit tests**".

**4. R-004 (dead cached-member spec).**  Remove the
`free_surface_godunov_state_` cached member and the `R_11`/`R_21`
compliance test criterion that references it.  Exact patch:

```diff
 ### Files to Modify
-- `dynamic/godunov_flux.hpp` — add `FreeSurfaceGodunov` declaration + a
-  cached `free_surface_godunov_state_` 9×9 matrix member.
-- `dynamic/godunov_flux.cpp` — add definition (~80 LOC) using the
-  eigenvector matrix `R` already built in the ctor; precompute the
-  Godunov-state projector at ctor time.
+- `dynamic/godunov_flux.hpp` — add `FreeSurfaceGodunov` declaration only.
+- `dynamic/godunov_flux.cpp` — add definition (~50 LOC) using `Zp_`
+  and `Zs_` already stored by the ctor.  No new data members: the
+  compliance block diag(-1/Zp, -1/Zs, -1/Zs) appears inline as
+  (invZp, invZs, invZs) multiplied into the velocity-traction
+  update — algebraically equivalent to a cached projector without
+  the dead-member spec-vs-implementation mismatch.
```

And in Phase 1 Acceptance Criteria:

```diff
- [ ] The ctor-built `R_11` and `R_21` compliance block match
- hand-derived `diag(-1/Zp, -1/Zs, -1/Zs)` bit-exactly (tested via
- `EXPECT_DOUBLE_EQ`).
+ [ ] A direct check
+     `FreeSurfaceGodunov(nor=(0,0,1), Q_self with SXX=1, rest=0)`
+     equals `FreeSurface(same inputs)` within 1 ULP per component.
+     (Inline invZp/invZs form is algebraically diag-compliance; no
+     separate member to inspect.)
```

Append REVIEW.md's `test_R004_free_surface_godunov_no_cached_member.cpp`
(lines 380-388) as a comment-only assertion in Phase 1's Acceptance
Criteria ("regression test: the Phase 1 patch must not add member
variables to `GodunovFlux`").

### Interfaces
No code interfaces changed in this phase — only plan-document text.

### Edge Cases to Handle
- If the v9.3.0 plan has been edited between the REVIEW.md run and
  this patch (e.g., someone added a paragraph to Phase 4), the
  diff-context search may fail.  Mitigation: use `Read` to re-verify
  the exact block text before each `Edit` call; if the context has
  drifted, abort that specific patch and note it in the review-
  incorporation log.

### Acceptance Criteria
- [ ] The 4 patches above are applied verbatim via `Edit` to the
      v9.3.0 plan.
- [ ] A `grep -n "sigma_yy = -sigma_n0" tpv102_debug_v9.3.0_debug_plan.md`
      returns NO matches after the patch (R-001 gone).
- [ ] A `grep -n "free_surface_godunov_state_" tpv102_debug_v9.3.0_debug_plan.md`
      returns NO matches after the patch (R-004 gone).
- [ ] A `grep -n "psi_at_entry\|MFEM_ASSERT.*psi_at_entry" tpv102_debug_v9.3.0_debug_plan.md`
      returns at least two matches (entry + exit, R-003).
- [ ] A `grep -n "subtracting dtau\|Q\[SXY, elem_dof(j)\] -= dtau" tpv102_debug_v9.3.0_debug_plan.md`
      returns at least one match (R-002 sign flip explicit).
- [ ] Four new "**Proposed unit tests**" bullets appear in Phases 1,
      3, 4, 5 acceptance criteria.
- [ ] Pandoc rebuild of the plan PDF (if the plan is published as PDF)
      reflects the edits.

### Dependencies
- Depends on: review of REVIEW.md (complete).
- Required by: Phase 2.

---

## Phase 2: Address MODERATE findings R-005..R-011

### Goal
After this phase, all 7 moderate-severity findings are folded into
the v9.3.0 plan via targeted textual patches.  Phase 7 is re-scoped
per R-011 (ADER + total coupling deferred to v9.3.1).

### Files to Modify
- `tpv102_debug_v9.3.0_debug_plan.md` — 7 patches as detailed below.
- `drivers/tpv102_driver.cpp` is NOT modified here; a NEW
  Phase 4 sub-bullet in the v9.3.0 plan specifies that TPV102 is
  migrated to total stress without adding a new dual-mode / fallback
  CLI.

### Detailed Requirements

**5. R-005 (CLI flag silent-fallback).**  Replace Phase 2's driver
pseudocode with the case-insensitive + warning variant:

```diff
-std::string fs_bc_str = GetStringArg(argc, argv,
-                                     "--free-surface-bc", "gamma");
-FreeSurfaceBCMode fs_bc_mode = (fs_bc_str == "godunov")
-                                 ? FreeSurfaceBCMode::Godunov
-                                 : FreeSurfaceBCMode::Gamma;
-wave.SetFreeSurfaceBCMode(fs_bc_mode);
-if (rank == 0) {
-   std::cout << "Free-surface BC: " << fs_bc_str << "\n";
-}
+std::string fs_bc_str = GetStringArg(argc, argv,
+                                     "--free-surface-bc", "gamma");
+std::string fs_bc_lower = fs_bc_str;
+std::transform(fs_bc_lower.begin(), fs_bc_lower.end(),
+               fs_bc_lower.begin(),
+               [](unsigned char c){ return std::tolower(c); });
+FreeSurfaceBCMode fs_bc_mode;
+if (fs_bc_lower == "godunov") {
+   fs_bc_mode = FreeSurfaceBCMode::Godunov;
+} else if (fs_bc_lower == "gamma") {
+   fs_bc_mode = FreeSurfaceBCMode::Gamma;
+} else {
+   if (rank == 0) {
+      std::cerr << "WARNING: unknown --free-surface-bc=\""
+                << fs_bc_str << "\"; falling back to gamma.\n";
+   }
+   fs_bc_mode = FreeSurfaceBCMode::Gamma;
+   fs_bc_str  = "gamma";
+}
+wave.SetFreeSurfaceBCMode(fs_bc_mode);
+if (rank == 0) {
+   std::cout << "Free-surface BC: " << fs_bc_str
+             << " (effective: "
+             << (fs_bc_mode == FreeSurfaceBCMode::Godunov
+                  ? "godunov" : "gamma") << ")\n";
+}
```

Append REVIEW.md's `test_R005_free_surface_bc_flag_parse.cpp`
(lines 470-478) under Phase 2's Acceptance Criteria as a
"**Proposed unit test**".

**6. R-006 (Phase 4 staging without a shipped fallback CLI).**  Replace
the old dual-mode driver-plumbing text with a total-only migration
statement.  Exact patch:

```diff
 Vector Q;
-if (state_rep == StateRepresentation::Total) { ... }
+InitializeStateTotal(Q, ndof_total,
+                     TPV102Params::sigma_n, TPV102Params::tau_ini);
+// Phase 4 may land before Phase 5's nucleation-to-Q injection.  In
+// that intermediate state, the acceptance fixture suppresses
+// nucleation inside the harness (e.g. nuc_dtau = 0) rather than
+// adding a shipped `--no-nucleation` runtime flag.
```

Additionally, replace the old Phase 4 "Files to Modify" bullet:

```diff
- `drivers/tpv102_driver.cpp` — parse CLI flag; add total-mode
-  initializer; branch the `InitializeState(Q)` call.
+ `drivers/tpv102_driver.cpp` — remove the dual-mode state-
+  representation concept from the Phase 4 text.  Initialize TPV102 in
+  total stress unconditionally; if a locked fixture is needed before
+  Phase 5, describe it as a test-harness/setup condition, not a
+  shipped CLI flag.
```

Append a Phase 4 proposed unit test that exercises the locked-fault
harness under total-stress initialization (the old
`test_R006_total_mode_nucleation_abort.cpp` no longer applies once the
runtime fallback CLI is removed).

**7. R-007 (prose error: "pre-stresses cancel").**  Replace the Phase
3 Math block's erroneous "cancel" sentence:

```diff
-Note the sum term pre-stresses cancel across $\pm$ sides identically
-(symmetric fault), so $\sigma_{n}^{*,\text{tot}} =
-\sigma_{n}^{*,\text{fluc}} + 2\eta_{p}\sigma_{n,0}/Z_{p} =
-\sigma_{n}^{*,\text{fluc}} + \sigma_{n,0}$ (using $\eta_{p} = Z_{p}/2$).
+On a symmetric fault, pre-stress on both sides has the same value
+$\sigma_{n,0}$, so the sum term **adds** to $2\sigma_{n,0}/Z_p$ (it
+does NOT cancel).  The velocity-difference term is pre-stress-free,
+and with $\eta_p = Z_p/2$ (harmonic mean for homogeneous):
+$\sigma_{n}^{*,\text{tot}} =
+\sigma_{n}^{*,\text{fluc}} + \eta_p \cdot 2\sigma_{n,0}/Z_p =
+\sigma_{n}^{*,\text{fluc}} + \sigma_{n,0}$.
```

No new test case — covered implicitly by R-008 / R-009 (see R-007
rationale in REVIEW.md line 599).

**8. R-008 (too-loose tilt tolerance).**  Replace Phase 1's
Acceptance Criterion #3 (the 100 ULP allowance) with a uniform 10
ULP requirement at all tilt angles:

```diff
 - [ ] New `seas_test_free_surface_godunov`:
-  - On `nor = (0, 0, 1)` (axis-aligned z=0), for each of 50 random
-    `Q_self` states, `F_h_gamma` and `F_h_godunov` agree to 10 ULP.
-  - On a 10°-tilted normal, they agree to 50 ULP.
-  - On a 45°-tilted + non-axis-aligned `(t1, t2)` choice from
-    `BuildFrame`, they still agree within 100 ULP (γ-mirror's
-    frame-invariance on a flat surface is exact; Godunov is exact by
-    construction).
+  - On `nor = (0, 0, 1)`, 10 ULP per component for 50 random Q_self.
+  - On `nor` rotated by 10°, 30°, 45°, 75° from vertical: 10 ULP.
+    The equivalence is algebraically exact (let
+    Δ := Q_ghost − Q_god; then A⁻·Δ = 0 by construction of Q_god,
+    so the resulting flux is identical to the last FP bit); any
+    tilt-dependent divergence indicates a bug in BuildFrame's
+    orthonormalisation or in the Q_god projection.
```

Append REVIEW.md's `test_R008_fs_godunov_tilt_invariance.cpp`
(lines 663-677) to Phase 1 Acceptance Criteria as a "**Proposed unit
test**" — replacing the vague multi-tolerance line in the original
plan.

**9. R-009 (nucleation ±-side symmetry + nodal vs QP).**  Replace
Phase 5's "Detailed Requirements" §1 "Fault-QP-to-global-DOF map"
with the ±-side-explicit version from REVIEW.md lines 729-761:

```diff
-**1. Fault-QP-to-global-DOF map.**  Build once at init:
-```cpp
-std::vector<int> fault_qp_to_elem_;   // element index for each fault QP
-std::vector<int> fault_qp_to_elem_dof_;// local DOF index in that element
-```
+**1. Fault-QP-to-elements-DOF map.**  Build once at init.  At each
+fault QP there are TWO elements (+ and − sides), and BOTH need the
+nucleation injection (the trial traction in `EvaluateTotal` averages
+`Q+[SXY]/Zs + Q−[SXY]/Zs`, so a one-sided injection would halve the
+perturbation).  Store:
+```cpp
+std::vector<int> fault_qp_to_elem_plus_;       // + side elem id
+std::vector<int> fault_qp_to_elem_minus_;      // − side elem id
+std::vector<int> fault_qp_to_elem_dof_plus_;
+std::vector<int> fault_qp_to_elem_dof_minus_;
+```
+Partition-seam fault faces (shared between MPI ranks) have only one
+local side; the other side's injection is performed on the peer rank
+via the same ghost-exchange mechanism WaveOperator already uses
+(`dynamic/wave_operator.inl::ComputeSharedFaceFluxRHS`).
```

And extend the "GaussLobatto" requirement with a concrete startup
check:

```diff
-For the Phase 5 minimum, require **GaussLobatto basis + interpolatory
-nodal injection** and MFEM_VERIFY this at setup.
+For the Phase 5 minimum, require **GaussLobatto basis + interpolatory
+nodal injection** and MFEM_VERIFY this at setup.  Concretely, verify
+(at driver init, rank 0 aborts propagate via `MFEM_ABORT`):
+```cpp
+MFEM_VERIFY(L2_FECollection::GetBasisType() == BasisType::GaussLobatto,
+            "Phase 5 total-mode nucleation requires GaussLobatto "
+            "basis; got basis type " << static_cast<int>(
+                L2_FECollection::GetBasisType()));
+```
+For each fault QP, additionally verify that the chosen
+`elem_dof_{plus,minus}` has basis-function value ≥ 1 − 1e-12 at the
+QP's reference coordinate (interpolatory check).  Abort at setup if
+any QP fails; the L2-projection fallback is deferred to v9.3.1.
```

Append REVIEW.md's `test_R009_nucleation_total_both_sides.cpp`
(lines 765-773) to Phase 5 Acceptance Criteria under "**Proposed unit
tests**".

**10. R-010 (test semantics unclear — corner fixture).**  Replace
Phase 6's `test_tpv102_free_surface_bc_variants` description:

```diff
-`tests/unit/test_tpv102_free_surface_bc_variants.cpp` — run the same
-rig with γ vs godunov free surface, assert agreement to 10 ULP on a
-flat surface and ≤ 10 Pa on the corner fixture (tolerant of the
-actual corner-pumping divergence we HOPE to see eliminated).
+`tests/unit/test_tpv102_free_surface_bc_variants.cpp` — run the rig
+with godunov free surface.  Two distinct sub-tests:
+ 1. **Equivalence gate** (flat axis-aligned surface): assert γ and
+    godunov agree to 10 ULP per component, per face, for the whole
+    run.  If this fails, the Phase 1 Godunov projection is buggy
+    (independent of corner physics).
+ 2. **Corner-elimination test** (fault-surface corner fixture):
+    assert `|σ_yy at the fault-corner QP under godunov BC|` ≤ 10 Pa
+    at t = 0.5 s.  This is the H-V92-G1 elimination test — it does
+    NOT compare against γ (γ-mirror's corner pump is the bug we
+    want to eliminate; comparing Godunov to γ would pass trivially
+    if Godunov also has a small corner pump of its own).
+ 3. **Comparison datum** (informational, not an assert): print
+    γ-mirror's corner σ_yy at t = 0.5 s to the test log so future
+    debugging can see what the original bug magnitude was.
```

**11. R-011 (ADER + total workaround is under-specified).**  Replace
Phase 7's "Medium risk" bullet with an explicit scope deferral:

```diff
-Under ADER, F_h is computed from the time-integrated state I, not
-Q.  The fault flux's imposed-state construction expects input in
-the same representation as Q.  Under total + ADER, everything stays
-total — the `EvaluateADERTotal` variant (not in v9.3.0 scope; add
-to ADER plan Phase 5 as a follow-up) is a straightforward
-copy of Phase 3 with `Evaluate` → `EvaluateTotal`.  Until that is
-landed, the v9.3.0 Phase 7 test uses an explicit workaround:
-rescale I by 1/Δt inside `ComputeADERFaceFluxRHS` (ADER plan Phase 6)
-before calling `EvaluateTotal`, then rescale the output by Δt.
+Under ADER, F_h is computed from the time-integrated state I, not
+Q.  The fault flux's imposed-state construction expects input in
+the same representation as Q.  Under total + ADER, everything stays
+total, but `EvaluateADERTotal` (a copy of Phase 3 with `Evaluate`
+→ `EvaluateTotal` and I instead of Q) is a hard prerequisite.
+It is NOT in v9.3.0 scope.
+
+**Consequence:** Phase 7's `--time-integrator=ader` acceptance test on
+the migrated total-stress TPV102 path IS BLOCKED on the ADER plan's
+Phase 5 landing `EvaluateADERTotal`.  Phase 7 in v9.3.0 adds ONLY the
+cross-reference note to v9.3.0 Phase 4 in the ADER plan.  The
+`test_ader_total_vs_rk4_total` acceptance test is DEFERRED to v9.3.1.
```

Additionally, update the "Scope and timeline" table in the v9.3.0
plan:

```diff
-| 7 ADER wiring        |  60 |  40 | 0.5 day |
-| **Total**            | **~750** | **~320** | **~6.5 days** |
+| 7 ADER wiring (banner only) |  20 |  15 | 0.25 day |
+| **Total**                   | **~710** | **~295** | **~6.25 days** |
```

And the "Deferred / explicitly NOT in scope" section:

```diff
 - **`EvaluateADERTotal`.**  Properly a follow-up addition to the ADER
   plan Phase 5, not v9.3.0.  v9.3.0 Phase 7 uses an explicit rescale
   workaround.
+- **`test_ader_total_vs_rk4_total`.**  Moved from v9.3.0 Phase 7 to
+  v9.3.1, contingent on `EvaluateADERTotal`.
```

### Interfaces
No runtime-code interfaces changed here — still plan-document edits.

### Edge Cases to Handle
- The locked-fault harness in Phase 4 must be described so readers do
  not mistake it for a shipped TPV102 runtime mode.  Use wording like
  "nucleation suppressed inside the acceptance fixture" rather than a
  new user-facing CLI.

### Acceptance Criteria
- [ ] 7 patches applied to the v9.3.0 plan, one per MODERATE finding.
- [ ] `grep -n "pre-stresses cancel"` returns NO matches (R-007 gone).
- [ ] `grep -n "100 ULP"` in Phase 1 context returns NO matches (R-008).
- [ ] `grep -n "fault_qp_to_elem_plus_"` returns at least 1 match (R-009).
- [ ] `grep -n "state-representation" tpv102_debug_v9.3.0_debug_plan.md`
      returns 0 matches in the TPV102 driver-plumbing / CLI context
      after the edits (R-006 / total-only migration).
- [ ] Scope and timeline table shows `~710 / ~295 / ~6.25 days` (R-011).
- [ ] Deferred section lists `EvaluateADERTotal` AND
      `test_ader_total_vs_rk4_total` (R-011).
- [ ] 7 "**Proposed unit tests**" bullets exist across Phases 1, 2,
      4, 5 acceptance criteria (R-005, R-006, R-008, R-009;
      R-007/R-010/R-011 have no new tests).

### Dependencies
- Depends on: Phase 1 of this plan.
- Required by: Phase 3.

---

## Phase 3: Address LOW findings R-012..R-014

### Goal
After this phase, the three low-severity findings are addressed via
small clarifying edits.  The plan's prose becomes more internally
consistent and removes leftover dual-mode / enum assumptions that no
longer apply under the total-only migration decision.

### Files to Modify
- `tpv102_debug_v9.3.0_debug_plan.md` — 3 patches.

### Detailed Requirements

**12. R-012 (clarifying comment on zeroed DOFData pre-stress).**
Insert an explanatory comment at the zeroing loop in Phase 4 driver
pseudocode:

```diff
    for (int i = 0; i < num_fault_total; i++) {
       dof_data[i].sigma_n0 = 0.0;
       dof_data[i].tau1_0   = 0.0;
       dof_data[i].tau2_0   = 0.0;
+      // NOTE: sigma_n_corr / tau_i_corr are OUTPUT fields (set by
+      // InitializeFaultDOFs to the physical pre-stress, and
+      // overwritten on first EvaluateTotal call).  We do NOT zero
+      // them here so the initial ParaView snapshot at t=0 reflects
+      // the physical pre-stress in the migrated total-stress path.
+      // Any refactor that reads
+      // dof_data[i].{sigma_n0,tau*_0} later must confirm the field
+      // has been re-initialised for total stress.
    }
```

**13. R-013 (DOFData "frozen" constraint overstated).**  Relax the
Constraints-section promise from "append-only" to a more honest
statement acknowledging the conditional DIAG tail:

```diff
-- **`DOFData` struct layout** — append-only.  New fields
-  (`Q_total_init_*`) at end of struct; `alignof(DOFData)` must not change.
+- **`DOFData` struct layout** — append-only in the non-DIAG region.
+  `sizeof(DOFData)` already varies by build configuration because
+  the `#ifdef SEAS_DIAG_FAULT_FLUX` tail adds a `bool diag_print`.
+  Any v9.3.0 new field MUST go BEFORE the `#ifdef SEAS_DIAG_FAULT_FLUX`
+  block so release builds (no DIAG) retain the non-DIAG layout
+  unchanged.  `alignof(DOFData)` must not change.
```

**14. R-014 (enum placement) is superseded by total-only migration.**
Remove the need for a TPV102 `StateRepresentation` enum entirely.
Apply to Phase 4 "Files to Modify":

```diff
- `dynamic/seas_dynamic_operator.hpp` — add `enum class StateRepresentation`.
- `dynamic/wave_operator.hpp` — add `state_rep_` member + setter.
+- remove the `StateRepresentation` design entirely from the TPV102
+  plan text.  Phase 4 should describe unconditional total-stress
+  initialization plus unconditional `EvaluateTotal` dispatch for the
+  TPV102 path.
```

Also verify that no other reference in the v9.3.0 plan points readers
to `StateRepresentation`, `seas_dynamic_operator.hpp`, or
`wave.SetStateRepresentation` in the TPV102 migration context.

### Edge Cases to Handle
- If either Phase 2 or Phase 4 still references `StateRepresentation`
  in the existing plan text, the Phase 3 edit must remove those
  references too.

### Acceptance Criteria
- [ ] 3 low-priority patches applied.
- [ ] `grep -n "append-only" tpv102_debug_v9.3.0_debug_plan.md`
      matches BOTH the old promise location AND the replacement text
      (i.e., the old bullet has been replaced, not just commented
      out).
- [ ] `grep -n "StateRepresentation\|state-representation"`
      in the plan returns 0 matches in the TPV102 migration context.

### Dependencies
- Depends on: Phase 2.
- Required by: Phase 4.

---

## Phase 4: Add §"Review-incorporation log (rev 2)"

### Goal
After this phase, the v9.3.0 plan ends with a new top-level section
that lists every finding addressed, its resolution, and the line(s)
in the plan that changed.  This makes the plan self-documenting for
future reviewers.

### Files to Modify
- `tpv102_debug_v9.3.0_debug_plan.md` — append a new section at the
  end (after the existing "Quick-reference for implementation agent"
  section, before EOF).

### Detailed Requirements

**Insert at the end of the document** (the exact header + body
follows; the implementing agent should write this verbatim):

```markdown
---

## Review-incorporation log (rev 2, 2026-04-21)

The findings below are consumed from `REVIEW.md` (14 items; 4
CRITICAL, 7 MODERATE, 3 LOW).  Each row cites the finding ID, the
v9.3.0 plan section where the fix landed, and a short note.  The
implementing agent SHOULD NOT need to re-read REVIEW.md to execute
v9.3.0 — every load-bearing change from the review is now in the
plan itself.

| ID | Sev | Landed in | Resolution |
|---|---|---|---|
| R-001 | CRITICAL | Phase 4 "Files to Create" | Pre-stress signs corrected: `Q[SYY] = +sigma_n0`, `Q[SXY] = -tau_ini`.  Rotation check added. |
| R-002 | CRITICAL | Phase 5 "Projection step" | Sign of bulk-Q nucleation injection made explicit (`Q[SXY, node] -= dtau`). |
| R-003 | CRITICAL | Phase 3 `EvaluateTotal` body | Added `#ifndef NDEBUG` psi-invariant entry/exit + DIAG print block. |
| R-004 | CRITICAL | Phase 1 "Files to Modify" + AC | Removed dead `free_surface_godunov_state_` cached member; replaced `R_11`/`R_21` AC with an `invZp`/`invZs` inline-compliance equivalence test. |
| R-005 | MODERATE | Phase 2 driver pseudocode | Case-insensitive `--free-surface-bc` parsing + warning on unknown value; banner reports both requested + effective mode. |
| R-006 | MODERATE | Phase 4 driver pseudocode + Phase 4 files-to-modify bullet | Reframed Phase 4 as total-only migration with a locked-fault harness in tests, not a shipped `--no-nucleation` / `--state-representation` fallback. |
| R-007 | MODERATE | Phase 3 Math block | Corrected "pre-stresses cancel" → "pre-stresses add" in the trial-traction derivation. |
| R-008 | MODERATE | Phase 1 AC #3 | Uniform 10 ULP tolerance at all tilt angles; rationale (A⁻·Δ = 0) added. |
| R-009 | MODERATE | Phase 5 Detailed Requirements §1 | Per-fault-QP injection must touch BOTH + and − element DOFs.  GaussLobatto + interpolatory-nodal check required at setup. |
| R-010 | MODERATE | Phase 6 test_tpv102_free_surface_bc_variants spec | Split into equivalence gate (flat, vs γ), elimination test (corner, vs 0), informational γ readout. |
| R-011 | MODERATE | Phase 7 "Medium risk" + Scope table + Deferred section | ADER + total test deferred to v9.3.1; Phase 7 scope reduced to cross-reference only under the total-only TPV102 path. |
| R-012 | LOW | Phase 4 driver pseudocode | Added maintainer comment explaining why sigma_n_corr / tau_i_corr are NOT zeroed in the migrated total-stress path. |
| R-013 | LOW | Constraints "Interface constraints" | Relaxed "append-only" language to acknowledge the conditional DIAG tail. |
| R-014 | LOW | Phase 4 "Files to Modify" | Removed the need for a TPV102 `StateRepresentation` enum under the total-only migration decision. |

### Patches that add new plan-level unit tests
| Phase | Test file | Finding ID |
|---|---|---|
| 1 | `tests/unit/test_R004_free_surface_godunov_no_cached_member.cpp` | R-004 |
| 1 | `tests/unit/test_R008_fs_godunov_tilt_invariance.cpp` | R-008 |
| 2 | `tests/unit/test_R005_free_surface_bc_flag_parse.cpp` | R-005 |
| 3 | `tests/unit/test_R003_evaluatetotal_psi_invariant.cpp` | R-003 |
| 4 | `tests/unit/test_R001_initstate_total_sign.cpp` | R-001 |
| 4 | `tests/unit/test_R006_phase4_locked_fault_harness.cpp` | R-006 |
| 5 | `tests/unit/test_R002_nucleation_sign_match.cpp` | R-002 |
| 5 | `tests/unit/test_R009_nucleation_total_both_sides.cpp` | R-009 |

All 8 tests are plan-level stubs.  The `/code-implement` agent will
flesh them out during v9.3.0 implementation using the TEST_ASSERT /
TEST_NEAR pattern established in
`tests/unit/test_fault_surface_vtu_k4.cpp` and
`tests/unit/test_fault_basis_dip_strike_symmetry.cpp`.

### Open items from REVIEW.md "Unreviewed Areas"

The reviewer flagged four unreviewed areas (REVIEW.md lines
1043-1060).  These remain open for v9.3.0 but are tracked here:

1. **ADER plan not reviewed** — a separate REVIEW.md round on
   `tpv102_ader_time_integration_plan.md` is the natural follow-up.
   Out of scope for v9.3.0 unless ADER lands before v9.3.0 finishes.
2. **Codebase walkthrough I-section not cross-checked** — R-001's
   rotation analysis does cross-check I-03 by direct computation;
   the other I-items rely on the walkthrough's citations.
3. **Phase 6 tolerance not verified against numerical noise floor** —
   the 1 Pa / 1 nm / 1 nm/s tolerances may need tightening or
   loosening based on a pilot run.  Add a Phase 6 sub-step:
   pilot-run the 1-element rig, measure noise floor, adjust the
   acceptance tolerance.
4. **`ComputeSharedFaceFluxRHS` dispatch parity** — Phase 4's
   instruction "Same dispatch in `ComputeSharedFaceFluxRHS`"
   requires line-by-line verification against
   `dynamic/wave_operator.inl:1318` that the shared branch exposes
   the same Evaluate/EvaluateTotal switch point.  Implement-time
   check (not a plan edit).
```

### Edge Cases to Handle
- If any of the 14 patches from Phases 1–3 FAILED to apply (e.g.,
  plan text drifted), the review-incorporation log MUST note the
  failure rather than claiming "Landed in …".  Use a row marker
  "⚠ FAILED — [reason]" and escalate to the user.

### Acceptance Criteria
- [ ] New section `## Review-incorporation log (rev 2, 2026-04-21)`
      exists at the end of the plan document.
- [ ] Table lists all 14 findings with Landed-in references.
- [ ] Table of proposed unit tests lists 8 tests (matching the count
      of R-001..R-009 that have tests; R-007, R-010, R-011, R-012,
      R-013, R-014 have no new tests).
- [ ] Section "Open items from REVIEW.md 'Unreviewed Areas'" lists
      all 4 items from REVIEW.md lines 1043-1060.

### Dependencies
- Depends on: Phases 1, 2, 3.
- Required by: nothing (final phase).

---

## Testing Strategy

This plan's deliverable is a document edit, not code — so "testing"
means verifying that the edits landed correctly.  Per-phase acceptance
criteria above specify grep-based verification: specific strings
must be present or absent in the patched plan.

Additionally, the plan-level C++ unit tests (inserted into the
v9.3.0 plan by this round) are themselves test sketches, to be
fleshed out during v9.3.0 implementation.  They are NOT run here —
they appear in the plan as "**Proposed unit tests**" bullets under
each Phase's Acceptance Criteria.

## Risk Assessment

### High risk
- **Drift between REVIEW.md line references and current plan
  content.**  If the v9.3.0 plan has been edited since REVIEW.md was
  generated, the diff contexts in Phases 1–3 may no longer match.
  Mitigation: before every `Edit`, use `Read` to re-verify the
  exact "old_string" text.  If it has drifted, update the
  old_string to the current actual text and re-issue the edit.

### Medium risk
- **R-006 staging-language drift.**  Review text written for the old
  dual-mode design can accidentally reintroduce a user-facing fallback
  CLI into the plan.  Mitigation: keep the total-only migration note at
  the top of this review-incorporation plan and re-check the final
  Phase 4 wording after every edit.

- **R-013 constraint relaxation perceived as scope creep.**
  Reviewers may read "append-only in the non-DIAG region" as a
  loosening of the v9.3.0 data-layout promise.  Mitigation: the
  replacement text explicitly ties the relaxation to the existing
  `#ifdef SEAS_DIAG_FAULT_FLUX` tail — no new layout variability
  is introduced, only honest acknowledgment of what is already there.

### Low risk
- **Review-incorporation log duplication with REVIEW.md.**  The log
  table summarises REVIEW.md at a high level; if REVIEW.md is later
  superseded (by a rev-2 review), the log becomes stale.
  Mitigation: the log explicitly cites "REVIEW.md (14 items)" so a
  future reviewer knows the source version.

- **Markdown rendering of the diff blocks.**  The plan uses
  ```diff ... ``` code fences.  Some PDF builds (pandoc + xelatex)
  render these correctly; others do not.  Mitigation: verify
  rendering with `pandoc ... --pdf-engine=xelatex` post-patch (same
  command used to build the walkthrough PDF).

## Scope and timeline

| Phase | Scope | Plan-text LOC (added/modified) | Est. wall-clock |
|---|---|---:|---:|
| 1 | CRITICAL R-001..R-004 | ~150 / ~90 | 1 day |
| 2 | MODERATE R-005..R-011 | ~220 / ~120 | 1.5 day |
| 3 | LOW R-012..R-014      |  ~40 /  ~30 | 0.5 day |
| 4 | Review-incorporation log | ~80 / 0   | 0.5 day |
| **Total** | | **~490 / ~240** | **~3.5 days** |

## Quick-reference for implementation agent

- Every `Edit` call in Phases 1–3 uses a `diff`-style old_string →
  new_string pair matching a REVIEW.md "Suggested fix" block
  verbatim.  Do NOT paraphrase; copy-paste the review's proposed
  text.
- For each CRITICAL finding, the plan edit is blocking on
  verification (use `Read` to confirm the patch landed, then
  `Grep` to confirm the regression gate described in the Acceptance
  Criteria).
- Phase 4's review-incorporation log is a pure append — use `Edit`
  with the existing last line of the plan as `old_string` and the
  existing last line + the log content as `new_string`.
- Do NOT touch any source code (`*.cpp`, `*.hpp`, `Makefile`) in
  this plan.  Every edit is confined to the plan `.md` file.
- After all 4 phases land, run `pandoc tpv102_debug_v9.3.0_debug_plan.md
  -o tpv102_debug_v9.3.0_debug_plan.pdf --pdf-engine=xelatex
  -V geometry:margin=1in -V fontsize=10pt -V mainfont="Helvetica"
  -V monofont="Menlo" --syntax-highlighting=tango --toc --toc-depth=2`
  to produce a PDF of the rev-2 plan for distribution.
