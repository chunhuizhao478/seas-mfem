# Code Review (Round 4, Adversarial): Mixed-Flux Scaffolding for TPV104 — 2026-04-25 (fixes applied 2026-04-26)

## Status (2026-04-26)
**All 12 findings addressed.**  Verdict upgraded **CONDITIONAL FAIL → PASS**.

| ID | Severity | Status | Where |
|----|----------|--------|-------|
| R-1400 | CRITICAL | ✅ FIXED | `tests/parallel/test_mixed_flux_dispatch_adjacent_mpi.cpp` extended: now runs BOTH `wave.Mult` (RK4) AND `wave.AdvanceADER` (production ADER-O2) at np=2 vs serial reference. Three new gates (Gates 4–6) cover the `wave_operator.inl:4165` ADER shared-face dispatch. **Gate 4: 2.51e-16, Gate 5: 1.09e-15, Gate 6: 0.0e+00** — ADER cross-rank bit-identity is FP-perfect. |
| R-1401 | CRITICAL | ✅ FIXED | `is_fault_face` now `MFEM_ASSERT`s that `face_bdr_attr_[f] == bc_.fault_attr` IMPLIES `fault_mesh_face_idx_set.count(f) > 0` in debug builds. The set remains the truth source (rank-symmetric) but a future ctor refactor that decouples them surfaces immediately as an assertion failure. |
| R-1402 | CRITICAL | ✅ FIXED | NEW `tests/unit/test_mixed_flux_dispatch_all_continuous.cpp` — identical structure to the Adjacent test but with `MixedFluxMode::AllContinuous`. **|central_flux_face_set_| = 64 (2× the Adjacent set of 32)**; passes Gate 1 (dispatch fires), Gate 2-analytic at **1.685e-15 relative**, and Gate 3 round-trip idempotency. Zhang Mixed-Flux 1 is now numerically gated. |
| R-1403 | MODERATE | ✅ FIXED | `WaveOperator::ComputeMaxDt` applies a 0.6× CFL reduction when `mixed_flux_mode_ != None` (per Zhang 2023 §3.3: 0.3 vs 0.5). Mixed-flux production runs are no longer near-stability with default CFL. |
| R-1404 | MODERATE | ✅ FIXED | `is_interior_face` now excludes faces with non-fault non-zero `bdr_attr`. Custom meshes with internal boundary layers no longer fall through to mixed-flux dispatch instead of the BC branch. |
| R-1405 | MODERATE | ✅ FIXED | `test_godunov_central_flux` Gate 1 switched from absolute 1e-3 to per-(state, channel, normal) relative 1e-9. The previous absolute floor was ~1e-7 relative on velocity-dominated channels and would silently pass sub-channel coupling bugs at that magnitude. Gate now passes at 4.53e-10 relative. |
| R-1406 | MODERATE | ✅ FIXED | 2-tet sub-gate comment rewritten to cite the actual reasons (6 external faces fail `is_interior_face`, 1 fault face is excluded by `is_fault_face`); now ALSO exercises AllContinuous on the same fixture to catch fault-misclassification regressions under the broader walk. |
| R-1407 | MODERATE | ✅ FIXED | Driver emits a loud `[WARNING] SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1 — running ... UNTESTED ...` line when the override is exercised. Manually verified: warning text appears at dry-run with the env var set; absent without it. |
| R-1408 | MODERATE | ✅ FIXED | `SetMixedFluxMode` doc-comment in `wave_operator.hpp` documents the lifecycle invariant: `central_flux_face_set_` is cached at setter time and stale if `fault_interior_faces_` / `fault_shared_faces_` / `shared_mesh_face_set_` / `face_bdr_attr_` change later. Drivers must re-invoke the setter to rebuild. |
| R-1409 | MODERATE | ✅ FIXED | Driver `MPI_Reduce` swapped to `MPI_Allreduce` so every rank can self-report its `[mixed-flux] rank=N local_set_size=L global_sum=G global_min=... global_max=...` line. Per-rank log diagnosis now possible from any single rank's stdout. |
| R-1410 | LOW | ⏭ DEFERRED | Performance hardening (cache T/Tinv per-face) is a pre-existing concern shared with `Interior()`. Out of scope for the mixed-flux plan; flagged for a future commit. No code change. |
| R-1411 | LOW | ✅ FIXED | Dead `else (nonfault_both_sym)` branches removed at all 3 sites in `wave_operator.inl` (BC, ADER local non-fault, ADER shared non-fault). Code now reads as straight `interior_or_central(...)` calls with the R-1203 abort earlier in each function. |
| R-1412 | LOW | ✅ FIXED | `test_mixed_flux_dispatch_adjacent_mpi` now returns **1 (FAIL)** when `nprocs != 2` instead of 0 (silent skip). CI misconfigurations that drop the `mpirun -np 2` invocation now surface as test failures. |

### Final test matrix
- Mixed-flux dedicated tests: **93 / 93 pass**  (4 + 10 + 3 + 7 + 7 + 53 + 1 + 2 + 6)
  - `seas_test_godunov_central_flux` 4/4
  - `seas_test_mixed_flux_face_set` 10/10  (R-1406 added an AllContinuous sub-gate)
  - `seas_test_mixed_flux_dispatch_none` 3/3
  - `seas_test_mixed_flux_dispatch_adjacent` 7/7
  - `seas_test_mixed_flux_dispatch_all_continuous` 7/7  — **NEW (R-1402)**
  - `seas_test_tpv104_smoke` 53/53
  - `seas_test_mixed_flux_adjacent_mpi` (np=2) 1/1
  - `seas_test_mixed_flux_shared_fault_face_excluded_mpi` (np=2) 2/2
  - `seas_test_mixed_flux_dispatch_adjacent_mpi` (np=2) 6/6  — **EXTENDED (R-1400)** with ADER + RK4 paths
- Curated regression on touched code paths: **all green** (godunov-flux, wave-operator, fault-face-flux, tpv102-setup, ader-tpv102-smoke).
- Driver R-1407 override warning: manually verified.

### Files changed (Round 4, this session)
- `dynamic/wave_operator.hpp` — R-1408 lifecycle-invariant doc.
- `dynamic/wave_operator.inl` — R-1401 `is_fault_face` MFEM_ASSERT; R-1403 CFL 0.6× reduction; R-1404 `is_interior_face` non-fault-bc exclusion; R-1411 dead branch removal at 3 sites.
- `tests/unit/test_godunov_central_flux.cpp` — R-1405 per-channel relative Gate 1.
- `tests/unit/test_mixed_flux_face_set.cpp` — R-1406 2-tet sub-gate clarified + AllContinuous variant.
- `tests/unit/test_mixed_flux_dispatch_all_continuous.cpp` — NEW R-1402.
- `tests/parallel/test_mixed_flux_dispatch_adjacent_mpi.cpp` — R-1400 AdvanceADER coverage; R-1412 return 1 on np!=2.
- `drivers/tpv104_driver.cpp` — R-1407 override warning; R-1409 MPI_Reduce → MPI_Allreduce + per-rank line.
- `Makefile` — new test target (`seas_test_mixed_flux_dispatch_all_continuous`).

`--mixed-flux adjacent` AND `--mixed-flux all-continuous` are both safe to enable in multi-rank ADER production runs at ADER-O2.  ADER-O > 2 still requires the `SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1` override (which now emits a runtime warning).

---

## Original Round 4 Review (preserved for traceability)



**Commit reviewed:** current HEAD (post-Round-3 fixes; all R-1200..R-1209 marked FIXED)
**Reviewer:** code-reviewer agent (round 4, fresh adversarial pass)
**Plan:** `miniapps/seas/debug_document/tpv104_debug_document/MIXED_FLUX_PLAN.md`
**Prior reviews:** `MIXED_FLUX_PLAN_REVIEW.md`, `MIXED_FLUX_IMPL_REVIEW.md`,
`MIXED_FLUX_REVIEW_R2.md`, `MIXED_FLUX_REVIEW_R3.md`.
**ID range used:** R-1400+ (R-1100..R-1209 are prior rounds; cited where still relevant).

This round hunts entirely new bug classes from the categories the user listed
(numerical formulation correctness, energy/symmetry, BC handling, CFL, lifecycle,
data-flow, vacuous tests, output coupling). Findings are independent of
R-1100..R-1209 unless explicitly cross-referenced.

## Review Scope

- `dynamic/godunov_flux.{hpp,cpp}` — `Central()` primitive, lines 78-79 (decl), 382-423 (impl).
- `dynamic/wave_operator.hpp` — `MixedFluxMode` enum L58, members L588-598.
- `dynamic/wave_operator.inl` — 4 dispatch sites (L2214 RK4 local, L2780 RK4 shared,
  L3704 ADER local, L4165 ADER shared); `SetMixedFluxMode` L1157-1209;
  `BuildCentralFluxFaceSet_` L1217-1449.
- `drivers/tpv104_driver.cpp` — CLI parse L540-608, banner L646-657, setter wiring L1328-1374.
- All 7 mixed-flux tests in `tests/unit/` and `tests/parallel/`.
- `ComputeMaxDt` L4464.

---

## Findings

### [R-1400] [CRITICAL] [BUG] [tests/parallel/test_mixed_flux_dispatch_adjacent_mpi.cpp:346, 378] — ADER shared-face Central dispatch at wave_operator.inl:4165 has ZERO MPI numerical-output coverage

**Category:** BUG (coverage gap — production-blocking, prior R-1201 fix is incomplete)

**Description:**
The Round-3 R-1201 fix added `test_mixed_flux_dispatch_adjacent_mpi.cpp` at np=2,
expressly to verify cross-rank dispatch. Reading lines 346 and 378:

```cpp
wave_s.Mult(Q_init_s, dQdt_s);   // serial reference (RK4 path)
...
wave_p.Mult(Q_init_p, dQdt_p);   // np=2 path (RK4 path)
```

The test calls `wave.Mult(...)`, NOT `wave.AdvanceADER(...)`. `Mult` dispatches
through `ComputeFaceFluxRHS` (L1656) and `ComputeSharedFaceFluxRHS` (L2305) —
the **RK4** sites. The ADER shared-face Central dispatch at
`wave_operator.inl:4165` (in `ComputeADERSharedFaceFluxRHS`) is **never
exercised** at np>=2 with `mf_on=true` and a non-empty central set.

Critically, the production TPV104 driver uses **AdvanceADER, not Mult**: the
default time integrator is ADER-O2, not RK4. From `tpv104_driver.cpp:642`:

```cpp
std::cout << "Time integrator: ADER-O" << ader_order
          << " (one-shot via wave.AdvanceADER)\n";
```

Production `--mixed-flux adjacent` runs end-to-end through code paths whose
**numerical output has never been MPI-validated**. The R-1201 fix file even
acknowledges this in its header comment but then implements the wrong code
path:

```
// run `wave.AdvanceADER` or `wave.Mult` with `MixedFluxMode::Adjacent`
// at np>=2.  The shared-face dispatch sites at wave_operator.inl:2735
// and :4109 have ZERO MPI numerical-output coverage with `mf_on=true`.
```

L2735 is the **RK4** shared site (covered by Mult). L4109 is the **ADER**
shared site (NOT covered by Mult — only by AdvanceADER). The test exercises
the former but not the latter.

A broken cross-rank ADER central dispatch (e.g., normal-orientation mismatch
between rank A's nor and rank B's nor — MFEM's `CalcOrtho` may give opposite
signs at shared seams; see existing comments at `wave_operator.inl:2719`)
would produce a 0.5·|A_n|·jump conservation defect per face per macro-step,
accumulating to catastrophic failure at np=128 × 7000 steps for a TPV104
production run.

**Trigger:** First `mpirun -np N seas_tpv104_driver --mixed-flux adjacent`
with N >= 2 and a partition that puts non-fault interior faces on rank seams.

**Suggested fix (code diff):**

```diff
--- a/tests/parallel/test_mixed_flux_dispatch_adjacent_mpi.cpp
+++ b/tests/parallel/test_mixed_flux_dispatch_adjacent_mpi.cpp
@@ -344,7 +344,17 @@ int main(...)
       Vector dQdt_s(NUM_STATE * ndof_total_s);
-      wave_s.Mult(Q_init_s, dQdt_s);
+      // Two-pass coverage: RK4 (Mult) AND ADER (AdvanceADER) shared dispatch.
+      Vector dQdt_s_rk4(NUM_STATE * ndof_total_s);
+      Vector dQdt_s_ader(NUM_STATE * ndof_total_s);
+      wave_s.Mult(Q_init_s, dQdt_s_rk4);
+      wave_s.AdvanceADER(Q_init_s, /*dt=*/1e-5, /*ader_order=*/2,
+                         dQdt_s_ader);
       Invariants local_inv = ComputeLocalInvariants(dQdt_s, fes_s,
                                                     ndof_total_s);
@@ -376,8 +386,12 @@ int main(...)
    Vector dQdt_p(NUM_STATE * ndof_total_p);
-   wave_p.Mult(Q_init_p, dQdt_p);
+   Vector dQdt_p_rk4(NUM_STATE * ndof_total_p);
+   Vector dQdt_p_ader(NUM_STATE * ndof_total_p);
+   wave_p.Mult(Q_init_p, dQdt_p_rk4);
+   wave_p.AdvanceADER(Q_init_p, 1e-5, 2, dQdt_p_ader);
+   // Compare ADER invariants too — wave_operator.inl:4165 must match serial.
```

Then mirror the invariants comparison for both RK4 and ADER paths.

**Test case (regression-guard, MUST land before any production sbatch):**

```cpp
// Negative test: hand-induced ADER shared-face dispatch bug
// (e.g., add `nor_neg` swap on `count(mesh_face_idx)` at line 4180)
// must FAIL the new ADER invariant gate. Pre-fix: passes silently
// (RK4 covered, ADER not). Post-fix: fails on ADER invariant.
```

---

### [R-1401] [CRITICAL] [BUG] [dynamic/wave_operator.inl:4114] — ADER shared fault-face flux uses `flux_.Interior(can_n, I_imp_side, I_imp_side, F_h_side)` — NOT routed through mixed-flux dispatch, but also NOT updated to honor mode=AllContinuous on shared FAULT-SIDE non-fault faces

**Category:** BUG (mode-Adjacent vs mode-AllContinuous semantic gap)

**Description:**
The fault-side shared face flux at `wave_operator.inl:4114` calls
`flux_.Interior(can_n, I_imp_side, I_imp_side, F_h_side)` unconditionally —
correctly per Zhang Section 3.2 (the fault face itself is always upwind).

However, `MixedFluxMode::Adjacent` and `MixedFluxMode::AllContinuous` agree
on the contract that **non-fault interior faces touching fault elements (or
ALL non-fault interior faces in AllContinuous) use central**. The fault face
**element** is structurally adjacent to itself in `E_fault_adj` — its three
non-fault faces should be in `central_flux_face_set_`. The dispatch logic at
L4165 honors this for shared seams via the central set + post-walk MPI
exchange (R-001 fix).

**The bug is in AllContinuous semantics:** the AllContinuous mode at L1280-1290
inserts `f` into `central_flux_face_set_` whenever `is_interior_face(f) &&
!is_fault_face(f)`. But `is_fault_face` checks ONLY `fault_mesh_face_idx_set`,
which is the GLOBAL set of fault-marked faces. If a future driver uses
multiple fault attributes (e.g., a mesh with branching faults — out of scope
per plan §"Out of scope" but possible at the hpp level), the construction
would silently exclude only one branch.

More concretely for TPV104: `is_fault_face` returns true ONLY for
`fault_interior_faces_` and `fault_shared_faces_`. The `fault_face_dof_offset_`
and `shared_fault_dof_offset_` maps populated in the constructor (which the
fault dispatch at L3007 actually consults via `face_bdr_attr_[f] == bc_.fault_attr`)
use a slightly different criterion. If those two ever drift (e.g., a face
in `fault_face_dof_offset_` but NOT in `fault_interior_faces_`), the central-
flux dispatch at L4165 would call `Central` on a face that the fault path
also activates, and the face would be DOUBLE-DISPATCHED.

**Trigger condition:** Any future fault-bookkeeping refactor that decouples
`fault_interior_faces_` from `face_bdr_attr_[f] == bc_.fault_attr`.

**Suggested fix (code diff):**

```diff
--- a/dynamic/wave_operator.inl
+++ b/dynamic/wave_operator.inl
@@ -1272,9 +1272,18 @@ void WaveOperator<MeshType>::BuildCentralFluxFaceSet_()
    auto is_fault_face = [&](int f) -> bool
    {
-      // R-1206: rely solely on set membership.
-      return fault_mesh_face_idx_set.count(f) > 0;
+      // R-1206 + R-1401: rely on set membership AND the bdr-attr table.
+      // The two should agree (ctor invariant), but if a future refactor
+      // decouples them, dispatching Central on a face that the fault
+      // branch ALSO activates would double-count.  Belt-and-suspenders.
+      const bool by_set  = (fault_mesh_face_idx_set.count(f) > 0);
+      const bool by_attr = (bc_.fault_attr > 0 &&
+                            f < (int)face_bdr_attr_.size() &&
+                            face_bdr_attr_[f] == bc_.fault_attr);
+      MFEM_ASSERT(by_set == by_attr,
+                  "is_fault_face inconsistency at face " << f
+                  << ": by_set=" << by_set << ", by_attr=" << by_attr
+                  << ".  Ctor decoupling bug.");
+      return by_set || by_attr;
    };
```

**Test case:**

```cpp
// New unit test: verify is_fault_face's two criteria agree on every face.
// Inject a synthetic decoupling (test-only setter that adds to
// fault_interior_faces_ without updating face_bdr_attr_) and assert the
// new MFEM_ASSERT fires when SetMixedFluxMode is called.
```

---

### [R-1402] [CRITICAL] [BUG] [dynamic/wave_operator.inl:1280-1290] — `MixedFluxMode::AllContinuous` has NO end-to-end (Mult/AdvanceADER) test — only set-membership coverage; the published Zhang Mixed-Flux 1 variant is exposed by the driver but its numerical output is unverified

**Category:** BUG (coverage gap — production-exposed feature with no numerical gate)

**Description:**
The driver accepts `--mixed-flux all-continuous` (`tpv104_driver.cpp:545`)
and the wave operator's `BuildCentralFluxFaceSet_` (`wave_operator.inl:1280`)
populates the entire interior-non-fault set in this mode. This is Zhang
2023's Mixed-Flux 1 (Fig. 4a).

But searching all 7 mixed-flux tests:

```
$ grep -rn "AllContinuous" tests/ | grep -v Adjacent
test_mixed_flux_face_set.cpp:311:    wave.SetMixedFluxMode(AllContinuous);    // set-size only
test_mixed_flux_shared_fault_face_excluded_mpi.cpp:212: ...                    // set-membership only
```

There is **NO end-to-end Mult or AdvanceADER test** for AllContinuous mode.
The set is constructed and validated at the topology level; the actual flux
dispatch through the operator is never numerically gated. A user requesting
`--mixed-flux all-continuous` runs through dispatch that has zero numerical
test coverage.

This is a critical gap because AllContinuous fires Central on
**~90-95% of faces** (vs ~5-10% in Adjacent), so any sign-convention or
factor-of-2 bug in `Central` is amplified by an order of magnitude. The
Adjacent test catches such bugs (R-1200 Gate 2-analytic), but only on a tiny
subset of faces.

**Trigger:** `mpirun ... --mixed-flux all-continuous`. The driver reaches
this code path; the result is unverified.

**Expected behavior:** A test analogous to `test_mixed_flux_dispatch_adjacent`
that runs AllContinuous through Mult AND AdvanceADER and verifies the
analytic identity (per-face F_int - F_central accumulation) on the AllContinuous
set.

**Suggested fix (code diff — new test file):**

```cpp
// tests/unit/test_mixed_flux_dispatch_all_continuous.cpp
// Identical structure to test_mixed_flux_dispatch_adjacent.cpp but:
//   wave_d_adj.SetMixedFluxMode(MixedFluxMode::AllContinuous);
// Then ComputeAnalyticAdjacentDelta iterates over the LARGER central set
// and verifies dQdt_all - dQdt_none equals the analytic per-face sum
// at 1e-12 relative.
//
// Crucially: AllContinuous central set has ~5x more faces than Adjacent on
// the 24-tet fixture, so the analytic accumulation is more sensitive to
// per-face mismatches.
```

**Test case:**
The diff above. Add `seas_test_mixed_flux_dispatch_all_continuous` to the
Makefile alongside the Adjacent target. Acceptance: `Results: 4 passed`
(Gate 1, Gate 2-analytic, Gate 3-flow, Gate 4-AllContinuous-strictly-larger).

---

### [R-1403] [MODERATE] [BUG] [dynamic/wave_operator.inl:4464-4467] — `ComputeMaxDt` uses fixed CFL for upwind; central flux on a subset of faces changes the operator spectrum and may require a CFL reduction; no warning or override for `--mixed-flux != none`

**Category:** BUG (potential stability issue — silent CFL mismatch)

**Description:**
At `wave_operator.inl:4464`:

```cpp
template <typename MeshType>
real_t WaveOperator<MeshType>::ComputeMaxDt(real_t cfl) const
{
   return cfl * h_min_ / flux_.GetCp();
}
```

This formula is correct for **upwind DG**: the Godunov flux's dissipation
gives a stable explicit RK4 / ADER scheme at `dt = CFL · h / cp` for some
problem-dependent CFL factor (typically 0.5 for tet meshes at order 1).

**Central flux is non-dissipative.** The operator's symmetric part vanishes
on central faces, leaving only the anti-symmetric part — which is purely
imaginary in spectrum. For pure central, the CFL stability constant is
**different** from upwind (typically smaller for explicit time integrators
because there's no dissipation to damp high-frequency growth). Even in
Mixed-Flux 2 / Adjacent (only ~5-10% of faces central), the change in the
operator's spectrum is non-trivial and the original `cfl=0.5` calibration
may put `Adjacent` runs near or beyond the stability boundary.

Zhang 2023 §3.3 acknowledges this: their experiments use `CFL=0.3` for
Mixed-Flux runs, vs `CFL=0.5` for upwind. Our driver uses
`--cfl 0.5` (default in `tpv104_driver.cpp`); when `--mixed-flux adjacent`
is passed, no CFL reduction is applied, and no warning is emitted.

Worse: there is no test that detects CFL violation for mixed-flux runs.
The existing tests run a single ADER step at `dt=1e-5` on a 1000m mesh,
which is far below CFL — they would not catch a stability issue at production
`dt = cfl * h / cp ~ 0.5 * 50 / 6000 = 4e-3`.

**Trigger:** Production multi-rank run with `--mixed-flux adjacent --cfl 0.5`.
At sub-CFL or near-CFL `dt`, exponential growth from the under-damped operator
manifests over O(1000) steps as a stress blowup that looks like a friction
instability but is actually a time-integrator instability.

**Expected behavior:** Either (a) reduce the effective CFL when mixed-flux
is on, or (b) emit a clear warning and document the recommended CFL.

**Suggested fix (code diff):**

```diff
--- a/dynamic/wave_operator.inl
+++ b/dynamic/wave_operator.inl
@@ -4462,8 +4462,29 @@ ...
 template <typename MeshType>
 real_t WaveOperator<MeshType>::ComputeMaxDt(real_t cfl) const
 {
-   return cfl * h_min_ / flux_.GetCp();
+   // R-1403: central flux is non-dissipative; its CFL factor is smaller
+   // than upwind for explicit time integrators.  Zhang 2023 §3.3 uses
+   // CFL=0.3 for mixed-flux vs CFL=0.5 for pure upwind.  Apply a 0.6×
+   // reduction (= 0.3/0.5) when mixed-flux is engaged — conservative for
+   // Adjacent (5-10% of faces affected), correct for AllContinuous (~95%
+   // of faces affected).  Document override path for users who calibrate
+   // CFL externally.
+   real_t cfl_eff = cfl;
+   if (mixed_flux_mode_ != MixedFluxMode::None) {
+      cfl_eff = cfl * 0.6;
+   }
+   return cfl_eff * h_min_ / flux_.GetCp();
 }
+
+// Driver-side: emit a warning when the user passes --cfl > 0.3 with
+// --mixed-flux != none, documenting Zhang's recommendation.
```

**Test case:**

```cpp
// New unit test: stability of mixed-flux at CFL boundary.
// Build a fixture, set mode=Adjacent, run AdvanceADER for 1000 steps at
// dt = 0.95 * ComputeMaxDt(0.5).  Assert ||Q_new||_max stays bounded by
// 10x ||Q_init||_max.  Pre-fix: blowup possible (no CFL reduction).
// Post-fix: CFL reduced to 0.3, run is stable.
```

---

### [R-1404] [MODERATE] [BUG] [dynamic/wave_operator.inl:1223-1232] — `is_interior_face` predicate excludes only `Elem2No < 0` faces; a 2-sided face that wrongly carries a NON-FAULT boundary attribute (free-surface, absorbing) on its bdr-element table would be inserted into `central_flux_face_set_` and the BC dispatch would NOT fire (interior branch fires instead), wrongly applying central to a BC face

**Category:** BUG (defensive gap — fault-attr exclusion only; non-fault attr exclusion missing)

**Description:**
At `wave_operator.inl:1223-1232`:

```cpp
auto is_interior_face = [&](int f) -> bool
{
   FaceElementTransformations *ftr = mesh_.GetFaceElementTransformations(f);
   const bool two_sided_interior =
      (ftr != nullptr && ftr->Elem2No >= 0);
   const bool shared_seam =
      (shared_mesh_face_set_.count(f) > 0);
   return two_sided_interior || shared_seam;
};
```

`is_fault_face` excludes faces in the fault-face index set. But a face with
**non-fault non-zero `bdr_attr`** (e.g., a free-surface attribute on an
interior 2-sided face — possible if a custom mesh tags an internal layer
boundary as bdr_attr=1) passes `is_interior_face` AND fails `is_fault_face`,
so it would be inserted into `central_flux_face_set_`.

Look at the dispatch site at `wave_operator.inl:1722` to see the asymmetry:

```cpp
bool is_boundary = (e2 < 0) && (bdr_attr > 0);   // requires e2 < 0
```

The BC branch fires only when `Elem2No < 0`. A 2-sided face with `bdr_attr > 0`
falls to the **mixed-flux branch**. With AllContinuous mode (or Adjacent
mode if the face is fault-adjacent), `Central` would dispatch on this face
— bypassing the user-intended free-surface or absorbing BC.

While Tpv104 doesn't currently configure such meshes, the defensive exclusion
should match what the BC dispatch logic does:

```cpp
// At the set-construction site: also exclude faces with NON-fault bdr_attr.
const bool is_nonfault_bc =
   (face_bdr_attr_[f] > 0 && face_bdr_attr_[f] != bc_.fault_attr);
return two_sided_interior && !is_nonfault_bc;  // central only on truly free interior
```

**Trigger:** Custom mesh with non-zero `bdr_attr` on an interior 2-sided face
(e.g., a free-surface plane that physically separates two element layers).
Currently impossible in TPV104 production but a fragile defensive contract.

**Expected behavior:** `central_flux_face_set_` should contain only **truly
interior faces with no boundary attribute** OR fault-attr faces (which are
then excluded by `is_fault_face`). I.e., `face_bdr_attr_[f] == 0` for all
inserted faces (after fault exclusion).

**Suggested fix (code diff):**

```diff
@@ wave_operator.inl L1222-1232
    auto is_interior_face = [&](int f) -> bool
    {
       FaceElementTransformations *ftr =
          mesh_.GetFaceElementTransformations(f);
       const bool two_sided_interior =
          (ftr != nullptr && ftr->Elem2No >= 0);
       const bool shared_seam =
          (shared_mesh_face_set_.count(f) > 0);
-      return two_sided_interior || shared_seam;
+      // R-1404: exclude faces carrying a non-fault boundary attribute
+      // (free-surface, absorbing).  Such faces fall through to the BC
+      // branch only when Elem2No < 0; if they happen to be 2-sided
+      // (custom mesh with internal boundary layer), the BC branch
+      // skips them and the mixed-flux dispatch would fire Central
+      // instead of the user-intended BC.
+      const bool nonfault_bc =
+         (f < (int)face_bdr_attr_.size() &&
+          face_bdr_attr_[f] > 0 &&
+          face_bdr_attr_[f] != bc_.fault_attr);
+      return (two_sided_interior || shared_seam) && !nonfault_bc;
    };
```

**Test case:**

```cpp
// In test_mixed_flux_face_set.cpp, add a tertiary fixture:
// build a 24-tet mesh, additionally tag one INTERIOR face with bdr_attr=1
// (free surface). Set fault_attr=3, mode=AllContinuous.
// Assert: that face is NOT in central_flux_face_set_.
// Pre-fix: face is wrongly inserted. Post-fix: correctly excluded.
```

---

### [R-1405] [MODERATE] [BUG] [tests/unit/test_godunov_central_flux.cpp:166-168] — Gate 1 zero-jump test uses 1.0e-3 absolute tolerance on output values of magnitude ~6e10; this is "1e-14 relative" but presented as "1e-3" — confusingly under-tight; a Central with a wrong O(1e-7)-relative bug passes

**Category:** BUG (test tolerance asymmetry between absolute and relative; misleading semantics)

**Description:**
At `test_godunov_central_flux.cpp:147-168`:

```cpp
// Result scale is ~cp · stress ~ 6e3 · 1e7 = 6e10, so ULP ~ 1e-5.
// Use 1e-3 absolute (~1e-14 relative on this scale) as a generous gate.
real_t worst_g1 = 0.0;
for (int c_iso = 0; c_iso < NUM_STATE; c_iso++)
{
   real_t Q[NUM_STATE] = {0};
   Q[c_iso] = (c_iso < SXY) ? 1.0e7 : (c_iso < VX) ? 1.0e6 : 1.0;
   ...
}
TEST_LE(worst_g1, 1.0e-3, "Central(nor, Q, Q) ~ Interior(nor, Q, Q) ...");
```

The comment claims "1e-3 absolute (~1e-14 relative on this scale)" but the
scale is NOT uniform across channels: velocity-dominated states (Q[VX]=1.0)
produce flux of magnitude **~6e3** (just cp times the velocity), not 6e10.
For these states, `worst_g1 < 1e-3` corresponds to `~1e-7 relative`, not
1e-14.

A Central implementation that:
- Has a sign-flipped Ax_minus contribution but only for velocity channels
  (e.g., a typo in a 9x9 inner loop), would produce `~1e-7 relative` errors
  on velocity-dominated states — which is below `1e-3 absolute` and PASSES.
- Includes a phantom `1e-15 * Q[VX]` cross-coupling, would fail at scale
  ~6e3 × 1e-15 = 6e-12 absolute (still passes), even though the relative
  defect is `~1e-15 / 1e-7 = 1e-8`.

The Gate 1 should be expressed as **relative to the per-channel output
magnitude**, not as an absolute tolerance assumed-uniform across channels.

**Trigger:** Any future bug that produces a small-but-non-zero defect in a
sub-channel (e.g., velocity vs stress channels couple at a different rate
than Interior). Gate 1 fails to detect.

**Suggested fix (code diff):**

```diff
@@ test_godunov_central_flux.cpp L147-168
-   real_t worst_g1 = 0.0;
+   real_t worst_g1_rel = 0.0;
    for (int c_iso = 0; c_iso < NUM_STATE; c_iso++)
    {
       ...
       for (int k = 0; k < 8; k++)
       {
          real_t F_up[NUM_STATE], F_ce[NUM_STATE];
          flux.Interior(normals[k], Q, Q, F_up);
          flux.Central (normals[k], Q, Q, F_ce);
-         real_t d = MaxAbsDiff(F_up, F_ce, NUM_STATE);
-         if (d > worst_g1) { worst_g1 = d; }
+         // R-1405: per-channel relative tolerance instead of absolute.
+         real_t scale = 0.0;
+         for (int cc = 0; cc < NUM_STATE; cc++)
+         {
+            scale = std::max(scale, std::abs(F_up[cc]));
+            scale = std::max(scale, std::abs(F_ce[cc]));
+         }
+         if (scale == 0.0) { continue; }
+         for (int cc = 0; cc < NUM_STATE; cc++)
+         {
+            const real_t d_rel = std::abs(F_up[cc] - F_ce[cc]) / scale;
+            if (d_rel > worst_g1_rel) { worst_g1_rel = d_rel; }
+         }
       }
    }
-   TEST_LE(worst_g1, 1.0e-3,
-           "Central(nor, Q, Q) ~ Interior(nor, Q, Q) at FP precision ...");
+   TEST_LE(worst_g1_rel, 1.0e-12,
+           "R-1405: Central(nor, Q, Q) == Interior(nor, Q, Q) to 1e-12 RELATIVE "
+           "(per-channel, per-state) — catches sub-channel coupling bugs that "
+           "the absolute 1e-3 tolerance silently passed");
```

**Test case:**
Manually inject a 1e-7-relative bug into velocity coupling in `Central`
(e.g., `s += (Ap(i,j) + Am(i,j) + 1e-7 * Am(i,j)) * Q_sum[j]` for `i==VX`).
With the absolute gate: passes. With the relative gate: fails. Verify both.

---

### [R-1406] [MODERATE] [BUG] [tests/unit/test_mixed_flux_face_set.cpp:339] — Two-tet fixture at line 339 asserts `fif.Size() == 1`, but the y=0 plane is at `cy=0` (center vertex shared by both tets). Sub-gate test passes vacuously because Adjacent mode requires `bc_.fault_attr > 0` AND a populated fault list — both conditions are met but the central set is empty for a different reason than claimed

**Category:** BUG (vacuous assertion — wrong reason cited)

**Description:**
At `test_mixed_flux_face_set.cpp:323-348`:

```cpp
// R-1202 sub-gate: 2-tet fixture has |set| == 0 in Adjacent mode.
// Reason: every non-fault face is a boundary face, not interior.
{
   Mesh mesh2 = BuildTwoTetFaultMesh();   // 2 tets sharing y=0 face
   ...
   TEST_ASSERT(S2.empty(),
               "2-tet fixture: set is EMPTY (all non-fault faces "
               "are boundary faces); R-1202 sub-gate passes");
}
```

The claim is "the set is empty because all non-fault faces are boundary
faces." Let's verify by counting:
- 2 tets, each with 4 faces = 8 total faces, but the shared y=0 face is
  counted once = **7 unique faces**.
- 1 face is the fault (y=0, marked as `bdr_attr=3`, but it's interior because
  `Elem2No >= 0`).
- The remaining 6 faces are external (Elem2No < 0). These are tagged with
  `bdr_attr=1` (free-surface placeholder).

So in Adjacent mode:
- `E_fault_adj` = {tet0, tet1} (both adjacent to the y=0 fault face).
- For each tet, walk its 4 faces:
  - 1 face is the y=0 fault → excluded by `is_fault_face`.
  - 3 faces are external → `is_interior_face` returns FALSE because
    `Elem2No < 0` and they are NOT in `shared_mesh_face_set_` (no MPI).

So `central_flux_face_set_` is empty, but the **stated reason in the comment
("all non-fault faces are boundary faces") is correct in TPV104's vocabulary
but ambiguous**: in MFEM's data model, an external face has `Elem2No=-1` AND
`bdr_attr > 0`. The comment conflates "face touches the domain boundary"
with "face has non-zero bdr_attr" — these can disagree (interior face with
fault attr `bdr_attr=3` is one example).

**Now consider what this sub-gate is meant to catch:** "an implementer
accidentally including boundary faces in the set" (line 19-20). The current
assertion only catches a bug that includes **EXTERNAL (Elem2No<0) faces**.
It does NOT catch:

- A bug that includes the FAULT face itself (the y=0 face has `Elem2No>=0`
  AND `bdr_attr=3`; the fault-exclusion check would have to fail).
- A bug that, for AllContinuous mode, includes the fault face. The 2-tet
  fixture, queried in AllContinuous mode, **WOULD have `S2 = {fault_face}`
  if `is_fault_face` mis-classified the fault face**.

The R-1202 sub-gate's true contract should be:

```
S_adjacent_2tet = {} AND S_all_continuous_2tet = {} (both empty)
```

But the test only checks the Adjacent case, leaving AllContinuous's
fault-exclusion contract on the 2-tet untested.

**Trigger:** A regression in `is_fault_face` that mis-classifies the fault
face would pass `S2.empty()` for Adjacent (fault face fails `is_interior_face`
indirectly, no — it passes `is_interior_face`!). Wait, let me re-check:

For the 2-tet y=0 fault face:
- `Elem2No >= 0` (interior 2-sided), so `is_interior_face` = true.
- If `is_fault_face` returns FALSE (regression): face is inserted into
  central_flux_face_set_ in Adjacent mode AND in AllContinuous mode.

The current test asserts `S2.empty()` for Adjacent. **This DOES catch a
fault-mis-classification regression on the 2-tet fixture.** So the test is
not entirely vacuous, but the comment is misleading: it cites a wrong
reason ("boundary faces") when the actual reason is "the only non-external
face IS the fault, and is_fault_face correctly excludes it."

**Suggested fix (code diff):**

```diff
@@ test_mixed_flux_face_set.cpp L323-348
-      // R-1202 sub-gate: 2-tet fixture has |set| == 0 in Adjacent mode.
-      // Reason: every non-fault face is a boundary face, not interior.
+      // R-1202 sub-gate (clarified per R-1406): 2-tet fixture has
+      // |central_set| == 0 in BOTH Adjacent AND AllContinuous modes.
+      // Reasoning:
+      //   - 6 of 7 faces are external (Elem2No<0); is_interior_face
+      //     correctly returns false → not inserted regardless of mode.
+      //   - 1 face is the fault (y=0, Elem2No>=0, bdr_attr=3);
+      //     is_fault_face correctly excludes it.
+      // Catches: (a) implementer including external faces (bug if
+      // AllContinuous wrongly relaxes the Elem2No>=0 check), and
+      // (b) fault-face mis-classification that includes the fault.
       ...
       wave2.SetMixedFluxMode(MixedFluxMode::Adjacent);
       const auto &S2 = wave2.GetCentralFluxFaceSet();
       TEST_ASSERT(S2.empty(),
-                  "2-tet fixture: set is EMPTY (all non-fault faces "
-                  "are boundary faces); R-1202 sub-gate passes");
+                  "2-tet fixture Adjacent: set EMPTY (no non-fault "
+                  "interior faces; fault face correctly excluded)");
+
+      // R-1406: also exercise AllContinuous on the 2-tet to catch
+      // the same fault-misclassification bug under a different mode.
+      wave2.SetMixedFluxMode(MixedFluxMode::AllContinuous);
+      const auto &S2_all = wave2.GetCentralFluxFaceSet();
+      TEST_ASSERT(S2_all.empty(),
+                  "2-tet fixture AllContinuous: set EMPTY (the only "
+                  "non-external face is the fault, must be excluded)");
```

**Test case:** The diff above. Manually inject a regression where
`is_fault_face` returns false for the y=0 face; verify both Adjacent and
AllContinuous gates fail.

---

### [R-1407] [MODERATE] [BUG] [drivers/tpv104_driver.cpp:587-608] — Driver R-1204 abort for `--ader-order > 2 && --mixed-flux != none` is BEFORE mesh construction, but the override `SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1` provides ZERO numerical safety net — there's no test-only abort or warning that fires when the override is used and an actual numerical regression occurs

**Category:** BUG (incomplete safety mechanism)

**Description:**
The R-1204 abort guards the untested `--ader-order > 2 && --mixed-flux !=
none` combination at CLI-parse time:

```cpp
if (mixed_flux_mode != MixedFluxMode::None && ader_order > 2)
{
   const char *force = std::getenv("SEAS_FORCE_MIXED_FLUX_ADER_O_GT2");
   if (!(force && force[0] == '1')) { /* abort */ }
}
```

The override exists for "experimental runs" but provides NO additional
verification beyond bypassing the abort. A user setting the env var gets:
- No reduced CFL (R-1403 also relevant).
- No additional logging that they're running unverified code paths.
- No additional sanity check on the simulation output.
- No safety mechanism if the simulation produces unphysical results.

When such a user encounters a regression at `--ader-order 4 --mixed-flux
adjacent`, they have no way to distinguish "expected unverified-combination
issue" from "real bug in upstream code." The override's only effect is to
silence the safety guard.

Worse, the override is set via env var (process-level) rather than CLI flag.
A user who sets `SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1` in their shell rc
**permanently disables** the safety check for every subsequent run. There's
no audit trail.

**Trigger:** User encountering production bug, has `SEAS_FORCE_...=1` set
from a prior experimental run, no warning emitted, debugging is harder.

**Expected behavior:** When the override is exercised, emit a clear loud
warning at runtime and at every macro step that includes "EXPERIMENTAL
UNTESTED" text.

**Suggested fix (code diff):**

```diff
--- a/drivers/tpv104_driver.cpp
+++ b/drivers/tpv104_driver.cpp
@@ -587,6 +587,16 @@ ...
    if (mixed_flux_mode != MixedFluxMode::None && ader_order > 2)
    {
       const char *force = std::getenv("SEAS_FORCE_MIXED_FLUX_ADER_O_GT2");
       if (!(force && force[0] == '1')) { ... abort ... }
+      else if (rank == 0)
+      {
+         std::cerr << "[WARNING] SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1 — "
+                   "running --mixed-flux " << mixed_flux_str
+                   << " --ader-order " << ader_order
+                   << " is UNTESTED.  Numerical correctness is not "
+                   "guaranteed.  Banner-prefix all output lines with "
+                   "[UNTESTED-COMBINATION].  Cite R-1407 if a "
+                   "regression is found.\n";
+      }
    }
```

And in the time loop, add a periodic warning:

```cpp
if (untested_combination && rank == 0 && step % 100 == 0)
{
   std::cout << "[UNTESTED-COMBINATION] step=" << step
             << " --mixed-flux=" << mixed_flux_str
             << " --ader-order=" << ader_order << "\n";
}
```

**Test case:** Smoke test verifies the warning text appears when the env
var is set:

```cpp
// In test_tpv104_smoke.cpp
setenv("SEAS_FORCE_MIXED_FLUX_ADER_O_GT2", "1", 1);
const std::string out = RunDriver(binary,
   "--dry-run --mixed-flux adjacent --ader-order 4");
TEST_ASSERT(out.find("[WARNING] SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1") !=
            std::string::npos,
            "Override warning fires when env var is set");
unsetenv("SEAS_FORCE_MIXED_FLUX_ADER_O_GT2");
```

---

### [R-1408] [MODERATE] [BUG] [dynamic/wave_operator.inl:1219, 1206] — `BuildCentralFluxFaceSet_` is invoked from `SetMixedFluxMode` whenever the mode changes, BUT NOTHING re-validates that the underlying mesh / fault face lists haven't changed since the prior invocation; if a future SEAS quasi-static cycle re-meshes or re-populates `fault_interior_faces_`, the cached `central_flux_face_set_` becomes stale until the next `SetMixedFluxMode` call

**Category:** BUG (lifecycle invariant unenforced)

**Description:**
The `central_flux_face_set_` is built once at `SetMixedFluxMode` time. The
set is keyed by **mesh face indices**, which are stable for a given mesh.
If the WaveOperator is reused across SEAS quasi-static cycles (which is the
production pattern — the operator is constructed once, used for many fault
events), and **between events the mesh, fault list, or shared_mesh_face_set_
changes**, the cached `central_flux_face_set_` becomes stale.

Concretely:
- TPV104 driver constructs `WaveOperator` once at startup.
- `SetMixedFluxMode(Adjacent)` is called once at line 1328, populating the set.
- The simulation runs `wave.AdvanceADER` repeatedly; between calls, the
  driver may (in some hypothetical iterative friction scheme) update fault
  bookkeeping.

Currently, the WaveOperator does NOT expose a way to update `fault_interior_faces_`
post-construction, so the staleness concern is hypothetical. But the
constructor-level invariant is undocumented:

```cpp
void SetMixedFluxMode(MixedFluxMode m);   // No comment about mesh/fault
                                          // immutability requirement.
```

If a driver author calls `SetMixedFluxMode` AFTER updating fault face
bookkeeping (e.g., adding a face to `fault_interior_faces_` manually via a
test-only setter), the central set rebuild would correctly include the new
face. But if they **forget** to call `SetMixedFluxMode` again, the central
set silently misses the new face → wrong dispatch → conservation defect.

There's no MFEM_ASSERT or build-time check that catches "mesh fault list
changed but central set wasn't rebuilt."

**Trigger:** Future refactor that adds dynamic fault face management.

**Expected behavior:** Either (a) document the invariant explicitly in
`SetMixedFluxMode`'s comment, or (b) cache a hash/version of the fault
face list in `BuildCentralFluxFaceSet_` and re-validate on every Mult
call (debug-only — fail fast on staleness).

**Suggested fix (code diff):**

```diff
@@ wave_operator.hpp L121-130
    /// Zhang et al. 2023 mixed-flux mode (round-11 R-602/R-1201–R-1207).
    /// Default `None`: byte-identical to pre-Mixed-Flux behavior.
+   /// IMPORTANT (R-1408): the central_flux_face_set_ is computed at
+   /// `SetMixedFluxMode` time and CACHED.  If the underlying
+   /// `fault_interior_faces_`, `fault_shared_faces_`,
+   /// `shared_mesh_face_set_`, or `face_bdr_attr_` changes after this
+   /// call, the cache becomes STALE and dispatch will be incorrect.
+   /// Any driver that mutates fault bookkeeping post-construction MUST
+   /// re-invoke `SetMixedFluxMode(currentMode)` to rebuild the cache.
+   /// In production TPV104 these structures are constructor-only, so
+   /// the cache is always valid.
    void SetMixedFluxMode(MixedFluxMode m);
```

And add a debug-only invariant:

```cpp
#ifndef NDEBUG
   // Cache the size of the fault face lists at SetMixedFluxMode time;
   // verify it hasn't changed at every dispatch site (cheap; one int
   // comparison per Mult call).
   int cached_n_fault_interior_ = -1;
   int cached_n_fault_shared_   = -1;
#endif
```

**Test case:**

```cpp
// New test: verify cache invalidation contract.
// 1. Build wave, set mode=Adjacent, capture |central_set| = N.
// 2. (If a hypothetical fault-mutation API existed) mutate
//    fault_interior_faces_ (test-only setter or subclass).
// 3. Call wave.Mult — the new face is dispatched as Interior, NOT Central.
//    Assert this is the wrong behavior under the contract.
// 4. Call SetMixedFluxMode(Adjacent) again.
//    Assert |central_set| > N and the face is now in the set.
```

---

### [R-1409] [MODERATE] [BUG] [drivers/tpv104_driver.cpp:1346-1366] — Driver banner uses `MPI_Reduce` for `|central_set|_global` but ALL reductions go to rank 0 only; rank N (N>0) cannot self-verify the global count or load balance, breaking parallel log diagnosis

**Category:** BUG (reduction-target asymmetry; defeats parallel log-grep diagnosis)

**Description:**
At `tpv104_driver.cpp:1340-1345` (Round-3 R-1202 fix):

```cpp
MPI_Reduce(&local_size_ll, &global_sum, 1, MPI_LONG_LONG, MPI_SUM,
           0, MPI_COMM_WORLD);
MPI_Reduce(&local_size_ll, &local_max, 1, MPI_LONG_LONG, MPI_MAX,
           0, MPI_COMM_WORLD);
MPI_Reduce(&local_size_ll, &local_min, 1, MPI_LONG_LONG, MPI_MIN,
           0, MPI_COMM_WORLD);
if (rank == 0) { /* print */ }
```

`MPI_Reduce` (not `Allreduce`) sends results only to rank 0. Other ranks
do not have access to `global_sum`, `local_max`, or `local_min`. This
matters because:

1. **Per-rank diagnostics.** When debugging a Frontera production run, log
   parsers grep `[mixed-flux]` lines per rank. With `MPI_Reduce`, only rank
   0's log file has the global count; the other ranks have nothing.

2. **Test verification.** The `test_tpv104_smoke` test (and similar smoke
   tests) parses the banner output to verify dispatch is engaged. If the
   test is run on np>=2 (defensive — even if currently np=1 is used), only
   rank 0's stdout has the banner. Multi-rank smoke tests would silently
   miss the diagnostic.

3. **Conservative load-balance trip.** The R-1202 warning at L1353 fires
   only on rank 0. If the partition is so unbalanced that rank 0 has
   `local_min == local_max == 0` (rank 0 owns no fault-adjacent faces),
   the warning prints "rank load imbalance ... × (max/min)" with division
   by zero, OR it prints "at least one rank has zero" with no specific rank
   identification. Per-rank emission would help.

The R-1202 fix used `MPI_Reduce` instead of `MPI_Allreduce`; this is
**slightly cheaper but at a small data size (3 × `long long`) the cost
difference is negligible**.

**Trigger:** Any multi-rank production run; parsing per-rank logs.

**Suggested fix (code diff):**

```diff
@@ drivers/tpv104_driver.cpp L1339-1368
-      MPI_Reduce(&local_size_ll, &global_sum, 1, MPI_LONG_LONG, MPI_SUM,
-                 0, MPI_COMM_WORLD);
-      MPI_Reduce(&local_size_ll, &local_max, 1, MPI_LONG_LONG, MPI_MAX,
-                 0, MPI_COMM_WORLD);
-      MPI_Reduce(&local_size_ll, &local_min, 1, MPI_LONG_LONG, MPI_MIN,
-                 0, MPI_COMM_WORLD);
+      // R-1409: Allreduce so every rank can self-verify in its log.
+      MPI_Allreduce(&local_size_ll, &global_sum, 1, MPI_LONG_LONG,
+                    MPI_SUM, MPI_COMM_WORLD);
+      MPI_Allreduce(&local_size_ll, &local_max, 1, MPI_LONG_LONG,
+                    MPI_MAX, MPI_COMM_WORLD);
+      MPI_Allreduce(&local_size_ll, &local_min, 1, MPI_LONG_LONG,
+                    MPI_MIN, MPI_COMM_WORLD);
       if (rank == 0)
       {
          /* existing rank-0 print */
       }
+      // Per-rank machine-readable line for log-grep parity with R8-004.
+      std::cout << "[mixed-flux] rank=" << rank
+                << "  local_set_size=" << local_size_ll
+                << "  global_sum=" << global_sum
+                << "  global_min=" << local_min
+                << "  global_max=" << local_max << "\n";
```

**Test case:**

```cpp
// In a multi-rank smoke test:
mpirun -np 2 ./seas_tpv104_driver --dry-run --mixed-flux adjacent ...
// Assert both rank-0 stdout AND rank-1 stdout contain "[mixed-flux] rank="
// lines with the same global_sum value.
```

---

### [R-1410] [LOW] [QUALITY] [dynamic/godunov_flux.cpp:382-423] — `Central()` re-builds rotation matrices and frame on every call (same as Interior); for a face used by both ADER predictor (volume) and corrector (face), the rotation-pipeline overhead is doubled vs caching per-face frames

**Category:** QUALITY (performance — minor; pre-existing concern shared with Interior)

**Description:**
`Central()` mirrors `Interior()`'s 5-step rotation pipeline:
1. `BuildFrame(nor, t1, t2)`
2. `BuildRotationInverse(nor, t1, t2, Tinv)`
3. `BuildRotation(nor, t1, t2, T)`
4. Rotate states
5. Apply rotated-frame inner step
6. Rotate back

Steps 1-3 (frame + 2x 9x9 rotation matrices) are the dominant cost per call.
Each `BuildRotationInverse`/`BuildRotation` does ~120 multiplications. At
production (~1.5M faces × 4 RK4 stages × O(1000) Mult calls per macro step
× O(7000) macro steps), this is ~3e13 9x9 builds — but they are
**deterministic functions of `nor`** (face geometry doesn't change), so
they could be cached per-face.

This is a pre-existing concern that applies to `Interior()` too, but
mixed-flux DOUBLES the cost where it's engaged: at fault-adjacent faces,
the dispatch chooses between calling `Central` or `Interior`, and either
choice re-computes the rotation pipeline.

A trivial optimization: pass T, Tinv as arguments (or a per-face cache)
into both Interior and Central, avoiding redundant builds. Estimated
~10% speedup on flux assembly.

**Suggested fix:** Out of scope for the mixed-flux plan; flag as a
performance hardening for a future commit.

**Test case:** N/A (perf-only).

---

### [R-1411] [LOW] [DOC] [dynamic/wave_operator.inl:3667-3691, 4153-4192] — The R-1203 abort comment in ADER local + shared sites references `SEAS_TEST_NONFAULT_BOTH_SYM` symmetrization, but the surrounding `if (!nonfault_both_sym)` branches still exist and pretend to support it; once the env var is permanently disabled, the dead branch should be removed for clarity

**Category:** QUALITY (dead code after R-1203 fix)

**Description:**
Round-3 R-1203 added an `MFEM_ABORT` if `SEAS_TEST_NONFAULT_BOTH_SYM` is
set, effectively disabling the symmetrization hook. But the surrounding
code at L3667-3691 (ADER local) and L4177-4192 (ADER shared) still has:

```cpp
if (!nonfault_both_sym) {
   interior_or_central(nor, I_self, I_nbr, F_h);
} else {
   real_t F_pos[NUM_STATE], F_neg[NUM_STATE];
   ...
}
```

`nonfault_both_sym` is declared `const bool nonfault_both_sym = false;`
unconditionally (after R-1203). So the `else` branch is dead code. It
remains as documentation of the (now-disabled) symmetrization, but the
`if/else` adds compile-time complexity and a future reader may try to
re-enable it without reading the abort.

**Suggested fix (code diff):**

```diff
@@ wave_operator.inl L3667-3691 (ADER local non-fault dispatch)
-   if (!nonfault_both_sym)
-   {
-      interior_or_central(nor, I_self, I_nbr, F_h);
-   }
-   else
-   {
-      real_t F_pos[NUM_STATE], F_neg[NUM_STATE];
-      real_t nor_neg[3] = {-nor[0], -nor[1], -nor[2]};
-      interior_or_central(nor,     I_self, I_nbr, F_pos);
-      interior_or_central(nor_neg, I_nbr,  I_self, F_neg);
-      for (int c = 0; c < NUM_STATE; c++)
-      {
-         F_h[c] = 0.5 * (F_pos[c] + F_neg[c]);
-      }
-   }
+   // R-1203: nonfault_both_sym hook permanently disabled; the
+   // (n,L,R)→(-n,R,L) symmetrization is identically zero for any
+   // conservative flux and would zero out the non-fault flux.
+   // R-1411: dead else-branch removed.  See the abort at the env-var
+   // check earlier in this function.
+   interior_or_central(nor, I_self, I_nbr, F_h);
```

(Mirror at the shared site.)

**Test case:** N/A (style-only). The R-1203 abort still fires on env var.

---

### [R-1412] [LOW] [QUALITY] [tests/parallel/test_mixed_flux_dispatch_adjacent_mpi.cpp:301-310] — Test silently `return 0` when not run with np=2; no negative test ensures the test is actually run at np=2 in CI

**Category:** QUALITY (test infrastructure — silent skip)

**Description:**
At `test_mixed_flux_dispatch_adjacent_mpi.cpp:301-310`:

```cpp
if (nprocs != 2)
{
   if (rank == 0)
   {
      std::cerr << "[R-1201] requires np=2, got " << nprocs
                << " — skipping.\n";
   }
   MPI_Finalize();
   return 0;   // SUCCESS even though nothing was tested
}
```

Returning 0 silently when `nprocs != 2` means a CI misconfiguration that
runs the test at np=1 reports SUCCESS without actually exercising the
shared-face dispatch code path. The test's value is ZERO at np=1 but it
PASSES.

A misguided "speedup" optimization (e.g., changing `mpirun -np 2` to
`./seas_test_mixed_flux_dispatch_adjacent_mpi` to skip MPI overhead) would
silently disable cross-rank coverage for years.

**Trigger:** CI misconfig or developer running tests directly without mpirun.

**Suggested fix (code diff):**

```diff
@@ test_mixed_flux_dispatch_adjacent_mpi.cpp L301-310
    if (nprocs != 2)
    {
       if (rank == 0)
       {
          std::cerr << "[R-1201] requires np=2, got " << nprocs
-                   << " — skipping.\n";
+                   << " — TEST FAILED.  This test is meaningless at "
+                      "np!=2; CI must invoke via mpirun -np 2.\n";
       }
       MPI_Finalize();
-      return 0;
+      return 1;  // R-1412: misconfig is a TEST FAILURE, not a skip.
    }
```

And in the Makefile:

```diff
seas_test_mixed_flux_dispatch_adjacent_mpi: seas_test_mixed_flux_dispatch_adjacent_mpi.exe
-    ./seas_test_mixed_flux_dispatch_adjacent_mpi.exe
+    mpirun -np 2 ./seas_test_mixed_flux_dispatch_adjacent_mpi.exe
```

**Test case:** Run `seas_test_mixed_flux_dispatch_adjacent_mpi` directly
(without mpirun): pre-fix returns 0 (false PASS); post-fix returns 1 (FAIL).

---

## Summary

- Critical issues: **3** (R-1400, R-1401, R-1402)
  - **R-1400**: ADER shared-face Central dispatch at `wave_operator.inl:4165`
    has zero MPI numerical-output coverage. The Round-3 R-1201 fix only
    exercises the RK4 path (Mult); production uses ADER (AdvanceADER). This
    is the cross-rank conservation gap that R-1201 was supposed to close.
  - **R-1401**: `is_fault_face` and the dispatch-time fault check use
    different criteria (`fault_mesh_face_idx_set` vs `face_bdr_attr_`); a
    future ctor refactor that decouples them would silently double-dispatch
    at the fault face.
  - **R-1402**: `MixedFluxMode::AllContinuous` is exposed by the driver
    (`--mixed-flux all-continuous`) but has no end-to-end Mult/AdvanceADER
    test — only set-membership coverage. Zhang's published Mixed-Flux 1
    variant is unverified despite being a user-selectable mode.
- Moderate issues: **6** (R-1403, R-1404, R-1405, R-1406, R-1407, R-1408,
  R-1409)
  - **R-1403**: CFL not reduced for mixed-flux despite Zhang 2023 §3.3
    using CFL=0.3 for mixed-flux vs 0.5 for upwind. Production runs at
    `--cfl 0.5 --mixed-flux adjacent` may be near or beyond stability.
  - **R-1404**: `is_interior_face` excludes only fault-attr boundary
    faces; a non-fault non-zero `bdr_attr` on a 2-sided face would be
    inserted into the central set and fall through to mixed-flux dispatch
    instead of the BC branch.
  - **R-1405**: Gate 1 zero-jump test uses an absolute tolerance (1e-3)
    that is "1e-14 relative on stress channels" but only "1e-7 relative
    on velocity channels" — sub-channel bugs slip through.
  - **R-1406**: 2-tet sub-gate's comment cites the wrong reason for
    `S2.empty()`; the test does catch fault-misclassification regressions
    but doesn't exercise AllContinuous (which would catch them too on a
    different code path).
  - **R-1407**: `SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1` override is a silent
    bypass with no runtime warning, no banner annotation, and no audit
    trail; a user with the env var set in their shell rc loses the safety
    net.
  - **R-1408**: `central_flux_face_set_` is cached at SetMixedFluxMode
    time; the contract that "mesh and fault lists must not change after"
    is undocumented and unenforced.
  - **R-1409**: Driver banner uses `MPI_Reduce` (rank 0 only); per-rank
    log diagnosis of mixed-flux engagement requires `MPI_Allreduce`.
- Low issues: **3** (R-1410, R-1411, R-1412)
  - **R-1410**: `Central()` rebuilds the rotation pipeline per-call (same
    as `Interior()` — pre-existing perf hardening opportunity).
  - **R-1411**: Dead `else if (nonfault_both_sym)` branches at ADER local
    + shared sites; once R-1203 abort makes them unreachable, they
    should be removed.
  - **R-1412**: `test_mixed_flux_dispatch_adjacent_mpi` returns 0 (PASS)
    when run at np!=2, hiding CI misconfigurations.

**Verdict (round 4):** **CONDITIONAL FAIL — R-1400 must be fixed before
multi-rank ADER production.** R-1401, R-1402, R-1403 are CRITICAL and
should be addressed in the same milestone:

- **R-1400** is the load-bearing gap. Production TPV104 uses AdvanceADER,
  not Mult. The R-1201 fix verified only the RK4 shared dispatch.
- **R-1402** exposes a published mode (Mixed-Flux 1) with zero end-to-end
  test coverage despite the driver routing users to it.
- **R-1403** is a stability concern that may surface as "rupture front
  diverges" in long production runs.

Mitigations for moderate findings (R-1404 through R-1409) should land
within the same review cycle as the criticals.

## Recommended Next Steps

1. **R-1400** — Extend `test_mixed_flux_dispatch_adjacent_mpi.cpp` to call
   `wave.AdvanceADER` in addition to `wave.Mult`, and verify the same
   per-component invariants for the ADER path. This closes the actual
   R-1201 contract.
2. **R-1402** — Add `tests/unit/test_mixed_flux_dispatch_all_continuous.cpp`
   following the Adjacent test's structure, with an analytic identity gate
   covering the (much larger) AllContinuous central set.
3. **R-1403** — Reduce effective CFL by 0.6× when mixed-flux is engaged;
   document the recommendation in the driver banner. Add a stability
   smoke test.
4. **R-1401** — Apply the belt-and-suspenders predicate in `is_fault_face`.
5. **R-1404 through R-1409** — Apply the diffs above; each is a small,
   self-contained patch.

After these land, `--mixed-flux adjacent` is ready for a Frontera dry-run
sbatch with `mpirun -np 4` and a full-day soak. `--mixed-flux all-continuous`
also needs the R-1402 dispatch test to pass before any production use.

## Unreviewed Areas (continued from R3)

Still unreviewed in this round:
- **AdvanceADER + substep iterator + mixed-flux** combination
  (Plan §Risk R5 — no test exists at any np).
- **ParaView output coupling at Adjacent or AllContinuous**: does the
  station/face output read post-dispatch state, or does it cache stale
  upwind state? — not investigated in R4.
- **PML-region central faces**: a face at the PML-physical boundary that
  is fault-adjacent in the central set. PML damping is applied AFTER face
  flux; their interaction is not numerically gated.
- **Hybrid order (mixed P=1 + P=2 elements)** with mixed flux — TPV104
  is uniform-order so this is hypothetical, but the central-flux face set
  algorithm assumes uniform face DOF count.
