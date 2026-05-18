# Code Review (Round 5): Post-round-4-fix audit (2026-05-18)

## Review Scope

- **Plan:** `safs/project_7.0_alternative/document/spatial_dynamic_rupture_plan.md` (rev-3)
- **Prior reviews / fix reports consumed (rounds 1-4):**
  - `debug_document/spatial_dynamic_rupture_review_2026-05-18.md` (round 1)
  - `debug_document/spatial_dynamic_rupture_fix_2026-05-18.md`
  - `debug_document/spatial_dynamic_rupture_review_round2_2026-05-18.md`
  - `debug_document/spatial_dynamic_rupture_fix_round2_2026-05-18.md`
  - `debug_document/spatial_dynamic_rupture_review_round3_2026-05-18.md`
  - `debug_document/spatial_dynamic_rupture_fix_round3_2026-05-18.md`
  - `debug_document/spatial_dynamic_rupture_review_round4_2026-05-18.md`
  - `debug_document/spatial_dynamic_rupture_fix_round4_2026-05-18.md`
- **Files reviewed (post-round-4-fix state):**
  - `miniapps/seas/dynamic/wave_operator.hpp` (R-401 SetTime override + flag; R-407 abort message)
  - `miniapps/seas/dynamic/wave_operator.inl` (R-401 dispatch guards in interior + shared arms)
  - `miniapps/seas/dynamic/spatial_setup.hpp` (R-405 helper extraction)
  - `miniapps/seas/dynamic/godunov_flux_pool.cpp` (R-406 MFEM_VERIFY guard)
  - `miniapps/seas/fault/fault_geometry.hpp` (R-402 eager NaN-fill; R-404 docstring)
  - `miniapps/seas/tests/unit/test_phaseh_lsw_forced_rupture.cpp` (F-8 added)
  - `miniapps/seas/tests/unit/test_spatial_setup.cpp` (S-9 added)
  - `miniapps/seas/Makefile` (link wave_operator.o + deps into phaseh-lsw test)
- **Domain context consulted:**
  - `mfem/linalg/operator.hpp` (TimeDependentOperator::SetTime contract; `[[noreturn]]` on mfem_error)
  - `mfem/linalg/vector.cpp::Min,Max` (fmin/fmax NaN propagation semantics)
  - `mfem/general/error.hpp` (MFEM_ABORT / MFEM_VERIFY definitions)
  - `miniapps/seas/dynamic/seas_dynamic_operator.hpp::SetTime` (fan-out chain to WaveOperator)
- **Round-4 fix sweep:** 558/558 spatial + Phase H + Phase 5 tests pass (per round-4 fix report); F-8 verifies the new flag flips; S-9 verifies Print no longer crashes.

## Findings

---

### [R-501] [MODERATE] [test_phaseh_lsw_forced_rupture.cpp:F-8 vs wave_operator.inl dispatch arms] — Neither the interior nor shared-fault `MFEM_VERIFY` from R-401 is exercised by any unit test

**Category:** EDGE_CASE / QUALITY (regression risk)

**Description:**
R-401 added two `MFEM_VERIFY(time_was_set_ || fdata.T_forced_rupture >= 1.0e8, ...)` guards in `wave_operator.inl` (lines ~3642 and ~4594 — the interior fault path and the shared-fault fallback).  The intent is to catch a driver that wires `LSW_ForcedRupture` dispatch without calling `wave.SetTime(t)`.

F-8 only exercises the **flag mechanism** — it constructs a `WaveOperator`, asserts `TimeWasSet() == false`, calls `SetTime(0.0)`, asserts `TimeWasSet() == true`.  It does NOT exercise the dispatch arm at all; the WaveOperator is built with `bc.fault_attr = 0` (no fault) so the dispatch code path is never reached during Mult.

Concretely:
- If a future commit removes the `MFEM_VERIFY(time_was_set_ ...)` guards from wave_operator.inl (because they "fire too often" or "are noisy"), no test catches the regression — R-401's silent-forced-rupture failure mode returns.
- If a future commit changes the threshold from `>= 1.0e8` to `>= 1.0e9` (or any other off-by-one), no test catches it either.
- If a future refactor breaks the conjunction logic (e.g., swaps `||` for `&&`), no test catches it.

The guard is a safety mechanism whose value depends on its activation rate; without a test that fires the guard, we cannot prove it is wired in the right place.

**Trigger:**
Any of the above hypothetical regressions.

**Actual behavior:**
The guard exists but is dead code from the test suite's perspective.

**Expected behavior:**
A unit test that constructs a `WaveOperator` with at least one fault DOF, sets `LSW_ForcedRupture`, constructs a `DOFData` with `T_forced_rupture = 0` (active forced rupture), runs `wave.Mult(Q, dQdt)` WITHOUT calling `SetTime` first, and asserts the call aborts via the guard.

**Suggested fix:**
Add a regression test that exercises the guard.  The minimum viable shape uses a single-element mesh with one fault face:

```cpp
// tests/unit/test_phaseh_lsw_forced_rupture.cpp — add F-9
static void F_9_dispatch_guard_aborts_without_settime()
{
   // Build a minimal WaveOperator + FaultFaceFlux + FaultGeometry with
   // 1 fault face and 1 fault DOF.  Configure dof_data[0].T_forced_rupture
   // = 0 (forced rupture active).  Set FaultFrictionLaw::LSW_ForcedRupture.
   //
   // Run Mult inside a forked child — without calling SetTime first.
   // The MFEM_VERIFY in wave_operator.inl must fire; child aborts.
   const bool aborted = RunInChild([]() {
      // ... mesh + WaveOperator + fault setup ...
      DOFData& d = dof_data[0];
      d.T_forced_rupture = 0.0;     // active forced rupture
      d.t0_decay_forced  = 0.5;
      wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW_ForcedRupture);
      // DELIBERATELY skip wave.SetTime(t).
      Vector Q(wave.Height()); Q = 0.0;
      Vector dQdt;
      wave.Mult(Q, dQdt);     // guard MFEM_VERIFY must fire here
   });
   TEST_ASSERT(aborted, "guard must abort when SetTime() was never called");
}
```

Full setup (fault face, fault basis, DOFData wiring) is non-trivial scaffolding — but it pairs with the existing `RunInChild` pattern in this file, so the integration risk is low.  If the scaffolding is too heavy, the alternative is to add a minimal `WaveOperator` test fixture in a shared header that `test_phaseh_lsw_forced_rupture` and `test_spatial_setup` both consume.

**Test case:**
(See F-9 sketch above — the test IS the case.)

---

### [R-502] [LOW] [wave_operator.inl:3642 / 4594 — per-DOF MFEM_VERIFY hoisting] — Guard fires per-DOF inside the ADER hot loop; could hoist to per-Mult

**Category:** QUALITY (performance)

**Description:**
The R-401 guard:

```cpp
MFEM_VERIFY(time_was_set_ || fdata.T_forced_rupture >= 1.0e8, ...);
```

Lives inside the per-DOF dispatch block in `ComputeADERFaceFluxRHS` (and the shared-face equivalent).  For a SAFS production run with ~10⁵ fault DOFs and ~10⁴ macro steps, the guard runs ~10⁹ times per simulation.  MFEM_VERIFY's branch is cheap (single bool compare), but the per-DOF granularity is finer than needed: `time_was_set_` is a per-WaveOperator field (one value per Mult call), and `T_forced_rupture >= 1.0e8` could be checked once across the whole `dof_data_` vector at SetTime / SetFaultFrictionLaw time instead of per DOF per step.

**Trigger:**
Any SAFS dynamic-rupture production run with `LSW_ForcedRupture`.

**Suggested fix:**
Hoist the `time_was_set_` portion of the check out of the inner loop.  Inside `ComputeADERFaceFluxRHS`, compute once at the top:

```cpp
const bool time_check_needed = (fault_friction_law_ ==
                                 FaultFrictionLaw::LSW_ForcedRupture)
                                && !time_was_set_;
```

Then inside the dispatch arm, only check the per-DOF `T_forced_rupture` if `time_check_needed`:

```cpp
else if (fault_friction_law_ == FaultFrictionLaw::LSW_ForcedRupture)
{
   if (time_check_needed)
   {
      MFEM_VERIFY(fdata.T_forced_rupture >= 1.0e8, ...);
   }
   fault_flux_->EvaluateADER_LSW_ForcedRupture(...);
}
```

This drops the guard cost to a single check per non-forced-rupture DOF when forced rupture is active and SetTime was missed.

Note: the cleaner alternative — checking once that `any DOF has T_forced_rupture < 1e8` at SetFaultFrictionLaw or SetTime time — would require iterating `dof_data_` which the wave operator doesn't directly own.

---

### [R-503] [LOW] [fault_geometry.hpp:Print — eta_values_(0) for new BP5 ctor] — Print emits the BP5 SCALAR default eta, not the SAFS per-DOF eta

**Category:** QUALITY (misleading diagnostic)

**Description:**
The round-4 R-402 fix populates `eta_values_(i) = bp5_params_.eta()` for all i in the new BP5 ctor.  For SAFS curvilinear-fault runs with a heterogeneous `MaterialField`, the actual per-DOF eta is computed inside `seed_static_dof_fields` and lives in `DOFData::eta_p / eta_s` (per-DOF, heterogeneous).  `FaultGeometry::eta_values_` is a BP5-legacy SCALAR field.

`Print()` outputs:

```cpp
os << "  eta: " << eta_values_(0) / 1e6 << " MPa·s/m\n";
```

For a SAFS run, this shows `bp5_params_.eta()` — typically the BP5 default `(2670 · 3464²)/(2 · 3464) = 4.62 MPa·s/m`.  But the user's SAFS material's per-DOF eta could span 1 – 20 MPa·s/m depending on local Vs.  The diagnostic is misleading.

**Trigger:**
SAFS driver calls `geom.Print(...)` for startup logging.

**Actual behavior:**
Print shows BP5 default eta as if it were the SAFS eta.

**Expected behavior:**
Either (a) note in the Print output that the SAFS new-ctor path doesn't use a single eta scalar (eta is per-DOF in DOFData), or (b) emit NaN for eta_values_ on the new ctor (matching the loud-fail pattern for a/dc/V_init) so the user understands the value is not meaningful.

**Suggested fix:**
Option (a) — emit a clarifying line when on the new ctor.  Hard to detect that condition cleanly, so option (b) is simpler:

```diff
@@ fault_geometry.hpp new BP5 ctor body @@
   const real_t eta_const = bp5_params_.eta();
   for (int i = 0; i < N; ++i)
   {
      a_values_(i)            = k_nan;
-     eta_values_(i)          = eta_const;
+     // SAFS per-DOF eta lives in DOFData::eta_p/eta_s, not here.
+     // NaN-fill to fail loud if any caller assumes this is real.
+     eta_values_(i)          = k_nan;
      dc_values_(i)           = k_nan;
      ...
   }
```

Then `Print()`'s `eta_values_(0)` line emits `eta: nan MPa·s/m` for SAFS, which is observably wrong-as-expected.

Counter-argument: removing the `eta_const` value could surprise BP5 tests that re-use the new ctor for diagnostics.  However, the only BP5 tests in scope use the LEGACY ctor (which still computes eta_values_ correctly via ComputeBP5Params), so this is safe.

**Test case:**
(Diagnostic-only; the existing S-9 already proves Print doesn't crash.)

---

### [R-504] [LOW] [spatial_setup.hpp — IP-aware InitializeFaultDOFs_Spatial validator duplication] — Validator block duplicates the centroid overload's checks

**Category:** QUALITY (drift risk)

**Description:**
After R-405 extracted the per-DOF body into `copy_lsw_and_forced_rupture_fields`, the two `InitializeFaultDOFs_Spatial` overloads still share an ~12-line validator block at the top:

```cpp
MFEM_VERIFY(ndof >= 0, ...);
MFEM_VERIFY(dof_to_elem.Size() == ndof, ...);
MFEM_VERIFY(lsw.mu_s.Size() == ndof, ...);
MFEM_VERIFY(lsw.mu_d.Size() == ndof, ...);
MFEM_VERIFY(lsw.d_c.Size()  == ndof, ...);
MFEM_VERIFY(tau_pre.Size()      == 2 * ndof, ...);
MFEM_VERIFY(sigma_n_eff.Size()  == ndof, ...);
MFEM_VERIFY(T_forced_s.Size()   == ndof, ...);
MFEM_VERIFY(t0_decay_s.Size()   == ndof, ...);
```

A future addition (e.g., new per-DOF parameter) requires synchronised edits to both overloads.  Easy to drift.

**Suggested fix:**
Extract a `validate_per_dof_arrays` helper in the `internal` namespace:

```diff
@@ spatial_setup.hpp internal namespace @@
+inline void validate_lsw_per_dof_arrays(
+   int ndof,
+   const Array<int>& dof_to_elem,
+   const SlipWeakeningPerDOFParams& lsw,
+   const Vector& tau_pre,
+   const Vector& sigma_n_eff,
+   const Vector& T_forced_s,
+   const Vector& t0_decay_s)
+{
+   MFEM_VERIFY(ndof >= 0, ...);
+   MFEM_VERIFY(dof_to_elem.Size() == ndof, ...);
+   ... (the nine MFEM_VERIFYs)
+}
```

Both overloads then call it.  This is a refactor — confirm against the broader project policy that prefers verbose-explicit MFEM_VERIFY over helper extraction.

---

## Summary

- Critical issues: **0**
- Moderate issues: **1**  (R-501: dispatch-arm guard not unit-tested — regression risk)
- Low issues: **3**  (R-502 per-DOF guard hoist; R-503 Print eta misleading; R-504 validator duplication)
- Plan compliance: **PARTIAL** — Round-4 fixed all 7 round-4 findings as documented.  Plan phases H/4/6 still partially unimplemented per round-3 R-301/R-305/R-306/R-314/R-315 (acknowledged informational).
- Verdict: **PASS WITH FIXES** — no remaining CRITICAL bugs.  R-501 is the only MODERATE: a regression-risk concern (guard not directly tested) that the round-4 fix acknowledged could not be exercised without a real fault-DOF fixture; the gap should be closed before Phase 4 driver lands so the guard becomes load-bearing.  R-502/R-503/R-504 are quality polish that does not affect correctness.

Round-3's R-301 (flux pool dispatch dead) and R-302 (LSW_ForcedRupture dispatch arm) are still PARTIALLY addressed: R-302 is wired but the production driver path requires R-401's flag-flip discipline.  Phase 4 driver is the missing piece — once it lands and calls `wave.SetTime(t)` per macro step, the round-3/round-4 fixes will all combine into a working forced-rupture pipeline.

## Unreviewed Areas

- **Heterogeneous WaveOperator code** — Phase H.2/H.3/H.5 still unimplemented (R-301 abort).  Nothing new to review.
- **Phase 4 driver source** — still not present.  When it lands, every R-401 guard fix becomes load-bearing; a test like R-501's F-9 will catch regressions there.
- **Phase 6 sbatches + verification scripts** — only `safs_smoke_8N_400r_dev.sbatch` and `verify_constant_tensor_projection.py` exist; the plan-mandated `verify_spatial_dyn_smoke_safs.py` and `compare_spatial_dyn_velocity_models.py` remain missing.
- **MPI behaviour of the new BP5 ctor + dispatch guards** — round-4 fix did not run a `mpirun -np 4` regression sweep; `ComputeGatherInfo()` in the new ctor and `time_was_set_` propagation across ranks are untested in parallel.  The shared-fault fallback dispatch path uses the same `time_was_set_` flag (per-rank), which should naturally be consistent if each rank calls `SetTime(t)` at the same macro step — but a parallel regression would confirm.
- **TPV byte-exact full regression** — only `seas_test_tpv104_checkpoint` (164/164) and `seas_test_compute_safs_params` (13/13) were run in round 4; the full `make test-tpv104` / `make test-tpv205` / `make test-bp5-*` byte-exact gates listed in the plan §TPV Benchmark Isolation were not re-run.  Recommended before Phase H wiring lands or any commit on `wave_operator.inl` merges.
- **Round 1-3 carry-over findings** — confirmed passing in round-4's 558/558 sweep; not re-audited in round 5.
