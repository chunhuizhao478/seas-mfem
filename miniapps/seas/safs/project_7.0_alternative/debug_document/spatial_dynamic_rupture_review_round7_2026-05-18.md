# Code Review (Round 7): Post-round-6 audit + NucleationKind generalisation (2026-05-18)

## Review Scope

- **Plan:** `safs/project_7.0_alternative/document/spatial_dynamic_rupture_plan.md` (rev-3)
- **Prior reviews / fix reports consumed (rounds 1-6):** see `debug_document/spatial_dynamic_rupture_review_round{2..6}_2026-05-18.md` and `..._fix_round{1..6}_2026-05-18.md`.
- **Files reviewed (round-6 changes):**
  - `miniapps/seas/dynamic/tpv205_substep_iterator.{hpp,cpp}` (R-601 SetForcedRuptureMode + t_sub_abs threading)
  - `miniapps/seas/drivers/spatial_dyn_driver.cpp` (R-602 paraview mu_eff; R-603 last_completed_step; R-604 ref_normal cleanup; R-607 #ifdef cleanup; NucleationKind dispatch)
  - `miniapps/seas/spatial/code/spatial_friction.{hpp,cpp}` (NucleationKind enum, OverstressSpec, ResolveOverstress stub, parser support)
  - `miniapps/seas/tests/unit/test_phaseh_lsw_forced_rupture.cpp` (F-10 iterator regression)
  - `miniapps/seas/tests/unit/test_phaseh_wave_operator_constant_parity.cpp` (round-6 added by implementer; reviewed end-to-end)
- **Domain context:**
  - `miniapps/seas/CLAUDE.md`
  - `miniapps/seas/dynamic/tpv205_friction.hpp::LSWFrictionCoefficient_TPV205`
  - `miniapps/seas/spatial/code/spatial_friction.hpp::LSWFrictionCoefficient_ForcedRupture`
  - `miniapps/seas/dynamic/wave_operator.inl::ComputeADER*FluxRHS` (dispatch + substep gate)
- **Round-6 test sweep:** 583 / 583 pass; spatial_dyn_driver builds clean; TPV104/205 byte-exact regression intact (164/164).

## Findings

---

### [R-701] [CRITICAL] [spatial_friction.cpp:resolve_forced_impl + spatial_dyn_driver.cpp:1188 SetForcedRuptureMode] — When `nucleation.kind == Overstress`, `ResolveForcedRupture` still populates strength-reduction `T_forced(r)` values, and the iterator's `SetForcedRuptureMode(true)` is UNCONDITIONAL; once the overstress stub fills in, an overstress-mode run will silently double-apply strength reduction

**Category:** BUG (latent until ResolveOverstress stub fills in)

**Description:**
The round-6 NucleationKind dispatch adds a `kind` switch to `[nucleation]` so users can pick `strength_reduction` (TPV26/27, the default) OR `overstress` (TPV205, currently stubbed).  The driver wires the dispatch tag based on `kind`:

```cpp
// drivers/spatial_dyn_driver.cpp:815-819
const bool use_strength_reduction =
   (cfg.nucleation.kind == spatial::NucleationKind::StrengthReduction);
wave.SetFaultFrictionLaw(use_strength_reduction
                         ? FaultFrictionLaw::LSW_ForcedRupture
                         : FaultFrictionLaw::LSW);
```

But two downstream pieces are NOT gated on `kind`:

1. **`ResolveForcedRupture` ignores `nuc.kind`** (`spatial_friction.cpp:1103-1109`):

   ```cpp
   if (!nuc.enabled)
   {
      p.T_forced_s = 1.0e9;
      p.t0_decay_s = 0.0;
      return p;
   }
   ```

   The early-out triggers ONLY on `!nuc.enabled`.  When `nuc.enabled = true` AND `nuc.kind = Overstress`, the body proceeds to compute strength-reduction `T_forced(r)` per hypocenter distance — exactly what overstress mode does NOT want.

2. **`SetForcedRuptureMode(true)` is unconditional** (`drivers/spatial_dyn_driver.cpp:1188`):

   ```cpp
   Tpv205SubStepIterator substep_iterator(fault_flux);
   substep_iterator.SetForcedRuptureMode(true);     // ALWAYS true
   ```

   So even when the dispatch tag is `LSW` (overstress path), the iterator's per-substep `mu_eff` calls `LSWFrictionCoefficient_ForcedRupture` which READS `d.T_forced_rupture` and `d.t0_decay_forced` from `InitializeFaultDOFs_Spatial` output.  Those fields were populated from `ResolveForcedRupture`'s WRONG strength-reduction values (see #1).  Result: iterator applies strength reduction to a mode the user explicitly chose to be overstress-only.

The driver-comment at line 941-948 claims the resolver "writes the 'never forced' sentinel (T = 1e9 everywhere)" for overstress mode — **the resolver does NOT do this**; the comment is wrong.

This bug is LATENT today because `ResolveOverstress` aborts in its stub body, so a user setting `[nucleation].kind = "overstress"` hits the abort before any wrong physics fires.  But the moment the overstress stub fills in with real per-DOF `(Δτ, Δσ_n)` perturbations, the iterator's incorrect strength-reduction application becomes a silent wrong-physics bug — the WORST kind, because the dispatch tag was correctly set to `LSW` and the user thinks the run is pure-overstress.

**Trigger:**
TOML config with `[nucleation] enabled = true` AND `[nucleation].kind = "overstress"`.  Currently → ResolveOverstress aborts (safe).  Once the stub fills in → silent strength-reduction applied on top of overstress.

**Actual behavior (latent):**
Iterator computes mu_eff via `LSWFrictionCoefficient_ForcedRupture` with hypocenter-based T_forced, applying strength reduction in addition to the overstress perturbation — double-counts the nucleation mechanism.

**Expected behavior:**
- `ResolveForcedRupture` must write the 1e9 sentinel for any `nuc.kind != NucleationKind::StrengthReduction`, so the iterator's `T_forced_rupture` reads are no-ops in non-strength-reduction modes.
- Driver must gate `SetForcedRuptureMode(...)` on `use_strength_reduction` so the iterator uses plain LSW for overstress.

Both fixes are belt-and-suspenders — either alone prevents the bug.  Apply both for defence in depth.

**Suggested fix:**

Fix (a) — sentinel for non-StrengthReduction kinds:

```diff
@@ spatial_friction.cpp::resolve_forced_impl @@
-   if (!nuc.enabled)
+   if (!nuc.enabled
+       || nuc.kind != NucleationKind::StrengthReduction)
    {
-      // TPV205 byte-exact default: T = 1e9, t_0 = 0.
+      // TPV205 byte-exact default: T = 1e9, t_0 = 0.
+      // Also returned for non-StrengthReduction kinds (e.g. Overstress)
+      // so the iterator's forced-rupture path is a no-op in those
+      // modes — the overstress perturbation lives in
+      // DOFData::tau{1,2}_nuc / sigma_n_nuc and is consumed by plain
+      // LSW instead (driver dispatches FaultFrictionLaw::LSW for
+      // kind=Overstress).
       p.T_forced_s = 1.0e9;
       p.t0_decay_s = 0.0;
       return p;
    }
```

Fix (b) — driver gates iterator's forced-rupture mode:

```diff
@@ drivers/spatial_dyn_driver.cpp:1187-1188 @@
    Tpv205SubStepIterator substep_iterator(fault_flux);
-   substep_iterator.SetForcedRuptureMode(true);
+   // R-701: only enable the iterator's forced-rupture mu_eff path
+   // when strength-reduction is the chosen mechanism.  Overstress
+   // mode runs plain LSW; the per-DOF Δτ/Δσ_n perturbation lives in
+   // DOFData::tau{1,2}_nuc / sigma_n_nuc populated at init time.
+   substep_iterator.SetForcedRuptureMode(use_strength_reduction);
```

And update the now-incorrect driver comment at line 941-948 to match the fixed behavior:

```diff
@@ drivers/spatial_dyn_driver.cpp:941-948 @@
-   // -----------------------------------------------------------------
-   // 12. Per-DOF forced-rupture times (Phase 1 / D-4) — populated only
-   //     when the strength-reduction kind is selected.  For overstress
-   //     mode the resolver writes the "never forced" sentinel (T = 1e9
-   //     everywhere) so the iterator's forced-rupture mu_eff path is a
-   //     no-op even when LSW_ForcedRupture happens to be wired.
-   //     Round-6: ResolveOverstress runs in parallel (stub).
-   // -----------------------------------------------------------------
+   // -----------------------------------------------------------------
+   // 12. Per-DOF forced-rupture times (Phase 1 / D-4).  ResolveForced
+   //     Rupture writes T_forced(r) per hypocenter distance when
+   //     cfg.nucleation is enabled AND kind == StrengthReduction; for
+   //     any other kind it returns the 1e9 sentinel so the iterator's
+   //     forced-rupture mu_eff path is a no-op (the iterator is also
+   //     gated by SetForcedRuptureMode below, so this is belt-and-
+   //     suspenders).
+   // -----------------------------------------------------------------
```

**Test case:**
```cpp
// tests/unit/test_spatial_friction_resolver.cpp — add F-4
static void F_4_ResolveForcedRupture_overstress_kind_returns_sentinel()
{
   // Mirror F-2's hypocenter geometry but set nuc.kind = Overstress.
   // Pre-fix: T_forced(0) computed per hypocenter distance — non-1e9.
   // Post-fix: T_forced(i) == 1e9 for every DOF.
   NucleationSpec nuc;
   nuc.enabled        = true;
   nuc.kind           = NucleationKind::Overstress;  // NEW
   nuc.hypocenter_x_m = 0.0;  nuc.hypocenter_y_m = 0.0;  nuc.hypocenter_z_m = 0.0;
   nuc.r_crit_m       = 4000.0;
   nuc.t0_decay_s     = 0.5;

   // 1 DOF at the hypocenter.
   Vector dofs(3); dofs = 0.0;
   Array<int> elem(1);  elem = 0;
   auto mat = MaterialField::MakeConstant(32e9, 32e9, 2670.0);
   TinyMeshHolder mh;

   SpatialFrictionResolver R;
   auto p = R.ResolveForcedRupture(nuc, dofs, elem, mat, mh.mesh());
   TEST_NEAR(p.T_forced_s(0), 1.0e9, 0.0,
             "Overstress kind returns 1e9 sentinel even when enabled");
}
```

---

### [R-702] [MODERATE] [tpv205_substep_iterator.cpp::Advance,AdvanceWithSubStepStates + wave_operator.inl shared-fault dispatch] — Interior-fault and shared-fault forced-rupture paths use DIFFERENT times when `ader_order >= 2`

**Category:** EDGE_CASE / ASSUMPTION

**Description:**
The round-6 fix threads `t_sub_acc = t_macro_start + Σ_{o'<o} deltaT[o']` through `StepOneQP_` so the iterator's per-substep `mu_eff` reads a current absolute time.  Good.

But the WAVE OPERATOR'S inline dispatch arm (used for shared-fault QPs, line 4787+ in wave_operator.inl) calls:

```cpp
fault_flux_->EvaluateADER_LSW_ForcedRupture(
   fdata, ..., dt, GetTime(), ...);
```

`GetTime()` returns the WaveOperator's stored `t`, which the driver sets ONCE per macro step (line 1218 `wave.SetTime(t)`).  So the shared-fault path always evaluates `f_2(t_macro_start)` — a SINGLE time per macro step.

For `ader_order = 1` (one substep), the two paths use the same time `t_macro_start` — consistent.

For `ader_order >= 2` (multiple substeps), the iterator evaluates `f_2` at `t_macro_start`, `t_macro_start + dt/O`, …, `t_macro_start + (O-1)·dt/O` per substep — but the shared-fault path uses only `t_macro_start`.  Interior and shared QPs at the same instant see different effective friction coefficients.  The discrepancy is bounded by `(O-1)/O · dt` per macro step (≤ dt total), so for a typical SAFS dt ≈ 0.01s it's at most ~10 ms of time-evaluation lag at the rupture front.

This is a **bounded** physical inconsistency, not a divergent bug — both paths converge to the correct answer as `dt → 0`.  But it's undocumented and could surface as a sharp "discontinuity" in the slip-rate field at MPI partition seams when the rupture front crosses a seam and ader_order > 1.

**Trigger:**
SAFS run with `ader_order ≥ 2` and rupture nucleating near an MPI partition seam.

**Actual behavior:**
Interior fault QPs see `f_2(t_sub_abs)`; shared fault QPs see `f_2(t_macro_start)`.

**Expected behavior:**
Either (a) document the inconsistency as a known limitation of the substep iterator vs inline dispatch, or (b) make the shared-fault dispatch use a per-substep time too.  (b) requires plumbing substep time into `ComputeADERSharedFaceFluxRHS` — bigger change.  (a) is the immediate fix.

**Suggested fix (option a — document):**

```diff
@@ wave_operator.inl: shared-fault LSW_ForcedRupture arm @@
+                  // R-702 (round-7): SHARED-FAULT TIME GRANULARITY.
+                  // GetTime() here is the macro-step start time.  The
+                  // iterator path (interior fault QPs) evaluates f_2(t)
+                  // at per-substep times t_macro + Σdeltaτ; this shared
+                  // arm uses only t_macro.  For ader_order >= 2 the
+                  // two paths see DIFFERENT effective friction
+                  // coefficients at the same physical instant.  The
+                  // mismatch is bounded by dt per macro step and
+                  // vanishes as dt → 0; it may surface as a small
+                  // discontinuity in slip-rate at MPI partition seams
+                  // when the rupture front crosses them under
+                  // ader_order > 1.  Future fix: thread substep time
+                  // into ComputeADERSharedFaceFluxRHS.
                   VerifyForcedRuptureTimeReady(
                      time_was_set_, fdata.T_forced_rupture);
                   fault_flux_->EvaluateADER_LSW_ForcedRupture(...);
```

(Same comment at the interior-fault arm.)

**Test case:**
```cpp
// Drive a SAFS-style synthetic run with ader_order=2 and a rupture
// front crossing an MPI partition seam.  Diff slip-rate fields on
// either side of the seam at peak rupture-front passage; assert the
// diff is bounded by the documented dt-scale.
// (Integration test; defer until SAFS verification scripts land.)
```

---

### [R-703] [MODERATE] [tests/unit/test_spatial_friction_*.cpp] — NucleationKind enum + parser + `ResolveOverstress` stub are NOT exercised by any unit test

**Category:** EDGE_CASE / QUALITY (coverage)

**Description:**
Round-6 added:
- `enum class NucleationKind { StrengthReduction, Overstress }`
- `[nucleation].kind = "strength_reduction" | "overstress"` parser path
- `[nucleation.overstress]` sub-block parser path
- `SpatialFrictionResolver::ResolveOverstress(...)` stub

`grep` for `NucleationKind`, `Overstress`, or `kind.*nucleation` in `tests/unit/test_spatial_friction_*.cpp` returns ZERO hits.  No test:
- Parses a TOML with `[nucleation].kind = "strength_reduction"` and asserts `cfg.nucleation.kind == NucleationKind::StrengthReduction`.
- Parses a TOML with `[nucleation].kind = "overstress"` and asserts the field round-trips.
- Parses a TOML with `[nucleation].kind = "garbage"` and asserts the parser aborts cleanly.
- Calls `ResolveOverstress(...)` with kind=StrengthReduction (expect: empty vectors).
- Calls `ResolveOverstress(...)` with kind=Overstress (expect: MFEM_ABORT).
- Parses a TOML with `[nucleation.overstress].direction = 2` (out of range) and asserts the parser aborts.

A typo in the parser (e.g., `"strenght_reduction"` typo in a future refactor) would not be caught.  A regression in `ResolveOverstress` (e.g., the abort getting removed by accident before the body is filled) would silently return empty vectors instead of failing loud.

**Trigger:**
Future commit modifying the parser or resolver — no test fires.

**Suggested fix:**
Add a `T_21_nucleation_kind_parser_*` block in `test_spatial_friction_config.cpp` and a `F_4_overstress_stub` block in `test_spatial_friction_resolver.cpp`:

```cpp
// test_spatial_friction_config.cpp — add T-21..T-23
static void T_21_nucleation_kind_default_is_strength_reduction()
{
   std::string toml = MinimalLSWHeader() + MinimalLSWBlock();
   toml += "[nucleation]\nhypocenter_x_m=0\nhypocenter_y_m=0\n"
           "hypocenter_z_m=0\nr_crit_m=4000\nt0_decay_s=0.5\n";
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.nucleation.kind == NucleationKind::StrengthReduction,
               "[nucleation] without `kind` defaults to StrengthReduction");
   TEST_ASSERT(cfg.nucleation.enabled, "[nucleation] block sets enabled=true");
}

static void T_22_nucleation_kind_overstress_parses()
{
   std::string toml = MinimalLSWHeader() + MinimalLSWBlock();
   toml += "[nucleation]\nhypocenter_x_m=0\nhypocenter_y_m=0\n"
           "hypocenter_z_m=0\nr_crit_m=4000\nt0_decay_s=0.5\n"
           "kind=\"overstress\"\n";
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.nucleation.kind == NucleationKind::Overstress,
               "[nucleation].kind=\"overstress\" round-trips");
}

static void T_23_nucleation_kind_garbage_aborts()
{
   std::string toml = MinimalLSWHeader() + MinimalLSWBlock();
   toml += "[nucleation]\nhypocenter_x_m=0\nhypocenter_y_m=0\n"
           "hypocenter_z_m=0\nr_crit_m=4000\nt0_decay_s=0.5\n"
           "kind=\"strenght_reduction\"\n";   // typo
   TEST_ASSERT(ParseAbortsInChild(toml),
               "[nucleation].kind with typo must abort");
}
```

```cpp
// test_spatial_friction_resolver.cpp — add F-4 / F-5
static void F_4_ResolveOverstress_StrengthReduction_returns_empty()
{
   NucleationSpec nuc;
   nuc.enabled = true;
   nuc.kind    = NucleationKind::StrengthReduction;
   Vector dofs(6); dofs = 0.0;
   SpatialFrictionResolver R;
   auto p = R.ResolveOverstress(nuc, dofs);
   TEST_ASSERT(p.delta_tau1_pa.Size() == 0,
               "ResolveOverstress on StrengthReduction returns empty");
}
static void F_5_ResolveOverstress_Overstress_aborts()
{
   const bool aborted = RunInChild([]() {
      NucleationSpec nuc;
      nuc.enabled = true;
      nuc.kind    = NucleationKind::Overstress;
      Vector dofs(3); dofs = 0.0;
      SpatialFrictionResolver R;
      (void)R.ResolveOverstress(nuc, dofs);
   });
   TEST_ASSERT(aborted, "stub ResolveOverstress aborts on Overstress kind");
}
```

---

### [R-704] [LOW] [drivers/spatial_dyn_driver.cpp:941-948 — incorrect comment] — Driver comment claims `ResolveForcedRupture` writes the sentinel for overstress mode; the resolver does NOT

**Category:** QUALITY (misleading documentation)

**Description:**
The current comment block at lines 941-948 says:

```cpp
// 12. Per-DOF forced-rupture times (Phase 1 / D-4) — populated only
//     when the strength-reduction kind is selected.  For overstress
//     mode the resolver writes the "never forced" sentinel (T = 1e9
//     everywhere) so the iterator's forced-rupture mu_eff path is a
//     no-op even when LSW_ForcedRupture happens to be wired.
```

But `ResolveForcedRupture` (per R-701) does NOT check `kind` and DOES populate strength-reduction times for `enabled && kind=Overstress`.  The comment lies about behavior — a reader who trusts it will assume R-701 doesn't exist.

**Suggested fix:**
Subsumed by R-701's suggested fix — update the comment to match the (fixed) behavior.

---

### [R-705] [LOW] [tpv205_substep_iterator.cpp:226 + 364 — `t_sub_acc = t_macro_start` start-of-substep convention] — Iterator evaluates `f_2` at the START of each substep; midpoint or Gauss–Lobatto would be more accurate for ADER

**Category:** EDGE_CASE / POSSIBLE

**Description:**
`Advance` and `AdvanceWithSubStepStates` both initialise `t_sub_acc = t_macro_start` and the first substep evaluates `f_2(t_macro_start)`.  For ADER quadrature with equal-weight subdivision (`dt/O` per substep), the more accurate friction evaluation points would be the midpoints `(t_macro_start + (o+0.5)·dt/O)` per substep.  Using start-of-substep introduces a half-step time-lag in the friction evaluation that compounds with `O`.

The error is bounded but matters for time-critical nucleation: at the hypocenter DOF with T_forced=0, `f_2(t_macro_start)` = 0 on the first substep even when the macro step physically straddles the start of forced rupture; midpoint would give the correct ramp value.

For O=1 substep, this is a half-step lag.  Documented in the F-10 test comment ("ADER-O1 with single midpoint"), but the implementation uses start-of-substep.  Minor inconsistency.

**Suggested fix:**
Either (a) document the start-of-substep choice (and the upper-bound time-lag), or (b) shift to midpoint:

```diff
@@ tpv205_substep_iterator.cpp::Advance @@
-   real_t t_sub_acc = t_macro_start;
+   real_t t_sub_acc = t_macro_start;   // tracks substep-START time
    for (int o = 0; o < O; ++o)
    {
       const real_t dt_sub      = deltaT_[o];
       ...
-      StepOneQP_(..., t_sub_acc, ...);
+      // R-705: midpoint sample for the friction time so f_2 sees the
+      // mean substep time, not the start.  For T_forced exactly inside
+      // [t_sub_acc, t_sub_acc + dt_sub) this avoids a half-step lag.
+      StepOneQP_(..., t_sub_acc + 0.5 * dt_sub, ...);
       ...
       t_sub_acc += dt_sub;
    }
```

POSSIBLE — depends on what the plan §Phase H.6 considers "current sim time".  The plan text isn't fully unambiguous.  Flag as POSSIBLE; downgrade to LOW.

---

### [R-706] [LOW] [drivers/spatial_dyn_driver.cpp:961-963 — `(void)overstress`] — Resolver result intentionally unused; loud `[[maybe_unused]]` would be clearer than `(void)cast`

**Category:** QUALITY (minor)

**Description:**
The driver consumes the `ResolveOverstress` result only as `(void)overstress;` to silence the unused-variable warning, with a comment explaining the consumer is deferred:

```cpp
const spatial::OverstressPerDOFParams overstress =
   resolver.ResolveOverstress(cfg.nucleation, dof_coords_3d);
(void)overstress;   // consumed in the (deferred) overstress init
```

The `(void)` cast works but is a C-style idiom.  C++17 `[[maybe_unused]]` is cleaner and self-documenting.  Minor style preference.

**Suggested fix:**
```diff
-   const spatial::OverstressPerDOFParams overstress =
+   [[maybe_unused]] const spatial::OverstressPerDOFParams overstress =
       resolver.ResolveOverstress(cfg.nucleation, dof_coords_3d);
-   (void)overstress;   // consumed in the (deferred) overstress init
```

---

## Summary

- Critical issues: **1**  (R-701: latent overstress double-application — surfaces when ResolveOverstress stub fills in)
- Moderate issues: **2**  (R-702 interior/shared time-granularity mismatch at ader_order ≥ 2; R-703 NucleationKind/Overstress not exercised by any test)
- Low issues: **3**  (R-704 misleading comment subsumed by R-701; R-705 substep-start vs midpoint convention; R-706 `(void)` cast style)
- Plan compliance: **PARTIAL** — Phase H Stage 1 fully shipped + tested (constant-parity gate green); Phase 4 driver shipped + builds clean; nucleation generalised to {strength_reduction, overstress} per user request; overstress resolver is a stub (deferred per the round-6 fix report).
- Verdict: **PASS WITH FIXES** — no remaining CRITICAL bugs on the strength-reduction production path (the round-6 fix correctly wired SetForcedRuptureMode and the F-10 regression test proves the iterator now consumes T_forced).  R-701 is the only CRITICAL and is LATENT (the stub abort masks it today).  Land R-701's two-line fix (ResolveForcedRupture sentinel for non-StrengthReduction kinds + driver gating SetForcedRuptureMode on use_strength_reduction) before the overstress stub fills in so the next implementer doesn't inherit a silent wrong-physics trap.  R-702 / R-703 are coverage and consistency concerns that should land in the next pass.

Cumulative trajectory 1→7: rounds 1-6 found 25 total findings (3+0+3+2+0+1 CRITICAL); round 7 finds 6 (1 latent CRITICAL).  The pattern is convergent — each fix round patches one layer and reveals issues one layer deeper; round 7's findings are about FUTURE behavior (overstress) and minor consistency, not active bugs in the strength-reduction production path.

## Unreviewed Areas

- **Overstress resolver body** — stub.  When implemented, R-701 + a paired set of integration tests will need to be in place.
- **MPI behaviour of the SAFS driver under np > 1** — the driver builds clean and `--config` smoke runs, but a `mpirun -np 4 seas_spatial_dyn_driver --config <fixture>.toml` end-to-end run was not invoked in this round.  The R-702 shared-fault time-granularity concern is the most relevant parallel issue.
- **TPV byte-exact full regression** (`make test-tpv104` / `make test-tpv205` / `make test-bp5-*`) — the iterator's signature change (`StepOneQP_` gains a `t_sub_abs` parameter) and the Tpv205SubStepIterator's behavioural addition (forced_rupture_mode_) should not affect TPV205 byte-exact runs (TPV205 driver does NOT call SetForcedRuptureMode, so the flag stays false and the original mu_eff path is preserved).  Should be confirmed with a re-run of the full TPV byte-exact suite before any commit that touches the iterator merges.
- **Round 1-6 carry-over findings** — assumed still resolved per the cumulative fix reports; not re-audited.
- **Phase H Stage 2** — per-element flux dispatch in `wave_operator.inl` hot loop, full bi-material MPI exchange — still not implemented (round-3 R-301 abort + round-5 R-407 guidance) — out of scope for this review.
