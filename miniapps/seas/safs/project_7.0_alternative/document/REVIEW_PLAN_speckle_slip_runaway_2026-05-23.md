# Code Review: PLAN_speckle_slip_runaway_2026-05-23.md (plan audit) — 2026-05-23

## Review Scope
- **Plan reviewed:** `safs/project_7.0_alternative/document/PLAN_speckle_slip_runaway_2026-05-23.md`
- **Companion:** `REVIEW_speckle_slip_runaway_debug_2026-05-23.md` (the debug review this plan
  derives from).
- **Code verified against:** `dynamic/tpv205_substep_iterator.cpp` (`StepOneQP_`,
  `AdvanceWithSubStepStates`, `Advance`), `dynamic/tpv205_friction.hpp` (`SolveLSW_TPV205`),
  `dynamic/fault_face_flux.{hpp,cpp}` (`DOFData`, `cs_`, `WriteBackState`),
  `drivers/spatial_dyn_driver.cpp` + `drivers/tpv205_driver.cpp` (iterator call sites, monitor),
  `dynamic/wave_operator.inl` (`ComputeADERSharedFaceFluxRHS` inline shared-QP solve).
- **Domain context:** CLAUDE.md (no-magic-constants, byte-exact contract, "report don't work
  around"), the reconcile review chain.

This is a **plan audit**: findings are defects in the plan's logic / code-level claims / proposed
snippets that would make the implementation wrong, byte-breaking, or ineffective if applied as
written. The plan is well-structured and its mechanism analysis is correct; the issues below are
in the proposed *fixes*, especially the central tripwire.

## Findings

### [R-001] CRITICAL [Phase 2, slip-plausibility tripwire, plan lines 169-172] — the tripwire is mathematically incapable of firing: it compares slip against (a multiple of) peak-velocity × time, a bound slip can never exceed

**Category:** BUG (the central new safety check does nothing)

**Description:**
The proposed bound is `slip_bound = slip_factor · V_substep_peak_global · t + slip_floor`, with
`slip_factor` default **4.0** and `V_substep_peak_global` the running max over sub-step `V_abs`.
But slip is the time-integral of velocity:
`slip_i = |∫₀ᵗ V_i ds| ≤ ∫₀ᵗ |V_i| ds ≤ (max_s |V_i(s)|)·t ≤ V_substep_peak_global · t`.
That last inequality is an **identity-level upper bound** — and slip is accumulated *only* in
`StepOneQP_` from `s.V1/s.V2` (`tpv205_substep_iterator.cpp:122-123`; the inline shared-QP solve
does NOT accumulate slip, `fault_face_flux.cpp:815-824`), so `slip_rate_substep_max` (hence
`V_substep_peak_global`) captures exactly the velocities that produced the slip. Therefore for any
`slip_factor ≥ 1`:
```
slip_i ≤ V_substep_peak_global · t ≤ slip_factor · V_substep_peak_global · t < slip_bound
```
**The tripwire can never trip** — including on the actual `max_slip = 1e5 m` runaway, because the
velocities that integrated to 1e5 m make `4·V_peak·t ≫ 1e5 m`. The check is comparing slip against
its own kinematic ceiling; it is self-defeating, not a plausibility test.

**Trigger:** the runaway it is designed to catch (or any input), with the default `slip_factor=4`.

**Actual behavior:** `n_implausible` is always 0; the tripwire is a no-op.

**Expected behavior:** flag slip that is inconsistent with a **physical** velocity ceiling, not
with the observed (already-pathological) peak.

**Suggested fix:** bound against a physical velocity ceiling decoupled from the observed peak —
reuse R-002's `c_s`-derived ceiling, so the test becomes "slip exceeds what a physically-bounded
slip rate could produce in the elapsed time":
```diff
- const real_t slip_factor  = .../*env, default 4.0*/;
- const real_t slip_bound   = slip_factor * V_substep_peak_global * t + slip_floor;
+ // Physical ceiling: a slip rate cannot exceed ~k·c_s for sustained time.
+ // Bounding against the OBSERVED peak is circular (slip ≤ V_peak·t identically),
+ // so use the same c_s ceiling as the R-002 V-sanity bound.
+ const real_t v_phys_ceiling = vsane_kcs * flux_cs;          // k·c_s, k~3
+ const real_t slip_bound     = v_phys_ceiling * t + slip_floor;
```
Now `slip_i ≤ V_peak·t` no longer protects a runaway whose `V_peak ≫ k·c_s`: the 1e5 m DOF
violates `slip_i ≤ k·c_s·t` (≈ 3·3000·1 ≈ 9e3 m) and trips.

**Test case:**
```cpp
// Identity check that the AS-WRITTEN tripwire never fires:
//   given any V-history, slip_i = sum(V*dt) and V_peak = max(V),
//   assert slip_i <= V_peak * t_total  (so slip_factor>=1 => bound never exceeded).
// Then the FIXED tripwire DOES fire on the runaway:
//   V history with V=2e4 m/s for 5 sub-steps of dt=1e-4 -> slip≈10 m at t=5e-4;
//   extrapolated to t=1 s with sustained over-slip -> slip=1e5 m;
//   assert slip_i (1e5) > k*c_s*t (3*3000*1=9e3) -> trips; and
//   assert slip_i (1e5) < slip_factor*V_peak*t (4*2e4*1=8e4)?  NO: 1e5>8e4 here,
//   but with V_peak=1e5 (the true peak) 4*1e5*1=4e5 > 1e5 -> old bound DOES NOT trip.
```

---

### [R-002] MODERATE [Phase 2, `slip_floor`, plan lines 169 & 184-186] — recommended `slip_floor = V_plate·tfinal` is ~1e-9 m, far too small to prevent the t→0 false trips it is meant to suppress

**Category:** BUG

**Description:**
For a SAFS *dynamic* run `V_plate ~ 1e-9 m/s` and `tfinal ~ O(s)`, so `V_plate·tfinal ~ 1e-9 m` —
seven-plus orders below physical event slip (~m). As a floor it is effectively zero, so it does
**not** prevent the `t→0` divide-by-tiny false trips the plan invokes it for (plan line 191:
"`slip_floor` prevents false trips"). The plan offers a correct alternative in the same breath
("or a small physical absolute (e.g. one `d_c`)"), but lists the broken option first as the
derived choice. (Note: once R-001 is fixed, the floor only guards very early time, but it still
must be a physical absolute.)

**Trigger:** early simulation time with the floor set to `V_plate·tfinal`.

**Actual behavior:** floor ≈ 0 → at small `t` any nonzero elastic slip can exceed the bound
(after R-001's fix) → false trips; before R-001's fix it is moot because the tripwire never fires.

**Expected behavior:** a physical absolute floor on the order of the slip-weakening distance.

**Suggested fix:**
```diff
- const real_t slip_floor   = /* derive: V_plate * tfinal, or a small absolute m */;
+ // Physical floor: ~ one slip-weakening distance d_c, so a DOF that has barely
+ // started weakening never trips. NOT V_plate*tfinal (~1e-9 m on a dynamic run).
+ const real_t slip_floor   = dof_data[0].lsw_d_c;   // or a small multiple thereof
```

**Test case:**
```cpp
// At t=1e-6 s with V_substep_peak_global=0 (pre-rupture) and a DOF with slip=0.1 m:
//   floor=V_plate*tfinal (~1e-9) -> slip_bound~k*c_s*1e-6 + 1e-9 ~ 3e-3 -> 0.1>3e-3 FALSE TRIP;
//   floor=d_c (0.4)             -> slip_bound~3e-3 + 0.4 = 0.403 -> 0.1<0.403 OK (no trip).
```

---

### [R-003] MODERATE [Phase 1, Files to Modify, plan lines 80-81 & Constraints line 48-49] — false premise: TPV205 uses `AdvanceWithSubStepStates`, not the no-callback `Advance`; the byte-exact argument is built on the wrong mechanism

**Category:** DEVIATION / ASSUMPTION (incorrect justification)

**Description:**
The plan states it keeps TPV205 byte-exact by **not touching the no-callback `Advance` overload
(`:137`)** and confines edits to `AdvanceWithSubStepStates`. But **TPV205 itself calls
`AdvanceWithSubStepStates`** — verified `drivers/tpv205_driver.cpp:414`
(`iterator.AdvanceWithSubStepStates(...)`), same overload SAFS uses
(`spatial_dyn_driver.cpp:439`). So "not touching `Advance`" protects nothing for TPV205; TPV205
runs straight through the modified code. Byte-exactness actually holds *only* because the edits are
inert (new field unread by TPV205 output; trace/`c_s` checks env-gated) — a more fragile property
than the plan claims, and one a maintainer could unknowingly break (e.g., later reading
`slip_rate_substep_max` in a shared output path). The wrong premise should be corrected so the
real invariant ("the AdvanceWithSubStepStates edits are inert for TPV205") is the one stated and
guarded.

**Trigger:** reading the plan's mitigation as written and assuming TPV205 is protected by overload
separation.

**Actual behavior:** plan asserts a protection that does not exist; byte-exactness rests on an
unstated, more fragile property.

**Expected behavior:** state that TPV205 *does* use `AdvanceWithSubStepStates`, and that
byte-exactness depends on the new field being write-only/unread by any TPV* output path and the
diagnostics being env-gated — then guard that with the parity regression (which the plan already
lists).

**Suggested fix (plan text):**
```diff
- Do **not** touch the no-callback `Advance` overload (`:137`) or `StepOneQP_` — keeps the `Q̄` path byte-exact.
+ NOTE: TPV205 ALSO calls `AdvanceWithSubStepStates` (tpv205_driver.cpp:414), so this overload is
+ NOT TPV205-free. Byte-exactness holds because (i) `slip_rate_substep_max` is written but read by
+ no TPV* output path, and (ii) the trace/`c_s` checks are env-gated. `StepOneQP_` and the
+ no-callback `Advance` stay untouched only to minimize surface area. Guard via the parity test.
```

**Test case:** `test_tpv205_evaluate_ader_lsw_parity` + TPV205 driver smoke byte-identical with
`SEAS_*` unset (already in the plan; this finding just corrects the *reason* it passes).

---

### [R-004] MODERATE [Phase 1, `DOFData` field, plan lines 73-75] — adding `slip_rate_substep_max` unconditionally contradicts the existing documented invariant that `DOFData` layout stays byte-for-byte; checkpoint/restart compatibility is unverified

**Category:** ASSUMPTION / EDGE_CASE

**Description:**
`DOFData` already has a field (`diag_print`) deliberately wrapped in `#ifdef SEAS_DIAG_FAULT_FLUX`
*specifically* "so struct layout matches pre-change byte-for-byte" in production
(`fault_face_flux.hpp:99-106`). The plan adds `slip_rate_substep_max` **unconditionally**, changing
`sizeof(DOFData)` and the field layout for every build. That directly violates the invariant the
existing code goes out of its way to preserve. The plan asserts (line 58-59) "never written to any
TPV* output, so struct-value output is unchanged" — true for *field-by-field* writers
(`station_writer.WriteStep`, the reconcile payload), but it does **not** address (a) any
checkpoint/restart that serializes `DOFData` by layout (the driver runs with `--checkpoint-every`),
or (b) the reason the `diag_print` author cared about layout in the first place. No raw
`sizeof(DOFData)`/`memcpy(DOFData)` serialization was found in a grep, but the checkpoint path was
not exhaustively traced, and the existing `#ifdef` precedent shows layout *is* treated as load-bearing.

**Trigger:** restarting from a checkpoint written by a pre-change binary; or any code that relies
on `DOFData` layout/size.

**Actual behavior:** struct layout changes for all builds; checkpoint/restart compatibility and the
documented byte-for-byte-layout invariant are unverified.

**Expected behavior:** either confirm `DOFData` is never serialized by layout (and update/relax the
`diag_print` comment to say so), or gate the new field consistently with `diag_print`.

**Suggested fix:** add an explicit checkpoint round-trip check to the acceptance criteria, and
EITHER document that `DOFData` is field-serialized only (so layout is free to change — and then the
`diag_print` `#ifdef` is over-cautious) OR gate the field:
```diff
+ // Verify before landing: does the checkpoint/restart path serialize DOFData by
+ // layout? If yes, gate slip_rate_substep_max behind a macro like diag_print, or
+ // bump the checkpoint version. The existing diag_print #ifdef (fault_face_flux.hpp:99)
+ // exists precisely to keep DOFData layout byte-stable — do not silently break it.
```

**Test case:**
```cpp
// Checkpoint round-trip: write dof_data with a pre-change build, restart with the new
// build; assert all physical fields restore bit-identical (or that the checkpoint
// serializes named fields, not raw struct bytes).
```

---

### [R-005] MODERATE [POSSIBLE] [Phase 1 AC line 141 + Phase 2 R-004] — `slip_rate_substep_max ≥ slip_rate` can be violated on SHARED fault QPs, and the "honest V_max" misses the inline shared-QP solve

**Category:** EDGE_CASE / ASSUMPTION

**Description:**
`slip_rate_substep_max` is maintained only in `AdvanceWithSubStepStates` (the iterator). But for
**shared** fault QPs, `ComputeADERSharedFaceFluxRHS` runs a *second* friction solve inline with the
macro `dt` AFTER the iterator and overwrites `dof_data.slip_rate` via `WriteBackState`
(`wave_operator.inl` Pass 1; `fault_face_flux.cpp:307`). That macro-`dt` `V_abs` is independent of
the iterator's sub-step `V_abs` and can exceed `slip_rate_substep_max`. So:
- the Phase-1 acceptance criterion `slip_rate_substep_max ≥ slip_rate` **for every DOF every step**
  (plan line 141) can FAIL on shared QPs (the clean 2-tet unit fixture won't expose it; production
  will);
- the "honest" `V_substep_peak_global` (R-004) does not see the inline shared-QP solve's velocity,
  so it under-reports there — the very blind spot R-004 is meant to close, reopened for shared QPs.

**Trigger:** a shared fault QP whose inline macro-`dt` solve yields `V_abs >` the iterator's
sub-step max at that QP.

**Suggested fix:** either (a) restrict the AC to interior QPs / state it holds only off the rank
seam, or (b) also update `slip_rate_substep_max` from the inline solve in
`ComputeADERSharedFaceFluxRHS` (and reconcile it like `slip_rate` should be — note this ties back
to the slip_rate-reconcile fix). Minimum: weaken the AC and document the shared-QP gap.
```diff
- [ ] `dof_data[i].slip_rate_substep_max ≥ dof_data[i].slip_rate` for every DOF every step
+ [ ] `slip_rate_substep_max ≥ slip_rate` for every INTERIOR fault DOF every step.  On SHARED
+     QPs `slip_rate` is overwritten by the macro-dt inline solve in ComputeADERSharedFaceFluxRHS,
+     so the invariant need not hold there until that solve also updates slip_rate_substep_max.
```

**Test case:**
```cpp
// np=2 shared-fault fixture: drive the inline macro-dt solve to V_abs greater than the
// iterator sub-step max at the shared QP; assert (pre-fix) slip_rate > slip_rate_substep_max
// there, demonstrating the AC violation; (post-fix b) the inline solve updates the max.
```

---

### [R-006] LOW [POSSIBLE] [Phase 1 trace, plan lines 100-102] — trace reconstructs `mu_eff` with `LSWFrictionCoefficient_TPV205`, but the SAFS shared-fault solve uses the forced-rupture coefficient; the printed `tau_str`/`mu_eff` may not match the shared-QP physics

**Category:** QUALITY (trace fidelity) / POSSIBLE pre-existing inconsistency

**Description:**
The trace recomputes `mu_eff = LSWFrictionCoefficient_TPV205(...)` (plan line 101), which matches
the iterator's `StepOneQP_:100` — good for the iterator path. But the SAFS *shared-fault* solve is
`EvaluateADER_LSW_ForcedRupture` (`fault_face_flux.cpp:917` →
`LSWFrictionCoefficient_ForcedRupture`). If `InitializeFaultDOFs_Spatial` sets active forced-rupture
parameters (CLAUDE.md notes it overwrites `T_forced_rupture`/`t0_decay_forced` per-DOF), then on
shared QPs the **iterator** (plain TPV205 LSW) and the **inline solve** (forced-rupture LSW) use
different `mu(δ,t)` for the same QP — a pre-existing discrepancy worth confirming, and the trace's
`tau_str`/`mu_eff` would reflect only the iterator's variant. If SAFS keeps `T_forced_rupture` at
the sentinel (forced-rupture inert), this is moot.

**Trigger:** SAFS DOFData with active forced-rupture fields, traced on a shared QP.

**Suggested fix:** (a) confirm whether SAFS gradual_overstress leaves `T_forced_rupture` at the
sentinel; if not, note in the trace that `mu_eff` is the iterator's TPV205 value and may differ
from the shared-QP inline solve; (b) separately verify the iterator using
`LSWFrictionCoefficient_TPV205` (not the forced-rupture variant) is intended for SAFS — if not,
that mismatch is itself a candidate speckle contributor and belongs in the deferred follow-up.

**Test case:** n/a (diagnostic-fidelity / verification note).

---

### [R-007] LOW [Phase 1 trace, plan line 93] — `#ifndef NDEBUG_SEAS_SLIP` wrapper is pointless and misleading

**Category:** QUALITY

**Description:**
The trace is wrapped in `#ifndef NDEBUG_SEAS_SLIP` ("always compiled; runtime-gated below"). Since
nothing defines `NDEBUG_SEAS_SLIP`, the guard is always true — it adds no compile-time control and
its name evokes `NDEBUG` (which it is unrelated to), inviting confusion. The runtime `slip_diag`
gate already provides the control.

**Suggested fix:** drop the `#ifndef`/`#endif` and keep only the runtime `slip_diag` gate; or, if a
compile-out is wanted, use a clearly-named opt-in macro (`#ifdef SEAS_ENABLE_SLIP_TRACE`) consistent
with the project's other diagnostic macros.

**Test case:** n/a (style; no behavior change).

---

## Summary
- Critical issues: 1 (R-001 — the slip-plausibility tripwire cannot fire as specified; it tests a
  kinematic identity)
- Moderate issues: 4 (R-002 floor magnitude; R-003 false byte-exact premise; R-004 unconditional
  `DOFData` field vs the documented layout invariant; R-005 shared-QP AC violation + monitor blind spot)
- Low issues: 2 (R-006 trace `mu_eff` fidelity on shared QPs; R-007 dead `#ifndef`)
- Plan compliance: **PARTIAL** — the observability substrate (Phase 1 trace + `slip_rate_substep_max`)
  and the finiteness guard (Phase 3) are sound and byte-safe; the **slip-plausibility tripwire
  (R-003 of the plan / the headline "localize" deliverable) is logically broken (R-001)** and the
  byte-exact justification (R-003) and `DOFData` layout handling (R-004) need correcting.
- **Verdict: PASS WITH FIXES.** Must fix R-001 (tripwire bound vs `c_s`, not observed peak) and
  R-002 (physical floor) before the tripwire is worth shipping; correct R-003/R-004 before relying
  on the byte-exact claim; address R-005 so the AC and monitor are honest on shared QPs. The trace
  (Phase 1) and finiteness guard (Phase 3) are good to implement as-is.

## Unreviewed Areas
- Whether the checkpoint/restart path serializes `DOFData` by raw layout (R-004) — grep found no
  raw `sizeof/memcpy(DOFData)`, but the `--checkpoint-every` code path was not opened end-to-end.
- The `WriteStep`/`StationWriter` field set — assumed field-by-field (so the new field is inert);
  not opened.
- Whether `V_max_global` is consumed anywhere beyond the env-gated `[DIAG]` block (e.g., an abort
  threshold or adaptive `dt`) — if it is, Phase 2's switch of the onset trigger to the sub-step max
  would be more than cosmetic; not traced.
- The deferred R-005/R-006/R-007 (tension, within-step relaxation, resolution) are correctly out of
  scope here and not reviewed.
