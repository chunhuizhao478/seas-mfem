# Code Review: TPV102 persistent-nucleation fix — round 3 (fresh audit)

Date: 2026-04-22 (round 3 — post round-2 fix application)
Branch: feature/elasticity-inertia
Commit under review: working tree on top of 632ad03 with round-2 fixes
applied to R-001..R-004 from the prior REVIEW.md (uncommitted).

## Review Scope

- Plan: `miniapps/seas/debug_document/tpv102_debug_document/tpv102_nucleation_code_fix_2026-04-22.md`
  plus round-2 review feedback in the preceding `REVIEW.md`.
- Prior bug analysis:
  `miniapps/seas/debug_document/tpv102_debug_document/tpv102_nucleation_code_review_2026-04-22.md`.
- Files reviewed (uncommitted working tree):
  - `miniapps/seas/dynamic/fault_face_flux.hpp`            (DOFData + signatures)
  - `miniapps/seas/dynamic/fault_face_flux.cpp`            (EvaluateTotal Steps 1–5)
  - `miniapps/seas/dynamic/tpv102_setup_total.hpp`         (ZeroDOFData + ApplyNucleationTotalPrestress)
  - `miniapps/seas/drivers/tpv102_driver.cpp`              (5 nucleation call sites, setup block)
  - `miniapps/seas/tests/unit/test_persistent_nuc_prestress_channel.cpp` (A/B/C + new A8)
  - `miniapps/seas/Makefile`                               (target + dependency wiring)
  - `miniapps/seas/jobs/tpv102/tpv102_200m_p1_2.0s_400rank_persistent_nuc_fix.sbatch`
- Cross-referenced against SeisSol sources:
  - `SeisSol/src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h`
    (`adjustInitialStress` at L422–450, imposed-state update at L336–363).
  - `SeisSol/src/DynamicRupture/FrictionLaws/CpuImpl/RateAndState.h`
    (`updateFrictionAndSlip` at L196–256 — `totalTraction` scale for friction,
    `tractionResults` scale for the Riemann imposed state).
- Domain context:
  - `miniapps/seas/CLAUDE.md` (sign conventions, BP5 FaultBasis, total-Q contract).

Starting premise: "assume at least 3 bugs." Re-executed all three passes
from scratch (not just verification of R-001..R-004 closure).

---

## Pass 1 — Round-2 findings closure

The round-2 REVIEW.md identified four findings (R-001/CRITICAL,
R-002/MODERATE, R-003/MODERATE, R-004/LOW). All four are now addressed by
the working tree:

- **R-001** — `fault_face_flux.cpp:302-316` introduces `const real_t
  {sigma_n,tau1,tau2}_fric = {sigma_n,tau1,tau2}_trial + data.{…}_nuc`
  as local temporaries; `tau*_trial` is never mutated. L333–342 uses the
  `_fric` scale for `Theta` and `solver_.Solve(…, |sigma_n_fric|, …)`.
  L349–360 computes the slip decomposition numerator on `_fric` and the
  Riemann-side `tau*_corr = tau*_trial - eta_s·V*` on the TRIAL scale
  (matching SeisSol `tractionResults.traction* = faultStresses.traction*
  - etaS·slipRate*` at `CpuImpl/RateAndState.h:235-240`). L367 leaves
  `sigma_n_corr = sigma_n_trial` (trial scale, matching SeisSol
  `imposedState[N] = faultStresses.normalStress` at
  `FrictionSolverCommon.h:347`). L393–398 stores the TOTAL-scale
  `data.{sigma_n,tau1,tau2}_corr = {sigma_n,tau1,tau2}_corr + data.{…}_nuc`
  for user-facing output (matching fluctuation-path `Evaluate` at
  `fault_face_flux.cpp:209-211`).
- **R-002** — Test A8 at `test_persistent_nuc_prestress_channel.cpp:207-250`
  now inspects `|Qip[VZ] - Qim[VZ]|` vs `d_a8.slip_rate` with a 5% relative
  + 1e-9 m/s absolute tolerance. Under the pre-R-001 buggy code this
  would read ~5.2 m/s vs V_abs ≈ 0.22 m/s and fail; under the fix both
  equal V_abs. Correctly gates the mutation regression.
- **R-003** — `SetSymmetricBackgroundQ` at `test_persistent_nuc_prestress_channel.cpp:121-129`
  now writes `Qp[SXZ] = Qm[SXZ] = +tau_ini` (face-local, post-rotation).
  Test A1 assertion updated to `+tau_ini` at L158-159.
- **R-004** — `ZeroDOFDataPreStressTotal` at `tpv102_setup_total.hpp:113-132`
  now also resets `sigma_n_nuc`, `tau1_nuc`, `tau2_nuc` to 0.

SeisSol cross-check (confirmed by reading `FrictionSolverCommon.h` and
`CpuImpl/RateAndState.h` directly): the fixed MFEM path is byte-for-byte
the SeisSol semantics. `totalTraction` (SeisSol) = `tau*_fric` (MFEM)
for friction; `tractionResults.traction*` (SeisSol) = `tau*_corr` on
trial scale (MFEM) for the Riemann imposed-state update;
initialStressInFaultCS / initialPressure (SeisSol) = `tau*_nuc /
sigma_n_nuc` (MFEM) and are NEVER added into the imposed velocity jump.

## Pass 2 — Fresh bug hunt

Re-ran all three passes on the changed files, not checking against any
prior finding list.

---

## Findings

### [R-005] [LOW] [fault_face_flux.hpp:32-46 + tests/unit/test_persistent_nuc_prestress_channel.cpp:273-275] — `tau*_nuc` modification invariant not asserted inside `EvaluateTotal` (analog of the `tau*_0 == 0` assert is missing)

**Category:** ASSUMPTION (contract documentation gap)

**Description:**
`EvaluateTotal` carries an `MFEM_ASSERT(data.sigma_n0 == 0.0 && data.tau1_0 == 0.0
&& data.tau2_0 == 0.0, …)` guard at `fault_face_flux.cpp:266-274` so a total-Q
driver that forgets `ZeroDOFDataPreStressTotal` aborts on the first friction
call (debug builds). No symmetric guard exists for the new `tau*_nuc` channel.
If a non-TPV102 total-Q driver forgets to call `ZeroDOFDataPreStressTotal` OR
`ApplyNucleationTotalPrestress`, the `tau*_nuc` fields silently keep whatever
bits were in the freshly-default-constructed `DOFData` (currently 0.0 thanks
to in-class initializers at `fault_face_flux.hpp:45-46`) — correct today, but
the invariant "`tau*_nuc` is 0 unless the caller opts in" is held only by
constructor convention, not by a runtime check.

Test B's "pre-poison tau2_nuc with garbage" (L273-275) relies on this exact
invariant being loose — a future refactor that tightens it would break Test B.

**Trigger:** A hypothetical total-Q driver that does NOT use nucleation yet
reuses `DOFData` with stale contents (e.g., checkpoint restart that deserialises
only a subset of fields).

**Actual behavior:** Stale `tau*_nuc` silently influences the friction solve
and the user-facing `data.tau*_corr` output without any warning.

**Expected behavior:** Either (a) document explicitly that `tau*_nuc` may be
nonzero outside the ZeroDOFDataPreStressTotal contract, or (b) add a debug-
only assert for drivers that wish to enforce `tau*_nuc == 0` (e.g., BP5).

**Suggested fix:** Documentation-only. Add a comment to
`fault_face_flux.cpp:266-274`:

```diff
    MFEM_ASSERT(data.sigma_n0 == 0.0 &&
                data.tau1_0   == 0.0 &&
                data.tau2_0   == 0.0,
                "EvaluateTotal contract violated: DOFData pre-stress "
                "fields must be zero (sigma_n0=" << data.sigma_n0 <<
                ", tau1_0=" << data.tau1_0 <<
                ", tau2_0=" << data.tau2_0 <<
                ").  Total-stress drivers must call "
-               "ZeroDOFDataPreStressTotal after InitializeFaultDOFs.");
+               "ZeroDOFDataPreStressTotal after InitializeFaultDOFs.  "
+               "Note: the separate {sigma_n,tau1,tau2}_nuc channel is "
+               "INTENTIONALLY not asserted here; those fields may be "
+               "nonzero (carrying the persistent nucleation forcing) "
+               "and are added to the friction input only — see Step 1b "
+               "below.  Callers that want to enforce tau*_nuc == 0 must "
+               "add their own guard.");
 #endif
```

**Test case:** Not required for LOW.

---

### [R-006] [LOW] [drivers/tpv102_driver.cpp:1282-1295] — RK4 stage 3 has no `ApplyNucleationTotalPrestress` call; relies on stage 2 having set `tau*_nuc` at `t + dt/2`

**Category:** ASSUMPTION (undocumented coupling)

**Description:**
RK4 stages 2 and 3 share the time argument `t + dt/2`. Stage 2 calls
`ApplyNucleationTotalPrestress(dof_data, fault_coords, t + dt_step/2.0)`
at `drivers/tpv102_driver.cpp:1265`. Stage 3 at L1282–1295 has NO such call
and assumes `tau*_nuc` is still set to the stage-2 value.

Since `ApplyNucleationTotalPrestress` is an overwrite-at-time helper (not an
accumulator) and no other code path touches `dof_data[i].tau*_nuc` between
stages, this assumption is currently correct. But there is no comment
explaining the absence of the call, and any future edit that inserts a
`ZeroDOFDataPreStressTotal` or similar between stages 2 and 3 would silently
zero `tau*_nuc` for stage 3 — producing a subtle O(dt) nucleation-amplitude
error distributed across each step.

**Trigger:** Any refactor that touches `dof_data` between RK4 stages 2 and 3
without being aware of the stage-3 implicit read.

**Suggested fix:** Add a one-line comment at `drivers/tpv102_driver.cpp:1282`:

```diff
-      // RK4 stage 3: at time t + dt/2, psi = psi_n + (dt/2) * psi_k2
+      // RK4 stage 3: at time t + dt/2, psi = psi_n + (dt/2) * psi_k2.
+      // No ApplyNucleationTotalPrestress call needed — stage 3 shares
+      // the stage 2 time argument (t + dt/2), and dof_data[i].tau*_nuc
+      // is still set to NucleationPerturbation(..., t + dt/2) from the
+      // stage 2 call above.  Do NOT insert code that touches tau*_nuc
+      // between stage 2 and stage 3.
       add(Q, dt_step / 2.0, k2, Q_tmp);
```

**Test case:** Not required for LOW.

---

### [R-007] [LOW] [dynamic/tpv102_setup_total.hpp:716-770] — Legacy bulk-Q path (`ApplyNucleationTotal`, `FaultQPNodalMap`, `BuildFaultQPNodalMap`, `FaultQPNucleationState`, `VerifyFaultQPNodalMapUnique`) is retained "to keep R-002 / R-009 / R-I06-003 unit tests green" but no longer exercises any real code path

**Category:** QUALITY

**Description:**
The diff in `tpv102_setup_total.hpp` retains ~470 LoC of dead-code machinery
for the bulk-Q delta-inject path that the round-1 review showed is
fundamentally wrong. The driver no longer calls any of those APIs, and the
only remaining callers are the three unit tests
`test_R002_nucleation_sign_match.cpp`,
`test_R009_nucleation_total_both_sides.cpp`,
`test_R_I06_003_nucleation_shape_scaling.cpp` — which test per-call arithmetic
on a broken API. Keeping them green carries no real value: they verify
properties of code that is never invoked in production.

This risks a future contributor:
1. Reading the retained doc comment (L716–736) that describes how the
   bulk-Q injection works and assuming it is the current approach.
2. Adding a new total-Q driver that calls `ApplyNucleationTotal` because
   the function still exists, walking into the round-1 bug.
3. Spending time "fixing" the bulk-Q path when the persistent-prestress
   channel has superseded it.

**Trigger:** Any contributor browsing `tpv102_setup_total.hpp` without
reading the debug-document chain first.

**Suggested fix:** After a Frontera validation run on the
persistent-prestress branch confirms the fix works end-to-end, delete:
- `FaultQPNodalMap`, `FaultQPNucleationState`,
  `ApplyNucleationTotal(Vector &Q, …)`,
  `BuildFaultQPNodalMap`, `PickNearestNodalDOF`,
  `VerifyFaultQPNodalMapUnique`.
- The three unit tests above.
- The Makefile entries for those tests.

Until then, annotate the deprecated API at the top of the file:

```diff
+// -----------------------------------------------------------------------
+// DEPRECATED (2026-04-22): The bulk-Q delta-inject nucleation path below
+// (FaultQPNodalMap, FaultQPNucleationState, ApplyNucleationTotal,
+// BuildFaultQPNodalMap, PickNearestNodalDOF, VerifyFaultQPNodalMapUnique)
+// is superseded by ApplyNucleationTotalPrestress at the bottom of this
+// file.  The production TPV102 driver no longer calls any of these
+// symbols; they remain only to keep three legacy unit tests building.
+// DO NOT use these in new code.  See
+//   debug_document/tpv102_debug_document/
+//     tpv102_nucleation_code_review_2026-04-22.md  (bug analysis)
+//     tpv102_nucleation_code_fix_2026-04-22.md     (replacement)
+// -----------------------------------------------------------------------
 /// @brief Per-fault-QP nodal-injection map (I-06 Phase 5).
```

**Test case:** Not required for LOW.

---

### [R-008] [LOW] [tests/unit/test_persistent_nuc_prestress_channel.cpp:103-129] — `SetSymmetricBackgroundQ` sign comment correct, but the Q[VX], Q[VY], Q[VZ] are implicitly zero-initialised and the test depends on that; make it explicit

**Category:** QUALITY

**Description:**
`SetSymmetricBackgroundQ` loops `for (int c = 0; c < NUM_STATE; c++) { Qp[c]
= 0; Qm[c] = 0; }` then sets only `SXX` and `SXZ`. All other components
(including VX, VY, VZ, SYY, SZZ, SYZ) are left at zero. The friction-solve
correctness critically depends on VX=VY=VZ=0 (equilibrium fault), SYY=0
(no added normal stress in that direction), and SXY=0 (no dip shear).

This is implicitly correct but is load-bearing for Test A's baseline
assertions. A future edit that reorders the initialisation or swaps in a
different fixture builder must preserve the zero-init of the other components.

**Trigger:** Any edit to the fixture builder.

**Suggested fix:** Add explicit comments:

```diff
 void SetSymmetricBackgroundQ(real_t *Qp, real_t *Qm)
 {
-   for (int c = 0; c < NUM_STATE; c++) { Qp[c] = 0; Qm[c] = 0; }
+   // Zero every wave-state component first: equilibrium fault has
+   // VX=VY=VZ=0 (no bulk motion), SYY=SZZ=SYZ=0 (no auxiliary stress),
+   // and SXY=0 (no dip shear; TPV102 is pure strike-slip).  The only
+   // nonzero components are the fault-normal compression (SXX) and
+   // the along-strike pre-stress (SXZ) set below.
+   for (int c = 0; c < NUM_STATE; c++) { Qp[c] = 0; Qm[c] = 0; }
    // SXX is the fault-NORMAL stress in the rotated frame.
    Qp[SXX] = TPV102Params::sigma_n;  Qm[SXX] = TPV102Params::sigma_n;
    // SXZ = sigma_NT2 = +tau_ini in fault-local coordinates (sign flip
    // from global -tau_ini under the BP5 rotation; see header above).
    Qp[SXZ] = +TPV102Params::tau_ini; Qm[SXZ] = +TPV102Params::tau_ini;
 }
```

**Test case:** Not required for LOW.

---

## Pass 3 — Quality / maintainability

- **Variable naming** (`tau*_fric` vs `tau*_trial` vs `tau*_corr`) is distinct
  and accurate; each local has a comment tying it to the SeisSol analog.
  No finding.
- **Function responsibilities** in `EvaluateTotal` (trial → friction →
  decomposition → imposed state → output) are cleanly separated with
  headers for each step. No finding.
- **Error handling**: `MFEM_VERIFY` on size mismatches in
  `ApplyNucleationTotalPrestress`. `MFEM_ASSERT` on the
  `{sigma_n,tau1,tau2}_0 == 0` contract preserved. No finding.
- **Dead code**: Flagged as R-007.
- **Inline comment density**: the fix added ~50 LoC of explanatory comments
  to `EvaluateTotal`. Higher than the rest of the file, but justified given
  the round-1/round-2 debug cost.

---

## Summary

- Critical issues: **0**
- Moderate issues: **0**
- Low issues: **4** (R-005 through R-008 — all documentation / dead-code
  hygiene, none affect correctness).
- Plan compliance: **FULL** — R-001..R-004 from round 2 all addressed with
  matching code; SeisSol semantics confirmed by direct cross-reference to
  `CpuImpl/RateAndState.h:222-240` and
  `FrictionSolverCommon.h:336-363`.
- Verdict: **PASS** — the persistent-prestress channel fix is now
  mathematically consistent with SeisSol, the Riemann imposed-velocity
  overshoot is eliminated, and Test A8 gates regressions.
  **The Frontera sbatch `tpv102_200m_p1_2.0s_400rank_persistent_nuc_fix.sbatch`
  is safe to submit from a CORRECTNESS standpoint.** The four LOW findings
  are cleanups to apply post-run or in a follow-up commit; none of them
  block or affect the Frontera V_max outcome.

### Expected Frontera V_max trajectory (R-001 fix baseline)

For the 400-rank / tfinal=2.0 s run with the fix:
- t < 0.3 s: `V_max ≤ 1e-9 m/s` (locked fault, G(t) still small).
- t ∈ [0.3, 0.6 s]: exponential climb as G(t) ramps from ~0.4 → ~0.8;
  `V_max` reaches ~1e-4 m/s by t ≈ 0.5 s (DRDG3D reference: 5.5e-5).
- t ≈ 0.8 s: `V_max ~ 5e-2 m/s` (DRDG3D: 4.8e-2).
- t ≈ 1.0 s: `V_max > 0.2 m/s`, rupture breaks out of the hypocentral patch.
- t ∈ [1.0, 1.5 s]: propagating rupture; `V_max` reaches 1–3 m/s.
- 6/9 SCEC stations slipped by t = 2.0 s.

If the Frontera run fails to follow this trajectory, investigate a second
bug — the persistent-prestress channel correctness is now validated from
first principles + SeisSol cross-reference + Test A8 gating, so a remaining
mismatch points to either the ADER corrector, the shared-fault
communication path, or a mesh / CFL interaction.

## Unreviewed Areas

- `wave_operator.inl` paths (`ComputeADERTimeIntegrated`, `ComputeADERFaceFluxRHS`,
  `ComputeADERSharedFaceFluxRHS`): these are unchanged by the current fix
  and were reviewed indirectly in round 1 — the R-001 bulk-Q radiation
  issue was fixed at the call site (driver) rather than inside the wave
  operator. If the Frontera run exposes a new issue, `wave_operator.inl`
  is the next place to look.
- `ApplyNucleationTotal` (legacy bulk-Q path): retained per R-007 note
  but no longer driver-reachable. Not re-reviewed here because it no
  longer affects production.
- Shared-fault cross-rank DOFData consistency (`R-101` check in
  `wave_operator.inl`): not re-reviewed; the log on the prior-baseline
  run showed "shared: 0" for the 400-rank / 200 m mesh, so the shared-
  fault path is not exercised in the sbatch. A future finer mesh or
  different partitioning would need shared-fault re-validation of the
  `tau*_nuc` replication (both ranks independently call
  `ApplyNucleationTotalPrestress` on their own `dof_data`; since
  `fault_coords` for a shared QP is geometrically identical on both
  sides, both ranks produce the same `tau*_nuc` — this appears safe
  by construction but is not unit-tested here).
