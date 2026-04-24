# TPV102 Debug v9.4.0 — SeisSol-Style Split-Prestress to Eliminate Forced-Response Pepper Artifact

> **Author:** review agent, 2026-04-22.
> **Predecessors:**
> - `tpv102_debug_v9.3.0_debug_plan.md` — total-stress migration (I-06).
> - `tpv102_nucleation_code_review_2026-04-22.md` — round-1 bulk-Q radiation
>   bug identification (V_max plateau at 2e-7 m/s on Frontera job 7672719).
> - `tpv102_nucleation_code_fix_2026-04-22.md` — round-2 persistent-prestress
>   channel fix landed as commit `28284b9` (`tau*_nuc` fields in DOFData;
>   `EvaluateTotal` friction-input uses `tau*_fric = tau*_trial + tau*_nuc`;
>   Riemann-side `tau*_corr` stays on trial scale).
>
> **Status:** PLAN ONLY. No code changed by this document.
>
> **Scope:** Close out the second class of TPV102 artifact — the "pepper"
> cell-scale noise pattern that appears on the fault surface under rupture
> drive, survives the round-2 persistent-prestress fix, is absent in
> no-nucleation runs, and decays when a one-shot perturbation is applied.
> This plan proposes an architectural revert to SeisSol's split-prestress
> fault representation and provides concrete, file-level patch targets.
>
> **Origin:** User-supplied plan (session transcript 2026-04-22) based on
> a floating-point forced-response analysis of the fault-QP trace
> reconstruction path. This document reproduces the plan verbatim and
> augments it with an adversarial review (concerns, risks, suggested
> hardening) per the review playbook.

---

## 1. Root Cause Hypothesis (PROPOSED)

Quoted from the origin plan:

> The TPV102 pepper artifact is a **forced-response numerical error**, not a
> physical instability or an MPI mapping bug. In MFEM's fault-face path, the
> local trace state at each fault quadrature point is reconstructed from
> large total-stress background values using
> `Q_self(q) = Σ_i shape_i(q) · Q_i`. For a mathematically uniform background
> state, the three fault QPs should reconstruct identical values, but in
> floating point the shape-weight permutations differ across QPs, so the
> accumulation order changes and produces ULP-scale per-QP differences.
> Those tiny differences are re-injected every step, especially in the
> `sigma_n` channel, then mixed by the nonlinear friction solve into
> traction/slip-rate differences. A one-time perturbation decays, but the
> per-step reconstruction keeps re-forcing the system, so the DG fault/bulk
> coupling accumulates the effect into visible per-element pepper under
> rupture drive.

**Observational support (from prior runs):**

1. No-nucleation cases stay clean — consistent with "ULP seed × 0 rupture
   drive = no amplification".
2. One-shot stress perturbations decay — consistent with "no persistent
   re-injection source".
3. Rupture-driven runs develop pepper — consistent with "friction
   nonlinearity amplifies a continuously re-seeded tiny noise".

**Numerical magnitude estimate (my check):**
- `sigma_n = 1.2e8 Pa` on every bulk DOF under total-Q.
- `1 ULP(1.2e8) ≈ 1.2e8 · 2^-52 ≈ 2.7e-8 Pa`.
- Under a split-prestress representation, bulk Q carries only the
  dynamic fluctuation (max ~1e6 Pa stress drop during rupture).
- `1 ULP(1e6) ≈ 2.2e-10 Pa`.
- Ratio: ~10^5 reduction in seed noise amplitude.

This supports the plan's claim that bulk Q scale matters. It is,
however, a HYPOTHESIS — see §6 R-001 for the gating requirement before
committing a 4-commit refactor.

---

## 2. Target Contract (PROPOSED)

At the fault, use the SeisSol split representation:

- **Dynamic part:** comes from interpolated bulk traces `Q_plus / Q_minus`.
  Bulk Q carries **fluctuations only** (no static prestress, no nucleation).
- **Persistent fault storage (per-QP `DOFData`):**
  - `sigma_n0, tau1_0, tau2_0` — static background prestress.
  - `sigma_n_nuc, tau1_nuc, tau2_nuc` — time-varying nucleation forcing.
  - Together the MFEM equivalent of SeisSol's `initialStressInFaultCS`.
- **Friction sees:** `total traction = persistent fault stress + dynamic trial traction`,
  computed pointwise per QP. Bulk Q is never the source of the large
  background magnitude.
- **Riemann imposed state:** uses fluctuation-scale `tau*_corr = tau*_trial - eta_s·V*`
  so the velocity jump carries only `-eta_s·V` (matching SeisSol
  `tractionResults.traction*` at `CpuImpl/RateAndState.h:235-240`).
- **User output:** `data.tau*_corr = data.tau*_0 + data.tau*_nuc + tau*_corr`
  (TOTAL physical traction, for station / VTU writers).

---

## 3. Patch Plan (PROPOSED — 5 commits)

Referenced line numbers verified against working tree at commit `28284b9`.

### Commit 1 — `FaultFaceFlux::Evaluate()` reads `tau*_nuc`

**Files:**
- `miniapps/seas/dynamic/fault_face_flux.cpp` (`Evaluate` at line 70)
- `miniapps/seas/dynamic/fault_face_flux.hpp` (doc comment on the
  `Evaluate` / `EvaluateTotal` contract)

**Edits:**

1. In `Evaluate` (around `fault_face_flux.cpp:116-119`), extend the
   friction input to include the nucleation channel:
   ```diff
   -  real_t sigma_n_total = data.sigma_n0 + sigma_n_trial;
   -  real_t tau1_total = data.tau1_0 + tau1_trial;
   -  real_t tau2_total = data.tau2_0 + tau2_trial;
   +  real_t sigma_n_total = data.sigma_n0 + data.sigma_n_nuc + sigma_n_trial;
   +  real_t tau1_total    = data.tau1_0   + data.tau1_nuc    + tau1_trial;
   +  real_t tau2_total    = data.tau2_0   + data.tau2_nuc    + tau2_trial;
   ```
   Note: `tau*_corr` (TRIAL scale) stays on the Riemann path unchanged.

2. Extend the user-output `data.*_corr` storage
   (`fault_face_flux.cpp:209-211`):
   ```diff
   -  data.tau1_corr = data.tau1_0 + tau1_corr;
   -  data.tau2_corr = data.tau2_0 + tau2_corr;
   -  data.sigma_n_corr = data.sigma_n0 + sigma_n_corr;
   +  data.tau1_corr    = data.tau1_0   + data.tau1_nuc    + tau1_corr;
   +  data.tau2_corr    = data.tau2_0   + data.tau2_nuc    + tau2_corr;
   +  data.sigma_n_corr = data.sigma_n0 + data.sigma_n_nuc + sigma_n_corr;
   ```

3. Leave `EvaluateTotal` (`fault_face_flux.cpp:247`) structurally in
   place but annotate as experimental/legacy in the class docstring.

**Why first:** Isolates the physics contract change inside one function
without touching driver, `wave_operator.inl`, or init paths. Any caller
with `tau*_nuc == 0` (the default) sees byte-identical behavior to the
current `Evaluate`. BP5 and existing fluctuation-path tests should be
unaffected — **verify this explicitly** (see §5 R-003 below).

### Commit 2 — TPV102 driver reverts to fluctuation bulk Q

**Files:**
- `miniapps/seas/drivers/tpv102_driver.cpp`
- `miniapps/seas/dynamic/tpv102_setup.hpp`
- `miniapps/seas/dynamic/tpv102_setup_total.hpp`

**Edits:**
1. Remove production use of `InitializeStateTotal` (`tpv102_setup_total.hpp:61`).
   Initialize bulk Q as fluctuation state: `Q = 0` (then any nucleation-
   free dynamic perturbation).
2. Remove production use of `ZeroDOFDataPreStressTotal`
   (`tpv102_setup_total.hpp:113`).
3. `InitializeFaultDOFs` (`tpv102_setup.hpp:40`) retains its existing
   writes to `sigma_n0`, `tau1_0 = 0`, `tau2_0 = tau_ini`.
4. Rename `ApplyNucleationTotalPrestress` (`tpv102_setup_total.hpp:761`)
   to a non-total-specific name, e.g., `ApplyNucleationPrestress`, and
   keep it as the ONLY production nucleation entry point.
5. Delete (or gate under `#if 0`) the bulk-Q nucleation entry point
   `ApplyNucleationTotal` (`tpv102_setup_total.hpp:654`) so it cannot
   accidentally be re-used; see §5 R-005.

### Commit 3 — `WaveOperator` fault dispatch reverts to `Evaluate` / `EvaluateADER`

**Files:** `miniapps/seas/dynamic/wave_operator.inl`

Replace at the four verified call sites:

| Site | Line | Current | Replace with |
|------|------|---------|--------------|
| RK4 interior fault | 1130 | `EvaluateTotal` | `Evaluate` |
| RK4 shared fault   | 1615 | `EvaluateTotal` | `Evaluate` |
| ADER interior fault| 1957 | `EvaluateADERTotal` | `EvaluateADER` |
| ADER shared fault  | 2246 | `EvaluateADERTotal` | `EvaluateADER` |

Do **not** change: canonical frame reconstruction, shared-face
role/tiebreaker logic, per-QP loop structure. Also: remove or relax any
`has_bulk_bg_` guards that assume total-Q dispatch; they no longer hold.

### Commit 4 — Tests

**Priority:**
- `miniapps/seas/tests/parallel/test_shared_fault_dof_data_consistency.cpp`
- `miniapps/seas/tests/unit/test_persistent_nuc_prestress_channel.cpp`
- `miniapps/seas/tests/unit/test_fault_face_flux_total_vs_fluctuation.cpp`
  (confirmed exists)
- `miniapps/seas/tests/parallel/test_shared_fault_role_consistency.cpp`
  (confirmed exists)

**New gate:** "TPV102 split-prestress path" — `Q_bulk = 0`, prestress in
DOFData, nucleation in `tau2_nuc`, assert fault output matches expected
total traction.

### Commit 5 — OPTIONAL arithmetic-parity improvements

Precompute `T_can / Tinv_can` per QP; precompute face interpolation
matrices. Only after Commits 1–4 confirm pepper reduction.

---

## 4. Plan's Short Summary (verbatim)

> To match SeisSol closely enough to avoid this error class:
> - bulk fault traces must be dynamic/fluctuation only
> - prestress and nucleation must live in persistent fault storage
> - friction must add those to dynamic trial traction pointwise
> - do not interpolate the large background prestress from bulk Q every step

---

## 5. Review / Comments / Suggestions

Per the review playbook: assume at least three bugs. Below are 7
findings — concerns with the plan that should be resolved (or
consciously accepted) before the refactor lands. IDs are referenced by
the fix agent.

### [R-001] [CRITICAL] Root-cause hypothesis is unverified — gate the 5-commit refactor on a local pepper reproducer

**Category:** ASSUMPTION

**Description:**
The plan invests 5 commits / multi-day refactor on the hypothesis that
per-step floating-point re-injection on a `sigma_n ~ 1.2e8 Pa` bulk state
is the pepper seed. The ULP math is plausible (~10^5 seed-amplitude
reduction under split-prestress) but none of the following has been
independently confirmed:

1. The pepper observed on Frontera is ACTUALLY per-QP ULP-scale
   reconstruction noise (vs. a mesh-scale instability, an MPI ghost-state
   drift, or a friction-solver conditioning issue).
2. The friction nonlinearity amplification factor is large enough to
   make a 2.7e-8 Pa seed visible (the factor must be ~10^6 to reach kPa
   station output — empirically plausible under high `psi` but not
   measured).
3. A 10^5 reduction in seed AMPLITUDE actually brings pepper below the
   visibility threshold (the system may be fundamentally
   instability-dominated, in which case reducing seed amplitude only
   delays onset by a few e-folds).

**Trigger:** Any acceptance decision based on "the split-prestress
refactor should fix pepper."

**Expected evidence before Commit 1:** A standalone local reproducer —
short (≤ 0.2 s), ranks-1–4, small mesh — that replays the pepper
phenomenon in a way that lets us:
- Quantify the seed amplitude per step (probe `Q_self[SXX]` dispersion
  across the 3 QPs of one face with `Q_bulk` uniformly set to
  `sigma_n0`).
- Confirm the seed scales as `sigma_n × eps_mach` rather than something
  mesh- or psi-dependent.
- Measure the growth rate of per-QP dispersion over ~100–1000 steps
  under nontrivial rupture drive (so the amplification factor is
  observable).

**Suggested fix:** Before Commit 1, author
`tests/unit/test_pepper_reproducer.cpp` (serial, fast) that:

```cpp
// Pseudocode — the goal is to expose the per-QP ULP dispersion, NOT
// replicate the full Frontera trajectory.
void test_R001_pepper_seed_is_ulp_per_qp() {
    // Build a 2-element tet mesh with one fault face, order 1,
    // Dunavant 3-pt face QP rule (3 QPs).
    // Fill bulk Q[SXX] = sigma_n0 = 1.2e8 at every DOF.
    // Rotate to face-local at each QP.  Assert:
    //   (a) at least two of the three QP values differ (ULP dispersion
    //       exists)
    //   (b) max |Q_self[SXX, QP_i] - Q_self[SXX, QP_j]| < 100 * ULP(sigma_n0)
    //       (dispersion is bounded at the expected scale)
    //   (c) the dispersion amplitude scales linearly with sigma_n0
    //       (parameter sweep sigma_n0 ∈ {1e6, 1e7, 1e8, 1e9}).
}
```

If this test passes the ULP bound AND points (a–c) hold, R-001 is
closed and the refactor is justified. If any of (a–c) fails, the root
cause is elsewhere (e.g., the friction solver) and the plan as-written
will not fix pepper — possibly making it worse by changing ADER
semantics while leaving the true bug in place.

**Test case:** As above — a standalone reproducer is itself the test.

---

### [R-002] [MODERATE] Plan does not address the ADER time-integration path's treatment of nucleation; could re-introduce the round-1 V_max plateau bug

**Category:** DEVIATION

**Description:**
The plan's Phase 2 item 6 says for ADER:

> minimum parity step: update persistent nucleation before each
> `AdvanceADER(...)`
> strict SeisSol parity step: redesign ADER fault update so nucleation
> is applied at each time substep/quadrature node, like SeisSol's
> `timeIndex` loop in
> `SeisSol/src/DynamicRupture/FrictionLaws/CpuImpl/BaseFrictionLaw.h:98`.

The "minimum parity" variant is what commit `28284b9` already does (one
`ApplyNucleationTotalPrestress` call at `t + dt/2` per step) and the
round-1 review flagged it as O(dt²) acceptable. Going back to
fluctuation-Q via `EvaluateADER` does NOT automatically ensure the
nucleation picture is consistent — specifically, `EvaluateADER`
(`fault_face_flux.cpp:399-424`) time-averages the input `I±` to `Q_avg`
and dispatches through `Evaluate`. The friction solver then reads the
CURRENT `data.tau*_nuc`, which is a SNAPSHOT at the step-start
midpoint, NOT time-integrated over `[t, t+dt]`.

If the plan holds `data.tau*_nuc = NucleationPerturbation(t_mid)`
constant over the ADER predictor stages, this is formally identical to
the round-1 behavior — fine for correctness but does not match SeisSol's
`timeIndex`-loop semantics either. The "strict SeisSol parity" variant
would require threading a time-varying `tau*_nuc(t_k)` through the CK
predictor's Cauchy-Kovalevskaya expansion, which is a non-trivial
redesign not detailed in the plan.

**Trigger:** Implementation ambiguity — the plan leaves the choice between
"minimum parity" and "strict parity" to the implementer and does not
specify whether the pepper reproducer requires the strict variant.

**Suggested fix:** Explicitly commit to the minimum-parity variant in
this document (Phase 2, Commit 3 scope). Document that strict SeisSol
ADER-sub-step nucleation is deferred to v9.5.0. Add the following
clarification to the plan's Phase 2 item 6:

```
ADER: the MINIMUM-PARITY approach is the v9.4.0 target.  Call
ApplyNucleationPrestress(dof_data, fault_coords, t + dt/2) once per
ADER step before wave.AdvanceADER.  data.tau*_nuc is then CONSTANT
for the entire step's predictor/corrector; the CK expansion treats it
as static.  O(dt²) error in the nucleation ramp amplitude is
documented as accepted under v9.4.0.  Strict sub-step parity is
deferred to v9.5.0 and is NOT a prerequisite for removing pepper.
```

**Test case:**
```cpp
// test_R002_v940_ader_nuc_timing.cpp
// Run 100 ADER-O(2) steps at full nucleation amplitude.
// Probe: the time-averaged fault traction in ParaView matches the
// midpoint nucleation value, not the endpoint or the step-start value.
// Tolerance: 5% of nuc_dtau (O(dt^2) bound at CFL 0.5).
```

---

### [R-003] [MODERATE] Adding `tau*_nuc` reads to `Evaluate` changes the BP5 fault-flux contract; must verify BP5 initializes `tau*_nuc = 0` and no BP5 code path ever writes nonzero `tau*_nuc`

**Category:** ASSUMPTION

**Description:**
Commit 1 extends `Evaluate()` — which is BP5's production fault flux — to
read `data.tau*_nuc` into the friction input (and `data.*_corr` output).
BP5 never writes to `tau*_nuc`; with the default `0` field initializer
at `fault_face_flux.hpp:45-46`, BP5 should see zero. BUT:

1. Any BP5 code path that zero-copies or memcpy's DOFData from a
   previous run, a checkpoint, or an uninitialised buffer could carry
   stale bits into `tau*_nuc`. The default-ctor initialiser only fires
   on a FRESH object.
2. `InitializeFaultDOFs` (`tpv102_setup.hpp:40-104`) does NOT touch
   `tau*_nuc` — BP5's own init helper may have the same omission.

If a single `tau*_nuc` entry on a BP5 fault leaks nonzero, BP5 cycle
simulations silently change behavior. No test currently catches this.

**Trigger:** Any BP5 run executed after Commit 1 lands.

**Suggested fix:** In `InitializeFaultDOFs` (`tpv102_setup.hpp:46-48`),
explicitly zero the nucleation fields:

```diff
    for (int i = 0; i < ndof; i++)
    {
       DOFData &d = dof_data[i];
+      // Nucleation channel (TPV102-specific; BP5 never writes these
+      // but Evaluate now reads them — zero here so a reused DOFData
+      // vector with stale bits cannot leak a nonzero nucleation into
+      // BP5's friction solve).
+      d.sigma_n_nuc = 0.0;
+      d.tau1_nuc    = 0.0;
+      d.tau2_nuc    = 0.0;
       // Impedances (homogeneous)
```

Also: mirror this in BP5's analogous init helper (search for
`std::vector<DOFData>` initializers in `bp5/`).

**Test case:**
```cpp
// test_R003_bp5_tau_nuc_invariant.cpp
// Run one BP5 cycle with a deliberately-prepoisoned DOFData.tau2_nuc =
// 1e6 Pa on one QP.  Expected (after fix): BP5 zero-clears it, cycle
// recurrence unchanged.
// Before fix: BP5 friction sees a 1 MPa nucleation perturbation,
// recurrence interval shifts measurably.
```

---

### [R-004] [MODERATE] `has_bulk_bg_` / `SetAbsorbingBackground` contract change in Phase 3 needs a dedicated regression — the v9.3.0 migration documented a BC correctness bug the fluctuation path had

**Category:** DEVIATION (regression risk)

**Description:**
Round-6 v9.3.0 purpose change #2 added `SetAbsorbingBackground(Q_bg)` and
total-Q-aware BC variants (`AbsorbingTotal`, `FreeSurface*Total`, `PML`
damping toward `Q_bg`) specifically because the fluctuation-BC path had
an identified radiation/damping bug against a nonzero background
(`drivers/tpv102_driver.cpp:510-519` comment explicitly says this).

The plan's Phase 3 item 7 reverts this — for fluctuation Q, `Q_bg = 0`
is natural, but that says NOTHING about whether the pre-v9.3.0 BC
variants (`Absorbing` without `Total`, `FreeSurface` without `Total`,
etc.) are still callable and correct. The plan does NOT enumerate:
1. Whether the old fluctuation-BC variants STILL EXIST in
   `wave_operator.inl` (they may have been deleted in the v9.3.0
   migration — grep required).
2. Whether `has_bulk_bg_` checks are a hard `MFEM_VERIFY` (which would
   abort the run after commit 3) or soft `MFEM_ASSERT` (only in debug).
3. What the expected free-surface reflection coefficients are for the
   fluctuation BC on top of a nonzero DOFData pre-stress.

**Trigger:** Commit 3 land + Frontera rerun with the default
absorbing-BC setting.

**Suggested fix:** Before Commit 3, run:
```bash
grep -n "SetAbsorbingBackground\|has_bulk_bg_\|FreeSurface\|FreeSurfaceTotal\|AbsorbingTotal\|Absorbing\b" \
     miniapps/seas/dynamic/wave_operator.inl
```
and enumerate every call site. For each:
- Does it have a fluctuation-path analog?
- Is the analog correctness-tested against the v9.3.0 total-Q path?
- What's the test coverage for the fluctuation BC on a nonzero-prestress
  fault?

Add a "BC correctness regression" subtest to Commit 3 that drives an
absorbing boundary with a known outgoing P-wave amplitude and asserts
the reflection coefficient is `< 1e-3` (the standard absorbing-BC
benchmark) under BOTH total-Q and the restored fluctuation path.

**Test case:**
```cpp
// test_R004_absorbing_bc_fluctuation_path.cpp
// Build a 1D-like setup: long box, one end absorbing, inject a P-wave
// pulse at the other end.  Fluctuation path: Q_bg = 0, DOFData prestress
// nonzero.  Assert |reflection| < 1e-3 after the pulse has left.
// If this fails, Phase 3 needs a dedicated BC fix, not a "just revert".
```

---

### [R-005] [MODERATE] Phase 2 item 5 proposes deleting `ApplyNucleationTotal` (bulk-Q path), but 3 unit tests still exercise it — the plan's test-repair budget is understated

**Category:** DEVIATION (scope)

**Description:**
The plan's Phase 2 item 5 says "Do not use `tpv102_setup_total.hpp:654`
in production" — meaning `ApplyNucleationTotal(Vector &Q, ...)`. Plan
Commit 4 lists 3 priority test edits but does NOT enumerate:

- `tests/unit/test_R002_nucleation_sign_match.cpp`
- `tests/unit/test_R009_nucleation_total_both_sides.cpp`
- `tests/unit/test_R_I06_003_nucleation_shape_scaling.cpp`

All three call `ApplyNucleationTotal` directly. If that function is
deleted, those 3 tests fail to compile. If it is retained but marked
deprecated, those tests exercise dead code — wasted CI budget — and the
new contributor risk flagged in the round-3 R-007 remains.

**Trigger:** `make test-v93-regression` after Commit 2 or Commit 4.

**Suggested fix:** Commit 4 must also:

```diff
# Makefile
-TEST_R002_NUCLEATION_SIGN_SRC     = tests/unit/test_R002_nucleation_sign_match.cpp
-TEST_R009_NUCLEATION_BOTH_SIDES_SRC = tests/unit/test_R009_nucleation_total_both_sides.cpp
-TEST_R_I06_003_NUC_SHAPE_SCALING_SRC = tests/unit/test_R_I06_003_nucleation_shape_scaling.cpp
```

and `git rm` the three `.cpp` files, plus remove their targets from
`test` and `test-v93-regression` lists. Those tests verify
per-call arithmetic on a now-dead API and keeping them green is
counterproductive.

**Test case:** Compile-time (`make all test`) — with `ApplyNucleationTotal`
deleted, `make test` must pass. Without the test deletions above, it
fails at link time.

---

### [R-006] [LOW] Plan's Phase 4 "optional precomputed rotations" is also a floating-point noise source; worth front-loading into Commit 3

**Category:** QUALITY

**Description:**
The plan labels `T_can / Tinv_can` precomputation as "optional for
arithmetic parity" — but those rotations are INSIDE the per-step fault
loop at `wave_operator.inl:1087-1108` (RK4 interior) and equivalents for
shared / ADER. They are rebuilt from the SAME inputs each step, which is
deterministic run-to-run, but the compiler is free to reorder FP
operations in that reconstruction. If the canonical-frame rebuild itself
contributes ULP-scale drift (e.g., cross products computed differently
for + and − sides in the shared-face case), that is AN ADDITIONAL
per-step noise seed on top of the `Q_self(q) = Σ shape·Q` effect the
plan focuses on.

This is minor compared to the `sigma_n × eps_mach` source but nonzero,
and precomputing has the side benefit of eliminating it. I'd fold it
into Commit 3 rather than deferring to Commit 5.

**Suggested fix:** Annotate Phase 4 item 8 as "can be folded into
Commit 3 if the implementer can precompute rotations without breaking
the shared-face role/tiebreaker logic." If done, verify via byte-exact
output that `T_can` and `Tinv_can` are the same post-precompute as
they were per-iteration.

**Test case:** Not required for LOW.

---

### [R-007] [LOW] The plan's acceptance gate ("pepper reproducer") is not a formal test target — make it a CI gate

**Category:** QUALITY

**Description:**
Phase 5 item 12 mentions the pepper reproducer as "the real acceptance
test" but no concrete test-target name is proposed. Without a CI gate,
the next refactor could silently re-introduce pepper.

**Suggested fix:** Define a formal target:

```diff
# Makefile
+TEST_TPV102_PEPPER_SRC = tests/unit/test_tpv102_pepper_reproducer.cpp
+TEST_TPV102_PEPPER_OBJ = $(TEST_TPV102_PEPPER_SRC:.cpp=.o)
+seas_test_tpv102_pepper_reproducer: $(TEST_TPV102_PEPPER_OBJ) ...
+test-tpv102-pepper-reproducer: seas_test_tpv102_pepper_reproducer
+	./seas_test_tpv102_pepper_reproducer
```

and include `test-tpv102-pepper-reproducer` in both `test` and
`test-v93-regression` (or create `test-v94-regression`). The test must
FAIL on commit `632ad03` (pre-v9.3.0 total-Q migration — pepper baseline)
and PASS on HEAD after the v9.4.0 commits land.

**Test case:** This IS the test target; defined in the R-001 fix.

---

## 6. Open Questions / Deferred Items

1. **ADER sub-step nucleation**: deferred to v9.5.0 per R-002.
2. **Precomputed face interpolation matrices**: plan's Phase 4 item 9;
   deferred to v9.5.0.
3. **BP5 end-to-end re-validation after Commit 1**: not in scope here
   but must happen before merging commit 1 to master.
4. **Pepper reproducer design**: the hardest technical decision —
   exposing ULP-scale dispersion without running a full Frontera job is
   non-trivial. R-001 proposes a unit-scale reproducer but the mapping
   from "per-QP ULP seed" to "visible pepper" is through multi-step
   rupture; a longer serial reproducer (~100k steps on a 20-element
   mesh) may be necessary.

---

## 7. Commit Order (PROPOSED — aligned with the user's plan)

1. **Commit 1** (`fault_face_flux`): `Evaluate()` reads `tau*_nuc`; doc
   update. Gated by R-001 reproducer and R-003 BP5 invariant.
2. **Commit 2** (driver / setup): revert TPV102 init to fluctuation Q;
   rename `ApplyNucleationTotalPrestress` → `ApplyNucleationPrestress`;
   delete bulk-Q nucleation entry.
3. **Commit 3** (`wave_operator`): switch fault dispatch to
   `Evaluate` / `EvaluateADER`; relax `has_bulk_bg_` asserts. Gated by
   R-004 absorbing-BC regression.
4. **Commit 4** (tests): repair `test_persistent_nuc_prestress_channel`
   (`EvaluateTotal` → `Evaluate`), `test_fault_face_flux_total_vs_fluctuation`,
   both shared-fault parallel tests; delete the three dead bulk-Q
   nucleation tests (R-005); add the pepper reproducer (R-001, R-007).
5. **Commit 5** (optional parity): precomputed rotations / face
   kernels; deferred unless pepper is not closed by Commits 1–4.

---

## 8. Acceptance / Exit Criteria for v9.4.0

- `test-tpv102-pepper-reproducer` PASSES on HEAD (R-001, R-007).
- `make test-v93-regression` PASSES (all tests updated to match new contract).
- BP5 cycle recurrence within 1% of pre-v9.4.0 value (R-003).
- Frontera 400-rank 2.0 s rerun shows:
  - V_max trajectory within 15% of DRDG3D reference at t ∈ {0.2, 0.5, 0.8, 1.0, 1.5} s
    (looser than round-3 expectation — the split-prestress refactor may
    change wave-propagation phasing at the ULP level).
  - 6/9 SCEC fault stations register nontrivial slip by t = 2.0 s.
  - No pepper visible in the t = 2.0 s VTU fault surface.

---

## 9. What This Plan Does NOT Do

- Does not change the friction solver (Brent's method, regularised
  Dieterich-Ruina — all unchanged).
- Does not change the per-QP canonical frame or shared-face role logic
  (all v7.1.0 / v8.0.0 R-101 / R-801 fixes preserved).
- Does not revert the round-2 persistent-prestress channel semantics —
  only widens them to `Evaluate` (fluctuation path).
- Does not touch Frontera sbatch infrastructure beyond adding a new
  v9.4.0 sbatch when ready.

---

## 9a. R-001 Reproducer Results (2026-04-22, this session)

The R-001 acceptance reproducer
(`miniapps/seas/tests/unit/test_tpv102_pepper_reproducer.cpp`,
`make test-tpv102-pepper-reproducer`) was implemented, built, and
executed at the configuration used by the Frontera production run:
order-1 L2 + GaussLobatto on a tetrahedron face with Dunavant 3-point
face quadrature (`ndof_per_el = 4`, `nqp = 3`).

**Measured per-QP dispersion of `Q_self(q) = Σ_i shape_i(q) · C` for
uniform `Q_dofs = C`:**

| `C` (Pa) | min `Q_self` | max `Q_self` | `|Δ|` (Pa) | relative |
|----------|-------------:|-------------:|-----------:|---------:|
| 1.0e+06  | 1.000000000000000e+06 | 1.000000000000000e+06 | **0.000e+00** | 0.000e+00 |
| 1.0e+07  | 1.000000000000000e+07 | 1.000000000000000e+07 | **0.000e+00** | 0.000e+00 |
| 1.0e+08  | 9.999999999999999e+07 | 1.000000000000000e+08 | **1.490e-08** | 1.490e-16 |
| 1.0e+09  | 1.000000000000000e+09 | 1.000000000000000e+09 | **0.000e+00** | 0.000e+00 |

All three acceptance subtests PASS (`(a)` dispersion nonzero at at
least one `C`, `(b)` ULP-bounded, `(c)` linear scaling — trivially
satisfied since only one magnitude shows any dispersion), but the
**substantive finding is that the dispersion is NOT systematic**:

- At the TPV102 production scale `C ~ 1.2e8 Pa`, the dispersion is
  ~1 ULP (1.49e-8 Pa) across the 3 QPs.
- At the fluctuation-Q hypothetical scale `C ~ 1.0e6 Pa`, dispersion
  is **zero** (the split-prestress refactor would eliminate the seed
  amplitude, but only because the seed is effectively zero there — not
  because of a linear C-scaling reduction).
- At `C = 1e7 Pa` and `C = 1e9 Pa`, dispersion is also zero. The ULP
  mechanism only manifests at one specific magnitude, inconsistent with
  a robust forced-response seed.

**Interpretation for the v9.4.0 decision:**

The plan's ULP-scale systematic re-injection hypothesis is **weakly
supported but inconsistent** on the minimal fixture:

1. There IS per-QP dispersion at the total-Q scale (1.49e-8 Pa per
   step). Over ~3500 steps of a 1-second rupture ramp at
   `dt = 2.85e-4 s`, the cumulative seed worst-case is
   ~5.2e-5 Pa. For pepper to reach visibly ~kPa-scale station noise,
   the amplification factor would need to be ~10^7–10^8 through the
   DG + friction + nonlinear-rate-state coupling. That is POSSIBLE
   under high `psi` and high slip-rate conditions, but it is not
   self-evident.
2. Dispersion is NOT monotonically scaling with `C`. Three of four
   magnitudes produce bit-exact partition-of-unity results (dispersion
   = 0). The plan predicted linear scaling; reality shows FP
   quantization noise that happens only at specific magnitudes. This
   suggests the seed mechanism is NOT a reliable forced-response
   driver — it depends on the specific FP rounding of the shape-weight
   sum, which is only noisy at certain `C` values.
3. The SPLIT-PRESTRESS refactor would cut the seed from 1.49e-8 Pa per
   step (at C ~ 1.2e8 Pa) to 0 Pa per step (at C ~ 1e6 Pa). That IS
   the plan's claimed benefit, just driven by quantization floor
   collapse rather than by the linear-scaling argument.

**Revised recommendation (supersedes §10 below):**

Under the originally proposed R-001 gate, subtests (a)/(b)/(c) all
pass — formally, the refactor is unblocked. But the empirical finding
weakens the confidence of the plan's architectural argument:

- **The pepper seed is real but small** (1 ULP at C = 1.2e8 Pa).
- **The split-prestress refactor reduces the seed amplitude**, but the
  reduction is sporadic (dispersion goes to zero at C = 1e6 Pa
  because of FP quantization, not because of a clean
  proportional-to-C mechanism).
- **The plan should proceed** with Commit 1 (extend `Evaluate` to read
  `tau*_nuc`) — that commit is correctness-preserving and costs
  nothing even if the pepper root cause is elsewhere.
- **Before Commits 2–4** (full revert of total-Q), a SECOND gate
  should be added: run a short (1 s, order-1, 200 m mesh) local
  rupture simulation under the existing v9.3.0 + persistent-prestress
  code and probe whether per-QP pepper grows exponentially or
  linearly. Exponential growth under R-001's small seed would confirm
  the amplification hypothesis; linear growth would suggest the seed
  is too small to be the primary driver and another root cause is at
  play.

**Next action (suggested):** Add a second test
`test_tpv102_pepper_amplification.cpp` that runs ~100 ADER-O(2) steps
on the 2-tet fixture from `test_persistent_nuc_prestress_channel.cpp`,
under full-saturation nucleation, and probes the per-QP dispersion of
`data.tau2_corr` as a function of step number. Only with BOTH the
per-step seed test (R-001, closed here) AND the multi-step
amplification test passing does Commit 2 become justified.

## 10. Reviewer's Bottom Line (my opinion)

The plan is architecturally sound — it matches SeisSol's split-prestress
pattern, which is a well-validated reference, and the ULP scaling
argument gives a plausible quantitative handle on the pepper seed
(~10^5 amplitude reduction). The 5-commit structure is reasonable and
the first commit is appropriately small.

**But** — the core hypothesis (ULP-scale bulk-Q reconstruction noise is
the pepper seed) has not been directly measured. R-001 is the single
most important gate: **before Commit 1 lands, the pepper reproducer
test must demonstrate the seed mechanism, not just be designed to fail
under current code**. If R-001's reproducer does NOT show per-QP
dispersion scaling with `sigma_n`, the root cause is elsewhere and
this plan will not fix pepper.

The plan's other risks — BP5 invariant (R-003), absorbing-BC
regression (R-004), and the ADER timing ambiguity (R-002) — are
tractable with the specific hardening steps above.

**If R-001 closes, I recommend proceeding** with the 5-commit sequence
as written, with R-002 / R-003 / R-004 / R-005 / R-007 addressed inline
with each commit.

**If R-001 cannot close** (the ULP dispersion doesn't manifest or doesn't
scale with `sigma_n`), re-open the root-cause investigation: candidates
include (a) friction-solver conditioning at high psi + high sigma_n, (b)
a shared-fault DOFData-drift mechanism R-101 doesn't catch under
rupture, or (c) a canonical-frame rebuild FP-drift.

---

## 11. Guard Test for the Pepper Fix

Any candidate fix emerging from this plan MUST be validated against the
existing minimal pepper reproducer before it is considered landed:

**Test:** `tests/unit/test_adjacent_triangle_fault_uniformity.cpp`

**Build / run:**

```bash
make seas_test_adjacent_triangle_fault_uniformity
./seas_test_adjacent_triangle_fault_uniformity                                 # serial
$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_adjacent_triangle_fault_uniformity  # 4-rank
# or via Makefile targets:
make test-adjacent-triangle-fault-uniformity-serial
make test-adjacent-triangle-fault-uniformity-parallel
```

**Fixture:** 8 tetrahedra with a planar y=0 fault carrying 8 fault
triangles. Uniform TPV102 background state (`SYY = sigma_n_bg`,
`SXY = tau_ini`, zero dynamic fluctuations). No nucleation (rupture
drive off — this is the "forced-response" seed probe). Order-1 DG.
20 explicit ADER steps.

**Assertions (fault-QP spread after 20 steps):**

- `slip_rate` spread <= 1.0e-10
- `tau1_corr` spread <= 1.0e-10  (the channel that currently fails)
- `tau2_corr` spread <= 1.0e-10
- `sigma_n_corr` spread <= 1.0e-10

**Current status — 2026-04-23, serial, post-v9.4.0 Commits 1-3 +
REVIEW round-2 R-001 re-enablement (rupture drive restored per the
file-header §1 contract):**

| channel        | spread (Pa)   | tol     | status   |
|----------------|--------------:|--------:|----------|
| `slip_rate`    | **4.56e-05**  | 1e-10   | **FAIL** |
| `tau1_corr`    | **2.24**      | 1e-10   | **FAIL** |
| `tau2_corr`    | **1.74e-06**  | 1e-10   | **FAIL** |
| `sigma_n_corr` | **2.05e-06**  | 1e-10   | **FAIL** |

The v9.4.0 split-prestress refactor DOES NOT close this guard test
under rupture drive.  Reducing the bulk-Q seed amplitude (from
`~sigma_n * eps_mach ~= 3e-8` at C = 1.2e8 Pa to machine-epsilon at
C = 0) is necessary but not sufficient.  Rupture drive with uniform
`tau2_nuc = nuc_dtau` uncorks a per-tet shape-weight DOF asymmetry
in the wave operator that amplifies into visible pepper on all four
channels.  This matches the reviewer's R-005 analysis:
`wave.Mult(Q_bg, k)` shows a 3.85e-4 Pa/s conservation-law residual
(per-tet asymmetry) even with Q_bg zero; the nonlinear friction
feedback turns that residual into the observed pepper.

**Historical (pre-rupture-drive, silently disabled) status:**
Prior to the REVIEW round-2 R-001 fix, this test had
`tau2_nuc = 0` hardcoded and ran as an equilibrium conservation
check that passed at `tau1_corr` spread ≈ 1.93e-11.  The hardcoding
contradicted the file-header §1 declared contract and was removed
in the round-2 fix so the guard exercises rupture-drive pepper as
intended.

**Previous (pre-v9.4.0) status for reference, total-Q dispatch
at equilibrium:** `tau1_corr` spread 2.06e-10 Pa (failed by 2x);
other 3 channels ~1e-16 Pa (passed).  See REVIEW.md round-2 for
the full diagnostic chain.

**Gate for accepting a pepper fix:**

A candidate fix is accepted only if ALL FOUR channels of this test
drop below `1e-10` spread in BOTH the serial and 4-rank parallel runs,
under uniform nucleation drive (`tau2_nuc = nuc_dtau` on every fault
QP, as mandated by the file-header §1 contract).  Dropping below
ULP-of-background (~1e-7 Pa at `sigma_n = 1.2e8 Pa` / ~2e-10 Pa at
`tau = 7.5e7 Pa`) is ideal; the `1e-10` threshold is a looser
structural guard against pepper reappearing.

**Confirmed NOT sufficient as pepper fixes (2026-04-23):**

1. v9.4.0 split-prestress refactor (Commits 1-3) alone: drops bulk-Q
   seed amplitude to machine-epsilon but the rupture-driven guard
   still fails at `tau1_corr` ≈ 2.24 Pa.

   Source: per-tet shape-weight DOF asymmetry in the wave operator's
   face assembly (`wave_operator.inl` lift operator L(i,q) =
   w_q · shape1_q(i)).  Non-uniform L(i,q) across per-tet DOFs
   causes per-QP rhs deposition to differ under otherwise-identical
   QP-level friction outputs.  Nonlinear friction amplifies the
   per-step delta into visible pepper.

**Still-open candidate architectural fixes (beyond v9.4.0):**

- Per-face DR batching with per-QP physics (SeisSol
  `BaseFrictionLaw::evaluate` pattern): keeps per-QP friction physics
  intact but makes face-level rhs deposition symmetric by construction.
- Higher-order DG (p ≥ 2) enforces polynomial smoothness across DOFs
  within a tet; does not eliminate the asymmetry but reduces it.
- Lumped-mass L2-projection of face flux (averages per-QP contributions
  into a single face amplitude before redistribution).

Each of these is a separate commit sequence beyond v9.4.0's split-
prestress scope.  Documentation / implementation owed as v9.5.0+.

This test serves as the regression guard once any pepper fix lands —
it is cheap (8 tets, 20 steps, seconds to run) and must stay green
across all future DG / fault-flux / friction-solver refactors.  Until
it DOES stay green, no production (400-rank Frontera) rupture run
should be submitted — the fault-surface output would reproduce the
pepper regardless of Frontera-scale investment.

**Companion diagnostic tests** (do NOT gate the fix, but narrow the
analysis if the guard fails):

- `test_pepper_amplification_chain.cpp` — T1-T4 drill-down from
  `Q_self` per-QP drift to friction/bulk response.
- `test_face_lift_operator.cpp` — structural `L(i,q) = w_q · shape1_q(i)`
  asymmetry showing which DOFs are most sensitive to per-QP noise.
- `test_pepper_qp_perturbation_response.cpp` — prescribed-QP probe via
  the 3x3 face interpolation `S`-solve (isolates single-QP response).
- `test_pepper_nonlinear_growth.cpp` — delta-only L1 + eps-ladder +
  3-regime (OFF / WEAK / STRONG nucleation) multistep growth-rate
  measurement. Documented baseline: `gamma_dQ ~ 0.998` (DECAY, not
  amplification from a one-shot seed) — consistent with a per-step
  re-injection rather than eigenmode growth.
- `test_tpv102_pepper_reproducer.cpp` — R-001 acceptance reproducer
  (per-QP `Q_self` dispersion vs magnitude `C`), referenced in §9a.
