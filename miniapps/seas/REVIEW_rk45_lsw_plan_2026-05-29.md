# Code Review: PLAN_rk45_lsw_mixed_flux_2026-05-29.md (fresh adversarial review, 2026-05-29)

> Note: the prior Phase-13 review previously in this file is archived at
> `REVIEW_phase13_round2_2026-05-29.md`. This REVIEW.md now reviews the RK45+LSW plan.

## Review Scope
- Plan reviewed: `miniapps/seas/document/mixed_flux_dev/PLAN_rk45_lsw_mixed_flux_2026-05-29.md`
  (a planning document — no implementation code exists yet).
- Code cross-checked against the plan's claims:
  - `dynamic/fault_face_flux.{hpp,cpp}` — `Evaluate`, `EvaluateADER_LSW`, `WriteBackState`, `BuildImposedState`
  - `dynamic/tpv205_friction.hpp` — `LSWFrictionCoefficient_TPV205`, `SolveLSW_TPV205`
  - `dynamic/rk_time_stepper.{hpp,cpp}` — `AdvanceRKCoupled_Spatial` (RS sibling)
  - `dynamic/wave_operator.{hpp,inl}` — `FaultFrictionLaw` enum, Mult-path fault dispatch (`:2617` interior, `:3260` shared)
  - `dynamic/nucleation_method.hpp` — `StaticOverstress`/`InstantaneousOverstressCircular` `ApplyAbsolute`
  - `drivers/spatial_dyn_driver.cpp` — guards (`:771`, `:778`), tableau/dispatch (`:2568`, `:2750`), `slip_rate_substep_max` reset (`:2737`), `SetFaultFrictionLaw`/`SetMixedFluxMode`/`SetCflRkAware`
  - `spatial/code/spatial_friction.cpp:1101` — TOML `time_integrator` parse
- Domain context consulted: `miniapps/seas/CLAUDE.md` (sign/frame conventions, "no local full-mesh runs"), the companion `BUILD_mixed_flux_rk_dynamic_rupture_2026-05-28.md`, project memory (mixed-flux/RK scope, byte-exact contract).

**Overall:** the plan's factual claims (line numbers, signatures, the dispatch architecture, the LSW physics derivation, the no-op-nucleation and law-agnostic-CFL claims, TOML `time_integrator` parsing) all check out against the code. If the prose is followed *literally*, the resulting code is correct. The findings below are (a) one acceptance criterion that is numerically false as written, (b) two porting traps the plan under-specifies given that `/code-implement` + `/code-fix` work mechanically off the RS sibling, and (c) a coverage gap: the plan modifies a parallel-only code path but provides no test that exercises it. No CRITICAL findings (no specified step produces wrong results when followed exactly).

---

## Findings

### [R-001] [MODERATE] [Phase 1 acceptance criterion / test L1] — "bit-identical for any dt" is numerically false (round-trip ULP gap)

**Category:** BUG (test design — the specified acceptance test will fail and send the implementer chasing a non-bug)

**Description:**
The Phase 1 acceptance criterion and test L1 assert that
`EvaluateLSW(data, Q±)` writes DOFData **bit-identical** to
`EvaluateADER_LSW(data, Q±·dt, dt)` "for any `dt > 0`" / "for several `dt`".
This is false in floating point. `EvaluateADER_LSW` reconstructs
`Q̄ = I·(1/dt) = (Q·dt)·(1/dt)` (`fault_face_flux.cpp:788-792`). For a general
`dt` (e.g. `0.1`), `(Q·dt)·(1/dt) ≠ Q` at the ~1-ULP level, and that gap
propagates through `ComputeTrialTraction` → `SolveLSW_TPV205` → `WriteBackState`,
so the DOFData fields differ by a few ULP. The round-trip is exact **only** when
`dt` is a power of two (including `1.0`, `0.5`, `0.25`).

**Trigger:** L1 run with `dt = 0.1` (or any non-power-of-two), asserting exact
equality.

**Actual behavior (as the test would behave):** `TEST_TRUE(a == b)` fails on the
last 1–2 ULP for the trial/solve-derived fields, even though `EvaluateLSW` is
correct.

**Expected behavior:** the equivalence should be asserted bit-exactly only for
`dt ∈ {1.0, 0.5, 0.25}` (powers of two), and with a tight relative tolerance
(≤ 1e-12) for general `dt`.

**Suggested fix (edit the plan):**
```diff
 ### Acceptance Criteria   (Phase 1)
-- [ ] `EvaluateLSW(data, Q±)` writes `data.{V1,V2,slip_rate,tau1_corr,tau2_corr,
-      sigma_n_corr}` **bit-identical** to `EvaluateADER_LSW(data, Q±·dt, dt)` for
-      any `dt > 0` ... (Exact, not O(dt²): ...)
+- [ ] For `dt ∈ {1.0, 0.5, 0.25}` (powers of two — the I/dt round-trip is exact):
+      `EvaluateLSW(data, Q±)` writes `data.{V1,V2,slip_rate,tau1_corr,tau2_corr,
+      sigma_n_corr}` **bit-identical** to `EvaluateADER_LSW(data, Q±·dt, dt)`, and
+      `Q_imp_± == I_imp_±/dt` exactly.
+- [ ] For a general `dt` (e.g. 0.1): the same fields match within relative
+      tolerance 1e-12 (the `(Q·dt)·(1/dt)` reconstruction in EvaluateADER_LSW
+      carries a ~1-ULP round-trip gap that propagates through trial→solve).
```
And the L1 test-matrix row: replace "for several `dt`" with "for `dt ∈ {1.0,0.5,0.25}` (exact) and `dt=0.1` (rtol 1e-12)".

**Test case (C++ sketch, `tests/unit/test_lsw_rk_mixed_flux.cpp`):**
```cpp
// L1: EvaluateLSW vs EvaluateADER_LSW round-trip
DOFData d = MakeTPV205InitDOF();          // lsw_* + impedances set, slip1/2 nonzero
real_t Qp[NUM_STATE], Qm[NUM_STATE];      // representative fault-local states
for (real_t dt : {1.0, 0.5, 0.25}) {      // powers of two => exact
   DOFData a = d, b = d;
   real_t qip[NUM_STATE], qim[NUM_STATE], iip[NUM_STATE], iim[NUM_STATE], Ip[NUM_STATE], Im[NUM_STATE];
   for (int c=0;c<NUM_STATE;c++){ Ip[c]=Qp[c]*dt; Im[c]=Qm[c]*dt; }
   flux.EvaluateLSW(a, Qp, Qm, qip, qim);
   flux.EvaluateADER_LSW(b, Ip, Im, dt, iip, iim);
   TEST_TRUE(a.V1==b.V1 && a.V2==b.V2 && a.slip_rate==b.slip_rate
             && a.tau1_corr==b.tau1_corr && a.tau2_corr==b.tau2_corr
             && a.sigma_n_corr==b.sigma_n_corr, "L1 bit-exact (pow2 dt)");
   for (int c=0;c<NUM_STATE;c++)
      TEST_TRUE(qip[c]==iip[c]/dt, "L1 Q_imp == I_imp/dt (pow2 dt)");
}
// dt=0.1: tolerance-only
{ DOFData a=d,b=d; /* ... */ TEST_NEAR(a.V1, b.V1, 1e-12*std::abs(b.V1)); }
```

---

### [R-002] [MODERATE] [Phase 3 — AdvanceRKCoupledLSW_Spatial slip combine] — `+=` vs `=` porting trap not called out; the RS sibling the plan says to model on uses `+=`

**Category:** ASSUMPTION (the plan is correct but under-specifies the single most likely mechanical-port bug)

**Description:**
Phase 3 instructs the implementer to build `AdvanceRKCoupledLSW_Spatial`
"structurally parallel to `AdvanceRKCoupled_Spatial`". That RS sibling accumulates
slip with **`+=`** at `rk_time_stepper.hpp:283-284`:
```cpp
dof_data[m].slip1 += slip1_inc;   // RS: slip is NOT staged, so dof_data.slip1 == slip1_n here
```
In the RS path that is correct because RS `Evaluate` never reads slip, so
`dof_data.slip1` is unchanged from `slip1_n` at combine time. In the **LSW** path,
the stepper writes the stage-local slip into `dof_data` before every `Mult`
(required, so `EvaluateLSW` sees the staged δ), so at combine time
`dof_data[m].slip1` holds the **last stage's** staged value, not `slip1_n`. A
mechanical port that keeps `+=` therefore computes
`slip_last_stage + dt·Σ b_i V_k` — double-counting the partial stage sum — a silent
wrong-front-speed bug with no crash. The plan's Phase 3 step 2 *does* write the
combine as `slip1 = slip1_n[m] + …` (correct), but it never explicitly says "do
NOT keep the RS sibling's `+=`", and the L6 test as worded ("constant V")
under-specifies the regime needed to catch it.

**Trigger:** implementer copies `rk_time_stepper.hpp:283` (`+=`) into the LSW
stepper while also adding the stage-local slip write.

**Actual behavior:** final slip is corrupted (last-stage staged slip + full
b-weighted increment); rupture front accelerates; no abort.

**Expected behavior:** combine MUST be `dof_data[m].slip1 = slip1_n[m] + dt·Σ_i b_i V1_k[i][m]`
(absolute from the step-start snapshot), and the regression test must assert the
**absolute** post-step slip (not an increment) in a regime where μ(δ) is constant
(δ ≥ d_c) so the closed-form expected value is exact.

**Suggested fix (edit the plan):**
Add to Phase 3 step 3 (the "slip staging is the only structural delta" paragraph):
```diff
+ **Combine must overwrite, not accumulate.** Because the stage loop overwrites
+ `dof_data[m].slip1/slip2` every stage, at combine time those fields hold the LAST
+ stage's staged slip — NOT `slip*_n`.  The combine therefore writes the ABSOLUTE
+ value `dof_data[m].slip1 = slip1_n[m] + dt·Σ_i b_i V1_k[i][m]` (`=`, from the
+ snapshot).  Do NOT copy the RS sibling's `dof_data[m].slip1 += slip1_inc`
+ (rk_time_stepper.hpp:283): in RS slip is unstaged so `+=` is correct, but in LSW
+ `+=` double-counts the last stage's partial sum.
```
And strengthen L6:
```diff
- [ ] **Slip accumulation:** ... frozen, prescribed `V` ... slip*_new == slip*_n + V*·dt
+ - [ ] **Slip accumulation (absolute, not increment):** single-QP LSW QP driven so
+       that δ ≥ d_c (μ = μ_d constant) and Q± held fixed across the step ⇒ V is
+       constant; assert the ABSOLUTE post-step `dof_data.slip1` equals
+       `slip1_n + V1·dt` to round-off (Σb=1 makes RK4-of-a-constant exact).  This
+       fails loudly if the combine uses `+=` instead of `=`.
```

**Test case (C++ sketch):**
```cpp
// L6: absolute slip after one step, constant-mu regime (catches += bug)
DOFData d = MakeTPV205InitDOF();
d.lsw_d_c = 0.1; d.slip1 = 0.5; d.slip2 = 0.0;  // delta=0.5 >= d_c => mu=mu_d const
std::vector<DOFData> dof{d};
Vector Q = MakeFixedFaultState();               // held constant => constant V
real_t dt = 0.01, t = 0.0; Vector Qn(Q.Size());
const real_t slip1_n = dof[0].slip1;
AdvanceRKCoupledLSW_Spatial(wave_LSW, dof, Q, dt, t, Qn, MakeRK4Tableau(), nullptr);
const real_t V1 = dof[0].V1;                      // endpoint V1 (constant across step)
TEST_NEAR(dof[0].slip1, slip1_n + V1*dt, 1e-14*std::abs(slip1_n + V1*dt));
// A '+=' combine yields slip1 ~ slip1_n + partial + V1*dt  != expected.
```

---

### [R-003] [MODERATE] [POSSIBLE] [Phase 3/Phase 2 — stepper⇄fault_friction_law_ coupling] — asymmetric guard lets a mismatched (stepper, law-flag) pair silently run the wrong fault kernel

**Category:** BUG (robustness / silent-wrong on misuse)

**Description:**
After Phase 2, `WaveOperator::Mult` selects the fault kernel purely from
`fault_friction_law_` (LSW → `EvaluateLSW`, RS → `Evaluate`). The *stepper* and the
*flag* are two independent settings that MUST agree:
- LSW stepper + `RateAndState` flag → `Mult` runs RS `Evaluate` (reads `data.psi`/`data.b`,
  which LSW DOFData leaves at defaults → NaN/garbage) while the stepper stages slip.
- RS stepper + `LSW` flag → `Mult` runs `EvaluateLSW` while the stepper integrates ψ
  and never stages slip → δ frozen at init → wrong physics.

The plan adds a self-check **only** to the LSW stepper (Phase 3 step 2:
`MFEM_VERIFY(wave.GetFaultFrictionLaw() == LSW)`). The RS stepper
`AdvanceRKCoupled_Spatial` has **no** such guard (confirmed: no `GetFaultFrictionLaw`
reference anywhere in `rk_time_stepper.{hpp,cpp}`), so the RS-stepper-with-LSW-flag
mismatch is unguarded. In the production driver the two are wired from a single
`is_lsw` (and cross-checked at `spatial_dyn_driver.cpp:2518-2527`), so this bites
mainly a future driver refactor or a miswired unit test — but the plan introduces
the asymmetry, and the fix is free (additive `MFEM_VERIFY`, no numerics change,
byte-exact preserved).

**Trigger:** any caller (e.g. a new unit test, or a future code path) that invokes
`AdvanceRKCoupled_Spatial` on a `WaveOperator` whose `fault_friction_law_ == LSW`.

**Actual behavior:** no abort; `Mult` runs `EvaluateLSW` with a frozen δ; the RS
stepper integrates a ψ that the fault solve ignores → silently wrong results.

**Expected behavior:** abort with a clear message, symmetric to the LSW stepper's
guard.

**Suggested fix (edit the plan — add to Phase 3 requirements):**
```diff
+ **Symmetric guard on the RS stepper.** Because Phase 2 makes the Mult fault
+ kernel a pure function of `fault_friction_law_`, add the mirror self-check to the
+ EXISTING `AdvanceRKCoupled_Spatial` (RS):
+   MFEM_VERIFY(wave.GetFaultFrictionLaw() == FaultFrictionLaw::RateAndState,
+               "AdvanceRKCoupled_Spatial integrates psi, but the WaveOperator's "
+               "fault_friction_law_ is not RateAndState; Mult would run the wrong "
+               "fault kernel.  Use AdvanceRKCoupledLSW_Spatial for LSW.");
+ This is additive (an MFEM_VERIFY, no arithmetic change) so the byte-exact RS-RK
+ contract is preserved; it just makes the (stepper, flag) mismatch fail loud.
```

**Test case:**
```cpp
// R-003: mismatched stepper/flag aborts (death test)
WaveOperator<Mesh> wave(/*...*/);
wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW);
std::vector<DOFData> dof = MakeTPV205InitDOFs(/*one fault QP*/);
EXPECT_MFEM_ABORT( AdvanceRKCoupled_Spatial(wave, dof, rs_cfg, rs, Q, dt, t, Qn,
                                            MakeRK4Tableau(), nullptr) );
wave.SetFaultFrictionLaw(FaultFrictionLaw::RateAndState);
EXPECT_MFEM_ABORT( AdvanceRKCoupledLSW_Spatial(wave, dof, Q, dt, t, Qn,
                                               MakeRK4Tableau(), nullptr) );
```

---

### [R-004] [MODERATE] [Phase 2 — shared-fault Mult dispatch / Testing Strategy] — the plan modifies a parallel-only code path (`wave_operator.inl:3260`) but no listed test exercises it

**Category:** EDGE_CASE (coverage gap on a code path the plan explicitly edits)

**Description:**
Phase 2 step 2 edits the **shared-fault** Mult dispatch at
`wave_operator.inl:3260` (the `ComputeSharedFaceFluxRHS` path, reached only on
`ParMesh` with a fault crossing a rank seam). Every unit test in the plan's L1–L7
matrix is **serial** (`Mesh`, single rank) and exercises only the interior-fault
path (`:2617`). A mistake in the `:3260` edit — wrong variable, an inverted branch,
or simply forgetting to add the LSW arm there while adding it at `:2617` — would
compile, pass all of L1–L7, and silently run the **RS** `Evaluate` on LSW DOFData
across rank seams (→ NaN via `data.b`/`data.psi`, or a desynced rupture front).
The MPI shared-fault check appears only in the Risk section ("add an MPI test
before production"), not as a numbered acceptance criterion or phase deliverable,
so a mechanical implementer can mark Phase 2 "done" with the parallel arm untested.

**Trigger:** any np>1 run with a fault on a rank boundary (i.e. every production
TPV205 run); not caught by `make test` serial unit coverage.

**Actual behavior (if the `:3260` arm is wrong/missing):** shared-fault QPs run the
RS kernel on LSW data → NaN or cross-rank rupture desync; serial tests stay green.

**Expected behavior:** an MPI test asserts that the shared-fault LSW dispatch
matches the interior-fault LSW dispatch (a fault QP on a seam gets the same
V/τ_corr as the same QP would interior).

**Suggested fix (edit the plan):** promote the MPI test from Risk #2 to a Phase 2
acceptance criterion and add an L8 row:
```diff
 ### Acceptance Criteria   (Phase 2)
+- [ ] **MPI shared-fault parity:** a 2-rank run with the fault crossing the rank
+      seam, `fault_friction_law_ == LSW`, produces V/τ*_corr/σ_n_corr on the shared
+      QPs equal (to round-off) to the values the same QPs get on a single-rank run
+      (interior path).  Guards the `:3260` edit, which no serial test reaches.
```
```diff
 ## Testing Strategy  (test matrix)
+| L8 | 2 | MPI (np=2) fault-on-seam: shared-fault LSW dispatch (:3260) gives the
+         same fault observables as the serial interior path (:2617). Run via
+         `mpirun -np 2 seas_test_lsw_rk_shared_fault_mpi`. |
```

**Test case (MPI, mirrors existing `seas_test_fault_*_mpi` harness):**
```cpp
// L8: serial vs 2-rank shared-fault LSW parity at one seam QP
// rank-collective: build the same fault QP interior (serial ref) and on a seam
// (np=2); SetFaultFrictionLaw(LSW); one Mult; MPI_Allreduce the seam QP's
// {V1,V2,tau1_corr,tau2_corr,sigma_n_corr}; assert == serial reference to 1e-12.
```

---

### [R-005] [LOW] [Phase 5 — config route] — record that the config-only `time_integrator="rk45"` route is valid (avoid a CLI-only assumption)

**Category:** QUALITY (minor — prevents an unnecessary CLI-only assumption)

**Description:**
Phase 5 hedges between a TOML `[numerics] time_integrator = "rk45"` and CLI
pinning, recommending CLI. That recommendation is fine, but the plan does not state
that the TOML route is actually wired — it is (`spatial/code/spatial_friction.cpp:1101`,
`toml_str(n, "time_integrator", "ader")`, default ADER, same `ader|rk4|rk45`
vocabulary as the CLI). Worth recording so the implementer/reviewer does not later
"discover" the TOML key is unsupported (it is) or, conversely, add the key to a
config and assume it is a no-op.

**Trigger:** N/A (documentation clarity).

**Actual/Expected:** add one sentence confirming the TOML key is parsed.

**Suggested fix (edit the plan, Phase 5 Files to Create note):**
```diff
+ NOTE: `[numerics].time_integrator` IS a parsed TOML key
+ (spatial_friction.cpp:1101, default "ader"), so the config-only route is valid;
+ CLI pinning is recommended only for the TPV31-style "scheme printed in the
+ banner + asserted config-only knobs" parity, not because the TOML key is missing.
```

(No test — documentation only.)

---

## Summary
- Critical issues: 0
- Moderate issues: 4 (R-001, R-002, R-003, R-004)
- Low issues: 1 (R-005)
- Plan compliance: N/A (this is a plan, not an implementation). The plan's factual
  claims about the codebase are **accurate** (verified line numbers, signatures,
  dispatch architecture, no-op nucleation, law-agnostic CFL, TOML parsing). The
  findings are corrections to two acceptance criteria, two porting-trap
  clarifications, and one coverage gap.
- Verdict: **PASS WITH FIXES** — the plan is sound and implementable for the happy
  path; apply R-001 (false bit-exact criterion), R-002 (`+=`/`=` trap + stronger
  slip test), R-003 (symmetric RS-stepper guard), and R-004 (MPI shared-fault test)
  before `/code-implement`, so the implementer is not led into a failing test
  (R-001) or a silent-wrong port (R-002/R-003) and does not ship an untested
  parallel path (R-004).

## Unreviewed Areas
- **CFL numerical stability for RK+central-flux LSW at p1** (Phase 5): the plan
  defers calibration to a Frontera dev smoke ("starting calibration numbers"),
  consistent with the project's "no local full-mesh runs" rule. The actual stable
  `dt` cannot be reviewed statically and is an empirical sign-off, not a plan defect.
- **SCEC TPV205 physics correctness of the end-to-end run** (Phase 5 final AC): a
  post-merge Frontera task against the benchmark overlays; out of scope for a
  static review of the plan.
- **`EvaluateLSW_ForcedRupture` / TPV26-27** and **bimaterial LSW**: explicitly
  out of scope per the plan's "Out of scope" section; the `Mult`-path abort on
  `LSW_ForcedRupture` (Phase 2) is the correct guard and was reviewed.
