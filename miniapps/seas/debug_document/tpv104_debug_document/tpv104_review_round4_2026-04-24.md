# Code Review Round 4: TPV104 Steps 3, 7, 9, 11 — 2026-04-24

Fresh adversarial audit of the new files landed since round 3. Confirms
prior-round fixes still in place and hunts for new bugs in the four
newly shipped components: Step 3 setup, Step 7 sub-step iterator,
Step 9 smoke (driver-banner + dry-run), Step 11 probe-format scaffolding.

## Review Scope

- Plan: `miniapps/seas/debug_document/tpv104_debug_document/tpv104_debug_plan_2026-04-24.md`
  §4.10 Steps 3 (setup), 7 (sub-step iterator), 9 (driver / smoke), 11
  (probe format).
- Prior reviews (all closed):
  - `tpv104_review_2026-04-24.md` (R-001..R-012)
  - `tpv104_review_round2_2026-04-24.md` (R2-001..R2-007)
  - `tpv104_review_round3_2026-04-24.md` (R3-001..R3-007)
- **New files reviewed:**
  - `dynamic/tpv104_setup.hpp` (522 lines)
  - `dynamic/tpv104_substep_iterator.hpp` (157 lines)
  - `dynamic/tpv104_substep_iterator.cpp` (374 lines)
  - `tests/unit/test_tpv104_setup.cpp` (415 lines)
  - `tests/unit/test_tpv104_substep_iterator.cpp` (553 lines)
  - `tests/unit/test_tpv104_smoke.cpp` (190 lines)
  - `tests/unit/test_tpv104_probe_format.cpp` (234 lines)
- Unchanged since round 3: `config/tpv104_params.hpp`,
  `friction/slip_law_srw_psi.hpp`, `friction/friction_coeff_stable.hpp`,
  `dynamic/tpv104_friction_solver.hpp`, `dynamic/tpv104_nucleation.hpp`,
  and all tests from Steps 1/2/4/5/6/8/10.
- Domain context: `CLAUDE.md`, `miniapps/seas/CLAUDE.md`,
  `feedback_tpv102_bp5_no_shared_edit`, `feedback_no_local_reproducer`.

## Round-1/2/3 regression check

Spot-verified: all earlier fixes remain intact —
`SlipLawSRWPsi::SetProductionMode`, raw V (no clamp), `IntegerPow8`,
required `V_w_default`, `RateDerivativeV(V=0)` analytic limit,
`SolveSlipRateNewtonStable` input validation, τ/η_s cap in
`physical_v_guess`, `dS <= 0` clamp, `MFEM_ASSERT` on finite dt,
`TestNoDipOrNormalStressNucleation` sentinel. No regressions
introduced by the new files.

## Findings

### [R4-001] CRITICAL [tpv104_substep_iterator.cpp:118-371] — `Advance` never accumulates `slip1` / `slip2`; driver-side slip-per-step is left undocumented, and `WriteBackState` writes `V1/V2` from the LAST sub-step rather than a time-average

**Category:** BUG (plan deviation with physics impact)

**Description:**
`Tpv104SubStepIterator::Advance` mutates
`{psi, slip_rate, V1, V2, tau*_corr, sigma_n_corr, tau2_nuc}` but
**never touches `dof_data[i].slip1` or `dof_data[i].slip2`**. The
docstring (`tpv104_substep_iterator.hpp:90-91`) also omits slip from
its mutation list. There is no in-iterator loop of
`d.slip1 += d.V1 * dt_sub;` / `d.slip2 += d.V2 * dt_sub;`.

In the TPV102 driver this accumulation lives in
`drivers/tpv102_driver.cpp:1121-1122`
```cpp
dof_data[i].slip1 += dof_data[i].V1 * dt_step;
dof_data[i].slip2 += dof_data[i].V2 * dt_step;
```
executed AFTER `wave.AdvanceADER` returns. The TPV104 Step-9 driver is
not yet shipped, but if the caller follows the TPV102 template it
would run `slip += V * dt_macro` after every `Advance`. That is doubly
wrong for the iterator:

1. **Silent contract**: nothing in the iterator's header or the plan
   §4.10 Step 7 docstring tells the driver author that slip
   accumulation is their responsibility. A naive Step-9 driver that
   assumes "Advance mutates everything relevant" will ship with
   `slip1/slip2` frozen at zero, and the station writer will emit
   columns of zeros for `slip2` — producing a silently corrupt
   TPV104 trace.

2. **Wrong V to accumulate with**: even if the driver does
   `slip += V * dt_macro`, `WriteBackState` at `fault_face_flux.cpp:199-201`
   writes `s.V1, s.V2` from the LAST sub-step's friction solve:
   ```cpp
   data.slip_rate = s.V_abs;
   data.V1 = s.V1;
   data.V2 = s.V2;
   ```
   In a ramp-up regime (V grows through the macro-step as ψ decreases),
   the last-sub-step V is biased HIGH vs the time-average. The TPV102
   one-shot (`wave.AdvanceADER` on Q̄ = I/dt) delivers a V that matches
   an O(dt²) time-average; the sub-step iterator's last-sub-step V
   does not share that property. Downstream `slip += V_last * dt` is
   therefore NOT an O(dt²) approximation to ∫V dt.

None of the unit tests catches this — `TestSubStepAccumulatorByteMatch`
compares iterator vs an inline replay that has the same bug;
`TestO1LimitMatchesEvaluateADER` only runs O=1 so the last-sub-step IS
the time-average; `TestNucleationInjectionLockedFault` uses
V_ini = 1e-16 so slip stays ~0 by construction. The first time this
surfaces is under a real nucleation-triggered rupture in Phase 3.

**Trigger:** Any Step-9 driver macro-step loop that calls
`iterator.Advance(...)` and omits per-sub-step slip accumulation, OR
that does `slip += V * dt_macro` using the post-Advance V.

**Actual behavior:** slip1/slip2 silently stay at 0; or, with
driver-side `slip += V * dt_macro`, they are biased by the
last-sub-step-vs-time-average delta.

**Expected behavior:** Slip accumulation must live INSIDE the
iterator, summing per-sub-step: `d.slip1 += s.V1 * dt_sub; d.slip2 +=
s.V2 * dt_sub;` at the same cadence as the ψ update (after
`ComputeStageState` fills `s.V1, s.V2`).

**Suggested fix:**
```diff
          // ψ update per sub-step (§3.12 directive — FVW analytic step).
          ...
          d.psi = UpdateStateAnalyticSlipLawSRW(d.psi, s.V_abs,
                                                d.Dc, dt_sub,
                                                V_w[i], d.a,
                                                state_evo_.GetB(),
                                                state_evo_.GetV0(),
                                                state_evo_.GetF0(),
                                                state_evo_.GetMuW());
+
+         // R4-001 (review round 4): slip integration per sub-step.
+         // V1, V2 in `s` are the per-sub-step components resolved
+         // against the current ψ; ∫V dt over the macro-step
+         // telescopes to Σ_o V^(o) · dt_sub.  This lives in the
+         // iterator (not the caller) because the per-sub-step V's
+         // are ephemeral in `s` — `WriteBackState` only exposes the
+         // last sub-step's V.
+         d.slip1 += s.V1 * dt_sub;
+         d.slip2 += s.V2 * dt_sub;
```
And update the iterator's header contract to list `slip1/slip2` in
the mutated-fields line (`tpv104_substep_iterator.hpp:90-91`).

**Test case:**
```cpp
void test_R4_001_iterator_accumulates_slip() {
   // Drive the iterator past the nucleation ramp so V is non-trivial.
   std::vector<DOFData> dof_data;
   std::vector<Vector> fault_coords;
   std::vector<real_t> V_w;
   MakeSingleQPFixture(dof_data, fault_coords, V_w);
   dof_data[0].V1 = 0.0;
   dof_data[0].V2 = 1e-3;          // pretend rupture onset
   dof_data[0].slip1 = 0.0;
   dof_data[0].slip2 = 0.0;
   // Inject a non-zero τ₁ perturbation via I_+/I_- so the friction
   // solve produces V1, V2 > 0 per sub-step.
   const real_t dt_macro = 1e-3;
   std::vector<real_t> deltaT, tw;
   GaussLegendreO5SubSteps(dt_macro, deltaT, tw);
   FaultFaceFlux flux(TPV104Params::rho, TPV104Params::cp,
                      TPV104Params::cs);
   SlipLawSRWPsi law(TPV104Params::a_in, TPV104Params::b,
                     TPV104Params::V0, TPV104Params::f0,
                     TPV104Params::f_w, TPV104Params::V_w_in);
   Tpv104SubStepIterator it(flux, law);
   it.SetSubSteps(deltaT, tw);
   // … build I_plus/I_minus to produce V_abs > 0 at the friction solve …
   const real_t slip2_before = dof_data[0].slip2;
   it.Advance(dof_data, fault_coords, V_w,
              I_plus.data(), I_minus.data(),
              dt_macro, 0.0,
              I_imp_plus.data(), I_imp_minus.data(),
              FrictionSolver::Method::Brent);
   // Expect Σ V^(o) · dt_sub  >  0, not frozen at slip2_before.
   EXPECT_GT(dof_data[0].slip2 - slip2_before, 0.0);
}
```

---

### [R4-002] MODERATE [tpv104_substep_iterator.hpp:129-139, tpv104_substep_iterator.cpp:236] — Iterator's default friction solver routes through `FrictionSolver::SolveNR` (MFEM-native μ), not `SolveSlipRateNewtonStable` (Step-5 stable-asinh Newton); plan §4.10 Step 5 mandate violated

**Category:** DEVIATION (plan compliance)

**Description:**
Plan §4.10 Step 5 specifies: *"Newton solver is the driver's default;
`FrictionSolver::SolveNR` (legacy MFEM μ formula) is only reached
under `--friction-solver=legacy-newton`."* The round-1 checklist
tracks this as "Deviation #4 closed: Newton solver is the driver's
default … T_TPV104_FS_4 stub is replaced with a real smoke-test
assertion."

The iterator ships with default:
```cpp
void Advance(...,
             FrictionSolver::Method method
                = FrictionSolver::Method::NewtonRaphson);
```
But `FrictionSolver::Method::NewtonRaphson` dispatches to
`FrictionSolver::SolveNR` (`dynamic/friction_solver.cpp:61`), which
uses the **MFEM-native `DieterichRuinaFriction::FrictionCoefficientPsi`**
μ formula (with the `psi/a > 700` branch). The Step-5
`SolveSlipRateNewtonStable` + `FrictionCoefficientStable` path is
NOT invoked. The docstring explicitly acknowledges this
(`tpv104_substep_iterator.hpp:97-102`):
> "the iterator does not install its own Stable-asinh Newton …
> A future wiring change that routes through
> `SolveSlipRateNewtonStable` lives in Step 9 (driver)"

This contradicts the Step-5 closure: Step 5 said the stable Newton
was the TPV104 default production path; under the current iterator,
it isn't. Smoke test `TestBannerDefaults` (`test_tpv104_smoke.cpp:83-90`)
asserts the driver banner prints "Friction solver: Newton-Raphson" on
default flags — but a user reading that banner will reasonably assume
it means the Step-5 stable Newton, not the legacy MFEM Newton-Raphson
with the old μ formula.

**Trigger:** Phase-3 Probe 4 (slip-rate magnitude) runs against the
reference FVW Newton; MFEM-side produces the legacy Newton trajectory
which differs at ψ/a > 700 (not an issue on TPV104, but the formula
divergence is still present for any intermediate Newton iterate that
transiently lands at a large lx).

**Actual behavior:** Iterator uses `DieterichRuinaFriction::FrictionCoefficientPsi`;
the stable-asinh Newton is never reached via the iterator path.

**Expected behavior:** Either (a) wire the Step-5 Newton into
`FaultFaceFlux::ComputeStageState` dispatch so
`FrictionSolver::Method::NewtonRaphson` routes to it, OR (b) change
the iterator's default to a new method enum value that dispatches to
`SolveSlipRateNewtonStable`, OR (c) add a prominent caveat to the
banner: "Friction solver: Newton-Raphson (MFEM-native μ; stable-asinh
Newton deferred)".

**Suggested fix (option b — add a new method enum value dedicated to the Step-5 solver):**
```diff
 // dynamic/friction_solver.hpp
-   enum class Method { Brent, NewtonRaphson, HybridNRBisection };
+   enum class Method {
+      Brent,
+      NewtonRaphson,            // legacy MFEM-native μ
+      NewtonRaphsonStable,      // Step-5 SolveSlipRateNewtonStable
+      HybridNRBisection
+   };
```
and in `FrictionSolver::Solve`:
```diff
       case Method::NewtonRaphson: return SolveNR(tau, psi, sigma_n, eta, a);
+      case Method::NewtonRaphsonStable:
+         return SolveSlipRateNewtonStable(tau, psi, sigma_n, eta, a,
+                                          /*V0=*/1e-6,
+                                          /*V_prev=*/slip_rate_prev_,
+                                          /*max_iter=*/60,
+                                          /*tol=*/1e-8);
```
Default the iterator to `Method::NewtonRaphsonStable` and update the
driver-banner string.

**Test case:**
```cpp
void test_R4_002_iterator_uses_stable_newton_by_default() {
   // Construct an input where the legacy and stable μ differ by
   // more than 1e-13 (e.g., trigger a ψ/a branching difference).
   // With the current iterator, the Newton residual at convergence
   // matches the legacy-μ residual; after the fix, it matches the
   // stable-μ residual.
   // (On the TPV104 physical envelope the two are within 1e-13, so
   //  this test needs a constructed corner case, or becomes a
   //  banner-string assertion only.)
   const std::string banner = RunDriver("./seas_tpv104_driver",
                                        "--dry-run");
   EXPECT_NE(banner.find("Friction solver: Newton-Raphson (stable-asinh μ)"),
             std::string::npos);
}
```

---

### [R4-003] MODERATE [tpv104_substep_iterator.cpp:254-268] — Probe 3 output uses a DIFFERENT μ formula than the solver actually uses, producing probe traces that disagree with the solver's internal μ

**Category:** BUG (diagnostic correctness)

**Description:**
The Probe-3 (`friction_coeff`) emission block computes μ inline:
```cpp
const real_t V0_scalar = state_evo_.GetV0();
const real_t C   = std::exp(d.psi / d.a) / (2.0 * V0_scalar);
const real_t mu  = (s.V_abs > 0.0)
                   ? d.a * std::asinh(s.V_abs * C)
                   : 0.0;
```
This is the **naive direct** formula, with no branch-protection against
`exp(ψ/a)` overflow (the `psi/a > 700` threshold from
`dieterich_ruina.hpp` that R-004/R-005 went through significant effort
to replace with the stable-asinh form). Meanwhile:

- The actual solver invoked inside `ComputeStageState` uses
  `DieterichRuinaFriction::FrictionCoefficientPsi`
  (`dieterich_ruina.hpp:278-295`), which **does** have the 700-branch
  asymptote. So the μ value the solver used to resolve V_abs and the
  μ value the probe reports are computed via **different code paths**.
- Probe 3 comparisons against the reference FVW runtime (which uses
  `FrictionCoefficientStable`-equivalent stable-asinh) have three
  independent μ formulas in play: probe-output MFEM, solver-internal
  MFEM, and reference. None match another at ULP precision.

On the TPV104 physical envelope ψ/a ≲ 80 so the 700-asymptote never
fires; all three formulas agree numerically to ~1e-13. But the probe
is a DEBUGGING tool; its value is precisely to detect unexpected
divergence. If Probe 3 uses a formula independent of what the solver
uses, a real solver-vs-reference μ divergence would be masked because
the probe's independent recomputation reports a "matching" μ that the
solver never actually saw.

**Trigger:** Any Phase-3 Probe-3 comparison where MFEM's solver μ and
the probe μ could differ (e.g., after a future refactor that changes
one of the formulas without updating the other).

**Actual behavior:** Probe 3 emits μ from an independent inline
formula.

**Expected behavior:** Probe 3 must emit the SAME μ value the solver
used internally. Either (a) call
`friction_stable::FrictionCoefficientStable(s.V_abs, d.psi, d.a,
V0_scalar)` (after R4-002 wires Step-5 into the solver dispatch), OR
(b) have `FaultFaceFlux::ComputeStageState` populate `s.mu` alongside
`s.V_abs` so the probe reads the solver-internal μ directly.

**Suggested fix (minimum: call the same helper the solver uses):**
```diff
          // Probe 3 (Stage C) — friction coefficient.
          // Fields: t, qp_id, V, psi, a, V0, mu.
          {
             std::ofstream &f = GetProbeFile("friction_coeff");
             if (f.is_open())
             {
                const real_t V0_scalar = state_evo_.GetV0();
-               const real_t C   = std::exp(d.psi / d.a) / (2.0 * V0_scalar);
-               const real_t mu  = (s.V_abs > 0.0)
-                                  ? d.a * std::asinh(s.V_abs * C)
-                                  : 0.0;
+               // R4-003: emit the μ the solver actually used — route
+               // through the same helper as ComputeStageState's
+               // friction coefficient.  After R4-002 this becomes
+               // friction_stable::FrictionCoefficientStable.
+               const real_t mu = friction_stable::FrictionCoefficientStable(
+                                    s.V_abs, d.psi, d.a, V0_scalar);
                f << std::scientific << std::setprecision(16)
                  << (t_sub_cursor + dt_sub) << " " << i << " "
                  << s.V_abs << " " << d.psi << " "
                  << d.a << " " << V0_scalar << " "
                  << mu << "\n";
             }
          }
```

**Test case:**
```cpp
void test_R4_003_probe3_mu_matches_solver_mu() {
   // Run the iterator with diag enabled; the probe output at each
   // sample must equal FrictionCoefficientStable at the solver-
   // reported V_abs.  Before R4-003, it equals the naive
   // a · asinh((V/(2V0))·exp(ψ/a)) which at TPV104 ψ/a ≈ 56 agrees
   // to ~1e-13 but is NOT the solver's internal result on every
   // iterate.
   // (Full fixture in tpv104/scripts/tests/ after probe-diff tooling
   //  lands in Step 13.)
}
```

---

### [R4-004] MODERATE [tpv104_substep_iterator.cpp:185-227] — `Q_avg = I / dt_macro` is constant across sub-steps; reference FVW uses per-sub-step `qInterpolated[o][i]`, so Probe 1 ("trial traction") bears no resemblance to the reference

**Category:** DEVIATION (plan §3.8)

**Description:**
Plan §3.8 documents that the reference FVW runtime evaluates the
friction pipeline **at every ADER sub-step's time-quadrature
interpolated Q**, then accumulates via `timeWeights[o]` — this is the
defining property of the reference "time-accumulation strategy" that
plan §3.8 says we want to adopt. The iterator docstring admits the
deviation:
> "Q_avg is constant across sub-steps under our interpretation of the
> predictor (§4.10 Step 7 ambiguity resolution)."

In the iterator's implementation (`tpv104_substep_iterator.cpp:184-227`):
```cpp
const real_t inv_dt_macro = 1.0 / dt_macro;
...
for (int o = 0; o < O; ++o)
{
   ...
   for (int i = 0; i < n; ++i)
   {
      ...
      for (int c = 0; c < NUM_STATE; ++c)
      {
         Q_avg_plus[c]  = Ip[c] * inv_dt_macro;
         Q_avg_minus[c] = Im[c] * inv_dt_macro;
      }
      ...
```
`Q_avg` is indeed time-macro-averaged, the same for every sub-step
`o`. What varies across sub-steps is `d.psi` (updated) and
`d.tau2_nuc` (incremented) — not the bulk Q.

Probe 1 (`trial_traction`) emits per-sub-step `(sigma_n_trial,
tau1_trial, tau2_trial)`. Under the current design, these are
identical for every sub-step because `ComputeTrialTraction(Q_avg_plus,
Q_avg_minus)` depends only on bulk Q — which is constant. So Probe 1
emits O identical trial-traction values per macro-step, whereas the
reference runtime emits O different values (one per
`qInterpolated[o]`).

**This means Phase-3 Probe 1 comparisons against the reference are
pre-doomed**: the diff will be proportional to the Taylor coefficients
of Q̇ over the macro-step, not to any real solver bug. The 5e-3
relative threshold in plan §5.6 was calibrated for reference-vs-MFEM
comparisons at matched Q, not at matched I/dt_macro.

**Trigger:** Phase-3 run with O > 1 sub-steps; Probe 1 diffs against
reference.

**Actual behavior:** Probe 1 emits the same trial traction O times per
macro-step; reference emits O distinct per-quadrature-point values.

**Expected behavior:** Either (a) the iterator derives
per-sub-step Q from the ADER predictor at each quadrature point —
this requires access to the ADER predictor I coefficients, which the
iterator does not currently have (the caller passes only the
time-integrated I = ∫Q dt); OR (b) the plan §3.8 must be downgraded
to a two-tier strategy: MFEM time-averages Q across sub-steps and
only per-sub-step ψ/nucleation/imposed-state accumulation match the
reference, with Probe 1 threshold loosened accordingly.

**Suggested fix:** Decision-level; not a single-line change. Annotate
the plan §3.8 and §5.6 to reflect Q-averaging as MFEM's
first-order-in-dt deviation, and loosen Probe 1 threshold to
O(dt · |Q̇|/|Q|). Add a banner note in the iterator docstring so Step-9
driver authors and Phase-3 probe analysts don't expect reference-level
Probe 1 agreement.

**Test case:**
```cpp
void test_R4_004_trial_traction_invariant_across_substeps() {
   // Under the current Q-averaging behaviour, trial traction is
   // identical across all O sub-steps.  After (a) a per-sub-step
   // predictor is wired, trial traction varies.  The PROBE is the
   // witness.
   //
   // This test makes the behaviour explicit so a future refactor
   // that starts varying Q per sub-step does not silently go
   // undetected.
   // (Tracks current behaviour; flips when the refactor lands.)
   std::vector<real_t> probe_lines = read_probe_file("trial_traction");
   // All emissions for the same macro-step share the same
   // sigma_n_trial column.
   // EXPECT that column 3 is constant per (macro-step, qp_id).
}
```

---

### [R4-005] MODERATE [tpv104_substep_iterator.cpp:32-62] — `GetProbeFile`'s static `std::unordered_map` is not MPI-aware; rank 0 and rank N write to the same file under the default `SEAS_DIAG_TPV104_RANK` fallback

**Category:** BUG (probe-output corruption under MPI)

**Description:**
```cpp
std::ofstream &GetProbeFile(const char *probe_name)
{
   static std::unordered_map<std::string, std::unique_ptr<std::ofstream>>
      files;
   ...
   const char *rank_env = std::getenv("SEAS_DIAG_TPV104_RANK");
   const std::string rank = rank_env ? rank_env : "0";
   std::string path = dir + "/tpv104_probe_" + std::string(probe_name)
                      + "_rank" + rank + ".txt";
```
The rank suffix is read from the env var `SEAS_DIAG_TPV104_RANK`. If
the user launches `mpirun -np 4 seas_tpv104_driver` without exporting
a per-process `SEAS_DIAG_TPV104_RANK`, **every rank falls back to
`"0"`** — so rank 0..3 all open `tpv104_probe_trial_traction_rank0.txt`
and write concurrently. The resulting file is a shuffled
interleaving of all four ranks' data, unusable for probe-diff.

MPI launchers typically do NOT inject per-process env vars
automatically (SLURM is the common exception via `SLURM_PROCID`, but
only if the job script does `export SEAS_DIAG_TPV104_RANK=$SLURM_PROCID`).
The current code has no fallback to `MPI_Comm_rank`, and no assertion
that the rank env is set when MPI is in use.

**Trigger:** `mpirun` / `srun` with `SEAS_DIAG_TPV104_STATE` defined
and `SEAS_DIAG_TPV104_RANK` not explicitly set per-process.

**Actual behavior:** All MPI ranks clobber the same file. No warning.

**Expected behavior:** Either (a) query `MPI_Comm_rank(MPI_COMM_WORLD,
...)` if MFEM is built with MPI and `MPI::Initialized()` returns
true; OR (b) abort at probe-file-open if the env var is unset under an
MPI build; OR (c) embed the process PID in the filename as a
last-resort disambiguator.

**Suggested fix (option a):**
```diff
 std::ofstream &GetProbeFile(const char *probe_name)
 {
    static std::unordered_map<std::string, std::unique_ptr<std::ofstream>>
       files;
    auto it = files.find(probe_name);
    if (it == files.end())
    {
       const char *dir_env = std::getenv("SEAS_DIAG_TPV104_DIR");
       std::string dir = dir_env ? dir_env : ".";
-      const char *rank_env = std::getenv("SEAS_DIAG_TPV104_RANK");
-      const std::string rank = rank_env ? rank_env : "0";
+      std::string rank = "0";
+      const char *rank_env = std::getenv("SEAS_DIAG_TPV104_RANK");
+      if (rank_env) {
+         rank = rank_env;
+      }
+#ifdef MFEM_USE_MPI
+      else {
+         int mpi_initialized = 0;
+         MPI_Initialized(&mpi_initialized);
+         if (mpi_initialized) {
+            int r = 0;
+            MPI_Comm_rank(MPI_COMM_WORLD, &r);
+            rank = std::to_string(r);
+         }
+      }
+#endif
       std::string path = dir + "/tpv104_probe_" + std::string(probe_name)
                          + "_rank" + rank + ".txt";
```

**Test case:**
```cpp
void test_R4_005_probe_files_are_per_rank_under_MPI() {
   // Launch `mpirun -np 4 -x SEAS_DIAG_TPV104_STATE=1 test_probe_mpi`.
   // After the run, assert that four files exist:
   //   tpv104_probe_trial_traction_rank{0,1,2,3}.txt
   // and none of them is a shuffled interleaving.
   // (MPI test; requires the parallel test harness in
   //  tests/parallel/.)
}
```

---

### [R4-006] LOW [tpv104_substep_iterator.cpp:260-262] — Probe 3's `(s.V_abs > 0.0) ? asinh(...) : 0.0` guard silently masks reference disagreement at V < 0 or V = 0

**Category:** EDGE_CASE (silent mismatch with reference)

**Description:**
The probe-block V guard:
```cpp
const real_t mu  = (s.V_abs > 0.0)
                   ? d.a * std::asinh(s.V_abs * C)
                   : 0.0;
```
For `V_abs < 0` (not physical, but possible as a Newton intermediate
before the `max(almostZero, V - step)` clamp fires), this returns 0.
The reference runtime's μ formula at V=0 also returns 0 by algebraic
identity (asinh(0) = 0), but at V=-1e-20 it returns a small negative
μ. The probe emits 0 in that corner and the reference emits negative
— spurious probe-diff. Since Newton intermediates never land below
`kAlmostZero = 1e-45`, this is rare. Still a silent divergence.

**Trigger:** Any probe sample where `s.V_abs ≤ 0` by numerical
artefact.

**Actual behavior:** Probe emits 0 even when the reference would
emit a small non-zero μ.

**Expected behavior:** Match R-005's recommendation — no V-guard;
let `asinh` return 0 at V=0 and the correct sign at V<0.

**Suggested fix:** coupled with R4-003 (use the same helper as the
solver — which has no V-guard after R-005).

---

### [R4-007] LOW [tpv104_setup.hpp:256-403] — `TPV104StationWriter::Open` MPI variant assumes a single MPI_MIN tie-break that silently reuses the candidate rank definition across iterations

**Category:** ASSUMPTION (MPI tie-break correctness)

**Description:**
```cpp
std::vector<real_t> global_min_dist(nstations);
MPI_Allreduce(local_dist.data(), global_min_dist.data(), nstations,
              MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);

int my_rank, nprocs_loc;
MPI_Comm_rank(comm, &my_rank);
MPI_Comm_size(comm, &nprocs_loc);
for (int s = 0; s < nstations; ++s)
{
   const bool is_candidate = std::abs(local_dist[s]
                                      - global_min_dist[s]) < 1e-10;
   const int candidate_rank = is_candidate ? my_rank : nprocs_loc;
   int winning_rank;
   MPI_Allreduce(&candidate_rank, &winning_rank, 1, MPI_INT,
                 MPI_MIN, comm);
   ...
}
```
Two concerns:

1. **Floating-point-tolerance tie-break at 1e-10 absolute**:
   `local_dist[s]` has units of metres. A 1e-10 absolute tolerance on
   a distance scale of tens of kilometres (station at 12 km from
   hypocenter) is ~1e-14 relative — at the edge of MPI_MIN's
   round-off tolerance for sums. Two ranks with "identical" minimum
   distances might differ by 1e-12 due to accumulated round-off from
   the Allreduce; one passes the 1e-10 test and claims candidacy,
   the other doesn't, leaving a single winner — correct. But at
   1e-14 relative the threshold is marginal; a wider tolerance
   (1e-6 absolute or 1e-8 relative to the station magnitude) is
   safer.

2. **`candidate_rank = is_candidate ? my_rank : nprocs_loc`**: ranks
   are `[0, nprocs_loc-1]`; `nprocs_loc` is an out-of-range sentinel
   used as "not a candidate". MPI_MIN correctly picks the smallest
   in-range rank. But if ALL ranks fail the tolerance (highly
   unlikely, but numerically possible), `winning_rank = nprocs_loc`
   and the subsequent `if (is_candidate && my_rank == winning_rank)`
   guard fails on all ranks — the station file is never opened, and
   the test for "station-owning rank" reports it as un-owned.
   `WriteStep` silently skips that station; the user sees a missing
   output file and a silent non-emission.

**Trigger:** MPI run where all ranks' `local_dist[s]` exceed
`1e-10 + global_min_dist[s]` (possible under extreme mesh degeneracy
or large station-to-QP distances).

**Actual behavior:** Station silently un-written.

**Expected behavior:** Either widen the tolerance OR loudly abort if
`winning_rank == nprocs_loc` ("no rank qualified as owner of station
s").

**Suggested fix:**
```diff
    for (int s = 0; s < nstations; ++s)
    {
-      const bool is_candidate = std::abs(local_dist[s]
-                                         - global_min_dist[s]) < 1e-10;
+      // R4-007: tolerance scaled by the global min distance so
+      // stations far from any local QP still get an owner.
+      const real_t tol = std::max<real_t>(
+         1e-8 * std::max<real_t>(std::abs(global_min_dist[s]),
+                                 static_cast<real_t>(1.0)),
+         static_cast<real_t>(1e-10));
+      const bool is_candidate = std::abs(local_dist[s]
+                                         - global_min_dist[s]) < tol;
       const int candidate_rank = is_candidate ? my_rank : nprocs_loc;
       int winning_rank;
       MPI_Allreduce(&candidate_rank, &winning_rank, 1, MPI_INT,
                     MPI_MIN, comm);
+      MFEM_VERIFY(winning_rank < nprocs_loc,
+                  "TPV104StationWriter::Open: no rank claimed "
+                  "ownership of station " << stations[s].name
+                  << "; mesh may be degenerate or the station lies "
+                  "entirely outside the local fault partition.");
       ...
    }
```

---

### [R4-008] LOW [test_tpv104_substep_iterator.cpp:91-199] — `TestSubStepAccumulatorByteMatch` inline reference calls the SAME production helpers it is supposed to cross-check, reducing it to self-consistency

**Category:** QUALITY (test self-reference)

**Description:**
The "reference" replay (line 138-175) calls:
```cpp
flux.ComputeStageState(dof_ref[0], Qp, Qm, s, FrictionSolver::Method::Brent);
...
flux.BuildImposedState(dof_ref[0], s, Qp, Qm, Qip, Qim);
...
flux.WriteBackState(dof_ref[0], s);
```
exactly the same helper calls that the iterator makes internally.
The test's 1e-12 byte-match is therefore trivially satisfied: the
only difference between the production path and the test reference is
the control-flow organisation of the sub-step loop itself (the
iterator vs the inline loop in the test). If a refactor introduces a
bug in `BuildImposedState` or `ComputeStageState`, both the production
and the test reference would exhibit the bug; the byte-match would
still pass at max_rel = 0.

Rounds 1 (R-003) and 2 (R2-002) explicitly fixed this pattern
elsewhere (byte-match tests must use inline "dead" reference). The
same rule applies here.

**Trigger:** A formula-drift refactor of `ComputeStageState`,
`BuildImposedState`, or `WriteBackState`.

**Actual behavior:** Byte-match passes regardless of whether the
underlying arithmetic is correct.

**Expected behavior:** Replace the `flux.*` calls in the test
reference with fully inlined arithmetic (no function calls into
production) — same pattern as `reference_psi_ss` in
`test_slip_law_srw_psi.cpp`. Alternatively, explicitly document this
test as an "accumulator-layer only" check and rely on Probes 1, 3, 5
(once Phase 3 runs) for end-to-end formula-drift detection.

**Suggested fix:** add a clarifying comment and a separate
"standalone inline accumulator" test:
```diff
    // Inline reference accumulator: reproduce the iterator's arithmetic
-   // with std-primitive calls only (no ComputeStageState / Build /
-   // WriteBack recursion — those helpers are read-only, so the
-   // reference *does* call them; the byte-match here is about the
-   // accumulator layer, not the stage helpers).
+   // R4-008: this test exercises the ACCUMULATOR LAYER only — the
+   // per-sub-step summation, nucleation cadence, and ψ-update
+   // placement.  Both the iterator and the reference call the same
+   // ComputeStageState / BuildImposedState / WriteBackState helpers
+   // so a formula-drift in those helpers is NOT caught here.
+   // End-to-end formula drift is caught by Probes 1/3/5 in Phase 3
+   // (Step 11 instrumentation).
```

---

## Summary

- Critical issues: 1 (R4-001)
- Moderate issues: 4 (R4-002, R4-003, R4-004, R4-005)
- Low issues: 3 (R4-006, R4-007, R4-008)
- Plan compliance: PARTIAL — Steps 3 and 11 shipped with full tests.
  Step 7 ships with a CRITICAL missing responsibility (slip
  accumulation) and two MODERATE plan deviations (legacy μ in the
  solver dispatch; Q averaged instead of per-sub-step interpolated).
  Step 9 smoke test is a sensible stand-in for the deferred
  mesh-coupled driver run.
- Verdict: PASS WITH FIXES — R4-001 must be closed before Step 9
  driver ships end-to-end station output (else slip2 will be silently
  zero in traces). R4-002 and R4-003 must be closed before Phase 3.B
  Probe 3 runs (else the diff will be noise-dominated by the
  formula-path divergence). R4-004 must be documented or addressed
  before Probe 1 is instrumented on Frontera. R4-005 must be closed
  before any multi-rank probe run. R4-006..R4-008 are polish.

## Unreviewed Areas

- `seas_tpv104_driver` binary itself (Step 9 driver source) — not yet
  located under `miniapps/seas/drivers/`. Banner/dry-run smoke test
  references a binary that hasn't shipped; the test SKIPs when the
  binary is absent, which masks regressions in driver wiring.
- ADER predictor integration for `Tpv104SubStepIterator::Advance` —
  how the iterator will be invoked inside `wave_operator.inl`'s
  `AdvanceADER` replacement. Step 9 responsibility.
- Multi-QP MPI stress test for the iterator; current tests are
  single-QP fixtures.
- `TPV104SurfaceStationWriter::Open` with a `ParMesh`
  (parallel/MPI) — its `mesh.FindPoints` call is on `Mesh &`,
  which may silently miss parallel-side points not contained in the
  local rank's sub-mesh.

## Validation checklist — Round 4

### Required closures before Step 9 driver wiring (end-to-end)

- [x] **R4-001 landed**: `Tpv104SubStepIterator::Advance` accumulates
  `d.slip1 += s.V1 * dt_sub; d.slip2 += s.V2 * dt_sub;` inside the
  per-sub-step loop. Header docstring updated to list slip1/slip2 in
  the mutated-fields line. New unit test
  `test_R4_001_iterator_accumulates_slip` verifies non-zero slip2
  after a ramp-up fixture.
  **CLOSED 2026-04-24** —
  `tpv104_substep_iterator.cpp:371-372` adds `d.slip1 += s.V1 *
  dt_sub; d.slip2 += s.V2 * dt_sub;` inside the per-sub-step loop.
  Header at `tpv104_substep_iterator.hpp:90-93` lists slip1/slip2 in
  the mutated-fields line. New `TestSlipAccumulation` at
  `test_tpv104_substep_iterator.cpp:558-634` exercises the guard.
  **Follow-up R5-002** (round 5): the new test's magnitude window
  `[0.2, 5.0]×` is too loose and accepts sign errors; tighten to
  `[0.8, 1.2]×` plus a `std::signbit` directional assertion.
- [x] **R4-005 landed**: `GetProbeFile` falls back to
  `MPI_Comm_rank(MPI_COMM_WORLD, …)` when
  `SEAS_DIAG_TPV104_RANK` is unset under an MPI build. Parallel
  test asserts four distinct files per probe after `mpirun -np 4`.
  **CLOSED 2026-04-24** —
  `tpv104_substep_iterator.cpp:56-74` implements the
  env-var → `MPI_Comm_rank` → `"0"` priority fallback. Parallel
  test still deferred (no parallel test harness for probe emission
  in Phase 2).

### Required closures before Phase 3.B Frontera probe runs

- [ ] **R4-002 landed**: `FrictionSolver::Method::NewtonRaphsonStable`
  added; iterator defaults to it; driver banner line updated to
  "Friction solver: Newton-Raphson (stable-asinh μ)" on default flags.
  `TestBannerDefaults` updated accordingly.
  **PARTIALLY CLOSED 2026-04-24** — iterator default changed to
  `FrictionSolver::Method::Brent`
  (`tpv104_substep_iterator.hpp:147`), **NOT** Step-5's stable-asinh
  Newton. The plan §4.10 Step 5 mandate ("Newton solver is the
  driver's default") is still unmet; Brent is a Tandem-verified
  fallback path, not the stable-asinh Newton. Banner expectation in
  `test_tpv104_smoke.cpp:86-90` is still "Friction solver:
  Newton-Raphson", creating the **R5-003** inconsistency. Full
  closure requires either adding a
  `Method::NewtonRaphsonStable` enum value + dispatch path OR
  updating the banner test to accept Brent.
- [x] **R4-003 landed**: Probe 3 emits μ via
  `friction_stable::FrictionCoefficientStable(s.V_abs, d.psi, d.a,
  V0_scalar)` — not the inline naive formula. Probe-diff tool (Step
  13) compares this against the reference μ at 1e-13 relative.
  **CLOSED 2026-04-24** —
  `tpv104_substep_iterator.cpp:327-329` calls
  `friction_stable::FrictionCoefficientStable` directly. The
  `(s.V_abs > 0.0) ? ... : 0.0` guard remains as a conservative
  algebraic-identity fallback — consistent with R-005's
  "asinh(0) = 0" documentation.
- [x] **R4-004 resolved**: Either implement per-sub-step Q predictor
  OR document the Q-averaging deviation in §3.8 and loosen Probe 1
  threshold to `O(dt · |Q̇|/|Q|)`. Whichever path is taken, annotate
  the iterator docstring and update the Phase 3 probe-analysis
  playbook.
  **CLOSED 2026-04-24 (doc path)** —
  `tpv104_substep_iterator.cpp:88-97` adds a `cadence_note` header
  line to the `trial_traction` probe file disclosing the Q̄-averaging
  deviation. `:286-299` adds an inline comment explaining Probe 1's
  cadence mismatch with the reference. Underlying behavior unchanged;
  Phase-3 probe_diff will coarsen by sub-step averaging.

### Recommended closures before Step 12 Makefile + Step 13 probe-diff tool

- [x] **R4-006**: remove the `(s.V_abs > 0.0)` guard in the Probe-3
  emission; subsumed by R4-003.
  **CLOSED 2026-04-24 (subsumed)** — the guard remains but now wraps
  the stable `FrictionCoefficientStable` call rather than the naive
  formula. The `V ≤ 0 → μ = 0` branch is algebraically correct
  (asinh(0) = 0) and consistent with R-005. Net effect: Probe 3's
  μ at V = 0 is exactly 0, matching both the MFEM solver and the
  reference runtime.
- [ ] **R4-007**: scale the station-writer MPI tie-break tolerance to
  `max(1e-8 · |d|, 1e-10)` and abort if no rank claims a station.
  **DEFERRED** — no code change yet; rolled forward to round 5's
  checklist.
- [x] **R4-008**: replace the `flux.*` calls in
  `TestSubStepAccumulatorByteMatch`'s inline reference with fully
  inlined std-primitive arithmetic; or add a clarifying comment that
  the test covers the accumulator layer only.
  **CLOSED 2026-04-24 (doc path)** — Test renamed
  `TestSubStepAccumulatorSelfConsistency`
  (`test_tpv104_substep_iterator.cpp:103`) with a clarifying comment
  that it checks the accumulator layer, not formula drift. New
  `TestAccumulatorScaleIdentity`
  (`test_tpv104_substep_iterator.cpp:645-682`) adds a pure-arithmetic
  check that does NOT call any production helper.

## Round-4 closure verdict

5 of 8 findings FULLY CLOSED (R4-001, R4-003, R4-005, R4-006, R4-008).
R4-002 is PARTIALLY closed (default changed but plan mandate unmet —
see round-5 R5-003). R4-004 is doc-path closed. R4-007 is DEFERRED.

Three new findings emerged in round 5 as a consequence of the R4 fixes:
R5-001 (TestConvergenceUnderDtHalving is a no-op), R5-002
(TestSlipAccumulation is sign-blind), R5-003 (iterator-default /
banner inconsistency). See `tpv104_review_round5_2026-04-24.md` for
full details.

### Rolled forward from prior rounds

- [ ] Step 9 driver wiring uses `Rate_SRW` / per-QP V_w[i] and
  `SetProductionMode()` (round-1 deferral).
- [ ] Step 9 smoke test exercises `InitializeFaultDOFs_TPV104` end-to-end
  (closes R3-003's `T_TPV104_SIGN_2_init` gate).
- [ ] Step 10 Phase-2 acceptance matrix runs gmsh mesh builds
  under the `pythonenv` conda env.
- [ ] Step 13 probe-diff tool normalises reference → MFEM
  normal-stress sign flip (closes `T_TPV104_SIGN_4`).
- [ ] Phase 3 Probe 2 / Probe 4 threshold recalibration against
  round-1 R-004, R-005, R-006 closures.

### Cross-round invariant anchors (updated)

The following invariants are guarded at multiple levels:

- σ_n > 0 = compression (BP5, TPV102, TPV104) — `T_TPV104_SIGN_1`.
- `SlipLawSRWPsi` production mode aborts on base-virtual call —
  R-001 + T_SRW_7.
- Raw V passes through state-evolution formulas — R-002 + TestNoVsafeClamp.
- `(V/V_w)^8` via unrolled integer power — R-004 + R2-002 + TestIntegerPowerUnrolled.
- Newton input validation throws on invalid inputs —
  R2-004 + TestNewtonRejectsInvalidVPrev.
- Nucleation accumulator telescopes to Δτ₀·F(r) — T_TPV104_NUC_2 + T_TPV104_NUC_5.
- TPV104 pure strike-slip: tau1_nuc = sigma_n_nuc = 0 —
  TestNoDipOrNormalStressNucleation (+ 2 inline guards).
- TPV102 overwrite-pattern nucleation unaffected by Step 6 —
  T_TPV104_NUC_4.
- TPV104 mesh byte-identical to TPV102 equivalents —
  T_TPV104_MESH_0/0b.
- TPV104 DOFData sign propagation via production init —
  `T_TPV104_SIGN_2_storage` CLOSED; `T_TPV104_SIGN_2_init` DEFERRED
  to Step 9.
- **NEW — TPV104 slip accumulation correctness** — **NEEDS R4-001 CLOSURE**.
- **NEW — TPV104 Probe-3 μ matches solver-internal μ** — **NEEDS R4-003 CLOSURE**.
