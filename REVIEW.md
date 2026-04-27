# Code Review: TPV102 setup & workflow vs TPV104 (proven working) — 2026-04-27

## Review Scope
- Plan: `/Users/chunhuizhao/projects/seas-mfem/PLAN.md` (BP5 ParaView; not the active TPV102 plan — used only to confirm shared paraview API).
- Implementation report (commit log): `git log --oneline -- miniapps/seas/drivers/tpv102_driver.cpp` shows `ffffc35 TPV102 pepper + BP5 in-progress scope: debug dumps, build support` (HEAD), `28284b9 TPV102 nucleation: persistent-prestress channel`, `632ad03 TPV102 v9.3.0: ADER default + total-Q only + ADER+total-Q 400-rank sbatch pair`, with the current driver a fresh port of TPV104 per `Makefile:352–369` ("TPV102 driver rebuilt against the TPV104 code flow").
- Files reviewed (TPV102 and the matched TPV104 reference):
  - `miniapps/seas/drivers/tpv102_driver.cpp` (modified, 2181 lines) vs `drivers/tpv104_driver.cpp` (2231 lines)
  - `miniapps/seas/dynamic/tpv102_setup.hpp` (modified, 504 lines) vs `dynamic/tpv104_setup.hpp` (552 lines)
  - `miniapps/seas/dynamic/tpv102_friction_solver.hpp` (untracked, 48 lines) vs `dynamic/tpv104_friction_solver.hpp` (178 lines)
  - `miniapps/seas/dynamic/tpv102_nucleation.hpp` (untracked, 133 lines) vs `dynamic/tpv104_nucleation.hpp` (143 lines)
  - `miniapps/seas/dynamic/tpv102_substep_iterator.{hpp,cpp}` (untracked, 145 + 378 lines) vs `dynamic/tpv104_substep_iterator.{hpp,cpp}` (224 + 687 lines)
  - `miniapps/seas/config/tpv102_params.hpp` vs `config/tpv104_params.hpp`
  - `miniapps/seas/friction/state_evolution.hpp` (modified — adds `AgingLawPsi::GetB/GetV0/GetF0`)
  - `miniapps/seas/Makefile` (TPV102_HEADERS / TPV102_SHARED_OBJS / `seas_tpv102_driver` rule)
  - One representative TPV102 sbatch (`jobs/tpv102/tpv102_200m_p1_1.5s_400rank_dev.sbatch`) and the TPV102 setup unit test (`tests/unit/test_tpv102_setup.cpp`).
- Domain context consulted: `CLAUDE.md` (project), `miniapps/seas/CLAUDE.md` (SEAS — *Friction Solver*, *Sign Conventions*), `tpv102/TPV102_GUIDE.md`, `feedback_tpv102_bp5_no_shared_edit`, `feedback_dynamic_folder_editable_for_tpv104`.

The TPV102 path is structurally a faithful port of TPV104 — same call sequence, same WaveOperator API, same FaultFaceFlux dispatch, same paraview wiring. I did not find a smoking-gun runtime-fatal divergence on the default (one-shot, Brent, `--paraview-dt`) configuration that the active sbatch scripts exercise; the default path mirrors TPV104. The findings below are the substantive deviations between the two paths, prioritised by what would block a TPV102 substep / probe / parameter-sweep run.

## Findings

### [R-001] [MODERATE] [POSSIBLE] [drivers/tpv102_driver.cpp:521,780,1971–1976; dynamic/tpv102_substep_iterator.hpp:111–127] — `--fault-iterator substep` defaults to NewtonRaphsonStable, contradicting CLAUDE.md "Brent required" for TPV102

**Category:** BUG / ASSUMPTION

**Description:**
The TPV102 substep iterator path forwards `MapSolver(friction_solver)` (default `NewtonRaphsonStable`, see `tpv102_driver.cpp:521,780`) into `Tpv102SubStepIterator::AdvanceWithSubStepStates(... method)`, which in turn calls `FaultFaceFlux::ComputeStageState(d, Q+, Q-, s, NewtonRaphsonStable)` (`tpv102_substep_iterator.cpp:343`). Project-level `miniapps/seas/CLAUDE.md` (§ "Friction Solver") records that TPV102/BP5 require Brent because Newton fails when ψ/a is large — and TPV102's equilibrium ψ/a is large by spec.

For `a_vw = 0.008`, `tau_ini = 75 MPa`, `sigma_n = 120 MPa`, `V_ini = 1e-12`, `V0 = 1e-6`:
`ComputeInitialPsi(a)` at `tpv102_params.hpp:131–139` returns `psi ≈ 0.74` → `psi/a ≈ 92`. The legacy `Newton-Raphson` documented in `debug v1` failed at `psi/a` of similar magnitude. `NewtonRaphsonStable` uses `friction_stable::FrictionCoefficientStable` (asinh-form, better than legacy) but still applies an unbracketed Newton iterate seeded with `V_prev = 1e-12` (the very `V_lo` regime CLAUDE.md flags as the failure mode).

By contrast, TPV104 has `a_in = 0.01`, ψ_ini ≈ 0.564, `psi/a ≈ 56` — well within the stable-asinh's empirical regime. So the same iterator default that's safe for TPV104 is risky for TPV102.

Production runs default to `--fault-iterator one-shot` (Brent, hard-coded inside `EvaluateADER`), so this bug is dormant by default. It activates the moment a user passes `--fault-iterator substep` (or sets it via a sbatch).

**Trigger:**
Run with `--fault-iterator substep` on a TPV102 mesh at `t=0` (V_ini=1e-12, ψ ≈ 0.74). The Newton solve at the locked initial state with `V_prev = 1e-12` either fails to converge in 60 iterations or returns an unphysical iterate.

**Actual behavior:**
`SolveSlipRateNewtonStable` exits with `*has_converged = false` (silent — the convergence flag is captured but not asserted by the iterator). The `EvalStageState` carries this unconverged V and propagates into the imposed-state accumulator, eventually yielding NaN or non-physical V_max in the time loop.

**Expected behavior:**
The substep path on TPV102 should default to `Method::Brent` — the same dispatch the one-shot path uses today. Brent is bracketed and works at large ψ/a per `debug v1`.

**Suggested fix:**
Override the default at the TPV102 substep iterator declaration (smallest blast radius):
```diff
--- a/miniapps/seas/dynamic/tpv102_substep_iterator.hpp
+++ b/miniapps/seas/dynamic/tpv102_substep_iterator.hpp
@@
    void Advance(std::vector<DOFData> &dof_data,
                 const std::vector<Vector> &fault_coords,
                 const real_t *I_plus_flat,
                 const real_t *I_minus_flat,
                 real_t dt_macro,
                 real_t t_macro_start,
                 real_t *I_imp_plus_flat,
                 real_t *I_imp_minus_flat,
                 FrictionSolver::Method method
-                  = FrictionSolver::Method::NewtonRaphsonStable);
+                  = FrictionSolver::Method::Brent);
@@
    void AdvanceWithSubStepStates(
       std::vector<DOFData> &dof_data,
       const std::vector<Vector> &fault_coords,
       const std::vector<std::vector<real_t>> &Q_pointwise_plus_per_substep,
       const std::vector<std::vector<real_t>> &Q_pointwise_minus_per_substep,
       real_t dt_macro,
       real_t t_macro_start,
       real_t *I_imp_plus_flat,
       real_t *I_imp_minus_flat,
       FrictionSolver::Method method
-         = FrictionSolver::Method::NewtonRaphsonStable);
+         = FrictionSolver::Method::Brent);
```
Defaults live only in the header; no `.cpp` change is needed beyond a matching comment update at `tpv102_substep_iterator.cpp:84,224` documenting the new default.

This keeps the user's ability to override with `--friction-solver newton-stable` for cross-code diagnostics and remains symmetric with `tpv104_substep_iterator.hpp` (TPV104 retains its NewtonRaphsonStable default because its smaller ψ/a is empirically safe).

**Test case:**
```cpp
// tests/unit/test_tpv102_substep_friction_solver_default.cpp
TEST(R001_tpv102_substep_default_brent) {
   // Synthetic single-QP DOFData at TPV102 equilibrium (ψ/a ≈ 92).
   DOFData d{};
   d.Zp_plus = d.Zp_minus = TPV102Params::Zp;
   d.Zs_plus = d.Zs_minus = TPV102Params::Zs;
   d.eta_p   = TPV102Params::Zp / 2.0;
   d.eta_s   = TPV102Params::eta_s;
   d.sigma_n0 = TPV102Params::sigma_n;
   d.tau2_0   = TPV102Params::tau_ini;
   d.a   = TPV102Params::a_vw;     // 0.008
   d.Dc  = TPV102Params::Dc;
   d.psi = ComputeInitialPsi(d.a); // ≈ 0.74 → ψ/a ≈ 92
   d.slip_rate = TPV102Params::V_ini;          // 1e-12
   d.V2 = TPV102Params::V_ini;

   FaultFaceFlux flux(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   AgingLawPsi state_evo(TPV102Params::b, TPV102Params::V0, TPV102Params::f0);
   Tpv102SubStepIterator iter(flux, state_evo);
   iter.SetSubSteps({1e-3}, {1.0});  // single-sub-step, dt=1ms

   std::vector<DOFData> dof_data{d};
   std::vector<Vector> coords(1, Vector(3));
   coords[0] = 0.0; coords[0](2) = -7.5e3;

   std::vector<std::vector<real_t>> Qp(1, std::vector<real_t>(NUM_STATE, 0.0));
   std::vector<std::vector<real_t>> Qm(1, std::vector<real_t>(NUM_STATE, 0.0));
   std::vector<real_t> Iimp_p(NUM_STATE, 0.0), Iimp_m(NUM_STATE, 0.0);

   // Default method must be Brent and must converge.
   iter.AdvanceWithSubStepStates(dof_data, coords, Qp, Qm,
                                 1e-3, 0.0, Iimp_p.data(), Iimp_m.data());
   EXPECT_TRUE(std::isfinite(dof_data[0].slip_rate));
   EXPECT_LT(dof_data[0].slip_rate, 1e-6);  // still locked / quasi-locked
}
```
Run as ASSERT under both default and explicit `Method::Brent` to confirm bracketed convergence; the legacy default (NewtonRaphsonStable) should produce NaN or absurdly large `slip_rate` at this equilibrium.

---

### [R-002] [MODERATE] [drivers/tpv102_driver.cpp:435–449; dynamic/tpv102_substep_iterator.cpp (entire file)] — Driver advertises `SEAS_DIAG_TPV102_STATE = ON/OFF` but the iterator emits no probe trace

**Category:** DEVIATION (banner / build-info contract not implemented)

**Description:**
`tpv102_driver.cpp:435–449` writes `[BUILD] SEAS_DIAG_TPV102_STATE = ON|OFF` to stderr and `build_info.txt`, mirroring `tpv104_driver.cpp:434–448`. The TPV104 sibling actually honors the macro: `tpv104_substep_iterator.cpp:14–23,30–116,322–390,420–468` opens `tpv104_probe_<name>_rank<R>.txt` files and writes per-sub-step `trial_traction`, `friction_coeff`, `state_evolution`, `slip_rate`, and `corrected_imposed` traces, plus a static `Tpv104SubStepIterator::CloseAllProbeFiles()` invoked from the driver before `MPI_Finalize()`.

The TPV102 substep iterator (`tpv102_substep_iterator.cpp` start-to-end) has **no** `#ifdef SEAS_DIAG_TPV102_STATE` blocks at all — not the file-registry header, not any probe writes, not a `CloseAllProbeFiles` static. Building with `-DSEAS_DIAG_TPV102_STATE` therefore advertises probe diagnostics that never get emitted; Phase-3 cross-code probe-diff tooling (the same machinery TPV104 uses for SeisSol comparison per `tpv104_substep_iterator.cpp:104–108`) cannot run on TPV102.

**Trigger:**
`make seas_tpv102_driver SEAS_EXTRA_CPPFLAGS="-DSEAS_DIAG_TPV102_STATE"` and run with `--fault-iterator substep`. The banner says probe is ON; no `tpv102_probe_*_rank*.txt` files appear in the output dir.

**Actual behavior:**
Banner contract violated; users (and CI parsers) infer that probe data is being captured when in fact it is not. There is also no symmetric `Tpv102SubStepIterator::CloseAllProbeFiles()` for the driver to call before `MPI_Finalize()` — silently fine today (no probe files to close), but inconsistent with TPV104.

**Expected behavior:**
Either (i) implement the TPV102 probe-writer block as a near-mirror of `tpv104_substep_iterator.cpp:14–23,30–116,322–390,420–468` plus a static `Tpv102SubStepIterator::CloseAllProbeFiles()`, OR (ii) drop the build banner so it does not falsely claim a feature.

**Suggested fix (low-effort path — drop the banner, keep symmetry honest):**
```diff
--- a/miniapps/seas/drivers/tpv102_driver.cpp
+++ b/miniapps/seas/drivers/tpv102_driver.cpp
@@
    if (rank == 0)
    {
 #ifdef SEAS_DIAG_FAULT_FLUX
       const char *diag_fault_flux = "SEAS_DIAG_FAULT_FLUX = ON";
 #else
       const char *diag_fault_flux = "SEAS_DIAG_FAULT_FLUX = OFF";
 #endif
-#ifdef SEAS_DIAG_TPV102_STATE
-      const char *diag_tpv102 = "SEAS_DIAG_TPV102_STATE = ON";
-#else
-      const char *diag_tpv102 = "SEAS_DIAG_TPV102_STATE = OFF";
-#endif
       std::fprintf(stderr, "[BUILD] %s\n", diag_fault_flux);
-      std::fprintf(stderr, "[BUILD] %s\n", diag_tpv102);
       std::ofstream binfo("build_info.txt");
       if (binfo.is_open())
       {
          binfo << "[BUILD] " << diag_fault_flux << "\n";
-         binfo << "[BUILD] " << diag_tpv102 << "\n";
          binfo.close();
       }
    }
```

If the longer fix is wanted (full probe parity), it must add `static void Tpv102SubStepIterator::CloseAllProbeFiles();` plus the per-iteration `GetProbeFile(...)` writes mirroring `tpv104_substep_iterator.cpp:36–116` and the five Probe-1..Probe-5 emit blocks. The driver's normal exit and NaN-tripwire branches must then call `Tpv102SubStepIterator::CloseAllProbeFiles();` before `MPI_Finalize()` (compare `tpv104_driver.cpp:1914,2156,2177`).

**Test case:**
```cpp
// tests/unit/test_tpv102_diag_state_banner.cpp
TEST(R002_diag_state_banner_implies_probe_emit) {
   // Build with -DSEAS_DIAG_TPV102_STATE.  Run --fault-iterator substep.
   // Assert: at least one tpv102_probe_<name>_rank0.txt file exists in
   //         the output dir AND the banner reports
   //         "SEAS_DIAG_TPV102_STATE = ON".
   //
   // If the banner is dropped (low-effort fix), assert the symmetric
   // negative: no probe files AND no [BUILD] SEAS_DIAG_TPV102_STATE line.
}
```

---

### [R-003] [LOW] [config/tpv102_params.hpp:131–139] — `ComputeInitialPsi` uses naive `sinh(arg)` and is at the edge of double precision

**Category:** ASSUMPTION / numerical robustness

**Description:**
`tpv102_params.hpp:131–139`:
```cpp
inline real_t ComputeInitialPsi(real_t a)
{
   real_t arg = TPV102Params::tau_ini / (TPV102Params::sigma_n * a);
   real_t psi = a * std::log(2.0 * TPV102Params::V0 / TPV102Params::V_ini * std::sinh(arg));
   return psi;
}
```

For the production VW value `a_vw = 0.008` we have `arg = 75e6/(120e6·0.008) = 78.125`, so `std::sinh(78.125) ≈ 7.0e33`. That fits in `double` (max ≈ 1.8e308), so today this works.

The TPV104 sibling `ComputeInitialPsiTPV104` (`tpv104_params.hpp:210–219`) instead uses the numerically stable `log(x · sinh(c)) = |c| + log((x/2)·-sign(c)·expm1(-2|c|))` reformulation, which is finite for arbitrarily large `|c|`. Any future tightening of TPV102 (lower `a_vw`, higher `tau_ini`, lower `V_ini`) — e.g. a TPV102-style probe with `a = 0.003` would push `arg ≈ 208`, `sinh(208) ≈ 8e89` (still in range); `a ≈ 0.0024` pushes `arg ≈ 260`, `sinh(260) ≈ 4e112`; `a ≈ 0.0017` pushes `arg ≈ 367`, `sinh` overflows to `+inf`. The naive form is a few parameter perturbations away from `inf` propagating into `psi`.

The setup unit test (`tests/unit/test_tpv102_setup.cpp:101–121` — `TestInitialStateEquilibrium`) only checks the round-trip equilibrium for `a_vw = 0.008`, which masks this fragility.

**Trigger:**
A future TPV102-variant scenario reduces `a_vw` toward ≈ 0.0017 (arg ≥ 367 → sinh=+inf), or doubles `tau_ini`.

**Actual behavior:**
`std::sinh(arg)` returns `+inf`, `log(... · inf)` returns `+inf`, `psi = a · inf = inf`. Downstream `SolveSlipRateNewtonStable` (`tpv104_friction_solver.hpp:115–120`) throws `std::runtime_error` on non-finite `psi` — loud failure rather than silent NaN, which is an acceptable failure mode but unnecessary if the formula is stable.

**Expected behavior:**
Use the same logsinh-stable formulation as TPV104.

**Suggested fix:**
```diff
--- a/miniapps/seas/config/tpv102_params.hpp
+++ b/miniapps/seas/config/tpv102_params.hpp
@@
 inline real_t ComputeInitialPsi(real_t a)
 {
-   // From tau = sigma_n * a * asinh(V/(2*V0) * exp(psi/a)):
-   //   sinh(tau/(sigma_n*a)) = V/(2*V0) * exp(psi/a)
-   //   psi = a * ln(2*V0/V * sinh(tau/(sigma_n*a)))
-   real_t arg = TPV102Params::tau_ini / (TPV102Params::sigma_n * a);
-   real_t psi = a * std::log(2.0 * TPV102Params::V0 / TPV102Params::V_ini * std::sinh(arg));
-   return psi;
+   // Numerically stable logsinh form (matches ComputeInitialPsiTPV104 in
+   // tpv104_params.hpp:210–219).  Stable for arbitrarily large arg.
+   const real_t arg = TPV102Params::tau_ini / (TPV102Params::sigma_n * a);
+   const real_t x   = 2.0 * TPV102Params::V0 / TPV102Params::V_ini;
+   const real_t sign_c = (arg >= 0.0) ? 1.0 : -1.0;
+   const real_t absC   = std::abs(arg);
+   return a * (absC + std::log(x / 2.0 * -sign_c * std::expm1(-2.0 * absC)));
 }
```

**Test case:**
```cpp
// tests/unit/test_tpv102_initial_psi_logsinh.cpp
TEST(R003_compute_initial_psi_no_overflow) {
   // Synthetic — pretend tau_ini/sigma_n is normal but a is small,
   // pushing arg into the regime where naive sinh would overflow.
   // (Drive via a parameter sweep, since TPV102Params is constexpr.)
   const real_t arg_big = 367.0;  // sinh(367) overflows double
   const real_t a_eq    = 1.0;    // no-op a, just exercises arg

   const real_t x = 2.0 * 1e-6 / 1e-12;
   const real_t sign_c = 1.0;
   const real_t absC   = arg_big;
   const real_t psi = a_eq * (absC + std::log(x / 2.0 * -sign_c
                                              * std::expm1(-2.0 * absC)));
   EXPECT_TRUE(std::isfinite(psi));

   // Round-trip: tau_check = sigma * a * asinh(V/(2V0) * exp(psi/a))
   //           = sigma * a * arg_big   (modulo round-off).
   // Verify: psi/a + log(V/(2V0)) ≈ |arg_big| within 1e-12 abs.
   EXPECT_NEAR(psi / a_eq + std::log(1.0 / x), arg_big, 1e-12);
}
```
The test passes with the logsinh form and FAILS (non-finite or NaN) under the naive form.

---

### [R-004] [LOW] [drivers/tpv102_driver.cpp:114,162,221,272,512] — Stale comments referring to `EvaluateADERTotal` mislead readers

**Category:** QUALITY / documentation

**Description:**
Multiple driver comments and one banner string claim the production fault dispatch routes through `FaultFaceFlux::EvaluateADERTotal`. Audit:
- `tpv102_driver.cpp:114` — "wave.AdvanceADER -> FaultFaceFlux::EvaluateADERTotal (see the honest banner below)"
- `tpv102_driver.cpp:162` — "wave.AdvanceADER -> EvaluateADERTotal hard-codes the default Method::Brent argument"
- `tpv102_driver.cpp:221` — `BannerOf(DispatchedSolver::Brent)` literally returns "Brent (hard-coded via EvaluateADERTotal; --friction-solver flag IGNORED)"
- `tpv102_driver.cpp:272` and `512` — same misleading language.

Empirically (`grep -n "EvaluateADERTotal" miniapps/seas/dynamic/wave_operator.inl` returns zero matches), `wave.AdvanceADER` calls `fault_flux_->EvaluateADER(...)` at `wave_operator.inl:3614` — the fluctuation-Q path — never `EvaluateADERTotal`. The driver initializes `Q = 0` (`tpv102_setup.hpp:128–132 InitializeState`) so the fluctuation form is correct; it's the comment that is wrong. Banner output that mentions a function not on the call graph pollutes log diffs against TPV104 and seeds incorrect mental models.

**Trigger:**
Any reader (or future maintainer) who debugs the fault dispatch by `grep`ping for `EvaluateADERTotal` will find no matches in the production path and conclude the driver is broken.

**Actual behavior:**
Comments and banner text mention a function that is not on the call graph.

**Expected behavior:**
Comments should refer to `FaultFaceFlux::EvaluateADER` (fluctuation-Q variant). The misleading "EvaluateADERTotal" banner text must be replaced.

**Suggested fix:**
```diff
--- a/miniapps/seas/drivers/tpv102_driver.cpp
+++ b/miniapps/seas/drivers/tpv102_driver.cpp
@@
-// R7-001/R7-005 note: on the current driver path this value is kept only
-// for future iterator wiring.  The production time loop runs Brent via
-// wave.AdvanceADER -> FaultFaceFlux::EvaluateADERTotal (see the honest
-// banner below); the returned Method is not routed through that call.
+// R7-001/R7-005 note: on the current driver path this value is kept only
+// for future iterator wiring.  The production time loop runs Brent via
+// wave.AdvanceADER -> FaultFaceFlux::EvaluateADER (fluctuation-Q
+// dispatch — wave_operator.inl:3614); the returned Method is not routed
+// through that call.
@@
    case DispatchedSolver::Brent:
-      return "Brent (hard-coded via EvaluateADERTotal; "
+      return "Brent (hard-coded via EvaluateADER fluctuation-Q dispatch; "
              "--friction-solver flag IGNORED)";
```
Apply the same `s/EvaluateADERTotal/EvaluateADER/` rewrite at lines 162, 272, and 512 of `tpv102_driver.cpp`. Apply the matching fix to `tpv104_driver.cpp:113,161,219,270,511` (the same stale comments are inherited there).

**Test case:**
None required — documentation-only change. A grep-based regression guard is sufficient:
```bash
# expected output: zero matches once the cleanup lands
grep -n "EvaluateADERTotal" miniapps/seas/drivers/tpv102_driver.cpp \
                          miniapps/seas/drivers/tpv104_driver.cpp
```

---

### [R-005] [LOW] [drivers/tpv102_driver.cpp:1486–1493] — Substep deltaT seeded from initial `dt`, then rescaled per call — but driver claims "TPV102 uses fixed dt"

**Category:** QUALITY / latent assumption

**Description:**
The driver seeds the substep iterator's quadrature with `deltaT[o] = dt / O` once before the time loop (`tpv102_driver.cpp:1491`), then `AdvanceADERWithSubStep` (`tpv102_driver.cpp:328–334`) rescales `deltaT_scaled[o] = configured_deltaT[o] * (dt_step / configured_sum)` every macro-step. The comment at lines 1486–1490 claims:

> "Using `dt` (the auto-CFL initial value) here works because TPV102 uses fixed dt in the time loop."

That promise is violated by `tpv102_driver.cpp:1925`:
```cpp
real_t dt_step = std::min(dt, tfinal - t);
```
On the LAST step, `dt_step < dt` whenever `tfinal` is not an exact multiple of `dt`. The rescale at L328–334 protects `Σ deltaT[o] == dt_step`, so the iterator's per-call verify (`tpv102_substep_iterator.cpp:121–137`) won't trip — but the comment is wrong, and the rescale-each-call is actually load-bearing. This contradicts the rationale that motivated seeding the quadrature once.

This is the same pattern in TPV104 (`tpv104_driver.cpp:1493–1505`) — both drivers inherit the same misleading comment.

**Trigger:**
Run with `tfinal` not aligned to the auto-CFL `dt`. The last macro-step has `dt_step < dt` and `AdvanceADERWithSubStep` rescales the configured quadrature.

**Actual behavior:**
The rescale path is hit on the final step. No correctness issue (both ends of the verify hold). The misleading comment leads a reader to think the quadrature is bit-stable across the run.

**Expected behavior:**
Either drop the "TPV102 uses fixed dt" claim or honestly document that the rescale at L328–334 handles dt_step variability.

**Suggested fix:**
```diff
--- a/miniapps/seas/drivers/tpv102_driver.cpp
+++ b/miniapps/seas/drivers/tpv102_driver.cpp
@@
       // The Σ deltaT==dt_macro check inside Advance/AdvanceWithSubStepStates
-      // is RELATIVE so the same configuration handles every macro-step
-      // even though dt may vary slightly (it doesn't in TPV102, but the
-      // iterator is general).  Using `dt` (the auto-CFL initial value)
-      // here works because TPV102 uses fixed dt in the time loop.
+      // is RELATIVE; AdvanceADERWithSubStep rescales the configured
+      // deltaT to the actual dt_step at each call (see the dt_scale loop
+      // at L328–334).  Seeding with the auto-CFL `dt` is just a
+      // convenient initial scale.  The final macro-step (where
+      // dt_step = tfinal - t < dt) IS rescaled at runtime.
```

**Test case:**
None — comment-only change. A grep guard suffices:
```bash
grep -nE "uses fixed dt in the time loop" miniapps/seas/drivers/tpv102_driver.cpp \
                                          miniapps/seas/drivers/tpv104_driver.cpp
# expected: 0 matches after the fix.
```

---

### [R-006] [LOW] [POSSIBLE] [drivers/tpv102_driver.cpp:1808–1882; io/paraview_output.hpp:985–1006,1076–1083] — `pv_no_domain` branch advances regime via `CommitSchedule`; the other branch does not, leading to drifting regime state

**Category:** EDGE_CASE

**Description:**
The `paraview_write` lambda's branching (TPV102 driver L1863–1873) is:
```cpp
if (pv_no_domain) {
   pv_out->CommitSchedule(time);          // advances current_regime_ + last_write_time_ + last_v_max_
} else {
   pv_out->UpdateFaultFieldsBP5(...);
   pv_out->ForceSave(step_num, time);     // advances last_write_time_ only
}
pv_out->WriteFaultSurfaceVTU(...);        // unrelated to scheduling state
```

`ForceSaveImpl` (`paraview_output.hpp:1076–1083`) updates `last_write_time_` but NOT `current_regime_` or `last_v_max_`. Consequently, on the domain-PV-enabled branch, the V_max-driven regime never advances even though `PeekShouldWrite` consults `adaptive_.NextRegime(V_max, current_regime_)` for its `Interval(...)` decision (`paraview_output.hpp:940–948`). The schedule effectively freezes at whatever regime was returned by the last `Save`/`ShouldWrite`/`CommitSchedule`. For TPV102 this is dormant (the production sbatch passes `--paraview-dt`, which sets `fixed_dt > 0` and bypasses adaptive), but it IS reachable via the ParaViewOutput default schedule when the user passes only `--paraview` (no `--paraview-dt`, no `--paraview-every`) — and the driver's fallback at L1738 sets `output_every_n_steps = output_interval_for_pv` which routes through the `output_every_n_steps > 0` branch in `PeekShouldWrite` (`paraview_output.hpp:966–969`), bypassing the adaptive regime entirely. So the regime drift is **dormant** on TPV102 today but is a footgun for any future TPV102 sbatch that wants the V_max-adaptive cadence.

This is identical behaviour in TPV104 (`tpv104_driver.cpp:1875–1894`); not TPV102-specific, but it surfaces in this review because the user asked about workflow inconsistencies between the two drivers.

**Trigger:**
Pass only `--paraview` to either driver (no `--paraview-every`, no `--paraview-dt`, with `--no-domain-pv` OFF) AND configure the BP5 `AdaptiveSchedule` to non-step-based intervals. The first frame writes; subsequent frames consult a frozen regime.

**Actual behavior:**
`current_regime_` stays at 0 (interseismic), so `adaptive_.Interval(V_max, 0) = 1 yr`. No further frames write within a 12 s TPV102 run.

**Expected behavior:**
The `else` branch should also advance `current_regime_` / `last_v_max_` so the schedule stays consistent.

**Suggested fix:**
```diff
--- a/miniapps/seas/drivers/tpv102_driver.cpp
+++ b/miniapps/seas/drivers/tpv102_driver.cpp
@@
       if (pv_no_domain)
       {
          pv_out->CommitSchedule(time);
       }
       else
       {
          pv_out->UpdateFaultFieldsBP5(pv_local_slip, pv_local_slip_rate,
                                       pv_local_traction, pv_local_state,
                                       pv_local_normal_stress);
          pv_out->ForceSave(step_num, time);
+         // ForceSave only advances last_write_time_; the V_max-adaptive
+         // schedule additionally needs current_regime_ / last_v_max_
+         // advanced for the next PeekShouldWrite to use the correct
+         // regime interval (paraview_output.hpp:985–990).
+         pv_out->CommitSchedule(time, V_max);
       }
```
Apply the symmetric fix to `tpv104_driver.cpp:1875–1894`.

**Test case:**
```cpp
// tests/unit/test_paraview_regime_drift.cpp
TEST(R006_paraview_regime_advance_under_force_save) {
   ParaViewOutput<Mesh> pv("/tmp/test_pv", mesh, /*order=*/1);
   pv.GetSchedule().dt_coseismic   = 0.01;
   pv.GetSchedule().v_coseismic    = 1e-3;
   pv.GetSchedule().dt_interseismic = 1.0;
   // V_max ramps from 1e-9 to 1e-2; the adaptive regime should switch.
   EXPECT_TRUE (pv.PeekShouldWrite(0, 0.000, 1e-9));   // regime=0, interval=1s
   pv.ForceSave(0, 0.000);
   pv.CommitSchedule(0.000, 1e-9);                     // post-fix: applied here
   bool peek_v_co = pv.PeekShouldWrite(1, 0.005, 1e-2);
   // Pre-fix:  regime stays at 0, interval=1s, t<dt_out -> false.
   // Post-fix: regime moves to coseismic (V_max ≥ 1e-3), interval=0.01s -> true.
   EXPECT_TRUE(peek_v_co);
}
```

---

### [R-007] [LOW] [Makefile:1722,2756] — `seas_tpv102_driver` is compiled with `-DSEAS_USE_MPI` but `seas_tpv104_driver` is not — both rely on `MFEM_USE_MPI`, so the macro is dead

**Category:** QUALITY / build asymmetry

**Description:**
```
$(TPV102_DRIVER_OBJ): %.o: ...
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TPV104_DRIVER_OBJ): %.o: ...
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@
```

The `SEAS_USE_MPI` macro is referenced only by `io/parallel_benchmark_output.hpp` and `io/bp5_parallel_output.hpp` (BP5-specific). `tpv102_driver.cpp` does NOT use `SEAS_USE_MPI` — it gates everything on `MFEM_USE_MPI`. `grep -rn "SEAS_USE_MPI" miniapps/seas/drivers miniapps/seas/dynamic` returns zero matches in either driver. The `-DSEAS_USE_MPI` flag is dead weight on the TPV102 driver and the asymmetry vs TPV104 is silent — anyone bisecting build flags will be misled.

**Trigger:**
Build inconsistency search (e.g., `make print-FLAGS` or a CMake migration) trips on the asymmetric flag.

**Actual behavior:**
Same observable build, but the flag claims a feature that the driver does not consult.

**Expected behavior:**
Either remove `-DSEAS_USE_MPI` from the TPV102 driver rule (preferred — symmetry with TPV104) or add it to TPV104 with a comment explaining why both need it.

**Suggested fix:**
```diff
--- a/miniapps/seas/Makefile
+++ b/miniapps/seas/Makefile
@@
 $(TPV102_DRIVER_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV102_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
 	@mkdir -p $(@D)
-	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@
+	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@
```

**Test case:**
None — build-rule change. Confirm with:
```bash
grep -rn "SEAS_USE_MPI" miniapps/seas/drivers miniapps/seas/dynamic miniapps/seas/dynamic/*.inl
# expected: zero matches in TPV102/TPV104 paths.
```

---

## Summary
- Critical issues: 0
- Moderate issues: 2  ([R-001] friction-solver default on substep path; [R-002] dead `SEAS_DIAG_TPV102_STATE` banner)
- Low issues: 5  ([R-003] logsinh; [R-004] stale `EvaluateADERTotal`; [R-005] misleading "fixed dt" comment; [R-006] regime-drift on `ForceSave` branch; [R-007] dead `-DSEAS_USE_MPI`)
- Plan compliance: PARTIAL — TPV102 path is structurally a faithful TPV104 port; the substep-default and probe-banner contracts deviate.
- Verdict: PASS WITH FIXES — apply [R-001] and [R-002] before any TPV102 substep run; the rest are cleanups that prevent future regressions.

## Notes for the user

The user reported "TPV102 doesn't run." This review did **not** identify a runtime-fatal divergence on the *default* (one-shot, Brent, `--paraview-dt`) configuration that the active sbatch scripts use; the default path is a faithful TPV104 mirror. If the failure mode is reproducible only with `--fault-iterator substep`, [R-001] is almost certainly the cause (Newton-stable failing at TPV102's large ψ/a). If the failure mode is silent banner/contract mismatch in CI parsers, [R-002] is the cause. If the user has a concrete failure log (stack trace, NaN site, hang signature, sbatch output), please attach it to a follow-up so the next review can target the actual symptom rather than the structural diff.

## Unreviewed Areas
- `dynamic/wave_operator.inl` (3000+ lines) was treated as a black box — it is shared between BP5/TPV102/TPV104 and on the [C2] no-touch list per `feedback_dynamic_folder_editable_for_tpv104`. The review confirmed the fault-dispatch site (`wave_operator.inl:3614 EvaluateADER`) but did not audit the per-substep imposed-state side-channel beyond contract.
- `dynamic/fault_face_flux.{cpp,hpp}` — shared with BP5; treated as an API. Only its public surface (`ComputeStageState`, `BuildImposedState`, `WriteBackState`, `EvaluateADER`, `EvaluateADERTotal`) was matched against caller expectations.
- `mesh/tpv102_200m.msh` and other mesh artefacts — not in scope for a code review.
- The sbatch scripts beyond `tpv102_200m_p1_1.5s_400rank_dev.sbatch` (one was sampled for build/run-flag conventions).
- Unit tests `test_tpv102_total_locked_fault.cpp`, `test_tpv102_total_absorbing_equilibrium.cpp`, `test_tpv102_pepper_reproducer.cpp` — these still target the *legacy* total-Q path (`tpv102_setup_total.hpp`), which the new driver does not include. Whether these tests should be retired or ported to the fluctuation-Q + iterator flow is a separate scoping decision.
