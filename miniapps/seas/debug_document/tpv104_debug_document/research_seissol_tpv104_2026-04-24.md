# SeisSol vs MFEM — TPV104 Workflow Inconsistency Research Report

**Date:** 2026-04-24 (revised — reflects commits 4d2f9e8 + ffffc35: R7/R8 honest-dispatch driver + TPV102 pepper/BP5 support files)
**Objective:** Audit the *actual* TPV104 MFEM implementation (the `tpv104_*.hpp`/`.cpp` family, `SlipLawSRWPsi`, `FrictionCoefficientStable`, `Tpv104SubStepIterator`, `tpv104_driver.cpp`) against the SeisSol FL=103 reference, and document every inconsistency with real code, file paths, and line numbers.

**Critical scope update (commit 4d2f9e8)**: The production TPV104 driver runs under **R7-001 option (b)**: Brent + one-shot ADER + macro-step ψ + macro-step nucleation. The `Tpv104SubStepIterator`, `SlipLawSRWPsi::SetProductionMode`, and `SolveSlipRateNewtonStable` library code **exists and passes unit tests** but is **NOT wired into the driver's time loop** because doing so requires editing `dynamic/wave_operator.inl` (extreme-care, [C2] no-touch). The CLI flags `--friction-solver`, `--fault-iterator`, `--fric-law` are parsed and validated but do **not** reach the dispatched code path — they are banner-only. The `--verify-dispatch` flag emits `[dispatch]` lines disclosing the *actual* code path per rank.

This is surfaced explicitly in `drivers/tpv104_driver.cpp:104-108`:

```cpp
// Map the --friction-solver CLI name to FrictionSolver::Method (plan §4.10.X).
//
// R7-001/R7-005 note: on the current driver path this value is kept only
// for future iterator wiring.  The production time loop runs Brent via
// wave.AdvanceADER -> FaultFaceFlux::EvaluateADERTotal (see the honest
```

and again at `tpv104_driver.cpp:314-324`:

```cpp
   // R7-001/R7-003/R7-006: these flags are accepted so smoke tests and
   // sbatch scripts from the plan can pass CLI values, but on the current
   // driver path (one-shot wave.AdvanceADER) they have NO effect on the
   // dispatched solver / iterator / friction law.  Runtime is:
   //   - friction solver: Brent (hard-coded via EvaluateADERTotal)
   //   - fault iterator : one-shot (Tpv104SubStepIterator not wired
   //                      because doing so requires editing
   //                      wave_operator.inl, on the extreme-care list)
   //   - friction law   : slip-SRW via per-QP
   //                      UpdateStateAnalyticSlipLawSRW (per-macro-step)
```

This changes the category tags in §2: several items that would be CONSISTENT under the iterator path become DIFFERENT-EXPECTED or OPEN under option (b). Each affected entry is re-tagged below with an **R7-001(b) note**.

**MFEM TPV104 files audited** (under `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/`):

| File | Lines | Purpose |
|------|-------|---------|
| `config/tpv104_params.hpp` | 250 | SCEC parameter struct + spatial helpers + station table |
| `friction/slip_law_srw_psi.hpp` | 315 | FVW `SlipLawSRWPsi` class + `UpdateStateAnalyticSlipLawSRW` |
| `friction/friction_coeff_stable.hpp` | 159 | Stable `μ(V, ψ, a)` evaluator (replaces the `> 700` asymptotic branch) |
| `dynamic/tpv104_setup.hpp` | 541 | Fault-DOF init + V_w side-channel + station writers |
| `dynamic/tpv104_nucleation.hpp` | 143 | `SmoothStepIncrement_TPV104` + cumulative accumulator |
| `dynamic/tpv104_friction_solver.hpp` | 178 | Newton-Raphson solver (wraps stable μ) |
| `dynamic/tpv104_substep_iterator.hpp/.cpp` | 179 + 496 | Per-sub-step pipeline (friction + ψ + nucleation + imposed state) |
| `drivers/tpv104_driver.cpp` | 455 | Driver main — dispatches to `Tpv104SubStepIterator` |

**SeisSol reference** (under `/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/`):

* `FrictionLaws/CpuImpl/FastVelocityWeakeningLaw.h`
* `FrictionLaws/CpuImpl/RateAndState.h`
* `FrictionLaws/CpuImpl/BaseFrictionLaw.h`
* `FrictionLaws/FrictionSolverCommon.h`
* `FrictionLaws/RateAndStateCommon.h`
* `Initializer/BaseDRInitializer.cpp`
* `Initializer/RateAndStateInitializer.cpp`
* `Numerical/GaussianNucleationFunction.h`

**Conventions**:

* **CONSISTENT** — algorithm matches; Phase 3 probe agrees to ULP/1e-12.
* **CONSISTENT-WITH-DEVIATION** — algorithm matches but emission cadence or similar non-algorithmic detail differs; the comparator must compensate.
* **DIFFERENT-EXPECTED** — architectural difference that produces O(dt²) or bounded drift; probe thresholds must accept the delta.
* **DIFFERENT-REQUIRES-MAPPING** — permutation/sign remap applied by the probe-diff tool, no code change needed.
* **OPEN** — identified deviation that Phase 3 must flag to the user for decision.

---

## 1. Pipeline overview — what each code does per ADER macro-step

### SeisSol (`DynamicRupture/FrictionLaws/CpuImpl/BaseFrictionLaw.h:96-181`)

Per ADER macro-step, per fault cell:

1. `precomputeStressFromQInterpolated` — build `(σ_n, τ_1, τ_2)` at all `misc::TimeSteps` sub-steps from `qInterpolated±[o]`.
2. For `o = 0 .. TimeSteps-1`:
   1. `adjustInitialStress` — add `ΔS · smoothStepIncrement(t_end, dt, t0)` to `initialStressInFaultCS` (cumulative mutation).
   2. `updateFrictionAndSlip` — (calcInitialVariables, updateStateVariableIterative, calcSlipRateAndTraction): 2 outer × up-to-60 Newton iterations, analytic ψ update per iteration.
3. `postcomputeImposedStateFromNewStress` — accumulate `imposedState±` as `Σ_o timeWeights[o] · (normalStress, tractions)`.
4. `resampleStateVar` — project Δψ onto a lower-order basis (post-hook).

### MFEM TPV104 — PRODUCTION PATH (R7-001 option (b), `drivers/tpv104_driver.cpp:718-771`)

The production time loop is a **macro-step-only** pipeline. Per ADER macro-step:

1. `ApplyNucleationIncremental_TPV104(dof_data, fault_coords, t + dt_step, dt_step)` — ONE cumulative nucleation increment per macro-step (not per sub-step). `tpv104_driver.cpp:726-730`.
2. `wave.AdvanceADER(Q, dt_step, ader_order, Q_new)` — one-shot ADER predictor-corrector. Internally dispatches to `FaultFaceFlux::EvaluateADERTotal` which hard-codes `FrictionSolver::Method::Brent` and calls `Evaluate` ONCE on `Q̄ = I_±/dt_macro`. Writes V1/V2/slip_rate/tau*_corr/sigma_n_corr back onto DOFData. `tpv104_driver.cpp:732-739`.
3. Per-QP loop: `UpdateStateAnalyticSlipLawSRW` ONCE with `V = dof_data[i].slip_rate` (ADER-macro-averaged V), `dt = dt_step`; accumulate `slip1 += V1·dt_step`, `slip2 += V2·dt_step`. `tpv104_driver.cpp:756-771`.

Explicit R7-007 disclosure at `tpv104_driver.cpp:746-755`:

```cpp
      // R7-007 disclosure: plan §3.12 mandates PER-SUB-STEP ψ integration
      // ("follow exactly SeisSol did").  On this driver path ψ is
      // integrated ONCE per macro-step with V = dof_data[i].slip_rate
      // (ADER-averaged V over the whole dt_step).  The per-sub-step
      // cadence requires wiring Tpv104SubStepIterator into the time
      // loop (blocked by R7-001 option a), which in turn requires
      // exposing per-sub-step I± from wave_operator.inl (extreme-care
      // no-touch).  Under a rapidly-changing V the macro-step ψ deviates
      // from the per-sub-step result by O(dt_macro²); Phase-3 probe 2
      // against SeisSol will quantify the gap.
```

### MFEM TPV104 — LIBRARY PATH (shipped, not wired) (`dynamic/tpv104_substep_iterator.cpp:191-493`)

`Tpv104SubStepIterator::Advance` implements the per-sub-step pipeline that *would* run if the extreme-care wiring were permitted. It is linked, unit-tested (`test_tpv104_substep_iterator.cpp`), and invoked by probe-format tests only. Per sub-step it does:

1. `ApplyNucleationIncremental_TPV104` at sub-step endpoint (`tpv104_substep_iterator.cpp:295-296`).
2. Per QP: `ComputeStageState` on `Q̄ = I_±/dt_macro` (reused every sub-step — see I-05).
3. Per QP: accumulate `slip1 += V1·dt_sub`, `slip2 += V2·dt_sub`.
4. Per QP: `UpdateStateAnalyticSlipLawSRW` with per-QP `V_w[i]` + `a[i]`.
5. Per QP: `BuildImposedState` on `Q̄`.
6. Accumulate `I_imp_± += timeWeights[o] · dt_macro · Q_imp_±^{(o)}`.
7. On last sub-step: `WriteBackState` to DOFData.

### Key architectural differences (driver-path = production)

1. **Nucleation cadence**: SeisSol adds `smoothStepIncrement` per sub-step of the ADER macro-step; MFEM R7-001(b) driver adds ONCE per macro-step (`tpv104_driver.cpp:728-730`). For `dt_macro ≲ dt_ref_sub × O` the cumulative totals agree at end of ramp, but the per-probe snapshots differ at sub-step granularity.
2. **ψ cadence**: SeisSol runs `updateStateVariable` inside the per-sub-step friction loop (analytic exponential step per sub-step, with `numberStateVariableUpdates = 2` outer iterations using Kaneko averaging); MFEM R7-001(b) calls `UpdateStateAnalyticSlipLawSRW` ONCE per macro-step with the ADER-averaged V.
3. **Friction-solver dispatch**: SeisSol SIMD Newton (ψ re-updated inside the outer 2×iteration Kaneko loop); MFEM R7-001(b) dispatches **Brent** via `FaultFaceFlux::EvaluateADERTotal` → `Evaluate` → `FrictionSolver::Solve(..., Method::Brent)` (the CLI `--friction-solver newton-stable` is IGNORED on this path).
4. **Trial-traction cadence**: SeisSol uses `qInterpolated[o]` per sub-step; MFEM's `EvaluateADERTotal` uses `Q̄ = I_±/dt_macro` once.
5. **Analytic ψ formula**: `ψ_ss + (ψ_0 - ψ_ss) · exp(-V·dt/L)` byte-matches SeisSol when called with matched (V, ψ_0, dt) inputs — MFEM's `UpdateStateAnalyticSlipLawSRW` passes T_SRW_5 at 1e-13 relative.

---

## 2. Inconsistency inventory

**Category column reflects the production R7-001 option (b) driver path.** The parenthesised "library" category is what the tag *would* be if `Tpv104SubStepIterator` were wired in — provided for future Phase 2 scope.

| # | Title | Category (R7-001(b) driver) | Category (library if wired) |
|---|-------|----------------------------|------------------------------|
| I-01 | Trial traction precompute | CONSISTENT | CONSISTENT |
| I-02 | Friction coefficient stable asinh-exp (library only — NOT on driver path) | DIFFERENT-EXPECTED | CONSISTENT |
| I-03 | State-evolution ODE + analytic integrator (FVW) | CONSISTENT | CONSISTENT |
| I-04 | Nucleation cumulative accumulator cadence | DIFFERENT-EXPECTED (macro-step) | CONSISTENT (per-sub-step) |
| I-05 | Trial traction sub-step cadence — Q̄ reused vs qInterpolated[o] | CONSISTENT-WITH-DEVIATION | CONSISTENT-WITH-DEVIATION |
| I-06 | Friction solver — Brent (driver) vs Newton (library) vs Newton-SIMD (SeisSol) | DIFFERENT-EXPECTED | DIFFERENT-EXPECTED |
| I-07 | Slip-rate decomposition — parallel to total traction | CONSISTENT | CONSISTENT |
| I-08 | Corrected traction (trial scale) to Riemann | CONSISTENT | CONSISTENT |
| I-09 | Imposed-state time-weighted accumulation | DIFFERENT-EXPECTED (single-shot) | CONSISTENT (O-weighted sum) |
| I-10 | Pre-stress storage — side-channel (both) | CONSISTENT | CONSISTENT |
| I-11 | Tangent-frame convention (t1/t2) | DIFFERENT-REQUIRES-MAPPING | DIFFERENT-REQUIRES-MAPPING |
| I-12 | Normal-stress sign convention | DIFFERENT-REQUIRES-MAPPING | DIFFERENT-REQUIRES-MAPPING |
| I-13 | Initial ψ inversion | CONSISTENT | CONSISTENT |
| I-14 | `f_LV` `max(0, ...)` clamp sign | OPEN | OPEN |
| I-15 | Default solver dispatch — Brent (driver) vs Newton-stable (library) vs Newton-SIMD (SeisSol) | DIFFERENT-EXPECTED | CONSISTENT |
| I-16 | State-variable resampling (postHook) | OPEN | OPEN |
| I-17 | ψ integration cadence — macro-step vs per-sub-step | DIFFERENT-EXPECTED | CONSISTENT |
| I-18 | Parameter values | CONSISTENT | CONSISTENT |
| I-19 | Mesh / ADER order / station coordinates | DIFFERENT-EXPECTED | DIFFERENT-EXPECTED |
| I-20 | CLI flags parsed but not routed (R7-001(b)) | DIFFERENT-EXPECTED (banner-only) | N/A |

---

## 3. Inconsistency entries

### I-01 Trial traction precompute [CONSISTENT]

**What SeisSol does** — `DynamicRupture/FrictionLaws/FrictionSolverCommon.h:180-193`:

```cpp
    for (auto index = Range::Start; index < Range::End; index += Range::Step) {
      auto i{startLoopIndex + index};
      VariableIndexing<RangeExecutor<Type>::Exec>::index(faultStresses.normalStress, o, i) =
          etaP * (qIMinus[o][U][i] - qIPlus[o][U][i] + qIPlus[o][N][i] * invZp +
                  qIMinus[o][N][i] * invZpNeig);

      VariableIndexing<RangeExecutor<Type>::Exec>::index(faultStresses.traction1, o, i) =
          etaS * (qIMinus[o][V][i] - qIPlus[o][V][i] + qIPlus[o][T1][i] * invZs +
                  qIMinus[o][T1][i] * invZsNeig);

      VariableIndexing<RangeExecutor<Type>::Exec>::index(faultStresses.traction2, o, i) =
          etaS * (qIMinus[o][W][i] - qIPlus[o][W][i] + qIPlus[o][T2][i] * invZs +
                  qIMinus[o][T2][i] * invZsNeig);
    }
```

**What MFEM does** — reused verbatim from `dynamic/fault_face_flux.cpp:40-66` (shared component, not edited for TPV104). The TPV104 sub-step iterator calls it via `FaultFaceFlux::ComputeStageState` at `dynamic/tpv104_substep_iterator.cpp:319-320`:

```cpp
         EvalStageState s;
         flux_.ComputeStageState(d, Q_avg_plus, Q_avg_minus, s, method);
```

Impedance fields per TPV104 QP — `dynamic/tpv104_setup.hpp:100-106`:

```cpp
      // Impedances (homogeneous half-space — ρ, cp, cs from TPV104Params).
      d.Zp_plus  = TPV104Params::Zp;
      d.Zp_minus = TPV104Params::Zp;
      d.Zs_plus  = TPV104Params::Zs;
      d.Zs_minus = TPV104Params::Zs;
      d.eta_p    = TPV104Params::eta_p;
      d.eta_s    = TPV104Params::eta_s;
```

with `TPV104Params::eta_s = Zs / 2.0` (`config/tpv104_params.hpp:49`). For homogeneous material SeisSol's `etaS = Zs·Zs/(2·Zs) = Zs/2`, identical.

**Implication**: Byte-for-byte algebraic parity on trial-traction precompute. Phase 3 probe 1 (stage A) expects ULP agreement on a single-QP fixture. **No action.**

---

### I-02 Friction coefficient stable asinh-exp [DIFFERENT-EXPECTED on driver, CONSISTENT in library]

**R7-001(b) note**: `FrictionCoefficientStable` is exercised by `SolveSlipRateNewtonStable` and the `Tpv104SubStepIterator` — neither is on the production driver path. The production path calls `FaultFaceFlux::EvaluateADERTotal` → `Evaluate` → `friction/dieterich_ruina.hpp::FrictionCoefficientPsi` (the legacy MFEM evaluator with the `ψ/a > 700` asymptotic branch). For TPV104's envelope `ψ/a ≤ 80`, both evaluators give ULP-equivalent results.


**What SeisSol does** — `DynamicRupture/FrictionLaws/RateAndStateCommon.h:60-87`:

```cpp
template <typename T>
SEISSOL_HOSTDEVICE constexpr T arsinhexp(T x, T expLog, T exp) {
  constexpr T Switch = 10;
  constexpr T Threshold = 50;
  constexpr T Log2 = 0.69314718055994530943;
  int xexp{};
  (void)std::frexp(x, &xexp);

  if (expLog + std::max(xexp, 0) * Log2 > Switch || expLog >= Threshold) {
    if (expLog <= 0) {
      exp = 1 / exp;
    }
    const T xa = std::abs(x);
    const T xs = x >= 0 ? 1 : -1;
    return xs * (expLog + std::log(xa + std::sqrt(xa * xa + exp * exp)));
  } else {
    if (expLog > 0) {
      exp = 1 / exp;
    }
    const auto v = exp * x;
    return std::asinh(v);
  }
}
```

Used as `μ = a · arsinhexp(V/(2V₀), ψ/a, cExp)` in `FastVelocityWeakeningLaw.h:118-122`.

**What MFEM does** — `friction/friction_coeff_stable.hpp:74-101` is a **verbatim port**:

```cpp
inline real_t ArsinhExp(real_t x, real_t cExpLog, real_t cExp)
{
   int xexp = 0;
   (void)std::frexp(x, &xexp);

   if (cExpLog + std::max(xexp, 0) * kLog2 > kSwitch ||
       cExpLog >= kThreshold)
   {
      // Stable branch:
      //   asinh(x·e^c) = sign(x) · (c + log(|x| + sqrt(x² + e^{-2c}))).
      // When cExpLog ≤ 0, ComputeCExp returned exp(cExpLog) — we must
      // flip to exp(-cExpLog) for the formula.
      real_t exp_neg = cExp;
      if (cExpLog <= 0.0) { exp_neg = 1.0 / exp_neg; }
      const real_t xa = std::abs(x);
      const real_t xs = (x >= 0.0) ? 1.0 : -1.0;
      return xs * (cExpLog + std::log(xa + std::sqrt(xa * xa
                                                     + exp_neg * exp_neg)));
   }
   else
   {
      // Normal branch: v = exp(cExpLog) · x; return asinh(v).
      real_t exp_pos = cExp;
      if (cExpLog > 0.0) { exp_pos = 1.0 / exp_pos; }
      const real_t v = exp_pos * x;
      return std::asinh(v);
   }
}
```

with `kSwitch = 10.0`, `kThreshold = 50.0`, `kLog2 = 0.69314718055994530943` (lines 49-51) — identical constants. Wrapper:

```cpp
// friction/friction_coeff_stable.hpp:133-141
inline real_t FrictionCoefficientStable(real_t V, real_t psi, real_t a,
                                        real_t V0)
{
   const real_t cLin    = 0.5 / V0;                        // 1 / (2 V0)
   const real_t cExpLog = psi / a;
   const real_t cExp    = ComputeCExp(cExpLog);
   const real_t lx      = cLin * V;
   return a * ArsinhExp(lx, cExpLog, cExp);
}
```

matches SeisSol's `getMuDetails` + `updateMu` combo verbatim.

**Implication**: Identical formula, identical branch thresholds, identical `frexp`-based switch. ULP-level agreement expected on TPV104's (V, ψ, a) envelope. Probe 3 (friction coefficient) target tolerance 1e-12 relative. **No action.**

---

### I-03 State-evolution ODE + analytic integrator [CONSISTENT]

**What SeisSol does** — `DynamicRupture/FrictionLaws/CpuImpl/FastVelocityWeakeningLaw.h:54-77`:

```cpp
    const real lowVelocityFriction =
        std::max(static_cast<real>(0),
                 static_cast<real>(this->f0[faceIndex][pointIndex] -
                                   (this->b[faceIndex][pointIndex] - localA) *
                                       log(localSlipRate / this->drParameters->rsSr0)));
    const real steadyStateFrictionCoefficient =
        localMuW + (lowVelocityFriction - localMuW) /
                       std::pow(1.0 + misc::power<8, double>(localSlipRate / localSrW), 1.0 / 8.0);
    const real steadyStateStateVariable =
        localA * rs::logsinh(this->drParameters->rsSr0 / localSlipRate * 2,
                             steadyStateFrictionCoefficient / localA);

    // exact integration of dSV/dt DGL, assuming constant V over integration step

    const auto preexp1 = -localSlipRate * (timeIncrement / localSl0);
    const real exp1v = std::exp(preexp1);
    const real exp1m = -std::expm1(preexp1);
    const real localStateVariable = steadyStateStateVariable * exp1m + exp1v * stateVarReference;
```

with `rs::logsinh` at `RateAndStateCommon.h:146-152`:

```cpp
template <typename T>
SEISSOL_HOSTDEVICE constexpr T logsinh(T x, T c) {
  const T sign = c >= 0 ? 1 : -1;
  const T absC = std::abs(c);
  return absC + std::log(x / 2 * -sign * std::expm1(-2 * absC));
}
```

and `misc::power<8, double>(r)` is a template-unrolled `(((r·r)·r)·…)` integer power chain.

**What MFEM does** — `friction/slip_law_srw_psi.hpp:92-124` is a **verbatim port**:

```cpp
inline real_t UpdateStateAnalyticSlipLawSRW(real_t psi_old, real_t V,
                                            real_t L, real_t dt,
                                            real_t V_w, real_t a,
                                            real_t b, real_t V0,
                                            real_t f0, real_t muW)
{
   // (1) f_LV = max(0, f0 - (b - a) · log(V/V0))
   //     Sign convention: matches the canonical FVW reference and
   //     SCEC TPV104 (Noda & Lapusta 2013).  The plan's §4.2.1 eq (2a)
   //     has a sign typo; the byte-match directive resolves it in favour
   //     of the reference implementation.
   const real_t f_LV = std::max(static_cast<real_t>(0),
                                f0 - (b - a) * std::log(V / V0));

   // (2) f_ss = muW + (f_LV - muW) / (1 + (V/V_w)^8)^(1/8)
   //     R-004: integer 8th power via unrolled binary chain, not std::pow.
   const real_t V_over_Vw   = V / V_w;
   const real_t V_over_Vw_8 = IntegerPow8(V_over_Vw);
   const real_t denom       = std::pow(1.0 + V_over_Vw_8,
                                       static_cast<real_t>(1.0 / 8.0));
   const real_t f_ss        = muW + (f_LV - muW) / denom;

   // (3) ψ_ss = a · logsinh(2·V0/V, f_ss/a)
   const real_t psi_ss = a * LogSinhStable(2.0 * V0 / V, f_ss / a);

   // (4) Analytic step:  ψ(t+dt) = ψ_ss + (ψ_0 - ψ_ss)·exp(-V·dt/L).
   //     Expressed as ψ_ss·(1 - exp(preexp1)) + exp(preexp1)·ψ_0 to match
   //     the canonical reference's round-off profile.
   const real_t preexp1 = -V * (dt / L);
   const real_t exp1v   = std::exp(preexp1);
   const real_t exp1m   = -std::expm1(preexp1);
   return psi_ss * exp1m + exp1v * psi_old;
}
```

with `IntegerPow8` the binary-chain unrolling (lines 63-68):

```cpp
inline real_t IntegerPow8(real_t x)
{
   const real_t x2 = x * x;
   const real_t x4 = x2 * x2;
   return x4 * x4;
}
```

and `LogSinhStable` the port of SeisSol's `logsinh` (lines 51-56):

```cpp
inline real_t LogSinhStable(real_t x, real_t c)
{
   const real_t sign_c = (c >= 0.0) ? 1.0 : -1.0;
   const real_t absC   = std::abs(c);
   return absC + std::log(x / 2.0 * -sign_c * std::expm1(-2.0 * absC));
}
```

**Implication**: Byte-match to the SeisSol reference for identical inputs. T_SRW_5 unit test enforces 1e-13 relative on 50 random input tuples. Probe 2 (state evolution) expects ULP agreement when inputs agree.

**Note**: `V_over_Vw_8` uses integer `IntegerPow8` (exactly ULP-identical to `misc::power<8, double>`), but the outer `(1 + X)^{1/8}` call uses `std::pow` rather than SeisSol's `std::pow` — same library call; ULP agreement holds.

---

### I-04 Nucleation cumulative accumulator + SmoothStepIncrement [DIFFERENT-EXPECTED on driver, CONSISTENT in library]

**R7-001(b) note**: The production driver calls `ApplyNucleationIncremental_TPV104(..., t + dt_step, dt_step)` ONCE per macro-step at `tpv104_driver.cpp:726-730`, not per sub-step. SeisSol calls `adjustInitialStress` inside the per-sub-step loop of `BaseFrictionLaw::evaluate`. The cumulative smoothStep increment telescopes to the same full-ramp value at the end regardless of cadence, so end-state `tau2_nuc = Δτ₀·F(r)` matches; transient per-ADER-stage snapshots differ at sub-step granularity. The library iterator `tpv104_substep_iterator.cpp:295-296` does call per-sub-step.


**What SeisSol does** — `Numerical/GaussianNucleationFunction.h:22-39`:

```cpp
template <typename T>
SEISSOL_HOSTDEVICE inline T smoothStep(T currentTime, T t0) {
  if (currentTime <= 0) {
    return 0.0;
  } else if (currentTime < t0) {
    const T tau = currentTime - t0;
    return std::exp(tau * tau / (currentTime * (currentTime - 2.0 * t0)));
  } else {
    return 1.0;
  }
}

template <typename T>
SEISSOL_HOSTDEVICE inline T smoothStepIncrement(T currentTime, T dt, T t0) {
  return smoothStep<T>(currentTime, t0) - smoothStep<T>(currentTime - dt, t0);
}
```

and the per-sub-step call at `FrictionLaws/FrictionSolverCommon.h:432-450`:

```cpp
  if (fullUpdateTime <= t0 + s0 && fullUpdateTime >= s0) {
    const real gNuc =
        gaussianNucleationFunction::smoothStepIncrement<real>(fullUpdateTime - s0, dt, t0);

    for (auto index = Range::Start; index < Range::End; index += Range::Step) {
      auto pointIndex{startIndex + index};
      for (unsigned i = 0; i < 6; i++) {
        initialStressInFaultCS[i][pointIndex] += nucleationStressInFaultCS[i][pointIndex] * gNuc;
      }
      initialPressure[pointIndex] += nucleationPressure[pointIndex] * gNuc;
    }
  }
```

**What MFEM does** — `dynamic/tpv104_nucleation.hpp:40-68` ports the smooth-step pair verbatim:

```cpp
inline real_t SmoothStep_TPV104(real_t current_time, real_t t0)
{
   if (current_time <= 0.0) { return 0.0; }
   if (current_time < t0)
   {
      const real_t tau = current_time - t0;
      return std::exp(tau * tau / (current_time * (current_time - 2.0 * t0)));
   }
   return 1.0;
}

inline real_t SmoothStepIncrement_TPV104(real_t current_time, real_t dt,
                                         real_t t0)
{
   MFEM_ASSERT(std::isfinite(dt) && dt > 0.0,
               "SmoothStepIncrement_TPV104: dt must be finite and "
               "positive; got dt = " << dt);
   return SmoothStep_TPV104(current_time, t0)
          - SmoothStep_TPV104(current_time - dt, t0);
}
```

and the per-sub-step accumulator at `dynamic/tpv104_nucleation.hpp:104-138`:

```cpp
inline void ApplyNucleationIncremental_TPV104(
   std::vector<DOFData> &dof_data,
   const std::vector<Vector> &fault_coords,
   real_t t_substep_end,
   real_t dt_substep)
{
   const int n = static_cast<int>(dof_data.size());
   MFEM_VERIFY(static_cast<int>(fault_coords.size()) >= n,
               "ApplyNucleationIncremental_TPV104: fault_coords size "
               << fault_coords.size() << " < dof_data size " << n);

   const real_t dS = SmoothStepIncrement_TPV104(t_substep_end, dt_substep,
                                                TPV104Params::nuc_T);
   // R3-004 (review round 3): clamp non-positive increments to 0.
   // `smoothStep` is monotonically non-decreasing so mathematically
   // dS ≥ 0, but floating-point round-off at the ramp boundaries can
   // produce a sub-ULP negative dS (two close-to-1 values subtracted
   // via different code paths).  `dS == 0.0` exact-equality would miss
   // this; `dS <= 0.0` both fast-paths the before/after-ramp branches
   // and guards against sub-ULP negative drift into the accumulator.
   if (dS <= 0.0) { return; }

   for (int i = 0; i < n; ++i)
   {
      const real_t along_strike = fault_coords[i](0);
      const real_t down_dip     = std::abs(fault_coords[i](2));
      const real_t dx = along_strike - TPV104Params::hypo_along_strike;
      const real_t dz = down_dip     - TPV104Params::hypo_down_dip;
      const real_t r  = std::sqrt(dx * dx + dz * dz);
      const real_t F  = NucleationSpatial_TPV104(r);
      dof_data[i].tau2_nuc += dS * F * TPV104Params::nuc_dtau;
      // tau1_nuc and sigma_n_nuc are intentionally not updated
      // (pure strike-slip; no normal-stress nucleation in TPV104).
   }
}
```

Called per sub-step at `dynamic/tpv104_substep_iterator.cpp:295-296`:

```cpp
      ApplyNucleationIncremental_TPV104(dof_data, fault_coords,
                                        t_substep_end, dt_substep);
```

**Implication**: Identical cumulative smoothStep-increment architecture. Both codes produce `tau2_nuc = Δτ₀·F(r)·g(t)` at the end of the ramp by telescoping the increments. This is the item where the **TPV102 overwrite-per-macro-step path** diverged from SeisSol; TPV104 fixed it by ADDING per sub-step. **No action.**

**Minor difference**: SeisSol adds the increment to `initialStressInFaultCS[3]` (the tangent1 shear — strike in SeisSol's convention, see I-11); MFEM adds to `data.tau2_nuc` (strike in MFEM's convention). The underlying physical field (strike-shear traction increment) is the same.

---

### I-05 Trial traction sub-step cadence — Q̄ reused vs `qInterpolated[o]` [CONSISTENT-WITH-DEVIATION]

**What SeisSol does** — `DynamicRupture/FrictionLaws/CpuImpl/BaseFrictionLaw.h:72-77` runs `precomputeStressFromQInterpolated` ONCE outside the sub-step loop:

```cpp
        common::precomputeStressFromQInterpolated(faultStresses,
                                                  impAndEta[ltsFace],
                                                  impedanceMatrices[ltsFace],
                                                  qInterpolatedPlus[ltsFace],
                                                  qInterpolatedMinus[ltsFace],
                                                  etaPDamp);
```

but it fills `faultStresses.normalStress[o][i]` for **all** sub-steps `o ∈ [0, misc::TimeSteps)` by reading `qInterpolated±[o]` — SeisSol precomputes Q at each ADER time-quadrature point separately and thus has O distinct trial tractions.

**What MFEM does** — `dynamic/tpv104_substep_iterator.cpp:275-311` time-averages the predictor ONCE and reuses it at every sub-step:

```cpp
   // Cache per-QP time-averaged bulk state once — Q_avg is constant
   // across sub-steps under our interpretation of the predictor (§4.10
   // Step 7 ambiguity resolution).  Stack-allocated per-DOF working
   // buffers are tiny (9 real_t each); a single heap buffer of size
   // NUM_STATE keeps the per-sub-step loop allocation-free.
   real_t Q_avg_plus[NUM_STATE];
   real_t Q_avg_minus[NUM_STATE];
   // ...
   for (int o = 0; o < O; ++o)
   {
      // ...
      for (int i = 0; i < n; ++i)
      {
         DOFData      &d  = dof_data[i];
         const real_t *Ip = I_plus_flat  + static_cast<ptrdiff_t>(i) * NUM_STATE;
         const real_t *Im = I_minus_flat + static_cast<ptrdiff_t>(i) * NUM_STATE;

         // Time-average predictor -> pointwise Q̄ for friction solve.
         for (int c = 0; c < NUM_STATE; ++c)
         {
            Q_avg_plus[c]  = Ip[c] * inv_dt_macro;
            Q_avg_minus[c] = Im[c] * inv_dt_macro;
         }

         // Friction stage chain.
         EvalStageState s;
         flux_.ComputeStageState(d, Q_avg_plus, Q_avg_minus, s, method);
```

So MFEM evaluates the trial traction on **the same Q̄ at every sub-step** — the trial traction values `s.sigma_n_trial, s.tau1_trial, s.tau2_trial` are identical across the O sub-steps of one macro-step.

This is known and disclosed at `dynamic/tpv104_substep_iterator.cpp:102-109`:

```cpp
         if (probe_name == std::string("trial_traction"))
         {
            (*f) << "# cadence_note: MFEM emits O identical rows per "
                    "macro-step (Q̄ = I/dt_macro used for every "
                    "sub-step); reference emits O distinct rows. "
                    "Step-13 probe_diff must coarsen by sub-step "
                    "averaging or mask before comparing.\n";
         }
```

and again at lines 326-338:

```cpp
         // R4-004 deviation disclosure (review round 4): the MFEM
         // iterator evaluates the trial traction at Q̄ = I_±/dt_macro
         // for every sub-step o (see §4.10 Step 7 ambiguity resolution
         // in the header docstring).  The reference runtime evaluates
         // at per-sub-step qInterpolated[o].  Probe-1 rows emitted by
         // MFEM are therefore IDENTICAL across the O sub-steps of one
         // macro-step, while the reference emits O distinct rows.
         // The Step-13 Python probe-diff must mask (or coarsen by
         // sub-step-averaging) this channel when comparing across
         // codes until the iterator routes the per-sub-step ADER
         // predictor through the trial traction call (a Step-9
         // driver-level extension that requires editing the ADER
         // predictor interface — currently blocked by [C2]).
```

**Implication**: The integrated imposed state agrees to O(dt²) (Simpson-rule argument on the Lipschitz-smooth residual — same argument used for MFEM's single-shot ADER in TPV102). But the **per-sub-step probe-1 dumps** are **not** comparable row-for-row. The `tpv104_column_map.py` probe-diff tool must:

1. Either **mask** probe-1 rows at sub-step cadence, or
2. **Average** SeisSol's O rows per macro-step before diffing against MFEM's single row.

The macro-step-integrated imposed-state comparison (probe 5) remains valid.

**Action**: Probe-diff script must implement the mask/coarsen. The plan §4.10 Step 9 notes that routing per-sub-step Q through the predictor would require editing `fault_face_flux.hpp` (Extreme Care, [C2]) — **deferred**.

---

### I-06 Friction solver outer loop — 2× Kaneko vs single call [DIFFERENT-EXPECTED]

**What SeisSol does** — `DynamicRupture/FrictionLaws/CpuImpl/RateAndState.h:151-195`:

```cpp
  void updateStateVariableIterative(
      bool& hasConverged, /* ... */) {
    std::array<real, misc::NumPaddedPoints> testSlipRate{0};
    for (uint32_t j = 0; j < settings.numberStateVariableUpdates; j++) {
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
        // update state variable using sliprate from the previous time step
        localStateVariable[pointIndex] =
            static_cast<Derived*>(this)->updateStateVariable(pointIndex,
                                                             ltsFace,
                                                             stateVarReference[pointIndex],
                                                             this->deltaT[timeIndex],
                                                             localSlipRate[pointIndex]);
      }
      // ...
      hasConverged = this->invertSlipRateIterative(/* ... */);

      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
        // For the next SV update, use the mean slip rate between the initial guess and the one
        // found (Kaneko 2008, step 6)
        localSlipRate[pointIndex] = 0.5 * (this->slipRateMagnitude[ltsFace][pointIndex] +
                                           std::fabs(testSlipRate[pointIndex]));
        // ...
      }
    }
  }
```

So per sub-step: `numberStateVariableUpdates = 2` outer iterations × (analytic ψ update + Newton-on-V with averaging).

**What MFEM does** — `dynamic/tpv104_substep_iterator.cpp:319-320` then `:413-419`:

```cpp
         EvalStageState s;
         flux_.ComputeStageState(d, Q_avg_plus, Q_avg_minus, s, method);
         // ...
         d.psi = UpdateStateAnalyticSlipLawSRW(d.psi, s.V_abs,
                                               d.Dc, dt_sub,
                                               V_w[i], d.a,
                                               state_evo_.GetB(),
                                               state_evo_.GetV0(),
                                               state_evo_.GetF0(),
                                               state_evo_.GetMuW());
```

The Newton-on-V call inside `ComputeStageState` runs ONCE per sub-step, using the ψ carried in from the previous sub-step (or the init ψ at t=0). Then ψ is updated once analytically using the just-solved V.

**Implication**:

* Both codes converge on the same (V, ψ) pair at equilibrium; the deviation is in **transient dynamics** during ψ's rapid change.
* MFEM uses the `NewtonRaphsonStable` solver by default (`drivers/tpv104_driver.cpp:236` + `tpv104_substep_iterator.hpp:152-153`), with the same residual form as SeisSol's solver but inside `ComputeStageState` and using the stable-asinh μ.
* For Lipschitz-smooth ODE behaviour in TPV104 (nucleation smoothly ramps V from 1e-16 to ~5 m/s over ~0.5 s) the Kaneko averaging is largely cosmetic — both codes converge to the same root.
* Phase 3 Probe 4 (V_abs) is expected to agree to 1e-8 relative, consistent with the 1e-8 `newtonTolerance` in both codes.

**Action**: None. Document the Kaneko-averaging difference in the probe-diff report.

---

### I-07 Slip-rate decomposition — parallel to total traction [CONSISTENT]

**What SeisSol does** — `DynamicRupture/FrictionLaws/CpuImpl/RateAndState.h:227-232`:

```cpp
      const auto divisor =
          strength + this->impAndEta[ltsFace].etaS * this->slipRateMagnitude[ltsFace][pointIndex];
      this->slipRate1[ltsFace][pointIndex] =
          this->slipRateMagnitude[ltsFace][pointIndex] * totalTraction1 / divisor;
      this->slipRate2[ltsFace][pointIndex] =
          this->slipRateMagnitude[ltsFace][pointIndex] * totalTraction2 / divisor;
```

**What MFEM does** — shared `FaultFaceFlux::CompleteFromVabs` at `dynamic/fault_face_flux.cpp:143-145`:

```cpp
      s.V1 = s.V_abs * (s.tau1_total) / (strength + data.eta_s * s.V_abs);
      s.V2 = s.V_abs * (s.tau2_total) / (strength + data.eta_s * s.V_abs);
```

where `tau_{1,2}_total = tau_{1,2}_0 + tau_{1,2}_nuc + tau_{1,2}_trial` from `fault_face_flux.cpp:97-99`. For TPV104 the persistent nucleation channel `data.tau2_nuc` is written by `ApplyNucleationIncremental_TPV104` (I-04) each sub-step.

**Implication**: Identical: `V_i = |V|·τ_i^tot/(strength + η_s·|V|)`. Both parallel to total traction. **No action.**

---

### I-08 Corrected traction (trial scale) to Riemann [CONSISTENT]

**What SeisSol does** — `RateAndState.h:235-240`:

```cpp
      tractionResults.traction1[timeIndex][pointIndex] =
          faultStresses.traction1[timeIndex][pointIndex] -
          this->impAndEta[ltsFace].etaS * this->slipRate1[ltsFace][pointIndex];
      tractionResults.traction2[timeIndex][pointIndex] =
          faultStresses.traction2[timeIndex][pointIndex] -
          this->impAndEta[ltsFace].etaS * this->slipRate2[ltsFace][pointIndex];
```

**What MFEM does** — `fault_face_flux.cpp:147-149`:

```cpp
      s.tau1_corr = s.tau1_trial - data.eta_s * s.V1;
      s.tau2_corr = s.tau2_trial - data.eta_s * s.V2;
```

Same trial-scale — no pre-stress/nucleation baseline re-added — so `tau_corr` fed to `BuildImposedState` is directly comparable to SeisSol's `tractionResults`. **No action.**

---

### I-09 Imposed-state time-weighted accumulation [DIFFERENT-EXPECTED on driver, CONSISTENT in library]

**R7-001(b) note**: The production driver runs `wave.AdvanceADER` → `FaultFaceFlux::EvaluateADERTotal` → `EvaluateTotal` on `Q̄ = I_±/dt_macro` ONCE per macro-step. The imposed state is a single `(1.0) · dt_macro · Q_imp_±` — equivalent to the library path only in the `O = 1` single-sub-step limit with `time_weights = {1.0}`. For SeisSol's genuine O-sub-step `Σ_o timeWeights[o] · Q_imp_±^{(o)}` accumulation, the MFEM production matches SeisSol to O(dt_macro²) by the Simpson-rule / Lipschitz-smooth-residual argument — same O(dt²) that TPV102's `EvaluateADER` (§3.8 of the debug plan) was justified under.


**What SeisSol does** — `FrictionSolverCommon.h:329-362` accumulates across `misc::TimeSteps` sub-steps using `timeWeights[o]` (see §1 overview).

**What MFEM does** — `dynamic/tpv104_substep_iterator.cpp:470-480`:

```cpp
         // Accumulate time-weighted imposed state in time-integrated
         // form: I_imp_± += timeWeights[o] · dt_macro · Q_imp_±^{(o)}.
         real_t *Iout_p = I_imp_plus_flat
                          + static_cast<ptrdiff_t>(i) * NUM_STATE;
         real_t *Iout_m = I_imp_minus_flat
                          + static_cast<ptrdiff_t>(i) * NUM_STATE;
         for (int c = 0; c < NUM_STATE; ++c)
         {
            Iout_p[c] += accum_scale * Q_imp_plus[c];
            Iout_m[c] += accum_scale * Q_imp_minus[c];
         }
```

with `accum_scale = weight · dt_macro` (line 289).

**Implication**: Same time-weighted accumulation pattern; both emit `I_imp_± = Σ_o w_o·dt·Q_imp_±^{(o)}`. The difference is in the inputs `Q_imp_±^{(o)}` — SeisSol computes them on `qInterpolated[o]` while MFEM uses `Q̄` for every `o` (see I-05). But the final `I_imp_±` passed back to the bulk integrator follows the same scaling convention. **No action**, given I-05 caveat.

---

### I-10 Pre-stress storage — side-channel (both) [CONSISTENT]

**What SeisSol does** — `initialStressInFaultCS[ltsFace][6][NumPaddedPoints]` holds the 6-component pre-stress in fault-local CS; bulk Q at t=0 is zero.

**What MFEM does** — `dynamic/tpv104_setup.hpp:112-123`:

```cpp
      d.sigma_n0 = TPV104Params::sigma_n;
      d.tau1_0   = 0.0;
      d.tau2_0   = TPV104Params::tau_ini;

      d.sigma_n_nuc = 0.0;
      d.tau1_nuc    = 0.0;
      d.tau2_nuc    = 0.0;
```

and `dynamic/tpv104_setup.hpp:167-174`:

```cpp
inline void InitializeState_TPV104(Vector &Q, int ndof_total)
{
   // ...
   Q.SetSize(NUM_STATE * ndof_total);
   Q = 0.0;
}
```

TPV104 is **fluctuation-Q only** — §3.10 directive is explicit at `dynamic/tpv104_setup.hpp:21-26`:

```cpp
// Fluctuation-Q-only (§3.10 directive): this header ships ONLY the
// mode-1 fluctuation-Q initialiser.  There is NO `InitializeStateTotal_TPV104`
// or `ApplyNucleationTotalPrestress_TPV104` — total-Q mode is dropped
// for TPV104 so bulk Q carries fluctuation only and pre-stress lives
// exclusively in the fault-CS side-channel (DOFData.sigma_n0, tau*_0).
// T_TPV104_SETUP_5 enforces this at link/grep time.
```

**Implication**: TPV104 matches SeisSol's single-mode architecture exactly: pre-stress in side-channel, bulk Q carries fluctuation. The dual-mode baggage from TPV102 is dropped. **No action.**

---

### I-11 Tangent-frame convention (t1/t2) [DIFFERENT-REQUIRES-MAPPING]

**What SeisSol does** — `BaseDRInitializer.cpp:198-202`:

```cpp
    VrtxCoords strike{};
    VrtxCoords dip{};
    misc::computeStrikeAndDipVectors(fault.normal, strike, dip);
    seissol::transformations::symmetricTensor2RotationMatrix(
        fault.normal, strike, dip, faultTractionToCartesianMatrixView, 0, 0);
```

and `:249-250`:

```cpp
    seissol::transformations::inverseSymmetricTensor2RotationMatrix(
        fault.normal, fault.tangent1, fault.tangent2, cartesianToFaultCSMatrixView, 0, 0);
```

so SeisSol's `tangent1 = strike`, `tangent2 = dip`. Friction-solver index convention (`Misc.h:164-177`): `T1 = XY = 3`, `T2 = XZ = 5`, so `traction1 = XY = strike-parallel`, `traction2 = XZ = dip-parallel`.

**What MFEM does** — `dynamic/tpv104_setup.hpp:28-35` documents the swap explicitly:

```cpp
// Coordinate frame ([C1]): pre-stress fields are written in the BP5 /
// Tandem canonical frame (`tangent1 = dip, tangent2 = strike`), mirroring
// `dynamic/tpv102_setup.hpp:62-78` under R-801 Option A.  TPV104 is pure
// strike-slip, so `tau2_0 = tau_ini` and `V2 = V_ini`; `tau1_0` and `V1`
// remain zero.  T_TPV104_SETUP_4 pins this layout; T_TPV104_SETUP_3
// independently verifies the `FaultBasis::ComputeOrientedFrame` output
// for the TPV104 reference normal + up vectors so a future BP5 edit
// that flips the convention would trip this gate immediately ([C2]).
```

and implementation at `dynamic/tpv104_setup.hpp:108-114`:

```cpp
      // Background pre-stress in the BP5 / Tandem canonical frame
      // (tangent1 = dip, tangent2 = strike).  TPV104 is pure strike-slip
      // so `tau1_0 = 0` (no dip) and `tau2_0 = tau_ini` (along-strike).
      // `sigma_n0 > 0` = compression (geology convention — §3.11).
      d.sigma_n0 = TPV104Params::sigma_n;
      d.tau1_0   = 0.0;
      d.tau2_0   = TPV104Params::tau_ini;
```

The station writer, `dynamic/tpv104_setup.hpp:369-379`, explicitly swaps columns so output matches the SCEC column order:

```cpp
         // SCEC column order:  t, h-slip, h-slip-rate, h-shear-stress,
         //                     v-slip, v-slip-rate, v-shear-stress,
         //                     n-stress, psi.
         // h = strike = component 2, v = dip = component 1.
         files_[s] << std::scientific << std::setprecision(10)
                   << t << " "
                   << d.slip2 << " "
                   << d.V2 << " "
                   << d.tau2_corr << " "
                   << d.slip1 << " "
                   << d.V1 << " "
                   << d.tau1_corr << " "
                   << d.sigma_n_corr << " "
                   << d.psi << "\n";
```

**Implication**:

* The tangent-frame labels are opposite (MFEM t1=dip, SeisSol t1=strike), but the **station file column order already compensates**: MFEM writes `slip2 / V2 / tau2_corr` (its strike component) into columns 2/3/4, matching SeisSol's "horizontal-slip" columns (strike on a vertical strike-slip fault).
* Phase 3 probe-diff tool `tpv104/scripts/tpv104_column_map.py` (filename in the file listing confirms this exists) handles the inner-pipeline swap.
* Probe 1 (trial traction) dump emits `(sigma_n_trial, tau1_trial, tau2_trial)` in MFEM's convention. The probe-diff must map MFEM `tau1 ↔ SeisSol traction2` and MFEM `tau2 ↔ SeisSol traction1` at every probe stage.

**Action**: Verify `tpv104_column_map.py` applies the swap per probe stage, not just at station output. **No algorithm change.**

---

### I-12 Normal-stress sign convention [DIFFERENT-REQUIRES-MAPPING]

**What SeisSol does** — `RateAndState.h:345-358` clips normalStress to ≤ 0:

```cpp
      normalStress[pointIndex] = std::min(static_cast<real>(0.0),
                                          faultStresses.normalStress[timeIndex][pointIndex] +
                                              this->initialStressInFaultCS[ltsFace][0][pointIndex] +
                                              /* ... */);
```

SeisSol convention: compression is **negative** `σ_n`.

**What MFEM does** — `config/tpv104_params.hpp:85`:

```cpp
   static constexpr real_t sigma_n = 120.0e6;       ///< Normal stress [Pa] (positive compression)
```

and `dynamic/tpv104_friction_solver.hpp:134-140`:

```cpp
   if (sigma_n <= 0.0)
   {
      // Fault in tension — no friction; V = tau / eta_s.
      if (iterations)    { *iterations = 0; }
      if (has_converged) { *has_converged = true; }
      return (eta_s > 0.0) ? (tau_abs / eta_s) : 0.0;
   }
   const real_t inv_eta_s = 1.0 / eta_s;
   const real_t sigma_n_abs = std::abs(sigma_n);
```

Positive-compression convention. The solver uses `|σ_n|` internally — algorithmic output is identical.

**Implication**: Raw probe dumps of `σ_n_trial` / `σ_n_total` / `σ_n_corr` will differ by a sign flip between codes. Phase 3 probe-diff must flip MFEM σ_n (or `abs()` both) before residual computation.

**Action**: `tpv104_column_map.py` applies `sign_fix = -1.0` to MFEM σ_n channels. **No algorithm change.**

---

### I-13 Initial ψ inversion [CONSISTENT]

**What SeisSol does** — `RateAndStateInitializer.cpp:141-163` for FVW:

```cpp
RateAndStateInitializer::StateAndFriction
    RateAndStateFastVelocityInitializer::computeInitialStateAndFriction(real traction1,
                                                                        real traction2,
                                                                        real pressure,
                                                                        real rsA,
                                                                        real /*rsB*/,
                                                                        real /*rsSl0*/,
                                                                        real rsSr0,
                                                                        real /*rsF0*/,
                                                                        real initialSlipRate) {
  StateAndFriction result{};
  const real absoluteTraction = misc::magnitude(traction1, traction2);
  const real tmp = std::abs(absoluteTraction / (rsA * pressure));
  result.stateVariable =
      rsA * std::log(2.0 * rsSr0 / initialSlipRate * (std::exp(tmp) - std::exp(-tmp)) / 2.0);
  // ...
}
```

i.e. `ψ_ini = a · log((2V₀/V_ini) · sinh(|τ|/(a·|σ_n|)))` with `sinh(x) = (e^x - e^-x)/2` inlined.

**What MFEM does** — `config/tpv104_params.hpp:208-217`:

```cpp
inline real_t ComputeInitialPsiTPV104(real_t a)
{
   const real_t arg = TPV104Params::tau_ini / (TPV104Params::sigma_n * a);
   const real_t x   = 2.0 * TPV104Params::V0 / TPV104Params::V_ini;
   const real_t sign_c = (arg >= 0.0) ? 1.0 : -1.0;
   const real_t absC   = std::abs(arg);
   // log(x · sinh(c)) via expm1 — stable for large |c|.
   return a * (absC + std::log(x / 2.0 * -sign_c
                               * std::expm1(-2.0 * absC)));
}
```

uses the logsinh-stable identity `log(x·sinh(c)) = |c| + log((x/2)·-sign(c)·expm1(-2|c|))`. SeisSol's form `rsA·log(x·(exp(tmp)-exp(-tmp))/2)` evaluates `sinh` directly — both expressions are mathematically equivalent; MFEM's choice is safer for very large `c` (TPV104 nominal `c = f/a = 33.3`, safely below overflow for either form).

**Numerical anchor**: With TPV104 params (`τ = 40 MPa, σ_n = 120 MPa, a = 0.01, V₀ = 1e-6, V_ini = 1e-16`) — `ψ_ini ≈ 5.6359184e-01`, matching the SeisSol reference trace at `tpv104_seisol_x2_0_x3_7.5.txt` row 1 column 9. Enforced by T_TPV104_P_1 unit test.

**Implication**: Formula identical to SeisSol FVW init (no `η·V_ini` radiation-damping subtraction, which is correct for FVW initialiser). **No action.**

---

### I-14 `f_LV` `max(0, ...)` clamp sign — plan typo resolved in favour of reference [OPEN — documented]

**What SeisSol does** — `FastVelocityWeakeningLaw.h:55-59`:

```cpp
    const real lowVelocityFriction =
        std::max(static_cast<real>(0),
                 static_cast<real>(this->f0[faceIndex][pointIndex] -
                                   (this->b[faceIndex][pointIndex] - localA) *
                                       log(localSlipRate / this->drParameters->rsSr0)));
```

i.e. `f_LV = max(0, f₀ - (b - a)·ln(V/V₀))`. **Note the minus sign** — for `b > a` (velocity weakening) and `V < V₀` (pre-nucleation), `ln(V/V₀)` is negative, so `-(b-a)·ln(V/V₀) > 0` and the clamp rarely activates.

**What MFEM does** — `friction/slip_law_srw_psi.hpp:98-104`:

```cpp
   // (1) f_LV = max(0, f0 - (b - a) · log(V/V0))
   //     Sign convention: matches the canonical FVW reference and
   //     SCEC TPV104 (Noda & Lapusta 2013).  The plan's §4.2.1 eq (2a)
   //     has a sign typo; the byte-match directive resolves it in favour
   //     of the reference implementation.
   const real_t f_LV = std::max(static_cast<real_t>(0),
                                f0 - (b - a) * std::log(V / V0));
```

**Implication**:

* The MFEM code byte-matches SeisSol — they are **consistent**.
* But the plan document (`tpv104_debug_plan_2026-04-24.md §4.2.1 eq 2a`) contains the *opposite* sign (`f₀ + (b-a)·ln(V/V₀)`). The MFEM comment notes this as a plan typo.
* **OPEN question**: Confirm with the user that the SCEC TPV104 PDF matches SeisSol's sign (which the MFEM code assumes). If the SCEC PDF states `f_LV = f₀ + (b-a)·ln(V/V₀)`, BOTH codes would be wrong and the plan would be correct — in that case Phase 2 needs to flip the sign in `slip_law_srw_psi.hpp:104`.

**Action**: Phase 3 probe 2 on (ψ_in → ψ_out) with identical input tuples should reveal whether SeisSol's sign is correct (test-fixture byte match) vs SCEC spec. If SeisSol is the ground truth, MFEM is correct; if SCEC PDF is the ground truth, the plan typo propagated into SeisSol's reference implementation and BOTH codes need the fix.

---

### I-15 Default solver dispatch — Brent (driver) vs Newton-stable (library) vs Newton-SIMD (SeisSol) [DIFFERENT-EXPECTED on driver, CONSISTENT in library]

**R7-001(b) note**: The default under the production driver path is **Brent in log10(V) space** via `friction/dieterich_ruina.hpp::SolveSlipRatePsi`, not `SolveSlipRateNewtonStable`. This is hard-coded by `FaultFaceFlux::EvaluateADERTotal` calling `FrictionSolver::Solve(..., FrictionSolver::Method::Brent)`. The `MapSolver` function + its `FrictionSolver::Method` enum + `kTpv104AlmostZero` Newton code are live but reachable only through the library `Tpv104SubStepIterator`. `drivers/tpv104_driver.cpp:406-414`:

```cpp
   const FrictionSolver::Method method = MapSolver(friction_solver);
   (void)method;  // R7-001 option (b): not routed through the time loop.
```

Both Brent and SeisSol's Newton converge on the same root of a Lipschitz-smooth residual to within `tol = 1e-8`; Probe 4 (V_abs) is expected to agree at 1e-8 relative. The `--verify-dispatch` CLI flag emits `[dispatch] rank=* solver=brent` lines per rank, making the actual dispatch observable in Frontera logs.


**What SeisSol does** — `RateAndState.h:299-340` runs a SIMD loop of Newton iterates over `NumPaddedPoints` QPs simultaneously:

```cpp
    for (uint32_t i = 0; i < settings.maxNumberSlipRateUpdates; i++) {
#pragma omp simd
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
        muF[pointIndex] =
            static_cast<Derived*>(this)->updateMu(pointIndex, slipRateTest[pointIndex], details);
        g[pointIndex] = -this->impAndEta[ltsFace].invEtaS *
                            (std::fabs(normalStress[pointIndex]) * muF[pointIndex] -
                             absoluteShearStress[pointIndex]) -
                        slipRateTest[pointIndex];
      }
      // ... (convergence check + derivative + update)
```

Settings: `maxNumberSlipRateUpdates = 60`, `newtonTolerance = 1e-8`, floor `max(rs::almostZero(), ...)` with `rs::almostZero() = 1e-45` (`RateAndStateCommon.h:22-30`).

**What MFEM does** — `dynamic/tpv104_friction_solver.hpp:71-173`:

```cpp
inline real_t SolveSlipRateNewtonStable(
   real_t tau_abs, real_t psi, real_t sigma_n, real_t eta_s,
   real_t a, real_t V0,
   real_t V_prev,
   int max_iter = 60,
   real_t tol   = 1e-8,
   /* ... */)
{
   /* ... validation / tension handling ... */
   const real_t inv_eta_s = 1.0 / eta_s;
   const real_t sigma_n_abs = std::abs(sigma_n);

   // R-006: first guess unclamped — match reference.
   real_t V = V_prev;

   int it = 0;
   real_t g = 0.0;
   for (it = 0; it < max_iter; ++it)
   {
      const real_t mu = friction_stable::FrictionCoefficientStable(
                           V, psi, a, V0);
      g = -inv_eta_s * (sigma_n_abs * mu - tau_abs) - V;

      if (std::abs(g) < tol)
      {
         /* ... */
         return V;
      }

      const real_t dmu = friction_stable::FrictionCoefficientStableDerivV(
                            V, psi, a, V0);
      const real_t dg  = -inv_eta_s * (sigma_n_abs * dmu) - 1.0;

      const real_t step = g / dg;
      V = std::max(kTpv104AlmostZero, V - step);
   }
```

with `kTpv104AlmostZero = 1e-45` (line 47), `max_iter = 60`, `tol = 1e-8`.

**Implication**: Same residual (`g = -(σ_n·μ - τ)/η_s - V`), same Newton update, same tolerance, same iteration budget, same floor. MFEM is a **single-QP** scalar version where SeisSol SIMD-vectorises over `NumPaddedPoints`; numerically identical. Probe 4 (V_abs) 1e-10 relative is realistic.

**Note**: The MFEM dispatch at `drivers/tpv104_driver.cpp:236` defaults to `FrictionSolver::Method::NewtonRaphsonStable` — this routes through `ComputeStageState(..., method)` to call `SolveNRStable` (inside `FrictionSolver::Solve`). The legacy Brent path (`Method::Brent`) and Tandem-style paths remain available but **are not used by default** for TPV104.

**Action**: None. **Consistent**.

---

### I-16 State-variable resampling (postHook) [OPEN]

**What SeisSol does** — `FastVelocityWeakeningLaw.h:144-164` projects `Δψ` onto a lower-order polynomial basis via the `resampleParameter` kernel, called in `postHook` after every macro-step.

**What MFEM does** — **no equivalent**. `ψ` is updated at the DOF points directly by `UpdateStateAnalyticSlipLawSRW` (`tpv104_substep_iterator.cpp:413-419`) and stored in `d.psi` with no post-hoc projection.

**Implication**:

* For smooth ψ fields (TPV104 initial condition + early nucleation), the resample is a near no-op. MFEM runs TPV104 at `p = 1` by default, so there are no high-order modes to suppress — the resample would be a no-op even if implemented.
* For sharp ψ fronts at higher polynomial order, SeisSol's resample suppresses high-order oscillations by a few percent per macro-step.
* **OPEN**: If Phase 3 probe 2 shows persistent ψ-drift at the rupture front that's not attributable to I-05/I-06, flag this as the likely cause and consider implementing `resampleParameter` in MFEM.

**Action**: Phase 3 probe 2 monitors. **No Phase 2 change.**

---

### I-17 ψ integration cadence — macro-step vs per-sub-step [DIFFERENT-EXPECTED on driver, CONSISTENT in library]

**What SeisSol does** — ψ is updated inside the per-sub-step friction loop at `BaseFrictionLaw.h:116-121` via `updateFrictionAndSlip` → `updateStateVariableIterative` (`RateAndState.h:151-195`) with `numberStateVariableUpdates = 2` Kaneko outer iterations. Each outer iteration runs `updateStateVariable` (analytic exponential step) followed by Newton on V using the averaged prior V guess.

**What MFEM R7-001(b) driver does** — `drivers/tpv104_driver.cpp:756-771`:

```cpp
      for (int i = 0; i < num_fault_total; ++i)
      {
         dof_data[i].psi = UpdateStateAnalyticSlipLawSRW(
            psi_n[i],
            dof_data[i].slip_rate,
            dof_data[i].Dc,
            dt_step,
            V_w[i],
            dof_data[i].a,
            TPV104Params::b,
            TPV104Params::V0,
            TPV104Params::f0,
            TPV104Params::f_w);
         dof_data[i].slip1 += dof_data[i].V1 * dt_step;
         dof_data[i].slip2 += dof_data[i].V2 * dt_step;
      }
```

ONE analytic step per macro-step with `V = dof_data[i].slip_rate` (the ADER-averaged V over the whole `dt_step`, as written by `EvaluateTotal`'s `WriteBackState`).

**Implication** (explicitly disclosed by the driver at `tpv104_driver.cpp:746-755`):

* End-of-macro-step ψ agrees with SeisSol at equilibrium because `UpdateStateAnalyticSlipLawSRW` is exact for constant V.
* During rapidly-changing V (nucleation ramp-up), ψ deviates from the per-sub-step result by O(dt_macro²).
* Phase 3 probe 2 (ψ in→out) will quantify the gap. Probe 2 must accept the O(dt_macro²) drift in the R7-001(b) comparison band — tighter bounds only apply if Phase 4 wires the iterator in.

### I-17b (library-only) `Σ deltaT == dt_macro` defensive check

`dynamic/tpv104_substep_iterator.cpp:245-261` enforces that configured sub-step sizes sum to `dt_macro`:

```cpp
   const real_t dtsum = std::accumulate(deltaT_.begin(), deltaT_.end(),
                                        static_cast<real_t>(0));
   // ...
   if (rel > sum_tol)
   {
      throw std::runtime_error(
         "Tpv104SubStepIterator::Advance: Σ deltaT[o] must equal "
         "dt_macro within " + std::to_string(sum_tol)
         + " rel; got Σ = " + std::to_string(dtsum)
         + ", dt_macro = " + std::to_string(dt_macro)
         + ", O = " + std::to_string(O_size));
   }
```

Defensive programming; SeisSol's cadence is internal and fixed by `misc::TimeSteps`. CONSISTENT behaviour; library code only.

---

### I-18 Parameter values [CONSISTENT]

**What SeisSol does**: Parameters read from a per-run input file (`parameter.par` + fault easi file). For canonical SCEC TPV104 per `docs/tpv104.rst` + SCEC validation PDF.

**What MFEM does** — `config/tpv104_params.hpp:41-98`:

```cpp
struct TPV104Params
{
   static constexpr real_t rho = 2670.0;
   static constexpr real_t cs  = 3464.0;
   static constexpr real_t cp  = 6000.0;
   // ...
   static constexpr real_t f0     = 0.6;
   static constexpr real_t V0     = 1.0e-6;
   static constexpr real_t b      = 0.014;
   static constexpr real_t L      = 0.4;
   static constexpr real_t f_w    = 0.1;
   static constexpr real_t a_in   = 0.01;
   static constexpr real_t da     = 0.01;
   static constexpr real_t a_out  = a_in + da;
   static constexpr real_t V_w_in  = 0.1;
   static constexpr real_t dV_w    = 0.9;
   static constexpr real_t V_w_out = V_w_in + dV_w;
   // ...
   static constexpr real_t sigma_n = 120.0e6;
   static constexpr real_t tau_ini = 40.0e6;
   static constexpr real_t V_ini   = 1.0e-16;
   // ...
   static constexpr real_t nuc_dtau = 45.0e6;
   static constexpr real_t nuc_T    = 1.0;
   // ...
   static constexpr real_t t_final  = 12.0;
};
```

with compile-time consistency checks at lines 102-110:

```cpp
static_assert(TPV104Params::a_out   == TPV104Params::a_in   + TPV104Params::da, /*...*/);
static_assert(TPV104Params::V_w_out == TPV104Params::V_w_in + TPV104Params::dV_w, /*...*/);
static_assert(TPV104Params::lambda == TPV104Params::rho * TPV104Params::cp * TPV104Params::cp
              - 2.0 * TPV104Params::mu, /*...*/);
static_assert(TPV104Params::eta_s  == TPV104Params::rho * TPV104Params::cs / 2.0, /*...*/);
```

Every value matches SCEC spec. **No action.**

---

### I-19 Mesh / ADER order / station coordinates [DIFFERENT-EXPECTED]

**What SeisSol does** per `tpv104_seisol_x2_0_x3_7.5.txt:1-11`:

```
# code=SeisSol (ADER-DG) without plasticity
# order of approximation in space and time= O5
# time_step= 4.735792e-03
# num_time_steps= 2534
```

TPV5 mesh, O = 5, p = 4.

**What MFEM does** — driver CLI defaults configurable, with TPV104-specific meshes at `tpv104/mesh/tpv104_{200m,500m,1000m}.geo` (Gmsh source). Station coordinates hard-coded in `config/tpv104_params.hpp:235-245`:

```cpp
inline constexpr StationTPV104 kStationsTPV104[9] = {
   {    0.0,  3.0e3,  "x2_0_x3_3"     },
   {    0.0,  7.5e3,  "x2_0_x3_7.5"   },
   {    0.0, 12.0e3,  "x2_0_x3_12"    },
   {  9.0e3,  7.5e3,  "x2_9_x3_7.5"   },
   { 12.0e3,  3.0e3,  "x2_12_x3_3"    },
   { 12.0e3, 12.0e3,  "x2_12_x3_12"   },
   { -9.0e3,  7.5e3,  "x2_-9_x3_7.5"  },
   {-12.0e3,  3.0e3,  "x2_-12_x3_3"   },
   {-12.0e3, 12.0e3,  "x2_-12_x3_12"  },
};
```

matching SeisSol's nine canonical stations.

**Implication**: `dt` differs by ~2 orders of magnitude. Probe-diff tool must interpolate both traces onto a common time grid before residual computation. Station coordinates are identical on both codes — station lookup uses nearest-QP on each side.

**Action**: `tpv104/scripts/tpv104_column_map.py` already performs time-interpolation (per repository content).

---

### I-20 CLI flags parsed but not routed [DIFFERENT-EXPECTED]

**What SeisSol does** — N/A (SeisSol has no equivalent CLI toggle; friction-law class is picked by FL-number in the parameter file).

**What MFEM does** — `drivers/tpv104_driver.cpp:314-325` accepts `--friction-solver`, `--fault-iterator`, `--fric-law` but they do NOT reach the dispatched code path under R7-001(b):

```cpp
   // R7-001/R7-003/R7-006: these flags are accepted so smoke tests and
   // sbatch scripts from the plan can pass CLI values, but on the current
   // driver path (one-shot wave.AdvanceADER) they have NO effect on the
   // dispatched solver / iterator / friction law.  Runtime is:
   //   - friction solver: Brent (hard-coded via EvaluateADERTotal)
   //   - fault iterator : one-shot (Tpv104SubStepIterator not wired
   //                      because doing so requires editing
   //                      wave_operator.inl, on the extreme-care list)
   //   - friction law   : slip-SRW via per-QP
   //                      UpdateStateAnalyticSlipLawSRW (per-macro-step)
```

The `MapSolver` call is still invoked and validated — unknown values abort eagerly via `MFEM_ABORT` — so the CLI is syntactically live, just not routed. The **`--verify-dispatch` flag** (line 302) emits one line per rank describing the actually-dispatched solver/iterator/friction law so Frontera logs can be grepped to confirm.

**Implication**: A user or sbatch script passing `--friction-solver newton-stable` silently runs Brent instead. The banner + the `--verify-dispatch [dispatch]` line are the only protection against this being mistaken for a silent mismatch. Phase 3 probe-diff must read the `[dispatch]` lines from the MFEM log and propagate them into the diff report.

**Action**: Phase 3 probe-diff reads `[dispatch]` lines; any sbatch script that passes a non-Brent solver while `[dispatch]` reports `solver=brent` must warn in the report. When `Tpv104SubStepIterator` is eventually wired into the driver (requires editing `wave_operator.inl` — currently [C2] no-touch), `MapSolver` + the `method` value become live and this entry becomes CONSISTENT.

---

## 4. Summary

### 4.1 Production R7-001(b) driver path (what actually runs on Frontera today)

| Category | Count | Entries |
|----------|-------|---------|
| **CONSISTENT** | 8 | I-01, I-03, I-07, I-08, I-10, I-13, I-18 (+ I-14 assuming SeisSol is ground truth) |
| **CONSISTENT-WITH-DEVIATION** | 1 | I-05 (per-sub-step cadence masking) |
| **DIFFERENT-EXPECTED** | 7 | I-02 (stable-μ unused), I-04 (macro-step nucleation), I-06 (Brent vs Newton-SIMD), I-09 (single-shot imposed state), I-15 (Brent dispatch), I-17 (macro-step ψ), I-19 (mesh/order/dt), I-20 (CLI not routed) |
| **DIFFERENT-REQUIRES-MAPPING** | 2 | I-11 (t1/t2 swap), I-12 (σ_n sign) |
| **OPEN** | 2 | I-14 (plan typo vs SCEC PDF), I-16 (resample — monitor probe 2) |

### 4.2 If/when `Tpv104SubStepIterator` is wired in (future Phase 2+)

| Category | Count | Entries |
|----------|-------|---------|
| **CONSISTENT** | 13 | I-01, I-02, I-03, I-04, I-07, I-08, I-09, I-10, I-13, I-15, I-17, I-18 (+ I-14) |
| **CONSISTENT-WITH-DEVIATION** | 1 | I-05 |
| **DIFFERENT-EXPECTED** | 2 | I-06, I-19 |
| **DIFFERENT-REQUIRES-MAPPING** | 2 | I-11, I-12 |
| **OPEN** | 2 | I-14, I-16 |
| **N/A** | 1 | I-20 |

**Headline**: The MFEM TPV104 **library** code (`SlipLawSRWPsi`, `FrictionCoefficientStable`, `Tpv104SubStepIterator`, `SolveSlipRateNewtonStable`, `ApplyNucleationIncremental_TPV104`) is algorithmically consistent with SeisSol FL=103 on 13/20 items. However, the **production driver** dispatches only to Brent + `EvaluateADERTotal` + macro-step ψ + macro-step nucleation — a deliberately conservative R7-001 option (b) because the per-sub-step wiring requires editing `dynamic/wave_operator.inl` (extreme-care [C2] no-touch list). Under R7-001(b) the consistency drops to 8 CONSISTENT + 7 DIFFERENT-EXPECTED items.

Phase 3 probe-diff MUST:
1. Read `[dispatch]` lines from MFEM logs and propagate to the diff report (I-20).
2. Accept O(dt_macro²) drift on probe 2 (ψ) because MFEM updates ψ once per macro-step, not per sub-step (I-17).
3. Accept O(dt_macro) drift on probe 1 (nucleation amplitude snapshot) because MFEM injects once per macro-step (I-04).
4. Accept O(dt_macro²) drift on probe 5 (imposed state) because MFEM runs `EvaluateADERTotal` once on Q̄ (I-09).
5. Mask/coarsen the trial-traction probe (I-05).
6. Swap MFEM tau1/tau2 and V1/V2 against SeisSol traction1/traction2 and slipRate1/slipRate2 (I-11).
7. Flip sign on MFEM σ_n channels (I-12).
8. Interpolate both codes' traces onto a common time grid (I-19).
9. Monitor probe 2 drift for evidence of the resample gap (I-16).
10. Confirm sign of `f_LV` matches SCEC PDF (I-14) — if SCEC PDF is ground truth and differs, BOTH codes need the flip.

---

## 5. References

**SeisSol source** (`/Users/chunhuizhao/projects/SeisSol/src/`):

* `DynamicRupture/FrictionLaws/CpuImpl/FastVelocityWeakeningLaw.h`
* `DynamicRupture/FrictionLaws/CpuImpl/RateAndState.h`
* `DynamicRupture/FrictionLaws/CpuImpl/BaseFrictionLaw.h`
* `DynamicRupture/FrictionLaws/FrictionSolverCommon.h`
* `DynamicRupture/FrictionLaws/RateAndStateCommon.h`
* `DynamicRupture/Initializer/BaseDRInitializer.cpp`
* `DynamicRupture/Initializer/RateAndStateInitializer.cpp`
* `Numerical/GaussianNucleationFunction.h`
* `DynamicRupture/Misc.h`
* `docs/tpv104.rst`

**MFEM TPV104 source** (`/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/`):

* `config/tpv104_params.hpp` (parameters + spatial helpers + stations)
* `friction/slip_law_srw_psi.hpp` (FVW state evolution + analytic integrator)
* `friction/friction_coeff_stable.hpp` (stable `μ`)
* `dynamic/tpv104_setup.hpp` (fault-DOF init + V_w side-channel + station writers)
* `dynamic/tpv104_nucleation.hpp` (cumulative accumulator)
* `dynamic/tpv104_friction_solver.hpp` (Newton-Raphson solver)
* `dynamic/tpv104_substep_iterator.{hpp,cpp}` (sub-step pipeline)
* `drivers/tpv104_driver.cpp` (driver main)
* `tpv104/scripts/tpv104_column_map.py` (probe-diff + station mapping)

**Benchmark data**:

* `tpv104/benchmark_data/seisol/tpv104_seisol_x2_*_x3_*.txt` (9 station traces)

---

**Report generated**: 2026-04-24 (revised to reflect commits 4d2f9e8 + ffffc35).
**Audit method**: Verbatim code inspection of every cited file+line range on both sides.
**Commits audited**:
* `4d2f9e8` — "TPV104 miniapp: R7/R8 honest-dispatch driver + full Phase-2 scope" — adds the TPV104 driver, library, sbatch scripts, probe-diff tooling, 264 unit tests; explicitly discloses R7-001 option (b) dispatch (banner + `--verify-dispatch` + docstrings).
* `ffffc35` — "TPV102 pepper + BP5 in-progress scope: debug dumps, build support" — Makefile + `FaultFaceFlux::EvaluateTotalFaceAveraged` infrastructure; the TPV104 miniapp builds against this commit's Makefile changes.

**Next step**: Phase 3 execution on Frontera (requires user approval per `feedback_frontera_approval.md`).
* **R7-001(b) path**: 8/20 CONSISTENT, 7/20 DIFFERENT-EXPECTED items accepted as O(dt²) or architectural baseline; probe thresholds must be widened accordingly (see §4.1).
* **Phase 2+ scope** (future): wire `Tpv104SubStepIterator` into the time loop by editing `dynamic/wave_operator.inl` (requires lifting the extreme-care [C2] no-touch rule via a separate approval); this would collapse I-02, I-04, I-09, I-15, I-17, I-20 to CONSISTENT or N/A.
