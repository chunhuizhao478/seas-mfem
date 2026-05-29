# Implementation Plan: RK4/RK45 + Linear Slip-Weakening (LSW) + Scalar Mixed Flux

**Date:** 2026-05-29
**Branch / worktree:** `worktree-feature-rk45-slipweakening`
**Scope:** SEAS-MFEM dynamic-rupture miniapp (`miniapps/seas/`), spatial driver
`drivers/spatial_dyn_driver.cpp`.
**Parallels:** `document/mixed_flux_dev/BUILD_mixed_flux_rk_dynamic_rupture_2026-05-28.md`
(the *rate-and-state* RK build doc — this is its LSW sibling).
**Immediate target:** SCEC **TPV205** (homogeneous, scalar Godunov, linear
slip-weakening) at **p1**, run with **`--mixed-flux adjacent --time-integrator
rk45`**.

---

## Overview

Today the explicit Runge–Kutta integrators (`rk4`/`rk45`) run **rate-and-state
friction only**; linear slip-weakening (LSW) is hard-blocked. The reason is a
single missing fault kernel: the RK time loop drives `WaveOperator::Mult`, whose
fault face calls the **instantaneous** `FaultFaceFlux::Evaluate` (rate-and-state
Riemann solve); there is no instantaneous LSW solve on the `Mult` path — only the
ADER time-integrated `EvaluateADER_LSW` (I-form), welded to the ADER substep
iterator. This plan adds (1) an instantaneous `FaultFaceFlux::EvaluateLSW`,
(2) a `Mult`-path dispatch on the existing `FaultFrictionLaw` tag, (3) an
LSW-specialized coupled-RK macro-step that integrates `(Q, slip)` (LSW has **no**
ψ state variable — its friction coefficient depends on accumulated slip), and
(4) driver wiring that relaxes the LSW-blocking guard. Everything is additive
behind the existing `--time-integrator` / `--mixed-flux` CLI selectors; the
byte-exact ADER regression contract (`make test`) is preserved.

### Physics (what is being discretized)

Linear slip-weakening friction (SCEC TPV5 §7-11, `dynamic/tpv205_friction.hpp`):

```
μ(δ)        = μ_s − (μ_s − μ_d) · min(δ/d_c, 1),    δ = √(slip1² + slip2²)
τ_strength  = μ(δ) · max(σ_n_total, σ_n_floor) + C0      (C0 = cohesion, 0 for TPV205)
V_abs       = max(0, (|τ_total| − τ_strength) / η_s)
(V1, V2)    = V_abs · (τ1_total, τ2_total) / |τ_total|
τ*_corr     = τ*_trial − η_s · V*                         (trial-scale; pre-stress re-added on writeback)
```

The state variable for LSW is the **accumulated slip** `δ`, integrated by the
method of lines: `d(slip*)/dt = V*`. This is the structural contrast with
rate-and-state, where the stepped state is ψ and slip is a pure output
accumulator. **In LSW the slip is BOTH the stage input (μ depends on δ) AND the
integrated output** — so the RK stepper must write the stage-local slip into
`DOFData` *before* each `Mult`, exactly as the RS stepper writes the stage-local
ψ (`rk_time_stepper.hpp:230-240`).

Why RK (not ADER) for mixed flux: central flux is non-dissipative → fault-adjacent
modes sit on the imaginary axis; ADER-O2's stability function `|R(iy)|²=1+y⁴/4>1`
amplifies them, and at p1 the ADER order is locked to O2 (`O ≤ p+1`). An explicit
RK of order ≥ 3 has imaginary-axis coverage and is stable. (Full derivation: the
companion BUILD doc §1.)

---

## Constraints

- **Byte-exact ADER contract.** The RS-ADER and LSW-ADER paths
  (`Evaluate`, `EvaluateADER`, `EvaluateADER_LSW`, `AdvanceADERWithSubStep_Spatial`,
  `Tpv205SubStepIterator`) must not change behaviour. `make test` stays green;
  the TPV102/104/205-ADER station outputs stay bit-identical. New code is reached
  only via `time_integrator ∈ {rk4, rk45}` **and** `law = slip_weakening`.
- **Scalar interior flux only.** The RK path is a scalar-Godunov feature and is
  mutually exclusive with the matrix/bimaterial path (`InteriorFlux::Matrix`,
  `BimaterialWaveOperator`). The existing guard at
  `spatial_dyn_driver.cpp:778-783` (`is_rk ⇒ interior_flux == scalar`) is **kept**.
  TPV205 is homogeneous (scalar), so this is satisfied.
- **Homogeneous material.** `EvaluateLSW` inherits `EvaluateADER_LSW`'s
  `homog_ok` guard (per-side impedances equal within 1e-12 relative); a bimaterial
  LSW face aborts loudly. TPV205 is homogeneous.
- **Slip-stateless fault kernel.** `EvaluateLSW` may READ `data.slip1/slip2` (to
  form δ) but must NOT WRITE them — slip is integrated by the stepper. This
  mirrors (a) `Evaluate` being ψ-stateless (`fault_face_flux.cpp:432`
  `MFEM_ASSERT(data.psi == psi_at_entry)`) and (b) `EvaluateADER_LSW`'s R-001
  invariant ("`data.slip{1,2}` is NOT touched by this call",
  `fault_face_flux.cpp:837-846`). A new NDEBUG slip-invariance assert enforces it.
- **Sign / frame conventions (CLAUDE.md).** Component 1 = dip, component 2 =
  strike (BP5/Tandem `FaultBasis`). Slip-rate direction parallel to total
  tangential traction. `DOFData.{V1,slip1,tau1_corr}` = dip; `{V2,slip2,tau2_corr}`
  = strike. `SolveLSW_TPV205` already encodes these; do not re-derive.
- **Reference check.** LSW is a SCEC TPV5/TPV205 law (not a Tandem QD law); the
  authoritative physics oracle is `drivers/tpv205_driver.cpp` +
  `config/tpv205_params.hpp` + `dynamic/tpv205_friction.hpp`, already validated
  against SCEC. Do not invent new physics — `EvaluateLSW` REUSES
  `LSWFrictionCoefficient_TPV205` + `SolveLSW_TPV205` verbatim.
- **Forced rupture out of scope.** TPV205 disables forced rupture
  (`dummy_T_forced = 1e9`, driver :1618). The driver maps `is_lsw → LSW` (never
  `LSW_ForcedRupture`, :1162). The `Mult`-path dispatch must therefore handle
  `LSW` and `RateAndState`; `LSW_ForcedRupture` on the `Mult`/RK path **aborts**
  (documented follow-up — TPV26/27 forced-rupture-on-RK is a separate task).

### Source-of-truth call sites (current HEAD; verify before editing)

| What | File:line | Note |
|---|---|---|
| LSW-blocking guard | `drivers/spatial_dyn_driver.cpp:771-777` | relax in Phase 4 |
| scalar-flux RK guard | `drivers/spatial_dyn_driver.cpp:778-783` | **keep** |
| `is_rk` definition | `drivers/spatial_dyn_driver.cpp:769-770` | `≠ ADER` |
| `is_lsw` definition | `drivers/spatial_dyn_driver.cpp:754-755` | |
| RK rate_state-required guard | `drivers/spatial_dyn_driver.cpp:2571-2573` | gate to `!is_lsw` |
| tableau select | `drivers/spatial_dyn_driver.cpp:2568-2584` | RK4/RK45 |
| time-loop advance dispatch | `drivers/spatial_dyn_driver.cpp:2750-2766` | add LSW branch |
| RK stepper call | `drivers/spatial_dyn_driver.cpp:2757-2758` | |
| `SetFaultFrictionLaw(is_lsw?…)` | `drivers/spatial_dyn_driver.cpp:1162-1163` | already correct |
| `SetMixedFluxMode` | `drivers/spatial_dyn_driver.cpp:1176` | already correct |
| `SetCflRkAware(true)` | `drivers/spatial_dyn_driver.cpp:1833` | law-agnostic, reused |
| Mult interior-fault solve | `dynamic/wave_operator.inl:2617` | add LSW dispatch |
| Mult shared-fault solve | `dynamic/wave_operator.inl:3260` | add LSW dispatch |
| ADER interior-fault dispatch (pattern to mirror) | `dynamic/wave_operator.inl:3822-3846` | |
| ADER shared-fault dispatch (pattern to mirror) | `dynamic/wave_operator.inl:4796-4817` | |
| `FaultFrictionLaw` enum | `dynamic/wave_operator.hpp:82-87` | RS=0, LSW=1, LSW_FR=2 |
| `fault_friction_law_` member | `dynamic/wave_operator.hpp:762` | |
| `EvaluateADER_LSW` (derive from) | `dynamic/fault_face_flux.cpp:746-863` | |
| `Evaluate` (ψ-stateless pattern) | `dynamic/fault_face_flux.cpp:328-439` | |
| `WriteBackState` / `BuildImposedState` | `dynamic/fault_face_flux.cpp:271-323` | reuse |
| `LSWFrictionCoefficient_TPV205` / `SolveLSW_TPV205` | `dynamic/tpv205_friction.hpp:54,104` | reuse |
| coupled-RK macro-step (RS) | `dynamic/rk_time_stepper.hpp:182-317` | LSW sibling |
| CFL RK-aware factors | `dynamic/wave_operator.inl:5550-5599` | reused as-is |

---

## Phase 1 — Instantaneous `FaultFaceFlux::EvaluateLSW`

### Goal
A new instantaneous LSW fault Riemann solve exists that is the `dt→0` limit of
`EvaluateADER_LSW` and is bit-exact equivalent to it when fed `I± = Q± · dt`.

### Files to Modify
- `dynamic/fault_face_flux.hpp` — declare `EvaluateLSW`.
- `dynamic/fault_face_flux.cpp` — implement `EvaluateLSW`.

### Detailed Requirements

1. **Signature** (mirror `Evaluate`, drop `method` — LSW is closed-form, no
   Brent):
   ```cpp
   /// Instantaneous LSW counterpart to Evaluate (the dt→0 limit of
   /// EvaluateADER_LSW: no I-form, no dt).  Reads ONLY the LSW-native fields
   /// data.lsw_mu_s / lsw_mu_d / lsw_d_c / lsw_cohesion; data.a / psi / Dc are
   /// not consumed.  SLIP-STATELESS: reads data.slip1/slip2 to form δ, never
   /// writes them (the RK stepper integrates slip).  Writes
   /// data.{slip_rate,V1,V2,tau1_corr,tau2_corr,sigma_n_corr} (TOTAL traction
   /// per WriteBackState).
   void EvaluateLSW(DOFData &data,
                    const real_t *Q_plus, const real_t *Q_minus,
                    real_t *Q_imp_plus, real_t *Q_imp_minus) const;
   ```

2. **Body** — `EvaluateADER_LSW` (`fault_face_flux.cpp:746-863`) with the I/dt
   wrap removed (no Step 0 `Q̄ = I/dt`, no Step 8 `I_imp = dt·Q_imp`; operate
   directly on `Q_plus`/`Q_minus`):
   - Keep the all-zero-LSW-fields misuse guard (`:761-769`).
   - Keep the homogeneous-material `homog_ok` guard (`:771-784`).
   - `ComputeTrialTraction(data, Q_plus, Q_minus, s.sigma_n_trial, s.tau1_trial, s.tau2_trial)`.
   - total traction: `s.sigma_n_total = data.sigma_n0 + data.sigma_n_nuc + s.sigma_n_trial`
     (and τ1/τ2 likewise); `s.Theta = √(τ1_total² + τ2_total²)`.
   - `delta = √(data.slip1² + data.slip2²)`;
     `mu_eff = LSWFrictionCoefficient_TPV205(delta, data.lsw_mu_s, data.lsw_mu_d, data.lsw_d_c)`.
   - `SolveLSW_TPV205(s.tau1_trial, s.tau2_trial, s.tau1_total, s.tau2_total,
     s.sigma_n_total, data.eta_s, mu_eff, s.V_abs, s.V1, s.V2, s.tau1_corr,
     s.tau2_corr, SigmaNStrengthFloorForLSW(), data.lsw_cohesion)`.
   - `s.sigma_n_corr = s.sigma_n_trial`.
   - **No slip update** (R-001 invariant).
   - `BuildImposedState(data, s, Q_plus, Q_minus, Q_imp_plus, Q_imp_minus)`.
   - `WriteBackState(data, s)`.

3. **Slip-invariance guard** (mirror the ψ guard in `Evaluate:340-342,428-438`):
   ```cpp
   #ifndef NDEBUG
      const real_t slip1_at_entry = data.slip1, slip2_at_entry = data.slip2;
   #endif
   ... // body
   #ifndef NDEBUG
      MFEM_ASSERT(data.slip1 == slip1_at_entry && data.slip2 == slip2_at_entry,
                  "FaultFaceFlux::EvaluateLSW mutated data.slip{1,2}; the "
                  "coupled-RK-on-slip integrator assumes this kernel is "
                  "slip-stateless.");
   #endif
   ```

### Edge Cases
- **Strength barrier** (`lsw_mu_s ≥ ½·mu_s_barrier`): handled inside
  `LSWFrictionCoefficient_TPV205` + `SolveLSW_TPV205` (V=0, τ_corr=τ_trial). No
  special-casing in `EvaluateLSW`.
- **Tensile σ_n** (σ_n_total ≤ 0): handled by `SolveLSW_TPV205` floor logic
  (`max(σ_n, sigma_n_floor)`); default floor disabled → byte-exact `max(σ_n,0)`.
- **δ = 0 at rest**: `μ(0) = μ_s`; correct (no slip yet).
- **All-zero LSW fields** (RS DOFData mis-routed): the misuse guard aborts.

### Acceptance Criteria
- [ ] For `dt ∈ {1.0, 0.5, 0.25}` (powers of two — the I/dt round-trip is exact):
      `EvaluateLSW(data, Q±)` writes `data.{V1,V2,slip_rate,tau1_corr,tau2_corr,
      sigma_n_corr}` **bit-identical** to `EvaluateADER_LSW(data, Q±·dt, dt)`, and
      `Q_imp_± == I_imp_±/dt` exactly.
- [ ] For a general `dt` (e.g. 0.1): the same fields match within relative
      tolerance 1e-12. (`EvaluateADER_LSW` reconstructs `Q̄ = (Q·dt)·(1/dt)`,
      which equals `Q` bit-exactly only for power-of-two `dt`; for a general `dt`
      the ~1-ULP round-trip gap propagates through trial→solve. The DOFData writes
      are otherwise dt-independent — same closed-form solve on `Q̄`, only `I_imp`
      rescales.)
- [ ] `data.slip1/slip2` unchanged across the call (NDEBUG assert; also checked in
      release by an explicit test compare).
- [ ] `make test` still green (no existing target calls `EvaluateLSW`; additive).

### Dependencies
- Depends on: nothing (pure addition).
- Required by: Phase 2, Phase 3.

---

## Phase 2 — Dispatch `EvaluateLSW` on the `Mult` path

### Goal
`WaveOperator::Mult` routes the fault face to `EvaluateLSW` when
`fault_friction_law_ == LSW`, and to the unchanged `Evaluate` otherwise — on both
the interior-fault and shared-fault branches.

### Files to Modify
- `dynamic/wave_operator.inl` — interior-fault dispatch (~`:2617`), shared-fault
  dispatch (~`:3260`).

### Detailed Requirements

1. **Interior-fault `Mult` dispatch** — replace the unconditional
   `fault_flux_->Evaluate(fdata, Q_plus_local, Q_minus_local, Q_imp_plus, Q_imp_minus);`
   at `:2617` with a switch on `fault_friction_law_`, mirroring the ADER
   dispatch's structure (`:3822-3846`):
   ```cpp
   if (fault_friction_law_ == FaultFrictionLaw::LSW)
   {
      fault_flux_->EvaluateLSW(fdata, Q_plus_local, Q_minus_local,
                               Q_imp_plus, Q_imp_minus);
   }
   else if (fault_friction_law_ == FaultFrictionLaw::LSW_ForcedRupture)
   {
      MFEM_ABORT("WaveOperator::Mult: FaultFrictionLaw::LSW_ForcedRupture has no "
                 "instantaneous solve on the Mult/RK path (only EvaluateADER_LSW_"
                 "ForcedRupture exists, for ADER).  Use --time-integrator ader "
                 "for forced-rupture configs.");
   }
   else
   {
      fault_flux_->Evaluate(fdata, Q_plus_local, Q_minus_local,
                            Q_imp_plus, Q_imp_minus);
   }
   ```
2. **Shared-fault `Mult` dispatch** — identical switch at the shared-fault
   `Evaluate` call (`:3260`). Use the same `Q_plus_local`/`Q_minus_local` already
   in scope there.
3. **Do NOT touch** the ADER call sites (`:3822-3846`, `:4796-4817`,
   `EvaluateADER`/`EvaluateADER_LSW`). Those remain the ADER path.
4. Default (`RateAndState = 0`) preserves the exact current RS `Mult` behaviour.

### Edge Cases
- `fault_friction_law_` defaults to `RateAndState` (`wave_operator.hpp:762`), so
  any operator not explicitly set to LSW is unchanged. (The spatial driver sets it
  at :1162; the ADER LSW path for TPV205 already relies on it.)
- The interior fault path may run on a mesh with **no** fault faces — the switch is
  inside the existing per-fault-face loop, so zero faces ⇒ zero dispatch (unchanged).

### Acceptance Criteria
- [ ] With `fault_friction_law_ == LSW`, one `Mult` on a TPV205-init `dof_data`
      writes the same `data.{V1,V2,slip_rate,tau*_corr,sigma_n_corr}` as calling
      `EvaluateLSW` directly (interior-fault, serial).
- [ ] With `fault_friction_law_ == RateAndState` (default), `Mult` output is
      bit-identical to the pre-change code on a TPV102-init `dof_data`
      (regression: `seas_test_*` that exercise the RS `Mult` path stay green).
- [ ] `LSW_ForcedRupture` + `Mult` aborts with the documented message.
- [ ] **MPI shared-fault parity (guards the `:3260` edit):** a 2-rank run with the
      fault crossing the rank seam and `fault_friction_law_ == LSW` produces
      `V/τ*_corr/σ_n_corr` on the shared QPs equal (to round-off) to the values the
      same QPs get on a single-rank (interior-path) run. This is REQUIRED for Phase
      2 sign-off because the shared-fault dispatch (`wave_operator.inl:3260`) is a
      parallel-only path that **none** of the serial tests L1–L7 reach — a
      wrong/missing LSW arm there compiles, passes `make test`, and only fails at
      np>1 (i.e. every production run). See test L8.
- [ ] `make test` green.

### Dependencies
- Depends on: Phase 1.
- Required by: Phase 3 (the RK stepper drives `Mult`), Phase 4.

---

## Phase 3 — LSW coupled-RK macro-step `AdvanceRKCoupledLSW_Spatial`

### Goal
A coupled explicit-RK macro-step exists that integrates `(Q, slip)` with a Butcher
tableau, staging the slip into `DOFData` before each `Mult` so `EvaluateLSW` sees
the stage-local δ. No ψ, no `PsiRateEvaluator`, no rate-state config.

### Files to Modify
- `dynamic/rk_time_stepper.hpp` — add the templated `AdvanceRKCoupledLSW_Spatial`
  next to `AdvanceRKCoupled_Spatial`. (Header-resident because templated on
  `MeshType`, like the RS sibling.) Also add the symmetric friction-law guard to
  the existing `AdvanceRKCoupled_Spatial` (requirement 4 below) — additive
  `MFEM_VERIFY` only, byte-exact preserved.

### Detailed Requirements

1. **Signature** — like `AdvanceRKCoupled_Spatial` (`rk_time_stepper.hpp:182-192`)
   but **without** `rs_cfg` / `rs`:
   ```cpp
   /// Advance (Q, slip) by ONE macro-step dt_step with tableau `tab`, coupling the
   /// bulk wave field and the LSW fault slip through the SAME weights.  LSW has no
   /// ψ; the stepped fault state is the accumulated slip δ (μ depends on it), so
   /// the stage-local slip is written into dof_data BEFORE each wave.Mult — the
   /// LSW analogue of AdvanceRKCoupled_Spatial's stage-local ψ write.
   ///   slip*^(i) = slip*_n + dt·Σ_{j<i} a_ij V*_k[j]   (written before Mult)
   ///   wave.Mult(Q^(i), k_i)   — runs the slip-stateless EvaluateLSW at slip^(i)
   ///   capture V1_k[i]/V2_k[i] = dof_data[m].V1/V2
   /// Final combine:
   ///   Q_new      = Q + dt·Σ_i b_i k_i
   ///   slip*_new  = slip*_n + dt·Σ_i b_i V*_k[i]
   /// Requires wave.GetFaultFrictionLaw() == FaultFrictionLaw::LSW.
   template <typename MeshType>
   void AdvanceRKCoupledLSW_Spatial(WaveOperator<MeshType>&  wave,
                                    std::vector<DOFData>&    dof_data,
                                    const Vector&            Q,
                                    real_t                   dt_step,
                                    real_t                   t_step_start,
                                    Vector&                  Q_new,
                                    const RKTableau&         tab,
                                    INucleationMethod*       nuc);
   ```

2. **Body** — structurally parallel to `AdvanceRKCoupled_Spatial:194-317`, with ψ
   replaced by slip as the staged state:
   - `MFEM_VERIFY(dt_step > 0.0, ...)`; `ValidateTableau(tab)`.
   - **Self-check (required):** `MFEM_VERIFY(wave.GetFaultFrictionLaw() ==
     FaultFrictionLaw::LSW, ...)` so a mis-wired (stepper, flag) pair fails loud
     (see requirement 4 for why this must be symmetric with the RS stepper).
   - Snapshot per-DOF slip at step start:
     `slip1_n[m] = dof_data[m].slip1; slip2_n[m] = dof_data[m].slip2;`.
   - Per-stage buffers `V1_k[s][n]`, `V2_k[s][n]` (no `psi_k`).
   - For each stage `i`:
     - `Q_stage = Q + dt·Σ_{j<i} a_ij k_j` (same as RS).
     - **Stage-local slip written BEFORE Mult:**
       `dof_data[m].slip1 = slip1_n[m] + dt·Σ_{j<i} a_ij V1_k[j][m]` (and slip2).
     - `if (nuc) nuc->ApplyAbsolute(dof_data, t_step_start + tab.c[i]*dt_step);`
       (no-op for TPV205; preserved for generality, e.g. instantaneous overstress).
     - `wave.Mult(Q_stage, k[i]);` — `EvaluateLSW` reads the staged δ.
     - Capture `V1_k[i][m] = dof_data[m].V1; V2_k[i][m] = dof_data[m].V2;` and
       update `dof_data[m].slip_rate_substep_max = max(…, dof_data[m].slip_rate)`.
   - Final combine: `Q_new = Q + dt·Σ_i b_i k_i`; then per DOF
     `dof_data[m].slip1 = slip1_n[m] + dt·Σ_i b_i V1_k[i][m]` (and slip2).
   - **Endpoint re-evaluation** (non-FSAL only, exactly as RS at `:299-316`):
     the FSAL predicate `a[s-1][j] == b[j] ∀j`. For non-FSAL (RK4), set
     `dof_data[m].slip1/slip2 = slip*_new` (already done by the combine), apply
     `nuc->ApplyAbsolute(dof_data, t_step_start + dt_step)`, run one extra
     `wave.Mult(Q_new, k_endpoint)` to refresh `dof_data.{slip_rate,V1,V2,tau*_corr,
     sigma_n_corr}` to true-endpoint values, and re-reduce `slip_rate_substep_max`.
     **Critical:** the endpoint `Mult` must run AFTER the slip combine so
     `EvaluateLSW` sees `slip_new` (= δ at t+dt), not the last stage's predictor
     slip. (Slip-stateless ⇒ this re-eval does not corrupt the integrated slip.)

3. **Slip staging is the only structural delta vs RS.** In the RS stepper, slip is
   a pure output accumulator (RS `Evaluate` does not read slip), so RS writes the
   combined slip only at the end. In LSW, slip is *also* the stage input, so it
   must be written at the start of every stage. Getting this wrong (e.g. forgetting
   the stage-local slip write, or running the endpoint `Mult` before the slip
   combine) produces a stale-δ bug: μ lags the rupture by one macro-step and the
   front speed is wrong. This is the LSW analogue of the git:8461c67 R-V92-K01
   half-step-lag bug the RS stepper's endpoint re-eval fixed.

   **Combine must overwrite, not accumulate (porting trap).** Because the stage
   loop overwrites `dof_data[m].slip1/slip2` every stage, at combine time those
   fields hold the LAST stage's staged slip — NOT `slip*_n`. The combine therefore
   writes the ABSOLUTE value `dof_data[m].slip1 = slip1_n[m] + dt·Σ_i b_i V1_k[i][m]`
   (`=`, from the step-start snapshot). Do **NOT** copy the RS sibling's
   `dof_data[m].slip1 += slip1_inc` (`rk_time_stepper.hpp:283`): in the RS stepper
   slip is unstaged so `dof_data.slip1` still equals `slip1_n` at combine time and
   `+=` is correct; in the LSW stepper `+=` double-counts the last stage's partial
   staged sum, silently corrupting the rupture front. Test L6 asserts the ABSOLUTE
   post-step slip precisely to catch this.

4. **Symmetric friction-law guard on the RS stepper (required).** After Phase 2,
   `WaveOperator::Mult` selects the fault kernel purely from `fault_friction_law_`,
   so the *stepper* and the *flag* must agree or `Mult` silently runs the wrong
   kernel: an RS stepper on an `LSW`-flagged operator runs `EvaluateLSW` (δ frozen,
   ψ integrated but ignored → wrong physics), and an LSW stepper on a
   `RateAndState`-flagged operator runs `Evaluate` on LSW DOFData (`data.psi`/`data.b`
   defaults → NaN). The LSW stepper's self-check (requirement 2) covers one
   direction; the existing `AdvanceRKCoupled_Spatial` (RS) currently has **no** such
   check (confirmed: no `GetFaultFrictionLaw` reference in
   `rk_time_stepper.{hpp,cpp}`). Add the mirror guard at the top of
   `AdvanceRKCoupled_Spatial`:
   ```cpp
   MFEM_VERIFY(wave.GetFaultFrictionLaw() == FaultFrictionLaw::RateAndState,
               "AdvanceRKCoupled_Spatial integrates psi, but the WaveOperator's "
               "fault_friction_law_ is not RateAndState; Mult would run the wrong "
               "fault kernel.  Use AdvanceRKCoupledLSW_Spatial for slip_weakening.");
   ```
   This is additive (an `MFEM_VERIFY`, no arithmetic change), so the byte-exact
   RS-RK contract is preserved; it only makes the (stepper, flag) mismatch fail
   loud. In the production driver the pair is wired from a single `is_lsw` and
   cross-checked at `spatial_dyn_driver.cpp:2518-2527`, so this guard primarily
   protects future refactors and unit tests that construct operators directly.

### Edge Cases
- `n == 0` (no local fault DOFs, e.g. a rank with no fault): all per-DOF loops are
  empty; the bulk `Q` integration still runs (RK on a pure wave field). Endpoint
  re-eval guarded by `n > 0` like the RS path.
- `nuc == nullptr`: skip nucleation (frictionless / test). TPV205 passes a no-op
  `NoNucleation` (no `[nucleation]` block), so `nuc->ApplyAbsolute` is a safe no-op.
- Barrier QPs (`lsw_mu_s ≥ ½ barrier`): `EvaluateLSW` returns V=0, so
  `V*_k[i][m] = 0` and slip never accumulates there — correct lock.

### Acceptance Criteria
- [ ] **Bulk equivalence:** on a frictionless cube (empty `dof_data`, `nuc=nullptr`),
      `AdvanceRKCoupledLSW_Spatial(RK4)` produces the same `Q_new` as the proven
      `DoRK4Step` idiom and as `AdvanceRKCoupled_Spatial` — bit-for-bit on the
      stage `k`-vectors, `Q_new` to round-off. (Reuse the T3 pattern from
      `test_rk_time_stepper.cpp:231-298`.)
- [ ] **Slip accumulation (absolute, not increment):** a single-QP LSW problem
      driven so that δ ≥ d_c (μ = μ_d constant) and `Q±` held fixed across the step
      ⇒ V is constant; assert the **absolute** post-step `dof_data.slip1` equals
      `slip1_n + V1·dt` to round-off (the RK4 quadrature of a constant is exact via
      Σb=1). The δ≥d_c regime pins a constant μ so the closed-form expected value is
      exact, and asserting the absolute slip (not an increment) fails loudly if the
      combine uses `+=` instead of `=` (the R-002 porting trap).
- [ ] **Stage-local δ is honoured:** a test where μ(δ) changes appreciably within a
      step (δ crossing `d_c`) shows the stepped slip differs from a (wrong)
      no-staging variant — i.e. the stage slip write is exercised, not dead code.
- [ ] **FSAL predicate:** DP45 skips the endpoint re-eval; RK4 performs it (assert
      via the `last_stage_is_endpoint` logic, as `test_rk_time_stepper.cpp` T6 does).
- [ ] **Friction-law guard (both directions, death test):**
      `AdvanceRKCoupledLSW_Spatial` aborts on a `RateAndState`-flagged operator, and
      `AdvanceRKCoupled_Spatial` aborts on an `LSW`-flagged operator (requirement 4).
- [ ] `make test` green.

### Dependencies
- Depends on: Phase 1, Phase 2.
- Required by: Phase 4.

---

## Phase 4 — Driver wiring + guard relaxation

### Goal
`seas_spatial_dyn_driver` runs `law=slip_weakening` with `--time-integrator
rk4|rk45 --mixed-flux adjacent` end-to-end, selecting the LSW RK stepper, keeping
all other paths byte-exact.

### Files to Modify
- `drivers/spatial_dyn_driver.cpp`.

### Detailed Requirements

1. **Relax the LSW-blocking guard** (`:771-777`). Replace the
   `MFEM_VERIFY(!is_rk || !is_lsw, …)` with a guard that allows RK+LSW for the
   supported configuration and still rejects the unimplemented forced-rupture case:
   ```cpp
   // RK + LSW is supported for the scalar-Godunov, non-forced-rupture path
   // (instantaneous EvaluateLSW + AdvanceRKCoupledLSW_Spatial).  LSW forced
   // rupture has no instantaneous solve yet (only EvaluateADER_LSW_ForcedRupture).
   // is_lsw maps to FaultFrictionLaw::LSW (never LSW_ForcedRupture) at :1162, so
   // there is no forced-rupture path here today; assert that contract.
   ```
   (No `MFEM_VERIFY` needed against forced rupture beyond the `Mult`-path abort in
   Phase 2, because `is_lsw → LSW`. If a future config introduces
   `LSW_ForcedRupture`, the `Mult` abort in Phase 2 catches it.)
2. **Keep** the scalar-flux RK guard (`:778-783`) verbatim.
3. **Gate the RK rate_state-required guard** (`:2571-2573`) under `!is_lsw`:
   ```cpp
   if (is_rk && !is_lsw)
   {
      MFEM_VERIFY(cfg.rate_state.has_value(),
                  "--time-integrator rk4|rk45 (rate_state) requires a "
                  "[friction.rate_state] block (the RK ψ coupling reads it).");
   }
   ```
   The tableau selection (`:2574-2583`) stays for both laws (LSW uses the same
   tableaus). Move the `cfg.rate_state` read out of the unconditional path.
4. **Time-loop advance dispatch** (`:2750-2766`) — three-way branch:
   ```cpp
   if (is_rk && is_lsw)
   {
      AdvanceRKCoupledLSW_Spatial(wave, dof_data, Q, dt_step, t, Q_new,
                                  rk_tab, nuc.get());
   }
   else if (is_rk)   // rate_state
   {
      AdvanceRKCoupled_Spatial(wave, dof_data, *cfg.rate_state, rs,
                               Q, dt_step, t, Q_new, rk_tab, nuc.get());
   }
   else              // ADER (both laws)
   {
      AdvanceADERWithSubStep_Spatial(wave, substep_iterator, dof_data,
                                     fault_coords, Q, dt_step,
                                     cfg.numerics.ader_order, t, Q_new, nuc_cb);
   }
   ```
5. **Already correct, do not duplicate:** `SetFaultFrictionLaw(is_lsw ? LSW : RS)`
   (:1162), `SetMixedFluxMode(ParseMixedFlux(cfg.numerics.mixed_flux))` (:1176),
   `SetCflRkAware(true)` for `is_rk` (:1833). The CFL RK-aware factors
   (`wave_operator.inl:5550-5599`) are law-agnostic (keyed only on mixed-flux mode
   + the RK-aware flag), so LSW reuses them unchanged. The single-source-of-truth
   cross-check at `:2518-2527` (wave-op law == is_lsw) still holds.
6. **Banner:** the existing startup banner already prints `law`, `time integrator`,
   `mixed flux` (`:801-823`). No new banner field strictly required; optionally add
   a one-line `[time-integrator] LSW coupled-RK (slip-stateless EvaluateLSW)` note
   on the `is_rk && is_lsw` branch for run-log clarity.
7. **`--dry-run --verify-dispatch`:** confirm the dispatch report prints
   `law=slip_weakening`, `time integrator=rk45`, `mixed flux=adjacent`,
   `interior flux=scalar` and does **not** abort.

### Edge Cases
- `is_rk && is_lsw` but `interior_flux=matrix`: rejected by the kept scalar-flux
  guard (:778-783) before reaching the dispatch.
- `is_rk && is_lsw` but a `[friction.rate_state]` block is also present in the
  TOML: the parser already forbids a rate_state block when `law=slip_weakening`
  (R-002, see :1462-1464), so `cfg.rate_state` is `nullopt`; the gated guard (step
  3) does not read it. Confirm no other unconditional `*cfg.rate_state` deref on
  the LSW path (grep `cfg.rate_state` — the PsiRate path is RS-only).
- DP45 vs RK4 selection unchanged (`:2574-2577`).

### Acceptance Criteria
- [ ] `seas_spatial_dyn_driver --config tpv205/configs/tpv205_spatial.toml
      --time-integrator rk45 --mixed-flux adjacent --dry-run --verify-dispatch`
      exits 0 with the expected scheme banner.
- [ ] The same command with `--time-integrator rk45` and a rate_state config
      (TPV102) still runs the RS RK path (no regression).
- [ ] `--time-integrator ader` on TPV205 is byte-identical to pre-change (ADER LSW
      regression).
- [ ] Builds clean: `make seas_spatial_dyn_driver` (conda `mfem-dev`).
- [ ] `make test` green.

### Dependencies
- Depends on: Phase 1, 2, 3.
- Required by: Phase 5.

---

## Phase 5 — TPV205 config, sbatch job, and verification

### Goal
A runnable, documented p1 / scalar / mixed-flux / RK45 TPV205 job exists and is
verified to be stable (the configuration that runs away under ADER at p1).

### Files to Create
- `tpv205/configs/tpv205_spatial_rk45_mixedflux.toml` — copy of
  `tpv205_spatial.toml` with `[numerics] time_integrator = "rk45"` added (or set
  `mixed_flux = "adjacent"`, `interior_flux = "scalar"` — already so) and a header
  documenting the scheme. *(Alternative: keep one TOML and pin the scheme via CLI
  flags in the sbatch, matching the TPV31 job-1 pattern. Decide per the
  "config-only vs CLI-pinned" convention — the TPV31 jobs CLI-pin the scheme and
  assert config-only knobs. Recommend the CLI-pinned sbatch + assertions pattern
  for parity with `jobs/tpv31_spatial/`.)*
  NOTE: `[numerics].time_integrator` IS a parsed TOML key
  (`spatial/code/spatial_friction.cpp:1101`, `toml_str(n, "time_integrator", "ader")`,
  same `ader|rk4|rk45` vocabulary as the CLI), so the config-only route is valid —
  CLI pinning is recommended only for the TPV31-style "scheme printed in the banner +
  asserted config-only knobs" parity, NOT because the TOML key is unsupported.
- `jobs/tpv205/tpv205_p1_rk45_mixedflux.sbatch` — modelled on
  `jobs/tpv31_spatial/tpv31_job1_p1_aderO2_noflux_50m.sbatch`: repo-root walk-up,
  module loads, pre-flight assertions on config-only knobs (`order = 1`,
  `interior_flux = "scalar"`), `--dry-run --verify-dispatch` gate, then `ibrun`
  with `--time-integrator rk45 --mixed-flux adjacent`.

### Detailed Requirements
1. The TOML/CLI scheme MUST resolve to: `law=slip_weakening`, `order=1`,
   `interior_flux=scalar`, `mixed_flux=adjacent`, `time_integrator=rk45`.
2. The sbatch pre-flight asserts the config-only knobs (`order`, `interior_flux`)
   and runs `--dry-run --verify-dispatch` before the production `ibrun`, exactly as
   the TPV31 job does (so a drifted TOML fails loudly, not silently).
3. CFL: start from the BUILD doc's RK+central-flux guidance — DRDG3D's empirical
   `CFL ≈ 0.3` for `Adjacent`. The driver's RK-aware factor for `Adjacent` is `0.6`
   (`wave_operator.inl:5576`); the config `cfl = 0.25` × DG de-rating is the
   starting point. **These are starting calibration numbers** (per the BUILD doc
   §5.4 and the existing memory note); the first dev-queue smoke must report the
   empirical stable dt before any production wall is committed.
4. Output / restart / ParaView flags mirror the TPV31 job (`--paraview-fault-vtu`,
   `--checkpoint-every`, `--output-dir`).
5. **No local full-mesh run.** Per project rule, verify locally by compile + unit
   tests + `--dry-run --verify-dispatch` only; the full TPV205 mesh runs on
   Frontera.

### Acceptance Criteria
- [ ] `--dry-run --verify-dispatch` on the new config/CLI passes locally.
- [ ] sbatch pre-flight assertions pass on the committed TOML.
- [ ] (Frontera, manual) the run is **stable** to `t_final` with bounded `V_max`
      (the ADER-p1-central-flux blow-up does not occur), and on-fault rupture
      times / final slip match the SCEC TPV205 reference within tolerance
      (`tpv205` benchmark overlays). This is the physics sign-off; it is a
      post-merge Frontera task, not a local gate.

### Dependencies
- Depends on: Phase 4.
- Required by: nothing (delivery).

---

## Testing Strategy

A new unit-test target `seas_test_lsw_rk_mixed_flux` (registered in `Makefile`
mirroring `seas_test_rk_time_stepper`, lines 224-226 / 857 / 2148-2195, with the
same object dependency list as `seas_test_rk4_conservation`:
`WAVE_OPERATOR_OBJ GODUNOV_FLUX_BIMATERIAL_OBJ PRECOMPUTED_FACE_FLUXES_OBJ
GODUNOV_FLUX_OBJ PML_LAYER_OBJ FAULT_FACE_FLUX_OBJ FRICTION_SOLVER_OBJ
RK_TIME_STEPPER_OBJ`) covering:

| ID | Phase | What |
|---|---|---|
| L1 | 1 | `EvaluateLSW(Q±)` writes DOFData bit-identical to `EvaluateADER_LSW(Q±·dt, dt)` for `dt ∈ {1.0, 0.5, 0.25}` (powers of two ⇒ exact), and within rtol 1e-12 for `dt=0.1` (the `(Q·dt)·(1/dt)` round-trip carries a ~1-ULP gap); `Q_imp == I_imp/dt`. |
| L2 | 1 | `EvaluateLSW` leaves `data.slip1/slip2` unchanged (release-mode compare, not just the NDEBUG assert). |
| L3 | 1 | Strength-barrier QP (`lsw_mu_s = mu_s_barrier`) ⇒ V=0, τ_corr=τ_trial. Tensile σ_n ⇒ default free-slide (floor 0). |
| L4 | 2 | `WaveOperator` with `SetFaultFrictionLaw(LSW)`: one `Mult` on TPV205-init `dof_data` ⇒ fault DOFData fields == direct `EvaluateLSW`. With default `RateAndState`, `Mult` == pre-change RS output (regression). `LSW_ForcedRupture` ⇒ abort. |
| L5 | 3 | Bulk RK4 equivalence: `AdvanceRKCoupledLSW_Spatial` == `DoRK4Step` == `AdvanceRKCoupled_Spatial` on a frictionless cube (k-vectors bit-exact, Q_new to round-off). |
| L6 | 3 | Absolute slip accumulation in the δ≥d_c constant-μ regime (`Q±` fixed): `dof_data.slip1 == slip1_n + V1·dt` to round-off — fails if the combine uses `+=` not `=` (R-002 trap). Plus stage-local-δ honoured (μ(δ) crossing `d_c` within a step changes the result vs a no-staging control). |
| L7 | 3 | FSAL predicate: DP45 skips endpoint re-eval, RK4 performs it (mirror T6). |
| L8 | 2 | **MPI (np=2)** fault-on-seam: the shared-fault LSW dispatch (`:3260`) gives the same fault observables (`V1,V2,tau1_corr,tau2_corr,sigma_n_corr`) as the serial interior path (`:2617`) to round-off. Run via `mpirun -np 2 seas_test_lsw_rk_shared_fault_mpi`. Registered as an MPI test target (mirror an existing `seas_test_fault_*_mpi`). This is the ONLY coverage of the parallel `:3260` edit. |

Plus the existing regression gate: full `make test` stays green (the ADER and RS-RK
paths are untouched), and the TPV205-ADER station output is bit-identical
(byte-exact contract).

**Reference solutions:** L1 uses `EvaluateADER_LSW` as the oracle (same closed
form). L5 uses the proven `DoRK4Step` idiom (git:8461c67). L6 uses the analytic
RK4-of-a-constant identity (Σb=1). The physics-level oracle (Phase 5) is the SCEC
TPV205 benchmark, via `drivers/tpv205_driver.cpp` and the `tpv205` benchmark
overlays.

---

## Risk Assessment

1. **Stale-δ / half-step-lag bug (highest risk).** If the LSW stepper fails to
   write the stage-local slip before `Mult`, or runs the endpoint `Mult` before the
   slip combine, μ(δ) lags by a macro-step and the rupture speed is wrong but the
   run does not crash. **Detection:** L6 (stage-local-δ honoured) + Phase 5 rupture
   time vs SCEC. Mirror the RS stepper's structure exactly (`rk_time_stepper.hpp`
   ψ-staging at `:230-240` and endpoint re-eval at `:299-316`).
2. **Slip double-counting on shared faces (MPI).** `EvaluateADER_LSW`'s R-001
   warns that the ADER shared-fault fallback re-invokes the LSW solve after the
   iterator, so the kernel must not accumulate slip. `EvaluateLSW` is slip-stateless
   (Phase 1) and the RK stepper integrates slip exactly once in the combine, so the
   double-count cannot occur — **but** the `Mult` shared-fault path
   (`wave_operator.inl:3260`) must be confirmed to call `EvaluateLSW` once per
   shared QP per stage (not twice). Add an MPI test (fault crossing a rank seam)
   before production, analogous to the BUILD doc §6.5 RS requirement.
3. **CFL miscalibration.** The RK-aware central-flux factors (`0.6`/`0.7`) are
   interim placeholders. An over-large dt blows up; an over-small dt wastes wall.
   **Detection:** the first dev-queue smoke reports the empirical stable dt;
   `dt→0`/`NaN` is the regression signal (CLAUDE.md). Do not commit production wall
   until the smoke calibrates dt.
4. **Accidental ADER-path change.** Editing `wave_operator.inl` near the fault
   dispatch risks touching the ADER arms. **Detection:** the byte-exact `make test`
   gate + TPV205-ADER station diff. Keep edits surgical (only the two `Mult`-path
   `Evaluate` call sites).
5. **`cfg.rate_state` deref on the LSW path.** A missed unconditional
   `*cfg.rate_state` would crash an LSW-RK run (no rate_state block). **Detection:**
   grep `cfg.rate_state` on the RK path before merging; the dispatch dry-run on
   TPV205 exercises it.
6. **`fault_friction_law_` not set / mismatched.** If a driver path drives the LSW
   RK stepper without `SetFaultFrictionLaw(LSW)`, `Mult` would run the RS
   `Evaluate` on LSW DOFData (NaN via `data.b`/`data.psi`). **Detection:** the
   stepper's optional `GetFaultFrictionLaw() == LSW` self-check (Phase 3 step 2) +
   the existing :2518-2527 cross-check.

---

## Out of scope (documented follow-ups)

- LSW **forced rupture** (TPV26/27) on the RK/Mult path — needs an instantaneous
  `EvaluateLSW_ForcedRupture` (the `f_2(t)` factor). The `Mult`-path dispatch
  aborts on `LSW_ForcedRupture` (Phase 2).
- **Bimaterial** LSW + mixed flux — excluded by the scalar-flux guard;
  bimaterial central flux + RK is the separate TPV31-job-2 task (memory:
  `tpv31-mixed-flux-rk45-needs-matrix-path-cpp`).
- **Adaptive** RK45 step-size control (embedded `bhat` error) — the tableau carries
  `bhat`, but the loop is fixed-step (CFL-bound). Adaptive control is the BUILD doc
  §5.0 Option B, deferred.
- `AllContinuous` mixed-flux mode for LSW — supported by the dispatch but needs its
  own (tighter) CFL calibration; the target uses `Adjacent`.
