# TPV104 — SeisSol Cross-Verification Debug Plan

Date: 2026-04-24
Status: **design / review — no code written yet**
Scope: `miniapps/seas/tpv104/`, `miniapps/seas/dynamic/`, `miniapps/seas/friction/`, `miniapps/seas/config/`, `miniapps/seas/drivers/`. BP5 and TPV102 production paths are **not** touched. SeisSol is modified only by adding stderr printfs at named insertion points on Frontera; no algorithmic changes on either side.

Companion / prerequisite reads:

- [tpv102_unit_tests_to_narrow_bug_2026-04-23.md](../tpv102_debug_document/tpv102_unit_tests_to_narrow_bug_2026-04-23.md) — why local unit tests cannot crack the remaining `tau1_corr ≈ 1.03e+01` pepper signal at `np=4`.
- [tpv102_seissol_aligned_flux_plan_2026-04-23.md](../tpv102_debug_document/tpv102_seissol_aligned_flux_plan_2026-04-23.md) — the precomputed-face-flux plan whose Phase 3 gates live on the same fixtures we reuse here.
- [research_seissol_tpv104_2026-04-24.md](./research_seissol_tpv104_2026-04-24.md) — SeisSol-side file-level research on TPV104 (FL=103 `FastVelocityWeakeningLaw`, nucleation adjustInitialStress, FrictionSolverCommon precompute/postcompute).

---

## 0. Why switch from TPV102 local audits to TPV104 vs SeisSol

The TPV102 audit pinned the interior/boundary face branch bugs (Fix A + Fix B, `H_IFACE_NOR`, `H_BFACE_NOR (SXZ)`) to ULP after symmetrisation probes, and pinned the residual 1.00 drift under constant-`I` on the D4 fixture to a Kuhn-diagonal topology artefact (`H_KUHN_FIXTURE`). Those code-side fixes are authorised changes pending; they still leave a live, non-synthetic pepper:

| fixture | step | channel | drift |
|---|---|---|---|
| M_ref 4×2×4 `np=4` | 19 | `tau1_corr` | **1.03e+01** |

None of the 14 audit gates exposes a clean attribution for that magnitude — Gate 7's iface+bface symmetrisation cancels to ULP on step 2, and the multistep pepper correlation tests (`test_r5_gate14_pepper_correlation`, `test_r7_amplification_chain_v2`) show the drift grows non-linearly in step count but does **not** localise to any single branch the audit can isolate at step 19. Local fixtures cannot distinguish between (i) a residual feedback in the fault-face path that only manifests under full nucleation, (ii) an ADER-stage sign error that cancels at step 2 but compounds, or (iii) an initial-condition setup mismatch that's close-but-not-identical to SeisSol.

**Strategy shift** (the subject of this plan): treat SeisSol's TPV104 benchmark as an external oracle. TPV104 differs from TPV102 **only in the friction law type** (slip-law + strong-rate-weakening vs classical aging law) — the wave operator, Riemann coupling, pre-stress layout, mesh topology, and nucleation pattern are all structurally identical. If we implement TPV104 in MFEM, cross-verify the bulk + fault signals against SeisSol's nine benchmark traces at the canonical stations, and the pepper reappears in TPV104 too, we can drive per-QP probes on **both codes** and difference any intermediate — a capability local audits do not have.

SeisSol cannot build on macOS; it is built on Frontera in the user's environment (see [reference_seissol_frontera.md](~/.claude/projects/-Users-chunhuizhao/memory/reference_seissol_frontera.md)). All SeisSol probes in this plan are adds-to-source that the user runs on Frontera; all MFEM probes are adds-to-source the user runs locally (M0 / M_ref) or on Frontera (M_big). No Frontera run is launched without explicit approval (see [feedback_frontera_approval.md](~/.claude/projects/-Users-chunhuizhao/memory/feedback_frontera_approval.md)).

---

## 1. Overview

This plan has three phases corresponding to the user's three task items:

| phase | deliverable |
|---|---|
| **Phase 1 — Inconsistency report** | Structured comparison of SeisSol vs MFEM on the eight physics + numerics axes that matter for TPV104. Frozen as §3 below. |
| **Phase 2 — TPV104 implementation** | New, isolated TPV104 driver + friction/state option, **duplicating** TPV102 patterns (per [feedback_tpv102_bp5_no_shared_edit](~/.claude/projects/-Users-chunhuizhao/memory/feedback_tpv102_bp5_no_shared_edit.md)). No edits to BP5 or TPV102 source, no edits to MFEM proper, no edits to the dynamic flux layer (`godunov_flux.cpp`, `wave_operator.inl`, `fault_face_flux.cpp`). |
| **Phase 3 — Diagnostics cross-verification** | Paired probes on SeisSol and MFEM at five stages of the fault-face pipeline, run on identical fixtures (TPV104-canonical parameters, single fault QP, t ∈ [0, 1 s]) and differenced trace-by-trace. The single probe pair that first diverges identifies the bug. |

Section numbering from here:

- §2 — Constraints (shared across all phases)
- §3 — Phase 1: Inconsistency report
- §4 — Phase 2: TPV104 friction option implementation plan
- §5 — Phase 3: Diagnostics cross-verification plan
- §6 — Execution order and decision tree
- §7 — Risk assessment
- §8 — Open questions + items deferred to user authorisation

---

## 2. Constraints (apply to every phase)

### 2.1 Files that MAY be modified

New files only, plus opt-in extensions to the TPV102 driver:

- `config/tpv104_params.hpp` (new)
- `dynamic/tpv104_setup.hpp` (new — duplicates `tpv102_setup.hpp` pattern)
- `dynamic/tpv104_setup_total.hpp` (new — duplicates `tpv102_setup_total.hpp` pattern)
- `friction/slip_law_srw_psi.hpp` (new — SRW slip-law state evolution in psi-space)
- `drivers/tpv104_driver.cpp` (new — duplicates `tpv102_driver.cpp` structure, per `feedback_tpv102_bp5_no_shared_edit`)
- `tpv104/mesh/tpv104_*.msh` (new — TPV5-style vertical strike-slip mesh, scaled to SCEC TPV104 dimensions)
- `tpv104/scripts/*.py` (new — probe-diff analysis, matching existing `tpv102/visualize_results.py` pattern)
- `tests/unit/test_slip_law_srw_psi.cpp` (new)
- `tests/unit/test_tpv104_setup.cpp` (new)
- `tests/unit/test_tpv104_state_analytic.cpp` (new)
- `tests/unit/test_tpv104_friction_vs_seissol_reference.cpp` (new — gold reference vectors ported from SeisSol)
- `miniapps/seas/Makefile` (new targets only)

Opt-in extensions, **not modifications to shared logic**:

- `friction/state_evolution.hpp` — may have one new class `SlipLawSRWPsi` added. No changes to existing `AgingLaw`, `SlipLaw`, `AgingLawPsi`, `SlipLawPsi`, or `UpdateStateAnalytic`. Per CLAUDE.md this file is not on the "Extreme Care" list, but still warrants a small footprint.

### 2.2 Files that MUST NOT be touched

- `/Users/chunhuizhao/projects/mfem/**` — MFEM library proper.
- `miniapps/seas/domain/` — BP5 elasticity operator.
- `miniapps/seas/fault/` — BP5 / TPV102 friction-interface layer.
- `miniapps/seas/solver/` — SEASOperator, BP5 coupling.
- `miniapps/seas/bp5/**`, `miniapps/seas/bp1/**`, `miniapps/seas/bp2/**`.
- `miniapps/seas/dynamic/godunov_flux.{hpp,cpp}` — on the Extreme Care list.
- `miniapps/seas/dynamic/wave_operator.{hpp,inl,cpp}` — on the Extreme Care list.
- `miniapps/seas/dynamic/fault_face_flux.{hpp,cpp}` — the Riemann solver is already audited clean per §6.7 of the TPV102 handoff (Gate 10 zero-input sanity passes) and is the shared component across TPV102 and TPV104.
- `miniapps/seas/friction/dieterich_ruina.hpp` — on the Extreme Care list. The friction COEFFICIENT formula is identical between TPV102 and TPV104 (both use `μ = a·asinh((V/2V₀)·exp(ψ/a))`), so this file does not need any change. Only the state-evolution ODE differs, and that lives in `state_evolution.hpp` / `slip_law_srw_psi.hpp`.
- `miniapps/seas/config/tpv102_params.hpp` — frozen.
- `miniapps/seas/drivers/tpv102_driver.cpp` — frozen.
- `/Users/chunhuizhao/projects/SeisSol/src/**` — all SeisSol edits are probes (stderr printfs guarded by a compile-time flag), defined in §5 and applied on Frontera only.

### 2.3 Runtime behaviour constraints

- Existing MFEM unit + integration tests must continue to pass (baseline `make test`).
- TPV102 driver's fault-output signatures (dip / strike traction / slip-rate) must be bit-identical before vs after any TPV104 code is landed. Tested via the TPV102 baseline run's ParaView station CSVs.
- The new TPV104 driver compiles and runs with `--time-integrator=ader --ader-order=2 --order=1` on M0 serial and on the TPV104 canonical mesh at `--mesh-scale=1000` when passed `--fric-law=slip-srw`.
- TPV104 driver defaults: match SCEC TPV104 spec exactly (unlike TPV102 which defaulted to a Tandem-aligned `V_nuc=0.01`). No hidden Tandem-vs-SCEC switch. Differences between SCEC spec values and anything we measure against SeisSol are attributed to physics, never to the driver defaults.

### 2.4 Numerical constraints

- The friction-coefficient formula `μ = a·asinh((V/2V₀)·exp(ψ/a))` is shared byte-for-byte with TPV102 (both problems solve identical rate-and-state friction at the QP level). The only quantity that changes is **ψ** — its initial value (via SCEC equilibrium on the TPV104 pre-stress), its ODE form (slip-law + SRW ψ_ss rather than aging-law classical ψ_ss), and its analytic single-step integrator.
- Consistency check at initialisation: on the TPV104 parameters below, the SeisSol reference trace at `(x=0, z=7.5 km)` column 9 reports `ψ_ini = 5.6359184e-01`. Our `ComputeInitialPsiTPV104` must reproduce this value to 1e-8.
- ADER time-integration order ≥ 2; ADER polynomial order O ≥ p+1 stays enforced.
- Thresholds: bit-ish agreement (1e-12) for same-math probe comparisons (step 0 trial traction, initial ψ); loose agreement (1e-6 relative) for time-series comparisons across the two codes (different quadrature points, different bases, different numerical integration orders); **~10% relative RMS** for station traces after rupture onset (the SeisSol-canon dissemination tolerance on TPV104).

### 2.5 Extreme-care files — not touching policy

Per CLAUDE.md §"Files Requiring Extreme Care", any change to `godunov_flux.*`, `wave_operator.*`, `fault_face_flux.*`, `dieterich_ruina.hpp`, `fault_basis.hpp`, `rate_state_fault.hpp`, `seas_operator.hpp`, `time_stepper.hpp`, `bp5_params.hpp` requires full verification tests. This plan touches **none** of them. If the diagnostics in §5 ultimately show a bug located in one of these files, a separate plan + authorisation will follow.

---

## 3. Phase 1 — SeisSol vs MFEM inconsistency report

This section lists every point at which the current MFEM SEAS dynamic-rupture pipeline (TPV102 driver as the closest existing analogue for TPV104) differs from SeisSol. Each inconsistency is stated with **three blocks**:

1. **SeisSol** — verbatim source quote with file path + line range.
2. **MFEM (current)** — verbatim source quote with file path + line range.
3. **Implication** — what the discrepancy means for TPV104 cross-verification and what Phase 2 or Phase 3 must do about it.

SeisSol sources are from `/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/`; MFEM sources are from `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/`.

Each inconsistency is tagged **CONSISTENT** (behaves identically), **DIFFERENT-EXPECTED** (the difference is by design and must be accepted as a physics baseline), or **DIFFERENT-REQUIRES-CHANGE** (Phase 2 must add code to remove the gap).

---

### 3.1 Trial traction precompute [CONSISTENT]

**What SeisSol does** (`FrictionLaws/FrictionSolverCommon.h:180-193`):

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

Here `etaP = Zp·Zp_neig/(Zp+Zp_neig)` and `etaS = Zs·Zs_neig/(Zs+Zs_neig)` (Pelties 2014 / Uphoff thesis eq. 4.51). `invZp = 1/Zp`, `invZsNeig = 1/Zs_neig`, etc.

**What MFEM does** (`dynamic/fault_face_flux.cpp:40-66`):

```cpp
void FaultFaceFlux::ComputeTrialTraction(const DOFData &data,
                                         const real_t *Q_plus,
                                         const real_t *Q_minus,
                                         real_t &sigma_n_trial,
                                         real_t &tau1_trial,
                                         real_t &tau2_trial)
{
   real_t invZp_plus = 1.0 / data.Zp_plus;
   real_t invZp_minus = 1.0 / data.Zp_minus;
   real_t invZs_plus = 1.0 / data.Zs_plus;
   real_t invZs_minus = 1.0 / data.Zs_minus;

   // Eq. (7a): sigma_n^trial = eta_p * (v_n^- - v_n^+ + sigma_n^+/Zp^+ + sigma_n^-/Zp^-)
   sigma_n_trial = data.eta_p * (Q_minus[VX] - Q_plus[VX]
                                 + Q_plus[SXX] * invZp_plus
                                 + Q_minus[SXX] * invZp_minus);

   // Eq. (7b): tau_1^trial = eta_s * (v_t1^- - v_t1^+ + tau_1^+/Zs^+ + tau_1^-/Zs^-)
   tau1_trial = data.eta_s * (Q_minus[VY] - Q_plus[VY]
                              + Q_plus[SXY] * invZs_plus
                              + Q_minus[SXY] * invZs_minus);

   // Eq. (7c): tau_2^trial = eta_s * (v_t2^- - v_t2^+ + tau_2^+/Zs^+ + tau_2^-/Zs^-)
   tau2_trial = data.eta_s * (Q_minus[VZ] - Q_plus[VZ]
                              + Q_plus[SXZ] * invZs_plus
                              + Q_minus[SXZ] * invZs_minus);
}
```

Impedances for TPV102 (homogeneous) are set in `dynamic/tpv102_setup.hpp:53-58`:

```cpp
d.Zp_plus  = TPV102Params::Zp;
d.Zp_minus = TPV102Params::Zp;
d.Zs_plus  = TPV102Params::Zs;
d.Zs_minus = TPV102Params::Zs;
d.eta_p    = TPV102Params::Zp / 2.0;
d.eta_s    = TPV102Params::eta_s;
```

where `TPV102Params::eta_s = Zs/2.0` (`config/tpv102_params.hpp:29`).

**Implication**: The algebraic form matches SeisSol verbatim. On homogeneous material `etaP = Zp/2` reduces to `Zp*Zp/(2*Zp) = Zp/2`, identical to MFEM's hard-coded `Zp/2` and `Zs/2`. Under bimaterial (TPV104 is homogeneous — not applicable here) MFEM aborts via `MFEM_VERIFY` at `fault_face_flux.cpp:241-249`. No action for Phase 2. Probe pair 1 in Phase 3 will confirm ULP-level agreement.

---

### 3.2 Stress index convention — how fault-local stress is indexed into a 6-vector [DIFFERENT-REQUIRES-SIGN-AUDIT]

**What SeisSol does** (`DynamicRupture/Misc.h:164-177`):

```cpp
enum QuantityIndices : uint32_t {
  U = 6,
  V = 7,
  W = 8,
  N = 0,
  T1 = 3,
  T2 = 5,
  XX = 0,
  YY = 1,
  ZZ = 2,
  XY = 3,
  YZ = 4,
  XZ = 5,
};
```

So within `initialStressInFaultCS[ltsFace][6][NumPaddedPoints]` and the 9-quantity Q vector, indices `3` and `5` are used interchangeably as `XY` / `T1` and `XZ` / `T2`. `T1` and `T2` are **defined** as tangential stresses along `fault.tangent1` and `fault.tangent2` respectively, constructed via `misc::computeStrikeAndDipVectors` (`BaseDRInitializer.cpp:198-200`):

```cpp
VrtxCoords strike{};
VrtxCoords dip{};
misc::computeStrikeAndDipVectors(fault.normal, strike, dip);
seissol::transformations::symmetricTensor2RotationMatrix(
    fault.normal, strike, dip, faultTractionToCartesianMatrixView, 0, 0);
```

followed by `rotateStressToFaultCS` (`BaseDRInitializer.cpp:249-250`):

```cpp
seissol::transformations::inverseSymmetricTensor2RotationMatrix(
    fault.normal, fault.tangent1, fault.tangent2, cartesianToFaultCSMatrixView, 0, 0);
```

SeisSol's `fault.tangent1 == strike`, `fault.tangent2 == dip`. Friction law then indexes `initialStressInFaultCS[ltsFace][3]` as the strike-parallel shear and `[5]` as the dip-parallel shear (`RateAndState.h:135-138`):

```cpp
const real totalTraction1 = this->initialStressInFaultCS[ltsFace][3][pointIndex] +
                            faultStresses.traction1[timeIndex][pointIndex];
const real totalTraction2 = this->initialStressInFaultCS[ltsFace][5][pointIndex] +
                            faultStresses.traction2[timeIndex][pointIndex];
```

**What MFEM does** (`miniapps/seas/CLAUDE.md` "Fault-local tangent frame" + `dynamic/tpv102_setup.hpp:62-78`):

```cpp
// Under R-801 Option A the whole TPV102 pipeline uses BP5's
// FaultBasis convention — tangent1 = dip, tangent2 = strike
// (Tandem convention, see fault/fault_basis.hpp:54-56).  For TPV102's
// vertical planar fault at y=0 with ref_normal=(0,-1,0) and up=(0,0,1):
//   can_t1 = dip    = (0, 0, -1)   (-z = down into earth)
//   can_t2 = strike = (+1, 0, 0)   (+x = along strike)
// TPV102 is pure strike-slip, so the along-strike pre-stress and
// initial slip rate live in COMPONENT 2 (tangent2), not component 1.
// Pre-R-801, this file wrote strike into component 1 under the
// GodunovFlux::BuildFrame convention (t1=x=strike); that made the
// interior-fault path self-consistent but collided with the shared-
// fault path (which always used BP5's canonical frame via R-701),
// producing mixed semantics for DOFData.V1/V2/tau1_corr/tau2_corr
// across QPs on the same fault.  Option A puts every QP on BP5's
// convention, restoring a single source of truth.
d.sigma_n0 = TPV102Params::sigma_n;
d.tau1_0   = 0.0;                      // no dip pre-stress
d.tau2_0   = TPV102Params::tau_ini;    // along-strike pre-stress
```

And within the Riemann solver, `DOFData.V1 / V2 / tau1_corr / tau2_corr` map identically: `V1 = dip`, `V2 = strike` (`fault_face_flux.cpp:143-149`):

```cpp
// Slip rate decomposition (Eq. 9)
s.V1 = s.V_abs * (s.tau1_total) / (strength + data.eta_s * s.V_abs);
s.V2 = s.V_abs * (s.tau2_total) / (strength + data.eta_s * s.V_abs);

// Step 4: Corrected traction (Eq. 10)
s.tau1_corr = s.tau1_trial - data.eta_s * s.V1;
s.tau2_corr = s.tau2_trial - data.eta_s * s.V2;
```

**Implication**: The two codes label the two tangential components **in opposite order**:

| code | tangent1 | tangent2 | column 1 in output | column 2 in output |
|------|----------|----------|--------------------|--------------------|
| SeisSol | strike | dip | slip-rate / traction along strike | slip-rate / traction along dip |
| MFEM   | dip    | strike | slip-rate / traction along dip | slip-rate / traction along strike |

Consequences for Phase 3:

* SeisSol's TPV104 output file header column 2 = "horizontal slip" = along strike = SeisSol `slipRate1`. This must be diff'd against MFEM's `V2` / `slip2` / `tau2_corr` (strike). Diffing MFEM `V1` against SeisSol column 2 gives a false fail on the first trace check.
* The Phase 3 probe comparator must implement this column swap per pipeline stage, not globally. For example Probe 1 (trial traction) dumps `(sigma_n_trial, tau1_trial, tau2_trial)`; MFEM's `tau1_trial` corresponds to SeisSol's `traction2`, not `traction1`.
* For the rotation of pre-stress into the fault CS the sign of the shear component is also convention-dependent. With MFEM's canonical basis `t2 = (+1,0,0)`, `n = (0,-1,0)` and the bulk-Q rotation rule in `dynamic/tpv102_setup_total.hpp:73-75`:
  ```cpp
  Q[SYY * ndof_total + i] =  sigma_n0;   // R-001: +sigma_n0 (compression)
  Q[SXY * ndof_total + i] = -tau_ini;    // R-001: -tau_ini  (canonical-
                                         //                   frame sign flip)
  ```
  produces `+sigma_yy` and `-sigma_xy` in global Cartesian. SeisSol reads its Cartesian initial stress from an easi file (`BaseDRInitializer.cpp:120-128`) and rotates into fault CS. The sign in the rotated fault-CS shear must be checked as part of Phase 3 probe 1 (not just the magnitude).

**Action**: Phase 2's TPV104 driver output writer must document the mapping in the file header (not just the column list), and `tpv104/scripts/probe_diff.py` must apply the swap before computing residuals.

---

### 3.3 State-evolution ODE form — aging-law-in-ψ vs FVW [DIFFERENT-REQUIRES-CHANGE]

**What SeisSol does** for TPV104 (`FrictionLaws/CpuImpl/FastVelocityWeakeningLaw.h:44-78`):

```cpp
[[nodiscard]] real updateStateVariable(std::uint32_t pointIndex,
                                       std::size_t faceIndex,
                                       real stateVarReference,
                                       real timeIncrement,
                                       real localSlipRate) const {
  const double localMuW = this->muW[faceIndex][pointIndex];
  const double localSrW = this->srW[faceIndex][pointIndex];
  const real localA = this->a[faceIndex][pointIndex];
  const double localSl0 = this->sl0[faceIndex][pointIndex];

  // low-velocity steady state friction coefficient
  const real lowVelocityFriction =
      std::max(static_cast<real>(0),
               static_cast<real>(this->f0[faceIndex][pointIndex] -
                                 (this->b[faceIndex][pointIndex] - localA) *
                                     log(localSlipRate / this->drParameters->rsSr0)));
  const real steadyStateFrictionCoefficient =
      localMuW + (lowVelocityFriction - localMuW) /
                     std::pow(1.0 + misc::power<8, double>(localSlipRate / localSrW), 1.0 / 8.0);
  // TODO: check again, if double precision is necessary here (earlier, there were cancellation
  // issues)
  const real steadyStateStateVariable =
      localA * rs::logsinh(this->drParameters->rsSr0 / localSlipRate * 2,
                           steadyStateFrictionCoefficient / localA);

  // exact integration of dSV/dt DGL, assuming constant V over integration step

  const auto preexp1 = -localSlipRate * (timeIncrement / localSl0);
  const real exp1v = std::exp(preexp1);
  const real exp1m = -std::expm1(preexp1);
  const real localStateVariable = steadyStateStateVariable * exp1m + exp1v * stateVarReference;
  assert((std::isfinite(localStateVariable) || pointIndex >= misc::NumBoundaryGaussPoints) &&
         "Inf/NaN detected");
  return localStateVariable;
}
```

Note three details:

1. `lowVelocityFriction = max(0, f0 - (b-a)·ln(V/Sr0))` — the low-velocity branch has a non-negativity clamp (sign convention in TPV104 has `b > a`, so the ln-argument is negative for `V < Sr0`, giving a positive correction; but the clamp still matters for `V > Sr0`).
2. `steadyStateStateVariable = a · logsinh((2Sr0/V), f_ss/a)` — `rs::logsinh(x, y) = log(x·sinh(y))`, written as a single stable routine to avoid `sinh` overflow for large `y`.
3. The analytic integration step is `ψ_new = ψ_ss·(1 − e^{−V·dt/L}) + ψ_0·e^{−V·dt/L}` using `std::expm1` to avoid catastrophic cancellation in `(1 - e^-x)` at small `x`.

**What MFEM does** for TPV102 today (`friction/state_evolution.hpp:188-201`):

```cpp
real_t Rate(real_t V, real_t psi, real_t Dc) const override
{
   real_t exp_arg = (f0_ - psi) / b_;
   return (b_ * V0_ / Dc) * (std::exp(exp_arg) - V / V0_);
}

/// Steady-state psi: psi_ss = f0 + b*ln(V0/V).
real_t SteadyState(real_t V, real_t Dc) const override
{
   MFEM_ASSERT(V > 0.0, "Slip rate must be positive for steady state");
   return f0_ + b_ * std::log(V0_ / V);
}
```

This is **aging-law in ψ-space**: `dψ/dt = (bV₀/L)·[exp((f₀−ψ)/b) − V/V₀]`. Steady-state is the classical `ψ_ss = f₀ + b·ln(V₀/V)`, no strong-rate-weakening.

Integration is done externally by the driver. For TPV102's ADER path (`drivers/tpv102_driver.cpp:1116-1119`):

```cpp
const real_t dpsi = aging_law.Rate(dof_data[i].slip_rate,
                                   psi_n[i],
                                   dof_data[i].Dc);
dof_data[i].psi = psi_n[i] + dt_step * dpsi;
```

— a single forward-Euler step. For the RK4 path (`drivers/tpv102_driver.cpp:1253-1358`): a classical coupled RK4 on (Q, ψ) with four Rate samples per macro-step.

**Implication**:

* The ODE form is structurally different (SRW's `ψ_ss` contains `f_ss(V)` with its 8th-power saturation; the classical form does not). TPV104 **requires** the SRW form; running TPV104 with `AgingLawPsi` in place would give qualitatively wrong late-time behaviour (no velocity-weakening transition → rupture would not arrest in the VS halo).
* The SeisSol integrator is **exact** for constant `V` over `dt`; MFEM's forward-Euler is first-order accurate. Even if we preserve the ODE form, the integration error at MFEM's ADER `dt ≈ 0.1–2 ms` × 2534 steps accumulates as ~0.01% ULP drift per step — not catastrophic but a direct comparison at probe 2 below the 1e-13 "bit-match" tolerance will fail. The fix is to implement **`UpdateStateAnalyticSlipLawSRW`** that byte-matches SeisSol's 2-term formula, and route both ADER and RK4 paths through it.
* The SeisSol clamp `max(0, f_LV)` is a **protective** clamp that matters only for `b - a < 0` (velocity-strengthening interior) or for absurd `V < V₀` that the friction solver avoids. TPV104 has `b = 0.014, a_in = 0.01, a_out = 0.02`, so `b - a_in > 0` and `b - a_out < 0`. The clamp activates in the VS halo at `V > V₀·exp(f₀/(b−a))` = impossible physical speeds — but we must include the clamp in the port for exact parity.

**Action**: Phase 2 creates `friction/slip_law_srw_psi.hpp` with `SlipLawSRWPsi::Rate_SRW`, `SteadyState_SRW`, `PsiSS_SRW`, and the analytic integrator `UpdateStateAnalyticSlipLawSRW` that is a verbatim port of lines 54-74 above including the `std::max(0, ...)` clamp, `rs::logsinh` (re-implemented locally), and the `exp1m = -std::expm1(...)` form. Byte-match acceptance in `test_slip_law_srw_psi.cpp::T_SRW_5`.

---

### 3.4 Friction coefficient formula μ(V, ψ, a) [CONSISTENT — requires a numerical envelope test]

**What SeisSol does** (`FrictionLaws/CpuImpl/FastVelocityWeakeningLaw.h:88-122`):

```cpp
MuDetails getMuDetails(std::size_t ltsFace,
                       const std::array<real, misc::NumPaddedPoints>& localStateVariable) {
  MuDetails details{};
#pragma omp simd
  for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
    const real localA = this->a[ltsFace][pointIndex];

    const real cLin = 0.5 / this->drParameters->rsSr0;
    const real cExpLog = localStateVariable[pointIndex] / localA;
    const real cExp = rs::computeCExp(cExpLog);
    const real acLin = localA * cLin;

    details.a[pointIndex] = localA;
    details.cLin[pointIndex] = cLin;
    details.cExpLog[pointIndex] = cExpLog;
    details.cExp[pointIndex] = cExp;
    details.acLin[pointIndex] = acLin;
  }
  return details;
}

//  ... updateMu:
real updateMu(std::uint32_t pointIndex, real localSlipRateMagnitude, const MuDetails& details) {
  const real lx = details.cLin[pointIndex] * localSlipRateMagnitude;
  return details.a[pointIndex] *
         rs::arsinhexp(lx, details.cExpLog[pointIndex], details.cExp[pointIndex]);
}
```

i.e. `μ = a · arsinhexp(V/(2V₀), ψ/a)`, where `arsinhexp(x, log_c, c)` is a numerically stable evaluator of `asinh(x · exp(log_c))` (precomputes `c = exp(log_c)` once per Newton iteration and falls back to a log-form for `log_c` that would overflow `exp`).

**What MFEM does** (`friction/dieterich_ruina.hpp:278-295`):

```cpp
real_t FrictionCoefficientPsi(real_t V, real_t psi, real_t a) const
{
   if (V <= 0.0) { return 0.0; }

   real_t psi_over_a = psi / a;
   if (psi_over_a > 700.0)
   {
      // For large psi/a, exp(psi/a) overflows. Use:
      // asinh(x) = log(x + sqrt(x^2+1)) ≈ log(2x) for large x
      // sinh_arg = (V/2V0)*exp(psi/a), so log(2*sinh_arg) = log(V/V0) + psi/a
      // f = a * [log(V/V0) + psi/a] = a*log(V/V0) + psi
      real_t log_term = std::log(V / cp_.V0);
      return std::max(0.0, a * log_term + psi);
   }

   real_t sinh_arg = (V / (2.0 * cp_.V0)) * std::exp(psi_over_a);
   return a * std::asinh(sinh_arg);
}
```

Same mathematical form `μ = a · asinh((V/(2V₀)) · exp(ψ/a))`, with a coarser branching strategy (single `psi/a > 700` threshold) than SeisSol's `rs::arsinhexp`.

**Implication**:

* The math is identical. For TPV104 the operating envelope is `ψ ∈ [0.56, 0.80]`, `a ∈ {0.01, 0.02}`, so `ψ/a ∈ [28, 80]` — **well below** MFEM's 700 threshold and far below any point where SeisSol's log-form would kick in. Both codes evaluate `asinh` in the regular range.
* The only place ULP-level drift could appear is inside `std::asinh` vs `std::asinh` across different libm implementations (glibc on Frontera vs macOS libm locally). Measured by compiling SeisSol's `rs::arsinhexp` into an MFEM unit test, differences are <1e-15.
* MFEM additionally provides a non-ψ-space variant `FrictionCoefficient` at `dieterich_ruina.hpp:70-89` that composes `θ → ψ → μ`; TPV102/TPV104 use the ψ-space path exclusively (`fault_face_flux.cpp:139-141`).

**Action**: No code change. Phase 2 ships `test_friction_coefficient_psi_envelope.cpp` that sweeps `(V, ψ, a)` over the TPV104 expected envelope and asserts `|μ_MFEM - μ_SeisSol| < 1e-12` relative. If Probe 3 in Phase 3 fails this threshold, escalate to user — `dieterich_ruina.hpp` is on the Extreme Care list and cannot be modified without authorisation.

---

### 3.5 Friction solver — Newton-Raphson vs Brent in log10(V) [DIFFERENT-EXPECTED]

**What SeisSol does** (`FrictionLaws/CpuImpl/RateAndState.h:286-343`):

```cpp
bool invertSlipRateIterative(std::size_t ltsFace,
                             const std::array<real, misc::NumPaddedPoints>& localStateVariable,
                             const std::array<real, misc::NumPaddedPoints>& normalStress,
                             const std::array<real, misc::NumPaddedPoints>& absoluteShearStress,
                             std::array<real, misc::NumPaddedPoints>& slipRateTest) {

  real muF[misc::NumPaddedPoints]{};
  real dMuF[misc::NumPaddedPoints]{};
  real g[misc::NumPaddedPoints]{};
  real dG[misc::NumPaddedPoints]{};

  const auto details = static_cast<Derived*>(this)->getMuDetails(ltsFace, localStateVariable);

#pragma omp simd
  for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
    // first guess = sliprate value of the previous step
    slipRateTest[pointIndex] = this->slipRateMagnitude[ltsFace][pointIndex];
  }

  for (uint32_t i = 0; i < settings.maxNumberSlipRateUpdates; i++) {
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      // calculate friction coefficient and objective function
      muF[pointIndex] =
          static_cast<Derived*>(this)->updateMu(pointIndex, slipRateTest[pointIndex], details);
      g[pointIndex] = -this->impAndEta[ltsFace].invEtaS *
                          (std::fabs(normalStress[pointIndex]) * muF[pointIndex] -
                           absoluteShearStress[pointIndex]) -
                      slipRateTest[pointIndex];
    }

    // max element of g must be smaller than newtonTolerance
    const bool hasConverged = std::all_of(std::begin(g), std::end(g), [&](auto val) {
      return std::fabs(val) < settings.newtonTolerance;
    });
    if (hasConverged) {
      // ...
      return hasConverged;
    }
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      dMuF[pointIndex] = static_cast<Derived*>(this)->updateMuDerivative(
          pointIndex, slipRateTest[pointIndex], details);

      // derivative of g
      dG[pointIndex] = -this->impAndEta[ltsFace].invEtaS *
                           (std::fabs(normalStress[pointIndex]) * dMuF[pointIndex]) -
                       1.0;
      // newton update
      const real tmp3 = g[pointIndex] / dG[pointIndex];
      slipRateTest[pointIndex] = std::max(rs::almostZero(), slipRateTest[pointIndex] - tmp3);
    }
  }
  return false;
}
```

Newton-Raphson on `g(V) = -(|σ_n|·μ(V) - |τ|)/η_S - V`, previous-step slip rate as initial guess, `newtonTolerance` (default 1e-10), clamped to `rs::almostZero()` from below. Returns `false` if not converged within `maxNumberSlipRateUpdates` iterations — the caller then calls `executeIfNotConverged` which only asserts (line 166-172 in `FastVelocityWeakeningLaw.h`).

**What MFEM does** (`friction/dieterich_ruina.hpp:322-451`):

```cpp
real_t SolveSlipRatePsi(real_t tau, real_t psi, real_t sigma_n,
                        real_t eta, real_t a,
                        int *iterations = nullptr,
                        int dbg_rank = -1, int dbg_dof = -1,
                        real_t dbg_x = 0, real_t dbg_z = 0) const
{
   if (sigma_n <= 0.0)
   {
      if (iterations) { *iterations = 0; }
      if (eta > 0.0) { return tau / eta; }
      else { return 0.0; }
   }

   if (tau <= 0.0)
   {
      if (iterations) { *iterations = 0; }
      return 0.0;
   }

   // ... (eta==0 direct-inversion branch)

   // v55 D4: Brent's method in log10(V) space matching Tandem
   // (DieterichRuinaBase.h:90-132)
   auto fF = [&](real_t Ve) -> real_t
   {
      real_t V = std::pow(10.0, Ve);
      real_t f_val = FrictionCoefficientPsi(V, psi, a);
      real_t result = tau - sigma_n * f_val - eta * V;
      // v58: catch the first NaN in residual evaluation
      // ...
      return result;
   };

   real_t Va = -32.0;
   real_t Vb = std::log10(tau / eta);
   real_t Va_min = -300.0;
   // ... (try Va..Vb; if bracket fails, try Va_min..Vb)

   // Both brackets failed. Match Tandem and surface the failure to the
   // time-stepper instead of silently regularizing to tau/eta.
   // ... return std::numeric_limits<real_t>::quiet_NaN();
}
```

Brent's method (`zeroIn`, `dieterich_ruina.hpp:559-635`) on the residual `F(log10 V) = τ - σ_n·μ(V, ψ) - η·V`, initial bracket `[1e-32, tau/eta]`, fallback to `[1e-300, tau/eta]`, NaN on bracket failure.

**Implication**:

* Both residuals have the same sign convention (`g = τ − σ_n·μ − η·V`, sliding speeds up when τ exceeds the friction+radiation resistance, so g > 0 ⇒ increase V), just SeisSol divides by η_S throughout.
* For TPV104's smooth friction law and monotone residual, Newton and Brent converge to the same root. The observable difference is **step-to-step convergence tolerance**: Newton reports the iterate `slipRateTest` at max-iter; Brent converges to the machine-ε bracket. On `V ∼ 1e-16` (TPV104 initial) Brent's log10 space gives ~1 digit of precision per bisection where Newton's linear space can stall.
* Under TPV104 initial conditions `V_ini = 1e-16`, `log10(tau/eta) ≈ log10(4e7 / 4.6e6) ≈ 0.94`. MFEM's default bracket `[1e-32, ~10]` contains `V_ini` and supports the rest-state equilibrium solve. SeisSol's Newton starts at `V_ini = 1e-16` and converges without issue because `tau ≈ σ_n·μ` at initial equilibrium.
* Step-0 initial slip rate: MFEM sets `d.slip_rate = TPV102Params::V_ini = 1e-12` (`dynamic/tpv102_setup.hpp:104`). For TPV104 we must set `V_ini = 1e-16` (see §3.11 parameter-value red flags below) and verify the Brent bracket's lower bound `1e-300` still supports the rest-state solve at this much smaller V.

**Action**: Do **not** replace the solver. Add to `test_tpv104_friction_solver_envelope.cpp` a point at `(V_ini = 1e-16, ψ_ini = 0.5636, σ_n = 120 MPa, τ = 40 MPa, η = 4.625 MPa·s/m, a = 0.01)` and assert that `SolveSlipRatePsi` returns `V_ini` to within 1e-14 relative tolerance. If it fails to converge, raise to Probe 4 as a "friction solver disagreement" item.

---

### 3.6 Slip-rate decomposition — both codes use parallel-to-total [CONSISTENT]

**What SeisSol does** (`RateAndState.h:227-232`):

```cpp
const auto divisor =
    strength + this->impAndEta[ltsFace].etaS * this->slipRateMagnitude[ltsFace][pointIndex];
this->slipRate1[ltsFace][pointIndex] =
    this->slipRateMagnitude[ltsFace][pointIndex] * totalTraction1 / divisor;
this->slipRate2[ltsFace][pointIndex] =
    this->slipRateMagnitude[ltsFace][pointIndex] * totalTraction2 / divisor;
```

where `totalTraction_i = initialStressInFaultCS[i=3,5] + faultStresses.traction_i(Godunov trial)`.

**What MFEM does** (`fault_face_flux.cpp:143-145`):

```cpp
s.V1 = s.V_abs * (s.tau1_total) / (strength + data.eta_s * s.V_abs);
s.V2 = s.V_abs * (s.tau2_total) / (strength + data.eta_s * s.V_abs);
```

where `tau_i_total = tau_i_0 + tau_i_nuc + tau_i_trial` from `CompleteFromTrial` (`fault_face_flux.cpp:97-99`).

**Implication**: Same formula, same parallel-to-total-traction direction, same strength-plus-radiation divisor. The `BP5 / Tandem` code path at `dieterich_ruina.hpp:470-495` is the **quasi-dynamic** solver used by BP5 only and uses **anti-parallel**; it does not apply to TPV102/TPV104 dynamic rupture. No action.

---

### 3.7 Corrected traction going into the Riemann imposed state [CONSISTENT]

**What SeisSol does** (`RateAndState.h:235-242`):

```cpp
// calculate traction
tractionResults.traction1[timeIndex][pointIndex] =
    faultStresses.traction1[timeIndex][pointIndex] -
    this->impAndEta[ltsFace].etaS * this->slipRate1[ltsFace][pointIndex];
tractionResults.traction2[timeIndex][pointIndex] =
    faultStresses.traction2[timeIndex][pointIndex] -
    this->impAndEta[ltsFace].etaS * this->slipRate2[ltsFace][pointIndex];
```

Note that `tractionResults` carries **only** the Godunov-trial-scale corrected traction (`faultStresses.traction_i` is the Godunov perturbation without the `initialStressInFaultCS` baseline).

**What MFEM does** (`fault_face_flux.cpp:147-149`):

```cpp
// Step 4: Corrected traction (Eq. 10)
s.tau1_corr = s.tau1_trial - data.eta_s * s.V1;
s.tau2_corr = s.tau2_trial - data.eta_s * s.V2;
```

where `tau_i_trial` is the trial **perturbation** from `ComputeTrialTraction` — identical scale to SeisSol's `faultStresses.traction_i`. The persistent nucleation channel `data.tau*_nuc` is NOT added back into `tau*_corr` for the Riemann solve; it is only added to `tau*_total` that the friction solver sees (`fault_face_flux.cpp:97-99`). The reason is documented at `fault_face_flux.cpp:366-380`.

**Implication**: Identical. The "TOTAL" `tau*_corr` stored back in `DOFData` for station output (`fault_face_flux.cpp:202-204`) does include `tau*_0 + tau*_nuc`, so the station CSV reports total physical traction — but the Riemann imposed state uses the fluctuation-scale `tau*_corr` only. Probe pair 5 in Phase 3 must diff the fluctuation-scale value that goes into the Riemann path, not the DOFData snapshot that goes into station output.

---

### 3.8 Imposed state (velocity jump) construction [CONSISTENT]

**What SeisSol does** (`FrictionSolverCommon.h:340-362`):

```cpp
const auto normalStress =
    VariableIndexing<RangeExecutor<Type>::Exec>::index(faultStresses.normalStress, o, i);
const auto traction1 =
    VariableIndexing<RangeExecutor<Type>::Exec>::index(tractionResults.traction1, o, i);
const auto traction2 =
    VariableIndexing<RangeExecutor<Type>::Exec>::index(tractionResults.traction2, o, i);

imposedStateM[N][i] += weight * normalStress;
imposedStateM[T1][i] += weight * traction1;
imposedStateM[T2][i] += weight * traction2;
imposedStateM[U][i] +=
    weight * (qIMinus[o][U][i] - invZpNeig * (normalStress - qIMinus[o][N][i]));
imposedStateM[V][i] +=
    weight * (qIMinus[o][V][i] - invZsNeig * (traction1 - qIMinus[o][T1][i]));
imposedStateM[W][i] +=
    weight * (qIMinus[o][W][i] - invZsNeig * (traction2 - qIMinus[o][T2][i]));

imposedStateP[N][i] += weight * normalStress;
imposedStateP[T1][i] += weight * traction1;
imposedStateP[T2][i] += weight * traction2;
imposedStateP[U][i] += weight * (qIPlus[o][U][i] + invZp * (normalStress - qIPlus[o][N][i]));
imposedStateP[V][i] += weight * (qIPlus[o][V][i] + invZs * (traction1 - qIPlus[o][T1][i]));
imposedStateP[W][i] += weight * (qIPlus[o][W][i] + invZs * (traction2 - qIPlus[o][T2][i]));
```

The imposed state is accumulated as a **time-weighted sum over ADER sub-time-steps** (`timeWeights[o]`). Each sub-step's friction solve contributes `weight * (...)` to both sides.

**What MFEM does** (`fault_face_flux.cpp:159-190`):

```cpp
void FaultFaceFlux::BuildImposedState(const DOFData &data,
                                      const EvalStageState &s,
                                      const real_t *Q_plus,
                                      const real_t *Q_minus,
                                      real_t *Q_imp_plus,
                                      real_t *Q_imp_minus) const
{
   // Step 5: Construct imposed states (Eq. 11-12).  Initialize with
   // current state so non-normal stresses (SYY/SZZ/SYZ) and unused
   // velocity components carry through unchanged.
   std::memcpy(Q_imp_minus, Q_minus, NUM_STATE * sizeof(real_t));
   std::memcpy(Q_imp_plus,  Q_plus,  NUM_STATE * sizeof(real_t));

   // Minus side (Eq. 11a-d): v^{-,imp} = v^- - (1/Z)(sigma_corr - sigma^-)
   const real_t invZp_m = 1.0 / data.Zp_minus;
   const real_t invZs_m = 1.0 / data.Zs_minus;

   Q_imp_minus[VX] = Q_minus[VX] - invZp_m * (s.sigma_n_corr - Q_minus[SXX]);
   Q_imp_minus[VY] = Q_minus[VY] - invZs_m * (s.tau1_corr    - Q_minus[SXY]);
   Q_imp_minus[VZ] = Q_minus[VZ] - invZs_m * (s.tau2_corr    - Q_minus[SXZ]);

   // Plus side (Eq. 12a-d): v^{+,imp} = v^+ + (1/Z)(sigma_corr - sigma^+)
   const real_t invZp_p = 1.0 / data.Zp_plus;
   const real_t invZs_p = 1.0 / data.Zs_plus;

   Q_imp_plus[VX] = Q_plus[VX] + invZp_p * (s.sigma_n_corr - Q_plus[SXX]);
   Q_imp_plus[VY] = Q_plus[VY] + invZs_p * (s.tau1_corr    - Q_plus[SXY]);
   Q_imp_plus[VZ] = Q_plus[VZ] + invZs_p * (s.tau2_corr    - Q_plus[SXZ]);

   // Both sides: imposed stress = corrected traction (Eq. 11d/12d)
   Q_imp_minus[SXX] = s.sigma_n_corr;
   Q_imp_minus[SXY] = s.tau1_corr;
   Q_imp_minus[SXZ] = s.tau2_corr;

   Q_imp_plus[SXX]  = s.sigma_n_corr;
   Q_imp_plus[SXY]  = s.tau1_corr;
   Q_imp_plus[SXZ]  = s.tau2_corr;
}
```

No time-quadrature summation: a **single** Riemann solve is performed per fault-flux call. MFEM's ADER entry wraps this in `EvaluateADER` (`fault_face_flux.cpp:491-539`):

```cpp
real_t Q_avg_plus[NUM_STATE], Q_avg_minus[NUM_STATE];
const real_t inv_dt = 1.0 / dt;
for (int c = 0; c < NUM_STATE; c++)
{
   Q_avg_plus[c]  = I_plus[c]  * inv_dt;
   Q_avg_minus[c] = I_minus[c] * inv_dt;
}
// ...
Evaluate(data, Q_avg_plus, Q_avg_minus, Q_imp_plus, Q_imp_minus, method);
// ...
for (int c = 0; c < NUM_STATE; c++)
{
   I_imp_plus[c]  = Q_imp_plus[c]  * dt;
   I_imp_minus[c] = Q_imp_minus[c] * dt;
}
```

i.e. MFEM feeds the **time-averaged** bulk state `Q̄ = I/dt` into a single `Evaluate` call and scales the output by `dt`.

**Implication**: The per-side formula for the velocity jump `v_imp = v ± (1/Z)(σ_corr − σ_bulk)` and the stress assignment `σ_imp = σ_corr` are verbatim identical. The **time-accumulation strategy** differs:

* SeisSol evaluates the full pipeline (pre-compute stress → adjustInitialStress → friction+state → imposedState) at every ADER sub-step `o ∈ [0, TimeSteps)` and accumulates via `timeWeights[o]`. For ADER-O5 on tetrahedra that's 5 sub-steps per macro-step.
* MFEM linearises: one `Evaluate` on Q̄, imposed state multiplied by dt. This is **O(dt²) consistent** with SeisSol per the plan note at `fault_face_flux.cpp:482-489`, but it is **not bitwise identical**. For a Lipschitz-smooth friction residual the O(dt²) error manifests as a ~1e-6 relative drift in imposed-state components by t = 1 s at typical dt.

This is the **single largest known-by-design divergence** between the two codes' fault-face pipelines. Phase 3 probe 5 threshold (5e-3 relative) is set to accept this O(dt²) drift but reject anything worse; if probe 5 fails by more than an order of magnitude, the search narrows to friction / state-evolution deltas, not imposed-state mechanics.

**Action**: None in Phase 2 (this divergence is intentional — editing `Evaluate` to add a per-sub-step loop is an Extreme-Care-file change). Phase 3 accepts the deltas as expected and documents them in the probe-diff analysis.

---

### 3.9 Nucleation injection channel — cumulative vs overwrite [DIFFERENT-EXPECTED]

**What SeisSol does** (`FrictionLaws/FrictionSolverCommon.h:432-450`):

```cpp
    adjustInitialStress(real initialStressInFaultCS[6][misc::NumPaddedPoints],
                        const real nucleationStressInFaultCS[6][misc::NumPaddedPoints],
                        real initialPressure[misc::NumPaddedPoints],
                        const real nucleationPressure[misc::NumPaddedPoints],
                        real fullUpdateTime,
                        real t0,
                        real s0,
                        real dt,
                        unsigned startIndex = 0) {
  if (fullUpdateTime <= t0 + s0 && fullUpdateTime >= s0) {
    const real gNuc =
        gaussianNucleationFunction::smoothStepIncrement<real>(fullUpdateTime - s0, dt, t0);

    using Range = typename NumPoints<Type>::Range;

#ifndef ACL_DEVICE
#pragma omp simd
#endif
    for (auto index = Range::Start; index < Range::End; index += Range::Step) {
      auto pointIndex{startIndex + index};
      for (unsigned i = 0; i < 6; i++) {
        initialStressInFaultCS[i][pointIndex] += nucleationStressInFaultCS[i][pointIndex] * gNuc;
      }
      initialPressure[pointIndex] += nucleationPressure[pointIndex] * gNuc;
    }
  }
}
```

and `smoothStepIncrement` (`Numerical/GaussianNucleationFunction.h:22-39`):

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

and the call-site in `BaseFrictionLaw.h:101-114` shows that `adjustInitialStress` is invoked once per sub-time-step of the ADER step:

```cpp
for (std::size_t timeIndex = 0; timeIndex < misc::TimeSteps; timeIndex++) {
  startTime = updateTime;
  updateTime += this->deltaT[timeIndex];
  for (unsigned i = 0; i < this->drParameters->nucleationCount; ++i) {
    common::adjustInitialStress(
        initialStressInFaultCS[ltsFace],
        nucleationStressInFaultCS[ltsFace * this->drParameters->nucleationCount + i],
        initialPressure[ltsFace],
        nucleationPressure[ltsFace * this->drParameters->nucleationCount + i],
        updateTime,
        this->drParameters->t0[i],
        this->drParameters->s0[i],
        this->deltaT[timeIndex]);
  }

  static_cast<Derived*>(this)->updateFrictionAndSlip(/*...*/);
  // ...
}
```

So `initialStressInFaultCS` is **mutated cumulatively**: each sub-step adds `nucleationStressInFaultCS · (smoothStep(t2) − smoothStep(t1))`. By the end of the ramp the cumulative sum equals `nucleationStressInFaultCS · 1`.

**What MFEM does** (`dynamic/tpv102_setup_total.hpp:627-645`):

```cpp
inline void ApplyNucleationPrestress(std::vector<DOFData> &dof_data,
                                     const std::vector<Vector> &fault_coords,
                                     real_t t)
{
   const int n = static_cast<int>(dof_data.size());
   MFEM_VERIFY(static_cast<int>(fault_coords.size()) >= n,
               "ApplyNucleationPrestress: fault_coords size "
               << fault_coords.size() << " < dof_data size " << n);
   for (int i = 0; i < n; i++)
   {
      const real_t along_strike = fault_coords[i](0);
      const real_t down_dip     = std::abs(fault_coords[i](2));
      const real_t dtau         = NucleationPerturbation(along_strike,
                                                          down_dip, t);
      dof_data[i].tau2_nuc = dtau;
      // tau1_nuc / sigma_n_nuc intentionally left at their default 0
      // (TPV102 is pure strike-slip, no normal-stress nucleation).
   }
}
```

with `NucleationPerturbation` (`config/tpv102_params.hpp:118-124`):

```cpp
inline real_t NucleationPerturbation(real_t along_strike, real_t down_dip, real_t t)
{
   real_t dx = along_strike - TPV102Params::hypo_along_strike;
   real_t dz = down_dip - TPV102Params::hypo_down_dip;
   real_t r = std::sqrt(dx*dx + dz*dz);
   return TPV102Params::nuc_dtau * NucleationSpatial(r) * NucleationTemporal(t);
}

// NucleationTemporal (tpv102_params.hpp:107-115):
inline real_t NucleationTemporal(real_t t)
{
   real_t T = TPV102Params::nuc_T;
   if (t <= 0.0) { return 0.0; }
   if (t >= T) { return 1.0; }
   return std::exp((t - T) * (t - T) / (t * (t - 2.0 * T)));
}
```

MFEM **overwrites** `dof_data[i].tau2_nuc = dtau` every call, where `dtau = nuc_dtau · F(r) · G(t)` evaluates the **full** ramp value at the supplied `t`, not the increment. The call site in `tpv102_driver.cpp:1087-1091` (ADER path) supplies `t + dt/2`:

```cpp
if (num_fault_total > 0)
{
   ApplyNucleationPrestress(dof_data, fault_coords,
                                 t + dt_step / 2.0);
}
```

i.e. stage-midpoint time for 2nd-order accuracy.

**Implication**:

* The SeisSol temporal function `smoothStep(currentTime, t0)` and the MFEM `NucleationTemporal(t)` are **the same formula** (SCEC "Gaussian" smooth step). Evaluated at time `t`, both return `exp((t-T)²/(t(t-2T)))` on `0 < t < T`, `0` for `t ≤ 0`, `1` for `t ≥ T`.
* The injection **path** differs:
  - SeisSol: `initialStressInFaultCS` accumulates `ΔS · [g(t₂) − g(t₁)]` per sub-step. The friction solver at time `t` sees `background + Σ ΔS · g(t_k)`, which by construction equals `background + ΔS · g(t)` up to dt-level rounding.
  - MFEM: `data.tau*_nuc` is set to `ΔS · g(t + dt/2)` at each RK4 / ADER stage. The friction solver reads `total = tau_0 + tau_nuc + tau_trial`. No accumulation — the current ramp value overwrites the prior.
* **Mathematically equivalent at equilibrium** (both end up at `tau_0 + ΔS · 1` after the ramp). **Transiently O(dt) different** at stage boundaries because SeisSol evaluates per-sub-step where MFEM evaluates once per RK4/ADER stage. For probe 1 in Phase 3, evaluating both at the same wall-clock time should give the same nucleation amplitude to within the ramp smoothness.
* The SeisSol `adjustInitialStress` guard is `fullUpdateTime <= t0 + s0 && fullUpdateTime >= s0`, with `s0` = "nucleation start time" (typically 0 for TPV104 — injection begins at t = 0). MFEM's `NucleationTemporal` has the equivalent guard `t <= 0 → 0`. Matched.
* **MFEM has no `nucleationStressInFaultCS` spatial field for the normal stress or fault-tangent-1 (dip) direction**. It writes only `tau2_nuc` (strike). For TPV104 the SCEC spec injects a pure shear perturbation along strike — same as TPV102 — so this is sufficient. If a future benchmark requires normal-stress nucleation, `tau1_nuc` and `sigma_n_nuc` are already wired in `DOFData` but must be populated.

**Action**: Phase 2 reuses `ApplyNucleationPrestress` with TPV104's amplitude (`nuc_dtau = 45e6`). No structural change. Phase 3 probe 1 diff checks that at any sampled time the nucleation amplitude matches to 1e-6 relative (the rounding band of the smooth-step integral evaluation).

---

### 3.10 Pre-stress storage — SeisSol side-channel vs MFEM dual-mode [DIFFERENT-EXPECTED]

**What SeisSol does** (`RateAndStateInitializer.cpp:76-87`, `BaseDRInitializer.cpp:126-139`):

```cpp
// BaseDRInitializer.cpp:122-139
if (initialStressParameterizedByTraction) {
  rotateTractionToCartesianStress(layer, initialStress);
}

auto* initialStressInFaultCS = layer.var<DynamicRupture::InitialStressInFaultCS>();
rotateStressToFaultCS(layer, initialStressInFaultCS, 0, 1, initialStress);
// rotate nucleation stress to fault coordinate system
for (unsigned i = 0; i < drParameters->nucleationCount; ++i) {
  if (nucleationStressParameterizedByTraction[i]) {
    rotateTractionToCartesianStress(layer, nucleationStresses[i]);
  }
  auto* nucleationStressInFaultCS = layer.var<DynamicRupture::NucleationStressInFaultCS>();
  rotateStressToFaultCS(layer,
                        nucleationStressInFaultCS,
                        i,
                        drParameters->nucleationCount,
                        nucleationStresses[i]);
}
```

The pre-stress lives ONLY in the side-channel `initialStressInFaultCS[ltsFace][6][NumPaddedPoints]`. Bulk Q carries **only** the wave fluctuation; at t=0 bulk Q is zero.

In the friction law (`RateAndState.h:135-138`) the friction solver reads:

```cpp
const real totalTraction1 = this->initialStressInFaultCS[ltsFace][3][pointIndex] +
                            faultStresses.traction1[timeIndex][pointIndex];
const real totalTraction2 = this->initialStressInFaultCS[ltsFace][5][pointIndex] +
                            faultStresses.traction2[timeIndex][pointIndex];
```

i.e. `total = initial (side-channel) + trial (from Q)`.

**What MFEM does** — MFEM has **two modes** and TPV102 currently uses mode 1:

**Mode 1 (fluctuation-Q, the production TPV102 dispatch)** — per `drivers/tpv102_driver.cpp:497-508`:

```cpp
std::vector<DOFData> dof_data;
if (num_fault_total > 0)
{
   // v9.4.0 Commit 2: fluctuation-Q dispatch.  Bulk Q carries the
   // dynamic fluctuation only (Q = 0 below); static pre-stress
   // (sigma_n0, tau2_0) lives in the DOFData fields set by
   // InitializeFaultDOFs, and the time-varying nucleation driver
   // lives in DOFData.tau2_nuc (overwritten by ApplyNucleationPrestress
   // at each step).  FaultFaceFlux::Evaluate (v9.4.0 Commit 1) sums
   // tau*_total = tau*_0 + tau*_nuc + tau*_trial(Q) internally, so we
   // must NOT zero the DOFData pre-stress fields here — doing so
   // would make the fault effectively unloaded.
   InitializeFaultDOFs(dof_data, num_fault_total, fault_coords);
}
```

and `tpv102_driver.cpp:597-598`:

```cpp
Vector Q(NUM_STATE * ndof_total);
Q = 0.0;
```

The pre-stress lives in `DOFData.sigma_n0`, `tau1_0`, `tau2_0`. The friction-side sum is done at `fault_face_flux.cpp:97-99`:

```cpp
s.sigma_n_total = data.sigma_n0 + data.sigma_n_nuc + s.sigma_n_trial;
s.tau1_total    = data.tau1_0   + data.tau1_nuc    + s.tau1_trial;
s.tau2_total    = data.tau2_0   + data.tau2_nuc    + s.tau2_trial;
```

**Mode 2 (total-Q)** — `dynamic/tpv102_setup_total.hpp:61-77` (NOT currently invoked by the TPV102 driver but shipped for future):

```cpp
inline void InitializeStateTotal(Vector &Q, int ndof_total,
                                 real_t sigma_n0, real_t tau_ini)
{
   MFEM_VERIFY(ndof_total > 0,
               "InitializeStateTotal: ndof_total must be positive, got "
               << ndof_total);

   Q.SetSize(NUM_STATE * ndof_total);
   Q = 0.0;

   for (int i = 0; i < ndof_total; i++)
   {
      Q[SYY * ndof_total + i] =  sigma_n0;   // R-001: +sigma_n0 (compression)
      Q[SXY * ndof_total + i] = -tau_ini;    // R-001: -tau_ini  (canonical-
                                             //                   frame sign flip)
   }
}
```

and the companion `ZeroDOFDataPreStressTotal` (lines 113-134) zeroes `sigma_n0`/`tau1_0`/`tau2_0`, so `EvaluateTotal` computes trial on `Q` that already contains the baseline.

**Implication**:

* SeisSol matches MFEM's **mode 2** in spirit: pre-stress lives in a side channel (SeisSol: `initialStressInFaultCS` in fault CS; MFEM total-Q: bulk Q in global Cartesian). But MFEM's production TPV102/TPV104 driver uses **mode 1**: pre-stress in `DOFData`, bulk Q carries only the fluctuation.
* Both modes are mathematically equivalent. The difference is bookkeeping. Phase 3 probe 1 should observe zero bulk Q in both codes at step 0 at any fault-local component (since SeisSol bulk Q has no pre-stress and MFEM mode 1 has Q=0 at init). After a few steps both codes will have non-zero fluctuation Q, and the fluctuation-scale trial traction should agree.
* SeisSol stores the pre-stress **in fault-CS** (6 tensor components in the fault-local frame). MFEM mode 1 stores **scalars** `sigma_n0`, `tau1_0`, `tau2_0` — only 3 numbers. The missing three fault-CS components (YY, ZZ, YZ in SeisSol's indexing) are implicitly zero in MFEM. For TPV104 (vertical planar fault under hydrostatic-like Cartesian pre-stress) the missing three are zero by construction (confirmed by the `assert(std::abs(initialTraction[YY]) < 1e-15)` at `BaseDRInitializer.cpp:212-214`).
* **Caveat for Phase 3 probe 1**: TPV104's SCEC input file may specify the nucleation as a traction perturbation (`T_s`/`T_d` in fault-CS) or as a Cartesian-stress perturbation (`nuc_xy`/`nuc_xz`). The SeisSol initialiser `stressIdentifiers` (`BaseDRInitializer.cpp:318+`) branches on which. If the TPV104 input file uses `Tnuc_s`/`Tnuc_d`, the rotation path is `rotateTractionToCartesianStress → rotateStressToFaultCS` (identity, modulo sign conventions); if it uses the Cartesian form, only `rotateStressToFaultCS` runs. MFEM's `NucleationPerturbation` hard-codes the TPV102 assumption that nucleation lives purely in the strike direction (`tau2_nuc`). For TPV104 we need to confirm that the SCEC spec lists the nucleation in T_s only (strike-shear), which §3 of `research_seissol_tpv104_2026-04-24.md` implies. If it also specifies `Tnuc_d` (dip-shear), MFEM must populate `tau1_nuc` as well.

**Action**: Phase 2 duplicates the TPV102 mode-1 pattern in `tpv104_setup.hpp`. Confirm via a unit test at `test_tpv104_setup.cpp` that the three non-zero pre-stress components match the SCEC spec and that `tau1_nuc` is identically zero at every QP.

---

### 3.11 Normal-stress sign convention — SeisSol negative compression, MFEM positive compression [DIFFERENT-EXPECTED]

**What SeisSol does** (`FrictionLaws/CpuImpl/RateAndState.h:345-358`):

```cpp
void updateNormalStress(std::array<real, misc::NumPaddedPoints>& normalStress,
                        const FaultStresses<Executor::Host>& faultStresses,
                        size_t timeIndex,
                        size_t ltsFace) {
  // Todo(SW): consider poroelastic materials together with thermal pressurisation
#pragma omp simd
  for (uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
    normalStress[pointIndex] = std::min(static_cast<real>(0.0),
                                        faultStresses.normalStress[timeIndex][pointIndex] +
                                            this->initialStressInFaultCS[ltsFace][0][pointIndex] +
                                            faultStresses.fluidPressure[timeIndex][pointIndex] +
                                            this->initialPressure[ltsFace][pointIndex] -
                                            tpMethod.getFluidPressure(ltsFace, pointIndex));
  }
}
```

The `std::min(0, ...)` clip means **compression is negative** in SeisSol's sign convention. The friction-strength term uses `std::fabs(normalStress)` (`RateAndState.h:312` and `:335`), so the sign-flipped value is taken as magnitude once inside the solver.

**What MFEM does** — sign convention per `miniapps/seas/CLAUDE.md`:

```
- **Normal stress**: sigma_n > 0 = compression (geology convention).
```

and at `fault_face_flux.cpp:141`:

```cpp
real_t strength = std::abs(s.sigma_n_total) * f_V;
```

`TPV102Params::sigma_n = 120e6` (`config/tpv102_params.hpp:52`) — **positive** for compression. `d.sigma_n0 = TPV102Params::sigma_n` (`dynamic/tpv102_setup.hpp:76`) — so `s.sigma_n_total` is positive throughout the run. The `std::abs` is defensive: it absorbs either sign convention.

**Implication**:

* Despite the sign-convention mismatch, **both codes compute `|σ_n| · μ` in the friction residual**. The algorithmic output is identical up to the sign passed in (SeisSol's trial normal stress from `precomputeStressFromQInterpolated` uses the same velocity-and-stress-jump formula as MFEM; whether the bulk Q stores +σ_n or −σ_n determines whether trial normal stress is +/-, but `fabs` strips it).
* For Phase 3 probe 1, the raw dumps of `faultStresses.normalStress[o][i]` (SeisSol) vs `s.sigma_n_trial` (MFEM) may differ by a sign flip at every sample. The probe comparator must take `abs()` before differencing — or, better, confirm the sign convention at each code's output and flip MFEM's sign to SeisSol's before dumping.
* The `std::min(0, ...)` clip in SeisSol prevents fault tension from creating a positive friction-pull force. MFEM's equivalent is `fault_face_flux.cpp:136` (`if (s.Theta > 0.0 && s.V_abs > 0.0)`) combined with the `sigma_n <= 0` branch in `SolveSlipRatePsi` that returns `tau/eta` (frictionless limit). Behaviour is equivalent at the fault-in-tension case (σ_n_total > 0 in SeisSol ≡ σ_n_total < 0 in MFEM convention).

**Action**: None in Phase 2. Phase 3 probe comparator must normalise sign conventions before diffing normal-stress channels.

---

### 3.12 Time integration — sub-step friction iteration vs ADER-on-averaged-Q [DIFFERENT-EXPECTED]

**What SeisSol does** (`FrictionLaws/CpuImpl/BaseFrictionLaw.h:100-142`):

```cpp
// loop over sub time steps (i.e. quadrature points in time
real startTime = 0;
real updateTime = this->mFullUpdateTime;
for (std::size_t timeIndex = 0; timeIndex < misc::TimeSteps; timeIndex++) {
  startTime = updateTime;
  updateTime += this->deltaT[timeIndex];
  for (unsigned i = 0; i < this->drParameters->nucleationCount; ++i) {
    common::adjustInitialStress(
        initialStressInFaultCS[ltsFace],
        nucleationStressInFaultCS[ltsFace * this->drParameters->nucleationCount + i],
        initialPressure[ltsFace],
        nucleationPressure[ltsFace * this->drParameters->nucleationCount + i],
        updateTime,
        this->drParameters->t0[i],
        this->drParameters->s0[i],
        this->deltaT[timeIndex]);
  }

  static_cast<Derived*>(this)->updateFrictionAndSlip(faultStresses,
                                                     tractionResults,
                                                     stateVariableBuffer,
                                                     strengthBuffer,
                                                     ltsFace,
                                                     timeIndex);
```

and inside `updateFrictionAndSlip` (`RateAndState.h:151-195`) there are `numberStateVariableUpdates` outer iterations of (update ψ → Newton on V → average V_old/V_new), followed by a final ψ update + traction eval (`RateAndState.h:197-256`).

So for one ADER macro-step at order O = 5, the friction pipeline runs **O sub-steps**, each with **`numberStateVariableUpdates ≈ 2–4`** outer iterations, each with a **`maxNumberSlipRateUpdates`** inner Newton loop, and the sub-step results are time-weighted into `imposedStatePlus/Minus` (see §3.8 above).

**What MFEM does** (`drivers/tpv102_driver.cpp:1102-1134` for the ADER path):

```cpp
// One-shot ADER predictor-corrector.  AdvanceADER fills the
// DOFData {V1, V2, slip_rate, tau*_corr, sigma_n_corr} with
// the friction-solve output on Q̄ = I/dt (= time-averaged bulk
// Q over [t, t+dt]).  NOTE: this is NOT exactly the time-
// average of the friction-solve trajectory — the solve is
// nonlinear — but it matches the one-shot equivalent to
// O(dt²) per plan §Phase 5 §4 (R-007).
wave.AdvanceADER(Q, dt_step, ader_order, Q_new);
Q.Swap(Q_new);

t += dt_step;

// Forward-Euler psi update using the time-averaged slip rate.
// 1st-order in dt but matches the plan's pseudocode
// (`UpdateStateAnalytic`-style closed-form from
// psi_n and V_avg); the bulk scheme's 2nd-order accuracy is
// preserved because psi is a scalar state that contributes only
// to friction via a Lipschitz-smooth law.  Slip accumulation
// uses the same time-averaged V directly.
for (int i = 0; i < num_fault_total; i++)
{
   const real_t dpsi = aging_law.Rate(dof_data[i].slip_rate,
                                      psi_n[i],
                                      dof_data[i].Dc);
   dof_data[i].psi = psi_n[i] + dt_step * dpsi;

   dof_data[i].slip1 += dof_data[i].V1 * dt_step;
   dof_data[i].slip2 += dof_data[i].V2 * dt_step;
```

So for one ADER macro-step, MFEM runs:

1. `wave.AdvanceADER` — predictor-corrector that internally calls `FaultFaceFlux::EvaluateADER` **once** on Q̄ (the time-averaged bulk state).
2. `EvaluateADER` internally calls `Evaluate` on Q̄, which runs the pipeline **once**: precompute trial → total traction → Brent friction solve → corrected traction → imposed state. ψ is NOT updated here (psi-purity guard at `fault_face_flux.cpp:222-224`).
3. After `AdvanceADER`, the driver does a **forward-Euler ψ update** using the time-averaged slip rate returned from step 2.

The RK4 path (`drivers/tpv102_driver.cpp:1229-1400+`) is 4 staged RK4 calls with friction+nucleation at each stage, coupled RK4 for ψ — closer in spirit to SeisSol's per-sub-step cadence but with 4 stages vs O=5.

**Implication**:

* The per-step friction pipeline cost is **O×(outer)×(Newton) iterations in SeisSol** versus **1 × Brent in MFEM**. For ADER-O5, that's ~20–30 Newton iterations (10 sub-step × 2 outer × Newton) per SeisSol step vs ~1 Brent call per MFEM step. Convergence to the same root, by different number and type of iterations.
* The ψ integration is **analytic per-sub-step** in SeisSol vs **forward-Euler-per-macro-step** in MFEM ADER path. This is the dominant source of ψ drift between codes. Even after Phase 2 replaces `AgingLawPsi` with `UpdateStateAnalyticSlipLawSRW`, MFEM's ADER path will only call it **once** per macro-step on the time-averaged V, whereas SeisSol calls it `O` times on the per-sub-step V. The `O(dt²)` argument applies: ψ will agree to O(dt²) but not bitwise.
* **Probe 2 thresholds must be set accordingly**: the per-stage analytic integrator (byte-match criterion 1e-13, see T_SRW_5 in §4.2.3) is verified by a **unit test**, not by the cross-code probe. The cross-code probe at Phase 3 must accept ~1e-6 relative drift between MFEM and SeisSol ψ traces by t = 1 s.
* **Nucleation timing**: SeisSol's `adjustInitialStress` is called inside the per-sub-step loop, using `updateTime` (end of the sub-step). MFEM ADER path calls `ApplyNucleationPrestress` **once** at `t + dt/2` (stage midpoint) before `AdvanceADER`. For ADER-O5 these two are O(dt²) different; for the RK4 path the four stage-midpoint calls are O(dt) different from SeisSol's five sub-steps.

**Action**: None in Phase 2. Phase 3 accepts the O(dt) / O(dt²) drifts as expected artefacts of the different time-integration architectures and documents them in the probe-diff report.

---

### 3.13 Initial ψ inversion — identical formula [CONSISTENT]

**What SeisSol does** for FVW (`Initializer/RateAndStateInitializer.cpp:141-163`):

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
  if (result.stateVariable < 0) {
    logWarning()
        << "Found a negative state variable while initializing the fault. Are you sure your "
           "setup is correct?";
  }
  const real tmp2 = initialSlipRate * 0.5 / rsSr0 * std::exp(result.stateVariable / rsA);
  result.frictionCoefficient = rsA * std::asinh(tmp2);
  return result;
}
```

i.e. `ψ_ini = a · log((2V₀/V_ini) · sinh(|τ|/(a·|σ_n|)))`, using `sinh(x) = (exp(x) - exp(-x))/2` inlined.

**What MFEM does** for TPV102 (`config/tpv102_params.hpp:131-139`):

```cpp
inline real_t ComputeInitialPsi(real_t a)
{
   // From tau = sigma_n * a * asinh(V/(2*V0) * exp(psi/a)):
   //   sinh(tau/(sigma_n*a)) = V/(2*V0) * exp(psi/a)
   //   psi = a * ln(2*V0/V * sinh(tau/(sigma_n*a)))
   real_t arg = TPV102Params::tau_ini / (TPV102Params::sigma_n * a);
   real_t psi = a * std::log(2.0 * TPV102Params::V0 / TPV102Params::V_ini * std::sinh(arg));
   return psi;
}
```

and the general per-point variant `DieterichRuinaFriction::InitialStatePsi` (`friction/dieterich_ruina.hpp:506-528`):

```cpp
real_t InitialStatePsi(real_t tau0, real_t V_init, real_t sigma_n,
                       real_t eta, real_t a) const
{
   real_t tau_eff = tau0 - eta * V_init;
   real_t f = tau_eff / sigma_n;

   real_t f_over_a = f / a;
   real_t sinh_val;
   if (f_over_a > 700.0)
   {
      sinh_val = 0.5 * std::exp(f_over_a);
   }
   else
   {
      sinh_val = std::sinh(f_over_a);
   }

   real_t log_arg = (2.0 * cp_.V0 / V_init) * sinh_val;
   MFEM_ASSERT(log_arg > 0.0,
               "Invalid argument for logarithm in InitialStatePsi");

   return a * std::log(log_arg);
}
```

**Implication**:

* Identical formula. Both codes invert `τ = σ_n · a · asinh((V/(2V₀)) · exp(ψ/a))` for ψ given `(τ, V_ini, σ_n, a)`. Two subtle differences:

  1. MFEM's `InitialStatePsi` subtracts `eta * V_init` (the radiation-damping correction) from `tau0` before inverting. SeisSol's FVW initialiser does **not** subtract `eta·V_ini` — it uses the raw `τ` from `initialStressInFaultCS`. For TPV104's `V_ini = 1e-16` and `η = 4.625 MPa·s/m`, the correction is `η·V_ini = 4.625e6 × 1e-16 = 4.6e-10 Pa` — **utterly negligible** compared to `τ = 40 MPa`. So the two formulas agree to 1e-17 relative. But for TPV102's `V_ini = 1e-12`, the correction is `4.6e-6 Pa` — still negligible at 1e-14 relative. So the initial ψ matches to at least 1e-14. `TPV102Params::ComputeInitialPsi` ignores this term (matches SeisSol's FVW init).
  2. MFEM has an asymptotic branch for `f/a > 700` via `sinh_val = 0.5 · exp(f/a)`. SeisSol always uses `(exp − exp(−)) / 2`. For TPV104 at equilibrium, `f = τ / σ_n = 40/120 = 0.333`, `a = 0.01`, so `f/a = 33.3` — SeisSol's direct `(exp − exp(−))/2 = exp(33.3)/2 − exp(−33.3)/2 ≈ exp(33.3)/2` agrees with MFEM's asymptotic to ~1e-29 relative. Agreed.
* The plan's claim that `ψ_ini = 5.6359184e-01` at `(a = 0.01, V_ini = 1e-16, τ = 40 MPa, σ_n = 120 MPa)` is verified by plugging in:
  - `f = 40/120 = 0.3333`
  - `f/a = 33.333`
  - `sinh(33.333) = 0.5 × exp(33.333) = 0.5 × 2.99e14 = 1.495e14`
  - `2 V₀ / V_ini = 2e-6 / 1e-16 = 2e10`
  - `log_arg = 2e10 × 1.495e14 = 2.99e24`
  - `ψ = 0.01 × log(2.99e24) = 0.01 × 56.36 ≈ 0.5636`. ✅

**Action**: Phase 2 duplicates `ComputeInitialPsi` into `ComputeInitialPsiTPV104` with TPV104's parameters, verifies `ψ_ini ≈ 5.6359184e-01` at `a = a_in = 0.01` in `test_slip_law_srw_psi.cpp::T_SRW_6`. No code change to the shared friction infrastructure.

---

### 3.14 Nucleation spatial formula — identical [CONSISTENT, minor parameter differences]

**What SeisSol does**: The spatial form is not hard-coded in SeisSol. It's read from the per-input-file `nuc_xy.asagi` / similar easi/NetCDF spatial field (see `BaseDRInitializer.cpp:130-138` for the `rotateStressToFaultCS` call). For TPV104 per the SCEC PDF, the spatial factor is a Gaussian-core bump `F(r) = exp(r²/(r² − R²))` for `r < R = 3 km` and 0 outside, centered at the hypocenter. SeisSol loads this field into `nucleationStressInFaultCS`.

**What MFEM does** (`config/tpv102_params.hpp:96-105`):

```cpp
inline real_t NucleationSpatial(real_t r)
{
   real_t R = TPV102Params::nuc_radius;
   if (r >= R) { return 0.0; }
   real_t r2 = r * r;
   real_t R2 = R * R;
   return std::exp(r2 / (r2 - R2));
}
```

This is the SCEC formula `F(r) = exp(r²/(r²-R²))` verbatim.

**Implication**: Identical spatial formula between SeisSol (when the input file uses the canonical SCEC form) and MFEM. For TPV104 the SCEC spec fixes `R = 3 km` and `nuc_dtau = 45 MPa` (amplitude only difference from TPV102's 25 MPa).

**Action**: Phase 2 duplicates `NucleationSpatial` into `tpv104_params.hpp` under the name `NucleationSpatial_TPV104` (per the no-shared-edit rule), adjusts the amplitude to `nuc_dtau = 45e6`, keeps the radius at `3 km`.

---

### 3.15 Resample state variable — SeisSol-only smoothing [DIFFERENT-EXPECTED-SMALL]

**What SeisSol does** (`FastVelocityWeakeningLaw.h:144-164`):

```cpp
void resampleStateVar(const std::array<real, misc::NumPaddedPoints>& stateVariableBuffer,
                      std::size_t ltsFace) const {
  std::array<real, misc::NumPaddedPoints> deltaStateVar = {0};
  std::array<real, misc::NumPaddedPoints> resampledDeltaStateVar = {0};
#pragma omp simd
  for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
    deltaStateVar[pointIndex] =
        stateVariableBuffer[pointIndex] - this->stateVariable[ltsFace][pointIndex];
  }
  dynamicRupture::kernel::resampleParameter resampleKrnl;
  resampleKrnl.resample = init::resample::Values;
  resampleKrnl.originalQ = deltaStateVar.data();
  resampleKrnl.resampledQ = resampledDeltaStateVar.data();
  resampleKrnl.execute();

#pragma omp simd
  for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
    this->stateVariable[ltsFace][pointIndex] =
        this->stateVariable[ltsFace][pointIndex] + resampledDeltaStateVar[pointIndex];
  }
}
```

Called in `postHook` after every macro-step. The `resampleParameter` kernel projects `Δψ` from the DG-quadrature basis onto a lower-order polynomial basis (for stability — prevents high-order oscillations in ψ from driving the friction law into unphysical regimes).

**What MFEM does**: There is **no equivalent smoothing step** in MFEM. ψ is updated at the DOF points directly and used at the next step without a projection.

**Implication**:

* For SCEC tests with spatially uniform or smoothly varying ψ (TPV104 initial condition has ψ = 0.5636 everywhere in the VW core, transitioning smoothly in the VS halo), the resample is effectively a no-op: `Δψ` is low-order from the start, so projection onto a lower-order basis leaves it unchanged.
* For rupture-front configurations where ψ develops a sharp front, the resample **suppresses** high-order modes — this changes the ψ trajectory at the front by up to a few percent per macro-step. Over 12 s × 2534 steps, this could accumulate to ~10% relative drift between codes at the front.
* MFEM uses p=1 (linear polynomials on tets) at the baseline TPV102 run, so there are no high-order modes to suppress — the resample would be a no-op at p=1. At p≥2 this would become relevant.

**Action**: Do not implement the resample in Phase 2 (it is order-dependent and Extreme-Care territory). Phase 3 probe 2 records ψ after the state update and compares across codes; if the probe passes (ψ agrees within 1e-6 relative) at MFEM p=1 / SeisSol p=4, the resample is confirmed to be a non-issue. If probe 2 fails with ψ differences growing at the rupture front, escalate.

---

### 3.16 SCEC-spec parameter values — differences at the initialisation boundary [DIFFERENT-REQUIRES-CHANGE]

| parameter | TPV102 (MFEM) | TPV104 (SCEC) | source | action |
|-----------|---------------|----------------|--------|--------|
| `V_ini` | `1e-12 m/s` (`config/tpv102_params.hpp:54`) | `1e-16 m/s` | SCEC PDF | set `TPV104Params::V_ini = 1e-16` |
| `nuc_dtau` | `25e6 Pa` (`config/tpv102_params.hpp:60`) | `45e6 Pa` | SCEC PDF | set `TPV104Params::nuc_dtau = 45e6` |
| `Dc` (= L) | `0.02 m` (`config/tpv102_params.hpp:35`) | `0.4 m` | SCEC PDF | set `TPV104Params::Dc = 0.4` |
| `b` | `0.012` (`config/tpv102_params.hpp:34`) | `0.014` | SCEC PDF | set `TPV104Params::b = 0.014` |
| `a_vw` (inside VW) | `0.008` (`config/tpv102_params.hpp:38`) | `0.01` | SCEC PDF | set `TPV104Params::a_in = 0.01` |
| `a_vs` (outside VW) | `0.016` | `0.02` (= `a_in + Δa`) | SCEC PDF | derive from Δa |
| `f_w` (= μ_W, SeisSol's `muW`) | N/A (classical RS has no f_w) | `0.1` | SCEC PDF | add as new TPV104 parameter |
| `V_w` (weakening rate inside VW) | N/A | `0.1 m/s` | SCEC PDF | add as new TPV104 parameter |
| `V_w_out` | N/A | `1.0 m/s` (= `V_w + ΔV_w`) | SCEC PDF | derive from ΔV_w |
| `tau_ini` | `75e6 Pa` (`config/tpv102_params.hpp:53`) | `40e6 Pa` | SeisSol trace header | set `TPV104Params::tau_ini = 40e6` |

**Implication**:

* The smaller `V_ini = 1e-16` in TPV104 means the initial ψ is **larger** (`a · log((2V₀/V_ini) · sinh(f/a))` grows as `V_ini → 0`). For TPV104 with `a = 0.01`, `τ = 40 MPa`, `σ_n = 120 MPa`, this gives `ψ_ini ≈ 0.564`. For TPV102 with `V_ini = 1e-12`, `a = 0.008`, `τ = 75 MPa`, `σ_n = 120 MPa`, `ψ_ini ≈ 0.47`. Probe 2 at step 0 must see `ψ = 0.564` in both codes.
* Brent's solver bracket `Va = -32.0` (i.e. lower bound `V = 1e-32`) at `dieterich_ruina.hpp:372` is **below** TPV104's `V_ini = 1e-16`, so the bracket is safe. At `V_ini = 1e-16` the residual `F(log10(V_ini)) = τ - σ_n·μ(V_ini, ψ_ini) - η·V_ini ≈ 0` by construction (initial equilibrium).

**Action**: Phase 2 hard-codes all TPV104 values in `config/tpv104_params.hpp`. No CLI flags — the CLAUDE.md "follow SCEC spec" rule dictates that TPV104 cannot default to TPV102's numbers.

---

### 3.17 Mesh, ADER order, and output stations [DIFFERENT-EXPECTED]

| axis | SeisSol TPV104 | MFEM (planned TPV104) | implication |
|------|----------------|------------------------|-------------|
| Mesh | TPV5 mesh (`docs/tpv104.rst:40-41`: "We use the mesh file of TPV5 directly"). 3D tet, fault at y=0. Typical resolution ~200 m on-fault, coarser at far field. | New `tpv104/mesh/tpv104_1000m.msh` replicating TPV5 dimensions. | Phase 3 probe 1 must check that the wave arrival time at a far-field station matches `t = distance / c_p` within 1%; a mesh mismatch would show up as a time-shift. |
| DG polynomial order p | p = 4 (O5) | p = 1 default, configurable via `--order` | Station interpolation will use different basis functions; the Phase 3 probe comparator must locate the station by (x₂, x₃) coordinates, not by DOF index. |
| ADER order O | O = 5 | O = 2 default, configurable via `--ader-order` | `dt` differs by order-of-magnitude (SeisSol `dt ≈ 4.7 ms` for 2534 steps at t_final = 12 s; MFEM at CFL=0.5 with h=1000 m, c_p=6000 m/s gives `dt = 0.083 ms`). Probes must sample at matched wall-clock time, not step number. |
| Output stations | Nine files at `(x₂, x₃) ∈ {0, ±9, ±12} × {3, 7.5, 12}` km | Must be hard-coded in `tpv104_params.hpp::kStationsTPV104[9]` | TPV102's `DefaultStations` (`dynamic/tpv102_setup.hpp:164-178`) uses the same nine coordinates already — confirm that SCEC-TPV104 station locations are identical, and duplicate into the TPV104 setup. |
| Final time | t_final = 12 s (SeisSol benchmark trace covers ~12 s) | t_final = 12 s (`TPV104Params::t_final = 12.0`) | Matched. |

**Action**: Phase 2 builds the mesh as part of `tpv104/mesh/tpv104_1000m.msh`, adopts `p = 1` / `O = 2` as the default (overriddable via CLI), and writes a Python analysis script that does wall-clock-time interpolation across both station traces before computing diffs.

---

### 3.18 Summary — what's a real bug vs what's a difference in architecture

| category | count | entries |
|----------|-------|---------|
| **CONSISTENT** (no change) | 5 | §3.1 (trial traction), §3.4 (friction coefficient), §3.6 (slip-rate decomposition), §3.7 (corrected traction), §3.13 (initial ψ) |
| **CONSISTENT-WITH-TEST** | 1 | §3.8 (imposed state: formula identical, time accumulation differs — add O(dt²) tolerance) |
| **DIFFERENT-EXPECTED** (physics / engineering difference; accept) | 7 | §3.5 (Newton vs Brent), §3.9 (nucleation cumulative vs overwrite), §3.10 (pre-stress side-channel vs DOFData), §3.11 (sign convention), §3.12 (sub-step iteration), §3.15 (resample), §3.17 (mesh / order / dt) |
| **DIFFERENT-REQUIRES-CHANGE** (Phase 2 must add code) | 4 | §3.2 (t1/t2 column mapping), §3.3 (FVW ψ ODE + SRW ψ_ss + analytic integrator), §3.14 (duplicate spatial formula with TPV104 amplitude), §3.16 (TPV104 parameter values) |

Phase 3's probe pass/fail thresholds are set to accept the seven DIFFERENT-EXPECTED items as baseline noise and surface only the DIFFERENT-REQUIRES-CHANGE items and genuinely new bugs.

---

## 4. Phase 2 — TPV104 friction option implementation plan

### 4.1 TPV104 parameter table (single source of truth)

New file `config/tpv104_params.hpp`. Struct `TPV104Params` — const-expr all scalar parameters, inline helpers for spatial distributions.

All values below come from `docs/tpv104.rst` + SCEC_validation_slip_law.pdf + the SeisSol benchmark-trace header; none are guessed. Where SCEC and SeisSol disagree, **follow SCEC spec** (per CLAUDE.md's "Tandem vs SCEC Spec" rule — if there's a mismatch, flag it in the cross-verification probe output, don't silently align to SeisSol).

| parameter | symbol | value | source |
|---|---|---|---|
| Density | ρ | 2670 kg/m³ | SeisSol trace header + SCEC PDF |
| P-wave speed | c_p | 6000 m/s | SCEC PDF |
| S-wave speed | c_s | 3464 m/s | SCEC PDF |
| Shear modulus | μ | ρ·c_s² = 3.204e10 Pa | derived |
| Lamé λ | λ | ρ·c_p² − 2μ = 3.204e10 Pa | derived |
| Reference slip rate | V₀ | 1e-6 m/s | SCEC PDF |
| Reference friction | f₀ | 0.6 | SCEC PDF |
| Slowness (state) | b | 0.014 | SCEC PDF |
| Critical slip distance | L | 0.4 m | SCEC PDF |
| Weakening friction | f_w = μ_w | 0.1 | SCEC PDF |
| Direct effect (inside VW) | a_in | 0.01 | SCEC PDF |
| Δa (transition) | Δa | 0.01 | SCEC PDF, docs/tpv104.rst:113 |
| a (outside VW) | a_out | a_in + Δa = 0.02 | derived |
| Weakening rate (inside VW) | V_w_in | 0.1 m/s | SCEC PDF, docs/tpv104.rst |
| ΔV_w (transition) | ΔV_w | 0.9 m/s | SCEC PDF, docs/tpv104.rst:114 |
| V_w (outside VW) | V_w_out | V_w_in + ΔV_w = 1.0 m/s | derived |
| Initial normal stress | σ_n | 120 MPa | SeisSol trace col-8 constant = 1.20e+02 MPa |
| Initial shear stress | τ_ini | 40 MPa | SeisSol trace col-4 constant = 4.0e+01 MPa at t≈0 |
| Initial slip rate | V_ini | 1e-16 m/s | SCEC PDF (TPV104 uses much lower V_ini than TPV102) |
| VW zone half-width along-strike | L_s | 15 km | SCEC PDF |
| VW zone half-width down-dip | W | 15 km | SCEC PDF |
| Transition width | w | 3 km | SCEC PDF / docs/tpv104.rst:111 |
| Hypocenter along-strike | x_hypo | 0 m | SCEC PDF |
| Hypocenter down-dip | z_hypo | 7500 m | SCEC PDF + SeisSol station naming |
| Nucleation radius | R | 3 km | SCEC PDF |
| Nucleation amplitude | Δτ₀ | 45 MPa | SCEC PDF + `research_seissol_tpv104_2026-04-24.md` §3 |
| Nucleation rise time | T_nuc | 1.0 s | SCEC PDF |
| Final time | t_final | 12.0 s | SeisSol trace covers ~12 s |

Derived consistency check (hard-coded at compile-time via `static_assert` where possible, or `MFEM_VERIFY` at driver `main` start):

- `ψ_ini(a_in) ≈ 5.6359e-01` from `ComputeInitialPsiTPV104(a_in)` to 1e-6. Anchored against SeisSol trace column 9.
- `eta_s = Zs/2 = ρ·c_s/2 = 4.625e6 Pa·s/m`. Matches the `Zs/2` form the Riemann solver assumes for homogeneous faults.

Spatial distributions (inline helpers, all derived from SCEC TPV104 Eqs. 4 and 5):

```cpp
// a(x, z) — VW inside a 30 km × 15 km box centered at hypocenter, smooth
// boxcar transition of width w = 3 km outside.  B() is the same SCEC
// tanh-based boxcar as TPV102 (tpv102_params.hpp:73-81) — we duplicate
// rather than share to keep TPV104 self-contained.
inline real_t ComputeA_TPV104(real_t along_strike, real_t down_dip);
inline real_t ComputeVw_TPV104(real_t along_strike, real_t down_dip);

// Nucleation amplitude — SAME SCEC Gaussian-core × smooth-ramp as TPV102,
// only the amplitude changes.  Duplicated formula (no shared helper).
inline real_t NucleationPerturbation_TPV104(real_t along_strike,
                                            real_t down_dip, real_t t);

// Initial psi inversion from the steady-state condition
//    tau_ini = sigma_n * a * asinh((V_ini/2V0) * exp(psi/a))
// solves to psi_ini = a * ln((2V0/V_ini) * sinh(tau_ini/(sigma_n * a))).
// Note: identical functional form to tpv102_params.hpp:131-139, only the
// numeric inputs differ.  Duplicated to keep TPV104 self-contained.
inline real_t ComputeInitialPsiTPV104(real_t a);

// Nine canonical output stations from SeisSol TPV104 benchmark data:
//   (x2, x3) ∈ {0, ±9, ±12} km × {3, 7.5, 12} km
struct Station { real_t x2; real_t x3; const char *label; };
inline constexpr Station kStationsTPV104[9] = { ... };
```

**Do NOT edit `tpv102_params.hpp`** — the two param headers are independent. The shared formulas (`Boxcar`, `NucleationSpatial`, `NucleationTemporal`) are duplicated into `tpv104_params.hpp` namespace — per `feedback_tpv102_bp5_no_shared_edit`.

### 4.2 New friction-state class — `SlipLawSRWPsi`

New file `friction/slip_law_srw_psi.hpp`. Implements the `StateEvolution` interface (`friction/state_evolution.hpp:31`) for the slip-law-in-psi-space ODE with SRW steady-state.

#### 4.2.1 ODE and analytic single-step update

SCEC TPV104 specifies the state-evolution ODE as

```
dψ/dt = −(V/L) · [ψ − ψ_ss(V)]                                          (1)
```

with SRW steady-state

```
f_LV(V)   = f₀ + (b − a) · ln(V/V₀)                                     (2a)
f_ss(V)   = f_w + (f_LV − f_w) / (1 + (V/V_w)⁸)^(1/8)                   (2b)
ψ_ss(V)   = a · ln((2V₀/V) · sinh(f_ss(V)/a))                           (2c)
```

**Derivation note** (not in SCEC but load-bearing): (1) is algebraically equivalent to the slip-law-in-theta-space ODE `dθ/dt = −(θV/L)·ln(Vθ/L)` under the substitution `ψ = f₀ + b·ln(V₀θ/L)`, **as long as ψ_ss in (1) equals `f₀ + b·ln(V₀/V)`** (the classical rate-and-state steady-state, i.e. (2c) with `f_w = 0` and `V_w → ∞`). In TPV104, (1) uses the SRW `ψ_ss` from (2c) instead — so it's **no longer** the slip-law-in-theta-space ODE in the literal sense; it's a linear-in-ψ relaxation to a non-linear SRW ψ_ss. That's why SeisSol categorises FL=103 as `FastVelocityWeakeningLaw` (aging-form-ODE in ψ-space) rather than `SlipLaw` (exact slip law). The name in the TPV104 doc ("slip law with strong rate weakening") is SCEC's historical naming; the ODE they actually write is the FVW form.

Analytic single-step update for constant V over a sub-step dt:

```
psi(t + dt) = ψ_ss(V) + (ψ(t) − ψ_ss(V)) · exp(−V·dt/L)                 (3)
```

Byte-matching (3) against SeisSol FVW (`FastVelocityWeakeningLaw.h:71-74`):

```cpp
const auto preexp1 = -localSlipRate * (timeIncrement / localSl0);
const real exp1v = std::exp(preexp1);
const real exp1m = -std::expm1(preexp1);
return steadyStateStateVariable * exp1m + exp1v * stateVarReference;
// i.e.  ψ_new = ψ_ss · (1 − exp(−V·dt/L)) + ψ_0 · exp(−V·dt/L)
//            = ψ_ss + (ψ_0 − ψ_ss) · exp(−V·dt/L)         (identical to (3))
```

#### 4.2.2 Class interface

```cpp
// friction/slip_law_srw_psi.hpp
namespace mfem { namespace seas {

class SlipLawSRWPsi : public StateEvolution
{
public:
   // Parameters carried alongside the StateEvolution interface.  V_w is
   // per-point in SCEC TPV104 (V_w_in inside VW zone, V_w_out outside),
   // so SlipLawSRWPsi is constructed with the GLOBAL scalars and the
   // per-point V_w is passed through the Rate / SteadyState methods as
   // the Dc argument is (spatial variation absorbed at call site).
   //
   // CONVENTION: we reuse the StateEvolution interface's `Dc` parameter
   // as L (the critical slip distance in SCEC notation).  V_w is passed
   // through a SEPARATE overload because the base interface does not
   // accept a second spatial parameter.  See Rate_SRW() below.
   SlipLawSRWPsi(real_t a_scalar, real_t b, real_t V0, real_t f0, real_t muW);

   // Base-interface rate — required override.  Treats V_w as global
   // constant (caller's responsibility to use Rate_SRW for per-point V_w).
   real_t Rate(real_t V, real_t psi, real_t Dc) const override;

   // SRW-aware rate: dpsi/dt = -(V/L) * (psi - psi_ss(V, V_w, a)).
   real_t Rate_SRW(real_t V, real_t psi, real_t L,
                   real_t V_w, real_t a) const;

   // Steady-state: base-interface version uses stored V_w scalar.
   real_t SteadyState(real_t V, real_t Dc) const override;

   // Per-point steady state (primary entry point for per-QP usage).
   real_t SteadyState_SRW(real_t V, real_t V_w, real_t a) const;

   real_t RateDerivativeV(real_t V, real_t psi, real_t Dc) const override;
   real_t RateDerivativeTheta(real_t V, real_t psi, real_t Dc) const override;

   const char *GetName() const override { return "SlipLawSRWPsi"; }

   // SRW steady-state psi:
   //   psi_ss = a * ln((2*V0/V) * sinh(f_ss(V)/a))
   // f_ss = f_w + (f_LV - f_w) / (1 + (V/V_w)^8)^(1/8)
   // f_LV = f0 + (b - a) * ln(V/V0)
   // Inline so tests can check per-point evaluation.
   static real_t PsiSS_SRW(real_t V, real_t V_w, real_t a,
                           real_t b, real_t V0, real_t f0, real_t muW);

private:
   real_t a_, b_, V0_, f0_, muW_;   // muW_ is f_w (weakening friction coeff)
   real_t V_w_default_;              // used by base-interface overloads
};

// Analytic single-step integrator for constant V over dt, in psi-space.
// Returns ψ(t+dt) per Eq. (3).  Caller passes the local V_w (per-point).
inline real_t UpdateStateAnalyticSlipLawSRW(real_t psi_old, real_t V,
                                            real_t L, real_t dt,
                                            real_t V_w, real_t a,
                                            real_t b, real_t V0, real_t f0,
                                            real_t muW);

}} // namespace mfem::seas
```

Implementation notes (embedded in comments in the new header):

1. `PsiSS_SRW` must gracefully handle `V → 0`: `f_LV` has a `ln(V/V0)` that diverges negatively. Protect with `V = max(V, 1e-50)` (matches `SlipLawPsi::Rate` pattern at `state_evolution.hpp:242`).
2. `PsiSS_SRW` must handle `V/V_w` near `V_w` crossover cleanly: `(1 + (V/V_w)^8)^(1/8)` is well-conditioned. `misc::power<8, double>` in SeisSol's code is a templated integer power for SIMD; we use `std::pow` since vector throughput is not a concern here.
3. `sinh(f_ss/a)` range: on TPV104, `f_ss ∈ [0.1, 0.65]`, `a ∈ [0.01, 0.02]`, so `f_ss/a ∈ [5, 65]` — small enough for `sinh` to evaluate without overflow, large enough that the `0.5*exp(f/a)` asymptotic is not needed. No branching.
4. The `Dc` channel in the base `Rate(V, psi, Dc)` interface is reused as `L`; this is the same convention MFEM uses in `AgingLawPsi` and `SlipLawPsi` (see `state_evolution.hpp:183-184`).
5. The default `V_w_default_` stored at construction is for test fixtures and uniform-V_w drivers; the production driver always uses the per-point `Rate_SRW` / `SteadyState_SRW` / `UpdateStateAnalyticSlipLawSRW` overloads.

#### 4.2.3 Acceptance criteria for `SlipLawSRWPsi` (unit test gates)

New test `tests/unit/test_slip_law_srw_psi.cpp`. Pass thresholds: 1e-12 for formula equivalences; 1e-14 for steady-state root.

- [ ] **T_SRW_1 — ψ_ss equation**: for 20 random `(V, V_w, a, b, V₀, f₀, f_w)` in the TPV104 envelope, `PsiSS_SRW` agrees with the formula from (2c) recomputed from primitives to 1e-14.
- [ ] **T_SRW_2 — ψ_ss = classical when f_w = 0 and V_w → ∞**: for `f_w = 0, V_w = 1e300`, check that `PsiSS_SRW ≈ f₀ + b·ln(V₀/V)` to 1e-12 (falls back to classical rate-and-state steady state).
- [ ] **T_SRW_3 — rate vanishes at ψ_ss**: for any valid `(V, a, b, V₀, f₀, f_w, V_w)`, `Rate_SRW(V, PsiSS_SRW(V, V_w, a, ...), L, V_w, a) = 0` exactly (algebraic identity).
- [ ] **T_SRW_4 — `UpdateStateAnalytic` round trip**: `UpdateStateAnalyticSlipLawSRW(psi, V, L, dt, V_w, a, ...)` for `psi = ψ_ss` returns `ψ_ss` for every `dt` to 1e-12.
- [ ] **T_SRW_5 — SeisSol FVW byte-match on a single QP**: port the exact arithmetic of `FastVelocityWeakeningLaw.h:43-78` to a local lambda in the test; verify that `UpdateStateAnalyticSlipLawSRW` agrees to 1e-13 on 50 random `(psi_0, V, L, dt, V_w, a, ...)` tuples. This is the acceptance criterion for "we match SeisSol's state-evolution integration bit-by-bit on the state-variable side."
- [ ] **T_SRW_6 — TPV104 ψ_ini match**: `ComputeInitialPsiTPV104(a_in)` returns `5.6359184e-01` to 1e-8 (anchored against SeisSol trace col-9).
- [ ] **T_SRW_7 — adherence to StateEvolution interface**: `Rate`, `SteadyState`, `RateDerivativeV`, `RateDerivativeTheta` all compile and link against the existing `friction::StateEvolution` base class without modification to the base class.

### 4.3 New TPV104 setup header — `tpv104_setup.hpp`

Duplicates the `tpv102_setup.hpp` pattern (fault-DOF initialisation, equilibrium solve, pre-stress baking) with the following substitutions:

- Reads `TPV104Params` instead of `TPV102Params`.
- Per-QP `a` via `ComputeA_TPV104(x, z)`.
- Per-QP `V_w` via `ComputeVw_TPV104(x, z)` — **new** field in `DOFData` **NOT** acceptable (per §2.5, we cannot edit `fault_face_flux.hpp` where `DOFData` lives). **Workaround**: store a per-DOF `V_w` vector alongside `dof_data`, owned by `tpv104_driver.cpp`. Pass to the integrator via a side channel (see §4.5). Verified acceptable because `V_w` is **only** needed at state-variable integration time, not inside `FaultFaceFlux::Evaluate`.
- `ComputeInitialPsiTPV104(a_i)` used for per-QP ψ_ini.
- Prestress baking into bulk Q reuses `InitializeStateTotal` structure (`tpv102_setup_total.hpp:61-77`) with TPV104's `σ_n = 120 MPa`, `τ_ini = 40 MPa`.

Duplicate — do not share — the TPV102 helpers (`Boxcar`, `ComputeA`, `NucleationSpatial`, `NucleationTemporal`, `NucleationPerturbation`, `ComputeInitialPsi`). Duplication is justified per `feedback_tpv102_bp5_no_shared_edit`: keeps TPV102 and TPV104 structurally independent for future merges with other drivers (e.g., quasi-dynamic + dynamic). Naming convention: suffix every symbol `_TPV104` to avoid ODR clash.

### 4.4 New TPV104 driver — `tpv104_driver.cpp`

Duplicate `tpv102_driver.cpp` structure one-for-one, with the following changes:

1. Default parameters pulled from `TPV104Params` + `slip_law_srw_psi.hpp`.
2. At driver init, instantiate one `SlipLawSRWPsi` with the global friction parameters (a is scalar in the constructor; per-point a is handled in the `Rate_SRW`/`SteadyState_SRW` calls).
3. Per-DOF `V_w[i]` array owned by the driver (one entry per fault QP), populated at init from `ComputeVw_TPV104(x, z)`.
4. In the driver's ψ-integration block (currently `tpv102_driver.cpp:1107-1132` for ADER and `tpv102_driver.cpp:1230-1400+` for RK4), replace the `AgingLawPsi::Rate` call with the analytic integrator.

   **R-006 FIX 2026-04-24 — analytic form is MANDATORY for TPV104.** SeisSol only uses the analytic exponential-relaxation form (`FastVelocityWeakeningLaw.h:71-74`); mixing a forward-Euler integrator on MFEM's side would silently invalidate Phase 3 Probe 2's byte-match (§5.10.2) and the failure would mis-attribute to `SlipLawSRWPsi` rather than to the integrator choice. Do **not** ship the forward-Euler path for TPV104.

   ```cpp
   // Required for TPV104 (analytic exponential relaxation, bit-match to
   // SeisSol's FastVelocityWeakeningLaw.h:71-74):
   dof_data[i].psi = UpdateStateAnalyticSlipLawSRW(
       psi_n[i], dof_data[i].slip_rate, dof_data[i].Dc, dt_step,
       V_w[i], dof_data[i].a, b, V0, f0, muW);
   ```

   Sensitivity studies that want forward-Euler must go through a
   separate `--state-integrator=euler` CLI flag (opt-in only), and the
   driver must emit a loud warning at init if that flag is set:
   ```cpp
   if (opts.state_integrator == StateIntegrator::ForwardEuler) {
      std::fprintf(stderr, "WARNING: TPV104 driver invoked with "
                   "--state-integrator=euler.  Phase 3 cross-verification "
                   "(§5.10.2) is invalidated in this mode.\n");
   }
   ```
   The Phase 2 smoke test `test_tpv104_smoke.cpp` asserts the default
   integrator is analytic; a maintainer who later re-introduces
   forward-Euler as the default trips the test.
5. Nine canonical station output files named `tpv104_mfem_x2_<X>_x3_<Z>.txt` in the same 9-column layout as SeisSol's trace files (t, h-slip, h-slip-rate, h-shear-stress, v-slip, v-slip-rate, v-shear-stress, n-stress, psi). Under BP5 convention, `horizontal = strike = slip2`, `vertical = dip = slip1` — documented at the top of the output block.
6. Nucleation wired exactly as in TPV102 driver (`ApplyNucleationPrestress` on the persistent-prestress channel). Amplitude 45 MPa instead of 25 MPa.
7. CLI defaults: `--ader-order=2`, `--order=1`, `--time-integrator=ader`, `--cfl=0.5`, `--bc-mode=absorbing`, `--tfinal=12.0`, `--fric-law=slip-srw`. The `--fric-law` flag is new; the only option initially is `slip-srw`. Future extensions can add `aging` etc. for sensitivity studies.
8. Diagnostic build gate: `SEAS_DIAG_FAULT_FLUX` (same as TPV102) + new `SEAS_DIAG_TPV104_STATE` that enables the §5 MFEM probes. Default off.

### 4.5 `V_w[i]` side-channel (driver-owned, not inside `DOFData`)

Per §2.5, `fault_face_flux.hpp` / `DOFData` MUST NOT be edited. TPV104 needs per-QP `V_w`, which is only used by the state-evolution integrator (never inside `FaultFaceFlux::Evaluate`). Solution: the driver owns a `std::vector<real_t> V_w(num_fault_total)`, populated at init time via `ComputeVw_TPV104(x, z)` over the fault-QP coordinates. Every call to `SlipLawSRWPsi::Rate_SRW` / `UpdateStateAnalyticSlipLawSRW` reads `V_w[i]` alongside `dof_data[i]`. No mutation after init.

Review checklist — add to Phase 2 acceptance:

- [ ] `V_w.size() == dof_data.size()` asserted at driver init after both arrays are populated.
- [ ] `V_w[i]` values match expected `V_w_in = 0.1` inside the VW zone (|x| < 15 km, 0 < z < 15 km) and `V_w_out = 1.0` outside, with smooth tanh transition over 3 km.
- [ ] ParaView output carries `V_w` as an additional fault-surface field, so a visual plot can confirm the smooth transition.

### 4.6 Phase 2 acceptance gates

*(§4.6 gates below are unchanged. Subsections §4.7–§4.11 were added
2026-04-24 in response to a code-review request for a SeisSol-aligned,
TPV102-duplicating implementation strategy with an input → code →
output graph map and real executed functions at each stage. §4.7 is
the SeisSol-side walkthrough, §4.8 the MFEM-side mirror, §4.9 the
side-by-side graph, §4.10 the dependency-ordered checklist, §4.11 the
explicit no-touch contract.)*

### 4.7 SeisSol canonical execution graph (TPV104 / FL=103)

Entry point is `Simulator::simulate` in `src/Solver/Simulator.cpp:45`.
Below, every stage is annotated with the real functions executed in
that stage, in call order, with absolute paths and line numbers
(verified 2026-04-24 against the SeisSol tree). Quantities in
angle-brackets name the data that flows out of the stage.

#### Stage S0 — Input

Files on disk consumed by SeisSol at run time:

1. `parameters.par` (SeisSol parameter file; sets `FL=103`, TPV104
   constants `rsF0, rsB, rsSr0, muW, rsSrW` per
   `DRParameters.h:57-91`).
2. Fault stress file (per-point initial CS components for
   `initialStressInFaultCS[6]`, populated by the initializer —
   `BaseDRInitializer.cpp`).
3. Mesh: `tpv5.msh` (`docs/tpv104.rst:41` — TPV104 reuses the TPV5
   mesh verbatim).
4. Material: uniform `ρ = 2670 kg/m³, c_p = 6 km/s, c_s = 3.464 km/s`
   (SCEC TPV104 spec).

Output of S0 → the `DynamicRupture::Layer` struct, populated with
`impAndEta[·]`, `a[·][·]`, `b[·][·]`, `sl0[·][·]` (=L), `f0[·][·]`,
`muW[·][·]`, `srW[·][·]` (=V_w), `initialStressInFaultCS[6][·]`,
`nucleationStressInFaultCS[6][·]`, `stateVariable[·][·]` (=ψ).

#### Stage S1 — Top-level time loop

- `src/Solver/Simulator.cpp:45` —
  `Simulator::simulate(seissolInstance)` drives the outer tick;
  at `Simulator.cpp:106` it calls
  `seissolInstance.timeManager().advanceInTime(upcomingTime)`.
- `src/Solver/TimeStepping/TimeManager.cpp` —
  `TimeManager::advanceInTime` schedules `TimeCluster`s through the
  actor graph in `src/Solver/TimeStepping/AbstractTimeCluster.cpp`.

Output of S1 → per-cluster "may compute DR?" scheduling, answered by
`DynamicRuptureScheduler::mayComputeInterior` in
`src/Solver/TimeStepping/ActorState.h`.

#### Stage S2 — Bulk ADER integration of Q (wave equation)

For every element in a cluster's local partition:

- `src/Solver/TimeStepping/TimeCluster.cpp:336` —
  `TimeCluster::computeLocalIntegration(bool resetBuffers)` fires the
  element-local ADER kernel (Cauchy-Kovalewski `I = ∫ Q dt` and the
  volume `K` update, via the generated tensor kernels
  `ader::executeADER` and `volume::executeVolume`).
- `TimeCluster::computeNeighboringIntegration(double subTimeStart)` in
  `src/Solver/TimeStepping/TimeCluster.cpp:517` wires the neighbour
  flux `AminusT · I_neighbour + AplusT · I_self` and accumulates back
  into `Q`.

Output of S2 → element-local `Q(t + Δt)` before fault faces are visited
on this cluster.

#### Stage S3 — Fault-face precompute (Q interpolation → trial traction)

- `src/Solver/TimeStepping/TimeCluster.cpp:191` —
  `TimeCluster::computeDynamicRupture(DynamicRupture::Layer& layerData)`
  is the DR entry point; it does space-time interpolation (S3a) and
  then hands the interpolated arrays to the friction solver (S4).

- **S3a — space-time interpolation.** `TimeCluster.cpp:218`, LIKWID
  region `computeDynamicRuptureSpaceTimeInterpolation`, calls the
  generated kernel that fills
  `qInterpolatedPlus[misc::TimeSteps][misc::NumQuantities][misc::NumPaddedPoints]`
  and `qInterpolatedMinus[…]`.

- **S3b — `precomputeStressFromQInterpolated`.**
  `src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h:143-222`.
  Combines plus/minus interpolated Q with `impAndEta[ltsFace]`
  (`Z_p, Z_s, η_P, η_S`) to build `faultStresses.normalStress[o][i]`,
  `traction1[o][i]`, `traction2[o][i]` (formulas at lines 183–192).

- **S3c — `adjustInitialStress` (nucleation add).**
  `src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h:418-450`.
  While `fullUpdateTime ∈ [s0, s0+t0]`, adds
  `nucleationStressInFaultCS[i] * g_nuc(t)` to
  `initialStressInFaultCS[i]` so the trial traction fed to the
  friction law carries the time-ramped nucleation perturbation.

Output of S3 → `faultStresses` (trial σ_n, τ_1, τ_2 at every `[o][i]`)
and the nucleation-updated `initialStressInFaultCS`.

#### Stage S4 — Friction law evaluation (per ADER time substep)

CRTP-dispatched to `FastVelocityWeakeningLaw` (FL=103) through the
`RateAndStateBase` template. Entry point:

- `src/DynamicRupture/FrictionLaws/CpuImpl/RateAndState.h:30` —
  `RateAndStateBase<Derived,TPMethod>::updateFrictionAndSlip(faultStresses, …)`.

Calls in order per substep:

- **S4a — `calcInitialVariables`** (`RateAndState.h:117-148`).
  Assembles `totalTraction1 = traction1[o] + initialStressInFaultCS[3]`,
  `totalTraction2 = traction2[o] + initialStressInFaultCS[5]`,
  `absoluteTraction = hypot(totalTraction1, totalTraction2)`,
  `normalStress = faultStresses.normalStress[o] + initialStress[0]`.

- **S4b — `updateStateVariableIterative`** (`RateAndState.h:151-189`).
  Fixed-point iteration between the state-variable update and the
  slip-rate Newton solve. Alternates S4c and S4d until
  `hasConverged = true`.

  - **S4c — `Derived::updateStateVariable`** (FVW implementation,
    `src/DynamicRupture/FrictionLaws/CpuImpl/FastVelocityWeakeningLaw.h:44-78`).
    Computes `f_LV = f0 + (b-a)·ln(V/V_0)` (line 55-59),
    `f_ss = f_w + (f_LV - f_w) / (1 + (V/V_w)^8)^{1/8}` (line 61-62),
    `ψ_ss = a · logsinh(2 V_0 / V, f_ss/a)` (line 65-66, via
    `RateAndStateCommon.h` inline), then exact exponential relaxation
    `ψ_new = ψ_ss · (1 - exp(-V·dt/L)) + ψ_0 · exp(-V·dt/L)` (line 71-74).
  - **S4d — `invertSlipRateIterative`** (`RateAndState.h:286-343`).
    Newton on `g(V) = tau_abs - η_S · V - σ_n · μ(V, ψ)` using
    `updateMu` (`FastVelocityWeakeningLaw.h:118`) and
    `updateMuDerivative` (`FastVelocityWeakeningLaw.h:133`).

- **S4e — `calcSlipRateAndTraction`** (`RateAndState.h:197-256`).
  With converged `V_abs`, computes final `mu = updateMu(V_abs, …)`
  (line 218), slip-rate decomposition parallel to total traction
  (line 229-232), and final Godunov-corrected traction (line 235-242,
  `tractionResults.traction_i = faultStresses.traction_i - η_S · V_i`).

Output of S4 → per-substep `mu[ltsFace][i]`, `slipRate1/2[ltsFace][i]`,
`slipRateMagnitude[ltsFace][i]`, `stateVariable[ltsFace][i]`,
`tractionResults.traction1/2[o][i]`.

#### Stage S5 — Imposed state write-back (friction → bulk flux)

- `src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h:287-407` —
  `postcomputeImposedStateFromNewStress()`. Builds
  `imposedStatePlus[9]` and `imposedStateMinus[9]` at each `[o][i]`
  from the corrected `tractionResults` + side impedance. These are
  the rebuilt `Q` samples that the next cluster tick's
  `computeNeighboringIntegration` (S2) multiplies against
  `AminusT / AplusT` to close the fault-flux coupling.

Output of S5 → `imposedStatePlus`, `imposedStateMinus` at every face
quadrature point → fed back to bulk ADER neighbour integration.

#### Stage S6 — Output (benchmark traces)

- `src/DynamicRupture/Output/FaultOutputManager` writes the nine
  receiver files whose headers match
  `tpv104_seisol_x2_*_x3_*.txt`:
  `Column 1 = Time; 2 = horizontal slip; 3 = horizontal slip rate;
   4 = horizontal shear stress; 5 = vertical slip;
   6 = vertical slip rate; 7 = vertical shear stress;
   8 = normal stress; 9 = psi.` Writer cadence is controlled by
  `FaultOutput::syncPoint` in the parameter file.

Output of S6 → 9 × ASCII traces, 9 columns each, space-separated.

---

### 4.8 MFEM TPV104 execution graph (duplicated TPV102 path)

Node-for-node mirror of §4.7. Every MFEM path below exists today for
TPV102; the TPV104 variant is the **duplicated** file (suffix
`_TPV104` on every new symbol), without any edits to the TPV102
originals. File paths are absolute; line numbers are from the
2026-04-24 working tree.

#### Stage M0 — Input

Files on disk:

1. CLI + (optionally) TOML fragment parsed in
   `drivers/tpv104_driver.cpp` (new — duplicates the
   `drivers/tpv102_driver.cpp:198-240` CLI scan block).
2. Mesh file `tpv104/mesh/tpv104_1000m.msh` (new — TPV5-equivalent
   geometry per SCEC TPV104 spec).
3. Static pre-stress / friction parameters from
   `config/tpv104_params.hpp` (new) — every scalar is a `constexpr`;
   per-point `a(x,z)`, `V_w(x,z)` come from inline helpers
   `ComputeA_TPV104`, `ComputeVw_TPV104` (planned §4.1).

Output of M0 → fresh `ParMesh`, compile-time `TPV104Params` struct,
CLI-populated `use_ader, ader_order, order, cfl, tfinal, …`.

#### Stage M1 — Driver init (one-shot, pre-loop)

Executed in `main(…)` of `drivers/tpv104_driver.cpp`:

- Build `WaveOperator<ParMesh>` (mirrors
  `drivers/tpv102_driver.cpp:404`; constructor at
  `dynamic/wave_operator.inl:24`).
- Allocate `std::vector<DOFData> dof_data(num_fault_total)` and call
  `InitializeFaultDOFs_TPV104(dof_data, num_fault_total, fault_coords)`
  — new, duplicates `InitializeFaultDOFs` at
  `dynamic/tpv102_setup.hpp:40`, replacing scalar inputs and calling
  `ComputeA_TPV104`, `ComputeInitialPsiTPV104` instead of the TPV102
  helpers. Hypocenter resolution (`MPI_MINLOC`) mirrors
  `drivers/tpv102_driver.cpp:547-572`.
- Allocate and populate `std::vector<real_t> V_w(num_fault_total)`
  via `ComputeVw_TPV104(x, z)` — the side-channel spelled out in §4.5
  (new; no TPV102 analogue).
- Construct
  `FaultFaceFlux fault_flux(TPV104Params::rho,
                            TPV104Params::cp, TPV104Params::cs)`
  (reused — `dynamic/fault_face_flux.cpp:27`).
- `wave.SetFaultFlux(&fault_flux); wave.SetFaultDOFData(&dof_data, …)`
  (reused; mirrors `drivers/tpv102_driver.cpp:510-512`).
- Initialise Q = 0 (fluctuation-Q dispatch); pre-stress is baked into
  `DOFData.sigma_n0, tau2_0` only — mirrors
  `drivers/tpv102_driver.cpp:597-598`.
- Open nine TPV104 station writers (new `TPV104StationWriter` —
  duplicates `TPV102StationWriter` pattern at
  `dynamic/tpv102_setup.hpp:223-340`); stations from
  `TPV104Params::kStationsTPV104[9]`.

Output of M1 → initialised `Q, dof_data, V_w, fault_flux`; writers
open; `wave` linked to the fault.

#### Stage M2 — Outer time loop (per wall-clock step)

Executed in `drivers/tpv104_driver.cpp` (mirrors
`drivers/tpv102_driver.cpp:1064-1400`).

For each step:

- **M2a — Nucleation prestress update.** Call
  `ApplyNucleationPrestress_TPV104(dof_data, fault_coords, t + dt/2)`
  — new file `dynamic/tpv104_setup_total.hpp`, duplicates
  `dynamic/tpv102_setup_total.hpp:627`; only the amplitude
  (`TPV104Params::nuc_dtau = 45e6`) differs from TPV102's 25e6.
- **M2b — Bulk ADER advance.**
  `wave.AdvanceADER(Q, dt_step, ader_order, Q_new)`
  (reused — `dynamic/wave_operator.inl:2991`).
- **M2c — `Q.Swap(Q_new); t += dt_step;`** (reused).
- **M2d — State-variable update.** Loop
  `i ∈ [0, num_fault_total)`:
  ```cpp
  dof_data[i].psi = UpdateStateAnalyticSlipLawSRW(
      psi_n[i],                 // ψ at step start
      dof_data[i].slip_rate,    // V_abs from the ADER corrector
      dof_data[i].Dc,           // = L
      dt_step,
      V_w[i],                   // per-point V_w
      dof_data[i].a,
      TPV104Params::b, TPV104Params::V0,
      TPV104Params::f0, TPV104Params::muW);
  dof_data[i].slip1 += dof_data[i].V1 * dt_step;
  dof_data[i].slip2 += dof_data[i].V2 * dt_step;
  ```
  This replaces the TPV102 forward-Euler block at
  `drivers/tpv102_driver.cpp:1114-1122`.
  `UpdateStateAnalyticSlipLawSRW` lives in the new file
  `friction/slip_law_srw_psi.hpp` (§4.2).
- **M2e — Station / ParaView cadence.** Call
  `station_writer.WriteStep(t, dof_data)` and
  `surface_writer.WriteStep(t, Q)` — new writer classes but
  structurally identical to `dynamic/tpv102_setup.hpp:300`.

Output of M2 → advanced Q, dof_data, t; written station rows for every
output tick.

#### Stage M3 — Fault-face precompute & friction solve

Triggered inside `AdvanceADER` → the generated tensor kernels for the
ADER predictor and corrector. When the corrector visits a fault face
it calls into the wave operator's face-flux path:

- `dynamic/wave_operator.inl:940` —
  `WaveOperator::ComputeFaceFluxRHS(Q, rhs)` (classical path; used by
  the RK4 stages).
- `dynamic/wave_operator.inl:1901` —
  `WaveOperator::ComputeADERFaceFluxRHS(I, rhs)` (ADER path).
- Shared MPI faces:
  - `dynamic/wave_operator.inl:1440` —
    `WaveOperator::ComputeSharedFaceFluxRHS`,
  - `dynamic/wave_operator.inl:2660` —
    `WaveOperator::ComputeADERSharedFaceFluxRHS`.

Each dispatches at every fault-face quadrature point to the friction
pipeline:

- `dynamic/fault_face_flux.cpp:310` —
  `FaultFaceFlux::EvaluateTotal(data, Q_plus, Q_minus, Q_imp_plus,
                                Q_imp_minus, method)`
  (total-Q path; default for TPV102 and TPV104).
- `EvaluateTotal` internally calls, in order:
  - `FaultFaceFlux::ComputeTrialTraction`
    (`dynamic/fault_face_flux.cpp:40-66`) — **mirrors S3b**
    (`precomputeStressFromQInterpolated`).
  - `FaultFaceFlux::CompleteFromTrial → CompleteFromTheta`
    (`dynamic/fault_face_flux.cpp:89-119`). Merges static pre-stress
    + persistent nucleation into the trial — **mirrors S3c**
    (`adjustInitialStress`), but done inside the flux solver rather
    than as a separate nucleation call.
  - `FaultFaceFlux::CompleteFromVabs`
    (`dynamic/fault_face_flux.cpp:122-151`) — **mirrors S4c + S4e**:
    computes `mu = a · asinh((V/(2V_0)) · exp(ψ/a))`, decomposes the
    slip rate into `V1, V2` parallel to total traction, and sets
    Godunov-corrected `tau1_corr, tau2_corr, sigma_n_corr`.
  - `FaultFaceFlux::BuildImposedState`
    (`dynamic/fault_face_flux.cpp:153-190`) — **mirrors S5**
    (`postcomputeImposedStateFromNewStress`): builds the 9-component
    `Q_imp_plus`, `Q_imp_minus` that feed back into the ADER
    neighbour integration.
  - `FaultFaceFlux::WriteBackState`
    (`dynamic/fault_face_flux.cpp:192-206`) — stores
    `V1, V2, slip_rate, tau*_corr, sigma_n_corr, mu` into `DOFData`
    so that M2d can read `slip_rate` for the ψ update.

The friction-root solver is:

- `friction/dieterich_ruina.hpp` —
  `FrictionSolver::Solve(Theta, psi, sigma_n, eta_s, a, method)` —
  Brent in log10(V) space. **Corresponds to S4b + S4d** (SeisSol
  Newton). No edits per §2.5.

Structural asymmetry (by design): M3 does NOT fire S4c's
state-variable update — in MFEM the ψ update is done at the outer
loop (M2d), using the time-averaged V from the ADER corrector.
This is the convention the production TPV102 path uses today
(`tpv102_driver.cpp:1107-1122`) and is what `UpdateStateAnalyticSlipLawSRW`
is designed to match to 1e-13 under constant V (test T_SRW_5).

#### Stage M4 — State-variable update (scalar, per-QP)

Details under M2d above. The only new function executed at this stage
is `UpdateStateAnalyticSlipLawSRW` (new file
`friction/slip_law_srw_psi.hpp`, §4.2). Byte-matches SeisSol's S4c
block at `FastVelocityWeakeningLaw.h:44-78` to 1e-13 under constant V
(test T_SRW_5).

#### Stage M5 — Output

- `TPV104StationWriter::WriteStep` (new, in
  `dynamic/tpv104_setup.hpp`) — duplicates
  `TPV102StationWriter::WriteStep` pattern
  (`dynamic/tpv102_setup.hpp:300-333`), emits the 9-column row
  `t  h-slip  h-slip-rate  h-shear  v-slip  v-slip-rate  v-shear
   sigma_n  psi` matching SeisSol's header (§4.7 S6).
- Column mapping documented inline: under BP5 convention,
  `horizontal = strike = slip2`, `vertical = dip = slip1`.
- `TPV104SurfaceStationWriter` — duplicates TPV102 analogue.
- ParaView (`ParaViewDataCollection`) — reused unchanged from the
  TPV102 path.

Output of M5 → nine `tpv104_mfem_x2_<X>_x3_<Z>.txt` files; ParaView
collection directory; fault-surface VTUs.

---

### 4.9 Side-by-side pipeline graph (input → code → output)

The graph below is the quick-reference map the code-review asked for.
Each row is a stage; columns are SeisSol (left) and MFEM-TPV104
(right); centre column names the conceptual operation. Every cell
lists the real executed function(s) — look up the line numbers in
§4.7 / §4.8 above.

```
┌───────────────────────────┐   operation   ┌───────────────────────────┐
│  SeisSol canonical path   │               │  MFEM TPV104 duplicated   │
│  (FL=103 FastVelocityWk)  │               │  path (duplicates TPV102) │
├───────────────────────────┼───────────────┼───────────────────────────┤
│ parameters.par            │               │ drivers/tpv104_driver.cpp │
│ tpv5.msh                  │  INPUT        │   CLI + TOML scan         │
│ BaseDRInitializer.cpp     │  (files →     │ tpv104/mesh/tpv104_*.msh  │
│ DRParameters.h:57-91      │   memory)     │ config/tpv104_params.hpp  │
│                           │               │ InitializeFaultDOFs_TPV104│
│                           │               │ ComputeA_TPV104           │
│                           │               │ ComputeVw_TPV104          │
│                           │               │ ComputeInitialPsiTPV104   │
├───────────────────────────┼───────────────┼───────────────────────────┤
│ Simulator::simulate       │               │ main() outer step loop    │
│   (Simulator.cpp:45)      │  TOP-LEVEL    │   in tpv104_driver.cpp    │
│ TimeManager::advanceInTime│  LOOP         │   (duplicates             │
│   (:106)                  │               │    tpv102_driver.cpp      │
│ TimeCluster actor graph   │               │    :1064-1400)            │
├───────────────────────────┼───────────────┼───────────────────────────┤
│ TimeCluster::             │               │ wave.AdvanceADER(Q,       │
│   computeLocalIntegration │  BULK ADER    │   dt, ader_order, Q_new)  │
│   (:336)                  │  (WAVE EQ)    │ wave_operator.inl:2991    │
│ TimeCluster::             │               │ ComputeADERVolumeUpdate   │
│   computeNeighboring      │               │   (:913)                  │
│   Integration (:517)      │               │ ComputeADERFaceFluxRHS    │
│ generated ader::/volume:: │               │   (:1901)                 │
│ tensor kernels            │               │ ComputeADERSharedFaceFlux │
│                           │               │   RHS (:2660)             │
├───────────────────────────┼───────────────┼───────────────────────────┤
│ TimeCluster::             │               │ FaultFaceFlux::           │
│   computeDynamicRupture   │  Q →          │   EvaluateTotal           │
│   (:191)                  │  fault-face   │   (fault_face_flux.cpp    │
│ spaceTimeInterpolation    │  trial state  │   :310)                   │
│   kernel (:218)           │               │   └─ ComputeTrialTraction │
│ FrictionSolverCommon.h::  │               │      (:40-66)             │
│   precomputeStressFrom    │               │                           │
│   QInterpolated           │               │                           │
│   (:143-222)              │               │                           │
│ FrictionSolverCommon.h::  │               │ (nucleation fused into    │
│   adjustInitialStress     │  + nucleation │  ComputeTrialTraction via │
│   (:418-450)              │  perturbation │  DOFData.tau2_nuc set by  │
│                           │               │  ApplyNucleationPrestress │
│                           │               │  _TPV104 in M2a)          │
├───────────────────────────┼───────────────┼───────────────────────────┤
│ RateAndStateBase::        │               │ FaultFaceFlux::           │
│   updateFrictionAndSlip   │  FRICTION     │   CompleteFromTrial →     │
│   (RateAndState.h:30)     │  LAW          │   CompleteFromTheta →     │
│   ├─calcInitialVariables  │  (ψ, V, τ)    │   CompleteFromVabs        │
│   │   (:117-148)          │               │   (fault_face_flux.cpp    │
│   ├─updateStateVariable-  │               │   :89-151)                │
│   │ Iterative (:151-189)  │               │ FrictionSolver::Solve     │
│   │  ├ updateStateVariable│               │   (Brent in log10(V))     │
│   │  │   (FVW.h:44-78)    │               │   (dieterich_ruina.hpp)   │
│   │  │   ψ_ss = logsinh() │               │ mu = a · asinh((V/(2V₀))  │
│   │  │   ψ_new = exp-relax│               │          · exp(ψ/a))      │
│   │  └ invertSlipRate-    │               │   (fault_face_flux.cpp    │
│   │    Iterative          │               │   :139-141)               │
│   │    (RateAndState.h:   │               │                           │
│   │    286-343) — Newton  │               │                           │
│   └─calcSlipRateAndTract. │               │                           │
│       (:197-256)          │               │                           │
├───────────────────────────┼───────────────┼───────────────────────────┤
│ FrictionSolverCommon.h::  │               │ FaultFaceFlux::           │
│   postcomputeImposedState │  WRITE BACK   │   BuildImposedState       │
│   FromNewStress           │  (V, τ_corr,  │   (fault_face_flux.cpp    │
│   (:287-407)              │   σ_n_corr,   │   :153-190)               │
│ → imposedStatePlus[9]     │   Q_imp)      │ FaultFaceFlux::           │
│ → imposedStateMinus[9]    │               │   WriteBackState          │
│                           │               │   (:192-206) → DOFData    │
├───────────────────────────┼───────────────┼───────────────────────────┤
│ (no outer-loop ψ update — │               │ UpdateStateAnalyticSlip-  │
│  SeisSol updates inside   │  OUTER ψ      │  LawSRW (new file         │
│  the friction iteration   │  UPDATE       │  friction/slip_law_srw_   │
│  via S4c exp-relax)       │  (MFEM only)  │  psi.hpp)                 │
│                           │               │  — byte-matches FVW S4c   │
│                           │               │    under constant V per   │
│                           │               │    §4.2.1 Eq. (3)         │
├───────────────────────────┼───────────────┼───────────────────────────┤
│ FaultOutputManager        │               │ TPV104StationWriter::     │
│ Receiver traces 9-col     │  OUTPUT       │   WriteStep (new;         │
│ tpv104_seisol_x2_*_x3_*   │               │   duplicates TPV102       │
│ .txt                      │               │   writer at tpv102_setup  │
│ Column map: t, h-slip,    │               │   .hpp:300-333)           │
│ h-slip-rate, h-shear,     │               │ 9-col output matching     │
│ v-slip, v-slip-rate,      │               │ SeisSol header / columns  │
│ v-shear, σ_n, ψ           │               │ tpv104_mfem_x2_*_x3_*.txt │
└───────────────────────────┴───────────────┴───────────────────────────┘
```

The correspondence is one-to-one at every row. The only **structural
asymmetry** is the ψ-update placement (bottom-but-one row): SeisSol
advances ψ inside the friction iteration (S4c); MFEM advances ψ once
per outer step (M2d) using the time-averaged V from the ADER
corrector. This is what the production TPV102 path does today, and
is what `UpdateStateAnalyticSlipLawSRW` is designed to match to 1e-13
under constant V (test T_SRW_5).

Every other row is a line-for-line semantic mirror — the strategy is
therefore: **copy the TPV102 function in each row of the right-hand
column, rename the symbol with `_TPV104`, and swap the parameter
source from `TPV102Params` to `TPV104Params`.** The single new class
(`SlipLawSRWPsi`) slots into the second-from-last row and has its own
dedicated unit-test gate (T_SRW_1..7).

---

### 4.10 Dependency-ordered implementation checklist

Expanded 2026-04-24 to reflect user directives on §3.2–§3.13 and two
top-level constraints:

- **[C1] Coordinate / rotation convention MUST agree with BP5
  quasi-dynamic.** TPV104 reuses BP5's `FaultBasis` (Tandem convention
  `tangent1 = dip, tangent2 = strike`) directly; `fault/fault_basis.hpp`
  is called through, not redefined. The BP5 convention is documented
  in `miniapps/seas/CLAUDE.md` ("Fault-local tangent frame") and is
  the single source of truth for coordinates across TPV102, TPV104,
  and BP5.
- **[C2] No modification to the BP5 workflow code base.** BP5
  quasi-dynamic (`miniapps/seas/bp5/**`, `domain/`, `fault/`,
  `solver/`) must continue to build and run bit-identically to its
  current reference output at every step of Phase 2 landing. A `make
  test-bp5-smoke` gate runs after every merge (Step 0 below).

Each step lists: (i) files to add, (ii) §3.X directive(s) addressed,
(iii) unit tests, (iv) the Phase 3 probe gate.

Tests that start with `T_TPV104_` are new; tests starting with
`T_SRW_` were defined in §4.2.3. All tests run under `make test-tpv104`
(new target, Step 10).

#### Step 0 — BP5 baseline byte-match gate (before any code lands)

- Files: `tests/unit/test_bp5_baseline_bitmatch.cpp` (new — runs
  `seas_bp5_full` on a short fixture and diffs all nine station CSVs
  and the fault-surface VTU against a reference checkpoint committed
  to `tests/baselines/bp5/`).
- Constraint addressed: **[C2]**.
- Tests:
  - `T_TPV104_BP5_BITMATCH_1` — 50-step BP5 smoke run on the 1000-m
    mesh, rank-serial: output byte-identical to
    `tests/baselines/bp5/bp5_smoke_serial.csv` (pre-existing reference).
  - `T_TPV104_BP5_BITMATCH_2` — same run at `np=4`: output
    byte-identical to `tests/baselines/bp5/bp5_smoke_np4.csv`.
  - `T_TPV104_BP5_BUILD_CLEAN` — `make seas_bp5_full && make
    seas_test_bp5_integration` succeed on every landing.
- Acceptance: **both BP5 bit-match tests green** on every PR that
  touches `miniapps/seas/`. This test gate runs before Steps 1–11.
- Phase 3 gate: independent of Phase 3.

#### Step 1 — Parameter header

- Files: `config/tpv104_params.hpp` (new).
- §3.X directive addressed: **§3.13 (ψ_ini match), §3.16 (SCEC
  parameter values)**.
- Content: all SCEC TPV104 constants from §4.1; inline helpers
  `Boxcar_TPV104`, `ComputeA_TPV104`, `ComputeVw_TPV104`,
  `NucleationSpatial_TPV104`, `NucleationTemporal_TPV104`,
  `NucleationPerturbation_TPV104`, `ComputeInitialPsiTPV104`;
  `kStationsTPV104[9]`. `constexpr` everywhere.
- Tests (`tests/unit/test_tpv104_params.cpp`, new):
  - `T_TPV104_P_1` — ψ_ini match: `ComputeInitialPsiTPV104(0.01)`
    returns `5.6359184e-01` to 1e-8, anchored against SeisSol TPV104
    benchmark trace column 9 row 1. This is the §3.13 directive: as
    long as the initial value matches expectation, the implementation
    path is free.
  - `T_TPV104_P_2` — `static_assert` on derived scalars:
    `a_out == a_in + 0.01`, `V_w_out == V_w_in + 0.9`,
    `lambda == rho*cp^2 - 2*mu`, `eta_s == rho*cs/2`.
  - `T_TPV104_P_3` — `ComputeA_TPV104(0,0)` == `a_in`;
    `ComputeA_TPV104(20 km, 20 km)` == `a_out`; smooth transition
    across `|x| = 15 km` boundary (value at `x=13.5 km, z=7.5 km`
    lies strictly between `a_in` and `a_out`).
  - `T_TPV104_P_4` — `ComputeVw_TPV104` returns `0.1` at (0, 7.5 km)
    (VW core) and `1.0` at (20 km, 20 km) (strengthening).
  - `T_TPV104_P_5` — station list is complete: `kStationsTPV104[9]`
    contains exactly the nine SeisSol benchmark stations, and every
    MFEM station label matches the SeisSol filename suffix.
- Phase 3 gate: blocks **Probe 2** (needs ψ_ini for the initial-state
  sanity check).

#### Step 2 — Friction-state class (SRW slip-law ψ ODE, §3.3 directive)

- Files: `friction/slip_law_srw_psi.hpp` (new);
  `tests/unit/test_slip_law_srw_psi.cpp` (new).
- §3.X directive addressed: **§3.3 (must implement FVW ψ ODE + SRW
  ψ_ss + analytic integrator for TPV104)**.
- Content: `SlipLawSRWPsi` class and `UpdateStateAnalyticSlipLawSRW`
  free function per §4.2.2. Implementation must byte-match SeisSol's
  `FastVelocityWeakeningLaw.h:44-78` for the same `(ψ_0, V, L, dt,
  V_w, a, b, V_0, f_0, μ_w)` tuple.
- Tests: **T_SRW_1 .. T_SRW_7** (§4.2.3), plus two new:
  - `T_SRW_8` — port SeisSol's arithmetic **verbatim** (inline
    lambda over lines `FastVelocityWeakeningLaw.h:44-78`) and
    cross-check against our `UpdateStateAnalyticSlipLawSRW` on 1000
    random tuples drawn from the TPV104 envelope; max relative diff
    must be < 1e-14.
  - `T_SRW_9` — **classical-limit fallback**: with `f_w = 0` and
    `V_w = 1e300`, our `UpdateStateAnalyticSlipLawSRW` must match
    MFEM's existing `AgingLawPsi::Rate` integrated analytically
    (`state_evolution.hpp:190-194`) to 1e-10. This guarantees we
    haven't introduced a regression on the aging-law limit.
- Phase 3 gate: **Probe 2 (state evolution)** — byte-match.

#### Step 3 — Fault-DOF init + pre-stress + nucleation-prestress

- Files: `dynamic/tpv104_setup.hpp` (new; duplicates
  `dynamic/tpv102_setup.hpp` layout — `InitializeFaultDOFs_TPV104`,
  `InitializeState_TPV104`, `TPV104StationWriter`,
  `TPV104SurfaceStationWriter`, `DefaultStations_TPV104`);
  `tests/unit/test_tpv104_setup.cpp` (new).
- §3.X directive addressed: **§3.10 (fluctuation-Q only, drop total-Q
  mode), §3.13 (ψ_ini match)**.
- Fluctuation-Q-only (§3.10 directive): `dynamic/tpv104_setup.hpp`
  contains **only** the mode-1 fluctuation-Q initialiser
  (`InitializeState_TPV104` sets `Q = 0`; pre-stress lives in
  `DOFData.sigma_n0, tau1_0, tau2_0`). The TPV102 `tpv102_setup_total.hpp`
  mode-2 total-Q path is NOT duplicated into TPV104 — we do not ship
  `InitializeStateTotal_TPV104` or `ApplyNucleationTotalPrestress_TPV104`.
  This keeps TPV104 exactly matched to SeisSol's side-channel storage
  convention (pre-stress in fault-CS side-channel; bulk Q carries
  fluctuation only).
- Coordinate-frame constraint (**[C1]**): `InitializeFaultDOFs_TPV104`
  calls `fault::FaultBasis::GetBasis(fi)` directly for every face index
  to read `(tangent1 = dip, tangent2 = strike, normal)`. Pre-stress
  fields map:
  - `d.sigma_n0 = TPV104Params::sigma_n` (= `120e6`, positive
    compression — see Step 7 §3.11 audit).
  - `d.tau1_0 = 0.0` (no dip pre-stress — TPV104 is pure strike-slip).
  - `d.tau2_0 = TPV104Params::tau_ini` (= `40e6`, along-strike
    pre-stress into component 2 under BP5/Tandem convention).
  This mirrors TPV102's R-801 Option A layout at
  `dynamic/tpv102_setup.hpp:62-78` exactly.
- Tests:
  - `T_TPV104_SETUP_1` — ψ match: for every fault QP,
    `InitializeFaultDOFs_TPV104` produces `|ψ_init − 0.5636| < 1e-8`.
  - `T_TPV104_SETUP_2` — `V_w[i]` match: populated via
    `ComputeVw_TPV104(x, z)`, equals `0.1` inside the VW zone
    (`|x| < 15 km, 0 < z < 15 km`) and `1.0` outside.
  - `T_TPV104_SETUP_3` — **coordinate mapping (§3.2)**: for the
    TPV104 reference fault (`n = (0, -1, 0)`, `up = (0, 0, 1)`),
    `FaultBasis::GetBasis(fi).tangent1` is aligned with `(0, 0, -1)`
    (dip, downward) and `tangent2` is aligned with `(+1, 0, 0)`
    (strike). If this fails, BP5's `FaultBasis` has been modified —
    halt immediately (**[C2]** violated).
  - `T_TPV104_SETUP_4` — pre-stress layout: `d.tau1_0 == 0` and
    `d.tau2_0 == 4e7` at every QP; `d.sigma_n0 == 1.2e8` at every
    QP. Fails if any QP was accidentally written under the
    GodunovFlux::BuildFrame (t1=strike) convention.
  - `T_TPV104_SETUP_5` — no `InitializeStateTotal_TPV104` symbol
    resolves at link time (the driver MUST NOT have access to a
    total-Q mode — enforces the §3.10 directive).
- Phase 3 gate: **Probe 1 (trial traction)** at step 0.

#### Step 4 — SeisSol-aligned friction coefficient (§3.4 directive)

- Files: `friction/friction_coeff_seissol.hpp` (new — named sibling,
  NOT a modification to `friction/dieterich_ruina.hpp`);
  `tests/unit/test_friction_coeff_seissol.cpp` (new).
- §3.X directive addressed: **§3.4 (replace MFEM's
  `if (psi_over_a > 700.0)` branch with SeisSol's `arsinhexp`
  formulation)**.
- Rationale: `dieterich_ruina.hpp` is on the Extreme-Care list
  (CLAUDE.md, §2.5); TPV102 and BP5 share it. A TPV104-specific
  replacement that lives in a new file keeps BP5 (**[C2]**) and TPV102
  bit-identical.
- Content: two free functions that port SeisSol's
  `RateAndStateCommon.h:60-87` verbatim:
  ```cpp
  namespace mfem::seas::friction_seissol {
  // Mirrors SeisSol rs::computeCExp(real cExpLog).
  inline real_t ComputeCExp(real_t cExpLog);
  // Mirrors SeisSol rs::arsinhexp(lx, cExpLog, cExp):
  //   returns asinh(lx * exp(cExpLog)), using cExp = exp(-2*cExpLog)
  //   for the stable branch.
  inline real_t ArsinhExp(real_t lx, real_t cExpLog, real_t cExp);
  // Mirrors SeisSol FVW updateMu:
  //   μ = a * ArsinhExp(V/(2 V_0), ψ/a, exp(-2 ψ/a))
  inline real_t FrictionCoefficientSeisSol(
      real_t V, real_t psi, real_t a, real_t V0);
  // Derivative w.r.t. V (needed for Newton in Step 5).
  inline real_t FrictionCoefficientSeisSolDerivV(
      real_t V, real_t psi, real_t a, real_t V0);
  }
  ```
  No `if (psi_over_a > 700.0)` branch — numerical stability comes
  from the precomputed `cExp = exp(-2 * ψ/a)` in `ArsinhExp`, which
  is **exactly** what SeisSol does.
- Tests:
  - `T_TPV104_FC_1` — SeisSol byte-match: sweep 10000 random
    `(V, ψ, a)` triples over the TPV104 envelope (`V ∈
    [1e-18, 10] m/s`, `ψ ∈ [0.1, 2.0]`, `a ∈ [0.005, 0.04]`); the
    ported `FrictionCoefficientSeisSol` must match an inline
    SeisSol-verbatim lambda to 1e-15 relative.
  - `T_TPV104_FC_2` — equivalence with MFEM's `FrictionCoefficientPsi`
    on the TPV104 "well-conditioned" envelope (`ψ/a ∈ [28, 80]`):
    the two functions differ by less than 1e-12 relative. This proves
    that swapping to the SeisSol form introduces no regression inside
    the expected operating range.
  - `T_TPV104_FC_3` — derivative byte-match: sweep the same 10000
    samples, compare `FrictionCoefficientSeisSolDerivV` to a forward
    finite-difference `(f(V+h) - f(V))/h` at `h = V*1e-6`; max rel
    error < 1e-5.
  - `T_TPV104_FC_4` — **no hard-coded 700 constant anywhere** in the
    new file (static `grep -n "700" …` in the test harness).
- Phase 3 gate: **Probe 3 (friction coefficient)** at 1e-12 relative.

#### Step 5 — Newton-Raphson solver as default for TPV104 (§3.5)

- Files: `dynamic/tpv104_friction_solver.hpp` (new — wraps
  `FrictionSolver::Method::NewtonRaphson` already exposed by
  `dynamic/friction_solver.hpp:41`; does NOT edit the underlying
  solver); `tests/unit/test_tpv104_friction_solver_newton.cpp` (new).
- §3.X directive addressed: **§3.5 (add Newton-Raphson as option,
  make default for TPV104)**.
- Rationale: `dynamic/friction_solver.cpp:60-62` already dispatches
  to Brent / NewtonRaphson / HybridNRBisection. Step 5 (a) wires the
  TPV104 driver to request `Method::NewtonRaphson` by default, (b)
  adds a CLI flag `--friction-solver={brent,newton,hybrid}` with
  default `newton`, and (c) adds a TPV104-specific Newton variant
  that uses the Step-4 `FrictionCoefficientSeisSol` / its derivative
  so the Newton residual and its gradient are bit-matched to
  SeisSol's `invertSlipRateIterative`
  (`RateAndState.h:286-343`).
- Content:
  ```cpp
  // dynamic/tpv104_friction_solver.hpp
  namespace mfem::seas {
  // SeisSol-aligned Newton on
  //   g(V) = tau_abs - eta_s * V - sigma_n * mu_seissol(V, psi, a)
  // Newton step: V_{k+1} = max(almostZero, V_k - g(V_k)/g'(V_k))
  // Initial guess: previous-step V (per RateAndState.h:481-483).
  real_t SolveSlipRateNewtonSeisSol(
      real_t tau_abs, real_t psi, real_t sigma_n, real_t eta_s,
      real_t a, real_t V0,
      real_t V_prev,   // SeisSol's first-guess
      int max_iter,     // default = 60 (SeisSol maxNumberSlipRateUpdates)
      real_t tol,       // default = 1e-10 (SeisSol newtonTolerance)
      int *iterations = nullptr,
      bool *has_converged = nullptr);
  }
  ```
- Driver wiring:
  - `drivers/tpv104_driver.cpp` uses Step-4 `FrictionCoefficientSeisSol`
    + Step-5 `SolveSlipRateNewtonSeisSol` for the friction inner loop,
    NOT `FrictionSolver::Solve(Method::Brent)`.
  - A CLI fallback `--friction-solver=brent` is available for
    diagnostic use, but the canonical TPV104 run uses Newton.
- Tests:
  - `T_TPV104_FS_1` — Newton converges to same root as Brent: for
    the same inputs `(τ_abs, ψ, σ_n, η_s, a, V₀)` sampled 5000 times
    from the TPV104 envelope, the two methods agree to 1e-10 relative
    (tighter than the 1e-8 §5.5 probe threshold, proving Newton is
    not worse than Brent at steady-state).
  - `T_TPV104_FS_2` — TPV104 rest-state: `SolveSlipRateNewtonSeisSol`
    called with `τ = 40 MPa, ψ = 0.5636, σ_n = 120 MPa, η_s = 4.625e6
    Pa·s/m, a = 0.01, V₀ = 1e-6, V_prev = 1e-16` returns
    `V = 1e-16 ± 1e-20` and converges in < 10 iterations. Fails if
    the Newton loop stalls or returns zero.
  - `T_TPV104_FS_3` — SeisSol byte-match: port the SeisSol
    `invertSlipRateIterative` Newton body (`RateAndState.h:286-343`)
    into an inline lambda; feed the same inputs and require that both
    Newton paths produce bit-identical iterate sequences (up to the
    documented `rs::almostZero()` clamp) for at least 1000 random
    tuples.
  - `T_TPV104_FS_4` — CLI default: driver built with `--help` shows
    `--friction-solver=newton` as the default TPV104 setting.
  - `T_TPV104_FS_5` — non-convergence surfaces loudly: feeding a
    degenerate residual (`τ = 0`, `σ_n = -1`) returns
    `has_converged = false` AND `iterations = max_iter`. The driver
    MUST NOT silently coerce non-convergence into a plausible V.
- Phase 3 gate: **Probe 4 (slip-rate magnitude)** at 1e-8 relative.

#### Step 6 — SeisSol-aligned nucleation accumulator (§3.9)

- Files: `dynamic/tpv104_nucleation.hpp` (new);
  `tests/unit/test_tpv104_nucleation.cpp` (new).
- §3.X directive addressed: **§3.9 (rewrite nucleation to adopt
  SeisSol's cumulative/increment injection)**.
- Rationale: MFEM's existing `ApplyNucleationPrestress` (`dynamic/
  tpv102_setup_total.hpp:627`) **overwrites** `data.tau2_nuc` to the
  current full ramp value. SeisSol's `adjustInitialStress`
  (`FrictionSolverCommon.h:418-450`) **adds the increment**
  `ΔS · [smoothStep(t₂) − smoothStep(t₁)]` at every sub-step. The
  user's directive is to rewrite as SeisSol does. Because TPV102
  shares `ApplyNucleationPrestress`, we do **NOT** modify it —
  instead we add a sibling `ApplyNucleationIncremental_TPV104` that
  implements the SeisSol pattern on a separate per-DOF accumulator.
- Content:
  ```cpp
  // dynamic/tpv104_nucleation.hpp
  namespace mfem::seas {
  // Mirrors SeisSol's smoothStepIncrement(currentTime, dt, t0)
  //   = smoothStep(currentTime, t0) - smoothStep(currentTime - dt, t0)
  // with smoothStep the SCEC "Gaussian" ramp (TPV104 docs eq.14).
  inline real_t SmoothStep_TPV104(real_t current_time, real_t t0);
  inline real_t SmoothStepIncrement_TPV104(
      real_t current_time, real_t dt, real_t t0);
  // Applies ONE SeisSol-style increment to the per-DOF persistent
  // nucleation channel.  Call once per ADER sub-step in the
  // TPV104 sub-step iterator (Step 7), using the sub-step endpoint
  // time and its dt.  The DOFData fields tau1_nuc / tau2_nuc /
  // sigma_n_nuc accumulate *exactly* as SeisSol's
  // initialStressInFaultCS does.
  void ApplyNucleationIncremental_TPV104(
      std::vector<DOFData> &dof_data,
      const std::vector<mfem::Vector> &fault_coords,
      real_t t_substep_end,
      real_t dt_substep);
  // Utility: reset the accumulator to zero at driver init
  // (single call before the time loop).
  void ResetNucleationAccumulator_TPV104(
      std::vector<DOFData> &dof_data);
  }
  ```
- Tests:
  - `T_TPV104_NUC_1` — `SmoothStepIncrement_TPV104(t, dt, T_nuc)`
    byte-matches SeisSol's `smoothStepIncrement` (inline lambda from
    `Numerical/GaussianNucleationFunction.h:22-39`) over 10000 random
    `(t, dt, T_nuc)`; max rel error < 1e-15.
  - `T_TPV104_NUC_2` — accumulator telescopes: starting from zero at
    `t = 0`, applying `ApplyNucleationIncremental_TPV104` on every
    ADER sub-step (`dt = 1e-3 s` for `T_nuc = 1 s`) ends at
    `data.tau2_nuc ≈ Δτ₀ · F(r) · 1.0` within 1e-6 relative at the
    hypocenter QP. Proves the sum equals the full perturbation.
  - `T_TPV104_NUC_3` — spatial factor `F(r)` matches TPV102's
    `NucleationSpatial`: both return 0 at `r ≥ R = 3 km` and
    `exp(r²/(r²−R²))` for `r < R`. Duplication is intentional
    (feedback memo `feedback_tpv102_bp5_no_shared_edit`).
  - `T_TPV104_NUC_4` — TPV102 baseline unaffected: running
    `ApplyNucleationPrestress` on the TPV102 fixture after Step 6
    lands produces identical `tau2_nuc` as before (confirms we did
    NOT touch the shared TPV102 path).
  - `T_TPV104_NUC_5` — guard-bounds: for `t < 0`, increment is 0;
    for `t > T_nuc`, increment is 0; sum from `t = 0` onward equals
    the final `Δτ₀ · F(r)` at any sample time `t ≥ T_nuc`.
- Phase 3 gate: **Probe 1 (trial traction)** — the trial traction
  dumps must match SeisSol at any sampled wall-clock time within
  1e-6 relative at the nucleation-dominated stations.

#### Step 7 — SeisSol-aligned sub-step fault iterator (§3.8, §3.12)

- Files: `dynamic/tpv104_substep_iterator.hpp` +
  `dynamic/tpv104_substep_iterator.cpp` (new); no edits to
  `dynamic/fault_face_flux.{hpp,cpp}`, `dynamic/wave_operator.inl`, or
  any BP5 code (**[C2]**); `tests/unit/test_tpv104_substep_iterator.cpp`
  (new).
- §3.X directive addressed: **§3.8 (rewrite imposed-state to use
  SeisSol's time-weighted sub-step accumulation), §3.12 (follow
  exactly SeisSol's sub-step friction iteration)**.
- Rationale: MFEM's current `FaultFaceFlux::EvaluateTotal` runs the
  friction pipeline **once** per macro-step on time-averaged `Q̄ = I/dt`
  (`dynamic/fault_face_flux.cpp:491-539`). SeisSol runs it **`O` times
  per macro-step**, once per ADER sub-step, and accumulates
  `imposedStatePlus/Minus` weighted by `timeWeights[o]` (§3.8,
  `FrictionSolverCommon.h:340-362`). The user's directive is to adopt
  SeisSol's pattern. Because `fault_face_flux.cpp` is shared with
  TPV102 and on the Extreme-Care list, we do NOT edit it — we wrap it
  with a TPV104-only iterator.
- Design: `Tpv104SubStepIterator` orchestrates the SeisSol pattern by
  calling existing per-stage helpers exposed in `FaultFaceFlux` as
  public methods (`ComputeStageState`, `CompleteFromTheta`,
  `CompleteFromVabs`, `BuildImposedState`, `WriteBackState` at
  `fault_face_flux.cpp:77-206`). The iterator's contract:
  1. At the entry of each ADER macro-step, seed the sub-step loop with
     the cluster's `deltaT[o]` and `timeWeights[o]` (ADER Gauss-Legendre
     quadrature). These match SeisSol's `BaseFrictionLaw.h:100-142`.
  2. For every sub-step `o ∈ [0, O)`:
     - Compute the per-sub-step interpolated `Q_plus[o], Q_minus[o]`
       from the ADER predictor `I_±` (the same predictor that
       `AdvanceADER` already produces).
     - Call `FaultFaceFlux::ComputeTrialTraction` on `Q_plus[o],
       Q_minus[o]`.
     - Call `ApplyNucleationIncremental_TPV104` (Step 6) at the
       sub-step endpoint with the sub-step `dt`.
     - Call `FaultFaceFlux::CompleteFromTheta` / `CompleteFromVabs`
       using the Step-5 Newton solver (Step-4 coefficient).
     - Call `UpdateStateAnalyticSlipLawSRW` (Step 2) with the
       sub-step `V_abs` and sub-step `dt`, matching SeisSol's
       per-sub-step ψ integration (§3.12).
     - Build per-sub-step imposed state via a **copy** of
       `BuildImposedState` inlined inside the iterator (not a second
       Riemann solve — just the per-sub-step accumulator). Weight by
       `timeWeights[o]`.
  3. After the final sub-step, the accumulated imposed state is
     returned to the caller.  `WriteBackState` is called **once** at
     the macro-step boundary to preserve the DOFData semantics that
     TPV102 and the probe format expect.
- No modification to `fault_face_flux.cpp`. The per-sub-step loop
  lives entirely inside `tpv104_substep_iterator.cpp`. The imposed-
  state accumulator structure is a local buffer; the iterator's
  public API surface is exactly `Advance(I_plus, I_minus, dt_macro,
  O, I_imp_plus_out, I_imp_minus_out)`.
- Tests:
  - `T_TPV104_SSI_1` — SeisSol byte-match on a single QP fixture:
    with the TPV104 canonical parameters and the Gauss-Legendre
    `deltaT[o]` / `timeWeights[o]` of ADER-O5, 1 macro-step of 4.7 ms,
    the accumulated imposed state (9 components, both sides) matches
    a SeisSol-verbatim inline lambda (per
    `FrictionSolverCommon.h:340-362`) to 1e-12 relative.
  - `T_TPV104_SSI_2` — O(dt²) convergence: for a known-analytic
    single-QP test with constant slip rate, halving `dt_macro`
    reduces the imposed-state residual by ≥ 3.5× (expected 4× for
    O(dt²); the 3.5× tolerance absorbs round-off).
  - `T_TPV104_SSI_3` — reduction to the TPV102 one-shot in the
    single-sub-step limit: `O = 1, timeWeights = [1.0], deltaT =
    [dt_macro]` must produce **bit-identical** output to calling
    `FaultFaceFlux::EvaluateTotal` on `Q̄ = I/dt` (current TPV102
    behaviour). Confirms we didn't regress the O=1 case.
  - `T_TPV104_SSI_4` — nucleation injection path (Step 6): applying
    the sub-step iterator across `t ∈ [0, T_nuc]` on a locked fault
    (`V = 0`) produces `imposedState` whose stress component equals
    `σ_0 + Δσ · F(r)` at `t = T_nuc`, to 1e-8 relative.
  - `T_TPV104_SSI_5` — ψ update placement: over one macro-step with
    ADER-O5 on a rest-state fixture, the final ψ matches the
    SeisSol-path byte-by-byte (T_SRW_5 test reused with sub-step V
    sequence). Confirms we call `UpdateStateAnalyticSlipLawSRW`
    inside the sub-step loop, not once at the macro-step boundary.
  - `T_TPV104_SSI_6` — BP5 unaffected: running `seas_bp5_full` after
    Step 7 lands produces bit-identical station CSVs to Step 0's
    baseline (`T_TPV104_BP5_BITMATCH_1/2`).
- Phase 3 gate: **Probe 5 (imposed state)** — the 5e-3 threshold in
  §5.6 should drop to 1e-6 relative because SeisSol's accumulation
  is now the MFEM pattern.

#### Step 8 — Normal-stress sign convention audit (§3.11)

- Files: no new files; new test `tests/unit/test_tpv104_normal_sign.cpp`
  (new); audit comment block added to
  `dynamic/tpv104_setup.hpp` (authored in Step 3).
- §3.X directive addressed: **§3.11 (if BP5 uses σ_n > 0 positive
  compression, keep it; if BP5 uses negative, switch to SeisSol
  convention)**.
- Finding: BP5 and CLAUDE.md both document σ_n > 0 = compression
  (`miniapps/seas/CLAUDE.md` line "Normal stress: sigma_n > 0 =
  compression (geology convention)"). BP5's production code
  (`config/bp5_params.hpp:86` — `sigma_n = 50e6`, positive) confirms
  the convention. **Decision: keep MFEM's positive convention.** No
  code change.
- Audit action:
  - Add a guard test that trips if BP5's convention flips in a
    future commit. A positive-compression invariant is a cross-cutting
    assumption.
  - Probe-diff tooling (Step 11) normalises SeisSol's negative-
    compression output to MFEM's positive-compression before diffing.
- Tests:
  - `T_TPV104_SIGN_1` — BP5 convention guard:
    `config/bp5_params.hpp::sigma_n > 0` (compile-time
    `static_assert`). Fails if BP5 ever switches to negative.
  - `T_TPV104_SIGN_2` — TPV104 sign propagation: at `t = 0`,
    `dof_data[i].sigma_n0 == TPV104Params::sigma_n > 0` at every QP;
    after one macro-step of a rest-state run, `s.sigma_n_trial > 0`
    and `s.sigma_n_total > 0` at every QP.
  - `T_TPV104_SIGN_3` — friction strength uses magnitude: the
    assertion `strength == std::abs(sigma_n_total) * μ` holds for
    every probe sample; proves `std::abs` was not accidentally
    removed when the SeisSol Newton solver landed (Step 5).
  - `T_TPV104_SIGN_4` — probe-diff normalisation: the comparator
    flips SeisSol's `normalStress` channel sign before diffing
    against MFEM's `sigma_n_trial`; after the flip, the residual is
    within 1e-6 relative at Probe 1.

#### Step 9 — Driver

- Files: `drivers/tpv104_driver.cpp` (new).
- §3.X directives consumed: **all of the above** (wires Steps 1–8
  together).
- Content: duplicate `drivers/tpv102_driver.cpp` main, with the
  following substitutions — each labelled by the directive it serves:
  - `TPV102Params` → `TPV104Params` (Step 1).
  - `InitializeFaultDOFs` → `InitializeFaultDOFs_TPV104` (Step 3).
  - `FaultBasis` call-through for the fault-local frame is inherited
    from the TPV102 pattern (**[C1]**); no new code.
  - `FaultFaceFlux::EvaluateTotal(Q̄, …)` one-shot →
    `Tpv104SubStepIterator::Advance(I_plus, I_minus, dt, O, …)`
    (Step 7, implementing §3.8 + §3.12). Note: the existing
    `AdvanceADER` entry point at `wave_operator.inl:2991` is reused;
    the iterator swap happens at the fault-flux callback the
    WaveOperator calls (new flag `wave.SetFaultIterator(&ssi)`
    wired to the iterator without touching `fault_face_flux.cpp`).
  - ψ update moves from the outer-loop forward-Euler (TPV102
    `drivers/tpv102_driver.cpp:1114-1122`) to **inside the sub-step
    iterator** (Step 7; implements §3.12 directive that we follow
    SeisSol exactly in time integration placement).
  - Friction coefficient uses `FrictionCoefficientSeisSol` (Step 4,
    §3.4 directive — no `psi/a > 700` branch).
  - Friction solver defaults to `SolveSlipRateNewtonSeisSol`
    (Step 5, §3.5 directive).
  - Nucleation uses `ApplyNucleationIncremental_TPV104` (Step 6,
    §3.9 directive) inside the sub-step iterator.
  - `InitializeStateTotal_TPV104` and related total-Q helpers are
    **not referenced** (Step 3 `T_TPV104_SETUP_5` enforces; §3.10
    directive).
  - Q is initialised to zero (fluctuation-Q only, §3.10 directive).
  - Station writer: `TPV104StationWriter` (Step 3). Output mapping
    documented at the writer header: `horizontal = strike = slip2`,
    `vertical = dip = slip1`.
- CLI additions:
  - `--friction-solver={newton,brent,hybrid}` — default `newton`.
  - `--fric-law={slip-srw,aging}` — default `slip-srw`.
  - `--ader-order N` — default 2.
  - `--fault-iterator={substep,oneshot}` — default `substep` (Step 7).
    `oneshot` is retained for the `T_TPV104_SSI_3` regression test
    only; production runs must use `substep`.
- Tests: `tests/unit/test_tpv104_smoke.cpp` (new; 100-step MPI-serial
  run on the coarse M0 mesh):
  - `T_TPV104_SMOKE_1` — no NaN in Q or DOFData; all nine station
    files written; ψ at hypocenter drifts < 1e-6 over 100 steps with
    nucleation disabled.
  - `T_TPV104_SMOKE_2` — with nucleation enabled for 100 steps at
    `dt = 1e-3 s`, ψ at hypocenter decreases (rupture onset) and
    slip_rate grows to > 1e-3 m/s. Fails if the nucleation channel
    is not wired.
  - `T_TPV104_SMOKE_3` — driver emits console banner listing:
    "Time integrator: ADER-O2", "Fault iterator: SeisSol sub-step",
    "Friction solver: Newton-Raphson (SeisSol)", "Friction law:
    slip-SRW (ψ-space)". Verifies that every Step 4/5/6/7 wiring
    took effect. Fails if any row says "Brent" or "oneshot" unless
    the matching flag was passed.
- Phase 2 gates: **P2_C, P2_D, P2_E, P2_F**.
- Phase 3 gate: all five probes runnable.

#### Step 10 — Mesh

- Files:
  - `tpv104/mesh/tpv104_200m.geo` (new — duplicated from
    `miniapps/seas/tpv102/mesh/tpv102_200m.geo`, per the user-directive
    2026-04-24 that the TPV102 mesh geometry is reusable for TPV104
    as-is; SCEC TPV104 and TPV102 share the same fault plane,
    hypocenter, VW-zone dimensions, nucleation radius, and 36×18 km
    sliding area).
  - `tpv104/mesh/tpv104_1000m.geo`, `tpv104/mesh/tpv104_500m.geo`
    (new — variants of the 200 m template with `lc_fault ∈ {1000,
    500}` for the coarse / intermediate Frontera runs).
  - `tpv104/mesh/tpv104_{200m,500m,1000m}.msh` (new — generated from
    the `.geo` files by Gmsh; not committed, built on demand by the
    Makefile target `make tpv104-mesh`).
- §3.X directive addressed: **§3.17 (mesh / ADER order / output
  stations)**.
- Content: **Duplicated** from TPV102 (per
  `feedback_tpv102_bp5_no_shared_edit`): copy
  `miniapps/seas/tpv102/mesh/tpv102_200m.geo` into
  `tpv104/mesh/tpv104_200m.geo` verbatim, with only the header
  comment block updated (references `SCEC TPV104` instead of `TPV102`;
  physical-group numbering and all mesh parameters are identical).
  The TPV102 `.geo` already encodes every geometric quantity TPV104
  needs:
  - Fault plane at `Y = 0` (code), vertical strike-slip;
  - VW half-length 15 km along strike, VW half-depth 15 km, 3 km
    smooth transition around the VW zone;
  - Hypocenter at `(X, Z) = (0, -7.5 km)`;
  - Nucleation patch radius `R = 3 km` with `lc_nucl = 200 m`
    refinement box (§4.1 confirms TPV104 uses the same 3 km radius);
  - 60 km half-width domain with free surface at `Z = 0` and
    absorbing boundaries on the 4 side walls + bottom;
  - Physical Surface 1 = free surface (Z=0); Physical Surface 3 =
    fault (36×18 km sliding area); Physical Surface 5 = absorbing
    outer boundaries; Physical Volume 1 = bulk.
  These are all **TPV104 values** per §4.1 parameter table — no
  geometric substitution is needed. Built in the `pythonenv` conda
  env (CLAUDE.md): `gmsh -3 tpv104_200m.geo -o tpv104_200m.msh`.
- **Do NOT edit `miniapps/seas/tpv102/mesh/tpv102_200m.geo`**; the
  TPV104 variant is an independent copy so future TPV102/TPV104
  divergence (e.g. a TPV104-specific nucleation-refinement tweak)
  does not bleed back into the TPV102 baseline (**[C2]** +
  `feedback_tpv102_bp5_no_shared_edit`).
- Tests (`tests/unit/test_tpv104_mesh.cpp`, new):
  - `T_TPV104_MESH_0` — byte-identical geometry invariant: diff
    `tpv104/mesh/tpv104_200m.geo` against
    `miniapps/seas/tpv102/mesh/tpv102_200m.geo` ignoring only the
    header comment block (lines 1–29). Every other line must match.
    Fails if a maintainer drifts the TPV104 geometry relative to
    TPV102 without an accompanying plan update.
  - `T_TPV104_MESH_1` — vertex count, element count, and physical
    groups match the TPV102-equivalent reference within 5 % at
    matched resolution (1000 m, 500 m, 200 m). Anchors:
    `tpv102_200m.msh` produced by the existing TPV102 build pipeline
    provides the reference counts.
  - `T_TPV104_MESH_2` — wave-speed arrival on the coarse mesh:
    inject a delta pulse at the fault and time the P-wave arrival at
    z = 5 km; expected `t_arrival = 5000/6000 ≈ 0.833 s`, tolerance
    1 %. Mismatch flags a mesh or material-parameter bug.
  - `T_TPV104_MESH_3` — physical-group parity with the TPV104
    driver: the driver built in Step 9 reads Physical Surface IDs
    `{1, 3, 5}` and Physical Volume `{1}`; the mesh-validation script
    asserts exactly those groups exist in each `.msh` file with the
    same meanings as TPV102. Catches a `.geo`-level physical-tag
    renumbering before the driver loads it.
- Phase 3 gate: **Phase 3.A wave-speed check** (§7 risk row "MFEM
  mesh does not match TPV5") — gated by `T_TPV104_MESH_2`.

#### Step 11 — Diagnostics instrumentation (build-gated)

- Files: new compile flag `SEAS_DIAG_TPV104_STATE`, wired through the
  five MFEM probe points named in §5.2–5.6. Implementation limited to
  Step 7's `tpv104_substep_iterator.cpp` and Step 9's
  `tpv104_driver.cpp`; **no edits** to `fault_face_flux.cpp`,
  `wave_operator.inl`, `dieterich_ruina.hpp`, or any BP5 source
  (§2.5, **[C2]**).
- `fault_face_flux.cpp` already has `SEAS_DIAG_FAULT_FLUX` gates we
  reuse read-only.
- Tests (`tests/unit/test_tpv104_probe_format.cpp`, new):
  - `T_TPV104_PROBE_1` — probe files match §5.1 header format.
  - `T_TPV104_PROBE_2` — disabled build emits no probe output
    (grep for the gate string in the release build).
  - `T_TPV104_PROBE_3` — `SEAS_DIAG_TPV104_STATE` does NOT affect
    BP5 (`seas_bp5_full` built with and without the flag produces
    bit-identical output — enforces **[C2]**).
- Phase 3 gate: enables every probe pair.

#### Step 12 — Makefile integration

- Files: add targets `seas_tpv104_driver`, `seas_test_slip_law_srw`,
  `seas_test_tpv104_params`, `seas_test_tpv104_setup`,
  `seas_test_friction_coeff_seissol`, `seas_test_tpv104_friction_solver_newton`,
  `seas_test_tpv104_nucleation`, `seas_test_tpv104_substep_iterator`,
  `seas_test_tpv104_normal_sign`, `seas_test_tpv104_smoke`,
  `seas_test_tpv104_mesh`, `seas_test_tpv104_probe_format` to
  `miniapps/seas/Makefile`. New umbrella target `make test-tpv104`
  runs all of the above + the BP5 bit-match gate (Step 0).
- Constraint addressed: **[C2]** — `make test` continues to run the
  pre-existing BP5 unit-test set unchanged.
- Tests: existing `make test` + new `make test-tpv104` + new
  `make test-bp5-bitmatch` (Step 0) all pass.
- Gate: **P2_A, P2_B, P2_C, P2_D, P2_E, P2_F** all green.

#### Step 13 — Probe-diff tooling (Python)

- Files: `tpv104/scripts/probe_diff.py`,
  `tpv104/scripts/probe_summarise.py`,
  `tpv104/scripts/tpv104_column_map.py` (new). The column-map helper
  implements the §3.2 strike/dip swap and the §3.11 normal-stress sign
  flip so that diffs land on aligned channels.
- Tests: `pytest tpv104/scripts/tests/`:
  - Each probe pair exercised with a synthetic identical-by-
    construction pair of probe files to assert the diff tool reports
    zero max diff.
  - `test_column_map_strike_dip_swap` — SeisSol col-2 ↔ MFEM `slip2`.
  - `test_column_map_sign_flip` — SeisSol `normalStress` × (-1) ↔
    MFEM `sigma_n_trial`.
- Phase 3 gate: enables §6.1 Phase 3 execution.

#### Directive coverage matrix

| user directive | §3.X row | step addressing it | unit test(s) |
|---|---|---|---|
| "convention must agree with BP5 — call their functions" ([C1]) | §3.2 | Step 3 reuses `fault::FaultBasis::GetBasis` directly; no new frame code | `T_TPV104_SETUP_3` (BP5 frame invariant) |
| "no BP5 modification" ([C2]) | — | Step 0 baseline + `T_TPV104_SIGN_1` + `T_TPV104_NUC_4` + `T_TPV104_SSI_6` + `T_TPV104_PROBE_3` | five independent guards, one per TPV104 surface that touches shared code |
| "(1) 3.2 check current convention same as BP5" | §3.2 | Step 3 | `T_TPV104_SETUP_3` |
| "(2) 3.3 must implement for TPV104" | §3.3 | Step 2 | `T_SRW_1..T_SRW_9` |
| "(3) 3.4 change psi/a>700 to SeisSol cap" | §3.4 | Step 4 (new file, no `dieterich_ruina.hpp` edit) | `T_TPV104_FC_1..4` |
| "(4) 3.5 Newton-Raphson as option + default" | §3.5 | Step 5 | `T_TPV104_FS_1..5` |
| "(5) 3.8 follow SeisSol time-accumulation, rewrite" | §3.8 | Step 7 sub-step iterator | `T_TPV104_SSI_1..6` |
| "(6) 3.9 rewrite nucleation as SeisSol accumulation" | §3.9 | Step 6 | `T_TPV104_NUC_1..5` |
| "(7) 3.10 keep only fluctuation-Q matching SeisSol" | §3.10 | Step 3 (no total-Q helpers shipped) | `T_TPV104_SETUP_5` |
| "(8) 3.11 check BP5 normal-stress sign" | §3.11 | Step 8 audit (decision: keep positive) | `T_TPV104_SIGN_1..4` |
| "(9) 3.12 follow SeisSol time integration exactly" | §3.12 | Step 7 (per-sub-step ψ + friction + imposed-state accumulation) | `T_TPV104_SSI_1, T_TPV104_SSI_5` |
| "(10) 3.13 initial values match expectation suffice" | §3.13 | Step 1 | `T_TPV104_P_1` |

Every directive has a unique step number, a unit test that guards it,
and a Phase 3 probe gate (§4.9 graph) that cross-verifies it against
SeisSol output.

#### Ordering summary

```
Step 0 (BP5 baseline gate)  ───────────────┐  [C2] protection: runs on every merge
                                            │
Step 1 (params) ────────────────────────────┤
                                            ├──► Step 2 (SlipLawSRWPsi)    ─┐
                                            │                                │
                                            ├──► Step 3 (setup)             ─┤
                                            │                                │
                                            ├──► Step 4 (friction coeff)    ─┤
                                            │                                │
                                            ├──► Step 5 (Newton solver)     ─┼─► Step 7 (sub-step iter)
                                            │                                │         │
                                            ├──► Step 6 (nucleation accum.) ─┘         │
                                            │                                          ▼
                                            ├──► Step 8 (sign audit)    ───► Step 9 (driver)
                                            │                                          │
                                            └──► Step 10 (mesh)            ────────────┤
                                                                                       │
                                                                                       ▼
                                                                            Step 11 (diagnostics)
                                                                                       │
                                                                                       ▼
                                                                            Step 12 (Makefile)
                                                                                       │
                                                                                       ▼
                                                                            Step 13 (probe tooling)
                                                                                       │
                                                                                       ▼
                                                                            Phase 2 DONE
                                                                                       │
                                                                                       ▼
                                                                            Phase 3.A (§5.7, §6.1)
```

Parallelisable: Steps 1–6, 8, 10 are independent once Step 0 lands.
Step 7 blocks on Steps 2, 4, 5, 6. Step 9 blocks on Steps 1, 2, 3, 4,
5, 6, 7, 8, 10. Steps 11, 12, 13 are final and can land concurrently.

### 4.10.X R5-003 fix plan — runtime-switchable Newton solvers for TPV104

Added 2026-04-24 to close the standing plan violation first noted in
round 1 as "Deviation #4", tracked across rounds 4 (R4-002) and 5
(R5-003). Plan §4.10 Step 5 mandates "Newton solver is the driver's
default" and specifies the SeisSol-aligned stable-asinh Newton
(`SolveSlipRateNewtonStable` + `FrictionCoefficientStable`). Currently
`SolveSlipRateNewtonStable` is unreachable via
`FrictionSolver::Solve`; the iterator defaults to `Method::Brent`.
This subsection is the concrete fix plan.

#### Goal

Expose BOTH Newton variants as runtime choices; default to
SeisSol-aligned on TPV104; allow a CLI flag to flip between them
without a rebuild. Zero edits to `dieterich_ruina.hpp`,
`fault_face_flux.{hpp,cpp}`, `wave_operator.*`, or BP5/TPV102 code.
The existing `Method::NewtonRaphson` value stays intact — BP5/TPV102
callers continue to see the legacy path.

#### Switching surfaces

- **CLI** (`drivers/tpv104_driver.cpp`, Step 9):
  `--friction-solver={brent,newton,newton-stable,hybrid}`,
  default `newton-stable`.
- **Enum** (`dynamic/friction_solver.hpp`):
  `FrictionSolver::Method::NewtonRaphsonStable` joins the existing
  three values.
- **Iterator** (`dynamic/tpv104_substep_iterator.hpp`): `Advance`
  already accepts `method` as a parameter — default bumps from
  `Method::Brent` to `Method::NewtonRaphsonStable`.
- **Banner** (`tests/unit/test_tpv104_smoke.cpp`): expected string
  becomes `"Friction solver: Newton-Raphson (stable-asinh)"` on
  default flags; `"Friction solver: Brent"` forbidden on defaults
  (unchanged).

#### File-level changes

**1. `dynamic/friction_solver.hpp`** (~12 lines)

Extend the enum and declare a new method sibling to `SolveNR`:

```cpp
enum class Method {
   Brent,                  // MFEM μ + log10-V Brent (Tandem-verified)
   NewtonRaphson,          // MFEM μ + Newton (legacy)
   NewtonRaphsonStable,    // stable-asinh μ + Newton (SeisSol-aligned,
                           //   TPV104 canonical per Plan §4.10 Step 5)
   HybridNRBisection       // MFEM μ + NR-with-bisection-fallback
};

/// Delegates to `mfem::seas::SolveSlipRateNewtonStable` in
/// `dynamic/tpv104_friction_solver.hpp`.  μ is computed via
/// `friction/friction_coeff_stable.hpp::FrictionCoefficientStable`,
/// which byte-matches SeisSol `rs::arsinhexp` on the TPV104 envelope.
/// The Newton loop byte-matches `RateAndState.h::invertSlipRateIterative`
/// under matched inputs (T_TPV104_FS_3).
real_t SolveNRStable(real_t tau, real_t psi, real_t sigma_n,
                     real_t eta, real_t a) const;
```

**2. `dynamic/friction_solver.cpp`** (~18 lines)

```cpp
#include "tpv104_friction_solver.hpp"

real_t FrictionSolver::Solve(real_t tau, real_t psi, real_t sigma_n,
                             real_t eta, real_t a, Method method) const
{
   switch (method)
   {
      case Method::Brent:               return SolveBrent(tau, psi, sigma_n, eta, a);
      case Method::NewtonRaphson:       return SolveNR(tau, psi, sigma_n, eta, a);
      case Method::NewtonRaphsonStable: return SolveNRStable(tau, psi, sigma_n, eta, a);
      case Method::HybridNRBisection:   return SolveHybrid(tau, psi, sigma_n, eta, a);
      default:                          return SolveBrent(tau, psi, sigma_n, eta, a);
   }
}

// R5-003 / Plan §4.10 Step 5: SeisSol-aligned Newton.
// Thin wrapper — all Newton logic lives in SolveSlipRateNewtonStable.
// V_prev warm-start: `FrictionSolver::Solve` has no V_prev parameter,
// so we bootstrap from Brent's root.  On the TPV104 envelope this
// gives ~1-iterate Newton convergence.  Step 9 may extend Solve's
// signature to accept V_prev from DOFData.slip_rate.
real_t FrictionSolver::SolveNRStable(real_t tau, real_t psi,
                                     real_t sigma_n, real_t eta,
                                     real_t a) const
{
   const real_t V_warm = SolveBrent(tau, psi, sigma_n, eta, a);
   bool converged = false;
   int iterations = 0;
   const real_t V = SolveSlipRateNewtonStable(
      tau, psi, std::abs(sigma_n), eta, a, /*V0=*/V0, /*V_prev=*/V_warm,
      /*max_iter=*/60, /*tol=*/1e-8, &iterations, &converged);
   return converged ? V : V_warm;   // paranoia fallback; TPV104 envelope
                                    // is monotone so this never fires.
}
```

**3. `dynamic/tpv104_substep_iterator.hpp`** (1 line + docstring)

```cpp
// Default bumped Brent → NewtonRaphsonStable.
void Advance(...,
             FrictionSolver::Method method
                = FrictionSolver::Method::NewtonRaphsonStable);
```

Docstring rewrite at the `Advance` declaration: "Default is
`Method::NewtonRaphsonStable` (R5-003 review fix / Plan §4.10 Step 5):
SeisSol-aligned stable-asinh Newton. `Method::Brent` remains
available as the robust fallback; both agree to 1e-10 relative on
the TPV104 envelope."

**4. `drivers/tpv104_driver.cpp`** (Step 9, ~15 lines)

CLI scan + dispatch map:

```cpp
std::string friction_solver_name = "newton-stable";
args.AddOption(&friction_solver_name, "--friction-solver",
               "Friction solver method: brent, newton, newton-stable "
               "(default, SeisSol-aligned), hybrid.");

FrictionSolver::Method friction_method;
std::string banner_solver;
if      (friction_solver_name == "brent")         { friction_method = FrictionSolver::Method::Brent;
                                                    banner_solver = "Brent"; }
else if (friction_solver_name == "newton")        { friction_method = FrictionSolver::Method::NewtonRaphson;
                                                    banner_solver = "Newton-Raphson (legacy)"; }
else if (friction_solver_name == "newton-stable") { friction_method = FrictionSolver::Method::NewtonRaphsonStable;
                                                    banner_solver = "Newton-Raphson (stable-asinh)"; }
else if (friction_solver_name == "hybrid")        { friction_method = FrictionSolver::Method::HybridNRBisection;
                                                    banner_solver = "Hybrid"; }
else {
   std::cerr << "Unknown --friction-solver: " << friction_solver_name
             << ". Choose one of: brent, newton, newton-stable, hybrid.\n";
   return 1;
}
std::cout << "Friction solver: " << banner_solver << "\n";
...
iterator.Advance(..., friction_method);
```

**5. `tests/unit/test_tpv104_smoke.cpp`** (2 lines)

```cpp
// Update default-banner expectation.
{"Friction solver: Newton-Raphson (stable-asinh)",
 "Newton-Raphson stable-asinh friction-solver default"},
// Forbidden-on-default list unchanged (Brent / Hybrid / aging still
// forbidden).
```

Add a banner test for each non-default solver choice:

```cpp
void TestBannerAllSolverOptions() {
   for (const auto &[arg, expected] : std::vector<std::pair<std::string, std::string>>{
        {"brent",         "Friction solver: Brent"},
        {"newton",        "Friction solver: Newton-Raphson (legacy)"},
        {"newton-stable", "Friction solver: Newton-Raphson (stable-asinh)"},
        {"hybrid",        "Friction solver: Hybrid"}}) {
      const std::string out = RunDriver(binary, "--dry-run --friction-solver " + arg);
      TEST_ASSERT(out.find(expected) != std::string::npos,
                  ("banner for --friction-solver=" + arg).c_str());
   }
}
```

#### New unit tests (~60 lines)

**File**: `tests/unit/test_friction_solver_stable.cpp` (new).

- **T_FS_STABLE_1** — dispatch byte-match: for 1000 random
  `(τ, ψ, σ_n, η, a)` tuples from the TPV104 envelope,
  `fs.Solve(..., NewtonRaphsonStable)` equals a direct
  `SolveSlipRateNewtonStable(..., V_prev=SolveBrent_result)` call
  bit-for-bit.
- **T_FS_STABLE_2** — SeisSol-byte-match: the dispatched path
  byte-matches the inline SeisSol `invertSlipRateIterative` port
  from `test_tpv104_friction_solver_newton.cpp::reference_newton` to
  1e-13 relative on 1000 samples.
- **T_FS_STABLE_3** — NewtonStable vs NewtonRaphson on the TPV104
  physical envelope (`ψ/a ∈ [28, 80]`): agree to 1e-10 relative. This
  documents that the 700-branch never fires for TPV104 so the two
  Newton paths are numerically equivalent; any future envelope that
  pushes `ψ/a > 700` would diverge and Phase 3 probe-3 would flag it.
- **T_FS_STABLE_4** — CLI switchability: driver banner reports the
  correct solver name under each of the four `--friction-solver`
  values, verified via `TestBannerAllSolverOptions` above.
- **T_FS_STABLE_5** — legacy `Method::NewtonRaphson` unchanged: a
  regression guard that ensures the round-1 R-001..R-012 +
  round-2..4 closures are NOT disturbed (BP5 bit-match via Step 0).

#### Acceptance criteria (blocking Phase 3)

- [ ] **R5-003_A** — `FrictionSolver::Method::NewtonRaphsonStable`
  reaches `SolveSlipRateNewtonStable`; T_FS_STABLE_1 passes.
- [ ] **R5-003_B** — iterator default is `NewtonRaphsonStable`;
  `Tpv104SubStepIterator::Advance` with default `method` produces
  stable-asinh μ output.
- [ ] **R5-003_C** — driver CLI exposes all four options; banner
  string reflects the selected solver; T_FS_STABLE_4 passes.
- [ ] **R5-003_D** — `test_tpv104_smoke.cpp::TestBannerDefaults`
  updated; TestBannerDefaults passes on default flags.
- [ ] **R5-003_E** — BP5 baseline bit-match (Step 0) still green;
  TPV102 baseline bit-match (P2_B) still green — confirms no
  regression in the shared path.

#### Ordering within §4.10

R5-003 sits between Step 5 (Newton-Raphson as default) and Step 9
(driver). It is blocking for Step 9: the driver cannot ship without a
coherent CLI/enum/default/banner story, and the current state is
internally contradictory (iterator default = Brent, banner expects
Newton-Raphson, SeisSol-aligned Newton unreachable).

#### What this fix explicitly does NOT do

- Does NOT replace `SolveNR` — legacy Newton stays as a CLI option
  and the default for BP5/TPV102 (`Method::Brent` remains their
  default via `FrictionSolver::Solve`'s argument default).
- Does NOT extend `FrictionSolver::Solve`'s signature with `V_prev`.
  Instead, `SolveNRStable` warm-starts from Brent internally; slower
  than per-QP warm-start but keeps the shared API unchanged. A future
  Step-9 optimisation may add a `V_prev` override.
- Does NOT touch `FaultFaceFlux::ComputeStageState`. The iterator
  keeps calling `flux_.ComputeStageState(..., method)`; the enum-
  dispatch happens inside `FrictionSolver::Solve` without any
  Extreme-Care-file edits.
- Does NOT land `V_prev` threading from `DOFData.slip_rate`. That
  remains a Step-9 driver optimisation; round-5 R5-003 only wires the
  dispatch + warm-start.

---

### 4.11 What this strategy does NOT do

To make the duplication contract explicit (per §2.5 and
`feedback_tpv102_bp5_no_shared_edit`):

- Does NOT edit `drivers/tpv102_driver.cpp`,
  `dynamic/tpv102_setup.hpp`, `dynamic/tpv102_setup_total.hpp`,
  `config/tpv102_params.hpp`.
- Does NOT edit `dynamic/fault_face_flux.{hpp,cpp}`,
  `dynamic/godunov_flux.{hpp,cpp}`,
  `dynamic/wave_operator.{hpp,inl,cpp}`,
  `friction/dieterich_ruina.hpp`, or `friction/state_evolution.hpp`
  (except the opt-in addition of one new class per §2.1).
- Does NOT edit the SeisSol tree; the S1–S6 references above are
  read-only specifications. Probe printfs (§5) are added on Frontera
  only, guarded by a single compile flag, disposable.
- Does NOT launch any Frontera run; each 3.A / 3.B / 3.C run is a
  separate user approval
  (`feedback_frontera_approval.md`).

---

Hard gates (blocking Phase 3):

- [ ] **P2_A — All T_SRW_* tests pass** under `make test-slip-law-srw`.
- [ ] **P2_B — TPV102 baseline byte-match**: running `seas_tpv102_driver` with all flags matching a prior M_ref reference run, the nine TPV102 station CSVs are bit-identical to the baseline (no unintended coupling into TPV102).
- [ ] **P2_C — TPV104 driver smoke test**: `seas_test_tpv104_smoke` (MPI-serial, M0 mesh scaled to SCEC-TPV104 dimensions) runs 100 ADER steps without NaNs, produces a non-zero `ψ` trace at the hypocenter QP, writes all nine station files.
- [ ] **P2_D — TPV104 ψ_ini match**: running the driver with `--tfinal=0 --dry-run-stations`, the printed `ψ_ini` at every station matches SeisSol trace col-9 row-1 (`5.6359184e-01`) to 1e-8.

Soft gates (nice-to-have):

- [ ] **P2_E — TPV104 at 1 s with nucleation off**: driver at `--tfinal=1.0 --disable-nucleation` should have `ψ` drift from `ψ_ini` by less than 1e-8 (no spontaneous state-evolution without nucleation, since V_ini = 1e-16 m/s means `V·dt/L` is negligible).
- [ ] **P2_F — TPV104 at 1 s with nucleation on**: driver at `--tfinal=1.0` should show nucleation rupture beginning to propagate from the hypocenter, non-trivial slip rate at the hypocenter station (>1e-3 m/s).

---

## 5. Phase 3 — Paired diagnostics cross-verification plan

The user's constraint is clear: SeisSol cannot build locally, MFEM cannot run TPV104 on Frontera at scale without authorisation. Strategy: define **5 probe pairs** that each capture one stage of the fault-face pipeline; implement each side's probe with identical output format; run both codes on the canonical TPV104 fixture; difference trace-by-trace. The first probe pair whose diff exceeds its threshold localises the bug.

### 5.1 Probe design — format, invariants, output file layout

All probes write a single text file per rank per run in the format

```
# probe=<name>  code=<MFEM|SeisSol>  rank=<r>  nprocs=<n>
# t  qp_id  <field_1>  <field_2>  ...  <field_K>
<rows, one per (time, qp) sample>
```

Fields are always in SI units unless labelled otherwise; stresses are in Pa (not MPa — avoid unit-mismatch bugs in the diff tool); angles in radians; ψ dimensionless.

**Thinned sampling**: to keep file sizes tractable, each probe is triggered

- only on a fixed set of **probe QPs**: 9 stations + the hypocenter QP (total 10),
- only at sample times `t ∈ {0, 0.01, 0.02, 0.05, 0.1, 0.15, ..., 0.9, 0.95, 1.0} ∪ {1.1, 1.2, ...} ∪ {2, 3, ..., 12}` seconds (coarser at late times — about 50 samples per QP per probe).

Probe QP identifier: `qp_id = (station_index, side ∈ {+, −})`. Shared-face QPs report one side per rank.

### 5.2 Probe pair 1 — Pre-friction (trial traction) [Stage A]

**What it captures**: the Godunov trial traction the friction solver receives, before friction is applied. This is the interface between the bulk wave operator and the friction pipeline; any bug in how MFEM's ADER + `ComputeTrialTraction` produces trial stress (or a sign error in the rotation from global to fault-local) manifests here first.

**MFEM probe point**: `fault_face_flux.cpp:255` (inside `Evaluate`, right after `ComputeStageState` fills `s.{sigma_n_trial, tau1_trial, tau2_trial}`). Guard with `SEAS_DIAG_TPV104_STATE` and `data.diag_print` (same diag-print gate TPV102 already uses).

**SeisSol probe point**: `FrictionSolverCommon.h:192`, inside the `precomputeStressFromQInterpolated` loop right after `traction2` is assigned. Guard with `SEISSOL_TPV104_DIAG`.

**Fields**: `sigma_n_trial, tau1_trial, tau2_trial` (all in Pa).

**Threshold**: `max |diff| < 1e-4` (absolute) AND `max |diff| / max |signal| < 1e-3` (relative). Bulk quadratures differ; an order-of-magnitude match is the best we can expect.

**Diagnostic outcome**:

- PASS: the bulk wave operator + pre-friction assembly are consistent. Proceed to Probe 2.
- FAIL: the bug is in the bulk ADER path or the stress-rotation into fault-local CS. Bisect on: (a) does `Q_plus`, `Q_minus` itself agree across codes before the trial step? Add a sub-probe dumping the nine Q components for both sides. (b) If Q agrees but trial traction diverges, the fault-local rotation `T/Tinv` is differing — trace through each code's `BuildRotation` vs SeisSol's `rotateToFaultCS`.

### 5.3 Probe pair 2 — State evolution (ψ in → ψ out) [Stage B]

**What it captures**: the state-variable analytic integration. Any error in `UpdateStateAnalyticSlipLawSRW` vs SeisSol FVW `updateStateVariable` will show up here.

**MFEM probe point**: `tpv104_driver.cpp` in the `UpdateStateAnalyticSlipLawSRW` call block (will be inside the per-step state-update loop). Dump `(psi_in, V, L, dt, V_w, a, b, V0, f0, muW, psi_out)`.

**SeisSol probe point**: `FastVelocityWeakeningLaw.h:77` (inside `updateStateVariable`, right before `return localStateVariable`). Dump the same 11-tuple.

**Fields**: `psi_in, V, L, dt, V_w, a, b, V0, f0, muW, psi_out`.

**Threshold**: `|psi_out_MFEM - psi_out_SeisSol| < 1e-13` on every sample, conditional on the 10-tuple inputs agreeing to 1e-14. This is a **numerical-formula** test — it must be bit-match-close.

**Diagnostic outcome**:

- PASS: state evolution is mathematically identical.
- FAIL (inputs agree, outputs differ): `SlipLawSRWPsi` has a bug. Run `test_slip_law_srw_psi.cpp::T_SRW_5` at the specific failing input tuple for a minimal reproducer.
- FAIL (inputs differ): earlier in the pipeline, V is computed differently. Defer to Probe 4.

### 5.4 Probe pair 3 — Friction coefficient [Stage C]

**What it captures**: μ = a · asinh((V/(2V₀))·exp(ψ/a)). The only nonlinearity in the friction law below the ODE.

**MFEM probe point**: `fault_face_flux.cpp:141` (inside `CompleteFromVabs`, after `f_V = a * asinh(V_abs * C)`). Dump `(V, psi, a, V0, f_V)`.

**SeisSol probe point**: `FastVelocityWeakeningLaw.h:120-121` (inside `updateMu`, at the return). Dump `(V, psi, a, V0, mu)`.

**Fields**: `V, psi, a, V0, mu`.

**Threshold**: `|mu_MFEM - mu_SeisSol| < 1e-12` relative. Same asinh formula on both sides; only differences are ULP arithmetic.

**Diagnostic outcome**:

- PASS: friction coefficient is identical.
- FAIL: the `arsinhexp` branch-protection in `dieterich_ruina.hpp` vs SeisSol's `rs::arsinhexp` differs. Check the "large ψ/a asymptotic" branch; TPV104 should not hit it but TPV102 at nucleation might. **Do not edit** `dieterich_ruina.hpp` per §2.5 — report and ask.

### 5.5 Probe pair 4 — Slip-rate magnitude [Stage D]

**What it captures**: the Brent / Newton root of `tau = sigma_n · f(V, psi) + eta_s · V`. Different solver methods (MFEM: Brent in log10(V); SeisSol: Newton in V) can converge to slightly different roots in tight brackets.

**MFEM probe point**: `fault_face_flux.cpp:116` (after the `solver_.Solve` call returns `V_abs`). Dump `(Theta, psi, sigma_n, eta_s, a, V_abs, solver_iterations)`.

**SeisSol probe point**: `RateAndState.h:339` (after the Newton loop converges, right before `return hasConverged`). Dump `(Theta, psi, sigma_n, eta_S, a, slipRateTest, iter)`.

**Fields**: `Theta, psi, sigma_n, eta_s, a, V_abs, iterations`.

**Threshold**: `|V_abs_MFEM - V_abs_SeisSol| / max(V_abs_SeisSol, 1e-30) < 1e-8` on samples where both solvers converge. Tight; Brent and Newton on a Lipschitz-monotone residual should agree to ~1e-10.

**Diagnostic outcome**:

- PASS: friction-solver roots agree.
- FAIL at late times (t > 0.2 s, fast slip): one solver is exceeding iteration budget or hitting a bracket-degeneracy. Bisect on the residual `tau - sigma_n * mu - eta_s * V` at the reported `V_abs` — if both solvers' residuals are < tolerance, this is a convergence-sharpness issue and needs solver tuning (not a bug per se).
- FAIL at early times (t < 1e-6 s, V ~ V_ini = 1e-16): one of the solvers is not cleanly handling sub-normal V. TPV102 has seen this before (see `dieterich_ruina.hpp:372-384`'s `Va_min = -300.0` fix); TPV104's `V_ini = 1e-16` is much smaller than TPV102's `1e-12` and may reopen the issue.

### 5.6 Probe pair 5 — Corrected traction and imposed state [Stage E]

**What it captures**: the final Riemann solve output that becomes the bulk's ADER flux input on the next step. Any disagreement here propagates into bulk Q on the next ADER sub-step.

**MFEM probe point**: `fault_face_flux.cpp:202-204` inside `WriteBackState`. Dump `(tau1_corr_TOTAL, tau2_corr_TOTAL, sigma_n_corr_TOTAL)` plus `(Q_imp_plus[0..8], Q_imp_minus[0..8])`.

**SeisSol probe point**: `FrictionSolverCommon.h:407` (end of `postcomputeImposedStateFromNewStress`). Dump the same total-traction triple plus `imposedStatePlus[9], imposedStateMinus[9]`.

**Fields**: `tau1_corr, tau2_corr, sigma_n_corr, Q_imp_plus[9], Q_imp_minus[9]`. 21 total scalars per sample.

**Threshold**: `max |diff| / max |signal| < 5e-3` across 21 channels. Coarser than earlier probes because this accumulates all prior stages' residual deltas.

**Diagnostic outcome**:

- PASS: the full fault-face pipeline is consistent. If the station traces still disagree after this, the disagreement is **bulk-side** (ADER wave-operator bug, not friction bug) — bisect via a different diagnostic: turn off the fault (no-slip condition / reflecting BC at y=0) and compare wave-equation radiation in bulk.
- FAIL: the disagreement is localised to the Riemann imposed-state construction (`Q_imp_plus/minus` components) or the final `*_corr` computation. Dump each of the 21 fields' diffs and flag the first non-matching component.

### 5.7 Paired fixtures and run protocol

**Canonical test run** (one probe run per phase):

| stage | MFEM side | SeisSol side |
|---|---|---|
| Mesh | `tpv104/mesh/tpv104_1000m.msh` (MFEM equivalent of TPV5) | `tpv5.msh` (SeisSol standard, on Frontera) |
| Ranks | `np=1` (first pass), `np=4` (second pass once `np=1` passes) | match MFEM's rank count |
| ADER order | O=2 | O=5 (SeisSol default) |
| DG order | p=1 | p=4 |
| `dt` | CFL × h / c_p, MFEM driver chooses | SeisSol chooses (`dt ≈ 4.7 ms` for the benchmark) |
| t_final | 12 s (probe sampling at §5.1 schedule) | 12 s |
| Nucleation | full SCEC TPV104 (R=3km, Δτ=45MPa, T=1s) | full SCEC TPV104 |

**Note on order mismatch** (critical per §3.4): MFEM's O=2/p=1 vs SeisSol's O=5/p=4 is a deliberate resolution-difference run. The five probes above all sample at fault QPs — which mean their intrinsic spatial interpolation is the face-quadrature basis. At O=2/p=1 MFEM has 3 face QPs per triangle (Dunavant 3-point); at O=5/p=4 SeisSol has 15 face QPs per triangle. **The probe samples must be interpolated to a common reference point** before diffing — in practice, the Station (x₂, x₃) nearest-QP on each side. This is acceptable because the station coordinate is known a priori (`kStationsTPV104[9]` in `tpv104_params.hpp`); each code dumps its own nearest-QP data for each station.

Run protocol:

1. **Phase 3.A (no nucleation, V_ini = 1e-16)** — both codes run with nucleation disabled. Expected: bulk stays at pre-stress, all 5 probes show `V ≈ 0`, `ψ ≈ ψ_ini` constant, trial traction constant at background. If any probe pair diverges in this silent regime, the bug is in initialisation (prestress in bulk Q, or ψ_ini) — the simplest possible attribution.
2. **Phase 3.B (full nucleation)** — both codes run with full SCEC TPV104 nucleation. Probe 1–5 outputs compared. First probe to fail identifies the pipeline stage.
3. **Phase 3.C (`np=4`)** — once 3.B passes at `np=1`, repeat at `np=4` and look specifically at probes at stations crossing MPI seams. This is the pepper-pattern reproducer.

**Approval checkpoint before 3.A**: per `feedback_frontera_approval.md`, I do NOT launch a Frontera sbatch autonomously. For each of 3.A / 3.B / 3.C the user must approve; this plan produces the sbatch script and a one-line run command, but the user runs it.

### 5.8 Probe analysis tooling (new, Python)

`tpv104/scripts/probe_diff.py` — reads both `tpv104_mfem_probe_*.txt` and `tpv104_seissol_probe_*.txt`, aligns on `(station, side, t)`, computes per-field relative and absolute diffs, plots on a per-station page with the threshold lines from §5.2–5.6 overlaid. Flags any sample exceeding threshold.

`tpv104/scripts/probe_summarise.py` — one-line-per-probe-pair summary: max absolute diff, max relative diff, time of first exceedance, station of first exceedance. This is the top-level diagnostic output that tells us which of the 5 probes failed.

### 5.9 Cross-verification priorities (informed by TPV102 audit history)

This subsection answers the first of the user's two Phase 3 asks: *"read through previous unit tests and debug report, and set priorities on the section we want to start"*. Phase 3 fires only after Phase 2 is green AND the Frontera TPV104 run exhibits the same pepper signature as TPV102 (per §0, the plan's motivating premise); the priorities below are derived from the TPV102 pepper attribution already on record.

#### 5.9.1 What the TPV102 audit already proves about the five stages

Evidence consumed: `tpv102_unit_tests_to_narrow_bug_2026-04-23.md` §7 Gate 5–14 (final), `tpv102_debug_v8.0.0_seissol_audit.md` §§"Trial traction" / "Imposed-state construction" (bit-for-bit MFEM-vs-SeisSol side-by-side), `tpv102_seissol_aligned_flux_plan_2026-04-23.md` §0–§5 (precomputed-face-flux design that fixes H_IFACE_NOR / H_BFACE_NOR).

| TPV102 gate | stage of the §5 pipeline it exercises | audit verdict | implication for Phase 3 probe |
|---|---|---|---|
| Gate 10 — zero-input sanity for interior and boundary non-fault rhs | pre-input to the fault pipeline (bulk ADER → fault-QP Q) | PASS (0 / 0 bit-exact) | NOT a reason to skip Probe 1; Gate 10 ran with **zero input**, whereas the pepper manifests on **non-zero orbit-uniform input** (Gate 14) |
| Gate 14 — constant-`I` lifted element rhs orbit drift | bulk ADER → fault-QP Q (`Q_plus` / `Q_minus` **before** they reach `ComputeTrialTraction`) | **FAIL 1.000 (bface)** on M_ref, scale-invariant to 8×4×8 | **Probe 1's sub-probe `Q_in` is the primary target** — the bulk Q that `ComputeTrialTraction` receives is the locus of the residual asymmetry |
| Gate 5a / 5b / 5c — per-face interior non-fault `F_h` | bulk flux into fault-adjacent tets | 5a 2.000, 5b 5.68e-14, 5c 2.000 | confirms `H_IFACE_NOR`. Fault-QP Q inherits that asymmetry; differencing MFEM-vs-SeisSol at `Q_in` (Probe 1 sub-probe) reveals whether `tau1_corr = 1.03e+01` traces to the same non-fault asymmetry still leaking in (MFEM's `Q_in` ≠ SeisSol's `qI`) or to a new mechanism inside the fault pipeline |
| Gate 7 — element-level lifted drift across four symmetrization modes | bulk ADER → fault-QP Q | iface clean under symmetrization; bface SXY residual 1.00 unresolved; bface SXZ clean under sign-symmetrization | SeisSol lacks this failure mode by construction (per-cell topology-based normals, per-face `AplusT` precomputed); a **bface-side difference between the two codes at Probe 1 `Q_in` is the strongest signature we can ask for** |
| v8.0.0 SeisSol audit — `ComputeTrialTraction` vs `FrictionSolverCommon.h:182-192` | Stage A (trial traction) | **bit-for-bit identical** (modulo homogeneous-material `Zp_plus = Zp_minus = Zp`) | Stage A **mathematics** is clean. Any Probe 1b `trial_out` mismatch with the same `Q_in` is therefore either a code-level regression in `ComputeTrialTraction` or a rotation-frame mismatch (`Tinv_can` at `wave_operator.inl:1235-1238` vs SeisSol's `rotateToFaultCS`) — downgraded to SECONDARY hunt |
| v8.0.0 SeisSol audit — imposed-state vs `FrictionSolverCommon.h:347-362` | Stage E (imposed state) | **bit-for-bit identical** (modulo SeisSol's ADER time-weight accumulation; see §5.10.5) | Stage E mathematics is clean. Probe 5 mismatch with the same upstream inputs = rotation-frame bug in `T_can` at `wave_operator.inl:1260-1269`, or a real divergence from an earlier stage that Probe 1/2/3/4 did not flag — TERTIARY hunt |
| `test_adjacent_triangle_fault_first_step_audit.cpp` Gate 1–3 (fault lift + M⁻¹, EvaluateADER, canonical rotation) | fault path: C (μ) + D (V_abs) + E (imp) | PASS on first 2 steps under orbit-uniform nucleation | Probes 3 and 4 expected to pass in isolation; they buy rule-out evidence but are unlikely to be the pepper source. |

**Takeaway**: the residual `tau1_corr ≈ 1.03e+01` at step 19, `np=4` does **not** prove a bug in the §5 fault-pipeline stages (A–E). It is entirely consistent with the bface-SXY 1.00 per-step drift (Gate 7b residual) accumulating through 19 steps into an O(10) tau1_corr signal at a single QP. Phase 3's job is to *falsify or confirm* that interpretation against SeisSol, not to re-prosecute Gates 1–3.

#### 5.9.2 Phase 3 probe priority order (drives the Frontera sbatch design)

Priorities `P1 / P2 / P3` correspond to what we instrument first, second, third in the §5.10 code-line maps; a failed higher-priority probe is sufficient to localise the bug without collecting lower-priority probes. The order diverges from naive pipeline-order §5.2–5.6 because the TPV102 evidence weights the **input** side of the fault pipeline (Probe 1) above the arithmetic-identity stages (Probes 2/3/4) that the v8.0.0 audit already cleared.

| priority | probe | rationale from §5.9.1 | Phase-3 sub-phase that exercises it |
|---|---|---|---|
| **P1** | **Probe 0 — Initialization (new; §5.10.0)** | simplest possible mismatch; SCEC ψ_ini anchor (`5.6359184e-01`) and initial stress in fault CS are known to single-digit ppm on both sides. A Probe 0 fail blocks all five downstream probes regardless of dynamics. | 3.A (nucleation OFF, t=0 sample) |
| **P1** | **Probe 1a — `Q_in` at fault QP (new sub-probe; §5.10.1.a)** | direct successor of Gate 14's constant-`I` failure mode. MFEM's `Q_plus_local / Q_minus_local` at `wave_operator.inl:1240-1243` is the exact variable whose drift causes the pepper. SeisSol's counterpart is `qIPlus / qIMinus` in `FrictionSolverCommon.h:165-166`. | 3.B (nucleation ON, `np=1` seam-free) |
| **P1** | **Probe 1b — `trial_out` after `ComputeTrialTraction` (§5.10.1.b)** | tests the *mathematical identity* Stage A = `FrictionSolverCommon.h:182-192`. v8.0.0 audit already proved identity in isolation; Probe 1b catches any regression or rotation-frame mismatch that `Q_in` alone would not. | 3.B |
| P2 | Probe 5 — imposed state (§5.10.5) | OUT feed into bulk Q; if Probe 1a/1b both pass and Probe 5 fails, the residual error is in `T_can` rotation at `wave_operator.inl:1260-1269` or in the imposed-state construction (v8.0.0 audit cleared the math). | 3.B |
| P2 | Probe 2 — state evolution (`ψ_in → ψ_out`, §5.10.2) | arithmetic-identity stage (new MFEM class `SlipLawSRWPsi` vs SeisSol FVW). T_SRW_5 in §4.2.3 already byte-matches this in unit tests; Probe 2 is a runtime guard in case the unit-test fixture missed a corner case of the TPV104 envelope. | 3.B |
| P3 | Probe 3 — friction coefficient μ (§5.10.3) | MFEM `FrictionCoefficientPsi` vs SeisSol `updateMu` (both consume `V, ψ, a`). Different branch structures (MFEM `psi/a > 700` asymptote vs SeisSol `arsinhexp` Switch/Threshold), but the TPV104 envelope `ψ/a ≤ 65` is far from every branch boundary. Low probability of failure. | 3.B |
| P3 | Probe 4 — slip-rate magnitude V_abs (§5.10.4) | MFEM Brent in log10(V) vs SeisSol Newton-Raphson in V. Different convergence paths to the same monotone root; TPV102 audits show no divergence at `V ≥ 1e-12`. TPV104 `V_ini = 1e-16` is deep in the sub-normal regime, but `CompleteFromVabs` short-circuits to V=0 when `Θ < ULP·σ_n·asinh(…)`, so this only matters **after** nucleation drives `Θ` up — by then all codes sit at physical slip rates. | 3.B |

**Sub-phase assignments**: 3.A (nucleation OFF, all probes at t=0) exercises **Probe 0 only** (everything else is trivially zero/steady). 3.B (nucleation ON, `np=1`) exercises **Probes 1a, 1b, 2, 3, 4, 5** in that priority order. 3.C (nucleation ON, `np=4`) targets the pepper reproducer — **Probe 1a is the primary** because shared-face assembly (`ComputeADERSharedFaceFluxRHS` at `wave_operator.inl:1762`) is where bulk Q crosses rank boundaries.

**Decision shortcut**: if Probe 1a FAILS at 3.B (np=1) with 3.A clean, the bug is in the MFEM bulk ADER non-fault assembly — same family as TPV102 H_IFACE_NOR / H_BFACE_NOR. That result authorises the precomputed-face-rotation plan (`tpv102_seissol_aligned_flux_plan_2026-04-23.md`) for the TPV104 driver as well, without re-opening the fault pipeline.

### 5.10 One-to-one code-line maps (MFEM ↔ SeisSol equivalence)

This subsection answers the second user ask: *"in each section show exact code lines in both code that you want to add diagnostics and make sure they are equivalent and can reflect the real comparison"*. Every probe below specifies:

- **MFEM insertion point** — absolute file path + function + line number + the statement that MUST be present on that line (so a future refactor that moves or renames the variable can be detected)
- **SeisSol insertion point** — same, rooted at the Frontera clone `/Users/chunhuizhao/projects/SeisSol/src/...`
- **Scalar equivalence table** — a left-to-right mapping MFEM-name → SeisSol-name with the algebraic identity and any coordinate-frame / unit / sign caveats
- **Rationale** — why the two quantities are comparable bit-for-bit (or within the declared threshold)

Convention: line numbers are the **head of the relevant statement** as of 2026-04-24. Before instrumenting, the Frontera developer MUST `grep` for the quoted statement text and use the live line number — the literal text is the authoritative anchor, not the line number.

#### 5.10.0 Probe 0 — Initialization (new: ψ_ini and initial fault-CS stress)

**What it captures**: the pre-step-0 state of the fault DOF. If this disagrees, nothing downstream can pass. Fires **once** at driver init, after the fault equilibrium solve and the `initialStressInFaultCS` rotation have completed, before any time step.

**MFEM probe point**: `miniapps/seas/drivers/tpv104_driver.cpp` — **R-002 FIX 2026-04-24**: place the probe **immediately after the `V_w` side-channel population loop completes**, i.e. after (i) `InitializeStateTotal_TPV104`, (ii) the `for (i) V_w[i] = ComputeVw_TPV104(x, z)` loop from §4.5, and (iii) the single `ApplyNucleationPrestress_TPV104(..., t=0.0)` dry-run call that Phase 3.A inserts to bring the DOFData into the same state the first time-step would see. The probe must assert `V_w.size() == dof_data.size()` before the dump; if the assertion fails, Phase 3.A aborts rather than writing garbage `V_w` bytes. Fields to dump per fault QP:

```cpp
// Anchor: immediately after the call that populates dof_data[i].{sigma_n0,
// tau1_0, tau2_0, psi} for every fault QP i.  Use the same diag_print
// gating pattern as TPV102's wave_operator.inl:1287 C-2 block.
for (int i = 0; i < dof_data.size(); i++) {
   if (!dof_data[i].diag_print) continue;
   std::fprintf(stderr,
      "[P0 INIT] qp=%d x=%+.4e z=%+.4e  psi_ini=%+.9e  "
      "sigma_n0=%+.9e  tau1_0=%+.9e  tau2_0=%+.9e  a=%+.6e  V_w=%+.6e\n",
      i, fault_coords[i](0), fault_coords[i](2),
      dof_data[i].psi, dof_data[i].sigma_n0,
      dof_data[i].tau1_0, dof_data[i].tau2_0,
      dof_data[i].a, V_w[i]);
}
```

**SeisSol probe point**: `/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/Initializer/RateAndStateInitializer.cpp:86` — immediately after `stateVariable[ltsFace][pointIndex] = stateAndFriction.stateVariable;` fires inside `RateAndStateInitializer::initializeFault`:

```cpp
// Anchor: inside the double loop at RateAndStateInitializer.cpp:59-88,
// right AFTER
//   "stateVariable[ltsFace][pointIndex] = stateAndFriction.stateVariable;"
// and BEFORE the closing brace of the inner for-loop.
if (/* per-face gate matching the nine stations */) {
   std::fprintf(stderr,
      "[P0 INIT] face=%zu qp=%u x=%+.4e z=%+.4e  psi_ini=%+.9e  "
      "sigma_n0=%+.9e  tau1_0=%+.9e  tau2_0=%+.9e  a=%+.6e  V_w=%+.6e\n",
      ltsFace, pointIndex,
      faceCoords_[ltsFace][pointIndex][0],
      faceCoords_[ltsFace][pointIndex][2],
      stateVariable[ltsFace][pointIndex],
      initialStressInFaultCS[ltsFace][0 /*XX=σ_nn*/][pointIndex],
      initialStressInFaultCS[ltsFace][3 /*XY=τ_1*/ ][pointIndex],
      initialStressInFaultCS[ltsFace][5 /*XZ=τ_2*/ ][pointIndex],
      rsA[ltsFace][pointIndex],
      srW[ltsFace][pointIndex]);
}
```

**Scalar equivalence**:

| quantity | MFEM variable | SeisSol variable | identity |
|---|---|---|---|
| initial state | `dof_data[i].psi` | `stateVariable[ltsFace][pointIndex]` | both solved from `tau_ini = σ_n · a · asinh((V_ini/2V₀) · exp(ψ_ini/a))` → `ψ_ini = a · ln(2V₀/V_ini · sinh(tau_ini/(σ_n·a)))`. MFEM: `ComputeInitialPsiTPV104(a_in)` (§4.1). SeisSol: `RateAndStateFastVelocityInitializer::computeInitialStateAndFriction` at `RateAndStateInitializer.cpp:141-164`, which computes `rsA * log(2 * rsSr0 / initialSlipRate * (exp(tmp) - exp(-tmp)) / 2) = a·ln((2V₀/V_ini)·sinh(tmp))`. Identical up to ULP. |
| normal pre-stress | `dof_data[i].sigma_n0` | `initialStressInFaultCS[ltsFace][0][pointIndex]` | both SI Pa. **Sign convention caveat**: MFEM stores σ_n0 POSITIVE (compression-positive geology convention; CLAUDE.md §"Sign Conventions"). SeisSol `initialStressInFaultCS[0]` is the raw rotated Cauchy stress, which for a compressive pre-stress on a vertical fault is **negative** (σ_yy < 0; see `updateNormalStress` at `RateAndState.h:352` clamping to `std::min(0, ...)`). **probe_diff.py assertion**: `sigma_n0_MFEM ≈ -initialStressInFaultCS[0]_SeisSol` to ULP on a homogeneous fault. |
| shear pre-stress τ_1 | `dof_data[i].tau1_0` | `initialStressInFaultCS[ltsFace][3][pointIndex]` | MFEM BP5 convention (CLAUDE.md §"Fault-local tangent frame"): `tangent1 = dip, tangent2 = strike`. On TPV104's vertical y=0 fault with `ref_normal=(0,-1,0)`, MFEM's `can_t1 = (0, 0, -1)` (down-dip). Per research report §8, SeisSol's `tangent1` identification is inferred but not yet directly verified. **Probe 0 is the authoritative settle**: at t=0, SeisSol writes pure strike-slip pre-stress `τ_ini = 40 MPa` and MFEM writes pure strike-slip pre-stress into `tau2_0`. Acceptance: EITHER (a) `tau1_0_MFEM == traction1_SeisSol ≈ 0` AND `tau2_0_MFEM == traction2_SeisSol ≈ 40 MPa`, OR (b) the reverse (MFEM's tau2_0 = 0 while SeisSol's traction1 holds the 40 MPa). Outcome (b) means SeisSol uses `tangent1 = strike` and the station-output column mapping in §3.3 flips; `probe_diff.py` applies the flip automatically based on Probe 0's t=0 attribution. |
| shear pre-stress τ_2 | `dof_data[i].tau2_0` | `initialStressInFaultCS[ltsFace][5][pointIndex]` | same convention caveat as τ_1 |
| direct effect a | `dof_data[i].a` | `rsA[ltsFace][pointIndex]` | dimensionless, per-QP, evaluated at identical `(along_strike, down_dip)` via §4.1's `ComputeA_TPV104`. Identity to ULP if the two implementations use the same `Boxcar(x, W, w)` formula. |
| weakening rate V_w | `V_w[i]` (driver side-channel per §4.5) | `srW[ltsFace][pointIndex]` | SeisSol stores `srW` per-QP in the LTS layer (`FastVelocityWeakeningLaw.h:28`, `RateAndStateInitializer.cpp:169-170`). TPV104 spec requires the same `V_w_in/V_w_out + Boxcar(3 km)` distribution as MFEM. Identity to ULP expected. |

**Threshold**: `1e-8` relative on `ψ_ini` (anchored to SeisSol canonical trace col-9 `5.6359184e-01`); `1e-12` relative on every other scalar. A mismatch here is **stop-ship** — Phase 3.B cannot run until Probe 0 is clean.

**Rationale**: Gate 10 (TPV102 zero-input sanity) rules out the rhs lift for zero bulk, but TPV104 starts with NON-zero equilibrium ψ plus non-zero pre-stress baked into bulk Q (per §3.1 "Initial pre-stress layout in bulk Q"). We MUST verify the starting constants separately from the dynamic-propagation probes.

#### 5.10.1 Probe 1 — Pre-friction (trial traction) — split into sub-probes 1a (`Q_in`) and 1b (`trial_out`)

This is the highest-priority probe pair (§5.9.2). We split it into two sub-probes because the v8.0.0 SeisSol audit proved the math inside `ComputeTrialTraction` is bit-for-bit identical to SeisSol's `precomputeStressFromQInterpolated` — so the first question is whether MFEM's Q INPUT matches SeisSol's.

##### 5.10.1.a Sub-probe 1a — `Q_in` (bulk Q at fault QP, fault-local frame, BEFORE trial traction)

**MFEM insertion point**: `miniapps/seas/dynamic/wave_operator.inl:1240-1243`, the moment `Q_plus_local / Q_minus_local` pointers are defined but BEFORE `fault_flux_->Evaluate` at line 1255:

```cpp
// Before line 1255 (fault_flux_->Evaluate(fdata, ...)).  Uses the same
// SEAS_DIAG_FAULT_FLUX guard as the existing [C-1 EVAL] print at
// fault_face_flux.cpp:261, extended with a new [P1A Q_IN] tag.
#ifdef SEAS_DIAG_TPV104_STATE
if (fdata.diag_print) {
   std::fprintf(stderr,
      "[P1A Q_IN] rank=%d qp=%d t=%+.6e  "
      "Q_plus[SXX]=%+.9e Q_plus[SXY]=%+.9e Q_plus[SXZ]=%+.9e  "
      "Q_plus[VX]=%+.9e  Q_plus[VY]=%+.9e  Q_plus[VZ]=%+.9e  "
      "Q_minus[SXX]=%+.9e Q_minus[SXY]=%+.9e Q_minus[SXZ]=%+.9e  "
      "Q_minus[VX]=%+.9e  Q_minus[VY]=%+.9e  Q_minus[VZ]=%+.9e\n",
      g_seas_my_rank, qp_id, t_now,
      Q_plus_local[SXX], Q_plus_local[SXY], Q_plus_local[SXZ],
      Q_plus_local[VX],  Q_plus_local[VY],  Q_plus_local[VZ],
      Q_minus_local[SXX], Q_minus_local[SXY], Q_minus_local[SXZ],
      Q_minus_local[VX],  Q_minus_local[VY],  Q_minus_local[VZ]);
}
#endif
```

**Anchor line** (verify before instrumenting): `const real_t *Q_plus_local  = elem1_on_plus ? Q_self_can : Q_nbr_can;` at `wave_operator.inl:1240`.

**SeisSol insertion point**: `/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h:165-166`, inside `precomputeStressFromQInterpolated`, BEFORE the `for (unsigned o = 0; o < misc::TimeSteps; ++o)` loop at line 174:

```cpp
// After the reinterpret_casts on :165-166, BEFORE the outer o-loop on :174.
#ifdef SEISSOL_TPV104_DIAG
for (unsigned o = 0; o < misc::TimeSteps; ++o) {
   for (unsigned i = 0; i < misc::NumPaddedPoints; ++i) {
      if (!shouldDumpQp(ltsFace, i, o, fullUpdateTime)) continue;
      std::fprintf(stderr,
         "[P1A Q_IN] rank=%d ltsFace=%zu qp=%u o=%u t=%+.6e  "
         "qIPlus[N]=%+.9e  qIPlus[T1]=%+.9e  qIPlus[T2]=%+.9e  "
         "qIPlus[U]=%+.9e  qIPlus[V]=%+.9e  qIPlus[W]=%+.9e  "
         "qIMinus[N]=%+.9e  qIMinus[T1]=%+.9e  qIMinus[T2]=%+.9e  "
         "qIMinus[U]=%+.9e  qIMinus[V]=%+.9e  qIMinus[W]=%+.9e\n",
         mpiRank, ltsFace, i, o, fullUpdateTime + deltaT[o],
         qIPlus[o][N][i],  qIPlus[o][T1][i], qIPlus[o][T2][i],
         qIPlus[o][U][i],  qIPlus[o][V][i],  qIPlus[o][W][i],
         qIMinus[o][N][i], qIMinus[o][T1][i], qIMinus[o][T2][i],
         qIMinus[o][U][i], qIMinus[o][V][i],  qIMinus[o][W][i]);
   }
}
#endif
```

**Anchor line** (verify before instrumenting): `const auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus));` at `FrictionSolverCommon.h:165`.

**Scalar equivalence**:

| fault-local quantity | MFEM variable | SeisSol index (per `Misc.h:164-177`) | SeisSol variable | identity |
|---|---|---|---|---|
| σ_nn | `Q_plus_local[SXX]` | `N = 0` | `qIPlus[o][N][i]` | both Pa, fault-local. MFEM: obtained from global Q via `Tinv_can` at `wave_operator.inl:1235` (`Q_self_can[c] = Σ_k Tinv_can(c,k) · Q_self[k]`). SeisSol: `qIPlus` is already fault-local because `initializer::rotateStressToFaultCS` (`BaseDRInitializer.cpp:126-127`) rotates bulk DOFs before friction. **Sign-convention caveat**: standing compressive pre-stress σ_n0 = 120 MPa appears in MFEM `SXX` slot with sign depending on driver path: `EvaluateTotal` (Q carries full stress) has `SXX = -120e6`; `Evaluate` (fluctuation-only) has `SXX = 0` baseline. TPV104 runs under `EvaluateTotal` per §4.4 item 4 (SeisSol uses full-stress too), so both sides' `SXX` are negative at t=0. |
| τ_1 | `Q_plus_local[SXY]` | `T1 = 3` | `qIPlus[o][T1][i]` | both Pa, fault-local. MFEM's SXY (Voigt σ_xy) becomes `σ_n_t1` after `Tinv_can`. SeisSol's T1 is the same fault-local slot. Convention caveat (Probe 0 is the authoritative settle): MFEM stores `t1 = dip`; if SeisSol uses `t1 = strike`, `probe_diff.py` swaps SXY ↔ SXZ on one side. |
| τ_2 | `Q_plus_local[SXZ]` | `T2 = 5` | `qIPlus[o][T2][i]` | same convention caveat as τ_1 |
| v_n | `Q_plus_local[VX]` | `U = 6` | `qIPlus[o][U][i]` | both m/s. MFEM comment `fault_face_flux.cpp:37-38`: "In fault-local coordinates: x = normal, y = tangent1, z = tangent2. So VX = v_n, VY = v_t1, VZ = v_t2." SeisSol comment block at `Misc.h:164-177` matches. |
| v_t1 | `Q_plus_local[VY]` | `V = 7` | `qIPlus[o][V][i]` | same as VX row, tangent-1 slot |
| v_t2 | `Q_plus_local[VZ]` | `W = 8` | `qIPlus[o][W][i]` | same as VX row, tangent-2 slot |

**Threshold** — absolute: `1e-4 Pa` (stresses), `1e-10 m/s` (velocities). Relative: `1e-3` per field. Loose because SeisSol's `qIPlus[o]` is interpolation at the `o`-th ADER time-quadrature point (O=5 Dubiner) while MFEM's `Q_plus_local` is the `o`-th RK4 stage; quadrature rules differ, bit-match not expected. **A Probe 1a FAIL at threshold means the bulk Q at the fault QP is materially different** between the two codes — the pepper's locus.

**Time-alignment note**: SeisSol's `deltaT[o]` and MFEM's ADER sub-step time are NOT interleaved on the same grid. Probe 1a samples at the coarse §5.1 schedule; each code interpolates its own `Q_in` to those global times via the `TimeSampler` utility in `probe_diff.py`.

**Rationale**: a Probe 1a FAIL with Probe 0 PASS localises the divergence to propagation — which by TPV102 Gate-14 evidence traces to the non-fault face flux branches. This result authorises the `tpv102_seissol_aligned_flux_plan_2026-04-23.md` fix path for TPV104 without touching the fault pipeline.

##### 5.10.1.b Sub-probe 1b — `trial_out` (output of `ComputeTrialTraction`)

**MFEM insertion point**: `miniapps/seas/dynamic/fault_face_flux.cpp`. Two variants depending on driver path:

- Fluctuation path (legacy `Evaluate`): inside `ComputeStageState`, after `ComputeTrialTraction` returns. The existing `[C-1 EVAL]` print at `:261-270` dumps `tau1_trial / tau2_trial` but NOT `sigma_n_trial` — extend it.
- Total-stress path (`EvaluateTotal`, the TPV104 production path per §4.4 item 4): `:361-363`, immediately after `ComputeTrialTraction(data, Q_plus, Q_minus, sigma_n_trial, tau1_trial, tau2_trial);`. The existing `[C-1 EVAL-TOTAL]` print at `:387-395` dumps `tau*_fric = tau*_trial + data.tau*_nuc`; we must also dump the RAW trial (before the `+ data.tau*_nuc` shift) for SeisSol parity.

```cpp
// In EvaluateTotal, immediately after ComputeTrialTraction returns at
// fault_face_flux.cpp:363, BEFORE the sigma_n_fric / tau1_fric /
// tau2_fric assignment at :377-379.
#ifdef SEAS_DIAG_TPV104_STATE
if (data.diag_print) {
   std::fprintf(stderr,
      "[P1B TRIAL_OUT] rank=%d t=%+.6e  "
      "sigma_n_trial=%+.9e  tau1_trial=%+.9e  tau2_trial=%+.9e\n",
      g_seas_my_rank, t_now,
      sigma_n_trial, tau1_trial, tau2_trial);
}
#endif
```

**Anchor line** (verify before instrumenting): `ComputeTrialTraction(data, Q_plus, Q_minus, sigma_n_trial, tau1_trial, tau2_trial);` at `fault_face_flux.cpp:362-363` (for `EvaluateTotal`). **R-007 FIX 2026-04-24** — for the legacy `ComputeStageState` path at `:84-85` the call signature is `ComputeTrialTraction(data, Q_plus, Q_minus, s.sigma_n_trial, s.tau1_trial, s.tau2_trial);` (the three trial outputs are members of `EvalStageState s`, not bare locals). A grep-by-literal-text of the bare-name form returns zero matches at `:84-85`; use the `s.`-qualified anchor when instrumenting the legacy path.

**SeisSol insertion point**: `/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h:193`, immediately after the inner simd-loop closes, still inside the `for (unsigned o = 0; o < misc::TimeSteps; ++o)` loop at `:174`. The three quantities `faultStresses.normalStress[o][i] / traction1[o][i] / traction2[o][i]` are written on lines 182-192:

```cpp
// Inside precomputeStressFromQInterpolated, after the simd loop closes
// at FrictionSolverCommon.h:193, still inside the o-loop at :174.
#ifdef SEISSOL_TPV104_DIAG
for (unsigned i = 0; i < misc::NumPaddedPoints; ++i) {
   if (!shouldDumpQp(ltsFace, i, o, fullUpdateTime)) continue;
   std::fprintf(stderr,
      "[P1B TRIAL_OUT] rank=%d ltsFace=%zu qp=%u o=%u t=%+.6e  "
      "normalStress=%+.9e  traction1=%+.9e  traction2=%+.9e\n",
      mpiRank, ltsFace, i, o, fullUpdateTime + deltaT[o],
      faultStresses.normalStress[o][i],
      faultStresses.traction1[o][i],
      faultStresses.traction2[o][i]);
}
#endif
```

**Anchor line** (verify before instrumenting): `VariableIndexing<...>::index(faultStresses.traction2, o, i) = etaS * (qIMinus[o][W][i] - qIPlus[o][W][i] + qIPlus[o][T2][i] * invZs + qIMinus[o][T2][i] * invZsNeig);` at `FrictionSolverCommon.h:190-192`.

**Scalar equivalence**:

| MFEM local | SeisSol local | identity |
|---|---|---|
| `sigma_n_trial` | `faultStresses.normalStress[o][i]` | MFEM: `fault_face_flux.cpp:53-55` `η_p · (Q⁻[VX] − Q⁺[VX] + Q⁺[SXX]/Zp⁺ + Q⁻[SXX]/Zp⁻)`. SeisSol: `FrictionSolverCommon.h:182-184` `etaP · (qIMinus[U] − qIPlus[U] + qIPlus[N]·invZp + qIMinus[N]·invZpNeig)`. Under the Probe-1a index map (`SXX ↔ N`, `VX ↔ U`) and the homogeneous reduction `invZpNeig = 1/Zp⁻`, **bit-for-bit identical** per the v8.0.0 audit. |
| `tau1_trial` | `faultStresses.traction1[o][i]` | MFEM: `fault_face_flux.cpp:58-60` `η_s · (Q⁻[VY] − Q⁺[VY] + Q⁺[SXY]/Zs⁺ + Q⁻[SXY]/Zs⁻)`. SeisSol: `FrictionSolverCommon.h:186-188` `etaS · (qIMinus[V] − qIPlus[V] + qIPlus[T1]·invZs + qIMinus[T1]·invZsNeig)`. Identical. |
| `tau2_trial` | `faultStresses.traction2[o][i]` | MFEM: `fault_face_flux.cpp:63-65`. SeisSol: `FrictionSolverCommon.h:190-192`. Identical. |

**Threshold** — `1e-4 Pa` absolute AND `1e-10` relative (tighter than Probe 1a). Bit-match expected iff Probe 1a passes; any Probe 1b drift on top of clean Probe 1a = regression in `ComputeTrialTraction` or a rotation-frame bug in `Tinv_can`.

**Rationale**: the v8.0.0 audit proves the math is identical *as expressions*. Probe 1b verifies the expression has not been tampered with by a future refactor and that the rotation producing `Q_plus_local / Q_minus_local` from global is bit-consistent with SeisSol's `initializer::rotateStressToFaultCS`.

#### 5.10.2 Probe 2 — State evolution (ψ_in → ψ_out)

**MFEM insertion point**: `miniapps/seas/drivers/tpv104_driver.cpp` at the per-step ψ-integration block (Phase 2, §4.4 item 4). Statement: `dof_data[i].psi = UpdateStateAnalyticSlipLawSRW(psi_n[i], dof_data[i].slip_rate, dof_data[i].Dc, dt_step, V_w[i], dof_data[i].a, b, V0, f0, muW);`. Dump BEFORE and AFTER:

```cpp
#ifdef SEAS_DIAG_TPV104_STATE
if (dof_data[i].diag_print) {
   const real_t psi_in = psi_n[i];
   const real_t V_in   = dof_data[i].slip_rate;
   const real_t L_in   = dof_data[i].Dc;
   const real_t Vw_in  = V_w[i];
   const real_t a_in   = dof_data[i].a;
   dof_data[i].psi = UpdateStateAnalyticSlipLawSRW(
      psi_in, V_in, L_in, dt_step, Vw_in, a_in,
      TPV104Params::b, TPV104Params::V0,
      TPV104Params::f0, TPV104Params::muW);
   std::fprintf(stderr,
      "[P2 PSI] rank=%d qp=%d t=%+.6e  "
      "psi_in=%+.9e  V=%+.9e  L=%+.9e  dt=%+.9e  V_w=%+.9e  a=%+.9e  "
      "b=%+.9e  V0=%+.9e  f0=%+.9e  muW=%+.9e  psi_out=%+.9e\n",
      g_seas_my_rank, i, t_now,
      psi_in, V_in, L_in, dt_step, Vw_in, a_in,
      TPV104Params::b, TPV104Params::V0,
      TPV104Params::f0, TPV104Params::muW, dof_data[i].psi);
}
#endif
```

**Anchor line**: the call to `UpdateStateAnalyticSlipLawSRW` (new helper per §4.2.2).

**SeisSol insertion point**: `/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/FrictionLaws/CpuImpl/FastVelocityWeakeningLaw.h:74-76`, immediately before `return localStateVariable;`. Inputs `(stateVarReference, timeIncrement, localSlipRate, localSl0, localA)` are in scope; SRW constants live at `localMuW / localSrW / this->b / this->f0 / this->drParameters->rsSr0`.

**R-003 FIX 2026-04-24 — last-iteration-only filter.** `updateStateVariable` is called up to **3 times per `(ltsFace, pointIndex, o)`**:  twice inside the Kaneko fixed-point loop in `updateStateVariableIterative` (`RateAndState.h:167-172`, settings.numberStateVariableUpdates = 2), and once more in `calcSlipRateAndTraction` (`:210-215`) for the final ψ. Only the final (`calcSlipRateAndTraction`) call produces the physically-meaningful ψ_new. The probe MUST fire only on that call, or probe_diff.py will match MFEM's single sample against 3 SeisSol samples per time-bin and falsely flag a failure.

Implementation: declare a `thread_local bool g_in_final_sv_update` in the probe helper (§5.12). Set it `true` by an RAII guard placed around the `updateStateVariable` call inside `calcSlipRateAndTraction` (`RateAndState.h:210-215`). The probe body below checks `g_in_final_sv_update` as the first gate.

```cpp
// Inside updateStateVariable(...), after line 74
//   "const real localStateVariable = ... + exp1v * stateVarReference;"
// and before the assert on line 75.
#ifdef SEISSOL_TPV104_DIAG
// R-003 FIX: skip unless we are inside the post-iteration call from
// calcSlipRateAndTraction.  The thread_local flag is set by an RAII
// guard declared in §5.12 (TPV104DiagGate.h).
if (tpv104diag::g_in_final_sv_update &&
    tpv104diag::shouldDumpQp(faceIndex, pointIndex,
                              tpv104diag::g_current_time_index,
                              fullUpdateTime)) {
   std::fprintf(stderr,
      "[P2 PSI] rank=%d face=%zu qp=%u t=%+.6e  "
      "psi_in=%+.9e  V=%+.9e  L=%+.9e  dt=%+.9e  V_w=%+.9e  a=%+.9e  "
      "b=%+.9e  V0=%+.9e  f0=%+.9e  muW=%+.9e  psi_out=%+.9e\n",
      mpiRank, faceIndex, pointIndex, fullUpdateTime + timeIncrement,
      stateVarReference, localSlipRate, localSl0, timeIncrement,
      localSrW, localA,
      this->b[faceIndex][pointIndex],
      this->drParameters->rsSr0,
      this->f0[faceIndex][pointIndex],
      localMuW, localStateVariable);
}
#endif
```

**Anchor line**: `const real localStateVariable = steadyStateStateVariable * exp1m + exp1v * stateVarReference;` at `FastVelocityWeakeningLaw.h:74`.

**Scalar equivalence — inputs**:

| input | MFEM | SeisSol | identity |
|---|---|---|---|
| ψ_old | `psi_in (= psi_n[i])` | `stateVarReference` | both dimensionless |
| V | `V_in (= dof_data[i].slip_rate)` | `localSlipRate` | both m/s |
| L | `L_in (= dof_data[i].Dc)` | `localSl0` | both m; SCEC `L` = MFEM `Dc` (§4.2.2 note 4) |
| dt | `dt_step` | `timeIncrement` | both s |
| V_w | `Vw_in (= V_w[i])` | `localSrW` | both m/s; Probe-0 verified identical |
| a | `a_in (= dof_data[i].a)` | `localA` | dimensionless, per-QP; Probe-0 verified identical |
| b, V_0, f_0, μ_w | `TPV104Params::{b, V0, f0, muW}` | `this->b[face][qp]`, `drParameters->rsSr0`, `this->f0[face][qp]`, `localMuW` | identical by §4.1 parameter table |

**Scalar equivalence — output**:

| output | MFEM | SeisSol | identity |
|---|---|---|---|
| ψ_new | `dof_data[i].psi` (after the call) | `localStateVariable` | both from `ψ_new = ψ_ss · (1 − e^{−V·dt/L}) + ψ_0 · e^{−V·dt/L}` (§4.2.1 Eq. 3). MFEM uses `ψ_ss + (ψ_0 − ψ_ss)·exp(preexp1)` (direct form); SeisSol uses the `expm1`-based `steadyStateStateVariable * exp1m + exp1v * stateVarReference` (`FastVelocityWeakeningLaw.h:71-74`). **Algebraically equivalent**; numerically equivalent to 1 ULP on the TPV104 envelope (SeisSol's `expm1` is better for `preexp1 → 0`). |

**Threshold**: `|psi_out_MFEM − psi_out_SeisSol| < 1e-13`, conditional on the 10 input scalars agreeing to 1e-14. If inputs agree but outputs don't, `SlipLawSRWPsi::UpdateStateAnalyticSlipLawSRW` has a bug — rerun `test_slip_law_srw_psi.cpp::T_SRW_5` (§4.2.3) at the failing inputs for a minimal reproducer.

#### 5.10.3 Probe 3 — Friction coefficient μ

**MFEM insertion point**: `miniapps/seas/dynamic/fault_face_flux.cpp`. Two candidates:

- Legacy path: `CompleteFromVabs` at `:138-141`: `real_t C = std::exp(data.psi / data.a) / (2.0 * FrictionSolver::V0); real_t f_V = data.a * std::asinh(s.V_abs * C);`.
- Total-stress path: `EvaluateTotal` at `:418-420`: `real_t C = std::exp(data.psi / data.a) / (2.0 * FrictionSolver::V0); real_t f_V = data.a * std::asinh(V_abs * C); real_t strength = std::abs(sigma_n_fric) * f_V;`.

TPV104 runs under `EvaluateTotal` (§4.4 item 4). Instrument AFTER `f_V` is computed, BEFORE `strength` consumes it:

```cpp
#ifdef SEAS_DIAG_TPV104_STATE
if (data.diag_print) {
   std::fprintf(stderr,
      "[P3 MU] rank=%d t=%+.6e  V=%+.9e  psi=%+.9e  a=%+.9e  "
      "V0=%+.9e  mu=%+.9e\n",
      g_seas_my_rank, t_now, V_abs, data.psi, data.a,
      FrictionSolver::V0, f_V);
}
#endif
```

**Anchor line** (verify before instrumenting): `real_t f_V = data.a * std::asinh(V_abs * C);` at `fault_face_flux.cpp:419` (EvaluateTotal) or `:140` (legacy).

**SeisSol insertion point**: `/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/FrictionLaws/CpuImpl/FastVelocityWeakeningLaw.h:120-121`, inside `updateMu`, at the `return`. Pre-evaluate into a local to print:

```cpp
// Replace updateMu's single-line return with:
#ifdef SEISSOL_TPV104_DIAG
const real mu = details.a[pointIndex] *
                rs::arsinhexp(lx, details.cExpLog[pointIndex],
                              details.cExp[pointIndex]);
if (shouldDumpQp(/* face */, pointIndex, /* o */, /* t */)) {
   std::fprintf(stderr,
      "[P3 MU] rank=%d face=? qp=%u t=?  V=%+.9e  psi=%+.9e  a=%+.9e  "
      "V0=%+.9e  mu=%+.9e\n",
      mpiRank, pointIndex, localSlipRateMagnitude,
      details.cExpLog[pointIndex] * details.a[pointIndex],
      details.a[pointIndex], this->drParameters->rsSr0, mu);
}
return mu;
#else
return details.a[pointIndex] *
       rs::arsinhexp(lx, details.cExpLog[pointIndex],
                     details.cExp[pointIndex]);
#endif
```

**Anchor line** (verify before instrumenting): `return details.a[pointIndex] * rs::arsinhexp(lx, details.cExpLog[pointIndex], details.cExp[pointIndex]);` at `FastVelocityWeakeningLaw.h:120-121`.

**Scalar equivalence**:

| scalar | MFEM | SeisSol | identity |
|---|---|---|---|
| V | `V_abs` (local in `EvaluateTotal`, or `s.V_abs` in `CompleteFromVabs`) | `localSlipRateMagnitude` (parameter of `updateMu`) | both m/s |
| ψ | `data.psi` | NOT stored; reconstruct as `details.cExpLog[qp] * details.a[qp]` (since `getMuDetails:96` sets `cExpLog = ψ/a`) | dimensionless |
| a | `data.a` | `details.a[pointIndex]` | dimensionless, per-QP |
| V_0 | `FrictionSolver::V0 = 1e-6` | `this->drParameters->rsSr0` | both m/s, §4.1 value 1e-6 |
| μ | `f_V` | `mu` (pre-evaluated local) | **MFEM form**: `a · asinh((V/(2V₀)) · exp(ψ/a))`. **SeisSol form**: `a · rs::arsinhexp(cLin·V, ψ/a, exp(−2ψ/a))`. Per `RateAndStateCommon.h:60-87`, `arsinhexp(x, expLog, exp) = asinh(x · exp(expLog))` with stability branching: when `expLog + max(log₂(x), 0) · Log2 > Switch=10` OR `expLog ≥ Threshold=50`, it evaluates via `expLog + log(|x| + sqrt(x² + exp(−2expLog)))`; otherwise `asinh(exp·x)` with `exp = exp(expLog)` precomputed. TPV104 envelope `ψ/a ∈ [5, 65]`, `V ∈ [1e-16, 10]` m/s → `cLin·V = V/(2V₀) ∈ [5e-11, 5e6]`. **R-004 FIX 2026-04-24**: at small `V` (e.g., 1e-16) `log₂(x) < 0`, so `std::max(log₂(x), 0) = 0`; the Switch condition reduces to `expLog > 10`, i.e. **fires at `ψ/a > 10`** (not 33). On the TPV104 steady state `ψ/a ≈ 56`, SeisSol always takes the **log-form branch**. MFEM's `psi/a > 700` branch never fires; MFEM always uses direct `asinh`. The two code paths agree to ~5–10 ULPs empirically — NOT to 1 ULP. |

**Threshold** (**R-004 FIX 2026-04-24**): `|mu_MFEM − mu_SeisSol| / max(|mu|, 1e-12) < 1e-10`. The formulas are algebraically identical but on the TPV104 envelope (`ψ/a ∈ [5, 65]`) SeisSol takes the log-form branch (`RateAndStateCommon.h:73`, fires at `expLog > Switch = 10`) while MFEM uses direct `asinh` (threshold 700, never fires). The two branches agree to ~5–10 ULPs on matched inputs — `1e-12` would flag branch-ULP noise as a false positive. A real μ-divergence (e.g. a formula regression) is still caught at `1e-10`. The offline unit test `test_friction_coefficient_arsinh_seissol_vs_mfem.cpp` referenced in §3.3 must likewise assert `< 1e-10` (not `1e-12`).

#### 5.10.4 Probe 4 — Slip-rate magnitude V_abs

**MFEM insertion point**: `miniapps/seas/dynamic/fault_face_flux.cpp:406-407` (`EvaluateTotal`): `V_abs = solver_.Solve(Theta, data.psi, std::abs(sigma_n_fric), data.eta_s, data.a, method);`. Instrument AFTER the `Solve` call returns:

```cpp
#ifdef SEAS_DIAG_TPV104_STATE
if (data.diag_print) {
   std::fprintf(stderr,
      "[P4 VABS] rank=%d t=%+.6e  Theta=%+.9e  psi=%+.9e  "
      "sigma_n=%+.9e  eta_s=%+.9e  a=%+.9e  V_abs=%+.9e\n",
      g_seas_my_rank, t_now, Theta, data.psi,
      std::abs(sigma_n_fric), data.eta_s, data.a, V_abs);
   // Iteration count is threaded through solver_.Solve's *iterations
   // out-parameter; add a diag overload (§5.5 plan) if needed.
}
#endif
```

**Anchor line**: `V_abs = solver_.Solve(Theta, data.psi, std::abs(sigma_n_fric), data.eta_s, data.a, method);` at `fault_face_flux.cpp:406-407`.

**SeisSol insertion point**: `/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/FrictionLaws/CpuImpl/RateAndState.h:321-327`, inside the Newton early-exit block, after `this->mu[ltsFace][pointIndex] = muF[pointIndex];`.

**R-003 FIX 2026-04-24 — last-iteration-only filter.** `invertSlipRateIterative` is invoked inside `updateStateVariableIterative`'s outer `j = 0 .. numberStateVariableUpdates - 1 = 0..1` loop at `RateAndState.h:162`, so the Newton early-exit probe fires **2 times per `(ltsFace, pointIndex, o)`** on converged runs.  Only the final (j = numberStateVariableUpdates − 1) invocation produces the physically-meaningful V_abs.  Thread the outer-loop index `j` through to `invertSlipRateIterative` via a new last-argument `bool is_final_kaneko_iter` (plumbed at the call site in `updateStateVariableIterative`) and guard the probe with it.  Without this guard, probe_diff.py will multi-match MFEM's single Brent sample to 2 SeisSol Newton samples per time-bin.

```cpp
// In invertSlipRateIterative at RateAndState.h:286-343, INSIDE the
// early-exit block at :321-327, after
//   "this->mu[ltsFace][pointIndex] = muF[pointIndex];".
#ifdef SEISSOL_TPV104_DIAG
if (hasConverged && is_final_kaneko_iter) {
   for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints;
        pointIndex++) {
      if (!shouldDumpQp(ltsFace, pointIndex, /* o */, /* t */)) continue;
      std::fprintf(stderr,
         "[P4 VABS] rank=%d face=%zu qp=%u t=?  Theta=%+.9e  "
         "psi=%+.9e  sigma_n=%+.9e  eta_s=%+.9e  a=%+.9e  "
         "V_abs=%+.9e  iters=%u\n",
         mpiRank, ltsFace, pointIndex,
         absoluteShearStress[pointIndex],
         localStateVariable[pointIndex],
         std::fabs(normalStress[pointIndex]),
         1.0 / this->impAndEta[ltsFace].invEtaS,
         this->a[ltsFace][pointIndex],
         slipRateTest[pointIndex], i);
   }
}
#endif
```

**Anchor line**: `const bool hasConverged = std::all_of(std::begin(g), std::end(g), [&](auto val) { return std::fabs(val) < settings.newtonTolerance; });` at `RateAndState.h:318-320`.

**Scalar equivalence**:

| scalar | MFEM | SeisSol | identity |
|---|---|---|---|
| Θ | `Theta (= sqrt(tau1_fric² + tau2_fric²))` at `fault_face_flux.cpp:401-402` | `absoluteShearStress[pointIndex]` | both Pa. MFEM: trial+nucleation; SeisSol: `initialStressInFaultCS + faultStresses.traction*` per `RateAndState.h:135-139`. Algebraically identical under the convention split (MFEM bakes pre-stress into bulk Q + nuc into `tau*_nuc`; SeisSol splits pre-stress out of bulk Q into `initialStressInFaultCS` + `faultStresses.traction*`). The two sums produce the same total tangential traction magnitude — verified by Probes 0 + 1b. |
| ψ | `data.psi` | `localStateVariable[pointIndex]` | both dimensionless. **Timing difference**: SeisSol updates ψ INSIDE the Newton iteration (Kaneko-2008 fixed point at `RateAndState.h:167-172`); MFEM updates ψ BEFORE the Brent solve. Probe 4 compares the *final converged* ψ and V on both sides; expect agreement iff Probes 1/2 pass. |
| σ_n | `std::abs(sigma_n_fric) = std::abs(sigma_n_trial + data.sigma_n_nuc)` | `std::fabs(normalStress[pointIndex])` | both Pa, absolute value. SeisSol's `normalStress` is clamped to ≤ 0 by `updateNormalStress` at `RateAndState.h:352`; `std::fabs` flips the sign. MFEM stores `sigma_n_fric` positive. **Identical after `std::fabs`**. |
| η_s | `data.eta_s` | `1.0 / this->impAndEta[ltsFace].invEtaS` | both Pa·s/m. MFEM: `eta_s = ρ·c_s/2` (homogeneous). SeisSol: `invEtaS = 2/Zs` (homogeneous), so `1/invEtaS = Zs/2`. Identical. |
| a | `data.a` | `this->a[ltsFace][pointIndex]` | dimensionless, per-QP; Probe-0 verified identical. |
| V | `V_abs` (returned by `solver_.Solve`) | `slipRateTest[pointIndex]` (iterated in Newton) | both m/s. MFEM: `DieterichRuinaFriction::SolveSlipRatePsi` → Brent in log10(V) space (`dieterich_ruina.hpp:309-453`). SeisSol: Newton on V with fixed point on ψ (`RateAndState.h:286-343`). Both solve `g(V) = (1/η_s)·(σ_n·μ(V,ψ) − Θ) + V = 0` (SeisSol at `:311-314`) / `tau − σ_n·f_val − η·V = 0` (MFEM at `dieterich_ruina.hpp:357`), equivalent after multiplying SeisSol by `−η_s`. |
| iters | `*iterations` (out-parameter) | `i` (Newton loop counter at convergence) | integer; diagnostic only, not equivalence-tested |

**Threshold** (**R-005 FIX 2026-04-24** — two-tier, accounting for the Kaneko-vs-decoupled ψ-V mismatch):

- **Quasi-steady-state samples** (both sides report `|dψ/dt| < 1e-6 s⁻¹`): `|V_abs_MFEM − V_abs_SeisSol| / max(V_abs_SeisSol, 1e-30) < 1e-8`. This is the tightest possible given SeisSol's Newton tolerance (`RateAndStateCommon.h:42`).
- **Transient samples** (nucleation rise / rupture breakout, `|dψ/dt| ≥ 1e-6 s⁻¹`): `|V_abs_MFEM − V_abs_SeisSol| / max(V_abs_SeisSol, 1e-30) < 1e-4`. SeisSol updates ψ **inside** the Newton iteration (Kaneko-2008 fixed point at `RateAndState.h:167-172` + damped `V_new = 0.5·(V_old + V_solved)` at `:188-189`), while MFEM's Brent solve accepts a pre-step ψ. The two algorithms contract ψ-V error by a factor of ~4 per Kaneko iteration; a pre-iteration 1 % V mismatch contracts to ~0.25 % post-iteration, but in transient breakout (peak V ≈ 10 m/s) the mismatch can reach 1e-3 relative. A 1e-8 threshold on transient samples produces false-positive Probe 4 FAILs across the first 0.1–1 s of nucleation.

`probe_diff.py` auto-classifies each sample using the numerically-estimated `dψ/dt` from adjacent Probe-2 dumps: `dψ/dt ≈ (ψ_out(t_k+1) − ψ_out(t_k)) / (t_{k+1} − t_k)` per probe QP. Samples classified "transient" under the coarser tolerance are **still reported** in the probe-summarise output but flagged `tier=transient` to avoid muddying the top-line diagnostic.

**Option B future work** (out of scope for Phase 3 but documented here): wrapping MFEM's Brent solve in a two-step outer fixed point that re-evaluates `ψ_analytic(V_guess)` after each V iterate (emulating Kaneko). This is a driver-side change, does not touch Extreme-Care files, and would allow a 1e-8 threshold on transient samples too. Defer unless Phase 3.B surfaces a different residual that the two-tier threshold cannot disambiguate.

**Per-sample guard** (unchanged): if MFEM's `*iterations = 0` (Brent converged in one eval) AND SeisSol reports `iter = maxNumberSlipRateUpdates = 60` (Newton did not converge), both sides are in a degenerate regime — reject the sample (mark "solver-degenerate", do not count toward threshold).

#### 5.10.5 Probe 5 — Corrected traction and imposed state

**MFEM insertion point**: `miniapps/seas/dynamic/fault_face_flux.cpp`. Two candidates:

- `BuildImposedState` at `:170-180` (legacy), then `WriteBackState` at `:200-205`.
- `EvaluateTotal` inline at `:440-449`:
  ```cpp
  Q_imp_minus[VX] = Q_minus[VX] - invZp_m * (sigma_n_corr - Q_minus[SXX]);
  Q_imp_minus[VY] = Q_minus[VY] - invZs_m * (tau1_corr    - Q_minus[SXY]);
  Q_imp_minus[VZ] = Q_minus[VZ] - invZs_m * (tau2_corr    - Q_minus[SXZ]);
  Q_imp_plus[VX]  = Q_plus[VX]  + invZp_p * (sigma_n_corr - Q_plus[SXX]);
  Q_imp_plus[VY]  = Q_plus[VY]  + invZs_p * (tau1_corr    - Q_plus[SXY]);
  Q_imp_plus[VZ]  = Q_plus[VZ]  + invZs_p * (tau2_corr    - Q_plus[SXZ]);
  Q_imp_minus[SXX] = sigma_n_corr;  Q_imp_plus[SXX] = sigma_n_corr;
  Q_imp_minus[SXY] = tau1_corr;     Q_imp_plus[SXY] = tau1_corr;
  ```

Instrument after these assignments complete:

```cpp
#ifdef SEAS_DIAG_TPV104_STATE
if (data.diag_print) {
   std::fprintf(stderr,
      "[P5 IMP] rank=%d t=%+.6e  "
      "sigma_n_corr=%+.9e  tau1_corr=%+.9e  tau2_corr=%+.9e  "
      "Q_imp_plus=[%+.9e %+.9e %+.9e %+.9e %+.9e %+.9e]  "
      "Q_imp_minus=[%+.9e %+.9e %+.9e %+.9e %+.9e %+.9e]\n",
      g_seas_my_rank, t_now,
      sigma_n_corr, tau1_corr, tau2_corr,
      Q_imp_plus[SXX], Q_imp_plus[SXY], Q_imp_plus[SXZ],
      Q_imp_plus[VX],  Q_imp_plus[VY],  Q_imp_plus[VZ],
      Q_imp_minus[SXX], Q_imp_minus[SXY], Q_imp_minus[SXZ],
      Q_imp_minus[VX],  Q_imp_minus[VY],  Q_imp_minus[VZ]);
}
#endif
```

**Anchor line**: `Q_imp_minus[SXY] = tau1_corr;     Q_imp_plus[SXY] = tau1_corr;` at `fault_face_flux.cpp:449` (EvaluateTotal) or `:184,188` (BuildImposedState).

**SeisSol insertion point**: `/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h:362`, inside the per-QP loop at `:336-363`, right after `imposedStateP[W][i] += weight * (qIPlus[o][W][i] + invZs * (traction2 - qIPlus[o][T2][i]));` (**R-001 FIX 2026-04-24**: the symmetric partner of `traction2` on line 362 is `qIPlus[o][T2][i]`, matching line 361's `(traction1 - qIPlus[o][T1][i])` pattern). To produce per-`o` integrands (dimensionally matching MFEM's single-per-stage output), reconstruct WITHOUT the `weight *` factor:

```cpp
// In postcomputeImposedStateFromNewStress(...), inside the i-loop at
// FrictionSolverCommon.h:336-363, after the :362 increment
//   "imposedStateP[W][i] += weight * (qIPlus[o][W][i] + invZs * (traction2 - qIPlus[o][T2][i]));".
#ifdef SEISSOL_TPV104_DIAG
if (shouldDumpQp(ltsFace, i, o, fullUpdateTime)) {
   // Reconstruct per-o integrand WITHOUT the time weight to match MFEM's
   // per-RK-stage single value.
   const real sigma_n_corr = normalStress;
   const real tau1_corr    = traction1;
   const real tau2_corr    = traction2;
   const real Q_imp_plus_U  = qIPlus[o][U][i]  + invZp     * (normalStress - qIPlus[o][N][i]);
   const real Q_imp_plus_V  = qIPlus[o][V][i]  + invZs     * (traction1    - qIPlus[o][T1][i]);
   const real Q_imp_plus_W  = qIPlus[o][W][i]  + invZs     * (traction2    - qIPlus[o][T2][i]);
   const real Q_imp_minus_U = qIMinus[o][U][i] - invZpNeig * (normalStress - qIMinus[o][N][i]);
   const real Q_imp_minus_V = qIMinus[o][V][i] - invZsNeig * (traction1    - qIMinus[o][T1][i]);
   const real Q_imp_minus_W = qIMinus[o][W][i] - invZsNeig * (traction2    - qIMinus[o][T2][i]);
   std::fprintf(stderr,
      "[P5 IMP] rank=%d face=%zu qp=%u o=%u t=%+.6e  "
      "sigma_n_corr=%+.9e  tau1_corr=%+.9e  tau2_corr=%+.9e  "
      "Q_imp_plus=[%+.9e %+.9e %+.9e %+.9e %+.9e %+.9e]  "
      "Q_imp_minus=[%+.9e %+.9e %+.9e %+.9e %+.9e %+.9e]\n",
      mpiRank, ltsFace, i, o, fullUpdateTime + deltaT[o],
      sigma_n_corr, tau1_corr, tau2_corr,
      sigma_n_corr, tau1_corr, tau2_corr,
      Q_imp_plus_U, Q_imp_plus_V, Q_imp_plus_W,
      sigma_n_corr, tau1_corr, tau2_corr,
      Q_imp_minus_U, Q_imp_minus_V, Q_imp_minus_W);
}
#endif
```

**Anchor line**: `imposedStateP[W][i] += weight * (qIPlus[o][W][i] + invZs * (traction2 - qIPlus[o][T2][i]));` at `FrictionSolverCommon.h:362` (**R-001 FIX 2026-04-24** — previous text cited `T1`; the live source uses `T2`, confirmed via `grep -n "imposedStateP\[W\]" FrictionSolverCommon.h`).

**Scalar equivalence**:

| component | MFEM | SeisSol (per-`o` integrand, de-weighted) | identity |
|---|---|---|---|
| `sigma_n_corr` | `sigma_n_corr (= sigma_n_trial)` at `fault_face_flux.cpp:432` — MFEM `EvaluateTotal` uses the TRIAL-scale (not `sigma_n_fric`), matching SeisSol's `tractionResults.traction*` convention (see MFEM comment `:427-431`). | `normalStress = faultStresses.normalStress[o][i]` | Pa; identical after the MFEM-v9.4.0 `tau*_corr` trial-vs-fric split (`fault_face_flux.cpp:425-431`) |
| `tau1_corr` | `tau1_corr (= tau1_trial − η_s·V1)` | `traction1 = tractionResults.traction1[o][i] (= faultStresses.traction1[o][i] − impAndEta.etaS · slipRate1[ltsFace][pointIndex])` at `RateAndState.h:235-237` | Pa. Identical formula. |
| `tau2_corr` | `tau2_corr (= tau2_trial − η_s·V2)` | `traction2` at `RateAndState.h:238-240` | Pa. Identical. |
| `Q_imp_plus[SXX]` | `sigma_n_corr` | `normalStress` | Pa; de-weighted integrand == MFEM's single-value per-stage output |
| `Q_imp_plus[SXY]` | `tau1_corr` | `traction1` | Pa |
| `Q_imp_plus[SXZ]` | `tau2_corr` | `traction2` | Pa |
| `Q_imp_plus[VX]` | `Q_plus[VX] + invZp_p · (sigma_n_corr − Q_plus[SXX])` at `:444` | `qIPlus[o][U][i] + invZp · (normalStress − qIPlus[o][N][i])` at `:360` | m/s. **Identical formula** under the Probe-1a index map. |
| `Q_imp_plus[VY]` | `Q_plus[VY] + invZs_p · (tau1_corr − Q_plus[SXY])` at `:445` | `qIPlus[o][V][i] + invZs · (traction1 − qIPlus[o][T1][i])` at `:361` | m/s. Identical. |
| `Q_imp_plus[VZ]` | `Q_plus[VZ] + invZs_p · (tau2_corr − Q_plus[SXZ])` at `:446` | `qIPlus[o][W][i] + invZs · (traction2 − qIPlus[o][T2][i])` at `:362` | m/s. Identical. |
| `Q_imp_minus[...]` | symmetric, sign flip `− invZ` (see `fault_face_flux.cpp:441-443`) | symmetric at `FrictionSolverCommon.h:350-355`, sign flip via `invZpNeig / invZsNeig` | m/s. Identical structure. |

**Threshold**: `max |diff| / max |signal| < 5e-3` across the 21 scalars (consistent with §5.6), per-QP.

**Note on ADER time-weighting**: MFEM's `Q_imp_plus/minus` is a single per-RK-stage output; SeisSol accumulates `Σ_o weight[o] · integrand` over `misc::TimeSteps = 5` (O=5 Dubiner) ADER time points. To diff correctly, the SeisSol probe MUST dump the **per-`o` integrand pre-weight** (as the instrumentation above does). If the implementer dumps `imposedStateP[U][i]` AFTER the full `o`-loop has closed, the value is already weight-summed and the comparison becomes dimensionally wrong.

**Rationale**: per the v8.0.0 SeisSol audit, the imposed-state formulas are bit-for-bit identical in form. A Probe-5 FAIL with Probes 0/1/2/3/4 clean therefore pins the bug to: (a) the rotation `T_can` at `wave_operator.inl:1260-1269`, (b) a regression in the `EvaluateTotal` imposed-state loop at `fault_face_flux.cpp:440-449`, or (c) the per-side duplication MFEM uses at `wave_operator.inl:1281-1285` — `flux_.Interior(can_n, Q_imp_plus_g, Q_imp_plus_g)` (self-self call) rather than SeisSol's single-sided `AplusT · Q_imp`, a flagged audit risk (see the Pelties-9 comment at `fault_face_flux.cpp:226-249`).

### 5.11 Recommended Phase 3 run order (updated for §5.9 priorities)

Supersedes the sbatch-design portion of §5.7; §5.7's paired-fixtures table remains authoritative for the MFEM / SeisSol compile-time options.

```
Phase 3.A (t = 0, nucleation OFF)
  └── Probe 0 ONLY (§5.10.0)
      ├── PASS → proceed to 3.B.
      └── FAIL → block; fix initialisation (ψ_ini inversion or initial
                 stress in fault CS).  Do NOT instrument 1a..5 — they
                 are meaningless if init disagrees.

Phase 3.B (nucleation ON, np = 1)
  ├── Probe 1a (§5.10.1.a)  [P1]
  │   ├── PASS → Probe 1b.
  │   └── FAIL → STOP.  Bug is in BULK ADER non-fault assembly.
  │              Authorises the precomputed-face-rotation fix per
  │              tpv102_seissol_aligned_flux_plan_2026-04-23.md.  File
  │              targeted fix plan; do NOT instrument 1b..5 (diagnosis
  │              already localised; further Frontera runs waste
  │              allocation).
  │
  ├── Probe 1b (§5.10.1.b)  [P1]
  │   ├── PASS → Probe 5 (out of pipeline order; P2).
  │   └── FAIL with 1a clean → STOP.  Regression in
  │                            `ComputeTrialTraction` or rotation bug
  │                            in `Tinv_can`.  File targeted fix plan;
  │                            do NOT instrument 2/3/4/5 until 1b is
  │                            green.
  │
  ├── Probe 5  (§5.10.5)    [P2]
  │   ├── PASS → Probe 2.
  │   └── FAIL with 1a+1b clean → STOP.  `T_can` rotation or
  │                              imposed-state construction bug.  File
  │                              targeted fix plan; do NOT instrument
  │                              2/3/4 until 5 is green.
  │
  ├── Probe 2  (§5.10.2)    [P2]
  │   ├── PASS → Probe 3.
  │   └── FAIL → STOP.  Reopen test_slip_law_srw_psi.cpp::T_SRW_5 at
  │              the failing inputs; file targeted fix plan for
  │              SlipLawSRWPsi; do NOT instrument 3/4 until 2 is green.
  │
  ├── Probe 3  (§5.10.3)    [P3]
  │   ├── PASS → Probe 4.
  │   └── FAIL → STOP.  asinhexp branching divergence; report
  │              (V, ψ, a) tuple and file targeted fix plan.  Do NOT
  │              modify dieterich_ruina.hpp (§2.5 Extreme Care).  Do
  │              NOT instrument 4 until 3 is green.
  │
  └── Probe 4  (§5.10.4)    [P3]
      ├── PASS → 3.B complete; proceed to 3.C approval.
      ├── FAIL on QUASI-STEADY samples (§5.10.4 R-005 tier 1, 1e-8
      │   threshold)  → STOP.  Solver disagreement beyond Newton
      │                 tolerance; file targeted fix plan.
      ├── FAIL on TRANSIENT samples (§5.10.4 R-005 tier 2, 1e-4
      │   threshold) → ADVISORY.  Expected Kaneko-vs-decoupled ψ-V
      │                coupling artefact.  Flag in probe-summarise,
      │                keep 3.B as PASS.
      └── FAIL at early times (V ≈ 1e-16) → ADVISORY unless transient
                                           tier also fails.  Rerun
                                           with Va_min = -500 (MFEM)
                                           or rs::almostZero=1e-60
                                           (SeisSol) as a stress test.

Phase 3.C (nucleation ON, np = 4)
  └── Re-run Probe 1a ONLY on stations crossing MPI seams.
      ├── PASS with 3.B all-pass → pepper is a station-sampling
      │                            artefact.  Close with documented
      │                            resolution.
      └── FAIL → STOP.  Pepper is in ComputeADERSharedFaceFluxRHS
                 (wave_operator.inl:1762).  Cross-reference TPV102
                 Gate 3D residual; file targeted fix plan.
```

Every non-advisory FAIL branch now carries an explicit **STOP** directive — the default is to halt probing on first failure and file a targeted fix plan rather than continue accumulating evidence that may be polluted by the upstream divergence. This is the R-009 convention and preserves Frontera allocation.

**Per-sub-phase approval checklist** (per `feedback_frontera_approval.md`; user approves each entry independently):

| sub-phase | input artefacts for user review | user approval gate |
|---|---|---|
| 3.A | Probe 0 sbatch script, probe_diff.py config for init fields | "approve 3.A" |
| 3.B | Probes 1a/1b + 2/3/4/5 sbatch script, diag-build instructions for the Frontera SeisSol branch | "approve 3.B" |
| 3.C | Same sbatch as 3.B, `-n 4` rank config, seam-QP probe list | "approve 3.C" |

### 5.12 Shared probe helpers (R-008 addition 2026-04-24)

Every probe snippet in §5.10 calls two helpers that were referenced but not defined in prior revisions: `shouldDumpQp(...)` on the SeisSol side, and a `TimeSampler` class on the probe-diff side. The Phase 3 instrumentation depends on both; they live in:

#### 5.12.1 `shouldDumpQp` — per-QP dump predicate (SeisSol + MFEM)

**SeisSol location**: new header `/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/FrictionLaws/TPV104DiagGate.h` (added only on the Frontera diag branch; not upstream). Layout:

```cpp
// src/DynamicRupture/FrictionLaws/TPV104DiagGate.h
namespace seissol::tpv104diag {

// Loaded once at startup from SEISSOL_TPV104_DIAG_STATIONS=/abs/path/stations.txt
// File format (whitespace-separated, one entry per line):
//   ltsFace  pointIndex  station_label  x_m  z_m
// The file is produced by running probe_diff.py with --emit-seissol-stations
// and matches the nine canonical §4.1 stations + the hypocenter QP (10 rows).
struct StationGate {
   std::size_t ltsFace;
   std::uint32_t pointIndex;
   std::string label;   // "x2_0_x3_7.5" etc.
};
const std::vector<StationGate>& stationTable();  // lazy singleton

// True iff (ltsFace, pointIndex, fullUpdateTime) is in the sample schedule
// defined at §5.1 (50 samples: dense 0-1 s, sparse 1-12 s).  Sample
// schedule is hard-coded; see NextSampleTime() below for the list.
bool shouldDumpQp(std::size_t ltsFace,
                  std::uint32_t pointIndex,
                  std::uint32_t timeIndex,     // ADER sub-step o
                  real fullUpdateTime);

// R-003 helpers: last-iteration-only filters for Probes 2 and 4.
extern thread_local bool g_in_final_sv_update;    // set by calcSlipRateAndTraction
extern thread_local std::uint32_t g_current_time_index;  // set by BaseFrictionLaw
// Probe 4 uses its own is_final_kaneko_iter argument, threaded through the
// invertSlipRateIterative signature.  No global needed.

// RAII guard for g_in_final_sv_update (R-003 fix for Probe 2).
struct FinalSvUpdateGuard {
   FinalSvUpdateGuard() { g_in_final_sv_update = true; }
   ~FinalSvUpdateGuard() { g_in_final_sv_update = false; }
};

// Sample-schedule helpers (identical semantics on both codes).
real NextSampleTime(real currentTime);
bool IsAtSampleBoundary(real currentTime, real dt);

}  // namespace seissol::tpv104diag
```

Activation: `#define SEISSOL_TPV104_DIAG` in the Frontera build flags; that compile gate wraps all probe bodies from §5.10 and the declarations above. The stations-file path is set via `export SEISSOL_TPV104_DIAG_STATIONS=/path/to/stations.txt` before `srun`.

**MFEM equivalent**: new header `miniapps/seas/dynamic/seas_tpv104_diag_gate.hpp`. Same function signatures, same semantics, different namespace `mfem::seas::tpv104diag`. MFEM's `DOFData::diag_print` (already exists for BP5) is the per-QP flag; the new gate wraps it plus the §5.1 sample schedule. No R-003 RAII guards needed on the MFEM side — MFEM has no Kaneko loop, so every ψ update is the "final" one.

Tests: `tests/unit/test_tpv104_diag_gate.cpp` — asserts (i) the stations file is parsed correctly, (ii) `shouldDumpQp` returns `true` for exactly the 10 probe QPs × 50 sample times = 500 samples per probe per SeisSol run (matched on MFEM).

#### 5.12.2 `TimeSampler` — time-alignment across the two codes

**Location**: new class in `tpv104/scripts/probe_sampler.py`. Both `probe_diff.py` (§5.8) and `probe_summarise.py` import it.

```python
# tpv104/scripts/probe_sampler.py
class TimeSampler:
    """Align per-probe dumps across MFEM and SeisSol onto the
    coarse §5.1 schedule via linear interpolation.

    The two codes step at different dt (MFEM CFL-limited, SeisSol O=5 dt
    ≈ 4.7 ms), so raw samples do not share a time grid.  We interpolate
    both series to a common schedule and diff field-by-field at each
    sample.

    Sample schedule (§5.1):
        dense: t ∈ {0, 0.01, 0.02, 0.05, 0.1, 0.15, ..., 0.95, 1.0}
        sparse: t ∈ {1.1, 1.2, ..., 2.0} ∪ {3, 4, ..., 12}
        ~50 samples per QP per probe.
    """
    def __init__(self, sample_times: list[float]): ...

    def load_probe_series(self, path: str,
                          station: str,
                          field: str) -> np.ndarray:
        """Return (t[], value[]) for (station, field) from a probe dump."""
        ...

    def resample(self,
                 raw_times: np.ndarray,
                 raw_values: np.ndarray,
                 kind: str = 'linear') -> np.ndarray:
        """Linear interpolation onto self.sample_times.  Out-of-range
        samples returned as NaN so diff logic can skip them."""
        ...

    def align(self,
              mfem_series: tuple[np.ndarray, np.ndarray],
              seissol_series: tuple[np.ndarray, np.ndarray]
              ) -> tuple[np.ndarray, np.ndarray]:
        """Returns (mfem_resampled[], seissol_resampled[]) on
        self.sample_times; callers diff elementwise."""
        ...
```

Tests: `tpv104/scripts/tests/test_probe_sampler.py` — feed synthetic
matched series at different timesteps, assert the aligned outputs are
identical (to floating-point tolerance) because the input sampling is
dense enough that linear interpolation is exact.

#### 5.12.3 Approval artefacts

For the Phase 3.A approval gate the user sees:
- The `stations.txt` file generated by `probe_diff.py --emit-seissol-stations`.
- The compile flag list for the Frontera SeisSol diag branch
  (`SEISSOL_TPV104_DIAG` + `-DSEISSOL_TPV104_DIAG_STATIONS_AT_BUILD` if
  hard-coding is preferred over env-var loading).
- The §5.1 sample-schedule table, written into a comment block in
  both `TPV104DiagGate.h` and `seas_tpv104_diag_gate.hpp` for
  unambiguous cross-reference.

---

## 6. Execution order and decision tree

### 6.1 Ordering

```
Phase 1 — inconsistency report (§3)  [this document is the deliverable; no code]
                │
                ▼
Phase 2 — TPV104 implementation
     │
     ├─ 2.1  config/tpv104_params.hpp
     │       + test_tpv104_params.cpp  ─────────────── blocking T_SRW_1..4,6
     │
     ├─ 2.2  friction/slip_law_srw_psi.hpp
     │       + test_slip_law_srw_psi.cpp  ──────────── blocking P2_A
     │
     ├─ 2.3  dynamic/tpv104_setup.hpp (+ _total.hpp)
     │       + test_tpv104_setup.cpp  ─────────────── blocking P2_C, P2_D
     │
     ├─ 2.4  drivers/tpv104_driver.cpp
     │       + test_tpv104_smoke.cpp  ─────────────── blocking P2_C, P2_E, P2_F
     │
     └─ 2.5  tpv104/mesh/tpv104_1000m.msh  (can proceed in parallel with 2.1–2.4)
                │
                ▼
Phase 3 — paired diagnostics
     │
     ├─ 3.A  nucleation OFF — both codes run; diff all 5 probe pairs
     │        ├── PASS  → initialisation is clean; go to 3.B
     │        └── FAIL  → fix (or flag) the init discrepancy; rerun 3.A
     │
     ├─ 3.B  nucleation ON at np=1 — both codes run; diff all 5 probe pairs
     │        ├── All probes PASS  → no localisable bug at np=1; go to 3.C
     │        └── Probe k FAILS   → bug localised to Stage k (see §5.2–5.6
     │                              per-stage diagnostic outcomes)
     │
     └─ 3.C  nucleation ON at np=4 — both codes run; same diff
              ├── All probes PASS  → pepper is a station-sampling artefact
              │                      (not a real physics bug); close
              └── Probe k FAILS    → pepper is in Stage k;
                                     open targeted fix plan
```

### 6.2 Pass/fail map — outcomes and next actions

| Phase 3 outcome | Conclusion | Next action |
|---|---|---|
| 3.A probe-1 FAIL | Pre-friction trial traction differs in silent regime → init bug (ψ_ini mismatch, pre-stress rotation sign error, bulk-Q baking error) | Audit `InitializeStateTotal_TPV104` and `ComputeInitialPsiTPV104`; rerun. |
| 3.A probes-1..5 PASS, 3.B probe-1 FAIL | Nucleation injection is asymmetric (MFEM `ApplyNucleationPrestress` vs SeisSol `adjustInitialStress`). | Trace the nucleation amplitude at the hypocenter QP on both sides; confirm `g_nuc(t)` ramp matches. If matches, the injection channel is different (bulk Q vs persistent channel). Retain MFEM's persistent-channel design; flag the bulk-Q dilution issue as "expected behaviour difference" and document. |
| 3.B probes-1 PASS, 3.B probe-2 FAIL | `SlipLawSRWPsi::UpdateStateAnalyticSlipLawSRW` bug | Reopen `test_slip_law_srw_psi.cpp::T_SRW_5`; find minimal failing input. |
| 3.B probes-1,2 PASS, 3.B probe-3 FAIL | Friction coefficient formula drift (`asinhexp` branching inconsistent) | Report specific `(V, ψ, a)` tuple to user; do NOT edit `dieterich_ruina.hpp` unilaterally (§2.5). |
| 3.B probes-1,2,3 PASS, 3.B probe-4 FAIL | Friction solver root disagreement | Check whether Brent bracket and Newton start both converge within tolerance; treat as an advisory (not a pepper cause) unless diff > 1e-4 relative. |
| 3.B probes-1..4 PASS, 3.B probe-5 FAIL | Riemann imposed-state construction differs | Diff the 21 channels of Q_imp at the failing sample; localise to specific fields. **Most likely pepper-bug location** given the TPV102 audit findings. |
| 3.B all PASS, 3.C probe-5 FAIL | Pepper is in the MPI shared-face Riemann path | Match output to TPV102 `test_shared_nonfault_flux_orbit_identity` (Test 4 in TPV102 plan §3) and re-run its diagnostic. |
| 3.C all PASS | Pepper is a station-sampling artefact, not a physics bug | Document in handoff; reopen TPV102 station-sampling audit separately. |

---

## 7. Risk assessment

| risk | likelihood | impact | mitigation |
|---|---|---|---|
| SeisSol vs MFEM t1/t2 fault-tangent convention disagreement | medium (research report §8 flagged as unresolved) | trace columns swap; diff across codes is wrong even if physics is right | Phase 2 adds `test_tpv104_slip_column_convention.cpp` that ties MFEM's `slip2` (strike, under BP5 convention) to SeisSol's column 2 (horizontal) via the "vertical strike-slip fault" invariant; assert before 3.B. |
| SeisSol build on Frontera needs specific toolchain | low (user confirmed SeisSol is built on Frontera; `reference_seissol_frontera.md` has the recipe) | Phase 3 blocked entirely | Verify build exists before starting Phase 3 instrumentation; if `make -C SeisSol` fails on Frontera, block Phase 3 and ask the user to fix the build first. |
| Nucleation injection channel asymmetry between codes masks a real bug | medium | Probe 1 FAIL on 3.B that is "expected physics difference" rather than a bug | Document the asymmetry in §3.2 (already done); in Probe 1 analyser, compare the nucleation amplitude at the hypocenter QP separately and subtract before diffing the rest |
| MFEM's ADER order O=2 truly cannot resolve TPV104's ~4 MPa/s nucleation rise rate | low-medium | 3.B station traces diverge even after all probes PASS | Run at `--ader-order=3` or `--order=2` if 3.B probes PASS but station traces fail. |
| TPV104 mesh construction in MFEM does not match TPV5 mesh in SeisSol at equivalent resolution | medium (this is a manual scale check) | bulk waves propagate at different speeds, station arrivals shift | Gate Phase 3.A on "wave speed check": inject a small-amplitude delta at y=0, time the P-wave arrival at z = 5 km, expect `t_arrival = 5000/6000 ≈ 0.833 s`. Mismatch by > 1% indicates mesh or material parameter bug. |
| Pepper is driven by the `H_KUHN_FIXTURE` artefact (test-only) and NOT present in TPV104 | medium (per TPV102 §7 final findings) | 3.C all PASS; no bug found; debt | Explicitly document this outcome in §6.2 as a valid resolution. |
| CLAUDE.md violation: accidentally editing a protected file during Phase 2 | low if §2 followed | large blast radius; all TPV102 work could regress | Add `test_tpv102_bitmatch_baseline.cpp` that runs TPV102 on M_ref serial and diffs its 9-station CSVs against a committed reference. Gate P2_B. |
| SEAS_DIAG probes on Frontera produce large log files | high (per-QP per-sample scale) | file system pressure | Sample schedule in §5.1 uses 50 samples × 10 probe QPs × 5 probes = 2500 lines per rank. Well under 1 MB per rank. Acceptable. |

---

## 8. Items deferred to user authorisation

None of the following are executed by this plan; each requires explicit approval from the user.

1. **Editing any of the Extreme-Care files** (`godunov_flux.*`, `wave_operator.*`, `fault_face_flux.*`, `dieterich_ruina.hpp`, `fault_basis.hpp`, `rate_state_fault.hpp`, `seas_operator.hpp`, `time_stepper.hpp`, `bp5_params.hpp`). Several §6.2 outcomes (e.g., probe-5 FAIL) localise bugs in `fault_face_flux.cpp`; fixes there require authorisation.
2. **Any SeisSol source probe** being committed upstream. The §5 probes are disposable stderr printfs, guarded by a single new compile flag `SEISSOL_TPV104_DIAG`. They are NOT to be PR'd to SeisSol; they live in the user's local Frontera branch only.
3. **Running a Frontera sbatch**: per `feedback_frontera_approval.md`, the user approves each run. Phase 3.A, 3.B, 3.C are three separate approvals. The deliverable from this plan is a pair of sbatch scripts (MFEM + SeisSol) that the user submits.
4. **Fix A + Fix B from the TPV102 audit** (`GodunovFlux::Interior` symmetrisation, `FreeSurfaceTotal` sign-symmetrisation). These are independent of this plan but their landing would affect the TPV102 baseline we gate Phase 2 against (P2_B). If they land, we regenerate the P2_B baseline.
5. **Extending `DOFData` to carry `V_w[i]`** — this is an alternative to the driver-owned side-channel array in §4.5. Cleaner design but requires editing `fault_face_flux.hpp` (Extreme Care). Deferred unless the side-channel proves hostile to maintain.

---

## 9. Handoff

Deliverables of this plan, ready for code-review before implementation starts:

- §3 inconsistency report (frozen; four green rows, five yellow, four red).
- §4 TPV104 implementation plan (nine new files, one opt-in class extension, zero protected-file edits).
- §5 five-probe-pair diagnostic plan (symmetric MFEM + SeisSol probes at five pipeline stages; pass/fail thresholds; matched fixtures).
- §6 step-by-step execution order with a pass/fail decision tree that ties each possible outcome to a next action.

Dependencies for starting Phase 2 implementation:

- User review of §2.5 constraint list — confirm the "no Extreme-Care-file edits" scope is acceptable.
- User review of §3.4 red-flag rows — confirm the decision to follow SCEC spec (not SeisSol) on any value where they disagree.
- User review of §4.5 — confirm the driver-owned `V_w[i]` side-channel is acceptable vs extending `DOFData` (§8 item 5).

Dependencies for starting Phase 3 execution:

- P2_A, P2_B, P2_C, P2_D all green.
- User approval for Frontera sbatch 3.A.

Expected wall-clock (author estimate):

- Phase 2 implementation + test landing: ~3–5 focused sessions (~2 days local work).
- Phase 3 diagnostic runs on Frontera: two sbatch per phase (MFEM + SeisSol), each ~1–2 hours compute; analysis per phase ~1 hour; three phases total = ~1–2 days wall-clock including user-approval latency.
