# TPV102 Debug v9.3.0 — SeisSol-parity patch set: free-surface Godunov projection, total-stress migration, and ADER integration point

> **Author:** planning agent, 2026-04-21.
> **Predecessors:**
> - `tpv102_codebase_walkthrough.md` — §"MFEM vs SeisSol implementation
>   inconsistencies" (I-01 … I-07).
> - `tpv102_ader_time_integration_plan.md` — the 8-phase ADER plan for
>   divergence I-05.
> - `tpv102_debug_v9.2.0_debug_plan.md` §19 — REVIEW findings R-V92-E01..E08.
> - REVIEW.md rev 3 — concrete root-cause hypothesis R-007 (bulk σ_yy pump)
>   and R-008 (RK4-averaged output).
>
> **Scope.** Patch the three remaining MFEM ↔ SeisSol divergences flagged
> in the walkthrough's I-section:
>
> - **I-04 (this plan):** Free-surface BC — add a Godunov-projection
>   variant alongside the existing γ-mirror.  Opt-in via
>   `--free-surface-bc={gamma|godunov}`.  Default keeps γ-mirror.
> - **I-05 (this plan):** Integration points only.  The actual ADER
>   implementation is `tpv102_ader_time_integration_plan.md`.  This plan
>   documents where ADER wiring plugs into v9.3.0's total-stress mode
>   (I-06) — the two landings must be coordinated.
> - **I-06 (this plan):** State representation — migrate TPV102 from
>   fluctuation-Q to total-stress Q.  After Phase 4, total stress is
>   the only supported TPV102 representation; there is no user-facing
>   `--state-representation` fallback to fluctuation mode.
>
> Each of I-04, I-06 lands independently of the other and independently
> of ADER.  I-05 requires I-06 as a prerequisite (SeisSol's ADER uses
> total stress internally — porting ADER onto fluctuation Q requires an
> extra subtraction step inside the CK recursion that clutters the
> implementation).  Order: **I-04 → I-06 → I-05 (ADER)**.
>
> Motivation (one-liner per divergence):
> - I-04 addresses plan v9.2.0 §19 H-V92-G1 (free-surface corner σ_yy pump).
> - I-06 eliminates R-007's "fluctuation amplification" mechanism.
> - I-05 eliminates H-V92-K by construction.
>
> **NOT in scope** for v9.3.0: GPU kernels, heterogeneous material
> extension to the per-side Pelties-9 flux, BP5-path changes, LTS.

## Overview

Three patches land in v9.3.0.  I-04 remains opt-in via
`--free-surface-bc={gamma|godunov}` so the old γ-mirror path stays
available for controlled A/B checks.  I-06 is not opt-in: it replaces
TPV102's fluctuation-Q path with total-stress Q.  After v9.3.0 + the
ADER plan land, the combination `--time-integrator=ader
--free-surface-bc=godunov` gives a TPV102 execution path that matches
SeisSol's at the equation-block level, with total stress as the only
supported TPV102 state representation.

## Constraints

### Interface constraints (frozen)
- **`WaveOperator::Mult(Q, dQdt)`** signature — unchanged.  No
  user-facing state-representation branch is added; TPV102 simply uses
  total-stress Q after Phase 4.
- **`FaultFaceFlux::Evaluate(data, Q_plus, Q_minus, Q_imp_plus, Q_imp_minus, method)`** —
  unchanged.  A second entry point `EvaluateTotal` (Phase 3) is added
  alongside; TPV102 dispatch uses `EvaluateTotal` exclusively from
  Phase 4 onward.  The old path may remain only as a migration
  reference in unit tests until follow-up cleanup.
- **`DOFData` struct layout** — append-only.  New fields
  (`Q_total_init_*`) at end of struct; `alignof(DOFData)` must not change.
- **CLI defaults** — `seas_tpv102_driver` has no
  `--state-representation` flag in the final v9.3.0 design.  TPV102
  runs total-stress by default/only; free-surface BC still defaults to
  γ-mirror unless `--free-surface-bc=godunov` is requested.
- **BP5 drivers** — untouched.  The total-stress migration is TPV102-only.
- **Plan v9.2.0 §19 `_k4` VTU diagnostic** — still populated under all
  v9.3.0 TPV102 runs.  Under ADER the `_k4 − avg` delta is zero by
  construction (from the ADER plan §7).  Under RK4 + total the delta
  is still informative.

### Dependency constraints
- No new external libraries.  Eigen headers already ship with SeisSol;
  MFEM has its own `DenseMatrix` and `DenseMatrixInverse`, which is
  what Phase 1 uses.
- No change to MFEM core.

### Convention constraints
- CLI flag pattern: `--free-surface-bc {gamma|godunov}`.  Lowercase,
  hyphenated.  There is no `--state-representation` runtime flag in the
  final v9.3.0 design.
- All new unit tests follow the `tests/unit/test_*` TEST_ASSERT /
  TEST_NEAR pattern.

### Numerical constraints
- **I-04:** on an axis-aligned free surface, γ-mirror and
  Godunov-projection must agree to 10 ULP per component per face.
  This is the equivalence-at-flat-surface regression gate.
- **I-06:** on a locked-fault fixture (V=0), the migrated TPV102
  total-stress path must produce:
  - `data.V1, data.V2, data.slip_rate` all 0.
  - `data.sigma_n_corr, tau1_corr, tau2_corr` equal to the physical
    pre-stress.
  - Bulk Q carrying the constant pre-stress field in global
    coordinates (`Q[SYY] = +sigma_n0`, `Q[SXY] = -tau_ini` for TPV102).
  - A dedicated migration test may compare this path against the frozen
    fluctuation formulation on a small fixture, but fluctuation mode is
    not user-exposed after v9.3.0.
- **Conservation:** both Godunov-projection and the migrated
  total-stress path
  must preserve discrete mass conservation
  (`∫_V Q_new = ∫_V Q − Δt ∫_∂V F_h · n`) to machine precision on a
  periodic domain.

---

## Phase 1: Free-surface Godunov-projection BC (I-04 part A)

### Goal
After this phase, `GodunovFlux::FreeSurfaceGodunov(nor, Q_self, F_h)`
exists and is unit-tested.  Its output equals `GodunovFlux::FreeSurface`
(γ-mirror) to 10 ULP on a flat axis-aligned surface and to 50 ULP on a
10°-tilted surface.  No integration into the wave operator yet — Phase 2
wires it.

### Files to Modify
- `dynamic/godunov_flux.hpp` — add `FreeSurfaceGodunov` declaration + a
  cached `free_surface_godunov_state_` 9×9 matrix member.
- `dynamic/godunov_flux.cpp` — add definition (~80 LOC) using the
  eigenvector matrix `R` already built in the ctor; precompute the
  Godunov-state projector at ctor time.

### Files to Create
- `tests/unit/test_free_surface_godunov.cpp` — compares γ-mirror and
  Godunov-projection on axis-aligned + tilted + fault-surface-corner
  fixtures.

### Detailed Requirements

**1. Math.** SeisSol's `getTransposedFreeSurfaceGodunovState`
(`src/Model/Common.h:251-296`) builds a 9×9 "Godunov state" matrix
`Q_god` such that the free-surface BC `σ·n = 0` is enforced as a
characteristic projection onto the subspace `Null(P_traction)` where
`P_traction` selects the three traction components of the rotated
state.

For isotropic elasticity, `R_{3\times3}` = eigenvectors of `A_x`
restricted to columns 0, 1, 2 (the right-going P + 2 S modes).  Split
into the 3 traction rows `matR11` and 3 velocity rows `matR21`:
$$
\mathbf{R}_{11} = \mathbf{R}\big[\{0, 3, 5\},\ \{0, 1, 2\}\big],\qquad
\mathbf{R}_{21} = \mathbf{R}\big[\{6, 7, 8\},\ \{0, 1, 2\}\big].
$$
The 3×3 "compliance" block is
$$
\mathbf{S} = -\mathbf{R}_{21}\,\mathbf{R}_{11}^{-1}.
$$
The 9-component Godunov state in the rotated frame is then
$$
Q_{\text{god}}[\text{trac}] = 0,\qquad
Q_{\text{god}}[\text{vel}] = Q_{\text{self}}[\text{vel}] + \mathbf{S}\,(\mathbf{0} - Q_{\text{self}}[\text{trac}]).
$$
Equivalently, SeisSol's $Q_{\text{god}}$ is the state that satisfies
$\sigma\!\cdot\!\hat{n} = 0$ with velocity incremented by
$\mathbf{Z}^{-1}\cdot\sigma_{\text{self}}$ where the impedance vector
$\mathbf{Z} = \mathrm{diag}(Z_p, Z_s, Z_s)$ arises naturally from
inverting the traction block of $R$.

For $R$ built in `godunov_flux.cpp:57-83` the explicit result is:
$$
\mathbf{S} = \mathrm{diag}\Big(-\tfrac{1}{Z_p},\ -\tfrac{1}{Z_s},\ -\tfrac{1}{Z_s}\Big).
$$
i.e., $Q_{\text{god}}[v_n]  = Q_{\text{self}}[v_n]  + \tfrac{1}{Z_p}\,Q_{\text{self}}[\sigma_{nn}]$,
$Q_{\text{god}}[v_{t_i}] = Q_{\text{self}}[v_{t_i}] + \tfrac{1}{Z_s}\,Q_{\text{self}}[\sigma_{nt_i}]$,
and $Q_{\text{god}}[\sigma_{nn}] = Q_{\text{god}}[\sigma_{nt_i}] = 0$.

**This is mathematically the same as the γ-mirror result** for a
planar free surface (proof: Riemann between Q_self and γ·Q_self
produces exactly this Q_god at the interface).  The divergence with
γ-mirror appears only when `BuildFrame` — internal to the Riemann
solver — picks a (t1, t2) choice that differs from the ideal local
frame, and only at ≥ O(h) sensitivity around corners.

**2. Method signature.**
```cpp
/// Free-surface flux via Godunov characteristic projection.
/// Replaces the γ-mirror ghost-cell state with a direct construction
/// of the σ·n=0 imposed state, bypassing the frame-choice sensitivity
/// of BuildFrame at corner geometries.
///
/// @param[in]  nor     Unit outward face normal (3 components).
/// @param[in]  Q_self  State on the local side (9 components, global frame).
/// @param[out] F_h     Numerical flux (9 components, global frame).
void GodunovFlux::FreeSurfaceGodunov(const real_t *nor,
                                     const real_t *Q_self,
                                     real_t *F_h) const;
```

**3. Implementation.**
```cpp
void GodunovFlux::FreeSurfaceGodunov(const real_t *nor,
                                     const real_t *Q_self,
                                     real_t *F_h) const
{
   // 1. Build orthonormal frame (Gram-Schmidt matching Interior).
   real_t t1[3], t2[3];
   BuildFrame(nor, t1, t2);

   // 2. Rotate into face-local frame.
   DenseMatrix Tinv(NUM_STATE, NUM_STATE);
   DenseMatrix T   (NUM_STATE, NUM_STATE);
   BuildRotationInverse(nor, t1, t2, Tinv);
   BuildRotation       (nor, t1, t2, T   );

   real_t Q_rot[NUM_STATE];
   Tinv.Mult(Q_self, Q_rot);

   // 3. Apply Godunov projection directly (no ghost cell; no Riemann
   //    solve against γ-mirror).  Enforces σ·n = 0 exactly; velocity
   //    is incremented by Z^{-1} · traction_self.
   const real_t invZp = 1.0 / Zp_;
   const real_t invZs = 1.0 / Zs_;

   real_t Q_god_rot[NUM_STATE];
   std::memcpy(Q_god_rot, Q_rot, NUM_STATE * sizeof(real_t));
   Q_god_rot[SXX] = 0.0;   // sigma_nn  = 0 at free surface
   Q_god_rot[SXY] = 0.0;   // sigma_nt1 = 0
   Q_god_rot[SXZ] = 0.0;
   Q_god_rot[VX] += invZp * Q_rot[SXX];
   Q_god_rot[VY] += invZs * Q_rot[SXY];
   Q_god_rot[VZ] += invZs * Q_rot[SXZ];

   // 4. Apply A_x^+ to the (self, god) pair in the rotated frame.
   //    Since Q_god has zero traction by construction, only the
   //    velocity contribution of A_x^- survives, and it couples to
   //    the Q_self traction — identical to the γ-mirror result at
   //    flat surfaces but not sensitive to BuildFrame's (t1, t2)
   //    choice at corners.
   real_t F_rot[NUM_STATE];
   ApplySplitFlux(Q_rot, Q_god_rot, F_rot);

   // 5. Rotate back.
   T.Mult(F_rot, F_h);
}
```

**4. Edge cases.**
- `nor` not unit: up to the caller (matches `FreeSurface` contract).
- Degenerate `Zp_ = 0` or `Zs_ = 0`: rejected at ctor time by
  `MFEM_VERIFY(rho > 0 && cp > 0 && cs > 0)`.

### Interfaces
- Public `FreeSurfaceGodunov` on `GodunovFlux`.
- No new data members beyond the impedance fields already present.

### Acceptance Criteria
- [ ] `make seas_tpv102_driver` compiles.
- [ ] All existing `seas_test_godunov_*` tests pass.
- [ ] New `seas_test_free_surface_godunov`:
  - On `nor = (0, 0, 1)` (axis-aligned z=0), for each of 50 random
    `Q_self` states, `F_h_gamma` and `F_h_godunov` agree to 10 ULP.
  - On a 10°-tilted normal, they agree to 50 ULP.
  - On a 45°-tilted + non-axis-aligned `(t1, t2)` choice from
    `BuildFrame`, they still agree within 100 ULP (γ-mirror's
    frame-invariance on a flat surface is exact; Godunov is exact by
    construction).
- [ ] The ctor-built `R_11` and `R_21` compliance block match
  hand-derived `diag(-1/Zp, -1/Zs, -1/Zs)` bit-exactly (tested via
  `EXPECT_DOUBLE_EQ`).

### Dependencies
- Depends on: nothing new (uses the existing ctor-built `R`, `Zp_`, `Zs_`).
- Required by: Phase 2.

---

## Phase 2: Wire `FreeSurfaceGodunov` into the wave operator via `--free-surface-bc`

### Goal
After this phase, `seas_tpv102_driver --free-surface-bc=godunov` runs
the TPV102 setup with the Godunov-projection free-surface BC at every
z=0 face.  Default remains γ-mirror.

### Files to Modify
- `dynamic/wave_operator.hpp` — add `enum class FreeSurfaceBCMode { Gamma, Godunov }`
  and a `free_surface_bc_mode_` member.  Set via a new setter.
- `dynamic/wave_operator.inl` — in `ComputeFaceFluxRHS` (and
  `ComputeSharedFaceFluxRHS` if any shared boundary face is possible;
  not in standard TPV102 but defensible) branch on
  `free_surface_bc_mode_` at the `case FaceBC::FreeSurface` line.
- `drivers/tpv102_driver.cpp` — parse the new CLI flag; call
  `wave.SetFreeSurfaceBCMode(mode)`.

### Detailed Requirements

**1. `WaveOperator::SetFreeSurfaceBCMode`:**
```cpp
enum class FreeSurfaceBCMode { Gamma = 0, Godunov = 1 };

void SetFreeSurfaceBCMode(FreeSurfaceBCMode m)
{ free_surface_bc_mode_ = m; }
```
Default = `Gamma`.

**2. Dispatch in `ComputeFaceFluxRHS` (around line 652 of `wave_operator.inl`):**
```cpp
case FaceBC::FreeSurface:
   if (free_surface_bc_mode_ == FreeSurfaceBCMode::Godunov) {
      flux_.FreeSurfaceGodunov(nor, Q_self, F_h);
   } else {
      flux_.FreeSurface(nor, Q_self, F_h);
   }
   break;
```
Same pattern in `ComputeSharedFaceFluxRHS` free-surface branch (if
reached — rare in practice).

**3. CLI flag in `drivers/tpv102_driver.cpp`:**
```cpp
std::string fs_bc_str = GetStringArg(argc, argv,
                                     "--free-surface-bc", "gamma");
FreeSurfaceBCMode fs_bc_mode = (fs_bc_str == "godunov")
                                 ? FreeSurfaceBCMode::Godunov
                                 : FreeSurfaceBCMode::Gamma;
wave.SetFreeSurfaceBCMode(fs_bc_mode);
if (rank == 0) {
   std::cout << "Free-surface BC: " << fs_bc_str << "\n";
}
```

**4. Edge cases.**
- `--free-surface-bc=xyz` (unrecognised) ⇒ fall through to default `gamma`
  and print a warning on rank 0.  Do not abort (per CLAUDE.md "Do not
  silently fall back" — the warning satisfies the "report" clause
  without forcing a run-time abort).
- Under MPI: every rank reads the same CLI; no broadcast needed.

### Acceptance Criteria
- [ ] `--free-surface-bc=gamma` (default) produces byte-identical
  station output to rev 3 on `tpv102_1000m.msh` / 1 rank / 2 s — verified
  by byte-compare.
- [ ] `--free-surface-bc=godunov` runs to completion; station σ_n at
  `flt_0_7.5` remains within ±1 MPa of pre-stress for at least t ≤ 2 s
  (the early-rupture window BEFORE the first reflection corner effect).
- [ ] New integration test `seas_test_free_surface_godunov_driver`:
  end-to-end TPV102 1-element mesh, 10 RK4 steps at dt=1e-4;
  `max|Q[SYY]|` at the fault QP stays ≤ 1 μPa under both γ and godunov
  modes (locked fault + quiescent initial condition ⇒ no energy to
  radiate ⇒ σ_yy should not grow).

### Dependencies
- Depends on: Phase 1.
- Required by: nothing in v9.3.0.  Optional gate for the "corner pump"
  H-V92-G1 diagnostic.

---

## Phase 3: `FaultFaceFlux::EvaluateTotal` — total-stress friction pipeline (I-06 part A)

### Goal
After this phase, a new `FaultFaceFlux::EvaluateTotal` method exists
that implements Pelties eq. 7–12 on **total-stress inputs**, not
fluctuations.  `Evaluate` may remain temporarily for migration checks,
but TPV102 no longer dispatches to it after Phase 4.  Unit-tested for
equivalence with `Evaluate` on a well-posed fixture.

### Files to Modify
- `dynamic/fault_face_flux.hpp` — add `EvaluateTotal` declaration.  Add
  a new optional `DOFData` field `struct InitialStress { real_t sigma_nn, tau_nt1, tau_nt2; }`
  for compactness, OR reuse existing `sigma_n0, tau1_0, tau2_0` — the
  plan uses the existing fields (no ABI change).
- `dynamic/fault_face_flux.cpp` — add definition (~80 LOC).

### Files to Create
- `tests/unit/test_fault_face_flux_total_vs_fluctuation.cpp` — prove
  equivalence on a range of inputs.

### Detailed Requirements

**1. Math — mapping between fluctuation and total at the fault QP.**

Let $Q^{\pm}_{\text{fluc}}$ be the fluctuation state (MFEM's current
Q) and $Q^{\pm}_{\text{tot}} = Q^{\pm}_{\text{fluc}} + Q^{\pm}_{\text{pre}}$
be the total, where $Q^{\pm}_{\text{pre}}$ is the constant pre-stress
state with $\sigma_{nn,\text{pre}} = \sigma_{n,0}$,
$\sigma_{nt_{1},\text{pre}} = \tau_{1,0}$,
$\sigma_{nt_{2},\text{pre}} = \tau_{2,0}$, and zero velocity /
passive-mode components.  Pre-stress is assumed SAME on both sides for
a symmetric fault (TPV102 case).

Under the total formulation, Pelties eq. 7 becomes:
$$
\sigma_{n}^{*,\text{tot}}
= \eta_{p}\Big(v_{n}^{-} - v_{n}^{+}
   + \tfrac{\sigma_{nn,\text{tot}}^{+}}{Z_{p}^{+}} + \tfrac{\sigma_{nn,\text{tot}}^{-}}{Z_{p}^{-}}\Big).
$$
Note the sum term pre-stresses cancel across $\pm$ sides identically
(symmetric fault), so $\sigma_{n}^{*,\text{tot}} = \sigma_{n}^{*,\text{fluc}} + 2\eta_{p}\sigma_{n,0}/Z_{p} = \sigma_{n}^{*,\text{fluc}} + \sigma_{n,0}$
(using $\eta_{p} = Z_{p}/2$).  **Result:** the total trial is the
fluctuation trial PLUS the pre-stress, by construction.

The friction solve input $\Theta$ then becomes the total shear-traction
magnitude directly — identical to SeisSol's
`sqrt(totalTraction1^2 + totalTraction2^2)` at `RateAndState.h:222-225`
— eliminating the `tau_i_total = tau_i_0 + tau_i*` reconstruction step
in MFEM.

**2. Method signature.**
```cpp
/// Total-stress variant of Evaluate.  Expects Q_plus, Q_minus to
/// carry TOTAL stresses (pre-stress + fluctuation).  Output Q_imp_±
/// also carry TOTAL stresses.  Mathematically equivalent to
/// Evaluate(data, Q_fluc_±, ...) + pre-stress shift; see plan §3.
///
/// @param[in,out] data  Per-DOF state.  sigma_n_corr, tau1_corr,
///                      tau2_corr stored as TOTAL values (same as
///                      the fluctuation path already does).
/// @param[in]  Q_plus_tot  Total-stress state on + side.
/// @param[in]  Q_minus_tot Total-stress state on - side.
/// @param[out] Q_imp_plus_tot  Total-stress imposed state on + side.
/// @param[out] Q_imp_minus_tot Total-stress imposed state on - side.
/// @param[in]  method Friction solver choice.
void EvaluateTotal(DOFData &data,
                   const real_t *Q_plus_tot,
                   const real_t *Q_minus_tot,
                   real_t *Q_imp_plus_tot,
                   real_t *Q_imp_minus_tot,
                   FrictionSolver::Method method = FrictionSolver::Method::Brent) const;
```

**3. Implementation.**
```cpp
void FaultFaceFlux::EvaluateTotal(DOFData &data,
                                   const real_t *Q_plus, const real_t *Q_minus,
                                   real_t *Q_imp_plus, real_t *Q_imp_minus,
                                   FrictionSolver::Method method) const
{
   // Homogeneous check identical to Evaluate.
   auto homog_ok = [](real_t a, real_t b) {
      return std::abs(a - b) <= 1e-12 * std::max(std::abs(a), std::abs(b));
   };
   MFEM_VERIFY(homog_ok(data.Zp_plus, data.Zp_minus) &&
               homog_ok(data.Zs_plus, data.Zs_minus),
               "Bimaterial fault face: EvaluateTotal assumes homogeneous.");

   // Step 1: Trial traction IS the total traction (Pelties eq. 7 on
   //         total states — not on fluctuations + pre-stress).
   //         SeisSol matches: RateAndState.h:222-240.
   real_t sigma_n_trial, tau1_trial, tau2_trial;
   ComputeTrialTraction(data, Q_plus, Q_minus,
                        sigma_n_trial, tau1_trial, tau2_trial);

   // Step 2: Solve friction equation for |V| using total Theta.
   //         No `tau_i_0 + tau_i*` reconstruction needed.
   const real_t Theta = std::sqrt(tau1_trial * tau1_trial
                                 + tau2_trial * tau2_trial);
   real_t V_abs = 0.0;
   if (Theta > 0.0)
   {
      V_abs = solver_.Solve(Theta, data.psi, std::abs(sigma_n_trial),
                            data.eta_s, data.a, method);
   }

   // Step 3: Slip-rate decomposition (identical formula; input is total).
   real_t V1 = 0.0, V2 = 0.0;
   real_t tau1_corr = tau1_trial, tau2_corr = tau2_trial;
   if (Theta > 0.0 && V_abs > 0.0)
   {
      real_t C = std::exp(data.psi / data.a) / (2.0 * FrictionSolver::V0);
      real_t f_V = data.a * std::asinh(V_abs * C);
      real_t strength = std::abs(sigma_n_trial) * f_V;
      V1 = V_abs * tau1_trial / (strength + data.eta_s * V_abs);
      V2 = V_abs * tau2_trial / (strength + data.eta_s * V_abs);
      tau1_corr = tau1_trial - data.eta_s * V1;
      tau2_corr = tau2_trial - data.eta_s * V2;
   }

   // Step 4: Imposed states in TOTAL — note the subtraction pattern
   //         matches the fluctuation path verbatim because the
   //         differences (sigma_n_corr - Q_plus[SXX]) are fluctuation-
   //         invariant under the pre-stress shift.
   real_t sigma_n_corr = sigma_n_trial;
   std::memcpy(Q_imp_minus, Q_minus, NUM_STATE * sizeof(real_t));
   std::memcpy(Q_imp_plus,  Q_plus,  NUM_STATE * sizeof(real_t));

   const real_t invZp = 1.0 / data.Zp_plus;  // homogeneous
   const real_t invZs = 1.0 / data.Zs_plus;

   Q_imp_minus[VX] = Q_minus[VX] - invZp * (sigma_n_corr - Q_minus[SXX]);
   Q_imp_minus[VY] = Q_minus[VY] - invZs * (tau1_corr    - Q_minus[SXY]);
   Q_imp_minus[VZ] = Q_minus[VZ] - invZs * (tau2_corr    - Q_minus[SXZ]);
   Q_imp_plus[VX]  = Q_plus[VX]  + invZp * (sigma_n_corr - Q_plus[SXX]);
   Q_imp_plus[VY]  = Q_plus[VY]  + invZs * (tau1_corr    - Q_plus[SXY]);
   Q_imp_plus[VZ]  = Q_plus[VZ]  + invZs * (tau2_corr    - Q_plus[SXZ]);

   Q_imp_minus[SXX] = sigma_n_corr;  Q_imp_plus[SXX] = sigma_n_corr;
   Q_imp_minus[SXY] = tau1_corr;     Q_imp_plus[SXY] = tau1_corr;
   Q_imp_minus[SXZ] = tau2_corr;     Q_imp_plus[SXZ] = tau2_corr;

   // Step 5: Update DOFData — store TOTAL directly (no sigma_n0 add).
   data.slip_rate    = V_abs;
   data.V1           = V1;
   data.V2           = V2;
   data.tau1_corr    = tau1_corr;    // TOTAL (not tau1_0 + tau1_corr)
   data.tau2_corr    = tau2_corr;
   data.sigma_n_corr = sigma_n_corr;
}
```

**4. Edge cases.**
- `Theta = 0`: identical to fluctuation path (both `Evaluate` and
  `EvaluateTotal` short-circuit to `V_abs = 0, V_1 = V_2 = 0`).
- Bimaterial: rejected with `MFEM_VERIFY`.

### Acceptance Criteria
- [ ] Compiles; existing `seas_test_fault_face_flux*` unchanged.
- [ ] New `test_fault_face_flux_total_vs_fluctuation`:
  - For 1000 random (Q_plus_fluc, Q_minus_fluc, pre-stress, ψ, a)
    fixtures, compute `Q_imp_fluc = Evaluate(Q_fluc)` and
    `Q_imp_tot  = EvaluateTotal(Q_fluc + pre)`.  Assert
    `|Q_imp_fluc + pre - Q_imp_tot| ≤ 1e-10 · |pre|` per component.
  - Same for `data.V1, V2, slip_rate` (bit-identical modulo solver
    tolerance ~1e-8).
  - Same for `data.{sigma_n_corr, tau_i_corr}` — both paths store
    TOTAL, so they must match to solver tolerance.
- [ ] The test above is documented as a migration-reference gate only;
  it does not imply that fluctuation mode remains a supported TPV102
  runtime option after Phase 4.

### Dependencies
- Depends on: nothing.
- Required by: Phase 4.

---

## Phase 4: Total-stress initialization + TPV102 dispatch migration (I-06 part B)

### Goal
After this phase, TPV102 initializes Q in total stress unconditionally.
There is no `--state-representation` CLI, and TPV102 fault dispatch
uses `EvaluateTotal` exclusively.  The bulk state carries the
pre-stress field at every DOF (not just fault DOFs).

### Files to Modify
- `dynamic/wave_operator.inl::Mult` — replace the TPV102 fault-face
  dispatch so it calls `EvaluateTotal`.
- `dynamic/wave_operator.inl::ComputeSharedFaceFluxRHS` — same
  replacement for shared fault faces.
- `drivers/tpv102_driver.cpp` — initialize TPV102 in total stress
  unconditionally and remove any design mention of a dual
  fluctuation/total runtime mode.

### Files to Create
- `dynamic/tpv102_setup_total.hpp` — companion to `tpv102_setup.hpp`
  exposing `InitializeStateTotal(Q, ndof_total, pre_stress_at_dof)`
  which fills every DOF's Q[SXX/SYY/SZZ/SXY/…] with the rotated
  pre-stress tensor in GLOBAL coordinates.  For TPV102 (homogeneous
  half-space under a vertical fault) the pre-stress is:
  ```
  sigma_xx = 0,   sigma_yy = -sigma_n0,   sigma_zz = 0,
  sigma_xy = -tau_ini,   sigma_yz = 0,   sigma_xz = 0
  ```
  — i.e. compressive normal stress of 120 MPa on the y=const plane
  plus shear of 75 MPa on the x-y plane.

### Detailed Requirements

**1. `Mult` replacement.** At the fault-face `if (is_fault)` block:
```cpp
fault_flux_->EvaluateTotal(fdata, Q_plus_can, Q_minus_can,
                           Q_imp_plus_can, Q_imp_minus_can, method);
```
Same replacement in `ComputeSharedFaceFluxRHS`.

**2. `InitializeStateTotal`:**
```cpp
inline void InitializeStateTotal(Vector &Q, int ndof_total,
                                 real_t sigma_n0, real_t tau_ini)
{
   Q.SetSize(NUM_STATE * ndof_total);
   Q = 0.0;
   // Pre-stress: sigma_yy = -sigma_n0 (compression on y=const plane),
   //             sigma_xy = -tau_ini   (right-lateral strike-slip loading).
   // Sign: MFEM uses sigma_nn > 0 = compression.  Global sigma_yy
   //       at the fault y=0 equals sigma_nn under the canonical frame
   //       (see walkthrough I-01 note).  Hence sigma_yy_pre = +sigma_n0.
   // Shear: tau2_0 = tau_ini in canonical frame (strike component).
   //       Global sigma_xy from canonical sigma_nt2 via T_can:
   //       sigma_xy_global = -tau_nt2_local (Voigt rotation with
   //       can_n = (0,-1,0), can_t2 = (+1,0,0)).
   //       So sigma_xy_pre = -tau_ini.
   for (int i = 0; i < ndof_total; i++) {
      Q[SYY * ndof_total + i] =  sigma_n0;
      Q[SXY * ndof_total + i] = -tau_ini;
   }
}
```

**3. Driver plumbing.**
```cpp
Vector Q;
InitializeStateTotal(Q, ndof_total,
                     TPV102Params::sigma_n, TPV102Params::tau_ini);
// Pre-stress is already baked into Q — zero the DOFData pre-stress
// fields to avoid double-counting in EvaluateTotal.
for (int i = 0; i < num_fault_total; i++) {
   dof_data[i].sigma_n0 = 0.0;
   dof_data[i].tau1_0   = 0.0;
   dof_data[i].tau2_0   = 0.0;
}
```

**4. Nucleation compatibility.**  `ApplyNucleation` writes to
`dof_data[i].tau2_0`.  After the migration, `tau2_0` is zeroed and
must remain bookkeeping-only.  Nucleation therefore has to write to the
BULK Q at the fault QP, not to the DOFData pre-stress field.  Add an
`ApplyNucleationTotal`
variant:
```cpp
inline void ApplyNucleationTotal(Vector &Q, ...,
                                 const std::vector<Vector> &fault_coords,
                                 int ndof_total, real_t t)
{
   // Option A: per-fault-QP sigma_xy injection (identifies which
   //           global Q DOF corresponds to which fault QP via a
   //           fault_qp_to_global_dof_ map built once in the driver).
   // Phase 4 may land before the full nucleation-to-Q injection from
   // Phase 5.  In that case, the Phase 4-only tests use a locked-fault
   // fixture with nucleation amplitude set to zero inside the harness.
}
```
Phase 5 adds the proper nucleation-to-Q injection for the production
TPV102 path.

### Acceptance Criteria
- [ ] There is no `--state-representation` parser in the final TPV102
  driver path.
- [ ] A 1-element locked-fault fixture with nucleation suppressed in
  the test harness runs to completion over 100 RK4 steps and station
  σ_n stays at 120 MPa to 1 Pa.
- [ ] Both interior and shared TPV102 fault paths call `EvaluateTotal`
  and no TPV102 dispatch path calls `Evaluate`.
- [ ] Existing v9.1.0 / v9.2.0 regression tests all pass.

### Dependencies
- Depends on: Phase 3.
- Required by: Phase 5 (nucleation), Phase 6 (tests), Phase 7 (ADER integration).

---

## Phase 5: Nucleation via bulk-Q injection (I-06 part C)

### Goal
After this phase, the migrated total-stress TPV102 path works with the
real nucleation patch.  The perturbation is injected into the bulk Q at
the fault QPs instead of into `DOFData.tau2_0`.

### Files to Modify
- `dynamic/tpv102_setup_total.hpp` — add `ApplyNucleationTotal`.
- `drivers/tpv102_driver.cpp` — call `ApplyNucleationTotal` instead of
  `ApplyNucleation` in the TPV102 driver at each RK4 stage.

### Detailed Requirements

**1. Fault-QP-to-global-DOF map.**  Build once at init:
```cpp
// For each fault QP index i, store the (element, local_dof, global_DOF_offset)
// tuple so nucleation can write directly into Q.  Use the face-local
// basis evaluation at the QP — Σ_j shape_j(ip) · Q[SXY, j] must equal
// -tau_ini - dtau(x,z,t).  The simplest implementation uses the L2
// nodal basis on the face which admits a direct interpolation.
std::vector<int> fault_qp_to_elem_;   // element index for each fault QP
std::vector<int> fault_qp_to_elem_dof_;// local DOF index in that element
```
This requires inspecting the L2 basis at the fault QP; for GaussLobatto
the nodal basis is interpolatory at the Lagrange node closest to the
QP — choose that one.  For the general case, solve the 3×3 local
interpolation system.

**2. Projection step.**  At each RK4 stage call, for each fault QP i
with perturbation `dtau(i, t)`:
```cpp
// In the canonical fault-local frame, nucleation adds to tau_nt2
// (strike shear).  In global frame this maps to sigma_xy (see I-03
// note).  The injection must PRESERVE the existing Q value elsewhere
// — we're ADDING dtau to Q[SXY, ...] at this QP only.
//
// Because Q is stored per ELEMENT (not per QP), and the QP lies
// strictly INSIDE the element, we use the local L2 Lagrange basis:
//   sigma_xy(ip) = Σ_j shape_j(ip) * Q[SXY, elem_dof(j)]
// Adding dtau at the QP requires adding dtau * M_ref^-1 * shape to
// the element DOFs, where M_ref is the element face-to-element
// coupling mass.
//
// Simpler approach: use a nodal GaussLobatto basis where shape_j(ip)
// is delta_{j, nearest_node}.  Then injection is a single DOF update.
```
For the Phase 5 minimum, require **GaussLobatto basis + interpolatory
nodal injection** and MFEM_VERIFY this at setup.

**3. Alternative if the interpolatory pick is infeasible:**
L2 project the nucleation field `dtau(x, z, t) · ê_xy · ξ(x, z)` (with
ξ the spatial nucleation profile) onto the element DG space at every
RK4 stage using MFEM's `LinearForm` integration.  More expensive but
avoids QP-vs-node alignment concerns.

### Acceptance Criteria
- [ ] The default TPV102 run (no state-representation flag) with
  nucleation ON runs to completion on the 1000 m fixture, 2 s wall-time.
- [ ] Station σ_n at hypocenter matches a frozen fluctuation-reference
  migration baseline within solver tolerance (1 kPa) over t ∈ [0, 2 s].
- [ ] Slip_strike at hypocenter matches the same migration baseline to
  within 1 mm over t ∈ [0, 2 s].
- [ ] No regression in the v9.2.0 `_k4` VTU diagnostic.

### Dependencies
- Depends on: Phase 4.
- Required by: Phase 6 integration tests.

---

## Phase 6: Migration + coverage tests (I-06 validation)

### Goal
After this phase, a regression matrix covers `{gamma, godunov}` on a
1-element TPV102 fixture for 100 RK4 steps, and a dedicated migration
test checks the new total-stress path against a frozen fluctuation
reference on a small fixture.  Fluctuation remains a reference only,
not a supported TPV102 runtime mode.

### Files to Create
- `tests/unit/test_tpv102_total_fluctuation_equivalence.cpp` — 1-element
  TPV102 rig, V_ini initial condition, 100 RK4 steps at dt=1e-4,
  assertion: station σ_n, slip_strike, slip_rate_strike from the new
  total-stress path agree with a frozen fluctuation-reference fixture
  to 1 Pa / 1 nm / 1 nm/s respectively.
- `tests/unit/test_tpv102_free_surface_bc_variants.cpp` — run the same
  rig with γ vs godunov free surface, assert agreement to 10 ULP on a
  flat surface and ≤ 10 Pa on the corner fixture (tolerant of the
  actual corner-pumping divergence we HOPE to see eliminated).

### Detailed Requirements
Both cells of the 1×2 matrix `{gamma, godunov}` run to completion on
the minimum fixture.  Output comparison is performed in-test via direct
DOFData inspection plus, if paraview is enabled, a file-level diff of
the VTU.  The migration regression against the fluctuation reference is
kept separate so runtime support is not confused with validation.

### Acceptance Criteria
- [ ] Both `{gamma, godunov}` cells produce run-complete TPV102 rigs;
  no NaN / abort / CFL overrun.
- [ ] The 1×2 grid of (σ_n at hypocenter at t=1 s) values agrees within
  1 Pa across `{gamma, godunov}` on the flat-surface fixture.
- [ ] The total-vs-frozen-fluctuation migration regression agrees to
  1 Pa in σ_n and 1 nm in slip_strike at t=1 s.

### Dependencies
- Depends on: Phases 1, 4, 5.
- Required by: v9.3.0 release.

---

## Phase 7: Wire v9.3.0 into the ADER plan (I-05 integration point)

### Goal
After this phase, the ADER plan's `--time-integrator=ader` flag (from
`tpv102_ader_time_integration_plan.md` Phase 7) combines cleanly with
the total-stress-only TPV102 path from I-06.  The combined mode
reproduces SeisSol's outer-loop behaviour.

### Files to Modify
- `tpv102_ader_time_integration_plan.md` — add a cross-reference note
  to v9.3.0 Phase 4 as the required state mode for ADER in TPV102.  Do NOT
  duplicate Phase 4 content.

### Files to Create
- `tests/unit/test_ader_total_vs_rk4_total.cpp` — a 1-element TPV102
  rig, run both combinations for 10 time steps,
  assert station output matches to O(dt) (ADER-O(2) is 2nd-order;
  RK4 is 4th-order; difference scales as dt² at the leading order
  given identical physics).

### Detailed Requirements

**1. Compatibility statement.**

Under ADER, the `ComputeADERTimeIntegrated` predictor from
ADER Phase 3 operates on total-stress Q directly — no modification
needed to the CK recursion because it only uses spatial derivatives of
Q and the constant material Jacobians $A_{d}$.

**2. Station-output alignment.**  Under ADER, the step-end DOFData is
the one-shot friction-solve result at $t_{n+1}$ (ADER plan §7).  Under
RK4 + total, DOFData after the averaging block is the RK4-average over
the step.  Under ADER + total, DOFData is the exact end-of-step
quantity at $t_{n+1}$ — this is the expected SeisSol-parity output.

### Acceptance Criteria
- [ ] `--time-integrator=ader` runs on the TPV102 1000 m fixture for
  2 s wall-time, 1 rank, using the total-stress TPV102 path.
- [ ] Station σ_n at hypocenter matches the RK4+total baseline
  to within 100 Pa (0.1 ‰ of σ_n0) for t ∈ [0, 2 s].
- [ ] The `_k4` VTU diagnostic (from plan v9.2.0 §19) shows
  `max |normal_stress_k4 − normal_stress|` ≤ 1 Pa — confirming that
  ADER's one-shot DOFData does not carry the RK4 averaging artifact.

### Dependencies
- Depends on: Phase 4, Phase 6, ADER plan Phase 7.
- Required by: v9.3.0 final commit.

---

## Testing Strategy

### Per-phase gates
Each phase has its own unit test(s); the existing regression suite
(`seas_test_fault_surface_vtu_continuity`, `seas_test_fault_surface_vtu_k4`,
`seas_test_fault_basis_dip_strike_symmetry`, `seas_test_godunov_rotation_identity`)
must stay green after every phase.

### Aggregate targets
Add to `Makefile`:
```
test-v93-regression: seas_test_free_surface_godunov \
                     seas_test_fault_face_flux_total_vs_fluctuation \
                     seas_test_tpv102_total_fluctuation_equivalence \
                     seas_test_tpv102_free_surface_bc_variants \
                     seas_test_ader_total_vs_rk4_total
```

### Reference comparison (manual, post-landing)
After Phase 7, a local 2 s / 1 rank run with the total-stress TPV102
path and all SeisSol-parity options enabled should produce station
output whose σ_n stays within ±1 MPa of
pre-stress at `flt_0_7.5` and whose slip_dip stays ≤ 0.001 m.  If it
does, the three divergences I-04/05/06 are sufficient to reproduce
SeisSol behaviour and the remaining v9.2 symptoms (if any) fall into
I-07 (shared-fault MPI) or a residual mesh/nucleation concern.

## Risk Assessment

### High risk
- **Phase 4 locked-fault harness.**  Using a locked fixture before the
  full nucleation-to-Q injection lands can look deceptively healthy.
  Mitigation: Phase 4 acceptance includes a pulse test where a small
  initial velocity perturbation is injected manually and the decay rate
  compared to a known analytical Green's function.

- **Phase 5 nucleation-to-Q injection on non-nodal bases.**  If the DG
  basis is modal / non-interpolatory, the single-DOF injection is
  wrong.  Mitigation: require GaussLobatto at Phase 5 entry; defer
  the L2-projection path to a v9.3.1 add-on if a non-nodal basis is
  introduced.

### Medium risk
- **Phase 7 ADER + total interaction with fault-flux injection.**
  Under ADER, F_h is computed from the time-integrated state I, not
  Q.  The fault flux's imposed-state construction expects input in
  the same representation as Q.  Under total + ADER, everything stays
  total — the `EvaluateADERTotal` variant (not in v9.3.0 scope; add
  to ADER plan Phase 5 as a follow-up) is a straightforward
  copy of Phase 3 with `Evaluate` → `EvaluateTotal`.  Until that is
  landed, the v9.3.0 Phase 7 test uses an explicit workaround:
  rescale I by 1/Δt inside `ComputeADERFaceFluxRHS` (ADER plan Phase 6)
  before calling `EvaluateTotal`, then rescale the output by Δt.

### Low risk
- **CLI flag ambiguity.**  `--free-surface-bc=Godunov` (capitalised)
  vs `godunov` — resolve by `std::tolower`-ing the input.
- **BP5 regression.**  BP5 drivers do not parse any v9.3.0 flag, so
  no regression possible.  Verified by a `make test-bp5-smoke`
  after each phase lands.

## Scope and timeline

| Phase | LOC (new) | LOC (modified) | Est. wall-clock |
|---|---:|---:|---:|
| 1 Godunov FS BC      |  80 |  30 | 1 day   |
| 2 BC dispatch flag   |  40 |  60 | 0.5 day |
| 3 EvaluateTotal      |  80 |  20 | 1 day   |
| 4 Total-Q init+dispatch| 120 |  80 | 1.25 day |
| 5 Nucleation→bulk Q  | 100 |  50 | 1 day   |
| 6 Equivalence tests  | 250 |  30 | 1 day   |
| 7 ADER wiring        |  60 |  40 | 0.5 day |
| **Total**            | **~750** | **~320** | **~6.5 days** |

## Deferred / explicitly NOT in scope

- **Heterogeneous material at fault.**  Both `Evaluate` and
  `EvaluateTotal` assume `Z_p^+ = Z_p^-`; bimaterial is a v9.4 topic.
- **Godunov projection on a TILTED free surface with curvature.**
  Phase 1 Godunov projection is frame-invariant for flat faces; the
  face-normal is re-built at every QP from the element Jacobian.
  Curved faces (not in TPV102) may need a curvature correction —
  deferred.
- **LTS (Local Time Stepping).**  SeisSol's main production feature.
  Out of scope; use global Δt.
- **`EvaluateADERTotal`.**  Properly a follow-up addition to the ADER
  plan Phase 5, not v9.3.0.  v9.3.0 Phase 7 uses an explicit rescale
  workaround.
- **BP5 total-stress migration.**  BP5 does not benefit (quasi-dynamic
  SIPG does not have a fluctuation/total split in the same sense).

## Quick-reference for implementation agent

- Follow CLAUDE.md "Files Requiring Extreme Care" — `wave_operator.inl`,
  `fault_face_flux.cpp` are on it; any edit requires the full v9.2.0
  `test-v92-regression-gates` suite to stay green AFTER each phase.
- Do not revert any v9.0.0 / v9.1.0 / v9.2.0 fix.  The Godunov-projection
  path is additive; the TPV102 total-stress path REPLACES the
  fluctuation-Q TPV102 runtime path.
- Phase ordering is: Phase 1 → Phase 2 → Phase 3 → Phase 4 → Phase 5 → Phase 6 → Phase 7.
  Phases 1, 2, 3 can land in parallel if desired (disjoint files);
  Phases 4–7 must be sequential.
- Every phase produces a commit with tests included.  The v9.3.0 tag
  is applied ONLY after Phase 7 passes the full regression matrix.
