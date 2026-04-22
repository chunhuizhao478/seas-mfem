# Implementation Plan: ADER Time Integration for SEAS-MFEM TPV102 Driver

> **Author:** planning agent, 2026-04-21
> **Scope:** add ADER (Arbitrary-order DERivative) time integration as an
> opt-in alternative to the current explicit RK4 in `seas_tpv102_driver`,
> following SeisSol's Cauchy-Kovalevskaya (CK) formulation. RK4 remains
> the default.
> **Primary motivation:** SeisSol's ADER-DG naturally produces a
> one-step-updated DOFData — there is no "stage-averaged vs stage-4"
> ambiguity as in RK4 (cf. REVIEW R-V92-E02 and plan v9.2.0 §19 H-V92-K).
> ADER will **eliminate H-V92-K by construction** if H-V92-K is the
> source of the σ_n drift observed in job 7668434.
> **Secondary motivation:** matches the Pelties 2012 benchmark setup
> (Dumbser & Käser 2006, ADER-DG scheme for 3-D isotropic elastic
> waves), enabling a direct point-by-point code comparison with SeisSol.

## Overview

The SEAS-MFEM TPV102 driver currently uses a 4-stage explicit RK4 on the
wave state `Q` plus an ad-hoc RK4-weighted averaging of `DOFData`
scalar fields for output (`drivers/tpv102_driver.cpp:858-874`). This
plan replaces the outer time-stepper — when `--time-integrator=ader`
is passed — with an ADER one-step predictor-corrector scheme:

1. **Predictor (local, per element):** expand `Q(x, t+τ) = Σ τ^k/k! · ∂^k_t Q(x,t)`
   using the Cauchy-Kovalevskaya recursion `∂_t^{k+1} Q = -Σ_d A_d ∂_{x_d} (∂_t^k Q)`.
   Compute the time-integrated state `I(x) = ∫_0^Δt Q(x, t+τ) dτ`.
2. **Corrector (global):** update `Q_new = Q + M^{-1} (∫_V ∇·(A·I) - ∫_∂V F_h·I)`.
   The fault face flux uses a time-integrated friction solve:
   `Evaluate` is called once per time step on `I_+`, `I_-` and returns
   `I_imp_+`, `I_imp_-` — NOT per stage. No stage-averaging of DOFData.

The plan is phased so each phase compiles cleanly, keeps RK4 as default,
passes all existing regression tests, and is independently reviewable.

## Constraints

### Interface constraints (frozen)
- **`WaveOperator::Mult(Q, dQdt)`** — MUST remain semantically identical
  so RK4 path and any existing caller (BP5 driver, unit tests) is
  unchanged. New ADER methods are added alongside.
- **`FaultFaceFlux::Evaluate(data, Q_plus, Q_minus, Q_imp_plus, Q_imp_minus, method)`** —
  MUST remain bit-identical for RK4 callers. ADER adds a new method
  `EvaluateADER` with a distinct signature.
- **`DOFData` struct layout** — MUST NOT grow new required fields that
  are only valid under ADER (BP5 shares the struct). New optional
  fields are permitted under a `#ifdef` or as a tail-of-struct
  addition.
- **PVD / VTU output field names** — MUST remain compatible with
  existing ParaView scripts. ADER writes to the same fields
  (`sigma_n_corr`, `V1`, etc.); the `_k4` diagnostic from plan v9.2.0
  §19 may or may not be populated under ADER (document which).
- **CLI default** — `seas_tpv102_driver` defaults to RK4 unchanged.

### Dependency constraints
- No new external libraries. All ADER math uses MFEM's `DenseMatrix`,
  `Vector`, `IntegrationRule`, and existing
  `GodunovFlux::BuildJacobian` / `BuildRotation`.
- No change to MFEM itself. All new code lives in
  `miniapps/seas/dynamic/`.
- Must preserve `mfem-dev` conda environment build path (no CMake
  changes needed beyond Makefile entries).

### Convention constraints
- Naming: ADER-specific types prefixed with `ADER` (e.g.
  `ADEROrder`, `ADERPredictor`, `ADERIntegrator`). C++ files under
  `dynamic/` per existing layout.
- All new unit tests in `tests/unit/test_ader_*.cpp` with the
  `test_fault_*` test-runner pattern already used (`TEST_ASSERT` and
  `TEST_NEAR` macros, main returns non-zero on failure, one output
  directory per test under `miniapps/seas/`).
- Follow CLAUDE.md "Files Requiring Extreme Care" list — `wave_operator.inl`
  is on it; any edit requires the full verification suite.

### Numerical constraints
- **Order consistency**: ADER order `O` should satisfy `O >= p + 1`
  where `p` is the DG spatial polynomial order, to avoid order
  reduction. Default: `O = p + 1 = 2` for TPV102 P1 (matches P1 time
  accuracy). Support `O = 2, 3, 4` at minimum.
- **CFL**: ADER's linear stability limit matches RK4's for the
  elastic wave system: `Δt_CFL = C_ader · h_min / c_p` with
  `C_ader ≈ 1 / (2p+1)` (Dumbser-Käser 2006 §4, Table 1).
  For P1, `C_ader ≈ 1/3`, same order as RK4.
- **Conservation**: `∫_V Q_new` must equal `∫_V Q - Δt · ∫_∂V F_h(I)` to
  machine precision (mass conservation). Tested in Phase 8 on a
  periodic domain.
- **Equivalence at linear orders**: On a linear problem (no friction,
  just the bulk wave operator), ADER-O and RK4 must produce
  Q_new identical to O(Δt^min(O, 4)) — this is the convergence-order
  test.

## Phase 1: Spatial derivative infrastructure

### Goal
After this phase, `WaveOperator::ApplySpatialDerivative(int dir, const Vector &Q, Vector &dQ_dxdir)`
exists and is unit-tested. It computes the element-local strong-form
spatial derivative `∂_{x_dir} Q` on every DG element without fluxes —
the building block for the ADER CK recursion.

### Files to Create
- `tests/unit/test_wave_operator_spatial_derivative.cpp` — unit test
  for `ApplySpatialDerivative` using a polynomial-exactness manufactured
  solution.

### Files to Modify
- `dynamic/wave_operator.hpp` — add declaration of `ApplySpatialDerivative`.
- `dynamic/wave_operator.inl` — add definition using the stiffness
  matrix `K_d` already built in the ctor for the volume RHS.
- `Makefile` — register the new test.

### Detailed Requirements

**1. Method signature.** Add as a public const method of `WaveOperator<MeshType>`:
```cpp
/// Compute the strong-form element-local spatial derivative
/// `dQ_dxdir = ∂_{x_dir} Q` on every element, with NO inter-element
/// flux coupling (pure element-local operation).  Uses the stiffness
/// matrix `K_d` already assembled in the ctor's volume-integral path.
///
/// @param[in]  dir      Spatial direction: 0 (x), 1 (y), 2 (z).
/// @param[in]  Q        State vector of size NUM_STATE * ndof_total.
/// @param[out] dQ_dxdir Output spatial derivative, same size as Q.
///                       Must be allocated by the caller.
void ApplySpatialDerivative(int dir,
                             const Vector &Q,
                             Vector &dQ_dxdir) const;
```

**2. Implementation.** For each element `e`:
- Let `ndof = fe->GetDof()`, `dof_offset = e * ndof_per_el_`.
- Build the reference-element stiffness integrator for the requested
  direction (MFEM's `DiffusionIntegrator` or a hand-rolled
  quadrature loop over `ir` using `fe->CalcDShape(ip, dshape)`).
- For each of `NUM_STATE = 9` components `c`:
  - Extract `Q_c(e) = Q[c * ndof_total + dof_offset + (0..ndof)]`.
  - Compute `dQ_c = M_e^{-1} · K_d^e · Q_c(e)` where `K_d^e` is the
    element-local stiffness matrix `K_d^e[i,j] = ∫_e φ_i · ∂_{x_dir} φ_j dV`
    and `M_e^{-1}` is the element mass inverse already cached in
    `elem_mass_inv_[e]`.
  - Write result to `dQ_dxdir[c * ndof_total + dof_offset + i]`.

**3. Convention on `K_d^e`.** The existing `ComputeVolumeRHS` at
`wave_operator.inl:489+` already performs the strong-form volume
integral `M^{-1} · ∫_e φ · ∂_d (A·Q) dV`, factored through
`GodunovFlux::BuildJacobian(d, A)`. The stiffness matrix is implicit
(built inline from `fe->CalcDShape` + quadrature). Factor out the
pure-derivative kernel into a helper `BuildElementStiffnessMatrix(int dir, int e, DenseMatrix &K_d_e)`
that returns the `ndof_per_el_ × ndof_per_el_` matrix `K_d^e[i,j] = ∫_e φ_i · ∂_{x_dir} φ_j dV`
computed via the same quadrature rule MFEM uses in `ComputeVolumeRHS`.

**4. Edge cases to handle.**
- `dir ∉ {0, 1, 2}`: `MFEM_VERIFY(dir >= 0 && dir < 3, ...)` — fail loud.
- `Q.Size() != NUM_STATE * ndof_total_`: `MFEM_VERIFY` — fail loud.
- `dQ_dxdir.Size()` may be 0 on entry: call `dQ_dxdir.SetSize(NUM_STATE * ndof_total_)`.
- PML elements: include, but `ApplySpatialDerivative` must not apply
  PML damping (that is separate from the spatial derivative).

### Interfaces
- Public: `ApplySpatialDerivative`.
- Private helper: `BuildElementStiffnessMatrix(int dir, int e, DenseMatrix &K) const`
  (may be a file-local free function if preferred).

### Edge Cases to Handle
- Empty mesh (`ne_ == 0`): `dQ_dxdir` is zero-length after `SetSize(0)`. No-op.
- Parallel mesh: only element-local data used; no MPI communication.

### Acceptance Criteria
- [ ] Compiles under `make seas_tpv102_driver` and `make seas_test_parallel_wave_operator`.
- [ ] Existing `seas_test_wave_operator` passes unchanged (no regression).
- [ ] New test `seas_test_wave_operator_spatial_derivative` passes with
      `dQ/dx_i (x_j^n) = n · x_i^{n-1} · δ_{ij}` for n ∈ {1, 2, 3} polynomial
      exactness on a P2 reference tet (order = max(2, p)).
- [ ] Runs in < 1 s for a 4×4×4 tet mesh at order 1.

### Dependencies
- Depends on: nothing (uses existing `elem_mass_inv_` and `IntRules`).
- Required by: Phase 3.

---

## Phase 2: Reference-element "star" matrices

### Goal
After this phase, a function `GodunovFlux::GetReferenceStarMatrix(int dir) const`
returns the 9×9 matrix `A^*_d = ∂Q/∂t contribution from ∂_{x_d} Q` in the
REFERENCE element frame. This is the ADER analogue of SeisSol's
"star matrices" (ref-frame Jacobian), needed for the CK recursion.

### Files to Modify
- `dynamic/godunov_flux.hpp` — add `GetReferenceStarMatrix` declaration
  and a `std::array<DenseMatrix, 3> ref_star_` precomputed member.
- `dynamic/godunov_flux.cpp` — precompute the star matrices in the
  existing ctor (~20 LOC) using the existing `BuildJacobian(dir, A)` output.

### Detailed Requirements

**1. What is a "star matrix"?** For an isotropic linear elastic cell in
a STRAIGHT-SIDED element (MFEM's default for linear tets), the
reference-frame Jacobian is related to the global-frame Jacobian by
`A^*_d = T^{-T} · A_d` where `T` is the element Jacobian (constant per
cell for linear tets). For a homogeneous material, `A^*_d` is the
same for every element and can be precomputed.

**2. TPV102-specific simplification.** The TPV102 driver uses a
homogeneous isotropic material (`rho, cp, cs` constant) and a
simplex mesh (MFEM tet). The star matrices `A^*_0, A^*_1, A^*_2` can
be precomputed ONCE from `BuildJacobian(d, A_d)` for d ∈ {0, 1, 2}
and stored. The Jacobian transform `T^{-T}` is absorbed into the
per-element stiffness matrices in Phase 1 — so `A^*_d = A_d` (global
frame) is sufficient for TPV102 at this phase. (Follow-up plan to
handle bent elements / heterogeneous material lives under
"Deferred — out of scope for v1".)

**3. Method signature.** Add as `const DenseMatrix& GodunovFlux::GetReferenceStarMatrix(int dir) const`
returning `ref_star_[dir]`.

**4. Ctor precompute.** At the end of `GodunovFlux::GodunovFlux(rho, cp, cs)`:
```cpp
for (int d = 0; d < 3; ++d) { BuildJacobian(d, ref_star_[d]); }
```
`ref_star_` is a `std::array<DenseMatrix, 3>` member.

### Interfaces
- `const DenseMatrix& GetReferenceStarMatrix(int dir) const`.

### Edge Cases to Handle
- `dir ∉ {0, 1, 2}`: `MFEM_VERIFY`. Fail loud.

### Acceptance Criteria
- [ ] All existing `seas_test_godunov_*` tests pass unchanged.
- [ ] `GetReferenceStarMatrix(0)` returns the same matrix as
      `BuildJacobian(0, A)` — verified by a new assertion in
      `seas_test_godunov_flux`.

### Dependencies
- Depends on: nothing (pure addition).
- Required by: Phase 3.

---

## Phase 3: Cauchy-Kovalevskaya predictor

### Goal
After this phase, `WaveOperator::ComputeADERTimeIntegrated(Q, dt, order, I)`
computes `I = ∫_0^dt Q(t+τ) dτ` at every DOF via the CK recursion, with
no face flux coupling. Unit-tested against a manufactured plane-wave
solution with known `I`.

### Files to Modify
- `dynamic/wave_operator.hpp` — add `ComputeADERTimeIntegrated` declaration.
- `dynamic/wave_operator.inl` — add definition using Phase 1
  `ApplySpatialDerivative` and Phase 2 `GetReferenceStarMatrix`.
- `Makefile` — new test entry.

### Files to Create
- `tests/unit/test_ader_ck_predictor.cpp` — convergence test of
  `ComputeADERTimeIntegrated` on a known plane-wave analytical
  solution.

### Detailed Requirements

**1. Math.** The strong form is
`∂Q/∂t = -A_x ∂Q/∂x - A_y ∂Q/∂y - A_z ∂Q/∂z ≡ L(Q)`.
Define `D(k) = ∂^k Q/∂t^k`. Recursion:
- `D(0) = Q`.
- `D(k+1) = L(D(k)) = -Σ_d A_d · (∂_{x_d} D(k))` for k = 0, 1, ..., O-2.

Time-integrated state:
- `I = Σ_{k=0}^{O-1} (dt^{k+1}/(k+1)!) · D(k)`.

For O=2: `I = dt·Q + (dt²/2) · L(Q)`.
For O=3: `I = dt·Q + (dt²/2) · L(Q) + (dt³/6) · L(L(Q))`.

**2. Method signature.** Public method of `WaveOperator<MeshType>`:
```cpp
/// ADER Cauchy-Kovalevskaya time-integrated state.
/// Computes I = ∫_0^dt Q(t+τ) dτ via Taylor expansion and
/// recursive spatial derivatives.
/// @param[in]  Q       State at t, size NUM_STATE * ndof_total.
/// @param[in]  dt      Time step.
/// @param[in]  order   ADER order O in {2, 3, 4}.  O=p+1 recommended
///                     where p is the DG polynomial order.
/// @param[out] I       Time-integrated state, size NUM_STATE * ndof_total.
///                     Units: [Q] × time.
void ComputeADERTimeIntegrated(const Vector &Q,
                                real_t dt,
                                int order,
                                Vector &I) const;
```

**3. Implementation.**
```cpp
void WaveOperator<MeshType>::ComputeADERTimeIntegrated(
   const Vector &Q, real_t dt, int order, Vector &I) const
{
   MFEM_VERIFY(order >= 2 && order <= 4, "ADER order must be in {2, 3, 4}");
   MFEM_VERIFY(Q.Size() == NUM_STATE * ndof_total_, "Q size mismatch");
   I.SetSize(NUM_STATE * ndof_total_);

   // Scratch storage for D(k) and L(D(k)).  Two buffers ping-pong.
   Vector D_curr(Q), D_next(Q.Size()), dQ_dxd(Q.Size());
   real_t fac = dt;  // dt^{k+1} / (k+1)!
   I = 0.0;
   I.Add(fac, D_curr);   // I += fac * D(0)

   for (int k = 0; k < order - 1; ++k) {
      // Compute D(k+1) = L(D_curr) = -Σ_d A_d · ∂_{x_d} D_curr.
      D_next = 0.0;
      for (int d = 0; d < 3; ++d) {
         ApplySpatialDerivative(d, D_curr, dQ_dxd);
         const DenseMatrix &A_d = flux_.GetReferenceStarMatrix(d);
         ApplyJacobianPerDOF(A_d, dQ_dxd, D_next, /*sign=*/-1.0);
      }
      // Advance factorial factor: fac *= dt / (k+2)
      fac *= dt / static_cast<real_t>(k + 2);
      I.Add(fac, D_next);   // I += fac * D(k+1)
      std::swap(D_curr, D_next);
   }
}
```

`ApplyJacobianPerDOF(A, X, Y, sign)` is a new helper in the same file
that computes `Y[c, dof] += sign · Σ_{c'} A[c, c'] · X[c', dof]` for
every `dof`. It is the per-DOF application of a 9×9 Jacobian to a
9-component state vector.

**4. Edge cases.**
- `dt = 0`: `I = 0`.
- `order = 2`: no recursion iterations beyond k=0.
- `Q.Size() = 0`: degenerate, return with empty `I`.
- `PML` damping: ADER predictor does NOT include PML. PML is applied
  in the correction step (Phase 4/6). Document this.

### Interfaces
- `ComputeADERTimeIntegrated` public on `WaveOperator`.
- `ApplyJacobianPerDOF` private helper (or free function in the `.inl`).

### Acceptance Criteria
- [ ] `make seas_tpv102_driver` succeeds.
- [ ] All existing wave-operator tests pass.
- [ ] New test `seas_test_ader_ck_predictor`:
  - For a plane P-wave `Q(x, t) = Q0 · sin(k·x − c_p·t)` with
    known `I(x) = ∫_0^dt Q(x, t+τ) dτ = (Q0/(c_p·k)) · (cos(k·x − c_p·(t+dt)) − cos(k·x − c_p·t))`,
    verify `||I_ader − I_analytic||_L2 ≤ C · dt^order` at three
    `dt` values (convergence rate = `order`).
  - For order=2, `dt = 0.01`: relative L2 error ≤ 1e-4.
  - For order=3, `dt = 0.01`: relative L2 error ≤ 1e-6.
- [ ] Runs in < 5 s for a 8×8×4 tet mesh at order 1, O=3.

### Dependencies
- Depends on: Phases 1, 2.
- Required by: Phases 4, 6.

---

## Phase 4: ADER volume integral update

### Goal
After this phase, `WaveOperator::ComputeADERVolumeUpdate(I, rhs)`
computes the volume contribution to the ADER one-step update:
`rhs += ∫_V ∇ · (A · I) dV` — i.e. the integrand evaluated on the
time-integrated state `I` instead of the instantaneous `Q` as in
`ComputeVolumeRHS`.

### Files to Modify
- `dynamic/wave_operator.hpp` — declaration.
- `dynamic/wave_operator.inl` — definition (mirrors `ComputeVolumeRHS`
  but takes `I` instead of `Q`; outputs into `rhs` additively).

### Detailed Requirements

**1. Math.** The ADER corrector's volume term is the STRONG-FORM
DG volume integral applied to the time-integrated state `I`:
```
rhs_e += Σ_d ∫_{e} φ_i · A_d · ∂_{x_d} I  dV   for each element e, each component.
```
In matrix form per element: `rhs_e += Σ_d K_d^e · (A_d · I_e)`
where `K_d^e` is the element stiffness in direction `d` (Phase 1) and
`A_d` is the reference star matrix (Phase 2).

**2. Method signature.**
```cpp
/// ADER volume-integral contribution: rhs += Σ_d K_d · A_d · I.
/// Same additive contract as ComputeVolumeRHS but evaluated on I,
/// not Q.  Does NOT zero rhs first; caller must.
void ComputeADERVolumeUpdate(const Vector &I, Vector &rhs) const;
```

**3. Relation to existing `ComputeVolumeRHS`.** `ComputeVolumeRHS(Q, rhs)`
computes `rhs = Σ_d K_d · A_d · Q` (integrating over every element,
strong form). `ComputeADERVolumeUpdate(I, rhs)` does the SAME operation
on `I` with additive output. **Implementation may share a static helper
`ComputeVolumeStrongForm(const Vector &src, Vector &dst, bool additive)`
to avoid duplication.**

### Edge Cases to Handle
- `I.Size() != NUM_STATE * ndof_total_`: `MFEM_VERIFY`.
- `rhs.Size() == 0`: `rhs.SetSize(NUM_STATE * ndof_total_); rhs = 0;` then add.

### Acceptance Criteria
- [ ] `make seas_tpv102_driver` succeeds.
- [ ] For `I = Q · dt` (equivalent to explicit Euler predictor),
      `ComputeADERVolumeUpdate(I, rhs)` equals `dt · ComputeVolumeRHS(Q, rhs)`
      bit-exact — tested in `seas_test_ader_ck_predictor`.

### Dependencies
- Depends on: Phase 1 (spatial derivative infrastructure).
- Required by: Phase 6.

---

## Phase 5: Time-integrated fault face flux — `FaultFaceFlux::EvaluateADER`

### Goal
After this phase, `FaultFaceFlux::EvaluateADER` computes imposed
time-integrated states `I_imp_+`, `I_imp_-` from time-integrated inputs
`I_+`, `I_-` via a ONE-STEP friction solve on the time-averaged
effective traction. This replaces per-RK4-stage calls to `Evaluate`.
Unit-tested for equivalence with 4 RK4 stages at the `dt → 0` limit.

### Files to Modify
- `dynamic/fault_face_flux.hpp` — add `EvaluateADER` declaration.
- `dynamic/fault_face_flux.cpp` — add definition (~80 LOC).

### Files to Create
- `tests/unit/test_fault_face_flux_ader_equivalence.cpp` — unit test
  comparing `EvaluateADER` output to `Evaluate` output at `dt → 0`.

### Detailed Requirements

**1. Math.** Let `I_± = ∫_0^dt Q_±(τ) dτ`. Then the time-averaged state
is `Q̄_± = I_± / dt`. Trial tractions using Pelties 2012 eq. (7) on
`(Q̄_+, Q̄_-)` give time-averaged `sigma_n_trial, tau1_trial, tau2_trial`.
Solve the friction equation once using these averaged trial values and
the state `psi_n` at `t_n` — the result is `V_abs_avg`, the time-
averaged slip-rate magnitude over `[t_n, t_{n+1}]`. Corrected tractions
are `tau_corr = tau_total - eta_s · V_avg`. Imposed time-integrated
states are then built via Pelties eq. (11)-(12) on the averaged inputs
scaled by `dt`:
`I_imp_± = dt · Q̄_imp_± = dt · [Q̄_± ± Z^{-1} · (σ_corr − Q̄_±[S])]`.
Equivalently: `I_imp_± = I_± ± dt · Z^{-1} · (σ_corr − I_±[S]/dt)`.

**2. Method signature.**
```cpp
/// ADER-flavoured friction solve: evaluates imposed time-integrated
/// states I_imp± from time-integrated inputs I±.
/// @param[in,out] data  Per-DOF friction state.  `psi`, `slip_rate`,
///                       `V1`, `V2`, `tau1_corr`, `tau2_corr`,
///                       `sigma_n_corr` are UPDATED to the time-
///                       AVERAGED values over [t_n, t_n+dt].
/// @param[in]  I_plus   Time-integrated + side state.
/// @param[in]  I_minus  Time-integrated − side state.
/// @param[in]  dt       Time step (for normalisation).
/// @param[out] I_imp_plus  Imposed time-integrated + state.
/// @param[out] I_imp_minus Imposed time-integrated − state.
/// @param[in]  method   Friction solver choice (default Brent).
void EvaluateADER(DOFData &data,
                   const real_t *I_plus, const real_t *I_minus,
                   real_t dt,
                   real_t *I_imp_plus, real_t *I_imp_minus,
                   FrictionSolver::Method method = FrictionSolver::Method::Brent) const;
```

**3. Implementation sketch.**
```cpp
// 1. Compute time-averaged inputs.
real_t Q_avg_plus[NUM_STATE], Q_avg_minus[NUM_STATE];
const real_t inv_dt = 1.0 / dt;
for (int c = 0; c < NUM_STATE; ++c) {
   Q_avg_plus[c]  = I_plus[c]  * inv_dt;
   Q_avg_minus[c] = I_minus[c] * inv_dt;
}

// 2. Call existing Evaluate on averaged inputs.  Evaluate updates
//    data.{sigma_n_corr, tau1_corr, tau2_corr, V1, V2, slip_rate, psi}
//    to VALUES AT THE TIME-AVERAGED STATE — exactly what we want for
//    the ADER single-step DOFData output.
real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
Evaluate(data, Q_avg_plus, Q_avg_minus, Q_imp_plus, Q_imp_minus, method);

// 3. Rescale imposed states back to time-integrated.
for (int c = 0; c < NUM_STATE; ++c) {
   I_imp_plus[c]  = Q_imp_plus[c]  * dt;
   I_imp_minus[c] = Q_imp_minus[c] * dt;
}
```

**4. Semantics of DOFData after call.** After `EvaluateADER`:
- `data.slip_rate, V1, V2` = time-averaged slip rate & components
  over [t_n, t_n+dt].
- `data.sigma_n_corr, tau1_corr, tau2_corr` = time-averaged corrected
  tractions over [t_n, t_n+dt].
- `data.psi` = UNCHANGED (caller must call `UpdateStateAnalytic` with
  the returned `slip_rate`).

**5. Edge cases.**
- `dt ≤ 0`: `MFEM_VERIFY(dt > 0)` — fail loud.
- `I_plus == nullptr` / `I_minus == nullptr`: not caught (trusted
  interface; same contract as `Evaluate`).

### Acceptance Criteria
- [ ] `make seas_test_fault_face_flux_ader_equivalence` succeeds.
- [ ] For a locked fault (`V=0` everywhere), `EvaluateADER` produces
      `data.sigma_n_corr = sigma_n0` within 1 Pa (bit-exact if `I_+ = I_- = 0`).
- [ ] For a symmetric pure strike-slip state, `I_imp_+ - dt·Q̄_+ = dt · Z_s^{-1} · (tau_corr - Q̄_+[SXZ]) · ê_strike`
      matches the hand-derived formula to 1e-10 relative tolerance.
- [ ] In the `dt → 0` limit (smooth `Q_±(τ)` with bounded derivative),
      the difference `|EvaluateADER(data, I_±, dt) − Evaluate(data, Q_±(t_n+dt/2), dt)|`
      scales as `O(dt²)` — verified by 3-point convergence run.

### Dependencies
- Depends on: nothing new (uses existing `Evaluate`).
- Required by: Phase 6.

---

## Phase 6: ADER one-step driver in the wave operator

### Goal
After this phase, `WaveOperator::AdvanceADER(Q, dt, order, Q_new)`
performs one complete ADER step: predictor → volume correction →
face flux (interior + shared + fault) → mass inverse. This is the
wave-operator-level counterpart of a single RK4 step.

### Files to Modify
- `dynamic/wave_operator.hpp` — add `AdvanceADER` declaration.
- `dynamic/wave_operator.inl` — add definition (~150 LOC).
  Re-use existing `ComputeFaceFluxRHS` / `ComputeSharedFaceFluxRHS`
  logic, refactored to take `I` (time-integrated) instead of `Q`
  (instantaneous).

### Detailed Requirements

**1. High-level flow.**
```cpp
void WaveOperator<MeshType>::AdvanceADER(
   const Vector &Q, real_t dt, int order, Vector &Q_new) const
{
   // Predictor: I = ∫_0^dt Q(t+τ) dτ via CK recursion.
   Vector I(Q.Size());
   ComputeADERTimeIntegrated(Q, dt, order, I);

   // Corrector:
   //   Q_new = Q + M^{-1} [ Vol(I) - Face(I) ]
   Vector rhs(Q.Size()); rhs = 0.0;
   ComputeADERVolumeUpdate(I, rhs);       // rhs += ∫ K_d A_d I
   ComputeADERFaceFluxRHS(I, rhs);        // rhs -= ∫_∂V F_h · I
   if constexpr (IsParallelMesh<MeshType>::value) {
      ComputeADERSharedFaceFluxRHS(I, rhs);
   }
   if (pml_layer_) {
      ApplyPMLDamping(I, rhs);            // PML on time-integrated I, scaled by 1
   }

   ApplyMassInverse(rhs);

   Q_new.SetSize(Q.Size());
   add(Q, 1.0, rhs, Q_new);                // Q_new = Q + rhs
}
```

**2. `ComputeADERFaceFluxRHS(I, rhs)`.** Mirrors `ComputeFaceFluxRHS`
but takes `I` (time-integrated) and calls `FaultFaceFlux::EvaluateADER`
for fault faces. For non-fault interior faces and boundary faces, the
Godunov flux is computed on time-integrated states:
`F_h = A_n^+ · I_self + A_n^- · I_nbr` — same formula, different input.
For BC's: Absorbing, FreeSurface unchanged (take `I` instead of `Q`).
For Fault: call `fault_flux_->EvaluateADER(data, I_plus, I_minus, dt, I_imp_plus, I_imp_minus)`
then `flux_.Interior(can_n, I_imp_side, I_imp_side, F_h_side)`.

**3. `ComputeADERSharedFaceFluxRHS`**: mirror of the shared-face path.

**4. Non-duplication strategy.** To minimise copy-paste:
- Refactor `ComputeFaceFluxRHS` to call `ComputeFaceFluxRHSCore(Q, rhs, /*is_ader=*/false, /*dt=*/0.0)`.
- Add `ComputeADERFaceFluxRHS(I, rhs)` as `ComputeFaceFluxRHSCore(I, rhs, /*is_ader=*/true, dt)`.
- The only inner-loop branch is at fault faces: `is_ader ? EvaluateADER : Evaluate`.

**5. Edge cases.**
- `order < 2` or `order > 4`: `MFEM_VERIFY`.
- Zero fault faces on this rank: all ADER-specific branches short-circuit
  as they do under RK4.
- PML: applies to `I`, not `Q`, so the damping contribution is
  `-α · I` (since the damping integrated over `dt` is `α·I`).

### Acceptance Criteria
- [ ] On a linear-elastic plane-wave initial condition (no fault),
      `AdvanceADER(Q, dt, 2, Q_new_ader)` matches `RK4-like 4-stage` to
      `O(dt²)` relative error — tested in a new
      `seas_test_ader_linear_wave_equivalence.cpp`.
- [ ] Existing `seas_test_wave_operator`, `seas_test_parallel_wave_operator`,
      `seas_test_fault_face_flux_*` all pass (RK4 path unchanged).
- [ ] `make test` end-to-end passes.

### Dependencies
- Depends on: Phases 3, 4, 5.
- Required by: Phase 7.

---

## Phase 7: Driver integration — `--time-integrator=ader` CLI flag

### Goal
After this phase, `seas_tpv102_driver --time-integrator=ader --ader-order=2`
runs the production TPV102 setup using the ADER time stepper from
Phase 6. The RK4 path is untouched. All existing command-line options
work unchanged under the default RK4.

### Files to Modify
- `drivers/tpv102_driver.cpp` — add CLI flags, branch on the flag, call
  either the existing RK4 loop or a new ADER loop.

### Detailed Requirements

**1. CLI flags.** Add to `GetArgs`:
- `--time-integrator {rk4|ader}` — default `rk4`.
- `--ader-order N` — integer in {2, 3, 4}, default 2. Used only when
  `--time-integrator=ader`.

**2. Driver loop refactor.** At the top of the time-stepping loop,
branch on the flag:
```cpp
if (time_integrator == "ader") {
   // ADER path.
   for (int step = 0; step < nsteps; step++) {
      real_t dt_step = std::min(dt, tfinal - t);
      if (dt_step <= 0.0) break;

      // Save psi for sub-step state evolution.
      for (int i = 0; i < num_fault_total; i++) {
         psi_n[i] = dof_data[i].psi;
      }

      // Apply nucleation at the stage midpoint (ADER time-averages
      // over [t, t+dt], so evaluate nucleation at t + dt/2 for
      // 2nd-order time accuracy; for O > 2, see R-007 below).
      ApplyNucleation(dof_data, num_fault_total, fault_coords,
                      t + dt_step / 2.0);

      // Advance Q via ADER one-step predictor-corrector.
      Vector Q_new(Q.Size());
      wave.AdvanceADER(Q, dt_step, ader_order, Q_new);
      Q = Q_new;

      t += dt_step;

      // ADER returns time-averaged DOFData.V1/V2/slip_rate/tau*_corr/
      // sigma_n_corr directly — no RK4-averaging block needed.
      // Update psi using the time-averaged slip rate.
      for (int i = 0; i < num_fault_total; i++) {
         dof_data[i].psi = UpdateStateAnalytic(
            psi_n[i], dof_data[i].slip_rate, dof_data[i].Dc,
            dt_step, TPV102Params::f0, TPV102Params::b, TPV102Params::V0);
         dof_data[i].slip1 += dof_data[i].V1 * dt_step;
         dof_data[i].slip2 += dof_data[i].V2 * dt_step;
      }

      // V_max tracking + paraview_write as in the RK4 loop.
      // [existing code follows]
   }
} else {
   // RK4 path — UNCHANGED from current.
}
```

**3. `_k4` diagnostic compatibility.** Under ADER, there is no stage-4
vs averaged distinction — `dof_data` carries the time-averaged values
directly. Behaviour: under ADER, populate `pv_local_*_k4 = pv_local_*`
so the VTU `<field>_k4 - <field>` delta is ZERO by construction. This
documents in the output that ADER eliminates H-V92-K by formulation.

**4. Nucleation timing for higher order.** At `O=2`, evaluating
nucleation once at `t+dt/2` is 2nd-order-accurate. For `O ≥ 3`, the
nucleation perturbation should be time-integrated across the step
like the bulk state. Scope limitation: at this phase we only support
`O=2` and `O=3` in production (nucleation at midpoint), with the
caveat documented. `O=4` is available for smooth test problems with
no nucleation.

**5. Logging.** Print at run start:
```
Time integrator: ADER-O(N)  (was: RK4)
```
where `N` is `ader_order`.

### Acceptance Criteria
- [ ] `seas_tpv102_driver --time-integrator=rk4 …` produces
      byte-identical station output to the pre-change baseline (RK4
      path must not have changed).
- [ ] `seas_tpv102_driver --time-integrator=ader --ader-order=2 …`
      runs to completion on the TPV102 1000 m laptop mesh in 5 min,
      and produces a valid VTU with `normal_stress_k4 == normal_stress`
      (delta zero by construction under ADER — this confirms the
      H-V92-K elimination claim at the output level).
- [ ] Existing `make test-bp5-*` unaffected (BP5 does not use this
      driver).

### Dependencies
- Depends on: Phase 6.
- Required by: Phase 8.

---

## Phase 8: Tests — convergence, equivalence, and TPV102 smoke

### Goal
After this phase, the ADER path has five passing tests:
1. CK predictor convergence (Phase 3).
2. Volume integral consistency (Phase 4).
3. Fault face flux equivalence with RK4 at dt → 0 (Phase 5).
4. Wave operator one-step equivalence with RK4 on a linear plane wave
   (Phase 6).
5. TPV102 smoke — `--time-integrator=ader --ader-order=2` runs the
   existing TPV102 1-element fixture to t=0.1s without diverging.

### Files to Create
- `tests/unit/test_ader_tpv102_smoke.cpp` — 1-element TPV102 fixture
  initial condition, run 10 ADER steps at dt=1e-3, assert:
  - `|dof_data[0].sigma_n_corr - 120e6| < 1e3` (no σ_n drift, since no
    rupture; ADER should preserve equilibrium to machine precision).
  - `|dof_data[0].V1| < 1e-12` (no spurious dip slip rate).
  - `|dof_data[0].V2 - V_ini| < 1e-15` (strike slip rate at initial
    equilibrium).

### Files to Modify
- `Makefile` — register the 5 new tests and add an aggregate target
  `make test-ader` that runs them in order.

### Detailed Requirements
All tests follow the existing pattern:
- `TEST_ASSERT` / `TEST_NEAR` macros.
- `main()` returns non-zero iff any test fails.
- Output directory under `miniapps/seas/test_ader_*` or per-test.

### Acceptance Criteria
- [ ] `make test-ader` — all 5 tests PASS.
- [ ] Total runtime < 2 min on laptop.
- [ ] `make test` (existing unit-test aggregate) — no regression.

### Dependencies
- Depends on: Phases 3, 4, 5, 6, 7.
- Required by: nothing (final phase).

---

## Testing Strategy

### Per-phase acceptance
Each phase has its own test that must pass before the next phase
starts. Per CLAUDE.md's "Files Requiring Extreme Care" policy, after
Phases 6 and 7 we run `make test` end-to-end (all BP5 and TPV102
tests) to catch any cross-contamination of the `wave_operator.inl` or
`fault_face_flux.cpp` files that the RK4 path uses.

### Convergence validation (Phase 3)
Use a 3D plane-wave analytical solution:
```
Q(x, t) = Q0 · sin(k·x − c_p·t · n̂)
```
where `n̂` is a unit propagation direction, `c_p = √((λ+2μ)/ρ)`, `Q0`
chosen to satisfy `A · Q0 · k̂ = c_p · Q0` (right eigenvector of `A_n`).
Evaluate `I_analytic(x) = ∫_0^dt Q(x, t_n+τ) dτ = (Q0 / c_p·|k|) · [cos(k·x − c_p·t_n) − cos(k·x − c_p·(t_n+dt))]`.
Compute `I_ader = ComputeADERTimeIntegrated(Q(t_n), dt, order, I_ader)`.
`||I_ader − I_analytic||_L2 ≤ C · dt^order` for order ∈ {2, 3, 4},
measured at `dt = {0.1, 0.05, 0.025}` of CFL.

### Equivalence validation (Phase 5)
`EvaluateADER(data, I_+, I_-, dt)` vs `Evaluate(data, Q_midpt_+, Q_midpt_-)`
at `I_± = dt · Q_midpt_±`: results equal to 1e-10 relative tolerance
for a fault QP in locked / early-rupture / saturated regimes.

### TPV102 reference comparison (Phase 7, manual — not a gated test)
After Phase 7 passes its unit tests, a local 1000 m laptop run under
ADER should produce:
- Station `flt_0_7.5` σ_n bounded within 120 MPa ± 1 MPa for t ∈ [0, 2s].
  (If H-V92-K is the actual cause, ADER will eliminate the σ_n drift.)
- VTU `normal_stress_k4 - normal_stress == 0` per cell, confirming the
  one-step property.

Frontera 12 s run is deferred to a separate sbatch plan, user-gated.

## Risk Assessment

### High risk
- **Stiffness matrix assembly in Phase 1.** The existing
  `ComputeVolumeRHS` builds `K_d^e` inline. Factoring it out could
  alter the strong-form vs weak-form sign or quadrature order
  subtly. Mitigation: the unit test in Phase 1 compares
  `ApplySpatialDerivative` against `DiffusionIntegrator::AssembleElementMatrix`
  to a P2 polynomial test function.

- **CK recursion numerical stability for O ≥ 3.** Higher-order
  derivatives amplify mesh-aligned noise. Mitigation: test Phase 3 at
  O=3 AND O=4 on the plane-wave benchmark; if O=4 fails, document as
  a known limit and ship with O=2/O=3 only.

- **Fault-face ADER EvaluateADER correctness.** The friction solve
  is nonlinear; averaging the input then solving is not bit-identical
  to solving then averaging, but the two agree to `O(dt²)`. Mitigation:
  Phase 5 unit test verifies this convergence rate explicitly.

### Medium risk
- **Driver refactor in Phase 7 introducing a regression in the RK4
  path.** Mitigation: before any Phase 7 change, snapshot the RK4
  station output on a small fixture (1000 m, 2 s, 1 rank); after the
  refactor, byte-compare. Any difference is a regression.

- **The `_k4` diagnostic fields were added under plan v9.2.0 §19 and
  are expected to be zero-delta under ADER.** If a bug in ADER
  produces nonzero deltas, the user may mistakenly conclude H-V92-K
  still exists. Mitigation: document the ADER-zero-delta property
  prominently in Phase 7 and in the test assertions.

### Low risk
- **Memory footprint for O = 4.** One extra `D_next` buffer of size
  `NUM_STATE · ndof_total ≈ 9 · 10⁷` doubles ≈ 700 MB per rank at
  production. Mitigation: reuse scratch buffers ping-pong; peak two
  copies.

- **BP5 unaffected.** BP5 uses a different driver (seas_bp5_full,
  quasi-dynamic SIPG), not this wave operator's ADER path. No risk of
  BP5 regression from this plan.

## Scope and timeline

| Phase | LOC (new) | LOC (modified) | Est. wall-clock |
|---|---:|---:|---:|
| 1 Spatial derivative | 120 | 60 | 1 day |
| 2 Star matrices | 20 | 15 | 0.5 day |
| 3 CK predictor | 180 | 30 | 1.5 day |
| 4 Volume update | 80 | 30 | 0.5 day |
| 5 Fault ADER flux | 100 | 10 | 0.5 day |
| 6 ADER advance | 250 | 150 | 2 day |
| 7 Driver flag | 150 | 80 | 1 day |
| 8 Tests + Makefile | 400 | 40 | 1 day |
| **Total** | **~1300** | **~415** | **~8 days** |

## Deferred / out of scope

1. **Heterogeneous material star matrices** — TPV102 is homogeneous,
   so we use a single global star matrix per direction. Bent elements
   and heterogeneous materials require per-element star matrices;
   punted to a v2 plan if/when BP5 or a future benchmark needs it.
2. **ADER local time stepping (LTS)** — SeisSol's main production
   feature, enabling per-element Δt. Massive addition; out of scope
   for this plan. Use global Δt = min CFL across all elements.
3. **Higher-order nucleation time integration** — the plan evaluates
   nucleation at `t + dt/2` (O(dt²) accurate). For O ≥ 3, nucleation
   should be Taylor-expanded to match. Deferred.
4. **GPU kernels** — SeisSol's ADER has GPU/batched variants. SEAS-MFEM
   is CPU-only; no GPU work needed now.
5. **ADER-specific CFL tuning** — use the same CFL as RK4 (both are
   CK-class explicit). A dedicated CFL study is a follow-up.
6. **BP5 quasi-dynamic ADER path** — BP5 uses implicit DG + adaptive
   RK45; ADER is an explicit scheme for the dynamic path. BP5 does
   not benefit from this plan. Out of scope.

## Quick-reference for implementation agent

- Read CLAUDE.md "Files Requiring Extreme Care" — `wave_operator.inl`,
  `fault_face_flux.cpp` are on it.
- Follow v9.1.0 commit pattern: one phase per commit; tests added in
  the same commit as the corresponding production code.
- Keep the RK4 path byte-identical. Before modifying any existing
  function, snapshot the output of the pre-change driver on a tiny
  fixture and re-compare after each phase.
- The plan allows phases to land individually: Phase 1 + 2 + 3 can
  ship and be reviewed in isolation (no user-visible behaviour
  change); Phases 4-6 land together (full ADER path alive but not
  invokable); Phase 7 wires the CLI; Phase 8 gates the tests.
