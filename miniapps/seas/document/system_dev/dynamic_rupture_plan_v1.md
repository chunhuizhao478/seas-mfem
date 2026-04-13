# Dynamic Rupture Implementation Plan v1

**Date:** 2026-04-11
**Status:** Draft — for review
**Benchmark:** SCEC TPV102 (half-space, aging law rate-and-state friction)
**References:**
- SCEC_validation_ageing_law.pdf (TPV101/102 benchmark specification)
- Dumbser & Käser (2006), ADER-DG for 3D elastic waves
- de la Puente et al. (2009), fault Riemann solver for dynamic rupture
- Pelties et al. (2014), SeisSol fault coupling
- SeisSol source code `/Users/chunhuizhao/projects/SeisSol`

---

## Overview

Build a 3D dynamic rupture DG code within the SEAS-MFEM refactored system, verified against SCEC benchmark TPV102. The code solves the velocity-stress elastic wave equation with an explicit DG scheme, coupling seismic wave propagation to rate-and-state friction on an embedded fault interface via a fault Riemann solver. Three objectives: (1) standalone dynamic rupture verified with TPV102, (2) QD↔FD transfer architecture for future hybrid simulation, (3) GPU acceleration path using MFEM device infrastructure.

---

## TPV102 Benchmark Summary

### Problem Setup

A planar fault at z=0 in a half-space (y > 0), with a free surface at y=0.

| Parameter | Value |
|-----------|-------|
| ρ (density) | 2670 kg/m³ |
| cₛ (S-wave speed) | 3464 m/s |
| cₚ (P-wave speed) | 6000 m/s |
| μ (shear modulus) | ρ·cₛ² = 32.04 GPa |
| λ (first Lamé) | ρ·cₚ² − 2μ = 32.04 GPa |
| ν (Poisson's ratio) | 0.25 |

### Friction Law: Regularized R&S Aging

$$\tau = a \sigma_n \text{arcsinh}\left[\frac{V}{2V_0} \exp\left(\frac{f_0 + b \ln(V_0\theta/L)}{a}\right)\right]$$

$$\frac{d\theta}{dt} = 1 - \frac{V\theta}{L}$$

| Parameter | Value |
|-----------|-------|
| f₀ | 0.6 |
| V₀ | 10⁻⁶ m/s |
| b | 0.012 |
| L (slip distance) | 0.02 m |
| a (VW interior) | 0.008 (a−b = −0.004, velocity-weakening) |
| a (VS exterior) | 0.016 (a−b = +0.004, velocity-strengthening) |
| Δa₀ (VS overshoot) | 0.008 |
| σₙ (normal stress) | 120 MPa (compressive) |
| τᵢₙᵢ (initial shear) | 75 MPa |
| Vᵢₙᵢ | 10⁻¹² m/s |

### Fault Geometry

- VW rectangle: |x| < W, 0 < y < W, where W = 15 km (30 km × 15 km)
- Transition layer width: w = 3 km (C∞ smooth taper via boxcar function B)
- Hypocenter: (x₀, y₀) = (0, 7.5 km)

### Nucleation

Stress perturbation: Δτ(x,y,t) = Δτ₀ · F(r) · G(t)

| Parameter | Value |
|-----------|-------|
| Δτ₀ | 25 MPa |
| R (radius) | 3 km |
| T (rise time) | 1 s |
| F(r) | exp(r²/(r²−R²)) for r<R, 0 otherwise |
| G(t) | exp((t−T)²/(t(t−2T))) for 0<t<T, 1 for t≥T |

### Initial Conditions

- σ_zx = τᵢₙᵢ = 75 MPa everywhere
- σ_zz = −σₙ = −120 MPa everywhere
- vₓ = ±Vᵢₙᵢ/2 for z ≷ 0
- θᵢₙᵢ(x,y) computed from friction law self-consistency (Eq. 6 of spec)

### Output Requirements

- Simulation: t = 0 to 12 s
- 9 fault stations: time histories of slip, slip rate, traction, log₁₀(θ)
- 6 free-surface stations (TPV102 only): displacement and velocity
- Rupture front contours: time when V first exceeds 1 mm/s

---

## Constraints

- **Interface constraints:** `ConstitutiveModel`, `FaultBasis`, `FaultGeometry`, `BoundaryConfig`, `SEASConfig` interfaces are unchanged. New dynamic rupture code extends but does not modify them.
- **Dependency constraints:** Uses MFEM `L2_FECollection` for DG space, `ODESolver` for explicit RK, `ParMesh` for parallel. No new external dependencies.
- **Convention constraints:** Follow existing SEAS-MFEM module organization (one directory per concern: `dynamic/`, files under `miniapps/seas/`). Follow existing naming patterns.
- **Numerical constraints:** CFL-limited explicit time stepping. Conservation of energy at fault interface. Friction solver convergence to |g| < 10⁻⁸.

---

## Governing Equations

### Velocity-Stress System (9 unknowns)

The 3D isotropic elastic wave equation in first-order form:

**State vector:** Q = (σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz, v_x, v_y, v_z)ᵀ

**System:**
```
∂Q/∂t + A ∂Q/∂x + B ∂Q/∂y + C ∂Q/∂z = 0
```

where A, B, C are 9×9 Jacobian matrices. For direction x:

```
A = | 0  0  0  0  0  0  -(λ+2μ)  -λ      -λ     |
    | 0  0  0  0  0  0  -λ       -(λ+2μ)  -λ     |
    | 0  0  0  0  0  0  -λ       -λ       -(λ+2μ) |
    | 0  0  0  0  0  0  0        -μ       0       |
    | 0  0  0  0  0  0  0        0        -μ      |
    | 0  0  0  0  0  0  -μ       0        0       |
    | -1/ρ  0  0  0  0  0  0  0  0                 |
    | 0  0  0  -1/ρ  0  0  0  0  0                 |
    | 0  0  0  0  -1/ρ  0  0  0  0                 |
```

(B and C follow by cyclic permutation of coordinates.)

**Eigenvalues:** ±cₚ, ±cₛ (×2), 0 (×3) — independent of propagation direction for isotropic media.

### DG Weak Form

In each element Tₘ with polynomial basis Φₗ(ξ,η,ζ) of degree N:

```
Q_h^(m) = Σ_l Q̂_l^(m)(t) · Φ_l(ξ,η,ζ)
```

Weak form (after integration by parts):

```
∫_T Φ_k ∂Q/∂t dV + ∫_∂T Φ_k F^h dS − ∫_T (∂Φ_k/∂x_j) F_j(Q) dV = 0
```

where F^h is the numerical flux at element interfaces.

### Godunov Upwind Flux at Interior Faces

Transform to face-normal coordinate system via rotation T:

```
F^h = ½ T(A + |A|)T⁻¹ Q_interior + ½ T(A − |A|)T⁻¹ Q_exterior
```

where |A| = R·|Λ|·R⁻¹ with |Λ| = diag(cₚ, cₛ, cₛ, 0, 0, 0, cₛ, cₛ, cₚ).

Physical meaning: outgoing waves use interior state, incoming waves use exterior state.

### Absorbing BC

```
F^abs = ½ T(A + |A|)T⁻¹ Q_interior
```

Only outgoing waves — incoming set to zero.

### Free Surface BC

```
F^free = ½ T(A + |A|)T⁻¹ Q + ½ T(A − |A|)·Γ·T⁻¹ Q
```

where Γ = diag(−1, 1, 1, −1, 1, −1, 1, 1, 1) mirrors stress components.

### Fault Riemann Solver (from SeisSol, Uphoff eq. 4.51)

**Impedances:**
```
Zₚ = ρ · cₚ,    Zₛ = ρ · cₛ
ηₚ = Zₚ⁺·Zₚ⁻ / (Zₚ⁺ + Zₚ⁻)    (harmonic mean)
ηₛ = Zₛ⁺·Zₛ⁻ / (Zₛ⁺ + Zₛ⁻)
```

**Step 1: Trial traction** (locked-fault Godunov state):
```
σₙ_trial = ηₚ · (v_n⁻ − v_n⁺ + σₙ⁺/Zₚ⁺ + σₙ⁻/Zₚ⁻)
τ₁_trial = ηₛ · (v_t1⁻ − v_t1⁺ + τ₁⁺/Zₛ⁺ + τ₁⁻/Zₛ⁻)
τ₂_trial = ηₛ · (v_t2⁻ − v_t2⁺ + τ₂⁺/Zₛ⁺ + τ₂⁻/Zₛ⁻)
```

**Step 2: Friction solve** (Newton-Raphson for slip rate magnitude V̂):
```
g(V̂) = −(1/ηₛ) · (|σₙ| · μ(V̂, θ) − Θ) − V̂ = 0

where Θ = |τ_total| = √((τ₁₀ + τ₁_trial)² + (τ₂₀ + τ₂_trial)²)
```

**Step 3: Slip rate decomposition** into tangential components:
```
strength = |σₙ| · μ(V̂, θ)
V₁ = V̂ · (τ₁₀ + τ₁_trial) / (strength + ηₛ · V̂)
V₂ = V̂ · (τ₂₀ + τ₂_trial) / (strength + ηₛ · V̂)
```

**Step 4: Corrected traction:**
```
τ₁_corrected = τ₁_trial − ηₛ · V₁
τ₂_corrected = τ₂_trial − ηₛ · V₂
```

**Step 5: Imposed state** (Uphoff eq. 4.60):
```
Plus side:  v⁺_imposed = v⁺ + (1/Zₛ⁺) · (τ_corrected − τ⁺)
Minus side: v⁻_imposed = v⁻ − (1/Zₛ⁻) · (τ_corrected − τ⁻)
Both sides: σ_imposed = τ_corrected   (traction continuity)
```

### CFL Condition

```
Δt ≤ CFL · h_min / cₚ,    CFL ≈ 1/(2N+1)
```

where h_min is the minimum element diameter and N is the polynomial order.

---

## Phase 1: Velocity-Stress Wave Operator (No Fault)

### Goal
A working explicit DG solver for 3D elastic wave propagation in a homogeneous medium, with Godunov flux at interior faces, verified by plane wave convergence test.

### Files to Create

- `dynamic/wave_state.hpp` — Q = (σ, v) state vector container, 9-component layout, conversion utilities
- `dynamic/wave_flux.hpp` — Jacobian matrices A, B, C; rotation matrix T; Godunov upwind flux computation
- `dynamic/wave_operator.hpp` — `WaveOperator` class: DG assembly for velocity-stress system, implements `TimeDependentOperator::Mult()`
- `dynamic/wave_operator.cpp` — Implementation of volume + face integrals
- `tests/unit/test_wave_flux.cpp` — Unit tests for Godunov flux correctness
- `tests/unit/test_wave_operator.cpp` — Plane wave convergence test

### Detailed Requirements

1. **`WaveState` (9-component Q container)**

   ```cpp
   // dynamic/wave_state.hpp
   namespace mfem::seas::dynamic {

   /// Index map for velocity-stress state vector Q[9]
   enum QIndex : int {
      SXX = 0, SYY = 1, SZZ = 2,  // Normal stresses
      SXY = 3, SYZ = 4, SXZ = 5,  // Shear stresses
      VX = 6, VY = 7, VZ = 8      // Velocities
   };
   static constexpr int NUM_STATE = 9;

   } // namespace
   ```

2. **`GodunovFlux` — Upwind flux computation**

   ```cpp
   // dynamic/wave_flux.hpp
   class GodunovFlux {
   public:
      GodunovFlux(real_t lambda, real_t mu, real_t rho);

      // Compute flux at an interior face given states on both sides
      // normal: outward unit normal from element_self
      // Q_self[9], Q_nbr[9]: extrapolated boundary states
      // flux[9]: output Godunov flux vector
      void ComputeInteriorFlux(const real_t normal[3],
                                const real_t Q_self[9],
                                const real_t Q_nbr[9],
                                real_t flux[9]) const;

      // Compute absorbing BC flux (only outgoing waves)
      void ComputeAbsorbingFlux(const real_t normal[3],
                                 const real_t Q_self[9],
                                 real_t flux[9]) const;

      // Compute free-surface BC flux (stress reflection)
      void ComputeFreeSurfaceFlux(const real_t normal[3],
                                   const real_t Q_self[9],
                                   real_t flux[9]) const;

      real_t GetMaxWaveSpeed() const { return cp_; }

   private:
      real_t lambda_, mu_, rho_;
      real_t cp_, cs_;  // P-wave and S-wave speeds

      // Build 9x9 rotation matrix T and its inverse for given normal
      void BuildRotation(const real_t normal[3],
                          DenseMatrix &T, DenseMatrix &Tinv) const;

      // Compute (A + |A|)/2 and (A - |A|)/2 for x-direction Jacobian
      void BuildSplitFluxMatrices(DenseMatrix &Aplus,
                                   DenseMatrix &Aminus) const;
   };
   ```

   Implementation: For a given face normal n, rotate Q to face-aligned coordinates via T, apply A⁺ = (A+|A|)/2 to Q_self and A⁻ = (A−|A|)/2 to Q_nbr in the rotated frame, rotate back.

   For homogeneous material, the split flux matrices A⁺, A⁻ are constant and can be precomputed once globally. For heterogeneous material across faces, they must be computed per face using both sides' properties (Riemann problem with material discontinuity).

3. **`WaveOperator` — DG wave equation RHS**

   ```cpp
   // dynamic/wave_operator.hpp
   class WaveOperator : public TimeDependentOperator {
   public:
      WaveOperator(ParMesh &mesh, int order,
                   const ConstitutiveModel &model,
                   real_t rho,
                   const BoundaryConfig &bdr_config);

      // Compute RHS: y = M⁻¹ · (volume_integral + face_flux_integral)
      // Called by explicit ODE solver
      void Mult(const Vector &Q, Vector &dQdt) const override;

      // CFL-limited time step
      real_t GetMaxDt() const;

      // Access FE space
      ParFiniteElementSpace &GetFESpace() { return *fes_; }
      const ParFiniteElementSpace &GetFESpace() const { return *fes_; }

   private:
      ParMesh &mesh_;
      int order_;
      real_t lambda_, mu_, rho_, cp_, cs_;
      std::unique_ptr<L2_FECollection> fec_;
      std::unique_ptr<ParFiniteElementSpace> fes_;

      // Precomputed element mass inverse (block-diagonal for DG)
      // Stored as dense matrix per element
      std::vector<DenseMatrix> elem_mass_inv_;

      // Precomputed Godunov flux
      GodunovFlux flux_;

      // BC classification
      BoundaryConfig bdr_config_;
      enum class FaceBC { Interior, Absorbing, FreeSurface, Fault };
      std::vector<FaceBC> face_bc_type_;

      // Assembly
      void AssembleElementMassInverse();
      void ComputeVolumeIntegral(const Vector &Q, Vector &rhs) const;
      void ComputeFaceFluxIntegral(const Vector &Q, Vector &rhs) const;
   };
   ```

   **FE space:** `L2_FECollection(order, 3)` with `vdim = 9` (9-component vector DG space). Each DOF stores all 9 components of Q.

   **Mass matrix:** Block-diagonal in DG — one dense matrix per element. Precomputed and inverted at setup. For orthogonal basis functions (Legendre), M is diagonal.

   **Volume integral:** For each element, compute ∫_T (∂Φ_k/∂x_j) · F_j(Q) dV using quadrature. F_j(Q) is the physical flux in direction j, computed from the Jacobian matrix A, B, or C applied to Q.

   **Face flux integral:** For each face (interior, boundary, shared), evaluate the Godunov flux F^h and accumulate ∫_F Φ_k · F^h dS for both adjacent elements. For shared faces in parallel, need ghost data exchange (MFEM's `ParFiniteElementSpace` handles this via `ExchangeFaceNbrData()`).

4. **Time stepping:** Use MFEM's `RK4Solver` (4th-order explicit Runge-Kutta) initially. CFL condition:
   ```
   Δt = CFL_safety · h_min / cₚ / (2·order + 1)
   ```
   with `CFL_safety = 0.5`.

### Edge Cases
- Degenerate elements (h_min → 0): CFL guard `MFEM_VERIFY(h_min > 0)`
- Machine-zero states: flux computation must handle Q = 0 without division by zero
- Parallel shared faces: ensure ghost exchange before face flux assembly

### Acceptance Criteria
- [ ] Plane wave propagation in homogeneous cube: P-wave and S-wave travel at correct speeds (error < 1% at 10 elements/wavelength)
- [ ] Convergence test: L2 error decreases at rate O(h^{N+1}) for polynomial order N = 1, 2, 3
- [ ] Energy conservation: total elastic energy constant (to machine precision) in closed domain with reflecting walls
- [ ] `make test` still passes (all existing QD tests)

### Dependencies
- Depends on: existing `ConstitutiveModel` (Phase 3), `BoundaryConfig` (Phase 4), MFEM `L2_FECollection`
- Required by: Phase 2, 3

---

## Phase 2: Boundary Conditions — Absorbing + Free Surface

### Goal
Wave operator handles absorbing and free-surface boundaries correctly, verified by reflection/absorption tests.

### Files to Create
- `tests/unit/test_wave_bc.cpp` — BC verification tests

### Files to Modify
- `dynamic/wave_operator.cpp` — Add BC-type dispatch in face flux assembly
- `dynamic/wave_flux.hpp` — Already has `ComputeAbsorbingFlux` and `ComputeFreeSurfaceFlux` from Phase 1

### Detailed Requirements

1. **BC classification at setup:** Parse `BoundaryConfig` to classify each boundary face:
   - `natural_attrs` → FreeSurface (traction-free)
   - `dirichlet_attrs` → Absorbing (non-reflecting)
   - `fault_attr` → Fault (Phase 3)

   For TPV102: free surface at y=0, absorbing on all other exterior boundaries.

2. **Absorbing BC implementation:** At absorbing faces, the flux is:
   ```
   F^abs = ½ T(A + |A|)T⁻¹ Q
   ```
   This is the outgoing-wave-only part of the Godunov flux. No neighbor state needed.

3. **Free-surface BC implementation:** At free-surface faces:
   ```
   F^free = ½ T(A + |A|)T⁻¹ Q + ½ T(A − |A|) · Γ · T⁻¹ Q
   ```
   where Γ = diag(−1, 1, 1, −1, 1, −1, 1, 1, 1) mirrors the face-normal and shear stress components.

   Physical meaning: virtual mirror element has same velocity but opposite face-normal/shear stresses → traction-free condition.

### Acceptance Criteria
- [ ] Free surface: plane P-wave hitting free surface reflects with correct amplitude and polarity (P→P and P→SV conversion)
- [ ] Absorbing BC: plane wave exits domain with < 1% reflected energy at normal incidence
- [ ] LOH.1 benchmark (optional): layer-over-halfspace with point source matches reference solution at surface stations

### Dependencies
- Depends on: Phase 1 (wave operator)
- Required by: Phase 4 (TPV102 needs free surface + absorbing)

---

## Phase 3: Fault Riemann Solver + R&S Aging Law Coupling

### Goal
Dynamic rupture on a planar fault with rate-and-state friction, computing slip rate from the Riemann problem and feeding corrected traction back into the wave operator.

### Files to Create
- `dynamic/fault_riemann.hpp` — Fault Riemann solver: trial traction → friction solve → imposed state
- `dynamic/fault_riemann.cpp` — Implementation
- `dynamic/aging_law_dr.hpp` — Aging law state update for dynamic rupture (analytic time integration)
- `tests/unit/test_fault_riemann.cpp` — Riemann solver unit tests

### Files to Modify
- `dynamic/wave_operator.hpp` — Add fault face handling, integrate with `FaultRiemann`
- `dynamic/wave_operator.cpp` — At fault faces, call `FaultRiemann` instead of standard Godunov flux

### Detailed Requirements

1. **`FaultRiemann` — Fault-face Riemann solver**

   ```cpp
   // dynamic/fault_riemann.hpp
   class FaultRiemann {
   public:
      struct FaultFaceData {
         real_t etaP, etaS, invEtaS;    // Impedance parameters
         real_t invZp_plus, invZp_minus;  // Inverse P-impedances
         real_t invZs_plus, invZs_minus;  // Inverse S-impedances

         // Initial stress (constant per DOF, set at t=0)
         real_t sigma_n0;    // Initial normal stress (positive = compression)
         real_t tau1_0;      // Initial shear traction component 1
         real_t tau2_0;      // Initial shear traction component 2

         // Friction parameters (spatially varying)
         real_t a, b, L, f0, V0;

         // State variable
         real_t theta;       // Current state variable value
         real_t slip;        // Accumulated slip magnitude
         real_t slip_rate;   // Current slip rate magnitude
      };

      struct FaultFaceResult {
         real_t sigma_n_corrected;
         real_t tau1_corrected, tau2_corrected;
         real_t V1, V2;              // Slip rate components
         real_t V_magnitude;         // Slip rate magnitude
         real_t theta_new;           // Updated state variable
      };

      FaultRiemann(real_t lambda, real_t mu, real_t rho);

      /// Solve the fault Riemann problem at one quadrature point.
      ///
      /// @param data   Per-DOF fault data (impedances, initial stress, friction params, state)
      /// @param Q_plus  State on plus side [9] (element owning face)
      /// @param Q_minus State on minus side [9] (neighbor element)
      /// @param dt      Current time step (for state variable update)
      /// @param result  Output: corrected traction, slip rate, new state
      void Solve(const FaultFaceData &data,
                 const real_t Q_plus[9],
                 const real_t Q_minus[9],
                 real_t dt,
                 FaultFaceResult &result) const;

      /// Compute imposed state for plus/minus sides from friction result.
      /// Used to set the Godunov flux at the fault face.
      void ComputeImposedState(const FaultFaceData &data,
                                const real_t Q_plus[9],
                                const real_t Q_minus[9],
                                const FaultFaceResult &result,
                                real_t imposed_plus[9],
                                real_t imposed_minus[9]) const;

   private:
      real_t lambda_, mu_, rho_;

      /// Newton-Raphson solve for slip rate magnitude.
      /// g(V̂) = −(1/ηₛ)(|σₙ|·μ(V̂,θ) − Θ) − V̂ = 0
      /// Returns slip rate magnitude.
      real_t SolveSlipRate(real_t sigma_n, real_t abs_traction,
                           real_t etaS, real_t invEtaS,
                           real_t a, real_t b, real_t L,
                           real_t f0, real_t V0,
                           real_t theta,
                           real_t V_guess) const;

      /// Regularized R&S friction coefficient.
      /// μ = a · arcsinh[V/(2V₀) · exp((f₀ + b·ln(V₀θ/L))/a)]
      static real_t FrictionCoeff(real_t V, real_t theta,
                                   real_t a, real_t b, real_t L,
                                   real_t f0, real_t V0);

      /// Derivative dμ/dV for Newton-Raphson.
      static real_t FrictionCoeffDeriv(real_t V, real_t theta,
                                        real_t a, real_t b, real_t L,
                                        real_t f0, real_t V0);
   };
   ```

2. **Newton-Raphson solver for slip rate** (mirroring SeisSol `invertSlipRateIterative`):

   Find V̂ such that g(V̂) = 0 where:
   ```
   g(V̂) = −(1/ηₛ)(|σₙ| · μ(V̂, θ) − Θ) − V̂
   g'(V̂) = −(1/ηₛ) · |σₙ| · dμ/dV − 1
   ```

   Newton update: V̂_{k+1} = max(ε, V̂_k − g_k/g'_k)
   Convergence: |g| < 10⁻⁸, max 60 iterations.

   The `arcsinh` evaluation must be numerically stable for very large arguments (when V·exp(c) overflows). Use: arcsinh(x) = ln(x + √(x²+1)) ≈ ln(2x) for large x.

3. **State variable update** — aging law analytic solution (mirroring SeisSol `AgingLaw.h`):
   ```cpp
   // dynamic/aging_law_dr.hpp
   inline real_t UpdateAgingLaw(real_t theta_old, real_t V, real_t L, real_t dt)
   {
      real_t x = -V * dt / L;
      real_t expx = std::exp(x);
      real_t expm1x = -std::expm1(x);  // = 1 - exp(x), accurate for small x
      return theta_old * expx + (L / V) * expm1x;
   }
   ```

   This is the analytic solution of dθ/dt = 1 − Vθ/L assuming V constant over dt.

4. **Integration with WaveOperator:** At fault faces during face flux assembly:
   - Call `FaultRiemann::Solve()` for each fault face DOF
   - Call `FaultRiemann::ComputeImposedState()` to get the imposed velocity/stress
   - Use imposed state as the "Godunov state" at the interface for both elements
   - The flux contribution to element m is: ∫_F Φ_k · (F_imposed − F_self) dS

5. **State variable iteration** (Kaneko et al. 2008): Per time step:
   1. Compute trial traction from current Q states
   2. First pass: solve for V with current θ, update θ analytically
   3. Second pass: solve for V with updated θ, average V
   This improves accuracy of the state variable update.

### Edge Cases
- Very small slip rate (V → 0): `arcsinh(V·C)` where C can be huge. Use stable `arsinhexp` implementation.
- Tensile normal stress (σₙ > 0): Clamp to zero friction: `σₙ_eff = min(σₙ, 0)`.
- Newton divergence: If not converged after 60 iterations, fall back to bisection on [0, Θ/ηₛ].

### Acceptance Criteria
- [ ] Locked fault (V=0 imposed): trial traction equals Godunov state, no slip
- [ ] Frictionless fault (μ=0): full stress drop, slip rate = Θ/ηₛ
- [ ] Known 1D solution: SH wave impinging on friction-controlled fault matches analytical Riemann solution
- [ ] Newton convergence: < 10 iterations for typical TPV102 parameters
- [ ] State variable: aging law θ(t) matches analytical solution for constant V

### Dependencies
- Depends on: Phase 1 (wave operator), existing `DieterichRuinaFriction` for reference
- Required by: Phase 4 (TPV102)

---

## Phase 4: TPV102 Full Benchmark Setup

### Goal
Complete TPV102 simulation producing SCEC-format output at all 9+6 stations, verifiable against community results.

### Files to Create
- `config/tpv102_params.hpp` — TPV102 parameter struct (all values from benchmark spec)
- `config/tpv102.toml` — TOML configuration file for TPV102
- `dynamic/tpv102_setup.hpp` — Initialization: spatial a(x,y), θ_ini(x,y), nucleation perturbation
- `drivers/tpv102_driver.cpp` — Main driver for standalone TPV102 simulation
- `tests/verification/tpv102_verification.cpp` — Verification test against reference data

### Detailed Requirements

1. **Mesh generation:** Half-space domain with fault at z=0.
   - Domain: [−50km, 50km] × [0, 30km] × [−50km, 50km] (x × y × z)
   - Fault plane: z=0, 0 < y < 30km (deep enough for VW+VS zones)
   - Free surface: y=0
   - Absorbing: all other boundaries
   - Resolution: ≈500m elements near fault, coarser away
   - Element type: hexahedra (compatible with existing infrastructure) or tetrahedra

   **Coordinate mapping from TPV102 spec to our system:**
   - TPV102: fault at z=0, free surface at y=0, y>0 is half-space
   - Our mesh: same convention (fault at z=0, y is depth from surface)
   - x = along-strike, y = depth (positive down from free surface), z = fault-normal

2. **`TPV102Params` struct:**
   ```cpp
   struct TPV102Params {
      // Material
      real_t rho = 2670.0, cs = 3464.0, cp = 6000.0;
      real_t mu() const { return rho * cs * cs; }
      real_t lambda() const { return rho * cp * cp - 2.0 * mu(); }

      // Friction
      real_t f0 = 0.6, V0 = 1e-6, b = 0.012, L = 0.02;
      real_t a_vw = 0.008, a_vs = 0.016, delta_a0 = 0.008;

      // Stress
      real_t sigma_n = 120e6, tau_ini = 75e6, V_ini = 1e-12;

      // Geometry
      real_t W = 15e3, w = 3e3;
      real_t x_hypo = 0.0, y_hypo = 7.5e3;

      // Nucleation
      real_t delta_tau0 = 25e6, R_nuc = 3e3, T_nuc = 1.0;

      // Simulation
      real_t t_final = 12.0;
   };
   ```

3. **Spatial parameter computation:**

   a(x,y) = 0.008 + Δa(x,y) where:
   ```
   Δa(x,y) = Δa₀ · [1 − B(x; W, w) · B(y − y₀; W/2, w)]
   ```

   B(x; W, w) = smooth C∞ boxcar function (Eq. 5 of spec):
   ```
   B(x; W, w) = 1                                          if |x| ≤ W
              = ½[1 + tanh(w/(|x|−W−w) + w/(|x|−W))]     if W < |x| < W+w
              = 0                                          if |x| ≥ W+w
   ```

   θ_ini(x,y): computed from Eq. 6 to satisfy friction law at t=0:
   ```
   θ_ini = (L/V₀) · exp[(a·ln(2·sinh(τ_ini/(a·σ_n))) − f₀ − a·ln(V_ini/V₀)) / b]
   ```

4. **Nucleation perturbation:** Applied as additional shear traction Δτ_x during simulation:
   ```
   Δτ(x,y,t) = Δτ₀ · F(r) · G(t)
   ```
   where r = √((x−x₀)² + (y−y₀)²). This is added to the fault traction at each time step for t < T.

5. **Station output:** SCEC format time series at 9 fault stations and 6 surface stations. Output interval: every time step or at fixed Δt = 0.01s.

### Acceptance Criteria
- [ ] Rupture nucleates at hypocenter and propagates bilaterally
- [ ] Rupture arrests at VW/VS boundary (within ±2 km)
- [ ] Peak slip rate at station (0, 7.5km) matches community median within 10%
- [ ] Rupture front arrival time at (0, 7.5km) matches community median within 5%
- [ ] Free surface displacement at station (0, 9km, z) shows expected P-wave and S-wave arrivals
- [ ] Simulation runs for full 12 seconds without instability

### Dependencies
- Depends on: Phase 1 (wave operator), Phase 2 (BCs), Phase 3 (fault Riemann)
- Required by: Phase 5 (QD↔FD transfer uses the same infrastructure)

---

## Phase 5: QD↔FD Transfer Architecture Design

### Goal
Design and implement the state conversion and snapshot protocol enabling future hybrid QD-dynamic simulations, as specified in refactoring plan v4 Section 12.1.

### Files to Create
- `dynamic/velocity_stress_state.hpp` — Container for Q = (σ, v) with conversion from displacement
- `dynamic/regime_switch.hpp` — QD→FD and FD→QD state transfer protocols
- `solver/seas_hybrid_operator.hpp` — Three-mode dispatcher skeleton (qd/dynamic/hybrid)
- `tests/unit/test_regime_switch.cpp` — Conversion accuracy tests

### Detailed Requirements

1. **`VelocityStressState` — Conversion between QD and FD representations:**

   ```cpp
   class VelocityStressState {
   public:
      /// Convert QD displacement field to velocity-stress field.
      /// σ = C : ∇u (element-local, no MPI needed)
      /// v = 0 (or smooth estimate from previous QD velocity)
      static void DisplacementToVelocityStress(
         const ParFiniteElementSpace &fes_disp,  // QD displacement space (H1 or DG)
         const ParGridFunction &displacement,
         const ConstitutiveModel &model,
         ParFiniteElementSpace &fes_wave,          // FD wave space (L2, 9-component)
         ParGridFunction &Q);

      /// Convert velocity-stress field back to displacement.
      /// u = ∫ v dt (accumulated during FD phase)
      /// Or: solve QD system with updated slip BC.
      static void VelocityStressToDisplacement(
         const ParFiniteElementSpace &fes_wave,
         const ParGridFunction &Q,
         ParFiniteElementSpace &fes_disp,
         ParGridFunction &displacement);

      /// Extract slip and state variable from FD fault data.
      static void ExtractFaultState(
         const FaultRiemann &fault_solver,
         Vector &slip,
         Vector &psi);
   };
   ```

2. **`RegimeSwitch` — Snapshot/restore protocol:**

   ```cpp
   struct QDSnapshot {
      Vector displacement;      // QD displacement field
      Vector slip;              // Fault slip vector
      Vector psi;               // State variable (theta or psi)
      real_t time;              // Snapshot time
   };

   class RegimeSwitch {
   public:
      /// Freeze QD state and convert to FD initial condition.
      void QDtoFD(const SEASQuasiDynamicOperator &qd_op,
                  const Vector &qd_state,
                  WaveOperator &wave_op,
                  Vector &fd_state,
                  QDSnapshot &snapshot);

      /// Restore QD state from snapshot + accumulated FD slip.
      void FDtoQD(const QDSnapshot &snapshot,
                  const Vector &fd_accumulated_slip,
                  const Vector &fd_final_psi,
                  SEASQuasiDynamicOperator &qd_op,
                  Vector &qd_state);
   };
   ```

3. **Hybrid operator skeleton** (functional implementation deferred to future phase):

   ```cpp
   // solver/seas_hybrid_operator.hpp
   class SEASHybridOperator : public TimeDependentOperator {
   public:
      enum class Regime { QD, FD, Transition };

      SEASHybridOperator(SEASQuasiDynamicOperator *qd_op,
                          WaveOperator *wave_op,
                          FaultRiemann *fault_solver,
                          real_t V_activate = 1e-3,
                          real_t V_deactivate = 1e-6);

      void Mult(const Vector &state, Vector &rate) const override;
      Regime GetCurrentRegime() const { return regime_; }

   private:
      Regime regime_ = Regime::QD;
      real_t V_activate_, V_deactivate_;
      // ... pointers to operators
   };
   ```

### Acceptance Criteria
- [ ] QD→FD conversion: stress field σ = C:∇u matches QD traction to machine precision
- [ ] Round-trip: QD→FD→QD with zero FD evolution returns original QD state (L2 < 1e-14)
- [ ] Fault state continuity: slip and θ are preserved across QD→FD transition
- [ ] `SEASHybridOperator` compiles and dispatches to QD mode (FD mode is stub that aborts with "not yet implemented")

### Dependencies
- Depends on: Phase 1-4 (FD infrastructure), existing QD infrastructure
- Required by: future hybrid simulation work

---

## Phase 6: GPU Acceleration Exploration

### Goal
Profile the CPU dynamic rupture code, identify bottlenecks, and implement GPU-accelerated versions of critical kernels using MFEM's device infrastructure.

### Files to Create
- `dynamic/wave_operator_gpu.cpp` — GPU-accelerated volume + face kernels
- `tests/unit/test_wave_gpu.cpp` — GPU vs CPU equivalence test
- `scripts/profile_tpv102.sh` — Profiling script

### Detailed Requirements

1. **Profiling:** Run TPV102 on CPU, profile with `perf`/`gprof` to identify bottlenecks. Expected:
   - ~60-70% in face flux computation (Godunov + fault Riemann)
   - ~20-30% in volume integral
   - ~5-10% in mass inverse application + time stepping overhead

2. **MFEM device infrastructure to use:**
   - `mfem::forall(N, [=] MFEM_HOST_DEVICE (int i) { ... })` for portable kernels
   - `Device::Configure("cuda")` or `Device::Configure("hip")` at startup
   - `Memory<real_t>` with `UseDevice()` for GPU-resident arrays
   - PA (Partial Assembly) mode for face integrals: precompute geometric data, apply matrix-free

3. **GPU kernels to implement:**

   a. **Element mass inverse application** (embarrassingly parallel per element):
   ```cpp
   // Per element: Q_new = M_inv * RHS
   MFEM_FORALL(e, num_elements, {
      // Load element RHS, multiply by precomputed M_inv, store
   });
   ```

   b. **Volume integral** (per element):
   ```cpp
   MFEM_FORALL(e, num_elements, {
      // Load Q for element e
      // Compute F_j(Q) at quadrature points
      // Contract with ∂Φ/∂x_j (precomputed)
      // Accumulate into RHS
   });
   ```

   c. **Face flux** (per face):
   ```cpp
   MFEM_FORALL(f, num_faces, {
      // Load Q from both sides at face quadrature points
      // Compute Godunov flux (or fault Riemann solve)
      // Accumulate into RHS for both elements
   });
   ```

   d. **Fault Riemann solver** (per fault face DOF):
   ```cpp
   MFEM_FORALL(i, num_fault_dofs, {
      // Load Q_plus, Q_minus at DOF i
      // Compute trial traction
      // Newton-Raphson for slip rate (branch-free version)
      // Update state variable
      // Compute imposed state
   });
   ```

4. **Performance targets:**
   - GPU speedup: ≥10× over single-core CPU for TPV102 at order 4
   - GPU vs CPU equivalence: L2 difference < 10⁻¹² (FP precision)

### Acceptance Criteria
- [ ] Profiling report documents CPU bottleneck breakdown
- [ ] GPU wave operator produces identical results to CPU (L2 < 10⁻¹²)
- [ ] GPU achieves ≥5× speedup over single-core CPU for TPV102
- [ ] Code compiles and runs correctly with `Device("cpu")` when no GPU available

### Dependencies
- Depends on: Phase 1-4 (working CPU code), MFEM compiled with CUDA/HIP support
- Required by: nothing (final optimization phase)

---

## Testing Strategy

### Per-Phase Testing

| Phase | Test | Type | Criterion |
|-------|------|------|-----------|
| 1 | Plane wave propagation | Unit | Speed matches cₚ, cₛ within 1% |
| 1 | Order convergence | Unit | L2 error ∝ h^{N+1} |
| 1 | Energy conservation | Unit | ΔE/E < 10⁻¹² in closed domain |
| 2 | Free surface reflection | Unit | P→P, P→SV amplitudes correct |
| 2 | Absorbing boundary | Unit | < 1% reflected energy |
| 3 | Locked fault | Unit | No slip, traction = Godunov state |
| 3 | Frictionless fault | Unit | Full stress drop |
| 3 | 1D SH wave + fault | Verification | Matches analytical solution |
| 4 | TPV102 nucleation | Verification | Rupture starts at t < 2s |
| 4 | TPV102 full | Verification | Matches community within 10% |
| 5 | QD↔FD round-trip | Unit | L2 < 10⁻¹⁴ |
| 6 | GPU equivalence | Unit | CPU vs GPU L2 < 10⁻¹² |

### Reference Data

- SeisSol TPV102 results (generate from SeisSol code at `/Users/chunhuizhao/projects/SeisSol`)
- SCEC community results (download from https://strike.scec.org/cvws/results.html)
- Analytical solutions for plane waves, 1D fault problems

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| Godunov flux implementation error | Medium | High | Phase 1 convergence test catches immediately; compare with SeisSol flux at same face |
| Fault Riemann solver Newton divergence | Medium | High | Bisection fallback; test with extreme parameters (V→0, σ_n→0) |
| CFL instability at fault faces | Medium | High | Conservative CFL safety factor; monitor energy balance |
| Mesh resolution insufficient for TPV102 | Low | Medium | Start with coarse mesh for smoke test, refine to 500m for verification |
| Parallel shared-face ordering for fault | Medium | Medium | Reuse existing `FaultBasis`/`FaultGeometry` infrastructure proven for BP5 |
| GPU memory overflow | Low | Medium | Start with small meshes; profile memory usage |
| MFEM L2 space performance for 9-component system | Medium | Medium | Benchmark against separate scalar spaces; optimize if needed |
| QD→FD stress field mismatch from DG vs H1 spaces | High | High | Phase 5 test catches this; may need L2 projection step |

### Known Tricky Areas

1. **Sign conventions at fault:** Our existing code uses `sign = (nor(0) > 0) ? -1 : 1`. The velocity-stress formulation's plus/minus convention must be consistent with this.

2. **Coordinate mapping:** TPV102 uses (x, y, z) where y=depth, z=fault-normal. Our BP5 code uses Y=fault-normal, Z=depth. Need explicit mapping in `tpv102_setup.hpp`.

3. **Friction solver stability:** The regularized arcsinh R&S law has numerical issues when V·exp(c) overflows double precision. SeisSol uses a stable `arsinhexp` function — must replicate.

4. **Mass matrix inversion:** For non-orthogonal basis on deformed hexahedra, M is not diagonal. Must use actual dense inverse per element (O(ndof³) per element, but ndof is small).

5. **Parallel fault face communication:** In the velocity-stress formulation, both sides of the fault are separate elements. If they're on different MPI ranks, need ghost exchange before Riemann solve. Reuse existing `FaultScatter` infrastructure.

---

## Project Layout After All Phases

```
miniapps/seas/
├── dynamic/                           NEW (all phases)
│   ├── wave_state.hpp                 Phase 1: Q = (σ,v) container
│   ├── wave_flux.hpp                  Phase 1: Godunov flux
│   ├── wave_operator.hpp              Phase 1: DG wave operator
│   ├── wave_operator.cpp              Phase 1: implementation
│   ├── wave_operator_gpu.cpp          Phase 6: GPU kernels
│   ├── fault_riemann.hpp              Phase 3: fault Riemann solver
│   ├── fault_riemann.cpp              Phase 3: implementation
│   ├── aging_law_dr.hpp               Phase 3: aging law update
│   ├── tpv102_setup.hpp               Phase 4: TPV102 initialization
│   ├── velocity_stress_state.hpp      Phase 5: QD↔FD conversion
│   └── regime_switch.hpp              Phase 5: snapshot/restore
├── config/
│   ├── tpv102_params.hpp              Phase 4: TPV102 parameters
│   └── tpv102.toml                    Phase 4: TOML config
├── drivers/
│   └── tpv102_driver.cpp              Phase 4: standalone driver
├── solver/
│   └── seas_hybrid_operator.hpp       Phase 5: three-mode dispatcher
└── tests/
    ├── unit/
    │   ├── test_wave_flux.cpp         Phase 1
    │   ├── test_wave_operator.cpp     Phase 1
    │   ├── test_wave_bc.cpp           Phase 2
    │   ├── test_fault_riemann.cpp     Phase 3
    │   ├── test_regime_switch.cpp     Phase 5
    │   └── test_wave_gpu.cpp          Phase 6
    └── verification/
        └── tpv102_verification.cpp    Phase 4
```

**New files: ~20. New LOC estimate: ~5,000–8,000.**

---

## Execution Order and Dependencies

```
Phase 1 (Wave Operator) ────────────────► Phase 2 (BCs) ──► Phase 4 (TPV102)
                                                    ▲              │
Phase 3 (Fault Riemann) ───────────────────────────┘              │
                                                                    ▼
                                              Phase 5 (QD↔FD Transfer)
                                                                    │
                                              Phase 6 (GPU Accel) ◄─┘
```

Phase 1 and Phase 3 can proceed in parallel. Phase 2 and Phase 3 must both complete before Phase 4. Phase 5 requires Phase 4 (proven FD code) plus existing QD code. Phase 6 requires Phase 4 (working CPU code to compare against).
