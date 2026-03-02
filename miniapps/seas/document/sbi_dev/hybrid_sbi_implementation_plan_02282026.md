# Hybrid FE-SBI Implementation Plan — Near-Fault FEM + Far-Field SBI

**Prepared:** 2026-02-28
**Prerequisite:** Pure SBI (`AntiplaneSBIOperator`) must be verified against BP1 FEM
results before starting this document. See `sbionly_implementation_plan_02282026.md`.

**Reference paper:** Abdelmeguid et al. (2019), JGR — "A Novel Hybrid FE-SBI Scheme
for Modeling Earthquake Cycles: Application to Rate and State Faults with LVZ."
Algorithm 1 and Eqs. 9, 16–17, 34–35 are the implementation targets.

**Scope:** Implement `AntiplaneFEMSBIOperator` — a hybrid operator that uses:
- A **small near-fault FEM strip** (width ±W_s ≈ 1–5 km) for near-fault heterogeneity
- **SBI at virtual boundaries** (x = ±W_s) as exact DtN boundary conditions

This enables: Low-velocity zones (LVZ), complex fault geometry, or any near-fault
structure that the homogeneous SBI cannot represent.

---

## Part A: Physics of the Hybrid Coupling

### A.1 Domain Decomposition

The full antiplane half-space problem is split into three regions:

```
x: -∞           -W_s      0       +W_s          +∞
   ←── SBI ────── | ← FEM strip → | ────── SBI ──→
      left half     near-fault       right half
      space (S⁻)   FEM domain       space (S⁺)
      x < -W_s      [-W_s, W_s]     x > W_s
```

- **FEM strip** `[-W_s, W_s] × [-Lz, 0]`: Solved with DG-BR2 (same as current code).
  Contains the fault at x=0, and any near-fault heterogeneity (LVZ, damage zone).
- **Right far-field** `x > W_s`: Homogeneous elastic half-space.
  Represented exactly by the SBI Dirichlet-to-Neumann (DtN) map at x = W_s.
- **Left far-field** `x < -W_s`: Same as right, by symmetry.

The SBI at each virtual boundary converts the FEM displacement on that boundary
into the traction that the far-field half-space exerts on the FEM domain.

### A.2 Dirichlet-to-Neumann (DtN) Map at the Virtual Boundary

For the right far-field (x > W_s), homogeneous elastic half-space with the FEM domain
to its left:

**Problem**: Given `u(W_s, z) = u_D(z)` on the boundary x = W_s,
find the traction `T_x(z) = μ · ∂u/∂x |_{x=W_s}` (outward normal = +x̂ into far field).

**Solution**: The displacement satisfying ∇²u = 0 in x > W_s with
`u → 0` as x → ∞ and `u(W_s, z)` periodic in z with period λ is:

```
u(x, z) = Σ_k  Û_k · exp(-|q_k|·(x - W_s)) · e^(i·q_k·z)
```

where `q_k = 2π·k/λ` and `Û_k = FFT(u_D)[k]`.

The traction on the FEM domain from the right half-space (outward normal is -x̂ from
the right half-space perspective, so traction = +μ·∂u/∂x from FEM's outward +x̂):

```
T_right(z) = +μ · ∂u/∂x|_{x=W_s^+}
           = +μ · Σ_k  (-|q_k|) · Û_k · e^(i·q_k·z)
           = -μ · iFFT( |q| · FFT(u_D) )
```

**DtN kernel for right boundary**:
```
T_right(z) = -μ · iFFT( |q| · FFT(u(W_s, z)) )
```

Similarly for the left boundary (x = -W_s), outward normal = -x̂:
```
T_left(z)  = +μ · iFFT( |q| · FFT(u(-W_s, z)) )
```

The **sign difference** reflects the fact that both far-field half-spaces resist
deformation of the FEM domain: positive u_D on the right boundary causes a negative
(inward) traction from the right half-space, and vice versa.

### A.3 Comparison: DtN Kernel vs Fault SBI Kernel

| Kernel | Formula | Factor | Context |
|---|---|---|---|
| **Fault SBI** (pure SBI) | `τ̂(q) = π·μ·|q|·ŝ(q)` | π | Fault dislocation in half-space; both sides contribute |
| **Virtual boundary DtN** | `T̂(q) = -μ·|q|·û(q)` | 1 (no π) | One-sided half-space; standard DtN map |

The factor of π in the fault kernel comes from integrating the half-space Green's function
over both sides of the fault (the ±-side contributions add via the method of images).
The DtN map for a one-sided half-space does not have this factor.

### A.4 Tectonic Loading in the Hybrid Scheme

Tectonic loading is applied via a **creep zone on the fault**, identical to the pure SBI.
The far-field half-spaces outside ±W_s are NOT driven by far-field Dirichlet BCs — instead,
the creep zone below Wf generates the background loading through the fault SBI kernel
(now embedded in the FEM traction computation via the fault's slip jump BC).

This means:
- FEM domain has **natural BCs** (zero traction) at the virtual boundaries x = ±W_s,
  PLUS the **SBI DtN traction** from the far-field half-spaces.
- Below-fault creep zone V = Vp is still represented on the fault (inside the FEM domain).
- The SBI DtN traction at x = ±W_s provides the elastic interaction correction from the
  far-field, making the small FEM strip accurate despite its truncation.

For the FEM-only approach (current bdrload code), this far-field interaction is
approximated by Dirichlet: `u(±Lx) = ±Vp/2·t`. For the hybrid scheme, the SBI provides
the **exact** far-field elastic response at x = ±W_s without any domain truncation error.

### A.5 Quasi-Dynamic Coupling and Predictor-Corrector

At each time step, the coupling between FEM and SBI is explicit with one corrector step
(Algorithm 1 of the paper). For quasi-dynamic SEAS, this is:

```
Given: state at time t_n (u^n, slip^n, theta^n)
Target: state at time t_{n+1} = t_n + dt

Predictor step:
  1. Extrapolate u at virtual boundaries: u_bdr^* = 2*u^n - u^{n-1}  (linear extrapolation)
  2. Compute SBI traction from u_bdr^*: T_sbi^* = DtN(u_bdr^*)
  3. Apply T_sbi^* as Neumann BC; advance slip^n → slip^*
  4. Solve FEM: K·u^* = f(slip^*, T_sbi^*)
  5. Extract u_bdr^** = u^*(x=±W_s)

Corrector step:
  6. T_sbi^{n+1} = DtN(0.5*(u_bdr^* + u_bdr^**))  [averaged prediction]
  7. Solve FEM again: K·u^{n+1} = f(slip^{n+1}, T_sbi^{n+1})

Note: For quasi-dynamic SEAS with RK45, the "corrector" is built into the RK stages.
The predictor-corrector is simplified to: use u from the previous RK stage to compute
the SBI traction for the current stage. The RK45 error estimator handles accuracy.
```

---

## Part B: Architecture

### B.1 Layer Diagram

```
SEASQuasiDynamicOperator<Mesh>      [solver/seas_operator.hpp — UNCHANGED]
  │
  ├── DomainOperator<Mesh>*
  │    ├── AntiplaneSBIOperator      [pure SBI — Phase 1, already planned]
  │    ├── AntiplaneBdrLoadOperator  [FEM + Dirichlet loading — current]
  │    └── AntiplaneFEMSBIOperator   [NEW — FEM strip + SBI DtN at boundaries]
  │          has: AntiplaneDomainOperator  (near-fault FEM — reuse existing)
  │          has: SBIDtNKernel (right)     [NEW — DtN map for right half-space]
  │          has: SBIDtNKernel (left)      [NEW — DtN map for left half-space]
  │          has: SBICreepZone             [reuse from sbionly]
  │
  └── RateStateFaultOperator<Mesh>*  [UNCHANGED]
```

Key design: `AntiplaneFEMSBIOperator` **wraps** the existing `AntiplaneDomainOperator`
(or `AntiplaneBdrLoadOperator`) and intercepts the `Solve()` call to:
1. Extract virtual boundary displacement from the last FEM solution
2. Compute SBI DtN tractions at the virtual boundaries
3. Apply those tractions as Neumann BCs in the FEM RHS
4. Call the underlying FEM `Solve()`

The FEM operator itself is unchanged. The SBI acts as a time-dependent Neumann
BC supplier that replaces the far-field Dirichlet approximation.

### B.2 Mesh Change: Small Strip vs Large Domain

| Parameter | Current FEM (bdrload) | Hybrid FE-SBI |
|---|---|---|
| Domain width | ±80 km (160 km total) | ±W_s ≈ ±2–5 km (4–10 km total) |
| Domain depth | 80 km | 80 km (unchanged) |
| Elements at 25m | ~963k | ~12k–30k (40× smaller) |
| MUMPS cost | O(963k^1.5) | O(25k^1.5) ≈ 350× faster |
| Far-field accuracy | Truncation error | Exact (SBI DtN) |

The small strip uses the same `.geo` file (`bp1_selfsimilar.geo`) with a
much smaller domain parameter D:

```bash
# Current (80 km):
gmsh -2 bp1/mesh/bp1_selfsimilar.geo -o bp1/mesh/bp1_ss_25m.msh -setnumber h 0.025 -setnumber D 80

# Hybrid FE-SBI (5 km strip):
gmsh -2 bp1/mesh/bp1_selfsimilar.geo \
     -o bp1/mesh/bp1_hybrid_5km_25m.msh \
     -setnumber h 0.025 -setnumber D 5 -setnumber W 4 -setnumber Zf 45 -setnumber r 1.0
```

Setting `D = W_s = 5` and `W = 4` (uniform zone extends nearly to the boundary).
Setting `r = 1.0` (no coarsening needed since the domain is small).
Result: ~12k elements at 25m resolution.

---

## Part C: Implementation Plan

### Phase 2a: `domain/sbi_dtn_kernel.hpp` — DtN Kernel for Virtual Boundary

**Purpose**: Computes the SBI Dirichlet-to-Neumann map at a vertical virtual boundary.
Very similar to `SBIKernel` from the pure SBI, but with different sign and without the
factor of π. Can share the FFT code from `SBIKernel` via a common base or utility.

**Class**: `SBIDtNKernel`

```cpp
class SBIDtNKernel {
public:
    // N_bdr: number of boundary DOFs along z-direction
    // dz: node spacing (m), mu: shear modulus (Pa)
    // side: +1 for right boundary (traction = -mu*|q|*u), -1 for left (+mu*|q|*u)
    SBIDtNKernel(int N_bdr, real_t dz, real_t mu, int side);
    ~SBIDtNKernel();

    // traction[i] = side * (-mu) * iFFT(|q_k| * FFT(u_bdr)[k])
    void ComputeTraction(const double* u_bdr, double* traction) const;

    int N()   const { return N_bdr_; }
    real_t Dz() const { return dz_; }

private:
    int    N_bdr_;
    real_t dz_, mu_;
    int    side_;             // +1 or -1

    // Precomputed: kernel_[k] = mu * |q_k|  (no factor of pi, no factor of 2)
    std::vector<double> kernel_;

    // FFTW3 r2c/c2r plans
    double       *buf_r_;
    fftw_complex *buf_c_;
    fftw_plan plan_r2c_, plan_c2r_;
};
```

**Key differences from `SBIKernel`**:
- No factor of π: `kernel_[k] = mu * |q_k|` (not `pi * mu * |q_k|`)
- Sign: `traction = side * (-mu) * iFFT(|q| * FFT(u_bdr))`
  - Right boundary (side=+1): `T = -mu * iFFT(|q| * FFT(u))`
  - Left boundary (side=-1): `T = +mu * iFFT(|q| * FFT(u))`
- The FFTW3 code is identical to `SBIKernel`; consider a shared `FFTWBuffer` utility

**Design note**: `SBIDtNKernel` and `SBIKernel` differ only in the kernel prefactor.
Consider refactoring to a shared `SBIFFTBase` class with the FFTW3 boilerplate, and
`SBIKernel` / `SBIDtNKernel` as thin wrappers. This is optional — keep separate for
clarity if preferred.

---

### Phase 2b: Virtual Boundary DOF Extraction

To apply the DtN traction, we need:
1. **Extraction**: get the displacement values at x = +W_s (and x = -W_s) from the FEM
   grid function `u` — these are the boundary DOF values along z.
2. **Application**: apply the computed SBI traction as a Neumann BC in the FEM RHS.

**Extraction** uses MFEM's `GridFunction::GetValues()` or direct DOF access:

```cpp
// Extract boundary displacement at FARFIELD_RIGHT (x = +W_s)
// Returns a Vector of length N_bdr sorted by z-coordinate
void ExtractBoundaryDisp(const GridFuncType& u,
                          int bdr_attr,         // FARFIELD_RIGHT = 2
                          Vector& u_bdr) const;
```

Implementation: iterate over mesh boundary faces with `bdr_attr == FARFIELD_RIGHT`,
get the DOF indices, sort by z-coordinate, fill `u_bdr`.

For DG L2 space: boundary DOFs are the face-interior DOFs on the boundary faces.
For each boundary face, evaluate `u` at the face center quadrature point to get
one displacement value per face. The result is N_bdr values sorted by z.

**Application** of SBI traction as Neumann BC uses `BoundaryLFIntegrator`:

```cpp
// SBI DtN traction coefficient: wraps the precomputed traction array
class SBIDtNTractionCoefficient : public Coefficient {
public:
    void SetTraction(const Vector& traction, const Vector& z_coords);
    real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override;
    // Interpolates stored traction at ip.z using linear interpolation
private:
    const Vector* traction_;
    const Vector* z_coords_;
};
```

Applied to the FEM linear form:
```cpp
SBIDtNTractionCoefficient right_traction_coeff;
right_traction_coeff.SetTraction(t_sbi_right, z_right);

// Add Neumann BC contribution at FARFIELD_RIGHT
b.AddBdrFaceIntegrator(
    new BoundaryLFIntegrator(right_traction_coeff),
    farfield_right_marker_);
// (same for left boundary with left_traction_coeff)
```

**Key difference from current bdrload approach**: The current code applies a
Dirichlet BC via `DGDirichletLFIntegrator` (which adds both RHS and stiffness
penalty terms). The SBI applies a **pure Neumann BC** via `BoundaryLFIntegrator`
(no stiffness modification — the stiffness matrix stays constant as in `AntiplaneDomainOperator`).

---

### Phase 2c: `domain/antiplane_femsbi_operator.hpp` — Hybrid Domain Operator

**Purpose**: Combines the near-fault FEM domain with SBI DtN tractions at virtual
boundaries. Wraps the existing `AntiplaneDomainOperator` (which already handles the
fault slip BC, stiffness assembly, and MUMPS factorization) and adds SBI at the
outer boundaries.

**Class**: `AntiplaneFEMSBIOperator`

```cpp
template <typename MeshType = mfem::Mesh>
class AntiplaneFEMSBIOperator : public DomainOperator<MeshType>
{
public:
    struct Params {
        real_t Ws;       // Virtual boundary half-width (m), e.g. 5e3
        real_t Wf;       // Fault rate-state depth (m)
        real_t Lf;       // Total creep fault depth (m)
        real_t dz_sbi;   // SBI boundary node spacing (m, should match FEM dz_fault)
        int    N_bdr;    // Number of DtN boundary nodes (= Lz / dz_sbi, power of 2)
        real_t mu;       // Shear modulus (Pa)
        real_t cs;       // Shear wave speed (m/s)
        real_t Vp;       // Plate velocity (m/s)
    };

    AntiplaneFEMSBIOperator(MeshType& mesh, int order, const Params& p);

    void SetupFaultInfo() override;

    // Solve sequence:
    // 1. Compute SBI DtN traction from previous u at virtual boundaries
    // 2. Apply SBI traction as Neumann BC in FEM RHS
    // 3. Call internal FEM solve (MUMPS with small-strip stiffness matrix)
    void Solve(real_t t, const Vector& slip, GridFuncType& u) override;

    // Same as AntiplaneDomainOperator::ComputeTraction (unchanged)
    void ComputeTraction(const GridFuncType& u, const Vector& slip,
                          Vector& traction) override;

    // Delegate to internal FEM operator
    int GetNumFaultDOFs() const override;
    void GetFaultDepths(Vector& depths) const override;
    const Array<int>& GetFaultDOFs() const override;
    int NumSlipComponents() const override { return 1; }
    int Dimension() const override { return 2; }

private:
    Params p_;

    // Internal FEM operator (small-strip mesh, natural BCs at virtual boundaries)
    AntiplaneDomainOperator<MeshType> fem_op_;

    // SBI DtN kernels for left and right virtual boundaries
    SBIDtNKernel dtn_right_;  // side = +1, applied at FARFIELD_RIGHT
    SBIDtNKernel dtn_left_;   // side = -1, applied at FARFIELD_LEFT

    // SBI creep zone (same as pure SBI — fault-based tectonic loading)
    SBICreepZone creep_zone_;

    // Boundary DOF coordinates for interpolation
    Vector z_right_, z_left_;     // z-coords of DOFs on each virtual boundary
    int    N_bdr_right_, N_bdr_left_;

    // Work buffers
    Vector u_bdr_right_, u_bdr_left_;      // boundary displacement (extracted from u)
    std::vector<double> t_sbi_right_, t_sbi_left_;  // SBI traction at boundary
    GridFuncType u_prev_;                  // u from previous Solve call (for DtN input)

    // Coefficient objects for applying SBI traction as Neumann BC
    SBIDtNTractionCoefficient right_traction_coeff_, left_traction_coeff_;

    real_t t_prev_;
};
```

**`SetupFaultInfo()` sequence**:
```
1. Call fem_op_.SetupFaultInfo()   (identifies fault faces, assembles stiffness)
2. Find boundary DOFs on FARFIELD_RIGHT and FARFIELD_LEFT
3. Extract z-coordinates, sort by z, store in z_right_, z_left_
4. Construct dtn_right_(N_bdr, dz_sbi, mu, side=+1)
5. Construct dtn_left_(N_bdr, dz_sbi, mu, side=-1)
6. Construct creep_zone_(N_rsf, N_creep, N_sbi_creep, Vp)
   (fault-based tectonic loading, as in pure SBI)
7. Initialize u_prev_ = 0 (no previous solution)
```

**`Solve(t, slip, u)` sequence**:
```
1. Advance creep zone: creep_zone_.Advance(t - t_prev_)
   (this updates the fault-based tectonic loading, same as pure SBI)

2. Extract virtual boundary displacement from u_prev_:
   ExtractBoundaryDisp(u_prev_, FARFIELD_RIGHT, u_bdr_right_)
   ExtractBoundaryDisp(u_prev_, FARFIELD_LEFT,  u_bdr_left_)

3. Compute SBI DtN traction:
   dtn_right_.ComputeTraction(u_bdr_right_.data(), t_sbi_right_.data())
   dtn_left_.ComputeTraction(u_bdr_left_.data(),  t_sbi_left_.data())

4. Update Neumann BC coefficients:
   right_traction_coeff_.SetTraction(t_sbi_right_, z_right_)
   left_traction_coeff_.SetTraction(t_sbi_left_,  z_left_)

5. Rebuild FEM RHS with SBI Neumann BCs:
   (The stiffness matrix K is unchanged — only the RHS changes)
   fem_op_.RebuildRHSWithNeumannBC(slip, right_traction_coeff_, left_traction_coeff_, rhs)

6. Solve: K · u = rhs  (MUMPS backsubstitution, matrix already factored)

7. Store result: u_prev_ = u;   t_prev_ = t
```

**Critical detail: RHS modification**:
The internal `AntiplaneDomainOperator` currently builds the RHS as:
```cpp
rhs = -K_slip * slip    // slip BC contribution
```
For the hybrid operator, the RHS is extended to:
```cpp
rhs = -K_slip * slip                                  // slip BC (unchanged)
    + boundary_integral(T_sbi_right, FARFIELD_RIGHT)  // right SBI Neumann
    + boundary_integral(T_sbi_left,  FARFIELD_LEFT)   // left SBI Neumann
```

This requires exposing a method on `AntiplaneDomainOperator` to accept additional
RHS contributions, OR reimplementing the RHS assembly in the hybrid operator.
**Recommended approach**: Add a `SetExternalNeumannBC(LinearForm& extra_terms)` method
to `AntiplaneDomainOperator` that is added to the RHS at Solve time. This is a
minimal, backward-compatible addition.

---

### Phase 2d: Hybrid Driver — `tests/verification/bp1_hybrid_sbi.cpp`

```cpp
// Hybrid FE-SBI driver (replaces AntiplaneBdrLoadOperator with AntiplaneFEMSBIOperator)

AntiplaneFEMSBIOperator<>::Params hybrid_p;
hybrid_p.Ws       = 5.0e3;     // 5 km virtual boundary half-width
hybrid_p.Wf       = params.Wf; // 15 km
hybrid_p.Lf       = 60.0e3;    // 60 km (fault + creep)
hybrid_p.dz_sbi   = 25.0;      // 25 m (match FEM mesh at z-boundary)
hybrid_p.N_bdr    = 4096;      // boundary DtN nodes (must be ≥ Lz/dz = 3200)
hybrid_p.mu       = params.mu;
hybrid_p.cs       = params.cs;
hybrid_p.Vp       = params.Vp;

// Small-strip mesh (5 km half-width instead of 80 km):
// gmsh -2 bp1_selfsimilar.geo -setnumber D 5 -setnumber W 4 -setnumber h 0.025
auto mesh = LoadGmshMesh("bp1/mesh/bp1_hybrid_5km_25m.msh");

auto domain = std::make_unique<AntiplaneFEMSBIOperator<>>(mesh, order, hybrid_p);
domain->SetupFaultInfo();

// All other code identical to bp1_bdrload.cpp:
// FaultGeometry, RateStateFaultOperator, SEASQuasiDynamicOperator, BenchmarkOutput, ...
```

**New Makefile target**:
```makefile
BP1_HYBRID_SBI_SRC = tests/verification/bp1_hybrid_sbi.cpp
BP1_HYBRID_SBI_OBJ = $(BP1_HYBRID_SBI_SRC:.cpp=.o)

seas_bp1_hybrid_sbi: $(BP1_HYBRID_SBI_OBJ)
    $(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(BP1_HYBRID_SBI_OBJ) $(MFEM_LIBS) -lfftw3
```

**New mesh generation**:
```makefile
bp1/mesh/bp1_hybrid_5km_25m.msh:
    gmsh -2 bp1/mesh/bp1_selfsimilar.geo \
         -o bp1/mesh/bp1_hybrid_5km_25m.msh \
         -setnumber h 0.025 -setnumber D 5 \
         -setnumber W 4 -setnumber Zf 45 -setnumber r 1.0
```

---

### Phase 2e: LVZ Test Case

The key motivation for the hybrid scheme is the Low-Velocity Zone (LVZ). For LVZ:

```
Near-fault strip [-W_s, W_s]: μ_LVZ < μ (e.g. μ_LVZ / μ = 0.4, 0.6, 0.8)
Far-field half-spaces: μ (homogeneous, handled by SBI DtN)
```

The FEM operator handles heterogeneous μ via the `MatrixConstantCoefficient`
in the `DiffusionIntegrator`. The heterogeneous shear modulus is set using a
`FunctionCoefficient` that returns μ_LVZ for |x| ≤ W_s and μ elsewhere.

```cpp
// Spatially varying shear modulus for LVZ:
class LVZCoefficient : public Coefficient {
public:
    LVZCoefficient(real_t mu_lvz, real_t mu_far) : mu_lvz_(mu_lvz), mu_(mu_far) {}
    real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override {
        // All points in the FEM domain are in the LVZ (|x| < Ws)
        return mu_lvz_;
    }
private:
    real_t mu_lvz_, mu_;
};
```

The SBI DtN kernel always uses the far-field μ (homogeneous), while the FEM stiffness
uses μ_LVZ inside the strip. This is the correct physical model.

---

## Part D: Verification Plan

### D.1 Unit Test: DtN Kernel Correctness

**Test 1**: Constant u_bdr → zero traction
```
u_bdr[i] = 1.0  for all i
Expected: T = 0 everywhere (|q_0| = 0, DC term)
TEST_NEAR(max|T|, 0.0, 1e-10, "constant displacement → zero DtN traction")
```

**Test 2**: Single-mode DtN
```
u_bdr[i] = A * sin(2π * i / N_bdr)  (wavenumber k=1)
Expected: T[i] = ∓ mu * (2π/λ) * A * sin(2π * i / N_bdr)
  (- for right boundary, + for left boundary)
TEST_NEAR(max |T - T_expected|, 0.0, 1e-8, "k=1 DtN traction amplitude")
```

**Test 3**: Decaying solution consistency
```
Take a known Laplace solution: u(x, z) = A * e^{-|q|x} * sin(q*z)
At x = W_s: u_bdr[i] = A * e^{-|q|*Ws} * sin(q*z_i)
Expected traction at right boundary: T = -mu * (-|q|) * u_bdr = mu * |q| * u_bdr
Verify DtN gives this value.
```

**Test 4**: Anti-symmetry (left vs right)
```
Apply same u_bdr to left and right kernels
Expected: T_right = -T_left  (opposite signs)
TEST_NEAR(T_right + T_left, 0.0, 1e-10, "left-right DtN anti-symmetry")
```

### D.2 Benchmark 1: Hybrid vs Pure SBI for Homogeneous Case

**Goal**: For homogeneous material (μ_LVZ = μ), the hybrid FE-SBI must match the
pure SBI result exactly (up to mesh/SBI resolution differences).

```bash
# Pure SBI (reference)
./seas_bp1_sbi --output-dir bp1/results_sbi_25m

# Hybrid FE-SBI (5 km strip, same 25m fault resolution)
mpirun -np 4 ./seas_bp1_hybrid_sbi \
    --mesh bp1/mesh/bp1_hybrid_5km_25m.msh \
    --output-dir bp1/results_hybridsbi_5km_25m
```

Expected agreement:
- First earthquake time: < 1% difference
- Recurrence interval: < 1% difference
- Interseismic V_max: < 2% difference

Any difference > 2% indicates an error in the SBI DtN traction sign, prefactor,
or coupling algorithm.

### D.3 Benchmark 2: Hybrid vs FEM-bdrload for Homogeneous Case

Compare the hybrid scheme against the current large-domain FEM:

```bash
# Current FEM baseline (80 km domain)
mpirun -np 300 ./seas_bp1_bdrload \
    --mesh bp1/mesh/bp1_ss_25m.msh \
    --output-dir bp1/results_ss_25m

# Hybrid FE-SBI (5 km strip, exact far-field)
mpirun -np 4 ./seas_bp1_hybrid_sbi \
    --mesh bp1/mesh/bp1_hybrid_5km_25m.msh \
    --output-dir bp1/results_hybridsbi_5km_25m
```

Expected:
- Hybrid should MATCH the FEM within ~1–2%
- The hybrid has ~350× fewer elements → comparable to pure SBI speed for the FEM portion
- Total cost: small FEM solve (O(25k^1.5)) + two FFTs (O(N log N)) << large FEM (O(963k^1.5))

### D.4 Benchmark 3: LVZ — Comparison with Literature

Run BP1 with LVZ (μ_LVZ / μ = 0.4, 0.6, 0.8) and compare with Figure 4 of the paper:

```
For each μ_ratio in [0.4, 0.6, 0.8]:
    Modify FEM stiffness to use mu * mu_ratio inside the strip
    Run 3000-year simulation
    Compare earthquake recurrence interval vs paper Fig. 4
    Compare slip rate evolution at probe depths vs paper Fig. 5
```

This is the primary validation that the hybrid scheme is correct for heterogeneous media.
The pure SBI cannot run this test; only the hybrid can.

### D.5 Virtual Boundary Width Study

Test that results are insensitive to the choice of W_s (virtual boundary location):

```
Run hybrid FE-SBI with W_s = 1 km, 2 km, 5 km, 10 km (homogeneous case)
Expected: first earthquake time changes by < 0.5% as W_s varies
If strong dependence on W_s: suspect DtN prefactor error or coupling bug
```

---

## Part E: Files Summary

### New files (Phase 2)

| File | Purpose | Dependencies |
|---|---|---|
| `domain/sbi_dtn_kernel.hpp` | DtN map FFT kernel for virtual boundary | FFTW3 |
| `domain/antiplane_femsbi_operator.hpp` | Hybrid FEM+SBI domain operator | `sbi_dtn_kernel.hpp`, `sbi_creep_zone.hpp`, `antiplane_operator.hpp` |
| `tests/unit/test_sbi_dtn_kernel.cpp` | Unit tests for DtN kernel | FFTW3 |
| `tests/unit/test_femsbi_operator.cpp` | Unit tests for hybrid operator | MFEM + FFTW3 |
| `tests/verification/bp1_hybrid_sbi.cpp` | Hybrid driver for BP1 + LVZ | all SEAS headers + FFTW3 |
| `bp1/mesh/bp1_hybrid_5km_25m.msh` | Small-strip mesh (5 km × 80 km, 25m) | Gmsh |

### Existing files modified (minimal)

| File | Change |
|---|---|
| `domain/antiplane_operator.hpp` | Add `SetExternalNeumannBC()` method (additive, backward-compatible) |
| `Makefile` | Add `seas_bp1_hybrid_sbi` target, mesh generation rule |

### Files requiring zero changes

Same as pure SBI: all solver, fault, friction, io, and common headers.

---

## Part F: Implementation Order

Execute after pure SBI (Phase 1) is fully verified:

1. **`sbi_dtn_kernel.hpp`** + `test_sbi_dtn_kernel.cpp` → verify DtN formula (Tests 1–4)
2. **Boundary extraction** (`ExtractBoundaryDisp` utility) → test on existing small mesh
3. **`antiplane_operator.hpp`** minimal addition: `SetExternalNeumannBC()` method
4. **`antiplane_femsbi_operator.hpp`** → `test_femsbi_operator.cpp` (Tests 9–11 equivalent)
5. **`bp1_hybrid_5km_25m.msh`** generation → verify mesh has correct physical groups
6. **`bp1_hybrid_sbi.cpp`** driver → homogeneous test (Benchmark 1 vs pure SBI)
7. **Benchmark 2**: compare to large-domain FEM
8. **LVZ test** (Benchmark 3): modify μ_LVZ, compare to paper Figs. 4–5

---

## Part G: Key Risks and Mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| DtN prefactor wrong (missing factor of π or 2) | High | Unit Test 3 (decaying solution) directly checks the prefactor; also run W_s sensitivity study |
| DtN sign error (traction direction reversed) | Medium | Test 4 (left-right anti-symmetry); if sign wrong, earthquake nucleates on wrong side |
| Virtual boundary DOF extraction wrong z-ordering | Medium | Print z_right_ array and verify it matches the mesh boundary face z-coords |
| SBI traction not updating after first step (u_prev_ stuck) | Medium | Add diagnostic: print max|T_sbi| per step; should grow from zero initially |
| Mesh has wrong boundary attributes after D change | Low | Print mesh boundary attribute counts after loading; check FARFIELD_LEFT/RIGHT attributes |
| LVZ case: μ discontinuity at x=0 (fault) causes artifacts | Low | The DG formulation handles μ-jumps across element boundaries naturally via the diffusion coefficient |
| Coupling instability (explicit SBI traction at virtual boundary) | Low-Medium | For quasi-dynamic SEAS, explicit coupling is stable when dt << W_s/cs (wave transit time); monitor residual in corrector step |

---

## Part H: Connection to Three-Document Plan

This document is the third in a series:

1. **`mfem_seas_system_updateplan_02282026.md`**
   Preparatory refactoring (TOML input, fault physical groups, unit tests, build guide)
   → Must be done first to establish clean infrastructure

2. **`sbionly_implementation_plan_02282026.md`**
   Pure SBI quasi-dynamic operator — no FEM, homogeneous only
   → Verifies core FFT physics; provides performance baseline

3. **This document** (`hybrid_sbi_implementation_plan_02282026.md`)
   Hybrid FE-SBI — small FEM strip + SBI DtN at virtual boundaries
   → Enables LVZ, complex near-fault geometry; requires pure SBI verified first

The three documents form a complete path from the current state (large-domain FEM)
to the full hybrid method described in Abdelmeguid et al. (2019).
