# Spectral Boundary Integral (SBI) Implementation Plan — Pure SBI, Antiplane

**Prepared:** 2026-02-28
**Scope:** Implement `AntiplaneSBIOperator` — a pure SBI domain operator for quasi-dynamic
antiplane earthquake cycle simulation. No FEM domain; only the fault is discretized.
FFT library: FFTW3.

**Reference paper:** Abdelmeguid et al. (2019), JGR — "A Novel Hybrid Finite
Element-Spectral Boundary Integral Scheme for Modeling Earthquake Cycles:
Application to Rate and State Faults With Low-Velocity Zones."

**Sequel:** Hybrid FE-SBI (for low-velocity zones) is a separate document.

---

## Part A: Physics and Mathematical Formulation

### A.1 Why SBI Instead of FEM?

The current MFEM SEAS miniapp solves the quasi-dynamic antiplane equilibrium equation
using a **Discontinuous Galerkin FEM** domain operator over a large 2D box (±80 km × 80 km).

| Property | FEM (current) | SBI (proposed) |
|---|---|---|
| Domain | 2D box, ±80 km | 1D fault line only |
| Cost per step | O(N^1.5) — MUMPS | O(N log N) — FFT |
| Heterogeneity | Arbitrary | Homogeneous half-space only |
| Accuracy | Finite domain truncation error | Exact half-space Green's function |
| Mesh | Required (large, 963k+ elements for 25m) | Not needed |
| Fault DOFs at 25m | ~600 (extracted from 2D mesh) | ~600 (1D directly) |

For a **homogeneous half-space** (BP1, BP2), the SBI gives the exact traction from
the half-space Green's function evaluated in Fourier space. No mesh is needed.
The speedup over FEM is 2–3 orders of magnitude per time step.

### A.2 Quasi-dynamic SBI Kernel

For a 2D antiplane fault in a homogeneous half-space, periodic with period λ, the
shear traction from a given slip distribution in Fourier space is (Eq. 22 of paper,
quasi-dynamic limit):

```
τ̂(q, t) = π·μ·|q|·ŝ(q, t)
```

where:
- `q = 2π·k/λ`  discrete wavenumber, k = -N/2, …, N/2 (FFTW convention)
- `ŝ(q)`  = FFT of the slip distribution along the fault (1D array)
- `μ` = shear modulus

In real space: `τ_SBI(z, t) = iFFT( π·μ·|q| · FFT(s(z, t)) )`

This is **instantaneous** — no time-history convolution is needed in quasi-dynamic mode.
The full dynamic kernel (Eq. 20–21 in paper) would add a convolution over slip history;
that is deferred to a future dynamic extension.

Key properties of the kernel:
- `q_0 = 0` → `kernel[0] = 0` (constant slip produces zero stress gradient — correct)
- Monotonically increasing with wavenumber → short-wavelength slip produces large traction
- Symmetric: `kernel[-k] = kernel[k]` (real-valued traction)

### A.3 Tectonic Loading via Creep Zone

To impose far-field plate loading without a FEM domain, the fault is extended below
the rate-state zone (depth Wf) to a total depth Lf, with prescribed slip rate V_creep = Vp.
The SBI kernel computes the background stress at the seismogenic zone from this creep.

```
z ∈ [0,   Wf]:   rate-and-state friction   (N_rsf nodes, ~600 at 25m)
z ∈ [Wf,  Lf]:   prescribed V = Vp         (N_creep nodes)
SBI period:  λ = N_sbi · dz_sbi ≥ 2·Lf   (zero-padded to avoid aliasing)
```

Recommended parameters for BP1:
- `Wf = 15 km`, `Lf = 60 km` (= 4×Wf for safety), `λ = 120 km`
- `dz_sbi = 25 m` → `N_rsf = 600`, `N_creep = 1800`, `N_sbi = 4096` (next power of 2 ≥ 2·2400)

This approach is standard in quasi-dynamic SBI codes (QDYN, Tandem SBI mode).

### A.4 Radiation Damping

Identical to the FEM approach. At each fault node:
```
τ_net = τ_SBI + τ_0 - η·V
```
where `η = μ/(2·cs)`. No change to `RateStateFaultOperator::ComputeRHS()`.

---

## Part B: Architecture — Where SBI Fits

### B.1 Layer diagram

The SBI operator is a **drop-in replacement** for the FEM domain operator. Every layer
above it (SEAS coupling, fault physics, time stepping, output) is unchanged.

```
SEASQuasiDynamicOperator<Mesh>     [solver/seas_operator.hpp  — UNCHANGED]
  │
  ├── DomainOperator<Mesh>*         [currently: AntiplaneBdrLoadOperator]
  │    └── AntiplaneSBIOperator     [NEW — replaces FEM domain entirely]
  │          has: SBIKernel          [NEW — FFT engine, FFTW3]
  │          has: SBICreepZone       [NEW — below-Wf slip manager]
  │
  └── RateStateFaultOperator<Mesh>* [fault/rate_state_fault.hpp — UNCHANGED]
        has: FaultGeometry           [UNCHANGED]
        has: DieterichRuinaFriction  [UNCHANGED]

DormandPrinceRK45                  [solver/time_stepper.hpp   — UNCHANGED]
BenchmarkOutput / ParallelBenchmarkOutput  [io/  — UNCHANGED]
```

### B.2 Interface compatibility

`AntiplaneSBIOperator` implements all pure-virtual methods of `DomainOperator<Mesh>`:

| Method | SBI implementation |
|---|---|
| `Solve(t, slip, u)` | FFT of full slip (RSF+creep) → store traction; `u` unused |
| `ComputeTraction(u, slip, traction)` | Return stored traction from last `Solve()` call |
| `GetNumFaultDOFs()` | Returns `N_rsf` |
| `GetFaultDepths(depths)` | Returns `[0, dz, 2·dz, ..., Wf]` |
| `GetFaultDOFs()` | Returns `[0, 1, ..., N_rsf-1]` |
| `NumSlipComponents()` | Returns `1` (antiplane) |
| `Dimension()` | Returns `2` |

`SEASQuasiDynamicOperator`, `RateStateFaultOperator`, and all output classes are
called through `DomainOperator*` virtual dispatch and need **zero changes**.

### B.3 The dummy mesh issue

`SEASQuasiDynamicOperator<Mesh>` passes a `GridFuncType& u` (an MFEM grid function)
to `Solve()`. For SBI there is no 2D displacement field.

**Mitigation**: Construct a 1D segment mesh of `N_rsf` segments at construction time.
Allocate a scalar L2 grid function on it as a compatibility shim. Its values are
never read — it is passed by reference and ignored inside `AntiplaneSBIOperator::Solve()`.

This avoids any changes to `SEASQuasiDynamicOperator`. Long-term clean-up
(a `DomainOperator<void>` specialization) is deferred to after verification.

---

## Part C: Implementation Plan

### Phase 1a: `domain/sbi_kernel.hpp` — FFT Traction Engine

**Purpose**: Encapsulates all FFTW3 calls. Completely independent of MFEM and SEAS.
Testable in isolation.

**Class**: `SBIKernel`

```cpp
class SBIKernel {
public:
    // N_sbi: total grid points (power of 2, must equal full SBI domain size)
    // dz: grid spacing (m), mu: shear modulus (Pa), lambda = N_sbi * dz
    SBIKernel(int N_sbi, real_t dz, real_t mu);
    ~SBIKernel();

    // Core: traction[i] = iFFT( pi*mu*|q_k| * FFT(slip)[k] )
    // Both arrays have length N_sbi (real-valued)
    void ComputeTraction(const double* slip, double* traction) const;

    int    N()      const { return N_sbi_; }
    real_t Dz()     const { return dz_; }
    real_t Lambda() const { return lambda_; }   // = N_sbi * dz

private:
    int    N_sbi_;
    real_t dz_, mu_, lambda_;

    // Precomputed: kernel_[k] = pi * mu * |q_k|, FFTW half-complex layout
    std::vector<double> kernel_;

    // FFTW3 plans (real-to-real, using r2c/c2r for efficiency)
    double       *buf_r_;    // real buffer, length N_sbi
    fftw_complex *buf_c_;    // complex buffer, length N_sbi/2+1
    fftw_plan plan_r2c_;     // forward: real → complex
    fftw_plan plan_c2r_;     // inverse: complex → real
};
```

**Key implementation details**:

```
Constructor:
  lambda = N_sbi * dz
  for k = 0, 1, ..., N_sbi/2:
    q_k = 2*pi*k / lambda
    kernel_[k] = pi * mu * q_k       (q_k is already non-negative in r2c layout)
  kernel_[0] = 0.0                   (DC component: constant slip → zero traction)
  Create FFTW plans: FFTW_MEASURE

ComputeTraction(slip, traction):
  1. Copy slip → buf_r_
  2. Execute plan_r2c_: buf_r_ → buf_c_   (forward FFT, complex output)
  3. for k = 0..N_sbi/2: buf_c_[k] *= kernel_[k]   (multiply by |q| kernel)
  4. Execute plan_c2r_: buf_c_ → buf_r_   (inverse FFT, real output)
  5. Normalize: buf_r_[i] /= N_sbi       (FFTW unnormalized)
  6. Copy buf_r_ → traction
```

FFTW r2c/c2r pair is used because both slip and traction are real-valued.
This is ~2× faster than the complex-to-complex pair and avoids imaginary part handling.

**Build**: link `-lfftw3`. Header includes `<fftw3.h>`.
On TACC: `module load fftw3` before building.

---

### Phase 1b: `domain/sbi_creep_zone.hpp` — Creep Zone Manager

**Purpose**: Maintains the monotonically-increasing slip in the below-Wf creep zone,
assembles the full slip array (RSF + creep + zero padding) for the SBI kernel, and
extracts the RSF-zone traction from the full output.

**Class**: `SBICreepZone`

```cpp
class SBICreepZone {
public:
    // N_rsf:   number of RSF-zone nodes (z in [0, Wf])
    // N_creep: number of creep-zone nodes (z in [Wf, Lf])
    // N_sbi:   total SBI grid size (power of 2, ≥ N_rsf + N_creep)
    // Vp:      plate velocity (m/s)
    SBICreepZone(int N_rsf, int N_creep, int N_sbi, real_t Vp);

    // Advance creep: slip_creep[i] += Vp * dt  (called once per RK stage)
    void Advance(real_t dt);

    // Assemble: full_slip = [slip_rsf | slip_creep | 0...0], length = N_sbi
    // full_slip must be pre-allocated to N_sbi
    void AssembleFullSlip(const mfem::Vector& slip_rsf, double* full_slip) const;

    // Extract: traction_rsf[i] = full_traction[i]  for i = 0..N_rsf-1
    void ExtractRSFTraction(const double* full_traction,
                             mfem::Vector& traction_rsf) const;

    real_t TimeElapsed() const { return t_; }

private:
    int    N_rsf_, N_creep_, N_sbi_;
    real_t Vp_, t_;
    std::vector<double> slip_creep_;   // length N_creep, monotonically increasing
};
```

**Key implementation details**:
- `AssembleFullSlip`: copies `slip_rsf[0..N_rsf-1]` → `full_slip[0..N_rsf-1]`, then
  copies `slip_creep[0..N_creep-1]` → `full_slip[N_rsf..N_rsf+N_creep-1]`, then
  sets `full_slip[N_rsf+N_creep..N_sbi-1] = 0` (zero padding)
- `Advance`: `t_ += dt`, `slip_creep_[i] += Vp * dt` for all i
- `ExtractRSFTraction`: `traction_rsf[i] = full_traction[i]` for i = 0..N_rsf-1

The zero padding ensures the SBI periodicity `λ = N_sbi * dz` is large enough that
image faults (at z = λ, 2λ, …) do not significantly affect the traction at z ∈ [0, Wf].

---

### Phase 1c: `domain/antiplane_sbi_operator.hpp` — Domain Operator Subclass

**Purpose**: Implements the `DomainOperator<MeshType>` interface using `SBIKernel`
and `SBICreepZone`. The MFEM mesh parameter is accepted for type compatibility but
not used for computation; all physics is 1D on the fault.

**Class**: `AntiplaneSBIOperator`

```cpp
template <typename MeshType = mfem::Mesh>
class AntiplaneSBIOperator : public DomainOperator<MeshType>
{
public:
    // Physical parameters (all SI units)
    struct Params {
        real_t Wf;      // Rate-state fault depth (m), e.g. 15e3
        real_t Lf;      // Total fault length incl. creep zone (m), e.g. 60e3
        real_t dz_sbi;  // Node spacing along fault (m), e.g. 25.0
        int    N_sbi;   // Total SBI domain grid size (power of 2), e.g. 4096
        real_t mu;      // Shear modulus (Pa)
        real_t cs;      // Shear wave speed (m/s)
        real_t Vp;      // Plate velocity (m/s)
    };

    explicit AntiplaneSBIOperator(const Params& p);

    // Must be called after construction, before Solve/ComputeTraction
    void SetupFaultInfo();

    // ---- DomainOperator<MeshType> interface ----

    // 1. Advance creep zone by dt = t - t_prev_
    // 2. Assemble full slip (RSF + creep + padding)
    // 3. Call SBIKernel to compute full traction
    // 4. Extract RSF traction and store in traction_
    // (GridFuncType& u is accepted but unused — dummy for interface compat)
    void Solve(real_t t, const mfem::Vector& slip,
               GridFuncType& u) override;

    // Return stored traction_ from last Solve() call
    void ComputeTraction(const GridFuncType& /*u*/,
                          const mfem::Vector& /*slip_bc*/,
                          mfem::Vector& traction) override;

    int                  GetNumFaultDOFs() const override;
    void                 GetFaultDepths(mfem::Vector& depths) const override;
    const mfem::Array<int>& GetFaultDOFs() const override;
    int                  NumSlipComponents() const override { return 1; }
    int                  Dimension()         const override { return 2; }

private:
    Params p_;
    int    N_rsf_, N_creep_;
    real_t eta_;            // = mu / (2 * cs), radiation damping

    SBIKernel    sbi_kernel_;
    SBICreepZone creep_zone_;

    // 1D dummy mesh and grid function (for interface compatibility only)
    std::unique_ptr<mfem::Mesh>             dummy_mesh_;
    std::unique_ptr<mfem::FiniteElementCollection> dummy_fec_;
    std::unique_ptr<mfem::FiniteElementSpace>     dummy_fes_;

    mfem::Vector       traction_;         // RSF-zone traction, set in Solve()
    mfem::Vector       depths_;           // z-coords of RSF nodes: [0, dz, 2dz, ..., Wf]
    mfem::Array<int>   fault_dof_ids_;    // [0, 1, ..., N_rsf-1]

    std::vector<double> full_slip_;       // work buffer, length N_sbi
    std::vector<double> full_traction_;   // work buffer, length N_sbi

    real_t t_prev_;   // time at end of last Solve() call
};
```

**`SetupFaultInfo()` sequence**:
```
1. N_rsf   = (int) round(Wf / dz_sbi)
2. N_creep = (int) round((Lf - Wf) / dz_sbi)
3. Verify N_rsf + N_creep ≤ N_sbi / 2   (aliasing check)
4. Construct SBIKernel(N_sbi, dz_sbi, mu)
5. Construct SBICreepZone(N_rsf, N_creep, N_sbi, Vp)
6. Build depths_: depths_[i] = i * dz_sbi  for i = 0..N_rsf-1
7. Build fault_dof_ids_: [0, 1, ..., N_rsf-1]
8. Allocate traction_(N_rsf), full_slip_(N_sbi), full_traction_(N_sbi)
9. Build dummy_mesh_: 1D segment mesh with N_rsf elements (for GridFuncType compat)
10. t_prev_ = 0.0
```

**`Solve(t, slip, u)` sequence** (called 6× per RK step by DormandPrinceRK45):
```
dt = t - t_prev_
if dt > 0: creep_zone_.Advance(dt)
creep_zone_.AssembleFullSlip(slip, full_slip_.data())
sbi_kernel_.ComputeTraction(full_slip_.data(), full_traction_.data())
creep_zone_.ExtractRSFTraction(full_traction_.data(), traction_)
t_prev_ = t
```

**Pre-stress initialization**: `RateStateFaultOperator::Init()` calls `ComputeTraction()`
at t=0 with zero slip to get initial traction, which is zero for SBI (zero slip, zero creep).
Pre-stress `tau_pre` is then set to balance `sigma_n * f(V_init, theta_ss) + eta * V_init`.
This is handled entirely by `RateStateFaultOperator::Init()` — no SBI-specific code needed.

---

### Phase 1d: Driver — `tests/verification/bp1_sbi.cpp`

Follows `tests/verification/bp1_bdrload.cpp` structure verbatim, with only the
domain operator construction changed:

```cpp
// Replace:
// auto domain = std::make_unique<AntiplaneBdrLoadOperator<ParMesh>>(mesh, order, params);

// With:
AntiplaneSBIOperator<>::Params sbi_p;
sbi_p.Wf      = params.Wf;          // 15 km
sbi_p.Lf      = 4.0 * params.Wf;    // 60 km
sbi_p.dz_sbi  = 25.0;               // 25 m
sbi_p.N_sbi   = 4096;
sbi_p.mu      = params.mu;
sbi_p.cs      = params.cs;
sbi_p.Vp      = params.Vp;

auto domain = std::make_unique<AntiplaneSBIOperator<>>(sbi_p);
domain->SetupFaultInfo();
```

Everything downstream — `FaultGeometry`, `RateStateFaultOperator`, `SEASQuasiDynamicOperator`,
`DormandPrinceRK45`, `BenchmarkOutput` — is constructed identically to `bp1_bdrload.cpp`.
No other changes.

**Makefile additions**:
```makefile
BP1_SBI_SRC = tests/verification/bp1_sbi.cpp
BP1_SBI_OBJ = $(BP1_SBI_SRC:.cpp=.o)

seas_bp1_sbi: $(BP1_SBI_OBJ)
    $(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(BP1_SBI_OBJ) $(MFEM_LIBS) -lfftw3

$(BP1_SBI_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
    @mkdir -p $(@D)
    $(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-bp1-sbi: seas_bp1_sbi
    mpirun -np 1 ./seas_bp1_sbi --output-dir bp1/results_sbi_test
```

---

## Part D: Verification Plan

### D.1 Unit Tests (`tests/unit/test_sbi_kernel.cpp`)

All tests follow the existing `TEST_ASSERT` / `TEST_NEAR` framework.

#### Test 1 — DC component: constant slip → zero traction
```
slip[i] = 1.0  for all i
Expected: traction[i] = 0 everywhere  (|q_0| = 0)
TEST_NEAR(max|traction|, 0.0, 1e-10, "constant slip → zero SBI traction")
```

#### Test 2 — Single sinusoidal mode k=1
```
slip[i] = A * sin(2π * i / N_sbi)
Expected amplitude: π * μ * (2π/λ) * A
For each i:
  TEST_NEAR(traction[i], expected_amp * sin(2π*i/N_sbi), 1e-8, "k=1 traction amplitude")
```

#### Test 3 — Superposition of two modes (k=1 and k=5)
```
slip[i] = A1*sin(2π*1*i/N) + A5*sin(2π*5*i/N)
Expected: traction = π*μ*(|q_1|*A1*sin(q_1*z) + |q_5|*A5*sin(q_5*z))
FFT the output traction and check amplitude at each mode independently.
```

#### Test 4 — Skew-symmetry: sin input → sin output (no phase shift)
```
If slip is purely a sin wave at wavenumber k, traction must also be sin at k
(no mixing between modes — this verifies the kernel is diagonal in Fourier space)
```

#### Test 5 — Eshelby estimate for uniform finite dislocation
```
slip[i] = s0  for i in [0, N_rsf-1],  0 elsewhere
Center traction (at i = N_rsf/2) should be approximately:
  Δτ ≈ μ * s0 / (π * Wf)     (antiplane Eshelby estimate)
TEST_NEAR(traction[N_rsf/2], mu*s0/(M_PI*Wf), 0.15 * mu*s0/(M_PI*Wf),
          "Eshelby estimate within 15% (periodicity correction expected)")
```

### D.2 Unit Tests (`tests/unit/test_sbi_creep_zone.cpp`)

#### Test 6 — Creep zone advance
```
creep.Advance(1e7);   // 1e7 s ≈ 0.3 yr
creep.AssembleFullSlip(zeros_N_rsf, full_slip);
Expected: full_slip[N_rsf + i] = Vp * 1e7 = 0.01 m  for i = 0..N_creep-1
          full_slip[j] = 0.0  for j ≥ N_rsf + N_creep  (zero padding)
```

#### Test 7 — RSF traction extraction
```
Set full_traction = [0,1,2,...,N_sbi-1] (sentinel values)
creep.ExtractRSFTraction(full_traction, traction_rsf)
TEST_ASSERT(traction_rsf[i] == i  for i = 0..N_rsf-1)
```

#### Test 8 — Elapsed time tracking
```
creep.Advance(100.0);
creep.Advance(200.0);
TEST_NEAR(creep.TimeElapsed(), 300.0, 1e-10, "elapsed time accumulates")
```

### D.3 Unit Tests (`tests/unit/test_sbi_operator.cpp`)

#### Test 9 — SetupFaultInfo populates correct sizes
```
op.SetupFaultInfo();
TEST_ASSERT(op.GetNumFaultDOFs() == expected_N_rsf)
Vector depths; op.GetFaultDepths(depths);
TEST_ASSERT(depths.Size() == expected_N_rsf)
TEST_NEAR(depths[0], 0.0, 1e-10, "first depth = 0")
TEST_NEAR(depths[N_rsf-1], Wf - dz_sbi, 1e-6, "last depth near Wf")
```

#### Test 10 — Solve with zero slip returns zero traction
```
Vector slip(N_rsf);  slip = 0.0;
GridFuncType u;  // dummy
op.Solve(0.0, slip, u);
Vector traction(N_rsf);
op.ComputeTraction(u, slip, traction);
TEST_NEAR(traction.Norml2(), 0.0, 1e-10, "zero slip → zero SBI traction at t=0")
```

#### Test 11 — Traction grows as creep advances
```
Vector slip(N_rsf);  slip = 0.0;   // RSF zone stays locked
GridFuncType u;
op.Solve(1e8, slip, u);   // 3 years of creep at Vp
op.ComputeTraction(u, slip, traction);
// Creep slip = Vp * 1e8 = 0.1 m → should produce positive traction in RSF zone
TEST_ASSERT(traction.Max() > 0.0, "creep slip → positive loading traction")
```

### D.4 Benchmark 1: BP1 SBI vs BP1 FEM-bdrload

**Command**:
```bash
# SBI run (serial, ~seconds per step instead of minutes)
./seas_bp1_sbi --output-dir bp1/results_sbi_25m

# FEM baseline (25m self-similar mesh)
mpirun -np 300 ./seas_bp1_bdrload \
    --mesh bp1/mesh/bp1_ss_25m.msh \
    --output-dir bp1/results_ss_25m
```

**Comparison metrics** (at all 15 BP1 probe depths):

| Metric | Target tolerance |
|---|---|
| Time of first earthquake | < 2% relative difference |
| Peak coseismic V_max | Same order of magnitude |
| Interseismic V_max | < 5% relative difference |
| Earthquake recurrence interval | < 2% relative difference |
| Pre-stress τ_0 at 7.5 km depth | < 0.1% (purely parameter-dependent, not method-dependent) |

Expected discrepancies and their causes:
- FEM truncates domain at ±80 km; SBI is exact but periodic with λ = 120 km
- Different fault node locations (FEM: mesh-based, SBI: uniform 25 m)
- Differences < 1–2% are expected and acceptable for the homogeneous case

### D.5 Benchmark 2: SBI Grid Convergence

Run BP1-SBI at multiple resolutions and compare to SCEC reference data:

```
dz_sbi = 200m → N_rsf = 75,   N_sbi = 512
dz_sbi = 100m → N_rsf = 150,  N_sbi = 1024
dz_sbi = 50m  → N_rsf = 300,  N_sbi = 2048
dz_sbi = 25m  → N_rsf = 600,  N_sbi = 4096
dz_sbi = 12.5m→ N_rsf = 1200, N_sbi = 8192
```

Compare first earthquake time and recurrence interval vs SCEC BP1 reference.
Plot `|error| vs dz_sbi` on log-log scale; expect slope ≈ 1 (midpoint quadrature).

### D.6 Benchmark 3: Performance Comparison

```bash
# Time 1000 steps of SBI (serial)
time ./seas_bp1_sbi --max-steps 1000 --output-dir /dev/null

# Time 1000 steps of FEM (parallel, 25m mesh)
time mpirun -np 300 ./seas_bp1_bdrload --mesh bp1/mesh/bp1_ss_25m.msh --max-steps 1000
```

Expected: SBI is **100–1000× faster** per step (FFT O(N log N) vs MUMPS O(N^1.5)).

---

## Part E: Files Summary

### New files (to be created)

| File | Purpose | External dependency |
|---|---|---|
| `domain/sbi_kernel.hpp` | FFT-based SBI traction kernel | FFTW3 (`-lfftw3`) |
| `domain/sbi_creep_zone.hpp` | Creep-zone slip manager | none |
| `domain/antiplane_sbi_operator.hpp` | `DomainOperator` subclass | `sbi_kernel.hpp`, `sbi_creep_zone.hpp` |
| `tests/unit/test_sbi_kernel.cpp` | Unit tests for FFT kernel | FFTW3 |
| `tests/unit/test_sbi_creep_zone.cpp` | Unit tests for creep zone | none |
| `tests/unit/test_sbi_operator.cpp` | Unit tests for full operator | FFTW3 |
| `tests/verification/bp1_sbi.cpp` | BP1 verification driver | FFTW3 + all SEAS headers |

### Existing files modified (minimal)

| File | Change |
|---|---|
| `Makefile` | Add `seas_bp1_sbi` target and `test-sbi-*` targets; link `-lfftw3` |

### Files requiring zero changes

```
domain/domain_operator.hpp        fault/rate_state_fault.hpp
fault/fault_geometry.hpp          fault/fault_nodes.hpp
solver/seas_operator.hpp          solver/time_stepper.hpp
solver/seas_bdrload_operator.hpp  io/benchmark_output.hpp
io/parallel_benchmark_output.hpp  io/checkpoint.hpp
friction/dieterich_ruina.hpp      friction/state_evolution.hpp
config/bp1_params.hpp             config/bp2_params.hpp
common/seas_types.hpp             common/mpi_context.hpp
```

---

## Part F: Implementation Order

Execute in this order to enable incremental testing at each step:

1. **`sbi_kernel.hpp`** → build `test_sbi_kernel.cpp` → run tests (Tests 1–5)
   - Zero SEAS/MFEM dependencies: can test with a standalone compile + FFTW3
2. **`sbi_creep_zone.hpp`** → build `test_sbi_creep_zone.cpp` → run tests (Tests 6–8)
   - Only MFEM `Vector` dependency; no FFTW3
3. **`antiplane_sbi_operator.hpp`** → build `test_sbi_operator.cpp` → run tests (Tests 9–11)
   - First full integration with MFEM interface
4. **`bp1_sbi.cpp`** → build `seas_bp1_sbi` → short run (10 years, quick check)
5. **Full BP1 SBI run** (3000 yr) → compare to FEM baseline (Benchmark 1)
6. **Grid convergence study** (Benchmark 3: 5 resolutions)
7. **Performance timing** (Benchmark 2)

---

## Part G: Key Risks and Mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| Dummy mesh incompatibility with `SEASQuasiDynamicOperator` | Medium | Build with serial `Mesh` only first; test thoroughly before adding `ParMesh` variant |
| Periodicity aliasing (image faults affect RSF zone) | Low | Use λ ≥ 4×Wf; validate with Eshelby unit test (Test 5) |
| FFTW3 not available on build machine | Low | Add `#ifdef SEAS_USE_FFTW3` guard; allow compile without SBI if FFTW3 absent |
| Creep zone dt accumulation error in RK45 stages | Medium | Use `t_elapsed` from RK stage time directly, not accumulated `dt`; compare to FEM at interseismic phase |
| Pre-stress initialization returns zero traction (t=0, zero slip, zero creep) | Expected | `RateStateFaultOperator::Init()` handles this correctly — see Phase 1c notes |

---

## Part H: Future Extension — Hybrid FE-SBI

Once the pure SBI is verified, the hybrid FE-SBI extension (for LVZ problems from
the paper) adds:

1. A **small FEM domain** (near-fault strip, ±W_s ≈ 5 km) using `AntiplaneDomainOperator`
2. **SBI traction at virtual boundaries** (x = ±W_s) applied as Neumann BC via
   `BoundaryLFIntegrator` on the FEM domain's FARFIELD_LEFT/RIGHT boundaries
3. Displacement from FEM at virtual boundaries fed back to SBI kernel each step
4. Coupling iteration (predictor-corrector from Algorithm 1 of paper)

This is documented separately in `hybrid_sbi_implementation_plan_02282026.md`
(to be written after pure SBI is verified).
