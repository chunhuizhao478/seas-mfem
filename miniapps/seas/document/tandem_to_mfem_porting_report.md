# Tandem to MFEM Porting Analysis Report

## Executive Summary

This report analyzes the challenges and strategies for porting earthquake simulation functionality from **Tandem** to **MFEM**, focusing on SCEC SEAS benchmark problems (BP1, BP2). Both codebases implement Discontinuous Galerkin (DG) methods for elasticity, but with fundamentally different architectural philosophies.

**Key Finding**: MFEM provides excellent foundational infrastructure (DG integrators, time stepping, mesh handling), but lacks the specialized earthquake simulation components present in Tandem. The porting effort is substantial but feasible, with the main challenges being:
1. Fault interface treatment and friction law integration
2. Domain-fault coupling architecture
3. Quasi-dynamic vs. fully-dynamic mode switching
4. Configuration system differences

---

## 1. Architectural Comparison

### 1.1 Overall Design Philosophy

| Aspect | Tandem | MFEM |
|--------|--------|------|
| **Design Pattern** | Template-heavy, compile-time polymorphism | Runtime polymorphism, virtual interfaces |
| **Compile-time Config** | Polynomial degree, dimension fixed at compile | Runtime configurable |
| **DG Implementation** | Custom local operator pattern | BilinearFormIntegrator hierarchy |
| **Solver Integration** | PETSc-native | Pluggable (hypre, PETSc, native) |
| **Parallelism** | MPI + custom communication | MPI with ParMesh, ParFESpace |

### 1.2 Code Organization

**Tandem Structure:**
```
tandem/
├── src/           # Core library (basis, mesh, quadrature, parallel)
├── app/           # Application layer (operators, forms, kernels)
│   ├── localoperator/   # Elasticity, Poisson, Friction
│   ├── form/            # DG operators, coupling adapters
│   └── tandem/          # SEAS-specific drivers
└── examples/      # Benchmark configurations (BP1, BP3, BP5, BP6)
```

**MFEM Structure:**
```
mfem/
├── fem/           # Finite element (integrators, forms, spaces)
├── linalg/        # Linear algebra, ODE solvers
├── mesh/          # Mesh handling
├── miniapps/      # Applications
│   └── seas/      # Target location for porting
│       └── bp2/   # Benchmark data exists
└── examples/      # Reference DG examples (ex9, ex17)
```

---

## 2. DG Method Implementation Comparison

### 2.1 Tandem DG Framework

**Core Classes:**
- `DGOperator<LocalOp>` (src/form/DGOperator.h) - Template wrapper
- `DGCurvilinearCommon<D>` (src/form/DGCurvilinearCommon.h) - Curvilinear base
- `AbstractDGOperator<D>` (src/form/AbstractDGOperator.h) - Interface

**Assembly Pattern:**
```cpp
// Tandem uses SFINAE-detected methods on local operators
template<typename LocalOp>
class DGOperator {
    void assemble() {
        // Volume integrals
        if constexpr (has_assemble_volume<LocalOp>)
            localOp.assemble_volume(...);
        // Skeleton (interior face) integrals
        if constexpr (has_assemble_skeleton<LocalOp>)
            localOp.assemble_skeleton(...);
        // Boundary integrals
        if constexpr (has_assemble_boundary<LocalOp>)
            localOp.assemble_boundary(...);
    }
};
```

**Supported DG Methods:**
- Interior Penalty (IP/SIPG)
- Bassi-Rebay 2 (BR2)

### 2.2 MFEM DG Framework

**Core Classes:**
- `DGTraceIntegrator` (fem/bilininteg.hpp:3375) - Advection
- `DGDiffusionIntegrator` (fem/bilininteg.hpp:3503) - Diffusion (SIPG/NIPG)
- `DGDiffusionBR2Integrator` (fem/bilininteg.hpp:3591) - BR2 method
- `DGElasticityIntegrator` (fem/bilininteg.hpp:3691) - Elasticity

**Assembly Pattern:**
```cpp
// MFEM uses runtime polymorphism with BilinearForm
BilinearForm a(&fes);
a.AddDomainIntegrator(new ElasticityIntegrator(lambda, mu));
a.AddInteriorFaceIntegrator(
    new DGElasticityIntegrator(lambda, mu, alpha, kappa));
a.AddBdrFaceIntegrator(
    new DGElasticityIntegrator(lambda, mu, alpha, kappa), bdr_marker);
a.Assemble();
```

### 2.3 Porting Challenge: DG Operator Mapping

| Tandem Component | MFEM Equivalent | Gap Analysis |
|------------------|-----------------|--------------|
| `Elasticity` local op | `DGElasticityIntegrator` | Equivalent functionality |
| `Poisson` local op | `DGDiffusionIntegrator` | Equivalent functionality |
| `assemble_skeleton()` | `AddInteriorFaceIntegrator()` | Different API, same concept |
| `assemble_boundary()` | `AddBdrFaceIntegrator()` | Different API, same concept |
| Curvilinear geometry | `IsoparametricTransformation` | MFEM handles internally |

**Challenge Level: LOW** - DG methods are well-matched.

---

## 3. Fault Interface Treatment

### 3.1 Tandem Fault Architecture

This is the **most complex** and **most critical** component to port.

**Key Classes:**
- `AbstractFrictionOperator` (app/form/AbstractFrictionOperator.h)
- `FrictionOperator<LocalOp>` (app/form/FrictionOperator.h)
- `AdapterOperator` (app/form/AdapterOperator.h)
- `RateAndStateBase` (app/localoperator/RateAndStateBase.h)
- `DieterichRuinaAgeing` (app/localoperator/DieterichRuinaAgeing.h)

**Fault Data Model:**
```cpp
// Tandem stores fault state per quadrature point
struct FaultState {
    // Per (D-1) dimensional fault element:
    real_t slip[TangentialComponents];  // D-1 components
    real_t psi;                         // State variable
    real_t V[TangentialComponents];     // Slip rate
    real_t tau[TangentialComponents];   // Shear stress
    real_t sn;                          // Normal stress
};
// Block size = numBasisFunctions * (TangentialComponents + 1)
```

**Domain-Fault Coupling via AdapterOperator:**
```cpp
// Traction computation: domain → fault
traction_at_fault = -stress_tensor * fault_normal;

// Slip velocity: fault → domain boundary condition
domain_bc = slip_rate / 2;  // Split for two-sided fault
```

### 3.2 MFEM Fault Treatment Options

MFEM does **not** have built-in fault handling. Options:

**Option A: Interior Face Integrator Approach**
```cpp
// Create custom integrator for fault faces
class RateStateFrictionIntegrator : public BilinearFormIntegrator {
    void AssembleFaceMatrix(const FiniteElement &el1,
                           const FiniteElement &el2,
                           FaceElementTransformations &Trans,
                           DenseMatrix &elmat) override;
    // Store friction state, compute traction, apply friction law
};
```

**Option B: Subdomain Interface Approach**
- Treat fault as internal boundary between subdomains
- Use `SubMesh` to extract fault surface
- Apply Neumann-type conditions with friction coupling

**Option C: Hybrid (Recommended)**
- Mark fault faces with special boundary attribute
- Use `AddInteriorFaceIntegrator` with fault-specific integrator
- Maintain separate fault `GridFunction` for state variables

### 3.3 Porting Challenge: Fault Interface

| Component | Porting Difficulty | Notes |
|-----------|-------------------|-------|
| Fault geometry tracking | Medium | Need to identify fault faces, store topology |
| Friction state storage | Medium | Separate vector, synced with domain |
| Traction computation | Medium | Face normal derivatives available |
| Slip rate ↔ BC coupling | High | Requires custom operator coupling |
| Rate-and-state law | Low | Pure math, directly portable |

**Challenge Level: HIGH** - Requires significant custom development.

---

## 4. Friction Law Implementation

### 4.1 Tandem Rate-and-State Friction

**Dieterich-Ruina Ageing Law** (app/localoperator/DieterichRuinaAgeing.h):

$$\tau = \sigma_n \cdot a \cdot \sinh^{-1}\left(\frac{V}{2V_0} \exp\left(\frac{f_0 + b\ln(V_0/V_{init})}{a}\right)\right)$$

$$\frac{d\psi}{dt} = -\frac{V}{L} \cdot \sinh^{-1}\left(\frac{V}{2V_0} \exp\left(\frac{\psi}{a}\right)\right)$$

**Key Parameters:**
| Parameter | Description | Typical Value (BP1) |
|-----------|-------------|---------------------|
| a | Direct effect coefficient | 0.010 (nucleation) to 0.025 (creep) |
| b | State evolution coefficient | 0.015 |
| L | Characteristic slip distance | 0.008 m |
| V0 | Reference velocity | 1e-6 m/s |
| f0 | Reference friction coefficient | 0.6 |
| η | Radiation damping | √(μρ)/2 |

### 4.2 MFEM Implementation Strategy

```cpp
// Proposed class structure
class DieterichRuinaAgeingLaw {
public:
    // Physics parameters (spatially variable)
    Coefficient *a_coeff, *b_coeff, *L_coeff;
    Coefficient *eta_coeff, *sn_pre_coeff, *tau_pre_coeff;

    // Constant parameters
    real_t V0, f0;

    // State evolution ODE
    void ComputeStateRHS(const Vector &state, const Vector &slip_rate,
                         Vector &dstate_dt);

    // Slip rate solver (Newton-Raphson)
    void ComputeSlipRate(const Vector &traction, const Vector &state,
                         Vector &slip_rate);

    // Effective traction with radiation damping
    void ComputeEffectiveTraction(const Vector &applied_traction,
                                  const Vector &slip_rate,
                                  Vector &effective_traction);
};
```

**Challenge Level: LOW** - Mathematical implementation is straightforward.

---

## 5. Time Integration Schemes

### 5.1 Tandem Time Integration

**PETSc TS Wrapper** (app/common/PetscTimeSolver.h):
- RK4, RK45 (adaptive), LSSP, implicit methods
- Supports multiple state vectors (1-3)
- Adaptive stepping with error control

**Quasi-Dynamic Mode:**
- Fault state evolves: dψ/dt = f(ψ, V)
- Domain: static equilibrium ∇·σ = 0 at each time step
- Implicit coupling: solve domain, update traction, update fault

**Fully-Dynamic Mode:**
- Three state vectors: velocity v, displacement u, fault state s
- dv/dt = (1/ρ)∇·σ(u) + f
- du/dt = v
- ds/dt = friction_rhs(s, τ(u))
- Explicit time stepping with CFL constraint

### 5.2 MFEM Time Integration

**Available Solvers** (linalg/ode.hpp):
- Explicit: ForwardEuler, RK2, RK3SSP, RK4, RK6, RK8
- Implicit: BackwardEuler, ImplicitMidpoint, SDIRK23/33/34
- IMEX: Additive RK methods

**TimeDependentOperator Interface:**
```cpp
class EarthquakeOperator : public TimeDependentOperator {
public:
    // For quasi-dynamic: only fault state evolves
    void Mult(const Vector &x, Vector &dxdt) const override;

    // For fully-dynamic: implicit solves
    void ImplicitSolve(real_t gamma, const Vector &x, Vector &k) override;

    // Internal state management
    void UpdateDomainSolution();  // Solve elasticity
    void UpdateTraction();        // Compute stress at fault
};
```

### 5.3 Porting Challenge: Time Integration

| Mode | Tandem | MFEM | Gap |
|------|--------|------|-----|
| Quasi-Dynamic | PETSc TS + static solve | ODESolver + BilinearForm::FormLinearSystem | Need coupling logic |
| Fully-Dynamic | PETSc TS explicit | RK4Solver | Direct mapping |
| Adaptive stepping | PETSc TS adaptive | Manual or custom | Needs implementation |

**Challenge Level: MEDIUM** - Coupling logic requires careful design.

---

## 6. Configuration and Input System

### 6.1 Tandem Configuration

**TOML Files** (bp1.toml):
```toml
final_time = 94608000000
mesh_file = "bp1.msh"
mode = "QD"
type = "poisson"
lib = "bp1.lua"
scenario = "bp1"
ref_normal = [-1, 0]
```

**Lua Scenario Scripts** (bp1.lua):
```lua
function BP1:a(x, y)
    local z = -y
    if z < self.H then return self.a0
    elseif z < self.H + self.h then
        return self.a0 + (self.amax - self.a0) * (z - self.H) / self.h
    else return self.amax end
end
```

### 6.2 MFEM Configuration Options

**Option A: Command-line arguments** (current MFEM style)
```cpp
OptionsParser args(argc, argv);
args.AddOption(&order, "-o", "--order", "Polynomial order");
args.AddOption(&dt, "-dt", "--time-step", "Time step size");
```

**Option B: JSON/TOML configuration file** (modern approach)
```json
{
    "mesh_file": "bp2.msh",
    "mode": "quasi_dynamic",
    "time": { "final": 3.15e10, "dt_init": 1.0 },
    "friction": {
        "type": "dieterich_ruina_ageing",
        "a_file": "a_distribution.gf",
        "b": 0.015, "L": 0.008
    }
}
```

**Option C: C++ Coefficient functions** (MFEM native)
```cpp
class RateStateFrictionA : public Coefficient {
    real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) {
        Vector x; T.Transform(ip, x);
        real_t z = -x(1);
        if (z < H) return a0;
        else if (z < H + h) return a0 + (amax - a0) * (z - H) / h;
        else return amax;
    }
};
```

**Challenge Level: LOW** - MFEM's Coefficient system is flexible.

---

## 7. Mesh and Geometry Handling

### 7.1 Mesh Format Differences

| Feature | Tandem | MFEM |
|---------|--------|------|
| Native format | GMSH (.msh) | MFEM mesh, GMSH, VTK, etc. |
| Mesh class | `GlobalSimplexMesh`, `LocalSimplexMesh` | `Mesh`, `ParMesh` |
| Face topology | `DGOperatorTopo`, `LocalFaces` | Built into `Mesh`, `FaceElementTransformations` |
| Curvilinear | Custom `Curvilinear` class | `IsoparametricTransformation` |

### 7.2 Fault Surface Identification

**Tandem:**
- Boundary conditions stored with BC enum (None, Natural, Fault, Dirichlet)
- `BoundaryMap` tracks fault facets
- Explicit fault surface mesh representation

**MFEM:**
```cpp
// Mark fault faces via boundary attributes
Array<int> fault_attr(mesh.bdr_attributes.Max());
fault_attr = 0;
fault_attr[FAULT_BOUNDARY_ID - 1] = 1;

// Or identify interior faces programmatically
for (int f = 0; f < mesh.GetNFaces(); f++) {
    int e1, e2;
    mesh.GetFaceElements(f, &e1, &e2);
    if (IsFaultFace(f)) {
        fault_faces.Append(f);
    }
}
```

**Challenge Level: LOW-MEDIUM** - MFEM mesh handling is comprehensive.

---

## 8. Parallel Implementation

### 8.1 Parallelism Comparison

| Aspect | Tandem | MFEM |
|--------|--------|------|
| Domain decomposition | METIS/ParMETIS | METIS (via ParMesh) |
| Communication | Custom Scatter, CommPattern | hypre IJ interface, HypreParVector |
| Parallel vectors | SparseBlockVector | HypreParVector, ParGridFunction |
| Parallel assembly | Custom ghost exchange | ParBilinearForm, parallel assembly |

### 8.2 MFEM Parallel Infrastructure

```cpp
// Parallel mesh and spaces
ParMesh pmesh(MPI_COMM_WORLD, mesh);
ParFiniteElementSpace pfes(&pmesh, &fec);

// Parallel forms
ParBilinearForm a(&pfes);
a.AddDomainIntegrator(new ElasticityIntegrator(lambda, mu));
a.AddInteriorFaceIntegrator(
    new DGElasticityIntegrator(lambda, mu, alpha, kappa));
a.Assemble();

// Parallel solve
HypreBoomerAMG prec;
CGSolver cg(MPI_COMM_WORLD);
cg.SetPreconditioner(prec);
cg.Mult(b, x);
```

**Challenge Level: LOW** - MFEM parallel infrastructure is mature.

---

## 9. Specific BP1/BP2 Benchmark Requirements

### 9.1 BP1 Benchmark (2D, Planar Fault)

**Problem Setup:**
- 2D anti-plane shear (scalar displacement)
- Vertical strike-slip fault at x=0
- Rate-and-state friction with depth-dependent a
- Quasi-dynamic approximation

**Tandem Implementation:**
- Uses Poisson equation (simplified 2D)
- Fault at domain boundary (not interior)
- Lua functions define spatial variation

**MFEM Implementation Needs:**
1. Scalar L2 DG space
2. Fault on boundary (simpler than interior)
3. Rate-and-state Neumann-type BC
4. Time integration coupling

### 9.2 BP2 Benchmark (3D, Planar Fault)

**Problem Setup:**
- 3D elasticity with vertical strike-slip fault
- Dipping fault geometry variations
- Full stress tensor treatment

**Additional Complexity:**
- Vector elasticity (3 components)
- Tangential slip in 2 directions
- More expensive computations

### 9.3 Implementation Priority

| Feature | BP1 | BP2 | Priority |
|---------|-----|-----|----------|
| Scalar DG | Required | - | High |
| Vector DG elasticity | - | Required | High |
| Rate-state friction | Required | Required | High |
| Fault boundary | Simple | Simple | High |
| Interior fault | Optional | Optional | Medium |
| 3D geometry | - | Required | Medium |
| Adaptive time step | Recommended | Required | Medium |

---

## 10. Implementation Roadmap

### Phase 1: Foundation (Weeks 1-2)

1. **Set up miniapp structure**
   ```
   miniapps/seas/
   ├── CMakeLists.txt
   ├── seas.cpp (main driver)
   ├── friction/
   │   ├── rate_state_friction.hpp
   │   └── dieterich_ruina.hpp
   └── operators/
       └── fault_interface_integrator.hpp
   ```

2. **Implement basic friction law**
   - Port DieterichRuinaAgeing mathematics
   - Create Coefficient-based parameter variation

3. **Implement fault boundary integrator**
   - For boundary faults (BP1 case)
   - Traction calculation
   - Slip rate boundary condition

### Phase 2: Quasi-Dynamic Solver (Weeks 3-4)

1. **Create EarthquakeOperator class**
   - Inherit from TimeDependentOperator
   - Manage domain-fault coupling

2. **Implement time stepping loop**
   - State evolution ODE
   - Static domain solve at each step
   - Traction update

3. **Verify against BP1/BP2 reference data**

### Phase 3: Fully-Dynamic Mode (Weeks 5-6)

1. **Add elastodynamic formulation**
   - Three-field state (u, v, s)
   - CFL time step calculation

2. **Optimize performance**
   - Partial assembly for DG
   - Matrix-free application if needed

### Phase 4: Advanced Features (Weeks 7-8)

1. **Interior fault treatment**
2. **Adaptive time stepping**
3. **Parallel scalability testing**
4. **Documentation and examples**

---

## 11. Risk Assessment

### High Risk Items

1. **Domain-fault coupling complexity**
   - Tandem's AdapterOperator is sophisticated
   - May require significant iteration to get right
   - Mitigation: Start with simpler boundary fault case

2. **Debugging physics**
   - Earthquake simulations are sensitive to parameters
   - Small errors can cause qualitatively wrong behavior
   - Mitigation: Extensive verification against benchmarks

### Medium Risk Items

3. **Performance parity**
   - Tandem uses compile-time optimization (kernel generation)
   - MFEM runtime approach may be slower
   - Mitigation: Use partial assembly, MFEM GPU support

4. **Configuration flexibility**
   - Tandem's Lua scripting is powerful
   - MFEM Coefficient system is less flexible
   - Mitigation: Support JSON/TOML config files

### Low Risk Items

5. **DG method implementation** - Well-matched between codes
6. **Time integration** - MFEM has comprehensive ODE solvers
7. **Parallel support** - MFEM ParMesh is mature

---

## 12. Recommendations

### Immediate Actions

1. **Start with BP1** (simpler 2D case)
2. **Use boundary fault** (avoid interior face complexity initially)
3. **Leverage ex17.cpp** as elasticity template
4. **Use PFF miniapp** as staggered solver template

### Architecture Decisions

1. **Friction as Coefficient** - Implement rate-state parameters as MFEM Coefficients
2. **Custom LinearFormIntegrator** - For fault traction contribution
3. **Separate fault state GridFunction** - Independent of domain solution
4. **JSON configuration** - Modern, readable configuration format

### Testing Strategy

1. **Unit tests** for friction law (compare to analytical)
2. **MMS verification** for DG elasticity
3. **BP1/BP2 validation** against reference data
4. **Scalability tests** for parallel performance

---

## 13. Conclusion

Porting Tandem's earthquake simulation capabilities to MFEM is **feasible but substantial**. The main challenges are:

1. **Fault interface treatment** (HIGH) - Requires custom integrators and coupling
2. **Domain-fault coupling** (HIGH) - Complex operator interaction
3. **Configuration system** (LOW) - MFEM is flexible
4. **DG methods** (LOW) - Good alignment between codes
5. **Time integration** (MEDIUM) - Logic porting, not algorithm

The effort is justified because:
- MFEM provides broader deployment options (GPU, diverse architectures)
- Better integration with LLNL/DOE ecosystem
- Active community and maintenance
- Existing DG infrastructure is solid

Estimated effort: **6-8 weeks** for core functionality, assuming familiarity with both codebases.

---

## Appendix A: Key File Mappings

| Tandem File | MFEM Equivalent | Port Required |
|-------------|-----------------|---------------|
| `src/form/DGOperator.h` | `BilinearForm` + integrators | No (use existing) |
| `app/localoperator/Elasticity.h` | `DGElasticityIntegrator` | No (use existing) |
| `app/localoperator/DieterichRuinaAgeing.h` | (new) `DieterichRuinaFriction` | Yes |
| `app/form/FrictionOperator.h` | (new) `FaultOperator` | Yes |
| `app/form/AdapterOperator.h` | (new) `DomainFaultCoupling` | Yes |
| `app/form/SeasQDOperator.h` | (new) `QuasiDynamicOperator` | Yes |
| `app/form/SeasFDOperator.h` | (new) `FullyDynamicOperator` | Yes |
| `app/common/PetscTimeSolver.h` | `ODESolver` hierarchy | No (use existing) |

## Appendix B: BP1 Parameter Reference

```cpp
// Physical parameters for BP1 benchmark
const real_t H = 15.0e3;        // Nucleation zone depth (m)
const real_t h = 3.0e3;         // Transition zone width (m)
const real_t a0 = 0.010;        // a in nucleation zone
const real_t amax = 0.025;      // a in creeping zone
const real_t b = 0.015;         // Rate-state b parameter
const real_t L = 0.008;         // Characteristic length (m)
const real_t V0 = 1.0e-6;       // Reference velocity (m/s)
const real_t f0 = 0.6;          // Reference friction
const real_t Vp = 1.0e-9;       // Plate velocity (m/s)
const real_t mu = 32.038e9;     // Shear modulus (Pa)
const real_t rho = 2670.0;      // Density (kg/m³)
const real_t eta = sqrt(mu*rho)/2.0;  // Radiation damping
const real_t sn_pre = 50.0e6;   // Normal pre-stress (Pa)
```

## Appendix C: Proposed Class Hierarchy

```
TimeDependentOperator
└── SEASOperator (base)
    ├── QuasiDynamicOperator
    │   ├── DGBilinearForm (elasticity)
    │   ├── FaultOperator
    │   │   └── RateStateFrictionLaw
    │   └── DomainFaultCoupling
    └── FullyDynamicOperator
        ├── DGBilinearForm (elastodynamics)
        ├── FaultOperator
        │   └── RateStateFrictionLaw
        └── DomainFaultCoupling

RateStateFrictionLaw (interface)
├── DieterichRuinaAgeing
├── DieterichRuinaSlip
└── SlipWeakening

DomainFaultCoupling
├── ComputeTraction() - domain → fault
└── ApplySlipRate() - fault → domain BC
```

---

*Report generated: February 5, 2026*
*Author: Analysis of Tandem (https://github.com/TEAR-ERC/tandem) and MFEM (https://mfem.org)*
