# BP5 (3D Full Elasticity) Implementation Plan

**Prepared:** 2026-03-01
**Scope:** Extension of the 2D antiplane MFEM SEAS miniapp to 3D full elasticity (BP5-QD benchmark).
**Prerequisite:** Phase 1a/1b (fault physical group detection, interface cleanup) — completed.

---

## Context

The SEAS miniapp currently supports BP1/BP2 (2D antiplane shear, scalar slip). BP5 is a 3D full elasticity benchmark with vector slip on a 2D fault plane. The recently completed Phase 1a/1b work added:
- `SEASBoundaryTags` + `FaultBoundaryData` for tag-based fault detection
- `DomainOperator` interface methods: `NumSlipComponents()`, `GetFaultCoords2D()`, `GetOffFaultDisplacement()`
- Template-generalized `SEASQuasiDynamicOperator<MeshType, DomainOpType>` and `SEASBdrLoadOperator<MeshType, DomainOpType>`

This plan extends the codebase to support BP5-QD by adding 3D elasticity, vector slip, 2D fault parameters, and BP5-specific output.

### Key BP5 differences from BP1/BP2

| Aspect | BP1/BP2 | BP5 |
|--------|---------|-----|
| Dimension | 2D | 3D |
| PDE | Laplace (antiplane) | Linear elasticity div(sigma)=0 |
| Displacement | Scalar u3 | Vector (u1, u2, u3) |
| Slip | Scalar s | Vector (s2, s3), no-opening [[u1]]=0 |
| Friction parameter a | 1D: a(z) | 2D: a(x2, x3) |
| Slip distance L | Constant | Spatially varying (0.14m / 0.13m nucleation) |
| Normal stress | Constant sigma_n | Constant (homogeneous material) |
| Radiation damping | eta*V scalar | eta*V vector (same scalar eta, applied per-component) |

### Tandem reference conventions (from bp5.lua)

- Tandem coords: x=along-strike, y=fault-normal, z=depth (upward positive)
- SCEC spec: x1=fault-normal, x2=along-strike, x3=depth (positive downward)
- `V_init` in nucleation zone: `(V_zero, 0.01)` (**not** 0.03 as in some docs)
- `tau_pre`: uses `Vi2` (second component) in asinh, returns `(-tau0*Vi1/Vi, -tau0*Vi2/Vi)`
- `boundary`: `u = (Vh, 0, 0)` where `Vh = +/-Vp*t/2` on +/-y walls
- `ref_normal = [0, -1, 0]` (fault normal in y direction)

---

## Phase 1: BP5 Parameters and Vector Friction

### Step 1.1: Create `config/bp5_params.hpp`

Model after `config/bp2_params.hpp`. Contains:

```cpp
struct BP5Params {
   // Material
   real_t rho = 2670.0;           // kg/m^3
   real_t cs = 3464.0;            // m/s
   real_t nu = 0.25;              // Poisson's ratio
   real_t mu;                     // = rho * cs^2 = 32.04 GPa
   real_t lambda;                 // = 2*nu*mu/(1-2*nu) = mu (for nu=0.25)

   // Fault geometry
   real_t Wf = 40.0e3;           // RSF fault width (depth) [m]
   real_t lf = 100.0e3;          // fault length (strike) [m]

   // Friction parameters
   real_t sigma_n = 25.0e6;      // normal stress [Pa]
   real_t a0 = 0.004;            // VW zone a
   real_t amax = 0.04;           // VS zone a
   real_t b0 = 0.03;             // constant b
   real_t L0 = 0.14;             // critical slip distance [m]
   real_t L_nuc = 0.13;          // nucleation zone L [m]
   real_t V0 = 1.0e-6;           // reference velocity
   real_t f0 = 0.6;              // reference friction
   real_t Vp = 1.0e-9;           // plate rate [m/s]
   real_t Vinit = 1.0e-9;        // initial velocity [m/s]

   // Geometry parameters for a(x2,x3)
   real_t hs = 2.0e3, ht = 2.0e3, H = 12.0e3;
   real_t l_vw = 60.0e3;         // VW length in strike
   real_t w_nuc = 12.0e3;        // nucleation zone width

   // 2D parameter functions
   real_t a_of_x2_x3(real_t x2, real_t x3) const;
   real_t L_of_x2_x3(real_t x2, real_t x3) const;
   void V_init_vec(real_t x2, real_t x3, real_t V[2]) const;
   void tau0_vec(real_t x2, real_t x3, real_t eta, real_t tau[2]) const;
   bool IsNucleationZone(real_t x2, real_t x3) const;
};
```

**2D `a(x2, x3)` function** (from BP5 spec Eq. 14):
```
a = a0,    if (hs + ht <= x3 <= hs + ht + H) AND (|x2| <= l/2)           [VW core]
a = amax,  if (x3 <= hs) OR (x3 >= hs + 2*ht + H) OR (|x2| >= l/2 + ht) [VS zones]
a = a0 + r*(amax - a0),  otherwise                                        [transition]
    where r = max(|x3 - hs - ht - H/2| - H/2, |x2| - l/2) / ht
```

**Nucleation zone**: `hs + ht <= x3 <= hs + ht + H` and `-l/2 <= x2 <= -l/2 + w`. In that zone `L = 0.13 m` and `V_init = (V_zero, 0.01)` (following Tandem's bp5.lua).

Reference: Tandem `bp5.lua` functions `a()`, `L()`, `Vinit()`, `tau_pre()`.

### Step 1.2: Add `SolveSlipRateVectorPsi()` to `friction/dieterich_ruina.hpp`

```cpp
void SolveSlipRateVectorPsi(const real_t tau_vec[2], real_t sigma_n,
                             real_t eta, real_t a, real_t psi,
                             real_t V_vec[2]) const
{
    real_t tau_abs = std::sqrt(tau_vec[0]*tau_vec[0] + tau_vec[1]*tau_vec[1]);
    real_t V_abs = SolveSlipRatePsi(tau_abs, psi, sigma_n, eta, a);
    if (tau_abs > 1e-30) {
        V_vec[0] = -(V_abs / tau_abs) * tau_vec[0];
        V_vec[1] = -(V_abs / tau_abs) * tau_vec[1];
    } else {
        V_vec[0] = 0.0;
        V_vec[1] = V_abs;
    }
}
```

This follows Tandem's `DieterichRuinaAgeing::slip_rate()`: solve scalar equation for `||V||`, then apply direction from `tau_hat`.

### Step 1.3: Unit tests for BP5 params and vector friction

- `tests/unit/test_bp5_params.cpp`: Verify `a(x2,x3)` transitions, `L()` nucleation zone, `tau0_vec` direction and magnitude
- `tests/unit/test_vector_friction.cpp`: Verify `SolveSlipRateVectorPsi` for known inputs, verify direction is anti-parallel to tau

---

## Phase 2: Generalize Fault Operator for Vector Slip

### Step 2.1: Add `SlipComponents` template parameter to `RateStateFaultOperator`

File: `fault/rate_state_fault.hpp`

```cpp
template <typename MeshType = Mesh, int SlipComponents = 1>
class RateStateFaultOperator
{
public:
   static constexpr int NumSlipComp = SlipComponents;
   static constexpr int StatePerNode = SlipComponents + 1; // slip(s) + psi
   // For SlipComponents=1: [slip, psi] per node (current layout)
   // For SlipComponents=2: [s2, s3, psi] per node
```

Changes required:
- `StatePerNode` becomes `SlipComponents + 1`
- `PreInit()`: initialize `s2=0, s3=0, psi=psi_ss`
- `Init()`: compute `tau_abs = ||(tau2+tau2_pre, tau3+tau3_pre)||`, solve scalar for V_abs
- `ComputeRHS()`: call `SolveSlipRateVectorPsi()` for `SlipComponents==2`, get `(V2,V3)`, set rates
- `GetSlip()`: extract both slip components for domain operator
- `GetSlipRate()`/`GetSlipRateMax()`: compute `||V||` from components

All changes must be backwards-compatible: `SlipComponents=1` (default) preserves existing scalar behavior identically.

### Step 2.2: Extend `FaultGeometry` for 2D coordinates

File: `fault/fault_geometry.hpp`

Add members:
- `coords_x2_`: along-strike coordinate per fault DOF
- `L_values_`: spatially varying critical slip distance
- `tau_pre_values_`: pre-stress vector (2 components per DOF for BP5)
- `V_init_values_`: initial velocity vector (2 components per DOF for BP5)

The `GatherToRoot`/`DeduplicateByDepth` logic generalizes to 2D deduplication by `(x2, x3)` distance.

### Step 2.3: Unit test for vector fault operator

- `tests/unit/test_vector_fault.cpp`: Create `RateStateFaultOperator<Mesh, 2>`, verify state layout, verify `ComputeRHS` with known 2-component traction

---

## Phase 3: 3D Mesh and Elasticity Domain Operator

This is the largest phase. It creates the 3D DG elasticity solver.

### Step 3.1: Create `bp5/mesh/bp5.geo`

Follow Tandem's `bp5.geo`:
- Two boxes split at fault plane (x1=0), with `BooleanFragments` to embed fault surface
- Nested rectangles for VW/VS/nucleation sub-regions on fault surface (for mesh refinement)
- Physical Surface tags following `SEASBoundaryTags` convention (tag 5 = fault)
- Domain: x1 in [-200, 200] km, x2 in [-100, 100] km, x3 in [0, 100] km
- Suggested element size: Delta_z = 1000 m

Reference: `/Users/chunhuizhao/projects/tandem/examples/tandem/3d/bp5.geo`

### Step 3.2: Create `domain/elasticity3d_operator.hpp`

Inherits `DomainOperator<MeshType>`. Key design:

```cpp
template <typename MeshType = Mesh>
class ElasticityDomainOperator : public DomainOperator<MeshType>
{
public:
   int NumComponents() const override { return 3; }
   int Dimension() const override { return 3; }
   int NumSlipComponents() const override { return 2; }

   void Solve(real_t time, const Vector &slip_bc, GridFuncType &u) override;
   void ComputeTraction(const GridFuncType &u, const Vector &slip_bc,
                         Vector &traction) override;
   void GetFaultCoords2D(Vector &x2, Vector &x3) const override;
   void GetOffFaultDisplacement(const std::vector<Vector> &pts,
                                 Vector &disp) const override;
```

**Implementation strategy: Use MFEM's built-in DG elasticity integrators (IP method first):**

- Volume term: `ElasticityIntegrator(lambda_coeff, mu_coeff)` — available in `fem/bilininteg.hpp`
- Interior face terms: `DGElasticityIntegrator(lambda_coeff, mu_coeff, alpha, kappa)` — available in `fem/bilininteg.hpp`
- Dirichlet BC: `DGElasticityDirichletLFIntegrator(...)` — available in `fem/lininteg.hpp`
- Free surface: Natural BC (zero contribution)

**FE space**: Vector-valued DG space with 3 components:
```cpp
auto fec = make_unique<DG_FECollection>(order, 3, BasisType::GaussLobatto);
auto fes = make_unique<FiniteElementSpace>(&mesh, fec.get(), 3, Ordering::byNODES);
```

### Step 3.3: Fault face identification for 3D

For 3D, a face is on the fault if:
- `|x1_centroid| < tol` AND `|x2_centroid| <= lf/2` AND `0 <= x3_centroid <= Wf`

Use tag-based detection (`FaultBoundaryData::FindFaultInteriorFaces`) when Gmsh Physical Surface(5) is available, with coordinate-based fallback.

### Step 3.4: Slip RHS assembly for 3D

The slip BC enters as a RHS contribution on fault interior faces. For each fault face:
- Prescribed jump: `[[u]] = (0, s2, s3)` (no-opening + tangential slip)
- Following the antiplane pattern in `antiplane_operator.hpp`, but extended to 3-component vector:
  - Consistency term: `- integral_e {{sigma(phi)}} . n * slip dA`
  - Penalty term: `+ integral_e kappa * phi . slip dA`

Reference pattern: `antiplane_operator.hpp:AssembleSlipRHS()` (lines ~400-500)

### Step 3.5: Traction computation for 3D

Extract `tau = {{sigma(u)}} . n_fault` on each fault face. For fault normal `n = (1,0,0)`:
- `tau_2 = {{sigma_21}}` = along-strike shear
- `tau_3 = {{sigma_31}}` = along-dip shear
- `sigma_n = -{{sigma_11}}` (constant for homogeneous material, verify but don't track dynamically)

Project onto fault DOFs using L2 mass-matrix weighted face integral (same pattern as antiplane traction).

### Step 3.6: Boundary conditions for 3D

The BP5 domain has six external faces. Following Tandem's `bp5.lua`:

- **Free surface (x3 = 0)**: Natural BC (zero traction) — no special treatment in DG
- **Bottom (x3 = Lz)**: Natural BC (zero traction)
- **Strike-parallel walls (x2 = +/-Ly)**: Dirichlet loading `u2 = +/-Vp*t/2`, `u1=0`, `u3=0` — applied weakly via `DGElasticityDirichletLFIntegrator`
- **Fault-normal walls (x1 = +/-Lx)**: Natural BC (zero traction) — push domain boundaries out until results converge

The stiffness matrix is assembled once (constant). Only the boundary RHS linear form (Dirichlet loading) needs updating each time step.

### Step 3.7: MMS verification test for 3D elasticity

Create manufactured solution with known displacement field, verify convergence of DG-IP method. Test traction extraction accuracy.

File: `tests/unit/test_elasticity3d.cpp`

---

## Phase 4: SEAS Coupling and Initial Conditions

### Step 4.1: BP5 SEAS operator instantiation

The generalized `SEASQuasiDynamicOperator` already accepts `DomainOpType`:
```cpp
using BP5Domain = ElasticityDomainOperator<ParMesh>;
using BP5SEAS = SEASQuasiDynamicOperator<ParMesh, BP5Domain>;
```

The `RateStateFaultOperator<MeshType, 2>` requires the SEAS operator to handle the 3-component state layout (`[s2, s3, psi]` per node). The `Mult()` method needs to:
1. Extract 2-component slip from state -> pass as `slip_bc` to domain `Solve()`
2. Receive 2-component traction from `ComputeTraction()`
3. Pass traction to fault `ComputeRHS()` which returns 3-component rate `[V2, V3, dpsi/dt]`

The SEAS operator's `Mult()` method queries `domain_->NumSlipComponents()` at runtime to handle both scalar and vector cases. Since `Mult()` is not in a hot inner loop, runtime dispatch is acceptable and avoids template proliferation.

### Step 4.2: Initialization for BP5

`SetInitialCondition()`:
1. `PreInit()`: set `s2=0, s3=0, psi=psi_ss(V_init)` per node
2. Domain solve with zero slip -> get initial traction `(tau2, tau3)`
3. `Init()`: verify stress balance with `tau_pre` — adjust `psi` if needed
4. Verify `V_max ~= V_init` or `V_nuc` in nucleation zone

### Step 4.3: Boundary condition time-dependence

The Dirichlet loading `u2 = +/-Vp*t/2` on the `+/-x2` walls changes with time. The stiffness matrix is assembled once (constant). The RHS linear form contribution from Dirichlet BC must be updated each time step by re-assembling only the boundary linear form with the new `u_D(t)`.

---

## Phase 5: BP5 Output and Driver

### Step 5.1: Create `io/bp5_benchmark_output.hpp`

Handles all 6 output types required by BP5 spec:
1. **On-fault time series** (10 stations, 7 fields each)
2. **Off-fault time series** (9 stations, 6 fields each)
3. **Global time series** (`global.dat`): V_max, moment rate
4. **Earthquake catalog** (`catalog.dat`)
5. **Slip/stress profiles** (8 files)
6. **Rupture time contour** (`rupture.dat`)

Probe stations from BP5 spec Table 2 (on-fault) and Table 3 (off-fault). Use MFEM `FindPoints()` for element location.

### Step 5.2: Create BP5 driver `bp5/bp5_verification.cpp`

Structure follows `tests/verification/bp1_verification_full.cpp`:
```
1. Load bp5.msh (Gmsh mesh)
2. Construct BP5Params
3. Construct ElasticityDomainOperator
4. Construct FaultGeometry with 2D params
5. Construct RateStateFaultOperator<MeshType, 2>
6. Construct SEASQuasiDynamicOperator<MeshType, ElasticityDomainOperator<MeshType>>
7. SetInitialCondition
8. Adaptive time stepping loop
9. BP5BenchmarkOutput at each step
```

### Step 5.3: Parallel driver `bp5/pbp5_verification.cpp`

Same as serial but with `ParMesh`, MPI context, parallel output. BP5 3D will likely require parallel execution from the start due to mesh size.

---

## Phase 6: Verification

1. **Unit tests pass**: All new tests (BP5 params, vector friction, vector fault, 3D elasticity, 3D MMS)
2. **Existing tests unchanged**: All BP1/BP2 tests continue to pass with zero changes
3. **MMS convergence**: 3D elasticity DG-IP shows expected convergence rate
4. **Static equilibrium**: Constant slip -> analytically verifiable traction
5. **Short simulation**: Run BP5 for ~10 years, verify V_max stays near V_init, nucleation zone shows higher V
6. **Long simulation**: Run to first earthquake, compare with Tandem reference output and other community codes

---

## Files Summary

### New files to create

| File | Purpose |
|------|---------|
| `config/bp5_params.hpp` | BP5 parameters, 2D a(x2,x3), L(x2,x3), tau0_vec |
| `domain/elasticity3d_operator.hpp` | 3D DG elasticity operator |
| `bp5/mesh/bp5.geo` | Gmsh geometry for BP5 3D domain |
| `io/bp5_benchmark_output.hpp` | BP5 output (6 output types) |
| `bp5/bp5_verification.cpp` | Serial BP5 driver |
| `bp5/pbp5_verification.cpp` | Parallel BP5 driver |
| `tests/unit/test_bp5_params.cpp` | BP5 parameter tests |
| `tests/unit/test_vector_friction.cpp` | Vector friction solver tests |
| `tests/unit/test_vector_fault.cpp` | Vector fault operator tests |
| `tests/unit/test_elasticity3d.cpp` | 3D elasticity MMS + traction tests |

### Existing files to modify

| File | Change |
|------|--------|
| `friction/dieterich_ruina.hpp` | Add `SolveSlipRateVectorPsi()` |
| `fault/rate_state_fault.hpp` | Add `SlipComponents` template param, vector state layout |
| `fault/fault_geometry.hpp` | Add 2D coords, L_values, tau_pre_values, V_init_values |
| `solver/seas_operator.hpp` | Update `Mult()` to handle NumSlipComponents > 1 |
| `Makefile` + `CMakeLists.txt` | Add all new targets |

### Key reference files

| File | Why |
|------|-----|
| `domain/antiplane_operator.hpp` | Blueprint for 3D operator (fault detection, slip RHS, traction) |
| Tandem `app/localoperator/Elasticity.h`/`.cpp` | Ground truth for 3D DG elasticity assembly |
| Tandem `examples/tandem/3d/bp5.lua` | Parameter functions, sign conventions |
| Tandem `examples/tandem/3d/bp5.geo` | 3D mesh geometry reference |
| MFEM `fem/bilininteg.hpp` | `DGElasticityIntegrator` API |
| MFEM `fem/lininteg.hpp` | `DGElasticityDirichletLFIntegrator` API |

---

## Key Design Decisions

### D.1: IP method first, BR2 later

Use MFEM's built-in `DGElasticityIntegrator` (IP/SIPG method) for Phase 3. This avoids writing a custom BR2 elasticity integrator and gets a working end-to-end system faster. BR2 can be added later as an option.

### D.2: Runtime dispatch for SlipComponents in SEAS operator

Rather than adding another template parameter to `SEASQuasiDynamicOperator`, query `domain_->NumSlipComponents()` at runtime in `Mult()`. This avoids template proliferation while keeping the hot inner loop (Newton iteration in friction) templated for performance.

### D.3: Constant normal stress

For BP5 with homogeneous material and no-opening condition, `sigma_n = 25 MPa` is constant for all time. The `ElasticityDomainOperator` will verify this but not track it dynamically. This simplifies the coupling significantly.

### D.4: Coordinate convention

Map from Tandem's (x=strike, y=normal, z=depth-up) to SCEC spec (x1=normal, x2=strike, x3=depth-down). All internal code uses the SCEC convention. The `.geo` file and output follow the SCEC convention.

---

## Implementation Order

**Phase 1** (BP5 params + vector friction) -> **Phase 2** (vector fault operator) -> **Phase 3** (3D mesh + elasticity) -> **Phase 4** (SEAS coupling) -> **Phase 5** (output + driver) -> **Phase 6** (verification)

Phases 1 and 2 can be developed and tested independently of Phase 3 (they only require the existing mesh for testing). Phase 3 is the largest and most critical phase. Phases 4 and 5 integrate everything.
