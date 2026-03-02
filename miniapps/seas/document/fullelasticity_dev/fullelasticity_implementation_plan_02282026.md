---

# BP5 SEAS Benchmark: Full Elasticity Implementation Plan

**Prepared:** 2026-02-28
**Scope:** Extension of the 2D antiplane MFEM SEAS miniapp to 3D full elasticity (BP5-QD benchmark). Covers BP5 physics specification, 5-phase implementation plan, and key design decisions.

**Prerequisite:** Before starting BP5 implementation, complete the preparatory refactoring described in `mfem_seas_system_updateplan_02282026.md` (input file infrastructure, fault physical group detection, interface cleanup). BP5 implementation assumes those changes are already in place.

---

## Part A: BP5 Problem Specification Summary

### A.1 Domain and Coordinate System

BP5 is a 3D half-space problem. The physical domain is:

```
(x1, x2, x3) ∈ (-∞, ∞) × (-∞, ∞) × (0, ∞)
```

with `x3` positive downward. The free surface lies at `x3 = 0`. For computation the half-space is truncated: Tandem's `bp5.geo` uses `x1 ∈ [-200, 200] km`, `x2 ∈ [-100, 100] km`, `x3 ∈ [0, 100] km` (all in km). The fault is a vertical plane at `x1 = 0`, extending `l_f = 100 km` in the `x2` direction and `W_f = 40 km` in the `x3` direction.

The displacement is a 3-component vector `u = (u1, u2, u3)`. The equilibrium equation is `div(sigma) = 0` (quasi-static for BP5-QD). Hooke's law reads `sigma_ij = K eps_kk delta_ij + 2 mu (eps_ij - 1/3 eps_kk delta_ij)` with `eps_ij = (du_i/dx_j + du_j/dx_i)/2`. Material parameters: `rho = 2670 kg/m³`, `cs = 3.464 km/s`, `nu = 0.25`, so `mu = rho * cs² ≈ 32.04 GPa` and `lambda = 2 nu mu / (1 - 2 nu) = 2 * 0.25 * mu / 0.5 = mu`.

### A.2 Boundary Conditions

**Free surface** (`x3 = 0`): Zero traction, `sigma_j3 = 0` for `j = 1, 2, 3`.

**Artificial far-field boundaries** (truncated half-space walls at `±x1`, `±x2`, `x3 = Lz`): The spec notes that numerical BCs will affect results and recommends pushing boundaries until results are independent of domain size. Tandem uses Dirichlet loading `u = ±(Vp/2)t` in the `x2`-direction on the `x2`-boundaries (the `boundary` function in `bp5.lua` applies `Vh = ±Vp*t/2` depending on `y > 1` vs `y < -1`), and zero traction (Natural) elsewhere. This is the same approach used for BP1/BP2.

### A.3 Fault Interface Conditions

At `x1 = 0`:
- **No-opening condition**: `u1(0+, x2, x3, t) = u1(0-, x2, x3, t)` — no normal displacement jump.
- **Slip vector**: `s_j(x2, x3, t) = u_j(0+) - u_j(0-)` for `j = 2, 3` (tangential components only).
- **Traction continuity**: The normal-direction traction `sigma_11` is equal and opposite across the fault; the tangential tractions `sigma_21` and `sigma_31` are continuous.
- These traction components define the normal stress `sigma = -sigma_11` (positive in compression) and the shear stress vector `tau = (sigma_21, sigma_31)`.

The fault is split into two zones by the `x2` and `x3` coordinates:

**RSF zone** (`|x2| <= l_f/2` and `0 <= x3 <= W_f`): Rate-and-state friction governs. The shear stress vector `tau = tau_0 + Delta_tau - eta * V` (quasi-dynamic), where `eta = mu/(2*cs)` is the shear-wave impedance (scalar, isotropic radiation damping), `V = (V2, V3)` is the slip velocity vector, and `Delta_tau` is the quasi-static stress transfer. The friction condition is `tau = F(V, theta) * V / ||V||` where `F = sigma_bar_n * f(||V||, theta)` and `f` is the regularized Dieterich-Ruina law.

**Creep zone** (outside RSF zone): Prescribed `V2 = Vp`, `V3 = 0` (plate-rate loading with no vertical component).

### A.4 Friction Law and State Evolution

The friction coefficient uses the regularized Dieterich-Ruina form (identical to BP1/BP2 in form):

```
f(V, theta) = a * asinh[ (V / 2V0) * exp((f0 + b * ln(V0 * theta / L)) / a) ]
```

with aging law `dtheta/dt = 1 - V * theta / L`.

The scalar pre-stress `tau_0` (a scalar applied in the direction `V/||V||`) satisfies:

```
tau_0 = sigma_bar_n * a * asinh[ (V_init / 2V0) * exp((f0 + b*ln(V0/V_init)) / a) ] + eta * V_init
```

The key BP5 parameters (Table 1 from spec):

| Parameter | Value |
|-----------|-------|
| `rho` | 2670 kg/m³ |
| `cs` | 3.464 km/s |
| `nu` | 0.25 |
| `a0` | 0.004 (VW zone) |
| `a_max` | 0.04 (VS zone) |
| `b0` | 0.03 (constant) |
| `sigma_bar_n` | 25 MPa |
| `L` | 0.14 m (0.13 m in nucleation zone) |
| `V_p` | 1e-9 m/s |
| `V_init` | 1e-9 m/s |
| `V_0` | 1e-6 m/s |
| `f_0` | 0.6 |
| `h_s` | 2 km (shallow VS zone width) |
| `h_t` | 2 km (VW-VS transition width) |
| `H` | 12 km (uniform VW region width in depth) |
| `l` | 60 km (length of uniform VW region in strike) |
| `W_f` | 40 km (RSF fault width in depth) |
| `l_f` | 100 km (fault length in strike) |
| `w` | 12 km (nucleation zone width) |
| `t_f` | 1800 years |
| `Delta z` | 1000 m (suggested) |

The 2D `a(x2, x3)` parameter function (Eq. 14 from spec, with `x3` down-positive):

```
a = a0,    if (h_s + h_t <= x3 <= h_s + h_t + H) AND (|x2| <= l/2)           [VW core]
a = a_max, if (x3 <= h_s) OR (x3 >= h_s + 2*h_t + H) OR (|x2| >= l/2 + h_t)  [VS zones]
a = a0 + r*(a_max - a0),  otherwise                                            [transition]
    where r = max(|x3 - h_s - h_t - H/2| - H/2, |x2| - l/2) / h_t
```

The nucleation zone occupies the corner: `h_s + h_t <= x3 <= h_s + h_t + H` and `-l/2 <= x2 <= -l/2 + w`. In that zone `L = 0.13 m`; elsewhere `L = 0.14 m`.

Initial conditions: `s_j(x2, x3, 0) = 0`, `theta(x2, x3, 0) = L / V_init`, uniform `V_init = V_p = 1e-9 m/s` everywhere except the nucleation zone where `V_2_init = 0.03 m/s` (note: Tandem uses `V_zero = 1e-20` for the `x3` component to avoid `log(0)` issues in output).

### A.5 Output Requirements

BP5 requires six output types:

1. **On-fault time series** at 10 stations (2D fault coordinates): fields `slip_2`, `slip_3`, `log10(slip_rate_2)`, `log10(slip_rate_3)`, `shear_stress_2`, `shear_stress_3`, `log10(state)`.

2. **Off-fault time series** at 9 stations (3D spatial positions on free surface and at depth): fields `disp_1`, `disp_2`, `disp_3`, `vel_1`, `vel_2`, `vel_3`.

3. **Source parameter time series** (`global.dat`): `V_max` and moment rate `M_t = integral mu V dA` over the VW patch area.

4. **Earthquake catalog** (`catalog.dat`): per-event `t_start`, `t_end`, `rup_area`, `avg_stress_start`, `avg_stress_end`, `avg_slip`. Event defined by `max(V) >= 1e-3 m/s`.

5. **Slip and stress evolution profiles** (8 files): along-dip at `x2 = 0` and along-strike at `x3 = 10 km`, for both `s_2`, `s_3`, `tau_2`, `tau_3`.

6. **Rupture time contour** (`rupture.dat`): `(x2, x3, t)` at which `V >= 1 mm/s` for first event.

### A.6 Key Differences from BP1/BP2

| Aspect | BP1/BP2 | BP5 |
|--------|---------|-----|
| Spatial dimension | 2D (`x1`, `x3`) | 3D (`x1`, `x2`, `x3`) |
| Displacement | Scalar `u3` (antiplane) | Vector `(u1, u2, u3)` |
| PDE | Laplace `∇²u3 = 0` | Linear elasticity `div(sigma) = 0` |
| Fault | 1D line at `x1=0` | 2D plane at `x1=0` |
| Slip | Scalar `s` | Vector `(s2, s3)` |
| Slip velocity | Scalar `V` | Vector `(V2, V3)`, friction uses `||V||` |
| Normal stress | Constant `sigma_n` | Computed from `sigma_11` (constant for homogeneous material) |
| Radiation damping | `eta = mu/(2cs)` scalar → scalar `tau` | `eta = mu/(2cs)` scalar but applied as `eta * V` vector |
| `a` parameter | 1D: `a(z)` | 2D: `a(x2, x3)` |
| `L` parameter | Constant | Spatially varying (nucleation zone) |
| Domain | 2D rectangle | 3D box with embedded fault surface |
| Output | 1D probe depth | 2D probe `(x2, x3)` coordinates |
| Off-fault output | None required | 9 off-fault stations (3D displacement + velocity) |

---

## Part B: Phased Implementation Plan

### Phase 1: Code Restructuring (Preparatory Refactoring)

**Goal:** Make the existing code ready to host a second domain operator without breaking BP1/BP2.

**Step 1.1: Generalize `SEASQuasiDynamicOperator` domain type.**

File: `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/solver/seas_operator.hpp`

Current line:
```cpp
using DomainOpType = AntiplaneDomainOperator<MeshType>;
```

Change to a second template parameter with default:
```cpp
template <typename MeshType = Mesh,
          typename DomainOpType = AntiplaneDomainOperator<MeshType>>
class SEASQuasiDynamicOperator : public TimeDependentOperator
```

This change is purely additive and backwards-compatible. All existing BP1/BP2 instantiations continue to work without change.

**Step 1.2: Add `GetFaultCoords2D` to `DomainOperator` base class.**

File: `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/domain/domain_operator.hpp`

Add virtual methods:
```cpp
/// Get 2D fault coordinates (x2, x3) at each fault DOF — needed for BP5
/// Default implementation returns (0, depth) using GetFaultDepths.
virtual void GetFaultCoords2D(Vector &coords_x2, Vector &coords_x3) const
{
    Vector depths;
    GetFaultDepths(depths);
    coords_x2.SetSize(depths.Size());
    coords_x2 = 0.0;
    coords_x3 = depths;   // For 2D problems x3 is just the depth
}

/// Get number of slip components per fault DOF (1 for antiplane, 2 for 3D)
virtual int NumSlipComponents() const { return 1; }
```

**Step 1.3: Generalize `RateStateFaultOperator` state layout.**

File: `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/fault/rate_state_fault.hpp`

Add a template parameter `int SlipComponents = 1`. For BP5, `SlipComponents = 2`.

```cpp
template <typename MeshType = Mesh, int SlipComponents = 1>
class RateStateFaultOperator
{
public:
    static constexpr int NumSlipComp = SlipComponents;
    static constexpr int StatePerNode = SlipComponents + 1;  // slip(s) + theta
    static constexpr int ThetaIndex = SlipComponents;
```

All `GetSlip`/`SetSlip` methods become aware of `NumSlipComp` DOFs. The scalar `slip_rate_` vector stays (stores `||V||`), but we add `slip_rate_vec_` storing all `SlipComponents` per node.

**Step 1.4: Add vector friction interface to `FrictionLaw`.**

File: `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/friction/friction_law.hpp`

```cpp
/// Solve for vector slip rate V = [V2, V3] given vector shear stress tau = [tau2, tau3].
/// Uses scalar V = ||V||, then applies direction: Vi = V * tau_hat / ||tau_hat||
/// Default implementation handles scalar case.
virtual void SolveSlipRateVector(const real_t tau_vec[], real_t sigma_n,
                                  real_t eta, real_t a,
                                  real_t theta, real_t V_out[],
                                  int ncomp) const;
```

This follows exactly Tandem's `DieterichRuinaAgeing::slip_rate()`, which:
1. Computes `tauAbs = norm(tau + tau_pre)`.
2. Solves the scalar equation `tauAbs = F(snAbs, V, psi) + eta * V` for scalar `V`.
3. Returns `Vi = -(V / tauAbs) * (tau + tau_pre)` as the vector.

**Step 1.5: Create `BP5Params` struct.**

File: `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/config/bp5_params.hpp`

Model after `bp2_params.hpp`. Key additions:
- 2D `a_of_x2_x3(real_t x2, real_t x3)` replacing `a_of_z(z)`
- `b0` (constant, unlike BP2 where `b` is fixed)
- Spatially varying `L_of_x2_x3(real_t x2, real_t x3)` (two values: 0.14 m and 0.13 m in nucleation zone)
- Nucleation zone geometry: `w`, `l`, `h_s`, `h_t`, `H` parameters
- `nu` (Poisson's ratio, needed for `lambda` computation)
- `tau_pre` as a 2-component vector (scalar `sigma_bar_n * f_ss`, then directed by `V_init_hat`)
- `V_init_vec[2]` (initial slip velocity vector: `(0, V_p)` everywhere except nucleation zone where `(V_zero, 0.03)`)
- `tau0_vec()` method returning a 2-component pre-stress vector

---

### Phase 2: 3D Mesh and BP5 Configuration

**Goal:** Produce a valid MFEM 3D mesh with the fault plane correctly identified as an interior surface, matching the BP5 geometry.

**Step 2.1: Create a BP5 Gmsh geometry file.**

File: `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/bp5/bp5.geo`

Follow Tandem's `bp5.geo` almost exactly. The key geometric elements are:
- Two boxes `Box(1)` and `Box(2)` for `y < 0` and `y > 0` halves of the domain (where `y` = Tandem's `x2`-axis, the along-strike direction). Tandem places `y ∈ [-100, 100]`, `x ∈ [-200, 200]`, `z ∈ [-100, 0]` (in km).
- A `Rectangle` surface for the fault at `x1 = 0` (Tandem's `x = 0`), with extent `[-l_f/2, l_f/2]` in strike and `[0, W_f]` in depth.
- Nested rectangles for the nucleation zone sub-regions (nuc1, nuc2, nuc3 for the transition zones and the nucleation patch).
- `BooleanFragments` to split the volumes and embed the fault surface.
- Physical Surface attributes: `1` = top/bottom (Dirichlet loading BCs); `3` = fault (interior surface); `5` = other exterior walls (Dirichlet loading or Natural BC depending on orientation).

The MFEM coordinate convention for BP5 should follow the spec: `x1` fault-normal (horizontal), `x2` along-strike (horizontal), `x3` depth (positive downward). In Tandem's Lua/geo, they use `x` = along-strike, `y` = fault-normal, `z` = depth-upward-positive, with the fault at `y = 0`. This requires a coordinate relabeling when translating from Tandem's geo to the spec's convention.

**Step 2.2: Create `BP5MeshGenerator` class.**

File: `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/domain/bp5_mesh.hpp`

Similar to `BP2MeshGenerator` but for 3D. Key struct:

```cpp
struct BP5BoundaryAttributes {
    static constexpr int FREE_SURFACE = 1;   // x3 = 0 (top) — Natural BC
    static constexpr int BOTTOM = 2;         // x3 = Lz — Natural or Dirichlet
    static constexpr int STRIKE_POS = 3;     // x2 = +L2 — Dirichlet plate loading
    static constexpr int STRIKE_NEG = 4;     // x2 = -L2 — Dirichlet plate loading
    static constexpr int NORMAL_POS = 5;     // x1 = +L1 — Natural BC
    static constexpr int NORMAL_NEG = 6;     // x1 = -L1 — Natural BC
};
```

Tandem's `bp5.lua` `boundary()` function applies `u2 = ±Vp/2 * t` on the `±x2` walls (strike-parallel faces), zero elsewhere. The plate loading approach for BP5 should follow this.

The `IsFaultFace3D` check: a face is on the fault if its centroid has `|x1| < tol`, `|x2| <= l_f/2`, and `0 <= x3 <= W_f`.

**Step 2.3: Mesh testing.**

Before building the full operator, verify the mesh: confirm correct face identification, correct element counts, correct boundary attribute assignments. A small `2 x 2 x 2` test mesh (two elements on each side of the fault) is sufficient for smoke tests.

---

### Phase 3: 3D Vector DG Elasticity Operator

**Goal:** Implement `ElasticityDomainOperator<MeshType>` satisfying the `DomainOperator` interface, solving the 3D linear elasticity problem with DG-BR2.

**Step 3.1: Understand what the operator must do.**

For BP5, the quasi-static equilibrium equation is:

```
∇ · sigma(u) = 0   in Omega
sigma(u) · n = 0   on free surface (Gamma_N)
u = u_D(t)         on Dirichlet boundaries (Gamma_D: ±x2 walls)
[[u_tangential]] = s   on fault (Gamma_F: x1 = 0)
[[u_normal]] = 0   on fault (no-opening condition)
```

The DG formulation for linear elasticity follows the same structure as for the Laplacian but with the stress tensor replacing the gradient:

```
B_h(u, v) = SUM_K integral_K sigma(u):grad(v) dV
           - SUM_e integral_e {{sigma(u)}} : [[v]] dA      (consistency)
           - SUM_e integral_e {{sigma(v)}} : [[u]] dA      (symmetry, SIPG)
           + SUM_e penalty_e integral_e [[u]] : [[v]] dA   (stability)
```

where `sigma(u) = lambda * div(u) * I + 2*mu * sym(grad(u))` (Voigt form), and `[[u]]` for the fault interface contains only the tangential components (the normal component is zero by the no-opening condition).

This is exactly the elasticity BR2 operator described in Tandem's `Elasticity.h`. The key observation from Tandem's implementation is:
- `NumQuantities = DomainDimension` (3 for 3D): the solution `u` is a 3-vector.
- The block size is `space.numBasisFunctions() * 3`.
- The stiffness matrix has size `(3 * ndof) × (3 * ndof)`.
- The slip BC is applied as a facet functional `fun_slip` that returns a 3-vector (with `u1 = 0` for the no-opening condition, `u2 = s2/2`, `u3 = s3/2` on the positive side).

**Step 3.2: Choose the MFEM implementation strategy.**

There are two options:

**Option A: Use MFEM's existing `ElasticityIntegrator` + `DGElasticityIntegrator`.**  
MFEM's built-in `ElasticityIntegrator` handles the volume term. MFEM does not provide a ready-made BR2 elasticity integrator, but does provide `DGElasticityIntegrator` for IP (SIPG/NIPG). If IP is acceptable (which it is — Tandem also supports IP for elasticity), this reduces implementation work significantly. The penalty parameter for 3D IP elasticity follows `kappa = p^2 * (2*mu + lambda)` scaled by face area/volume ratio.

**Option B: Implement a custom BR2 elasticity integrator.**  
Write `DGElasticityBR2Integrator` following the same lifting-operator pattern as the existing `BR2InteriorFaceIntegrator` in `dg_br2_integrator.hpp`, but generalized from scalar `K*grad(u)` to tensor `sigma(u)`. This is the more accurate and consistent approach (Tandem uses BR2 by default). The main complexity is that the "lifting" step now involves a rank-4 stiffness tensor rather than a scalar coefficient.

**Recommendation:** For Phase 3 start with Option A (IP method) to get a working end-to-end system, then implement Option B (BR2) in Phase 3b. This mirrors the existing antiplane code which also supports both IP and BR2 via the `DGMethod` enum.

**Step 3.3: MFEM vector FE space for 3D elasticity.**

```cpp
// Use vector-valued DG space
auto fec = std::make_unique<DG_FECollection>(order, 3, BasisType::GaussLobatto);
auto fes = std::make_unique<FiniteElementSpace>(&mesh, fec.get(), 3, Ordering::byNODES);
```

The `3` as the third argument to `FiniteElementSpace` creates a vector space with 3 components.

**Step 3.4: `ElasticityDomainOperator` class skeleton.**

File: `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/domain/elasticity3d_operator.hpp`

```cpp
template <typename MeshType = Mesh>
class ElasticityDomainOperator : public DomainOperator<MeshType>
{
public:
    ElasticityDomainOperator(MeshType &mesh, int order,
                              real_t lambda, real_t mu,
                              real_t Vp, real_t Wf_depth, real_t lf_strike,
                              DGMethod method = DGMethod::IP);

    int NumComponents() const override { return 3; }
    int Dimension() const override { return 3; }
    int NumSlipComponents() const override { return 2; }

    // Solve equilibrium: div(sigma(u)) = 0 with slip BC on fault
    // slip_bc has layout [s2_0, s3_0, s2_1, s3_1, ...] (2 components per fault DOF)
    void Solve(real_t time, const Vector &slip_bc, GridFuncType &u) override;

    // Compute traction vector (sigma · n) on fault interior faces
    // Returns [tau2_0, tau3_0, tau2_1, tau3_1, ...] and [sigma_n_0, ...]
    // For BP5 the normal stress is stored separately since friction needs it
    void ComputeTraction(const GridFuncType &u, const Vector &slip_bc,
                          Vector &traction) override;

    // Additional BP5 method: compute full traction (tangential + normal)
    void ComputeFullTraction(const GridFuncType &u, const Vector &slip_bc,
                              Vector &tau_vec, Vector &sigma_n_vec);

    // Off-fault displacement at prescribed points for output
    void GetOffFaultDisplacement(const std::vector<Vector> &points,
                                  Vector &displacements) const override;

    int GetNumFaultDOFs() const override { return num_fault_dofs_; }
    void GetFaultDepths(Vector &depths) const override;
    void GetFaultCoords2D(Vector &coords_x2, Vector &coords_x3) const override;
    const Array<int> &GetFaultDOFs() const override { return fault_dofs_; }

private:
    // Fault identification: (|x1| < tol) && (|x2| < lf/2) && (0 < x3 < Wf)
    bool IsFaultFace3D(int face) const;

    // Slip BC assembly: impose [[u]] = (0, s2, s3) on fault faces
    void AssembleSlipContribution(Vector &rhs, const Vector &slip_bc) const;

    real_t lambda_, mu_, Vp_, Wf_, lf_;
    // ... (similar private state to AntiplaneDomainOperator)
};
```

**Step 3.5: Stiffness matrix assembly.**

The stiffness matrix is assembled using:

```cpp
// Volume term: integral_Omega sigma(u):grad(v) dOmega
a_form_->AddDomainIntegrator(new ElasticityIntegrator(lambda_coeff, mu_coeff));

// Interior face terms (IP method): DGElasticityIntegrator
a_form_->AddInteriorFaceIntegrator(
    new DGElasticityIntegrator(lambda_coeff, mu_coeff, alpha=-1.0, kappa));

// Boundary face terms for Dirichlet BCs (loading walls)
// Free surface gets Natural BC (zero contribution)
// Loading walls at ±x2 get weak Dirichlet via DGElasticityIntegrator on boundary
```

The key difference from the antiplane case is that MFEM's `DGElasticityIntegrator` handles the interior face terms for all three displacement components simultaneously, assembling a `(3*ndof) × (3*ndof)` block matrix per face.

**Step 3.6: Fault slip RHS assembly.**

The slip boundary condition enters as a right-hand side contribution. For each fault interior face the DG residual contribution is:

```
F_slip = - {{sigma(v)}} : [[s]] n  +  penalty * [[v]] : [[s]] n
```

where `[[s]] = (0, s2, s3)` (the no-opening condition sets `[[s1]] = 0`). This requires iterating over fault interior faces and computing the face integral of the test function gradient dotted with the prescribed slip jump.

In practice, this is assembled by setting up a `LinearForm` with a custom face integrator for the fault faces that encodes the slip jump. The key formula, following the antiplane pattern:

```cpp
// For each fault face e:
// Add to global RHS:
//   - integral_e {{sigma(phi_i) n}} . slip ds  (for each test function phi_i)
//   + integral_e penalty * phi_i . slip ds     (for each test function on fault DOFs)
```

**Step 3.7: Traction computation.**

After solving for `u`, the traction vector on the fault is extracted as:

```cpp
// For each fault face e:
// tau_h(x2, x3) = {{sigma(u_h)}} . n_fault
// where n_fault = (1, 0, 0) or (-1, 0, 0) depending on side orientation
//
// Extract tangential components: tau = (sigma_21, sigma_31) from the (1) row
// The normal stress: sigma_n = -sigma_11 (positive in compression)
//
// Since both sides of a homogeneous medium contribute, and the jump is zero
// for the normal component, sigma_n is the same from both sides.
```

The traction is projected onto the fault DOFs using the same DOF-extraction logic as in `AntiplaneDomainOperator::ComputeTraction`.

---

### Phase 4: 3D Fault Geometry and Traction Decomposition

**Goal:** Extend `FaultGeometry` and `RateStateFaultOperator` to handle the 2D fault and vector slip.

**Step 4.1: Extend `FaultGeometry` to 2D.**

The class stores:
- `coords_x2_`, `coords_x3_`: Vector of `(x2, x3)` coordinates for each fault DOF.
- `a_values_`: 2D function `a(x2, x3)` evaluated at each DOF.
- `L_values_`: 2D function `L(x2, x3)` (spatially varying critical slip distance).
- `eta_values_`: Constant `mu/(2*cs)` (same scalar for all DOFs, all components).
- `sn_pre_values_`: Constant effective normal stress (25 MPa for BP5, since homogeneous).
- `tau_pre_values_`: Pre-stress vector `[tau2_pre, tau3_pre]` per DOF — computed from `BP5Params::tau0_vec()`.
- `V_init_values_`: Initial velocity vector `[V2_init, V3_init]` per DOF.

The `GatherToRoot` / `DeduplicateByDepth` logic generalizes to 2D deduplication using a `(x2, x3)` distance tolerance instead of 1D depth tolerance.

**Step 4.2: Extend `RateStateFaultOperator` for vector slip.**

With `SlipComponents = 2`, the state layout per fault DOF is `[s2, s3, psi]` (using Tandem's psi-space formulation for robustness).

The `PreInit` method initializes `s2 = 0`, `s3 = 0`, `psi = psi_ss` (steady-state at `V_init`).

The `Init` method:
1. Receives `traction` vector of size `2 * num_nodes` (tau2 and tau3 components).
2. For each node: computes `tau_abs = ||(tau2 + tau2_pre, tau3 + tau3_pre)||`.
3. Solves the scalar equation `tau_abs = sigma_n * F(V_abs, psi) + eta * V_abs` for `V_abs`.
4. Computes `psi_0` from initial condition.
5. Returns `V_max = max(V_abs)`.

The `ComputeRHS` method:
1. Receives `traction` vector of size `2 * num_nodes`.
2. For each node: calls `SolveSlipRateVector` to get `(V2, V3)`.
3. Sets `rate[s2_index] = V2`, `rate[s3_index] = V3`.
4. Sets `rate[psi_index] = state_evolution(V_abs, psi)`.

The `GetSlip` method now extracts both components: returns a vector of size `2 * num_nodes` with layout `[s2_0, s3_0, s2_1, s3_1, ...]`.

**Step 4.3: Vector friction solver.**

The `DieterichRuinaFriction` class needs a new method:

```cpp
// Solve for vector slip rate given vector shear stress
// Input: tau_vec[2] = (tau2, tau3) including pre-stress
//        sigma_n (positive in compression)
//        eta, a, psi
// Output: V_vec[2] = (V2, V3)
void SolveSlipRateVector(const real_t tau_vec[], real_t sigma_n,
                          real_t eta, real_t a, real_t psi,
                          real_t V_vec[]) const
{
    real_t tau_abs = std::sqrt(tau_vec[0]*tau_vec[0] + tau_vec[1]*tau_vec[1]);
    real_t V_abs = SolveSlipRatePsi(tau_abs, psi, sigma_n, eta, a);
    // Direction: V is anti-parallel to tau (fault motion opposes resistance)
    // Following Tandem: Vi = -(V/tauAbs) * tau
    if (tau_abs > 1e-30) {
        V_vec[0] = -(V_abs / tau_abs) * tau_vec[0];
        V_vec[1] = -(V_abs / tau_abs) * tau_vec[1];
    } else {
        V_vec[0] = 0.0;
        V_vec[1] = V_abs;  // Default to strike direction
    }
}
```

The sign convention requires careful attention. In Tandem, the `tau_pre` is included in the input `tau_vec` before computing `tauAbs`. The returned slip rate opposes the applied shear stress direction.

**Step 4.4: Normal stress handling.**

For BP5, since the material is homogeneous and there is no normal displacement jump, `sigma_n` is constant (25 MPa). The `ElasticityDomainOperator::ComputeFullTraction` will also return `sigma_n_vec` (one value per fault DOF), but for homogeneous material they will all equal the pre-stress value plus the quasi-static contribution (which is zero since slip on a fault between identical materials does not change normal stress). This is explicitly stated in the BP5 spec: "Since slip on a fault separating identical materials does not alter the normal traction, sigma_n remains constant." Therefore `sigma_n = sigma_bar_n = 25 MPa` everywhere on the fault for all time.

**Step 4.5: Radiation damping in 3D.**

For BP5-QD, the radiation damping term is `eta * V` where `eta = mu/(2*cs)` is a scalar (isotropic). The vector form of the quasi-dynamic constitutive relation is:

```
tau = tau_0 + Delta_tau - eta * V
```

where `tau`, `V` are 2-vectors. This means the radiation damping is applied component-wise with the same scalar `eta`. This is identical to Tandem's implementation (see `DieterichRuinaAgeing::Params::eta` and `tau_hat`). No coupling between strike and dip components through the damping term exists.

For the quasi-dynamic radiation damping in 3D, the correct treatment (as confirmed by BP5 spec Eq. 10) is that `eta = mu/(2cs)` is half the shear-wave impedance. This is a scalar that multiplies the full 2-component velocity vector independently per component. No P-wave coupling term appears in the quasi-dynamic approximation because inertia is entirely neglected.

---

### Phase 5: Integration and Testing

**Step 5.1: Create BP5 SEAS operator.**

File: `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/solver/seas_operator.hpp` (generalized) or a new typed alias in `bp5/bp5_verification.cpp`:

```cpp
// Typed alias for BP5
using BP5DomainOp = ElasticityDomainOperator<SEASMesh>;
using BP5FaultOp  = RateStateFaultOperator<SEASMesh, 2>;  // 2 slip components
using BP5SEASOp   = SEASQuasiDynamicOperator<SEASMesh, BP5DomainOp>;
```

The `SEASQuasiDynamicOperator` is already agnostic about the domain operator type once the second template parameter is introduced. No further changes are needed.

**Step 5.2: Create `BP5BenchmarkOutput` class.**

File: `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/io/bp5_benchmark_output.hpp`

Key responsibilities:
- Maintain 10 on-fault probe locations as `(x2, x3)` pairs.
- Maintain 9 off-fault probe locations as `(x1, x2, x3)` triples.
- Write per-step: on-fault time series (8 fields), off-fault time series (7 fields).
- Track `V_max` and moment rate for `global.dat`.
- Detect seismic events and populate `catalog.dat`.
- Maintain per-time-step slip and stress profiles for the 8 slip/stress evolution files.
- Write `rupture.dat` for first event.

The 2D probe interpolation (fault surface) can use MFEM's `GridFunction::GetValue(elem, ip)` after finding the element containing the probe point via `mesh.FindPoints()`.

**Step 5.3: Create BP5 driver program.**

File: `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/bp5/bp5_verification.cpp`

This follows the structure of `tests/verification/bp1_verification_full.cpp`:

```
1. Load/generate mesh from bp5.msh (or generate inline for tests)
2. Construct BP5Params
3. Construct ElasticityDomainOperator (3D elasticity, BP5 geometry)
4. Construct FaultGeometry3D (with 2D a(x2,x3), L(x2,x3))
5. Construct RateStateFaultOperator<Mesh, 2> (2 slip components)
6. Construct SEASQuasiDynamicOperator<Mesh, ElasticityDomainOperator<Mesh>>
7. SetInitialCondition
8. Time stepping loop with adaptive dt
9. BP5BenchmarkOutput at each step
10. Final output to disk
```

**Step 5.4: Testing strategy.**

The testing sequence mirrors the BP1/BP2 approach:

a. **MMS test for elasticity operator**: Construct a manufactured solution for 3D linear elasticity with a known displacement jump across `x1 = 0`. Verify second-order (for DG-IP) or third-order (for DG-BR2) convergence of displacement error in L2 norm.

b. **Traction extraction test**: With a known displacement field (e.g., from MMS), verify that `ComputeFullTraction` returns the correct `(tau2, tau3, sigma_n)`.

c. **Static equilibrium test**: Apply a constant slip `s2 = s0` (no `s3`) uniformly on the fault and verify that the resulting `tau2` matches the analytical formula for strike-slip on a half-space.

d. **Vector friction test**: With known `(tau2, tau3)`, verify that `SolveSlipRateVector` returns `V` with correct magnitude and direction.

e. **Full BP5 short simulation**: Run for 10 years and verify `V_max ≈ V_p` (pre-seismic), then nucleation around the initial perturbation.

---

## Part C: Key Design Decisions

### C.1 Template Strategy for Scalar to Vector Generalization

The current codebase uses `template <typename MeshType = Mesh>` for serial/parallel switching. The natural extension for scalar-to-vector generalization is to add a `int SlipComponents = 1` template parameter to `RateStateFaultOperator` and a `typename DomainOpType` parameter to `SEASQuasiDynamicOperator`.

The alternative — runtime polymorphism via virtual functions on the fault operator — would require boxing the state vector and adding indirection to the hot inner loop. Since `RateStateFaultOperator::ComputeRHS` is called at every ODE stage (potentially thousands of times per earthquake), compile-time specialization via templates is preferred.

However, to avoid template explosion, the recommendation is:
- Use `SlipComponents` as a template parameter only on `RateStateFaultOperator`.
- Expose the result as flat `Vector` objects (already done) so the SEAS operator, ODE solver, and I/O system all remain template-free.
- The domain operator type in `SEASQuasiDynamicOperator` should default to the antiplane operator for backwards compatibility, with the 3D type supplied explicitly for BP5.

This approach closely mirrors Tandem's design, where `RateAndState<Law>` is templated on the friction law, `RateAndStateBase` is parameterized by `TangentialComponents = DomainDimension - 1`, and the domain DG operator (`Elasticity` vs `Poisson`) is swapped by changing the `type` field in the configuration file.

### C.2 Fault Traction Decomposition (Normal Stress, Shear Stress Vector) in 3D

The full Cauchy stress on the fault face with outward normal `n = (1, 0, 0)` (normal to the fault plane `x1 = 0`, pointing into the positive side) gives a traction vector:

```
T = sigma . n = (sigma_11, sigma_21, sigma_31)
```

The three components have specific roles:
- `T_1 = sigma_11`: Normal traction (negative = compression). Normal stress `sigma_n = -T_1`.
- `T_2 = sigma_21`: Along-strike shear traction (positive = right-lateral).
- `T_3 = sigma_31`: Along-dip shear traction (positive = reverse on positive side).

For BP5, since the material is homogeneous and the no-opening condition holds, `sigma_11` is continuous across the fault and equals the pre-stress value. Therefore only the tangential components `(sigma_21, sigma_31)` need to be extracted and passed to the friction solver.

The traction is extracted using the DG numerical flux:

```
tau_h = {{sigma(u_h)}} . n_fault = 0.5 * (sigma(u_h)|_+ + sigma(u_h)|_-) . n_fault
```

where `sigma(u)|_±` are the stress tensors evaluated from the `±` sides of each fault face. This is the standard DG average, and is what `AntiplaneDomainOperator::ComputeTraction` already computes for the scalar case.

For the 3D elasticity operator, the stress tensor computation at a quadrature point `q` is:

```
sigma_ij(u_h) = lambda * div(u_h) * delta_ij + mu * (du_i/dx_j + du_j/dx_i)
```

The traction components needed are:
```
tau2 = sigma_21 = mu * (du_2/dx_1 + du_1/dx_2)
tau3 = sigma_31 = mu * (du_3/dx_1 + du_1/dx_3)
sigma_n = -sigma_11 = -(lambda * div(u) + 2*mu * du_1/dx_1)
```

In the DG average, for a face with normal `n = (1, 0, 0)`, only the `x1`-derivatives of `u1, u2, u3` and `x2, x3`-derivatives of `u1` contribute. Practically, `ComputeFullTraction` iterates over fault interior faces and for each face:

1. Evaluates the displacement gradient tensor from both `Elem1` and `Elem2` at face quadrature points.
2. Computes `sigma` on both sides.
3. Averages: `{{sigma}} = 0.5*(sigma^+ + sigma^-)`.
4. Extracts `T = {{sigma}} . (1, 0, 0)`.
5. Projects `(T_2, T_3)` onto fault DOFs using the same L2 projection as the scalar case (mass-matrix weighted face integral).

### C.3 How to Handle SH/P/SV Wave Coupling in Quasi-Dynamic Radiation Damping

This is the most subtle aspect of BP5. In a 3D strike-slip problem, the physical waves excited by fault slip are:
- **SH wave** (antiplane shear): generated by along-strike slip (`s2`).
- **P and SV waves** (in-plane): generated by along-dip slip (`s3`).

In the full dynamic problem, the radiation damping is anisotropic: the SH wave radiates from strike-slip with `eta_SH = mu/(2*cs)`, while the P/SV waves radiate from dip-slip with `eta_P = (lambda + 2*mu)/(2*cp)` for the P component and `eta_SV = mu/(2*cs)` for the SV component. However, in the quasi-dynamic approximation, only the dominant shear wave contribution is retained.

**BP5 uses the isotropic shear-wave approximation:** `eta = mu/(2*cs)` applied identically to both `V2` and `V3` components. The spec (Eq. 10) writes `tau = tau_0 + Delta_tau - eta * V` where `V = (V2, V3)`, meaning the same scalar `eta` multiplies each component. This is confirmed in Tandem's `DieterichRuinaAgeing`:

```cpp
struct Params {
    double eta;    // Single scalar, applied to both components
    ...
};
// In tau_hat: return tau + tau_pre + eta * V  (component-wise scalar multiply)
```

This simplification avoids the need to distinguish SH from P/SV modes in the quasi-dynamic framework. For the MFEM implementation, `eta = mu/(2*cs)` is a single scalar stored per fault node (constant for homogeneous medium), and it multiplies each component of the velocity vector independently. No cross-component coupling terms need to be implemented.

If in the future a more accurate anisotropic radiation damping is needed (for example, following Sanz-Alonso & Simpson 2024 or Lapusta et al. 2000), the `FaultGeometry` struct could be extended to store `eta_2` (along-strike) and `eta_3` (along-dip) separately, and the `SolveSlipRateVector` call would use anisotropic damping. This is a natural extension point.

### C.4 Boundary Conditions for the 3D Elastic Problem

The BP5 domain has six external faces. Following Tandem's `bp5.lua` and `bp5.geo`:

**Free surface (`x3 = 0`)**: Natural BC (zero traction). This is critical for BP5 because surface effects (surface waves, free-surface reflection) are physically important for the 3D problem. In DG, Natural BC requires no special treatment — the face integral simply vanishes for Neumann-zero.

**Bottom (`x3 = Lz`)**: Natural BC (zero traction). The domain is large enough that the bottom boundary has negligible influence.

**Strike-parallel walls (`x2 = ±Ly`)**: Dirichlet loading `u2 = ±Vp*t/2`, `u1 = 0`, `u3 = 0`. This applies the far-field plate motion as a Dirichlet condition on the `x2`-displacement at the along-strike boundaries. In DG, this is applied weakly using `DGElasticityIntegrator` on the boundary faces with the Dirichlet data. The time-dependence requires updating the linear form RHS at each time step (the stiffness matrix is assembled once and cached).

**Fault-normal walls (`x1 = ±Lx`)**: These are far from the fault. Following Tandem's `diri()` surface (Physical Surface 5), these are Dirichlet zero (zero displacement). However, more careful analysis suggests Natural BC (zero traction) may be preferable for a truncated domain, since zero displacement can reflect stress waves. The BP5 spec says the choice of truncation BC affects results; the recommendation is to use Natural BC (simpler, avoids spurious reflections in QD) and push the domain boundaries out until results converge.

**Fault interior surface (`x1 = 0`, RSF region)**: Not an external boundary but an interior interface. Treated via face integrators in the DG formulation, with the slip jump as the interface condition.

The boundary condition analysis is already well-documented for the 2D case in `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/document/tandem_boundary_condition_analysis.md`. For 3D BP5, the same principles apply: the loading is applied via Dirichlet BC on the `±x2` boundaries (far-field strike-parallel walls), mimicking the plate motion, while all other external boundaries use Natural BC.

### C.5 Pre-stress Initialization for the 3D Vector Case

For BP5, the pre-stress is a 2-vector `tau_0 = (tau2_0, tau3_0)`. The spec (Eq. 17) specifies:

```
tau_0 = tau^0 * V_init / ||V_init||
```

where `tau^0` is the scalar steady-state friction stress and `V_init / ||V_init||` provides the direction. Since the initial velocity is `V_init = (0, Vp)` outside the nucleation zone (only along-strike motion), the initial pre-stress is `tau_0 = (0, -tau^0)` (along-dip component zero, along-strike component equal to steady-state friction in the `-x2` direction — the sign is such that the fault is loaded in the right-lateral direction).

In the nucleation zone, `V_init = (Vzero, 0.03)`, so the direction is mostly in the `x2`-direction, with a tiny `x3` component from `Vzero`.

The `BP5Params::tau0_vec(x2, x3)` method must:
1. Compute `a = a(x2, x3)` and `L = L(x2, x3)`.
2. Compute `V_init = V_init(x2, x3)` (from nucleation zone test).
3. Compute `V_abs = ||V_init||`.
4. Compute scalar `tau_0 = sigma_n * a * asinh[(V_abs/(2*V0)) * exp((f0 + b*ln(V0/V_abs))/a)] + eta * V_abs`.
5. Return `tau_0_vec = -tau_0 * V_init / V_abs`.

The negative sign follows Tandem's convention where `tau_pre` is stored as the negative of the resisting stress (so that the physical stress is `tau + tau_pre`).

---

## Summary of New Files Required

The following new files need to be created for BP5 (no existing files need to be deleted or restructured, only extended):

| New File | Purpose |
|----------|---------|
| `config/bp5_params.hpp` | BP5 physical parameters, 2D a(x2,x3), L(x2,x3), tau_pre_vec |
| `domain/elasticity3d_operator.hpp` | 3D full elasticity DG operator, inherits DomainOperator |
| `domain/bp5_mesh.hpp` | 3D mesh generator, boundary attributes for BP5 |
| `bp5/bp5.geo` | Gmsh geometry for BP5 3D domain with embedded fault |
| `integrator/dg_elasticity_br2_integrator.hpp` | BR2 lifting for elasticity (Phase 3b, optional) |
| `io/bp5_benchmark_output.hpp` | BP5 output format (2D probes, off-fault stations, catalog) |
| `bp5/bp5_verification.cpp` | BP5 driver program |

The following existing files need to be modified (all additive, backwards-compatible):

| Existing File | Modification |
|---------------|-------------|
| `domain/domain_operator.hpp` | Add `NumSlipComponents()`, `GetFaultCoords2D()`, `GetOffFaultDisplacement()` with default implementations |
| `solver/seas_operator.hpp` | Add `DomainOpType` as second template parameter (with default = antiplane) |
| `fault/rate_state_fault.hpp` | Add `SlipComponents` template parameter (default = 1), generalize state layout |
| `friction/friction_law.hpp` | Add `SolveSlipRateVector()` virtual method with default scalar delegation |
| `friction/dieterich_ruina.hpp` | Implement `SolveSlipRateVector()` following Tandem's vector direction logic |
| `fault/fault_geometry.hpp` | Extend to store 2D coords, generalize GatherToRootDedup for 2D |

---

### Critical Files for Implementation

- `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/domain/domain_operator.hpp` — Abstract base class to extend with `NumSlipComponents()`, `GetFaultCoords2D()`, and `GetOffFaultDisplacement()` virtual methods; this file governs the interface that all three new components (3D domain operator, 3D fault geometry, BP5 output) must satisfy.

- `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/fault/rate_state_fault.hpp` — Core fault state operator that requires the `SlipComponents` template generalization and vector slip rate logic; the state layout change from `[slip, theta]` to `[s2, s3, psi]` per node is the single most invasive change in the entire codebase.

- `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/domain/antiplane_operator.hpp` — The primary reference pattern for implementing `ElasticityDomainOperator`; the fault-face identification, DOF extraction, slip RHS assembly, and traction computation patterns here are all directly adaptable to 3D, making this the most valuable blueprint file.

- `/Users/chunhuizhao/projects/tandem/app/localoperator/Elasticity.h` and `/Users/chunhuizhao/projects/tandem/app/localoperator/Elasticity.cpp` — The reference 3D elasticity DG operator in Tandem; the `assemble_skeleton` and `traction_skeleton` methods, the stress tensor computation, and the BR2 penalty formula here are the ground truth for the MFEM port.

- `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/solver/seas_operator.hpp` — The SEAS time integration driver that needs only the second template parameter addition to accommodate the 3D elasticity domain operator while remaining backwards-compatible with all BP1/BP2 code.

---
