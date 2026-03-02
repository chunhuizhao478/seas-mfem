# Fault Interface Implementation Guide: Tandem to MFEM

## Executive Summary

This document provides a detailed implementation roadmap for porting Tandem's fault interface treatment to MFEM. It includes:
1. Complete analysis of Tandem's fault coupling architecture
2. Mapping to MFEM's face integration infrastructure
3. Concrete code implementation suggestions with examples

---

# Part I: Tandem Fault Interface Architecture

## 1. Overview of Tandem's Coupling Architecture

Tandem implements fault interfaces through a **three-operator coupling pattern**:

```
┌─────────────────────────────────────────────────────────────────┐
│                     SeasQDOperator                              │
│  (Coordinates domain-fault coupling for quasi-dynamic mode)     │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│   ┌──────────────┐      ┌──────────────┐      ┌──────────────┐ │
│   │  DGOperator  │◄────►│AdapterOperator│◄────►│FrictionOper- │ │
│   │  (Elasticity)│      │   (Coupling)  │      │    ator      │ │
│   └──────────────┘      └──────────────┘      └──────────────┘ │
│         │                      │                      │         │
│         ▼                      ▼                      ▼         │
│   Displacement u         Traction T              State (S, ψ)   │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### Data Flow at Each Time Step

```
1. State → Slip BC
   FrictionOperator holds state (S, ψ)
         ↓
   AdapterOperator.slip_bc(state) → λ(fctNo, f_q)
         ↓
   Returns function that interpolates S to quadrature points

2. Slip BC → Displacement
   DGOperator.set_slip(slip_bc_function)
         ↓
   Assembles system with slip as Dirichlet BC on fault
         ↓
   LinearSolver.solve() → u

3. Displacement → Traction
   AdapterOperator.traction(u, T)
         ↓
   Computes σ(u)·n at fault quadrature points
         ↓
   Projects to fault nodes: T_nodal = M^{-1} ∫ φ_i T_q dA

4. Traction → State Rate
   FrictionOperator.rhs(T, state, dstate/dt)
         ↓
   For each fault node:
     V = slip_rate(σ_n, τ, ψ)    [Newton solver]
     dS/dt = V
     dψ/dt = state_evolution(V, ψ)
```

---

## 2. FrictionOperator: Detailed Implementation

### 2.1 Class Structure

```cpp
// File: app/form/FrictionOperator.h
template <typename LocalOperator>
class FrictionOperator : public AbstractFrictionOperator {
private:
    // Local operator containing friction law (RateAndState<DieterichRuinaAgeing>)
    std::unique_ptr<LocalOperator> lop_;

    // Topology reference for domain
    std::shared_ptr<DGOperatorTopo> topo_;

    // Fault-to-facet mapping
    std::shared_ptr<BoundaryMap> fault_map_;

    // Scratch memory for temporary computations
    Scratch<double> scratch_;

    // Maximum slip velocity (for CFL/output)
    double VMax_ = 0.0;

public:
    // Number of DOFs per fault element
    std::size_t block_size() const { return lop_->block_size(); }

    // Number of local fault elements
    std::size_t num_local_elements() const { return fault_map_->local_size(); }
};
```

### 2.2 State Vector Layout

```
State vector per fault element (block):
┌─────────────────────────────────────────────────────────────┐
│ Node 0: [S_x, S_y, ψ]  │ Node 1: [S_x, S_y, ψ] │ ... │     │
└─────────────────────────────────────────────────────────────┘

For 3D with P=1 triangular fault elements:
  - 3 nodes per element
  - 3 quantities per node (S_x, S_y, ψ)
  - Block size = 3 × 3 = 9

Memory layout (row-major):
  state[0] = S_x at node 0
  state[1] = S_y at node 0
  state[2] = ψ at node 0
  state[3] = S_x at node 1
  ...
```

### 2.3 Key Methods

#### `pre_init()` - Set Initial Slip

```cpp
void FrictionOperator::pre_init(BlockVector& state) {
    auto state_handle = state.begin_access_readwrite();

    for (std::size_t faultNo = 0; faultNo < num_local_elements(); ++faultNo) {
        auto state_block = state_handle.subtensor(slice{}, faultNo);
        lop_->pre_init(faultNo, state_block, scratch_);
    }

    state_handle.end_access();
}

// In RateAndState<Law>::pre_init():
void pre_init(std::size_t faultNo, Vector<double>& state, ...) {
    auto s_mat = reshape(state, nbf, NumQuantities);
    std::size_t index = faultNo * nbf;

    for (std::size_t node = 0; node < nbf; ++node) {
        // Set initial slip from parameters
        auto Sinit = law_.S_init(index + node);  // Returns array<double, D-1>
        for (std::size_t t = 0; t < TangentialComponents; ++t) {
            s_mat(node, t) = Sinit[t];
        }
        // Note: state variable ψ left uninitialized (computed in init())
    }
}
```

#### `init()` - Initialize State Variable from Traction

```cpp
void FrictionOperator::init(double time,
                            BlockVector const& traction,
                            BlockVector& state) {
    VMax_ = 0.0;
    auto traction_handle = traction.begin_access_readonly();
    auto state_handle = state.begin_access_readwrite();

    for (std::size_t faultNo = 0; faultNo < num_local_elements(); ++faultNo) {
        auto traction_block = traction_handle.subtensor(slice{}, faultNo);
        auto state_block = state_handle.subtensor(slice{}, faultNo);

        double local_VMax = lop_->init(time, faultNo,
                                        traction_block, state_block, scratch_);
        VMax_ = std::max(VMax_, local_VMax);
    }
}

// In RateAndState<Law>::init():
double init(double time, std::size_t faultNo,
            Vector<double const> const& traction,
            Vector<double>& state, ...) {

    auto s_mat = reshape(state, nbf, NumQuantities);
    auto t_mat = reshape(traction, nbf, DomainDimension);
    double VMax = 0.0;
    std::size_t index = faultNo * nbf;

    for (std::size_t node = 0; node < nbf; ++node) {
        // Extract stress components
        double sn = t_mat(node, 0);  // Normal stress (component 0)
        auto tau = get_tau(node, t_mat);  // Shear stress (components 1:D-1)

        // Compute initial state variable from equilibrium
        double psi = law_.psi_init(index + node, sn, tau);
        s_mat(node, PsiIndex) = psi;

        // Compute initial slip rate for output
        auto V = law_.slip_rate(index + node, sn, tau, psi);
        VMax = std::max(VMax, norm(V));
    }

    return VMax;
}
```

#### `rhs()` - Compute State Evolution Rates

```cpp
void FrictionOperator::rhs(double time,
                           BlockVector const& traction,
                           BlockVector const& state,
                           BlockVector& result) {
    VMax_ = 0.0;
    auto traction_handle = traction.begin_access_readonly();
    auto state_handle = state.begin_access_readonly();
    auto result_handle = result.begin_access_readwrite();

    for (std::size_t faultNo = 0; faultNo < num_local_elements(); ++faultNo) {
        auto traction_block = traction_handle.subtensor(slice{}, faultNo);
        auto state_block = state_handle.subtensor(slice{}, faultNo);
        auto result_block = result_handle.subtensor(slice{}, faultNo);

        double local_VMax = lop_->rhs(time, faultNo,
                                       traction_block, state_block,
                                       result_block, scratch_);
        VMax_ = std::max(VMax_, local_VMax);
    }
}

// In RateAndState<Law>::rhs():
double rhs(double time, std::size_t faultNo,
           Vector<double const> const& traction,
           Vector<double const>& state,
           Vector<double>& result, ...) {

    auto s_mat = reshape(state, nbf, NumQuantities);
    auto r_mat = reshape(result, nbf, NumQuantities);
    auto t_mat = reshape(traction, nbf, DomainDimension);
    double VMax = 0.0;
    std::size_t index = faultNo * nbf;

    for (std::size_t node = 0; node < nbf; ++node) {
        double sn = t_mat(node, 0);
        double psi = s_mat(node, PsiIndex);
        auto tau = get_tau(node, t_mat);

        // Core friction law: solve for slip rate
        auto Vi = law_.slip_rate(index + node, sn, tau, psi);
        double V = norm(Vi);
        VMax = std::max(VMax, V);

        // Store slip rate: dS/dt = V
        for (std::size_t t = 0; t < TangentialComponents; ++t) {
            r_mat(node, t) = Vi[t];
        }

        // Store state rate: dψ/dt
        r_mat(node, PsiIndex) = law_.state_rhs(index + node, V, psi);
    }

    return VMax;
}
```

---

## 3. DieterichRuinaAgeing Friction Law

### 3.1 Mathematical Formulation

**Governing Equation (Balance at fault)**:
```
τ_abs = σ_n · f(V, ψ) + η · V

where:
  τ_abs = |τ + τ_pre|           (absolute shear stress magnitude)
  σ_n = sn_pre - sn             (absolute normal stress, compression positive)
  f(V, ψ) = a · asinh((V/(2V₀)) · exp(ψ/a))   (friction coefficient)
  η = radiation damping coefficient
```

**State Evolution**:
```
dψ/dt = (b/L) · (V₀ · exp((f₀ - ψ)/b) - V)

At steady state (dψ/dt = 0):
  ψ_ss = f₀ - b · ln(V/V₀)
```

### 3.2 Slip Rate Solver (Newton-Raphson)

```cpp
auto slip_rate(std::size_t index, double sn,
               std::array<double, TangentialComponents> const& tau,
               double psi) const
    -> std::array<double, TangentialComponents>
{
    auto eta = p_[index].get<Eta>();
    auto tauAbsVec = tau + p_[index].get<TauPre>();
    double snAbs = -sn + p_[index].get<SnPre>();
    double tauAbs = norm(tauAbsVec);

    double V = 0.0;

    if (eta == 0.0) {
        // No radiation damping: closed-form solution
        V = Finv(index, snAbs, tauAbs, psi);
    } else if (snAbs <= 0.0) {
        // Fault in tension: viscous sliding
        V = tauAbs / eta;
    } else {
        // With radiation damping: solve R(V) = 0
        // R(V) = τ_abs - σ_n · f(V, ψ) - η · V = 0

        // Bracket the solution
        double a = 0.0;
        double b = tauAbs / eta;

        // Root finding using Brent's method
        auto residual = [this, &index, &snAbs, &tauAbs, &psi, &eta](double V) {
            return tauAbs - this->F(index, snAbs, V, psi) - eta * V;
        };
        V = zeroIn(a, b, residual);
    }

    // Return slip rate vector in shear stress direction
    return -(V / tauAbs) * tauAbsVec;
}

// Friction strength function
double F(std::size_t index, double snAbs, double V, double psi) const {
    auto a = p_[index].get<A>();
    double e = exp(psi / a);
    double f = a * asinh((V / (2.0 * cp_.V0)) * e);
    return snAbs * f;
}

// Inverse (no radiation damping)
double Finv(std::size_t index, double snAbs, double tauAbs, double psi) const {
    auto a = p_[index].get<A>();
    double r = tauAbs / snAbs;
    return cp_.V0 * (exp((r - psi) / a) - exp(-(r + psi) / a));
}
```

### 3.3 State Evolution

```cpp
double state_rhs(std::size_t index, double V, double psi) const {
    double L = p_[index].get<L>();
    // dψ/dt = (b/L) · (V₀ · exp((f₀ - ψ)/b) - V)
    return cp_.b * cp_.V0 / L * (exp((cp_.f0 - psi) / cp_.b) - V / cp_.V0);
}
```

### 3.4 Initial State Computation

```cpp
double psi_init(std::size_t index, double sn,
                std::array<double, TangentialComponents> const& tau) const {
    double snAbs = -sn + p_[index].get<SnPre>();
    double tauAbs = norm(tau + p_[index].get<TauPre>());
    auto Vi = norm(p_[index].get<Vinit>());

    if (Vi == 0.0) return cp_.f0;

    auto a = p_[index].get<A>();
    auto eta = p_[index].get<Eta>();

    // Solve: τ_abs = σ_n · f(V_init, ψ_init) + η · V_init
    // Rearrange: sinh((τ_abs - η·V)/(a·σ_n)) = (V/(2V₀)) · exp(ψ/a)
    // Therefore: ψ = a · ln((2V₀/V) · sinh(...))

    double s = sinh((tauAbs - eta * Vi) / (a * snAbs));
    double l = log((2.0 * cp_.V0 / Vi) * s);
    return a * l;
}
```

---

## 4. AdapterOperator: Domain-Fault Coupling

### 4.1 Class Structure

```cpp
// File: app/form/AdapterOperator.h
template <typename LocalOperator>
class AdapterOperator : public AbstractAdapterOperator {
private:
    std::shared_ptr<LocalOperator> adapted_lop_;  // Elasticity operator
    std::unique_ptr<Adapter<LocalOperator>> lop_; // Adapter local operator
    std::shared_ptr<DGOperatorTopo> topo_;
    std::shared_ptr<BoundaryMap> fault_map_;
    Scratch<double> scratch_;
};
```

### 4.2 Slip Boundary Condition

```cpp
auto slip_bc(BlockView const& state)
    -> std::function<void(std::size_t, Matrix<double>&, bool)>
{
    return [this, &state](std::size_t fctNo, Matrix<double>& f_q, bool) {
        // Reverse lookup: facet number → fault number
        auto faultNo = fault_map_->bndNo(fctNo);
        if (faultNo == BoundaryMap::INVALID) return;

        // Get state block for this fault element
        auto state_block = state.get_block(faultNo);

        // Interpolate slip to quadrature points
        lop_->slip(faultNo, state_block, f_q);
    };
}

// In Adapter<>::slip():
void slip(std::size_t faultNo, Vector<double const> const& state,
          Matrix<double>& f_q) {
    auto s_mat = reshape(state, nbf, NumQuantities);

    // f_q has shape [DomainDimension × num_quad_points]
    for (std::size_t q = 0; q < num_quad_points; ++q) {
        // Interpolate slip components using basis functions
        for (std::size_t d = 0; d < DomainDimension; ++d) {
            double slip_d = 0.0;
            for (std::size_t i = 0; i < nbf; ++i) {
                if (d < TangentialComponents) {
                    slip_d += basis_q(i, q) * s_mat(i, d);
                }
                // Normal component (d == D-1) is zero (no opening)
            }
            f_q(d, q) = slip_d;
        }

        // Transform from fault coordinates to domain coordinates
        // f_q[:, q] = FaultBasis · f_q[:, q]
        transform_to_domain(faultNo, q, f_q);
    }
}
```

### 4.3 Traction Computation

```cpp
void traction(BlockView const& displacement, BlockVector& result) {
    auto result_handle = result.begin_access_readwrite();

    for (std::size_t faultNo = 0; faultNo < num_local_elements(); ++faultNo) {
        auto fctNo = fault_map_->fctNo(faultNo);
        auto const& info = topo_->info(fctNo);

        // Get displacement from adjacent domain elements
        auto u0 = displacement.get_block(info.up[0]);
        auto u1 = displacement.get_block(info.up[1]);

        // Allocate temporary for traction at quadrature points
        Matrix<double> traction_q(DomainDimension, num_quad_points);

        // Compute traction using elasticity operator
        if (info.up[0] == info.up[1]) {
            // Boundary fault: single element
            adapted_lop_->traction_boundary(fctNo, info, u0, traction_q);
        } else {
            // Interior fault: two elements
            adapted_lop_->traction_skeleton(fctNo, info, u0, u1, traction_q);
        }

        // Transform to fault coordinates and project to nodes
        auto result_block = result_handle.subtensor(slice{}, faultNo);
        lop_->traction(faultNo, traction_q, result_block, scratch_);
    }
}

// In Adapter<>::traction():
void traction(std::size_t faultNo, Matrix<double> const& traction_q,
              Vector<double>& result, ...) {
    auto t_mat = reshape(result, nbf, DomainDimension);

    // For each component d:
    for (std::size_t d = 0; d < DomainDimension; ++d) {
        // Project quadrature values to nodes: t_nodal = M^{-1} ∫ φ_i t_q dA
        for (std::size_t i = 0; i < nbf; ++i) {
            double sum = 0.0;
            for (std::size_t q = 0; q < num_quad_points; ++q) {
                // Transform traction to fault coordinates
                double t_fault_d = 0.0;
                for (std::size_t k = 0; k < DomainDimension; ++k) {
                    t_fault_d += fault_basis_q(d, k, faultNo, q) * traction_q(k, q);
                }
                sum += weight_q(q) * basis_q(i, q) * t_fault_d * normal_length_q(q);
            }
            // Apply mass matrix inverse
            double t_nodal = 0.0;
            for (std::size_t j = 0; j < nbf; ++j) {
                t_nodal += mass_inv(i, j) * integral(j);
            }
            t_mat(i, d) = t_nodal;
        }
    }
}
```

### 4.4 Fault Basis Construction

```cpp
// Construct orthonormal basis aligned with fault
void prepare(std::size_t faultNo, FacetInfo const& info, ...) {
    // 1. Compute face normal at each quadrature point
    for (std::size_t q = 0; q < num_quad_points; ++q) {
        auto normal = compute_normal(faultNo, q);

        // 2. Check orientation against reference normal
        double dot = ref_normal_ · normal;
        bool flip = (dot < 0);
        if (flip) normal = -normal;

        // 3. Construct fault basis using up direction
        // Column 0: normal (outward from "-" side)
        // Columns 1 to D-1: tangent directions
        auto basis = construct_orthonormal_basis(normal, up_direction_);
        if (flip) basis = -basis;  // Flip all if needed

        // Store
        fault_basis_q[faultNo][q] = basis;
        sign_flipped_[faultNo][q] = flip;
    }

    // 4. Precompute mass matrix and inverse
    for (i, j in basis_functions) {
        M(i, j) = ∑_q weight(q) * normal_length(q) * φ_i(q) * φ_j(q);
    }
    M_inv = M.inverse();
}
```

---

## 5. SeasQDOperator: Time Integration Coupling

### 5.1 Overall Flow

```cpp
class SeasQDOperator {
public:
    void initial_condition(BlockVector& state) {
        // Step 1: Set initial slip
        friction_->pre_init(state);

        // Step 2: Scatter state to ghost partitions
        update_ghost_state(state);

        // Step 3: Solve domain with initial slip BC
        solve(0.0, make_state_view(state));

        // Step 4: Compute initial traction
        update_traction(make_state_view(state));

        // Step 5: Initialize state variable from traction
        friction_->init(0.0, traction_, state);
    }

    void rhs(double time, BlockVector const& state, BlockVector& result) {
        // Called by time integrator to get d(state)/dt

        // Step 1: Parallel communication
        update_ghost_state(state);

        // Step 2: Solve domain equilibrium
        solve(time, make_state_view(state));

        // Step 3: Compute traction from displacement
        update_traction(make_state_view(state));

        // Step 4: Compute friction rates
        friction_->rhs(time, traction_, state, result);
    }
};
```

### 5.2 Domain Solve with Slip BC

```cpp
void solve(double time, BlockView const& state_view) {
    // Set slip as boundary condition on fault
    dgop_->set_slip(adapter_->slip_bc(state_view));

    // Optional: set external Dirichlet BCs
    if (fun_boundary_) {
        dgop_->set_dirichlet((*fun_boundary_)(time));
    }

    // Assemble and solve: A·u = b
    linear_solver_.update_rhs(*dgop_);
    linear_solver_.solve();

    // Scatter solution for traction computation
    dgop_->set_slip(invalid_slip_bc());
    disp_scatter_.begin_scatter(linear_solver_.x(), disp_ghost_);
    disp_scatter_.wait_scatter();
}
```

---

# Part II: MFEM Implementation Guide

## 6. Mapping Tandem Concepts to MFEM

| Tandem Component | MFEM Equivalent | Implementation Strategy |
|------------------|-----------------|-------------------------|
| `DGOperator` | `BilinearForm` + Integrators | Use existing infrastructure |
| `FrictionOperator` | New `FaultOperator` class | Custom implementation |
| `AdapterOperator` | New `FaultCouplingOperator` | Custom face integrator |
| `BlockVector` | `BlockVector` or `Vector` | MFEM has BlockVector |
| `BoundaryMap` | Boundary attribute + face list | Use mesh attributes |
| `SeasQDOperator` | `TimeDependentOperator` subclass | Custom operator |
| `PetscTimeSolver` | `ODESolver` hierarchy | Use existing solvers |

---

## 7. Proposed MFEM Class Architecture

### 7.1 Class Diagram

```
TimeDependentOperator
└── SEASOperator
    ├── QuasiDynamicOperator    (implicit domain solve)
    └── FullyDynamicOperator    (explicit elastodynamics)

BilinearFormIntegrator
├── DGElasticityIntegrator     (existing, for domain)
└── FaultTractionIntegrator    (new, computes traction)

LinearFormIntegrator
└── FaultSlipLFIntegrator      (new, applies slip BC)

Coefficient
├── DieterichRuinaFriction     (new, friction coefficient)
└── FaultParameterCoefficient  (new, spatially varying params)

Custom Classes:
├── FaultOperator              (manages fault state evolution)
├── FaultCouplingOperator      (domain-fault interface)
├── FaultGeometry              (fault topology and basis)
└── RateStateFrictionLaw       (friction law implementation)
```

### 7.2 File Organization

```
miniapps/seas/
├── CMakeLists.txt
├── seas.cpp                    # Main driver
├── seas_solver.hpp             # SEASOperator implementation
│
├── friction/
│   ├── friction_law.hpp        # Base class for friction laws
│   ├── dieterich_ruina.hpp     # Dieterich-Ruina ageing law
│   └── slip_weakening.hpp      # Alternative friction law
│
├── operators/
│   ├── fault_operator.hpp      # FaultOperator (state evolution)
│   ├── fault_coupling.hpp      # FaultCouplingOperator
│   └── seas_operator.hpp       # QuasiDynamic/FullyDynamic
│
├── integrators/
│   ├── fault_traction_integrator.hpp
│   └── fault_slip_integrator.hpp
│
├── geometry/
│   └── fault_geometry.hpp      # Fault surface handling
│
└── config/
    └── seas_config.hpp         # Configuration parsing
```

---

## 8. Detailed Implementation: FaultOperator

### 8.1 Header File

```cpp
// File: miniapps/seas/operators/fault_operator.hpp

#ifndef MFEM_SEAS_FAULT_OPERATOR_HPP
#define MFEM_SEAS_FAULT_OPERATOR_HPP

#include "mfem.hpp"
#include "../friction/friction_law.hpp"
#include "../geometry/fault_geometry.hpp"

namespace mfem {
namespace seas {

/**
 * @brief Manages fault state evolution for SEAS simulations
 *
 * State vector layout per fault face (for 3D):
 *   [slip_strike, slip_dip, psi] at each integration point
 *
 * Traction vector layout:
 *   [sigma_n, tau_strike, tau_dip] at each integration point
 */
class FaultOperator {
public:
    /**
     * @param fault_fes  Finite element space on fault surface (L2 or DG)
     * @param friction   Friction law implementation
     * @param fault_geom Fault geometry information
     */
    FaultOperator(FiniteElementSpace *fault_fes,
                  RateStateFrictionLaw *friction,
                  FaultGeometry *fault_geom);

    ~FaultOperator();

    /// Number of state DOFs (slip components + state variable)
    int GetStateSize() const;

    /// Number of traction DOFs (normal + shear components)
    int GetTractionSize() const;

    /// Initialize slip to prescribed initial values
    void PreInit(Vector &state);

    /// Initialize state variable from initial traction
    /// @param traction Current traction at fault points
    /// @param state    State vector (slip already set, psi computed)
    void Init(real_t time, const Vector &traction, Vector &state);

    /// Compute state evolution rates: d(state)/dt
    /// @param traction Current traction at fault points
    /// @param state    Current state (slip, psi)
    /// @param rate     Output: d(slip)/dt, d(psi)/dt
    /// @return Maximum slip velocity (for CFL/output)
    real_t ComputeRHS(real_t time,
                      const Vector &traction,
                      const Vector &state,
                      Vector &rate);

    /// Get current maximum slip velocity
    real_t GetMaxSlipVelocity() const { return VMax_; }

    /// Project state to GridFunction for visualization
    void GetStateGridFunction(const Vector &state,
                              GridFunction &slip_gf,
                              GridFunction &psi_gf);

protected:
    FiniteElementSpace *fault_fes_;
    RateStateFrictionLaw *friction_;
    FaultGeometry *fault_geom_;

    int dim_;                    // Domain dimension
    int num_slip_components_;    // D-1
    int num_state_components_;   // D (slip + psi)
    int num_traction_components_;// D (normal + shear)

    real_t VMax_;                // Maximum slip velocity

    // Work arrays
    Vector slip_rate_;
    Vector state_rate_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_OPERATOR_HPP
```

### 8.2 Implementation

```cpp
// File: miniapps/seas/operators/fault_operator.cpp

#include "fault_operator.hpp"

namespace mfem {
namespace seas {

FaultOperator::FaultOperator(FiniteElementSpace *fault_fes,
                             RateStateFrictionLaw *friction,
                             FaultGeometry *fault_geom)
    : fault_fes_(fault_fes),
      friction_(friction),
      fault_geom_(fault_geom),
      VMax_(0.0)
{
    dim_ = fault_fes_->GetMesh()->Dimension() + 1;  // Fault is D-1 dimensional
    num_slip_components_ = dim_ - 1;
    num_state_components_ = dim_;  // slip + psi
    num_traction_components_ = dim_;  // normal + shear
}

int FaultOperator::GetStateSize() const
{
    return fault_fes_->GetTrueVSize() * num_state_components_;
}

int FaultOperator::GetTractionSize() const
{
    return fault_fes_->GetTrueVSize() * num_traction_components_;
}

void FaultOperator::PreInit(Vector &state)
{
    // Set initial slip from friction law parameters
    const int npts = fault_fes_->GetTrueVSize();

    for (int i = 0; i < npts; i++) {
        Vector coord(dim_);
        fault_geom_->GetPointCoordinates(i, coord);

        // Get initial slip from friction law
        Vector Sinit(num_slip_components_);
        friction_->GetInitialSlip(coord, Sinit);

        // Store in state vector
        for (int c = 0; c < num_slip_components_; c++) {
            state(i * num_state_components_ + c) = Sinit(c);
        }
        // Leave psi uninitialized (computed in Init)
    }
}

void FaultOperator::Init(real_t time, const Vector &traction, Vector &state)
{
    const int npts = fault_fes_->GetTrueVSize();
    VMax_ = 0.0;

    for (int i = 0; i < npts; i++) {
        // Extract traction components
        real_t sigma_n = traction(i * num_traction_components_ + 0);
        Vector tau(num_slip_components_);
        for (int c = 0; c < num_slip_components_; c++) {
            tau(c) = traction(i * num_traction_components_ + 1 + c);
        }

        // Get coordinates for spatially varying parameters
        Vector coord(dim_);
        fault_geom_->GetPointCoordinates(i, coord);

        // Compute initial state variable
        real_t psi = friction_->ComputeInitialState(coord, sigma_n, tau);
        state(i * num_state_components_ + num_slip_components_) = psi;

        // Compute initial slip rate for output
        Vector V(num_slip_components_);
        friction_->ComputeSlipRate(coord, sigma_n, tau, psi, V);
        VMax_ = std::max(VMax_, V.Norml2());
    }
}

real_t FaultOperator::ComputeRHS(real_t time,
                                  const Vector &traction,
                                  const Vector &state,
                                  Vector &rate)
{
    const int npts = fault_fes_->GetTrueVSize();
    VMax_ = 0.0;

    for (int i = 0; i < npts; i++) {
        // Extract current state
        Vector slip(num_slip_components_);
        for (int c = 0; c < num_slip_components_; c++) {
            slip(c) = state(i * num_state_components_ + c);
        }
        real_t psi = state(i * num_state_components_ + num_slip_components_);

        // Extract traction
        real_t sigma_n = traction(i * num_traction_components_ + 0);
        Vector tau(num_slip_components_);
        for (int c = 0; c < num_slip_components_; c++) {
            tau(c) = traction(i * num_traction_components_ + 1 + c);
        }

        // Get coordinates for spatially varying parameters
        Vector coord(dim_);
        fault_geom_->GetPointCoordinates(i, coord);

        // Compute slip rate from friction law
        Vector V(num_slip_components_);
        friction_->ComputeSlipRate(coord, sigma_n, tau, psi, V);
        real_t V_mag = V.Norml2();
        VMax_ = std::max(VMax_, V_mag);

        // Store slip rate: dS/dt = V
        for (int c = 0; c < num_slip_components_; c++) {
            rate(i * num_state_components_ + c) = V(c);
        }

        // Compute and store state rate: dψ/dt
        real_t psi_rate = friction_->ComputeStateRate(coord, V_mag, psi);
        rate(i * num_state_components_ + num_slip_components_) = psi_rate;
    }

    return VMax_;
}

} // namespace seas
} // namespace mfem
```

---

## 9. Detailed Implementation: RateStateFrictionLaw

### 9.1 Header File

```cpp
// File: miniapps/seas/friction/friction_law.hpp

#ifndef MFEM_SEAS_FRICTION_LAW_HPP
#define MFEM_SEAS_FRICTION_LAW_HPP

#include "mfem.hpp"

namespace mfem {
namespace seas {

/**
 * @brief Base class for rate-and-state friction laws
 */
class RateStateFrictionLaw {
public:
    virtual ~RateStateFrictionLaw() = default;

    /// Get initial slip at a point
    virtual void GetInitialSlip(const Vector &x, Vector &Sinit) = 0;

    /// Compute initial state variable from stress state
    virtual real_t ComputeInitialState(const Vector &x,
                                        real_t sigma_n,
                                        const Vector &tau) = 0;

    /// Compute slip rate from stress and state
    /// Solves: τ = σ_n · f(V, ψ) + η · V
    virtual void ComputeSlipRate(const Vector &x,
                                  real_t sigma_n,
                                  const Vector &tau,
                                  real_t psi,
                                  Vector &V) = 0;

    /// Compute state evolution rate
    virtual real_t ComputeStateRate(const Vector &x,
                                     real_t V_mag,
                                     real_t psi) = 0;
};

} // namespace seas
} // namespace mfem

#endif
```

### 9.2 Dieterich-Ruina Implementation

```cpp
// File: miniapps/seas/friction/dieterich_ruina.hpp

#ifndef MFEM_SEAS_DIETERICH_RUINA_HPP
#define MFEM_SEAS_DIETERICH_RUINA_HPP

#include "friction_law.hpp"

namespace mfem {
namespace seas {

/**
 * @brief Dieterich-Ruina ageing law friction
 *
 * Friction coefficient:
 *   f(V, ψ) = a · asinh((V/(2V₀)) · exp(ψ/a))
 *
 * State evolution (ageing law):
 *   dψ/dt = (b/L) · (V₀ · exp((f₀ - ψ)/b) - V)
 *
 * Stress balance with radiation damping:
 *   τ = σ_n · f(V, ψ) + η · V
 */
class DieterichRuinaAgeing : public RateStateFrictionLaw {
public:
    /// Constant parameters (uniform over fault)
    struct ConstantParams {
        real_t V0 = 1.0e-6;     // Reference velocity [m/s]
        real_t b = 0.015;       // State evolution parameter
        real_t f0 = 0.6;        // Reference friction coefficient
    };

    DieterichRuinaAgeing(const ConstantParams &cp);

    /// Set spatially varying parameter 'a'
    void SetACoefficient(Coefficient *a_coeff) { a_coeff_ = a_coeff; }

    /// Set radiation damping coefficient η
    void SetEtaCoefficient(Coefficient *eta_coeff) { eta_coeff_ = eta_coeff; }

    /// Set characteristic slip distance L
    void SetLCoefficient(Coefficient *L_coeff) { L_coeff_ = L_coeff; }

    /// Set normal pre-stress
    void SetNormalPreStress(Coefficient *sn_pre_coeff) { sn_pre_coeff_ = sn_pre_coeff; }

    /// Set shear pre-stress
    void SetShearPreStress(VectorCoefficient *tau_pre_coeff) { tau_pre_coeff_ = tau_pre_coeff; }

    /// Set initial slip velocity
    void SetInitialVelocity(VectorCoefficient *Vinit_coeff) { Vinit_coeff_ = Vinit_coeff; }

    /// Set initial slip
    void SetInitialSlip(VectorCoefficient *Sinit_coeff) { Sinit_coeff_ = Sinit_coeff; }

    // RateStateFrictionLaw interface
    void GetInitialSlip(const Vector &x, Vector &Sinit) override;
    real_t ComputeInitialState(const Vector &x, real_t sigma_n, const Vector &tau) override;
    void ComputeSlipRate(const Vector &x, real_t sigma_n, const Vector &tau,
                          real_t psi, Vector &V) override;
    real_t ComputeStateRate(const Vector &x, real_t V_mag, real_t psi) override;

protected:
    ConstantParams cp_;

    Coefficient *a_coeff_ = nullptr;
    Coefficient *eta_coeff_ = nullptr;
    Coefficient *L_coeff_ = nullptr;
    Coefficient *sn_pre_coeff_ = nullptr;
    VectorCoefficient *tau_pre_coeff_ = nullptr;
    VectorCoefficient *Vinit_coeff_ = nullptr;
    VectorCoefficient *Sinit_coeff_ = nullptr;

    /// Evaluate coefficient at point (handles null coefficients)
    real_t EvalCoeff(Coefficient *c, const Vector &x, real_t default_val);
    void EvalVectorCoeff(VectorCoefficient *c, const Vector &x, Vector &val);

    /// Friction strength function: F(V, ψ) = σ_n · f(V, ψ)
    real_t FrictionStrength(real_t a, real_t sigma_n, real_t V, real_t psi);

    /// Inverse friction (no radiation damping)
    real_t FrictionInverse(real_t a, real_t sigma_n, real_t tau_mag, real_t psi);

    /// Newton solver for slip rate with radiation damping
    real_t SolveSlipRate(real_t a, real_t eta, real_t sigma_n,
                          real_t tau_mag, real_t psi);
};

// Implementation
inline real_t DieterichRuinaAgeing::FrictionStrength(
    real_t a, real_t sigma_n, real_t V, real_t psi)
{
    real_t e = std::exp(psi / a);
    real_t f = a * std::asinh((V / (2.0 * cp_.V0)) * e);
    return sigma_n * f;
}

inline real_t DieterichRuinaAgeing::FrictionInverse(
    real_t a, real_t sigma_n, real_t tau_mag, real_t psi)
{
    real_t r = tau_mag / sigma_n;
    return cp_.V0 * (std::exp((r - psi) / a) - std::exp(-(r + psi) / a));
}

real_t DieterichRuinaAgeing::SolveSlipRate(
    real_t a, real_t eta, real_t sigma_n, real_t tau_mag, real_t psi)
{
    if (eta == 0.0) {
        return FrictionInverse(a, sigma_n, tau_mag, psi);
    }

    if (sigma_n <= 0.0) {
        // Fault in tension: viscous sliding
        return tau_mag / eta;
    }

    // Newton-Raphson solver for: R(V) = τ - F(V,ψ) - η·V = 0
    real_t V_lo = 0.0;
    real_t V_hi = tau_mag / eta;

    // Brent's method or bisection
    const int max_iter = 50;
    const real_t tol = 1.0e-12;

    for (int iter = 0; iter < max_iter; iter++) {
        real_t V_mid = 0.5 * (V_lo + V_hi);
        real_t R = tau_mag - FrictionStrength(a, sigma_n, V_mid, psi) - eta * V_mid;

        if (std::abs(R) < tol || (V_hi - V_lo) < tol * V_mid) {
            return V_mid;
        }

        if (R > 0) {
            V_lo = V_mid;
        } else {
            V_hi = V_mid;
        }
    }

    return 0.5 * (V_lo + V_hi);
}

void DieterichRuinaAgeing::ComputeSlipRate(
    const Vector &x, real_t sigma_n, const Vector &tau, real_t psi, Vector &V)
{
    int dim = tau.Size();

    // Get spatially varying parameters
    real_t a = EvalCoeff(a_coeff_, x, 0.01);
    real_t eta = EvalCoeff(eta_coeff_, x, 0.0);
    real_t sn_pre = EvalCoeff(sn_pre_coeff_, x, 0.0);

    Vector tau_pre(dim);
    EvalVectorCoeff(tau_pre_coeff_, x, tau_pre);

    // Compute absolute stresses
    real_t sigma_n_abs = sn_pre - sigma_n;  // Compression positive

    Vector tau_abs(dim);
    add(tau, tau_pre, tau_abs);
    real_t tau_mag = tau_abs.Norml2();

    // Solve for slip rate magnitude
    real_t V_mag = SolveSlipRate(a, eta, sigma_n_abs, tau_mag, psi);

    // Slip rate in shear stress direction
    if (tau_mag > 1.0e-14) {
        V.Set(-V_mag / tau_mag, tau_abs);
    } else {
        V = 0.0;
    }
}

real_t DieterichRuinaAgeing::ComputeStateRate(
    const Vector &x, real_t V_mag, real_t psi)
{
    real_t L = EvalCoeff(L_coeff_, x, 0.008);

    // dψ/dt = (b/L) · (V₀ · exp((f₀ - ψ)/b) - V)
    return (cp_.b / L) * (cp_.V0 * std::exp((cp_.f0 - psi) / cp_.b) - V_mag);
}

real_t DieterichRuinaAgeing::ComputeInitialState(
    const Vector &x, real_t sigma_n, const Vector &tau)
{
    int dim = tau.Size();

    real_t a = EvalCoeff(a_coeff_, x, 0.01);
    real_t eta = EvalCoeff(eta_coeff_, x, 0.0);
    real_t sn_pre = EvalCoeff(sn_pre_coeff_, x, 0.0);

    Vector tau_pre(dim), Vinit(dim);
    EvalVectorCoeff(tau_pre_coeff_, x, tau_pre);
    EvalVectorCoeff(Vinit_coeff_, x, Vinit);

    real_t sigma_n_abs = sn_pre - sigma_n;
    Vector tau_abs(dim);
    add(tau, tau_pre, tau_abs);
    real_t tau_mag = tau_abs.Norml2();
    real_t V_init = Vinit.Norml2();

    if (V_init < 1.0e-14) {
        return cp_.f0;
    }

    // Solve for ψ: τ = σ_n · f(V_init, ψ) + η · V_init
    // sinh((τ - η·V)/(a·σ_n)) = (V/(2V₀)) · exp(ψ/a)
    // ψ = a · ln((2V₀/V) · sinh((τ - η·V)/(a·σ_n)))

    real_t s = std::sinh((tau_mag - eta * V_init) / (a * sigma_n_abs));
    real_t l = std::log((2.0 * cp_.V0 / V_init) * s);
    return a * l;
}

} // namespace seas
} // namespace mfem

#endif
```

---

## 10. Detailed Implementation: Fault Traction Integrator

### 10.1 Face Integrator for Traction Computation

```cpp
// File: miniapps/seas/integrators/fault_traction_integrator.hpp

#ifndef MFEM_SEAS_FAULT_TRACTION_INTEGRATOR_HPP
#define MFEM_SEAS_FAULT_TRACTION_INTEGRATOR_HPP

#include "mfem.hpp"

namespace mfem {
namespace seas {

/**
 * @brief Computes traction on fault faces from domain displacement
 *
 * Given displacement field u on the domain, computes:
 *   T = σ(u) · n
 * where σ = λ tr(ε) I + 2μ ε is the Cauchy stress.
 *
 * Result is stored in a fault-based vector:
 *   [sigma_n, tau_1, tau_2, ...] at each integration point
 */
class FaultTractionIntegrator {
public:
    FaultTractionIntegrator(Coefficient &lambda, Coefficient &mu,
                            const Array<int> &fault_faces);

    /// Compute traction on all fault faces
    /// @param u_gf    Domain displacement GridFunction
    /// @param traction Output vector (sized for fault DOFs)
    void ComputeTraction(const GridFunction &u_gf, Vector &traction);

    /// Get number of traction DOFs
    int GetTractionSize() const;

protected:
    Coefficient *lambda_, *mu_;
    Array<int> fault_faces_;  // List of face indices that are faults

    Mesh *mesh_;
    FiniteElementSpace *fes_;

    // Cached data
    int dim_;
    int num_fault_points_;

    // Work arrays
    DenseMatrix dshape_, stress_, grad_u_;
    Vector shape_, nor_;
};

void FaultTractionIntegrator::ComputeTraction(
    const GridFunction &u_gf, Vector &traction)
{
    const FiniteElementSpace *fes = u_gf.FESpace();
    Mesh *mesh = fes->GetMesh();

    int traction_offset = 0;

    for (int fi = 0; fi < fault_faces_.Size(); fi++) {
        int face_idx = fault_faces_[fi];

        // Get face transformation
        FaceElementTransformations *Trans =
            mesh->GetInteriorFaceTransformations(face_idx);

        if (Trans == nullptr) {
            // Boundary face
            Trans = mesh->GetBdrFaceTransformations(face_idx);
        }

        // Get elements on both sides
        const FiniteElement &el1 = *fes->GetFE(Trans->Elem1No);
        const FiniteElement *el2 = (Trans->Elem2No >= 0) ?
            fes->GetFE(Trans->Elem2No) : nullptr;

        // Get element DOFs
        Array<int> vdofs1, vdofs2;
        fes->GetElementVDofs(Trans->Elem1No, vdofs1);
        if (el2) {
            fes->GetElementVDofs(Trans->Elem2No, vdofs2);
        }

        // Extract element displacement vectors
        Vector u1, u2;
        u_gf.GetSubVector(vdofs1, u1);
        if (el2) {
            u_gf.GetSubVector(vdofs2, u2);
        }

        // Quadrature rule for face
        int order = 2 * el1.GetOrder() + 1;
        const IntegrationRule &ir = IntRules.Get(Trans->GetGeometryType(), order);

        for (int q = 0; q < ir.GetNPoints(); q++) {
            const IntegrationPoint &ip = ir.IntPoint(q);

            // Map face point to element points
            Trans->SetAllIntPoints(&ip);
            const IntegrationPoint &eip1 = Trans->GetElement1IntPoint();

            // Compute stress at this point
            ComputeStressAtPoint(el1, *Trans->Elem1, u1, eip1, stress_);

            // Get face normal
            CalcOrtho(Trans->Jacobian(), nor_);
            nor_ /= nor_.Norml2();  // Unit normal

            // Traction: T = σ · n
            Vector T(dim_);
            stress_.Mult(nor_, T);

            // Store in output vector
            // Component 0: normal traction (T · n)
            traction(traction_offset++) = T * nor_;

            // Components 1 to dim-1: shear traction (tangential)
            // Project T onto tangent plane: T_shear = T - (T·n)n
            Vector T_shear(dim_);
            add(T, -(T * nor_), nor_, T_shear);

            // For 2D: single shear component
            // For 3D: two shear components (need tangent basis)
            if (dim_ == 2) {
                // Tangent = perpendicular to normal
                traction(traction_offset++) = T_shear(0) * nor_(1) - T_shear(1) * nor_(0);
            } else {
                // 3D: need to construct tangent basis
                // This requires fault geometry information
                // For now, store full tangential vector
                for (int d = 0; d < dim_ - 1; d++) {
                    traction(traction_offset++) = T_shear(d);  // Simplified
                }
            }
        }
    }
}

void FaultTractionIntegrator::ComputeStressAtPoint(
    const FiniteElement &el,
    ElementTransformation &Trans,
    const Vector &u_el,
    const IntegrationPoint &ip,
    DenseMatrix &stress)
{
    int ndof = el.GetDof();

    // Evaluate shape function gradients
    el.CalcDShape(ip, dshape_);

    // Transform to physical space
    DenseMatrix Jinv(dim_);
    CalcInverse(Trans.Jacobian(), Jinv);
    DenseMatrix dshape_phys(ndof, dim_);
    Mult(dshape_, Jinv, dshape_phys);

    // Compute displacement gradient
    grad_u_.SetSize(dim_, dim_);
    grad_u_ = 0.0;

    for (int i = 0; i < ndof; i++) {
        for (int d1 = 0; d1 < dim_; d1++) {
            for (int d2 = 0; d2 < dim_; d2++) {
                grad_u_(d1, d2) += u_el(i + d1 * ndof) * dshape_phys(i, d2);
            }
        }
    }

    // Compute strain: ε = (∇u + ∇u^T) / 2
    DenseMatrix strain(dim_);
    for (int i = 0; i < dim_; i++) {
        for (int j = 0; j < dim_; j++) {
            strain(i, j) = 0.5 * (grad_u_(i, j) + grad_u_(j, i));
        }
    }

    // Compute stress: σ = λ tr(ε) I + 2μ ε
    real_t lam = lambda_->Eval(Trans, ip);
    real_t mu = mu_->Eval(Trans, ip);
    real_t trace_eps = strain.Trace();

    stress.SetSize(dim_);
    for (int i = 0; i < dim_; i++) {
        for (int j = 0; j < dim_; j++) {
            stress(i, j) = 2.0 * mu * strain(i, j);
            if (i == j) {
                stress(i, j) += lam * trace_eps;
            }
        }
    }
}

} // namespace seas
} // namespace mfem

#endif
```

---

## 11. Detailed Implementation: QuasiDynamicOperator

### 11.1 Complete Quasi-Dynamic Solver

```cpp
// File: miniapps/seas/operators/seas_operator.hpp

#ifndef MFEM_SEAS_OPERATOR_HPP
#define MFEM_SEAS_OPERATOR_HPP

#include "mfem.hpp"
#include "fault_operator.hpp"
#include "../integrators/fault_traction_integrator.hpp"

namespace mfem {
namespace seas {

/**
 * @brief Quasi-dynamic SEAS operator
 *
 * Implements the time-dependent operator:
 *   d(state)/dt = RHS(t, state)
 *
 * where state = [slip, psi] and RHS involves:
 *   1. Solving static elasticity with slip BC
 *   2. Computing traction from displacement
 *   3. Computing slip rate and state rate from friction law
 */
class QuasiDynamicOperator : public TimeDependentOperator {
public:
    QuasiDynamicOperator(FiniteElementSpace *domain_fes,
                         FaultOperator *fault_op,
                         FaultTractionIntegrator *traction_integrator,
                         Coefficient *lambda,
                         Coefficient *mu,
                         const Array<int> &fault_bdr_attr);

    ~QuasiDynamicOperator();

    /// Initialize state at t=0
    void InitialCondition(Vector &state);

    /// Compute d(state)/dt given current state
    void Mult(const Vector &state, Vector &rate) const override;

    /// Get current displacement solution
    const GridFunction &GetDisplacement() const { return *u_gf_; }

    /// Get traction on fault
    const Vector &GetTraction() const { return traction_; }

protected:
    // Domain components
    FiniteElementSpace *domain_fes_;
    BilinearForm *a_form_;        // Elasticity stiffness
    LinearForm *b_form_;          // RHS (with slip BC contribution)
    GridFunction *u_gf_;          // Displacement solution

    // Fault components
    FaultOperator *fault_op_;
    FaultTractionIntegrator *traction_integrator_;

    // Material coefficients
    Coefficient *lambda_, *mu_;

    // Boundary attributes for fault
    Array<int> fault_bdr_attr_;

    // Linear solver
    CGSolver *solver_;
    GSSmoother *prec_;
    SparseMatrix *A_;

    // Work vectors
    mutable Vector traction_;
    mutable Vector b_;
    mutable Vector X_, B_;

    /// Solve domain problem with given slip BC
    void SolveDomain(const Vector &state) const;

    /// Update traction from current displacement
    void UpdateTraction() const;

    /// Apply slip as boundary condition
    void ApplySlipBC(const Vector &state, Vector &rhs) const;
};

QuasiDynamicOperator::QuasiDynamicOperator(
    FiniteElementSpace *domain_fes,
    FaultOperator *fault_op,
    FaultTractionIntegrator *traction_integrator,
    Coefficient *lambda,
    Coefficient *mu,
    const Array<int> &fault_bdr_attr)
    : TimeDependentOperator(fault_op->GetStateSize()),
      domain_fes_(domain_fes),
      fault_op_(fault_op),
      traction_integrator_(traction_integrator),
      lambda_(lambda),
      mu_(mu)
{
    fault_bdr_attr_ = fault_bdr_attr;

    // Create bilinear form for elasticity
    a_form_ = new BilinearForm(domain_fes_);
    a_form_->AddDomainIntegrator(new ElasticityIntegrator(*lambda_, *mu_));

    // Add DG interior face integrator
    real_t alpha = -1.0;  // SIPG
    real_t kappa = 10.0;  // Penalty (adjust based on mesh)
    a_form_->AddInteriorFaceIntegrator(
        new DGElasticityIntegrator(*lambda_, *mu_, alpha, kappa));

    // Add DG boundary face integrator on fault
    a_form_->AddBdrFaceIntegrator(
        new DGElasticityIntegrator(*lambda_, *mu_, alpha, kappa),
        fault_bdr_attr_);

    a_form_->Assemble();
    a_form_->Finalize();
    A_ = &a_form_->SpMat();

    // Create linear form for RHS
    b_form_ = new LinearForm(domain_fes_);

    // Displacement GridFunction
    u_gf_ = new GridFunction(domain_fes_);
    *u_gf_ = 0.0;

    // Linear solver
    prec_ = new GSSmoother(*A_);
    solver_ = new CGSolver();
    solver_->SetRelTol(1e-12);
    solver_->SetAbsTol(1e-14);
    solver_->SetMaxIter(1000);
    solver_->SetOperator(*A_);
    solver_->SetPreconditioner(*prec_);

    // Allocate work vectors
    traction_.SetSize(fault_op_->GetTractionSize());
    b_.SetSize(domain_fes_->GetTrueVSize());
    X_.SetSize(domain_fes_->GetTrueVSize());
    B_.SetSize(domain_fes_->GetTrueVSize());
}

void QuasiDynamicOperator::InitialCondition(Vector &state)
{
    // Step 1: Set initial slip
    fault_op_->PreInit(state);

    // Step 2: Solve domain with initial slip BC
    SolveDomain(state);

    // Step 3: Compute initial traction
    UpdateTraction();

    // Step 4: Initialize state variable
    fault_op_->Init(0.0, traction_, state);
}

void QuasiDynamicOperator::Mult(const Vector &state, Vector &rate) const
{
    // Step 1: Solve domain equilibrium with current slip BC
    SolveDomain(state);

    // Step 2: Compute traction from displacement
    UpdateTraction();

    // Step 3: Compute friction rates
    fault_op_->ComputeRHS(GetTime(), traction_, state, rate);
}

void QuasiDynamicOperator::SolveDomain(const Vector &state) const
{
    // Reset RHS
    b_ = 0.0;

    // Apply slip boundary condition contribution to RHS
    ApplySlipBC(state, b_);

    // Solve A*u = b
    B_ = b_;
    X_ = 0.0;
    solver_->Mult(B_, X_);

    // Copy solution to GridFunction
    u_gf_->SetFromTrueDofs(X_);
}

void QuasiDynamicOperator::UpdateTraction() const
{
    traction_integrator_->ComputeTraction(*u_gf_, traction_);
}

void QuasiDynamicOperator::ApplySlipBC(const Vector &state, Vector &rhs) const
{
    // Create slip coefficient from state
    // This requires extracting slip from state and creating a VectorCoefficient

    // For each fault boundary face, add DG Dirichlet contribution
    // using DGElasticityDirichletLFIntegrator pattern

    Mesh *mesh = domain_fes_->GetMesh();
    int dim = mesh->Dimension();

    // ... (implementation details depend on state layout and fault geometry)
    // Key: Use b_form_->AddBdrFaceIntegrator() pattern
}

} // namespace seas
} // namespace mfem

#endif
```

---

## 12. Implementation Roadmap Summary

### Phase 1: Core Infrastructure (Week 1-2)

| Task | Files | Dependencies |
|------|-------|--------------|
| 1.1 Create friction law base class | `friction/friction_law.hpp` | None |
| 1.2 Implement Dieterich-Ruina | `friction/dieterich_ruina.hpp` | 1.1 |
| 1.3 Create fault geometry class | `geometry/fault_geometry.hpp` | None |
| 1.4 Create FaultOperator | `operators/fault_operator.hpp` | 1.1, 1.3 |

### Phase 2: Domain-Fault Coupling (Week 3-4)

| Task | Files | Dependencies |
|------|-------|--------------|
| 2.1 Implement traction integrator | `integrators/fault_traction_integrator.hpp` | None |
| 2.2 Implement slip BC integrator | `integrators/fault_slip_integrator.hpp` | None |
| 2.3 Create QuasiDynamicOperator | `operators/seas_operator.hpp` | 2.1, 2.2, 1.4 |
| 2.4 Test with simple fault geometry | `tests/` | 2.3 |

### Phase 3: Verification (Week 5-6)

| Task | Files | Dependencies |
|------|-------|--------------|
| 3.1 Implement BP1 benchmark | `examples/bp1.cpp` | 2.3 |
| 3.2 Compare with reference data | `bp2/compare_results.py` | 3.1 |
| 3.3 Implement BP2 benchmark | `examples/bp2.cpp` | 3.1 |
| 3.4 Performance optimization | All | 3.3 |

### Critical Path

```
friction_law.hpp → dieterich_ruina.hpp → fault_operator.hpp
                                              ↓
fault_geometry.hpp ─────────────────────────→ seas_operator.hpp
                                              ↑
fault_traction_integrator.hpp ───────────────┘
```

---

## 13. Key Differences and Adaptations

| Tandem | MFEM Adaptation |
|--------|-----------------|
| BlockVector by element | Standard Vector with strided access |
| SFINAE-based dispatch | Virtual function overrides |
| Lua configuration | Coefficient-based parameters |
| PETSc linear solver | MFEM CGSolver or hypre |
| Custom scatter | MFEM ParGridFunction exchange |
| Template compile-time | Runtime polymorphism |

---

## 14. Testing Strategy

### Unit Tests

1. **Friction law**: Verify slip rate solver against analytical solutions
2. **State evolution**: Check steady-state behavior
3. **Traction computation**: Compare with analytical stress fields

### Integration Tests

1. **Simple shear**: Uniform fault with known analytical solution
2. **BP1 benchmark**: Compare slip evolution with reference
3. **BP2 benchmark**: Full 3D verification

### Performance Tests

1. **Scaling**: Measure time vs. problem size
2. **Solver iterations**: Track convergence
3. **Memory usage**: Monitor allocations

---

*Document Version: 1.0*
*Generated: February 5, 2026*
