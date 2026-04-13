# SEAS Miniapp Architecture

> Post-refactor (7 phases). Last updated: 2026-04-11.
>
> For detailed physics derivations and weak-form equations, see `CODEBASE_GUIDE.md`.
> For critical numerical constraints and debug history, see `CLAUDE.md`.

MFEM-based Discontinuous Galerkin code for SCEC SEAS (Sequences of Earthquakes and Aseismic Slip) benchmark problems. Couples a DG elasticity domain solver with rate-and-state friction on an embedded fault interface via a quasi-dynamic approximation. Supports 2D antiplane (BP1/BP2) and 3D full elasticity (BP5).

**Code size:** ~60K LOC C++ (excluding extern/toml11).

---

## Module Inventory

| Module | Directory | Key Files | Purpose |
|--------|-----------|-----------|---------|
| **Config** | `config/` | `seas_config.hpp`, `seas_config_parser.hpp`, `seas_config_bridge.hpp`, `bp{1,2,5}_params.hpp` | TOML config parsing, benchmark parameters, spatial functions a(x), Dc(x) |
| **Constitutive** | `constitutive/` | `constitutive_model.hpp`, `linear_elastic.hpp` | Material model interface + isotropic linear elastic |
| **Domain** | `domain/` | `domain_operator.hpp`, `elasticity_operator.hpp`, `antiplane_operator.hpp`, `boundary_config.hpp` | DG elasticity solvers (3D and antiplane), boundary conditions |
| **Fault** | `fault/` | `rate_state_fault.hpp`, `fault_geometry.hpp`, `fault_basis.hpp`, `face_quadrature.hpp` | Fault ODE operator, DOF management, coordinate transforms |
| **Friction** | `friction/` | `dieterich_ruina.hpp`, `friction_law.hpp`, `state_evolution.hpp` | Rate-and-state friction (Brent solver), aging/slip laws |
| **Integrator** | `integrator/` | `dg_elasticity_br2_integrator.hpp`, `dg_elasticity_ip_combined_integrator.hpp` | DG bilinear form integrators: BR2 and IP methods |
| **Solver** | `solver/` | `seas_operator.hpp`, `time_stepper.hpp`, `seas_bdrload_operator.hpp` | SEAS ODE coupling operator, Dormand-Prince RK45 |
| **I/O** | `io/` | `bp5_parallel_output.hpp`, `probe_output.hpp`, `checkpoint.hpp`, `paraview_output.hpp` | SCEC benchmark output, checkpointing, VTK visualization |
| **Common** | `common/` | `mpi_context.hpp`, `parallel_utils.hpp`, `fault_scatter.hpp`, `logging.hpp` | MPI wrappers, parallel utilities, logging |
| **Trace** | `trace/` | `face_trace_logger.hpp` | Face-level traction diagnostics |
| **Driver** | `drivers/` | `seas_driver.cpp` | Primary BP5 parallel driver (TOML-configured) |
| **Legacy** | root | `pseas.cpp` | Old BP2 driver (kept for reference) |
| **Tests** | `tests/` | `unit/` (38), `parallel/` (11), `verification/` (7) | Unit, MPI parallel, and full benchmark tests |

---

## Directory Layout

```
miniapps/seas/
├── drivers/
│   └── seas_driver.cpp             ← PRIMARY ENTRY POINT (TOML-based BP5 driver)
│
├── config/                         ← Configuration layer
│   ├── seas_config.hpp                 SEASConfig: top-level config aggregate
│   ├── seas_config_parser.hpp          TOML file parser
│   ├── seas_config_bridge.hpp          Config → domain/solver object bridge
│   ├── bp5_params.hpp                  BP5 parameters + spatial functions a(x2,x3), Dc(x2,x3)
│   ├── bp5_mesh_utils.hpp              BP5 mesh loading utilities
│   ├── bp2_params.hpp                  BP2 parameters + a(z)
│   └── bp1_params.hpp                  BP1 parameters
│
├── constitutive/                   ← Material models
│   ├── constitutive_model.hpp          Abstract: ComputeStress(), ComputeTangent()
│   └── linear_elastic.hpp              σ = λ tr(ε)I + 2με
│
├── domain/                         ← Domain PDE solvers
│   ├── domain_operator.hpp             Abstract base: Solve(), ComputeTraction()
│   ├── elasticity_operator.hpp         3D DG elasticity (BR2/IP), stiffness assembly, traction
│   ├── antiplane_operator.hpp          2D scalar DG (BP1/BP2)
│   ├── antiplane_bdrload_operator.hpp  Boundary-load variant (BP1)
│   ├── boundary_config.hpp             BC specification (Dirichlet/Natural/Fault attrs)
│   ├── domain_config.hpp               Solver tuning (penalty, BLR tol, face basis)
│   ├── bp2_mesh.hpp                    BP2 graded mesh generator
│   └── seas_boundary_tags.hpp          Gmsh physical surface tag definitions
│
├── integrator/                     ← DG bilinear form integrators
│   ├── dg_elasticity_br2_integrator.hpp        BR2 interior + boundary (3D vector)
│   ├── dg_elasticity_ip_combined_integrator.hpp IP interior + boundary + traction
│   ├── dg_elasticity_ip_penalty_integrator.hpp  IP penalty with full elasticity tensor
│   └── dg_br2_integrator.hpp                    BR2 scalar (antiplane)
│
├── friction/                       ← Friction laws
│   ├── friction_law.hpp                Abstract: FrictionCoefficient(), SolveSlipRate()
│   ├── dieterich_ruina.hpp             Regularized DR, Brent solver, psi-space
│   └── state_evolution.hpp             AgingLaw, SlipLaw, AgingLawPsi, SlipLawPsi
│
├── fault/                          ← Fault interface coupling
│   ├── fault_geometry.hpp              Fault DOF management, spatial params a(x), η, τ_pre
│   ├── rate_state_fault.hpp            Fault ODE operator: traction → (V, dψ/dt)
│   ├── fault_basis.hpp                 Per-face coordinate frames (normal, dip, strike)
│   ├── face_quadrature.hpp             Multi-DOF face basis for fault discretization
│   └── fault_nodes.hpp                 Fault node extraction and DOF mapping
│
├── solver/                         ← Top-level coupling & time integration
│   ├── seas_operator.hpp               SEAS quasi-dynamic ODE operator (domain + fault)
│   ├── time_stepper.hpp                Dormand-Prince RK45 + AdaptiveTimeStepper
│   └── seas_bdrload_operator.hpp       Boundary-load coupling (BP1)
│
├── io/                             ← Output & checkpointing
│   ├── bp5_parallel_output.hpp         Distributed SCEC fltst probe output
│   ├── bp5_benchmark_output.hpp        BP5-specific SCEC format
│   ├── benchmark_output.hpp            Serial benchmark output
│   ├── parallel_benchmark_output.hpp   Parallel output with gather + dedup
│   ├── probe_output.hpp                Single-station file writer
│   ├── paraview_output.hpp             VTK/ParaView diagnostic output
│   └── checkpoint.hpp                  Binary checkpoint/restart
│
├── common/                         ← Shared utilities
│   ├── mpi_context.hpp                 MPI RAII wrapper (GlobalMax, GlobalSum, etc.)
│   ├── parallel_utils.hpp              MPI scatter/gather helpers
│   ├── fault_scatter.hpp               Distributed fault data communication
│   ├── logging.hpp                     Rank-aware logging (SEAS_LOG)
│   ├── mpi_check.hpp                   MPI error checking macros
│   ├── mpi_tags.hpp                    MPI message tag definitions
│   └── seas_types.hpp                  Enums, typedefs (BCMode, ProblemType, etc.)
│
├── trace/
│   └── face_trace_logger.hpp           Per-face traction diagnostics
│
├── tests/
│   ├── unit/           (~38 tests)     Fast unit tests (`make test`, ~2 min)
│   ├── parallel/       (~11 tests)     MPI parallel tests
│   └── verification/   (7 tests)       BP1/BP2/BP5 full benchmark simulations
│
├── extern/toml11/                      TOML config parsing (header-only)
├── pseas.cpp                           Legacy BP2 driver
├── Makefile                            Primary build system
├── CMakeLists.txt                      CMake alternative
├── ARCHITECTURE.md                     This file
├── CLAUDE.md                           Critical numerical constraints
└── CODEBASE_GUIDE.md                   Physics-to-numerics walkthrough
```

---

## Class Hierarchy

### Domain Layer

```
ConstitutiveModel                           [constitutive/constitutive_model.hpp]
│  Abstract: ComputeStress(), ComputeTangent(), GetMaxWaveSpeed()
└── LinearElastic                           [constitutive/linear_elastic.hpp]
      σ = λ tr(ε)I + 2με

DomainOperator<MeshType>                    [domain/domain_operator.hpp]
│  Abstract: Solve(), ComputeTraction(), GetFaultDOFs()
│  MeshType = Mesh (serial) | ParMesh (parallel)
├── ElasticityDomainOperator<MeshType>      [domain/elasticity_operator.hpp]
│     3D vector elasticity, DG BR2 or IP, ConstitutiveModel-based
└── AntiplaneDomainOperator<MeshType>       [domain/antiplane_operator.hpp]
      2D scalar Laplace, DG IP or BR2
```

### Friction Layer

```
FrictionLaw                                 [friction/friction_law.hpp]
│  Abstract: FrictionCoefficient(), SolveSlipRate()
└── DieterichRuinaFriction                  [friction/dieterich_ruina.hpp]
      Regularized DR with Brent solver
      Psi-space: SolveSlipRatePsi(), SolveSlipRateVectorPsi()

StateEvolution                              [friction/state_evolution.hpp]
│  Abstract: Rate(), SteadyState()
├── AgingLaw         dθ/dt = 1 - Vθ/Dc
├── SlipLaw          dθ/dt = -Vθ/Dc · ln(Vθ/Dc)
├── AgingLawPsi      dψ/dt = (bV₀/Dc)[exp((f₀-ψ)/b) - V/V₀]
└── SlipLawPsi       dψ/dt = -(V/Dc)(ψ - ψ_ss)
```

### Fault Layer

```
FaultGeometry<MeshType>                     [fault/fault_geometry.hpp]
  Precomputes: depths, a(x), Dc(x), τ_pre(x), η, V_init(x)

RateStateFaultOperator<MeshType, SlipComp>  [fault/rate_state_fault.hpp]
  SlipComp = 1 (BP2 scalar) | 2 (BP5 vector)
  ComputeRHS(): traction → slip rate + state evolution

FaultBasis                                  [fault/fault_basis.hpp]
  Per-face: (normal, tangent_dip, tangent_strike)
  Embedding: (dip,strike) → (X,Y,Z)   Projection: (X,Y,Z) → (dip,strike,σ_n)

FaceQuadrature                              [fault/face_quadrature.hpp]
  Multi-DOF face basis functions for fault discretization
```

### Solver Layer

```
TimeDependentOperator (MFEM)
└── SEASQuasiDynamicOperator<MeshType, DomainOpType, FaultOpType>
      │                                     [solver/seas_operator.hpp]
      │  Mult(): domain solve + traction + friction → ODE RHS
      │  SetInitialCondition(): 4-phase init
      │
      ├── domain_: DomainOpType*            (non-owning)
      ├── fault_:  FaultOpType*             (non-owning)
      ├── mpi_ctx_: MPIContext*
      ├── u_gf_:    GridFunction            (displacement solution)
      ├── slip_, traction_, local_slip_, local_traction_
      └── elastic_sigma_n_: bool            (normal stress feedback, default ON)

DormandPrinceRK45                           [solver/time_stepper.hpp]
  7-stage embedded RK pair (5th/4th order), FSAL
  PI controller for dt adaptation
  V-guard: reject stages with V > 100·V_stage0
  MPI: parallel error reduction via GlobalMax
```

### Type Aliases (Template Instantiations)

```cpp
// Parallel BP5 — primary production path
using PBP5SEASOp = SEASQuasiDynamicOperator<
    ParMesh,
    ElasticityDomainOperator<ParMesh>,
    RateStateFaultOperator<ParMesh, 2>>;

// Serial BP2
using SBP2SEASOp = SEASQuasiDynamicOperator<
    Mesh,
    AntiplaneDomainOperator<Mesh>,
    RateStateFaultOperator<Mesh, 1>>;
```

---

## Configuration System

Parameters flow: **TOML file → CLI overrides → SEASConfig → bridge → domain/solver objects**

```
SEASConfig                                  [config/seas_config.hpp]
├── MeshConfig          file, scale (1000.0), order (1)
├── MaterialConfig      rho (2670), cs (3464), nu (0.25) → mu(), lambda()
├── FrictionConfig      V0, f0, b, L0, L_nuc, a0, amax, sigma_n
├── LoadingConfig       Vp, V_init, V_nuc, delta_tau_factor, smooth_nucleation
├── FaultGeomConfig     Wf, lf, hs, ht, H, l_vw, w_nuc
├── BoundaryTomlConfig  dirichlet_attrs, natural_attrs, fault_attr
├── SolverConfig        dg_method (IP/BR2), solver_type, penalty_factor, blr_tol
├── TimeConfig          t_final, atol (1e-7), rtol (1e-50), max_steps, tandem_time_stepping
├── OutputConfig        output_dir, prefix, write_every_step, ref_dir
└── SimulationConfig    mode (qd / dynamic / hybrid)
```

Parsed by `seas_config_parser.hpp`, bridged to runtime objects by `seas_config_bridge.hpp`.

CLI overrides: `--mesh`, `--order`, `--V-nuc`, `--delta-tau-factor`, `--tandem-time-stepping`, etc.

---

## Execution Flow: Full Simulation (seas_driver.cpp)

```
┌────────────────────────────────────────────────────────────────────┐
│ Stage 1: CONFIGURATION                              (lines 50-96) │
│   Parse TOML + CLI → SEASConfig → BP5Params, BoundaryConfig, etc. │
├────────────────────────────────────────────────────────────────────┤
│ Stage 2: MESH                                      (lines 109-154)│
│   Load Gmsh file → Mesh → ParMesh (MPI partition)                 │
├────────────────────────────────────────────────────────────────────┤
│ Stage 3: DOMAIN OPERATOR                           (lines 156-172)│
│   LinearElastic(λ,μ) → ElasticityDomainOperator(mesh, BCs, DG)   │
│   Assembles stiffness K once (cached, reused all steps)           │
├────────────────────────────────────────────────────────────────────┤
│ Stage 4: FAULT COMPONENTS                          (lines 174-188)│
│   FaultGeometry + DieterichRuinaFriction + AgingLawPsi            │
│   → RateStateFaultOperator<ParMesh, 2>                            │
├────────────────────────────────────────────────────────────────────┤
│ Stage 5: SEAS OPERATOR + INIT                      (lines 190-204)│
│   PBP5SEASOp(domain, fault) + elastic σ_n feedback                │
│   4-phase init: PreInit → Solve → Init(equilibrium) → Verify     │
├────────────────────────────────────────────────────────────────────┤
│ Stage 6: I/O SETUP                                 (lines 206-241)│
│   ParallelBP5BenchmarkOutput (distributed probe ownership)        │
│   ProbeOutput (global V_max log, root rank only)                  │
├────────────────────────────────────────────────────────────────────┤
│ Stage 7: TIME STEPPER                              (lines 252-302)│
│   DormandPrinceRK45: atol=1e-7, rtol=1e-50                       │
│   dt_init = min(0.01·L_nuc/V_nuc, CFL), V-guard=100              │
├────────────────────────────────────────────────────────────────────┤
│ Stage 8: TIME LOOP                                 (lines 304-401)│
│   while (t < t_final && step < max_steps):                        │
│     ode_solver.Step(seas_op, state, t, dt)                        │
│     earthquake detection: V_max > 1e-3 m/s enter, < 1e-6 exit    │
│     adaptive output scheduling                                     │
├────────────────────────────────────────────────────────────────────┤
│ Stage 9: FINAL OUTPUT                              (lines 403-420)│
│   Force write, summary: steps, rejections, earthquake count       │
└────────────────────────────────────────────────────────────────────┘
```

---

## Execution Flow: One ODE RHS Evaluation (Mult)

Called 7 times per RK45 step. This is the core physics cycle:

```
DormandPrinceRK45::Step(seas_op, state, t, dt)    [time_stepper.hpp]
  │
  │  For each of 7 RK stages:
  │
  └── seas_op.Mult(state, rate)                    [seas_operator.hpp]
        │
        ├── 1. fault_->GetSlip(state, slip_)
        │       Extract slip from interleaved state vector
        │
        ├── 2. domain_->ExpandOwnedToLocalFault(slip_, local_slip_)
        │       MPI: broadcast owned DOFs → ghost DOFs
        │
        ├── 3. domain_->Solve(t, local_slip_, u_gf_)
        │       Assemble RHS: b = b_slip + b_dirichlet
        │       Reuse cached stiffness K
        │       Solve K*u = b via MUMPS/CG+AMG
        │
        ├── 4. domain_->ComputeTraction(u_gf_, local_slip_, local_traction_)
        │       DG flux: τ = {{σ(u)·n}} - penalty·([[u]] - slip)
        │       Project (X,Y,Z) → (dip, strike, σ_n)
        │
        ├── 5. domain_->RestrictToOwnedFault(local_traction_, traction_)
        │       Select owned subset from full local traction
        │
        └── 6. fault_->ComputeRHS(traction_, state, rate)
                For each owned fault DOF:
                  τ_total = τ_pre + traction
                  V = SolveSlipRateVectorPsi(τ, ψ, σ_n, η, a)  [Brent]
                  dψ/dt = aging_law(|V|, ψ, Dc)
```

### Data Flow Diagram

```
  state = [slip_dip, slip_strike, ψ] × N_fault
    │
    ▼
┌──────────────────┐
│ 1. GetSlip()     │  Extract slip from interleaved state vector
└────────┬─────────┘
         │ slip (fault-local: dip, strike)
         ▼
┌──────────────────┐
│ 2. EmbedSlipQP() │  Transform (dip,strike) → (X,Y,Z) via FaultBasis
└────────┬─────────┘
         │ slip (global 3D)
         ▼
┌──────────────────────────────────────────┐
│ 3. ElasticityDomainOperator::Solve()     │
│    K u = b_slip + b_dirichlet            │
│    (K cached, only RHS reassembled)      │
└────────┬─────────────────────────────────┘
         │ u (displacement field)
         ▼
┌──────────────────────────────────────────┐
│ 4. ComputeTraction()                     │
│    τ = {{σ(u)·n}} - penalty·([[u]]-slip) │
│    Project (X,Y,Z) → (dip, strike, σ_n) │
└────────┬─────────────────────────────────┘
         │ traction (fault-local)
         ▼
┌──────────────────────────────────────────┐
│ 5. ComputeRHS()                          │
│    τ_total = τ_pre + traction            │
│    V = Brent solve: τ = σ_n·f(V,ψ)+η·V  │
│    dψ/dt = aging_law(|V|, ψ, Dc)        │
└────────┬─────────────────────────────────┘
         ▼
  rate = [V_dip, V_strike, dψ/dt] × N_fault
```

---

## Data Ownership

| Data | Owner | Shape | Notes |
|------|-------|-------|-------|
| Mesh (`ParMesh`) | `main()` | Distributed | MFEM standard partitioning |
| Displacement `u_gf_` | `SEASQuasiDynamicOperator` | Full local mesh DOFs | `GridFunction` on DG FE space |
| State vector | `main()` | `[N_owned * state_per_node]` | ODE state: only owned fault DOFs |
| Slip `slip_` | `SEASQuasiDynamicOperator` | `[N_owned * slip_comps]` | Extracted from state |
| Local slip `local_slip_` | `SEASQuasiDynamicOperator` | `[N_local * slip_comps]` | Owned + ghosts after expand |
| Traction `traction_` | `SEASQuasiDynamicOperator` | `[N_owned * slip_comps]` | After restrict-to-owned |
| Stiffness `cached_Ah_` | `ElasticityDomainOperator` | Sparse matrix | Assembled once, reused |
| Fault geometry | `FaultGeometry` | Per-DOF arrays | Depths, a(z), eta, Dc, tau_pre |
| Friction constants | `DieterichRuinaFriction` | Global constants | V0, f0, b (Dc is per-DOF in BP5) |

### State Vector Layout

**BP2** (antiplane, 1 slip component, `state_per_node = 2`):
```
state = [slip_0, theta_0, slip_1, theta_1, ..., slip_{N-1}, theta_{N-1}]
         (or psi_i if use_psi_=true)
```

**BP5** (3D, 2 slip components, `state_per_node = 3`):
```
state = [slip_dip_0, slip_strike_0, psi_0, slip_dip_1, slip_strike_1, psi_1, ...]
```

Traction follows the same interleaving but without the state variable component.

---

## DG Formulation

### Weak Form (3D Elasticity)

Find u in V_h such that for all v in V_h:

```
sum_K integral_K sigma(u):epsilon(v) dx
  - sum_F integral_F ( {{sigma(u)}} . [[v]] + {{sigma(v)}} . [[u]] ) ds
  + sum_F integral_F penalty * [[u]] . [[v]] ds
  = integral_{Gamma_f} ( ... slip terms ... ) ds
  + integral_{Gamma_D} ( ... Dirichlet loading ... ) ds
```

where `sigma = lambda*tr(epsilon)*I + 2*mu*epsilon`, `{{.}}` = average, `[[.]]` = jump.

### SIPG Penalty (IP Method)

```
penalty = kappa * (p+1)^2 / h * max(lambda + 2*mu)
```

with `kappa` a user-tunable constant. The IP path applies the full elasticity tensor in the penalty term (not just a scalar).

### BR2 Method

Bassi-Rebay 2 replaces the standard penalty with a lifting operator:

```
L([u]) = M_K^{-1} * integral_e phi * [u] . n ds
```

The BR2 penalty term is `eta_BR2 * integral_K C:L([u]) : L([v]) dx`, where `eta_BR2 = n_faces_per_element` (3 for 2D triangles, 4 for 3D tets, etc.).

### Fault Interface Flux

On fault faces, the jump is modified: `[[u]] - delta` where `delta` is the prescribed slip. The traction extracted for the friction law is:

```
tau = {{sigma(u)}} . n - penalty * ([[u]] - delta)
```

This DG-corrected traction feeds into the friction solver as the elastic loading.

---

## Friction Coupling

### Regularized Dieterich-Ruina Law

```
f(V, psi) = a * asinh[ (V / 2V0) * exp(psi / a) ]
```

where `psi = f0 + b * ln(V0 * theta / Dc)` (logarithmic state variable).

### Stress Balance (Quasi-Dynamic)

```
tau_total = sigma_n * f(V, psi) + eta * V
```

where `tau_total = tau0 + tau_elastic`, `eta = mu / (2*cs)` is the radiation damping coefficient.

### Slip Rate Solver

Solve `tau_total = sigma_n * f(V, psi) + eta * V` for V using Brent's method:
- Bracket: `[0, tau/eta]` (monotonic, unique solution guaranteed)
- BP5 vector: solve for `|V|` via Brent in `log10(V)` space, then decompose `V_vec = (|V| / |tau|) * tau_vec`

### State Evolution

Aging law (psi-space): `dpsi/dt = (b*V0/Dc) * [exp((f0-psi)/b) - V/V0]`

Steady state: `psi_ss = f0 + b * ln(V0/V)` when `dpsi/dt = 0`.

---

## Time Integration

Dormand-Prince RK45 (7-stage, 5th-order, FSAL) with adaptive step control:

- Error norm: `max_i |err_i / (atol + rtol * |y_i|)|` (L-infinity, default)
- PI controller: `dt_new = safety * dt * err_norm^(-1/5)`
- Parallel: `MPI_Allreduce(err_norm, MPI_MAX)` ensures consistent accept/reject
- V-guard: reject RK stages where `V_max > 100 * V_stage0` (prevents runaway amplification)
- Settings: `atol = 1e-7`, `rtol = 1e-50` (pure absolute), `dt_min = 1e-6 s`, `dt_max = 0.1 yr`

---

## Boundary Conditions (BP5)

Controlled by `BoundaryConfig` (configured via TOML or `BCMode` enum):

| Mode | Dirichlet Faces | Natural Faces | Usage |
|------|----------------|---------------|-------|
| `FarField` (default) | attrs 1-4 (X=+-Lx, Y=+-Ly) | attrs 5-6 (Z=0, Z=-Lz) | Matches Tandem |
| `XOnly` | attrs 1-2 (X=+-Lx only) | 3-6 | Antiplane-style |
| `AllDirichlet` | attrs 1-6 (all) | none | Legacy (incorrect) |

Dirichlet loading: `u_X = sgn(Y) * Vp * t / 2` (along-strike plate motion).

Tandem mesh tag convention: Physical Surface 1 = Natural (top/bottom), 3 = Fault, 5 = Dirichlet (far-field).

---

## MPI Communication Pattern

Most computation is rank-local. MPI happens at these points:

| When | Operation | Pattern |
|------|-----------|---------|
| Per RK stage | `ExpandOwnedToLocalFault()` | Scatter: owned slip → ghost DOFs |
| Per RK stage | `domain_->Solve()` | MUMPS distributed solve (internal MFEM/Hypre) |
| Per RK stage | `fault_->ComputeRHS()` | Purely local (no MPI) |
| Per RK step | Error norm reduction | `MPI_Allreduce(MAX)` — consistent accept/reject |
| Per RK step | Max slip rate | `MPI_Allreduce(MAX)` — earthquake detection, dt control |
| At I/O setup | Probe station ownership | `MPI_Allreduce(MIN)` on distances (computed once) |
| At output | Station writes | Each rank writes only its owned stations |

**Critical:** The error norm `MPI_Allreduce(MAX)` prevents rank divergence; without it, ranks disagree on accept/reject and deadlock.

---

## Coordinate System (BP5)

```
SCEC Convention:        Tandem/MFEM Mesh:
  x1 = fault-normal  →  Y (fault at Y=0)
  x2 = along-strike  →  X
  x3 = depth (+down)  →  -Z (Z=0 surface, Z<0 depth)
```

Fault-local frame (FaultBasis):
- `tangent1` = dip direction = (0, 0, +1) (downward into earth)
- `tangent2` = strike direction = (1, 0, 0) (along-strike, dominant slip)
- `normal` = (0, 1, 0) or (0, -1, 0) depending on face orientation

---

## Initialization Sequence (4-Phase)

1. **PreInit**: slip=0, psi = f0 + b*ln(V0/V_init) (steady-state placeholder)
2. **First domain solve**: K*u=b with zero slip, compute initial elastic traction
3. **Init**: Compute psi(0) from stress equilibrium: `tau0 + tau_elastic = sigma_n * f(V_init, psi) + eta * V_init`
4. **Verification re-solve**: Re-solve domain, re-compute traction, verify equilibrium error < 1e-6

---

## BP5 Spatial Parameters

| Parameter | VW Core | VS Zone | Transition | Nucleation Zone |
|-----------|---------|---------|------------|-----------------|
| a | 0.004 | 0.04 | linear blend | 0.004 |
| b | 0.03 | 0.03 | 0.03 | 0.03 |
| Dc | 0.14 m | 0.14 m | 0.14 m | 0.13 m |
| sigma_n | 25 MPa | 25 MPa | 25 MPa | 25 MPa |
| V_init | V_init | V_init | V_init | V_nuc (configurable) |

VW core: `hs+ht <= depth <= hs+ht+H` AND `|x2| <= l_vw/2` (depth 4-16 km, +/-30 km strike).

Nucleation zone: left edge of VW, `-l_vw/2 <= x2 <= -l_vw/2 + w_nuc`, same depth range.

---

## Key Physics Equations Summary

| Equation | Formula | Code Location |
|----------|---------|---------------|
| Momentum balance | ∇·σ(u) = 0 | `elasticity_operator.hpp` Solve() |
| Constitutive law | σ = λ tr(ε)I + 2με | `linear_elastic.hpp` ComputeStress() |
| Fault stress balance | τ = σ_n·f(V,ψ) + η·V | `dieterich_ruina.hpp` SolveSlipRatePsi() |
| Friction coefficient | f = a·asinh[(V/2V₀)·exp(ψ/a)] | `dieterich_ruina.hpp` FrictionCoefficientPsi() |
| Aging law (ψ-space) | dψ/dt = (bV₀/Dc)[exp((f₀-ψ)/b) - V/V₀] | `state_evolution.hpp` AgingLawPsi::Rate() |
| Radiation damping | η = μ/(2c_s) | `fault_geometry.hpp` constructor |
| DG traction | τ = {{σ·n}} - penalty·([[u]] - slip) | `elasticity_operator.hpp` ComputeTraction() |

---

## Building and Testing

```bash
conda activate mfem-dev              # Build environment

# Build
make all                              # Everything
make seas_driver                      # Primary BP5 driver

# Unit tests (fast, ~2 min)
make test

# Parallel tests
mpirun -np 4 seas_test_bp5_parallel_smoke
mpirun -np 8 seas_test_parallel_elasticity

# Full verification (long-running)
mpirun -np 8 seas_bp5_full --config bp5/config/bp5_1000m.toml
```

---

## Known Limitations

1. **Matrix assembly**: Stiffness matrix assembled and stored (not matrix-free). Memory-bound for large meshes at high polynomial order.
2. **Direct solver default**: MUMPS (or MUMPS_BLR) for BP5. Iterative solvers (CG+AMG) available but less robust for DG.
3. **Single fault plane**: Code assumes one planar fault at Y=0. Multi-fault or non-planar faults not supported.
4. **Quasi-dynamic only**: No fully dynamic wave propagation. Radiation damping (eta*V) approximates dynamic effects. (Note: `feature/elasticity-inertia` branch in progress.)
5. **DG noise on unstructured meshes**: BR2 traction can show O(h) noise at sharp element-size transitions. Smooth mesh grading recommended.
