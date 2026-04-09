# SEAS Miniapp Architecture

MFEM-based Discontinuous Galerkin code for SCEC SEAS (Sequences of Earthquakes and Aseismic Slip) benchmark problems. Couples a DG elasticity domain solver with rate-and-state friction on an embedded fault interface via a quasi-dynamic approximation. Supports 2D antiplane (BP1/BP2) and 3D full elasticity (BP5).

## Module Inventory

| Module | Directory | Purpose |
|--------|-----------|---------|
| **Config** | `config/` | Benchmark parameters: `bp1_params.hpp`, `bp2_params.hpp`, `bp5_params.hpp` |
| **Domain** | `domain/` | DG elasticity solvers: antiplane (scalar Laplace) and 3D vector elasticity |
| **Fault** | `fault/` | Fault geometry, basis transforms, rate-and-state operator, face quadrature |
| **Friction** | `friction/` | Friction laws (Dieterich-Ruina), state evolution (aging, slip, psi-space) |
| **Integrator** | `integrator/` | DG bilinear form integrators: BR2 and IP penalty for scalar and 3D elasticity |
| **Solver** | `solver/` | SEAS quasi-dynamic ODE operator, Dormand-Prince RK45 time stepper |
| **I/O** | `io/` | Benchmark output (SCEC format), ParaView VTK, checkpoint, probe sampling |
| **Common** | `common/` | MPI context, parallel utilities, type definitions |
| **Trace** | `trace/` | Face-level traction diagnostics (stress/correction decomposition) |
| **Driver** | `pseas.cpp` | Main parallel BP2 driver |
| **Tests** | `tests/` | Unit (24), verification (7), parallel (9) test programs |

## Class Hierarchy

```
TimeDependentOperator (MFEM)
  +-- SEASQuasiDynamicOperator<MeshType, DomainOpType, FaultOpType>
        |-- domain_: DomainOperator<MeshType>*
        |     +-- AntiplaneDomainOperator<MeshType>     (BP1/BP2: scalar)
        |     +-- ElasticityDomainOperator<MeshType>     (BP5: 3D vector)
        |-- fault_:  RateStateFaultOperator<MeshType, SlipComponents>*
        |     |-- geom_:      FaultGeometry<MeshType>
        |     |-- friction_:  FrictionLaw*
        |     |     +-- DieterichRuinaFriction
        |     |-- evolution_: StateEvolution*
        |     |     +-- AgingLaw / SlipLaw / AgingLawPsi
        |-- mpi_ctx_: MPIContext*
        +-- face_tracer_: FaceTraceLogger* (optional diagnostics)

DormandPrinceRK45 (custom adaptive RK45 time integrator)
  |-- mpi_ctx_: MPIContext* (parallel error reduction)

FaultBasis (fault coordinate transforms: global <-> dip/strike/normal)
FaceQuadrature (high-order fault DOF locations on faces)
```

## Execution Flow: One Time Step

```
DormandPrinceRK45::Step(seas_op, state, t, dt)
  |
  |-- For each RK stage k (7 stages, Dormand-Prince FSAL):
  |     |
  |     seas_op.Mult(state_k, rate_k)
  |       |
  |       |-- 1. fault_->GetSlip(state, slip_)
  |       |       Extract slip components from interleaved state vector
  |       |
  |       |-- 2. domain_->ExpandOwnedToLocalFault(slip_, local_slip_)
  |       |       MPI: broadcast owned DOFs to ghost DOFs on neighbor ranks
  |       |
  |       |-- 3. domain_->Solve(t, local_slip_, displacement_)
  |       |       Assemble RHS: slip contribution + Dirichlet loading
  |       |       Reuse stiffness matrix (assembled once)
  |       |       Solve K*u = b via MUMPS/CG+AMG/etc.
  |       |
  |       |-- 4. domain_->ComputeTraction(displacement_, local_slip_, local_traction_)
  |       |       Loop over fault interior + shared faces
  |       |       DG flux: tau = mu*{{grad(u).n}} - penalty*(jump(u) - slip)
  |       |       BP5: transform to fault-local (dip, strike, normal)
  |       |
  |       |-- 5. domain_->RestrictToOwnedFault(local_traction_, traction_)
  |       |       Select owned subset from full local traction
  |       |
  |       +-- 6. fault_->ComputeRHS(traction_, state, rate)
  |               For each owned fault DOF:
  |                 tau_total = tau0 + traction
  |                 V = friction_->SolveSlipRate(tau_total, psi, sigma_n, eta, a)
  |                 dpsi/dt = evolution_->Rate(V, psi, Dc)
  |               MPI_Allreduce(V_max, MPI_MAX) for global max
  |
  |-- Compute error estimate from embedded 4th-order solution
  |-- MPI_Allreduce(err_norm, MPI_MAX) for consistent accept/reject
  |-- Accept: advance t, state = y_next, k[0] = k[6] (FSAL)
  +-- Reject: shrink dt, retry
```

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

## Time Integration

Dormand-Prince RK45 (7-stage, 5th-order, FSAL) with adaptive step control:

- Error norm: `max_i |err_i / (atol + rtol * |y_i|)|` (L-infinity, default)
- PI controller: `dt_new = safety * dt * err_norm^(-1/5)`
- Parallel: `MPI_Allreduce(err_norm, MPI_MAX)` ensures consistent accept/reject
- Settings: `atol = 1e-7`, `rtol = 1e-50` (pure absolute), `dt_min = 1e-6 s`, `dt_max = 0.5 yr`

## Boundary Conditions (BP5)

Controlled by `BCMode` enum:

| Mode | Dirichlet Faces | Natural Faces | Usage |
|------|----------------|---------------|-------|
| `FarField` (default) | attrs 1-4 (X=+-Lx, Y=+-Ly) | attrs 5-6 (Z=0, Z=-Lz) | Matches Tandem |
| `XOnly` | attrs 1-2 (X=+-Lx only) | 3-6 | Antiplane-style |
| `AllDirichlet` | attrs 1-6 (all) | none | Legacy (incorrect) |

Dirichlet loading: `u_X = sgn(Y) * Vp * t / 2` (along-strike plate motion).

Tandem mesh tag convention: Physical Surface 1 = Natural (top/bottom), 3 = Fault, 5 = Dirichlet (far-field).

## MPI Communication Pattern

1. **ExpandOwnedToLocalFault**: Before domain solve, broadcast owned slip to ghost DOFs on neighbor ranks (all-to-all scatter at partition boundaries)
2. **Domain Solve**: Each rank solves local portion; MFEM parallel solvers handle internal communication
3. **ComputeTraction**: Each rank computes traction at its local fault faces (interior + shared)
4. **RestrictToOwnedFault**: Select owned DOF subset from full local traction
5. **ComputeRHS**: Each rank computes friction/state evolution for owned DOFs only
6. **GetMaxSlipRate**: `MPI_Allreduce(MPI_MAX)` for global V_max
7. **RK45 error norm**: `MPI_Allreduce(MPI_MAX)` for consistent accept/reject
8. **Output gather**: `MPI_Gatherv` to root + deduplication at partition boundaries

## Coordinate System (BP5)

```
SCEC Convention:        Tandem/MFEM Mesh:
  x1 = fault-normal  ->  Y (fault at Y=0)
  x2 = along-strike  ->  X
  x3 = depth (+down)  ->  -Z (Z=0 surface, Z<0 depth)
```

Fault-local frame (FaultBasis):
- `tangent1` = dip direction = (0, 0, +1) (downward into earth)
- `tangent2` = strike direction = (1, 0, 0) (along-strike, dominant slip)
- `normal` = (0, 1, 0) or (0, -1, 0) depending on face orientation

## Initialization Sequence (4-Phase)

1. **PreInit**: slip=0, psi = f0 + b*ln(V0/V_init) (steady-state placeholder)
2. **First domain solve**: K*u=b with zero slip, compute initial elastic traction
3. **Init**: Compute psi(0) from stress equilibrium: `tau0 + tau_elastic = sigma_n * f(V_init, psi) + eta * V_init`
4. **Verification re-solve**: Re-solve domain, re-compute traction, verify equilibrium error < 1e-6

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

## File Map

```
miniapps/seas/
  pseas.cpp                              Main parallel BP2 driver
  Makefile                               Build: make all / make test
  CMakeLists.txt                         CMake alternative
  ARCHITECTURE.md                        This file
  CLAUDE.md                              AI assistant instructions
  CODEBASE_GUIDE.md                      Physics-to-numerics walkthrough

  config/
    bp1_params.hpp                       BP1 benchmark parameters
    bp2_params.hpp                       BP2 benchmark parameters + mesh generator
    bp5_params.hpp                       BP5 parameters, spatial functions a(x2,x3), Dc(x2,x3)

  domain/
    domain_operator.hpp                  Abstract domain solver interface
    antiplane_operator.hpp               2D scalar DG (BP1/BP2), BR2 + IP
    antiplane_bdrload_operator.hpp       Variant with boundary-only loading
    elasticity_operator.hpp              3D vector DG elasticity (BP5), BR2 + IP
    bp2_mesh.hpp                         BP2 graded mesh generator (BP2MeshGenerator)
    seas_boundary_tags.hpp               Gmsh physical surface tag definitions

  fault/
    fault_geometry.hpp                   Fault DOF management, depth-dependent parameters
    fault_basis.hpp                      FaultBasis: global <-> (dip, strike, normal) transforms
    fault_nodes.hpp                      Fault node identification
    rate_state_fault.hpp                 RateStateFaultOperator: slip rate + state evolution
    face_quadrature.hpp                  FaceQuadrature: high-order DOF locations on faces

  friction/
    friction_law.hpp                     Abstract friction law interface
    dieterich_ruina.hpp                  Regularized Dieterich-Ruina, Brent solver, psi-space
    state_evolution.hpp                  AgingLaw, SlipLaw, AgingLawPsi

  integrator/
    dg_br2_integrator.hpp               Scalar BR2 DG integrator (BP1/BP2)
    dg_elasticity_br2_integrator.hpp     3D elasticity BR2 integrator (BP5)
    dg_elasticity_ip_penalty_integrator.hpp      IP penalty with full elasticity tensor
    dg_elasticity_ip_combined_integrator.hpp      Combined IP integrator

  solver/
    seas_operator.hpp                    SEASQuasiDynamicOperator: couples domain + fault
    seas_bdrload_operator.hpp            Variant with boundary-load driving
    time_stepper.hpp                     DormandPrinceRK45: adaptive RK45, parallel error

  io/
    benchmark_output.hpp                 SCEC format serial output
    parallel_benchmark_output.hpp        Parallel output with gather + dedup
    bp5_benchmark_output.hpp             BP5-specific SCEC output
    bp5_parallel_output.hpp              BP5 parallel output
    paraview_output.hpp                  ParaView VTK export
    probe_output.hpp                     Point probe sampling
    checkpoint.hpp                       Checkpoint/restart

  common/
    seas_types.hpp                       Enums, typedefs
    mpi_context.hpp                      MPIContext: Allreduce, Gatherv wrappers
    parallel_utils.hpp                   Parallel utility functions

  trace/
    face_trace_logger.hpp                Per-face traction decomposition diagnostics

  tests/
    unit/                                24 unit tests (friction, DG, fault, I/O, etc.)
    verification/                         7 verification tests (BP1/BP2/BP5 full runs)
    parallel/                             9 parallel tests (consistency, MPI, scaling)

  bp1/mesh/                              BP1 Gmsh geometries and meshes
  bp2/mesh/                              BP2 Gmsh geometries and meshes
  bp5/mesh/                              BP5 Gmsh geometries and meshes (250m-4000m)
  bp5/frontera/                          Frontera HPC SLURM scripts
  jobs/bp5/                              130+ SLURM job scripts (various configs)

  debug_document/bp5_debug_document/     Debug history: v1-v62 (critical domain knowledge)
```

## Known Limitations

1. **Matrix assembly**: Stiffness matrix assembled and stored (not matrix-free). Memory-bound for large meshes at high polynomial order.
2. **Direct solver default**: MUMPS (or MUMPS_BLR) for BP5. Iterative solvers (CG+AMG) available but less robust for DG.
3. **Single fault plane**: Code assumes one planar fault at Y=0. Multi-fault or non-planar faults not supported.
4. **Quasi-dynamic only**: No fully dynamic wave propagation. Radiation damping (eta*V) approximates dynamic effects.
5. **DG noise on unstructured meshes**: BR2 traction can show O(h) noise at sharp element-size transitions. Smooth mesh grading recommended.
