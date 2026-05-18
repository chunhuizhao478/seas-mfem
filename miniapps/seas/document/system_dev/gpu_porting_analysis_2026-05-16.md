# GPU Porting Analysis — SEAS-MFEM Paraview Branch

**Date:** 2026-05-16
**Author:** code-explore agent
**Scope:** Enable GPU (CUDA / HIP) support for the SEAS-MFEM-paraview branch,
covering dynamic rupture (TPV102, TPV104, TPV205) and quasi-dynamic (BP5)
benchmarks.
**Status:** Exploration + task plan. No code changes proposed yet.

---

## 1. Overview

The `seas-mfem-paraview` tree contains two physics regimes:

- **Dynamic rupture (TPV102 / TPV104 / TPV205)** — explicit, velocity–stress
  DG with Godunov + ADER-DG, custom hand-rolled assembly
  ([`dynamic/wave_operator.hpp`](miniapps/seas/dynamic/wave_operator.hpp),
  [`wave_operator.inl`](miniapps/seas/dynamic/wave_operator.inl), 5k+ lines).
- **Quasi-dynamic (BP5)** — implicit DG elasticity solved each ODE-RHS via
  MUMPS / CG-AMG, coupled to a rate-and-state friction ODE
  ([`domain/elasticity_operator.hpp`](miniapps/seas/domain/elasticity_operator.hpp),
  [`fault/rate_state_fault.hpp`](miniapps/seas/fault/rate_state_fault.hpp)).

Both regimes today run **CPU-only**. There is *no* `Device::Configure(...)`,
no `MFEM_FORALL`, no `Vector::UseDevice` anywhere in the SEAS miniapp; every
hot loop reads/writes via raw `Vector::GetData()` (host-only pointers). The
Phase 6 of [`dynamic_rupture_plan_v4.md`](miniapps/seas/document/system_dev/dynamic_rupture_plan_v4.md)
sketches a GPU plan but has not been started.

The MFEM library itself does ship a maturing GPU stack (CUDA, HIP, RAJA,
OCCA, CEED), built around `Device::Configure`, `MFEM_FORALL`,
`Vector::UseDevice`, and a per-integrator Partial Assembly (PA) /
Element Assembly (EA) path. The catch is that **MFEM's PA path currently
covers only tensor-product elements (hex), and the SEAS miniapp meshes are
all tetrahedral** (see `bp5/mesh/bp5.geo`, `tpv102/mesh/tpv102_*m.geo`,
`tpv205/mesh/tpv2053d_100m.geo`). On tets, GPU acceleration in MFEM has to
go through explicit `mfem::forall` kernels written by the application, or
through the `dFEM` (auto-diff) miniapp infrastructure new in MFEM 4.9.

---

## 2. Scope

| Target benchmark | Regime | Hot operator | Driver |
|---|---|---|---|
| TPV102 | Dynamic, RS friction | `WaveOperator::Mult` / `AdvanceADER` | `drivers/tpv102_driver.cpp` |
| TPV104 | Dynamic, FVW friction (stable asinh) | same + `Tpv104SubStepIterator` | `drivers/tpv104_driver.cpp` |
| TPV205 | Dynamic, LSW (linear slip-weakening) | same + `EvaluateADER_LSW` | `drivers/tpv205_driver.cpp` |
| BP5    | Quasi-dynamic, RS friction | `ElasticityDomainOperator::Solve` + `ComputeTraction` + `RateStateFaultOperator::ComputeRHS` | `drivers/seas_driver.cpp` |

**Approximate code size to touch:**
- `dynamic/` — ~16 200 LOC (12 large `.hpp`/`.inl`/`.cpp`).
- `domain/elasticity_operator*` — ~5 000 LOC.
- `integrator/dg_elasticity_*` — ~2 000 LOC of custom DG face matrix assembly.
- `fault/rate_state_fault.hpp` — ~1 000 LOC.

---

## 3. Architecture (what each code path actually does)

### 3.1 Dynamic rupture pipeline

`WaveOperator<MeshType>` ([`dynamic/wave_operator.hpp`](miniapps/seas/dynamic/wave_operator.hpp))
solves the 3D velocity–stress system

$$\partial_t \mathbf{Q} + A_x \partial_x \mathbf{Q} + A_y \partial_y \mathbf{Q} + A_z \partial_z \mathbf{Q} = 0$$

with $\mathbf{Q} \in \mathbb{R}^9$ = (σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz, v_x, v_y, v_z).
Storage is **component-major**: `Q[c * ndof_total_ + e * ndof_per_el_ + i]`.

Two stepping paths live in the same class:

1. **Classical explicit RK** — `Mult(Q, dQdt)` at [wave_operator.inl:684](miniapps/seas/dynamic/wave_operator.inl):
   ```
   ComputeVolumeRHS(Q, dQdt)                 // raw per-element loop
   ComputeFaceFluxRHS(Q, dQdt)               // raw per-face loop + Godunov + fault
   ComputeSharedFaceFluxRHS(Q, dQdt)         // MPI ghost faces (parallel only)
   ApplyPMLDamping(Q, dQdt)                  // optional sponge layer
   ApplyMassInverse(dQdt)                    // per-element DenseMatrix::Mult
   ```
2. **ADER-DG predictor–corrector** — `AdvanceADER(Q, dt, order, Q_new)` at
   [wave_operator.inl:4717](miniapps/seas/dynamic/wave_operator.inl):
   ```
   ComputeADERTimeIntegrated(Q, dt, order, I)  // Cauchy-Kovalevskaya recursion
   ComputeADERVolumeUpdate(I, rhs)             // strong-form volume on I
   ComputeADERFaceFluxRHS(I, dt, rhs)          // interior + boundary + fault
   ComputeADERSharedFaceFluxRHS(I, dt, rhs)    // shared faces
   ApplyMassInverse(rhs); Q_new = Q + rhs
   ```

Fault face flux is closed by `FaultFaceFlux::EvaluateADER` /
`EvaluateADER_LSW`, which runs a Brent or Newton root-find per fault QP
(rate-and-state) or a closed-form solve (LSW). The friction solver already
ships a dual-method API (`FrictionSolver::Method` =
{Brent, NewtonRaphson, NewtonRaphsonStable, HybridNRBisection}) where
**`Newton*` variants are explicitly labeled "for GPU"** in
[`dynamic/friction_solver.hpp:26`](miniapps/seas/dynamic/friction_solver.hpp) —
the API was designed with GPU in mind even though the dispatch site is
still CPU-only.

### 3.2 Quasi-dynamic (BP5) pipeline

The driver's `Mult(state, rate)` ODE step at
[`solver/seas_operator.hpp:318`](miniapps/seas/solver/seas_operator.hpp) calls,
per RK stage:

1. `domain_->Solve(t, slip, u)` — assembles RHS, solves the cached
   `HypreParMatrix` `K` via MUMPS, SuperLU, STRUMPACK, or HypreBoomerAMG+CG
   (selectable through TOML `[solver].solver_type`).
2. `domain_->ComputeTraction(u, slip, traction, normal_traction)` — DG flux
   recovery, written as **hand-rolled per-face per-QP loops on the host** in
   [`integrator/dg_elasticity_ip_combined_integrator.hpp`](miniapps/seas/integrator/dg_elasticity_ip_combined_integrator.hpp)
   (`ComputeTractionAtQuadPoints`) + `FaceQuadrature::GalerkinProject`.
3. `fault_->ComputeRHS(traction, state, rate)` — per-DOF friction root-find
   loop in [`fault/rate_state_fault.hpp:354`](miniapps/seas/fault/rate_state_fault.hpp).

Stiffness `K` is **cached** (assembled once) — so the per-step cost is
dominated by RHS rebuild, sparse solve, traction recovery, and the friction
loop.

---

## 4. MFEM GPU surface area

| MFEM capability | File | Status for SEAS |
|---|---|---|
| Backend abstraction `Device::Configure("cuda"|"hip"|"cpu")` | [`general/device.hpp`](general/device.hpp) | Not invoked anywhere in `miniapps/seas/`. Need to add to drivers. |
| `MFEM_FORALL` / `mfem::forall` kernel macros | [`general/forall.hpp`](general/forall.hpp) | Not used. Volume / face loops are plain `for(int e...)`. |
| `Vector::UseDevice(true)`, `Read/Write/ReadWrite`, `HostRead/HostWrite` | [`linalg/vector.hpp`](linalg/vector.hpp), [`general/mem_manager.hpp`](general/mem_manager.hpp) | Not used. All access via raw `GetData()` (host-only). |
| `BilinearForm::SetAssemblyLevel(AssemblyLevel::PARTIAL)` PA path | [`fem/bilinearform.hpp`](fem/bilinearform.hpp) | Could in principle be used for BP5, but our `DGElasticityIPCombinedIntegrator` / `DGElasticityBR2Integrator` only implement `AssembleFaceMatrix` (full-matrix path). No `AssemblePA` / `AddMultPA`. |
| Built-in `ElasticityIntegrator::AssemblePA` | [`fem/integ/bilininteg_elasticity_pa.cpp`](fem/integ/bilininteg_elasticity_pa.cpp) | Available — but currently **tensor-product (hex) only** and only assembles the volume term, not the SIPG / BR2 face terms we need. |
| `DGDiffusionIntegrator` PA | [`fem/integ/bilininteg_dgdiffusion_pa.cpp`](fem/integ/bilininteg_dgdiffusion_pa.cpp) | Scalar, hex-only, "no derivatives on faces" — doesn't match our 3D vector SIPG. |
| `DGTraceIntegrator` PA | [`fem/integ/bilininteg_dgtrace_pa.cpp`](fem/integ/bilininteg_dgtrace_pa.cpp) | Scalar trace — usable for the antiplane (BP2) scalar path, not for full 3D elasticity. |
| `L2FaceRestriction` (face DOF scatter/gather) | [`fem/restriction.hpp`](fem/restriction.hpp) | Needed for any PA-style face kernel. Not used by us. |
| Hypre on GPU (`HypreUsingGPU()`) — AMG / CG / GMRES move to device | [`linalg/hypre.cpp`](linalg/hypre.cpp) | Requires HYPRE ≥ 2.31.0 (runtime selectable CPU/GPU) per [`INSTALL:647-651`](INSTALL). MUMPS/SuperLU stay CPU. |

**Bottom line:** for the dynamic-rupture path we cannot lean on MFEM's PA
machinery at all (custom 9-component flux + Godunov + ADER are not MFEM
integrators). We must write our own device kernels via `mfem::forall`. For
BP5 we have a mixed picture: the **linear solve** can move to GPU via
HYPRE+CUDA (BoomerAMG-on-GPU + CG-on-GPU), but the **DG assembly,
traction recovery, and friction loop** are still custom CPU code and
need device equivalents.

---

## 5. GPU readiness audit (what blocks porting today)

| Pattern (host-only) | Where | Impact |
|---|---|---|
| `const real_t *p = Q.GetData();` then loop on host | [wave_operator.inl:715,774,863,978,1637,2086,2787,3305,…](miniapps/seas/dynamic/wave_operator.inl) (~30 sites) | Every hot loop reads host pointers. Replace with `Q.Read()` + `mfem::forall`. |
| `for (int e = 0; e < ne_; e++) { ... cached_per_element data ... }` | `ComputeVolumeRHS`, `ApplyMassInverse`, `AssembleElementMassInverse`, `ApplySpatialDerivative`, ADER volume | These are the embarrassingly-parallel kernels — ideal for `forall(e, ne, ...)`. |
| `std::vector<DenseMatrix> elem_mass_inv_` | [wave_operator.hpp:615](miniapps/seas/dynamic/wave_operator.hpp), [wave_operator.inl:4928](miniapps/seas/dynamic/wave_operator.inl) | Per-element dense matrices stored as a host `std::vector<DenseMatrix>`. To use on device, pack into a flat `Vector` or `DenseTensor` with `UseDevice(true)`. |
| Per-face loop with `for (int f = 0; f < mesh_.GetNumFaces(); f++) { mesh_.GetInteriorFaceTransformations(f); ... }` | `ComputeFaceFluxRHS` (~600 LOC), `ComputeSharedFaceFluxRHS` | Requires `L2FaceRestriction` + precomputed per-face geometric factors on device. Hardest single kernel. |
| MFEM `BilinearForm + ElasticityIntegrator + DGElasticityIPCombinedIntegrator + …->Assemble(0)` | [elasticity_operator_assembly.inl:23-86](miniapps/seas/domain/elasticity_operator_assembly.inl) | Full assembly into `SparseMatrix`/`HypreParMatrix`. Stays on CPU unless we (a) write PA kernels for our custom DG integrators, or (b) accept that assembly is CPU but the solve moves to GPU via HYPRE. |
| Sparse solver = MUMPS / SuperLU / STRUMPACK | [elasticity_operator_assembly.inl:292-330](miniapps/seas/domain/elasticity_operator_assembly.inl) | These have no GPU path in MFEM. To go GPU on BP5's linear solve, must switch to `HypreBoomerAMG + CGSolver` (already wired at line 374, 388, 394). |
| Brent / NR friction root-find per fault DOF/QP | [fault/rate_state_fault.hpp:375-440](miniapps/seas/fault/rate_state_fault.hpp), [dynamic/fault_face_flux.cpp `EvaluateADER`](miniapps/seas/dynamic/fault_face_flux.cpp) | Embarrassingly parallel per DOF/QP. Branch divergence is the main concern; the codebase already has Newton-Raphson and HybridNRBisection variants intended for GPU ([`friction_solver.hpp:79`](miniapps/seas/dynamic/friction_solver.hpp)). |
| `std::map<int, int> fault_face_dof_offset_`, `std::vector<bool> shared_fault_elem1_on_plus_`, … | [wave_operator.hpp:721-746](miniapps/seas/dynamic/wave_operator.hpp) | Host-side topology bookkeeping read inside the per-face kernels. Need to flatten to `Array<int>` (which can move to device) before kernel launch. |
| ParGridFunction::ExchangeFaceNbrData() per substep | [wave_operator.inl:147,1830,2787,4760](miniapps/seas/dynamic/wave_operator.inl) | MPI ghost exchange. With CUDA-aware MPI + `Memory::DEVICE` this can stay device-resident; without it, hits host↔device transfer in the hot loop. |
| Mesh is **tetrahedral** (Geometry::TETRAHEDRON) | [wave_operator.inl:71](miniapps/seas/dynamic/wave_operator.inl), all `*.geo` files | MFEM's PA path is currently hex-only. Blocks any drop-in PA use. |

---

## 6. Call graph (where GPU work would go)

```
seas_tpv102_driver / tpv104_driver / tpv205_driver  ─┐
                                                     ├─ Device::Configure("cuda"|"hip")   [add at startup]
                                                     │
   ODESolver::Step(WaveOperator, Q, t, dt)            │
     └─ WaveOperator::Mult(Q, dQdt)   OR             │
        WaveOperator::AdvanceADER(Q, dt, order, Qn)  │
           ├─ ComputeVolumeRHS         ◄── forall(e, ne, ...)
           ├─ ComputeFaceFluxRHS       ◄── forall(f, nf_interior, ...) +
           │                                forall(b, nf_boundary, ...) +
           │                                fault branch: forall(q, nq_fault, ...)
           ├─ ComputeSharedFaceFluxRHS ◄── ghost exchange + forall(qf, nq_shared, ...)
           ├─ ApplyPMLDamping          ◄── forall(e, ne, ...)
           ├─ ApplyMassInverse         ◄── forall(e, ne, ...)
           └─ (ADER) ApplySpatialDerivative + ComputeADERTimeIntegrated CK recursion
                                       ◄── forall(e, ne, ...) × O(order) iterations

seas_driver  (BP5 QD)                                 │
   PETSc TS / DormandPrince ─ seas_op.Mult(state, rate)
     ├─ domain_->Solve(t, slip, u)                   │
     │    ├─ AssembleSlipContribution  (CPU)         │
     │    ├─ AssembleDirichletLoading  (CPU)         │
     │    └─ solver_->Mult(B, X)       ◄── HYPRE+CUDA path (AMG-GPU + CG-GPU)
     ├─ domain_->ComputeTraction(u, slip, T)         │
     │    └─ DGElasticityIPCombined::ComputeTractionAtQuadPoints
     │                                 ◄── forall(f, nf_fault, ...) custom kernel
     └─ fault_->ComputeRHS(T, state, rate)
          └─ for i in [0, num_fault_dofs)
               FrictionSolver::SolveNR / SolveNRStable
                                       ◄── forall(i, num_fault_dofs, ...) (NR, fixed iters)
```

---

## 7. Detailed walkthrough of the hardest kernels

### 7.1 `WaveOperator::ComputeVolumeRHS` — easy

[wave_operator.inl:772-835](miniapps/seas/dynamic/wave_operator.inl)

- **Inputs:** `Q` (9 × `ndof_total_` doubles), per-element FE / Jacobian.
- **Math:** Strong-form DG volume integral
  $\text{rhs}_{c,i,e} \mathrel{+}= \sum_q w_q \sum_d (\partial_{x_d}\phi_i)(x_q) [A_d Q(x_q)]_c$.
- **Algorithm:** Triple nested loop (e → q → c). At each QP, contract shape
  function, evaluate `Q` at QP (sum over `ndof`), apply 3 dense 9×9 Jacobians
  `Ax_`, `Ay_`, `Az_`, scatter back to RHS.
- **Why it's easy to port:** purely element-local, fixed-shape arithmetic.
  `forall(e, ne, ...)` with shape/dshape evaluated through `DofToQuad`
  tables and Jacobians precomputed via `GeometricFactors`.

### 7.2 `WaveOperator::ApplyMassInverse` — easy

[wave_operator.inl:4897-4922](miniapps/seas/dynamic/wave_operator.inl)

- 9 small `DenseMatrix::Mult` calls per element. Replace `elem_mass_inv_`
  (a `std::vector<DenseMatrix>`) with a flat `Vector` packed as
  `[e * ndof² + i*ndof + j]` and `UseDevice(true)`, then a forall over
  `(e, ndof²)` does the matvec.

### 7.3 `WaveOperator::ComputeFaceFluxRHS` — hard

[wave_operator.inl:2084 onward](miniapps/seas/dynamic/wave_operator.inl) (~600 LOC)

- Per interior face: classify (regular / fault / boundary / shared seam),
  fetch `Elem1`, `Elem2` Jacobians + normals, evaluate `Q^-`, `Q^+`,
  call `GodunovFlux::Evaluate` or `FaultFaceFlux::Evaluate`, scatter
  contributions back to both element RHSs.
- **Branching:** 4 face categories (regular interior, central-flux mixed,
  fault interior, fault shared), and several sub-modes
  (`MixedFluxMode`, `UsePrecomputedFaceFluxes`, `FreeSurfaceBCMode`,
  `FaultFrictionLaw`). To make a single GPU kernel work, **separate the
  faces into typed arrays at setup time** and launch one `forall` per
  category — analogous to how MFEM's PA path uses different
  `L2FaceRestriction` instances per face type.
- The fault branch contains a root-finder (`FaultFaceFlux::EvaluateADER`)
  which has variable iteration count — use `SolveNRStable` (fixed iters)
  to keep warp divergence small.

### 7.4 `WaveOperator::AdvanceADER` — moderate, but recursive

[wave_operator.inl:4717-4827](miniapps/seas/dynamic/wave_operator.inl)

- Cauchy-Kovalevskaya recursion: `O-1` outer iterations, each calling
  `ApplySpatialDerivative` per direction `d ∈ {0,1,2}` (= 3 forall launches
  per CK step), plus a final volume + face flux update on the integrated
  state `I`.
- Storage cost: 5–8 element-sized scratch buffers (already cached as
  `ader_*_buf_` per R-1501). When ported, declare these with
  `UseDevice(true)` once at construction.

### 7.5 `RateStateFaultOperator::ComputeRHS` (BP5) — easy

[fault/rate_state_fault.hpp:354](miniapps/seas/fault/rate_state_fault.hpp)

- For each fault DOF: extract `tau` and `psi`, call
  `DieterichRuinaFriction::SolveSlipRateVectorPsi` (Brent), produce
  `[V_dip, V_strike, dpsi/dt]`. Independent across DOFs.
- The Newton variants already exist in `friction_solver.hpp` and the
  bookkeeping (per-DOF arrays of `a`, `eta`, `Dc`, `tau_pre`, `V_init`) is
  already in flat `Vector`s on `FaultGeometry`. Mostly a `UseDevice(true)`
  flip + a `forall(i, num_fault_dofs_, ...)` call.

### 7.6 BP5 stiffness / RHS / traction recovery — partly hard

- **`AssembleStiffness`** ([elasticity_operator_assembly.inl:23-86](miniapps/seas/domain/elasticity_operator_assembly.inl))
  uses full assembly into a `HypreParMatrix`. Assembly stays on CPU
  unless we write PA kernels for `DGElasticityIPCombinedIntegrator` /
  `DGElasticityBR2Integrator` (large effort, ~weeks; needs tensor-product
  faces in MFEM PA — not currently there for vector elasticity).
- **`AssembleSlipContribution`** ([assembly.inl:420](miniapps/seas/domain/elasticity_operator_assembly.inl))
  and **`AssembleDirichletLoading`** ([assembly.inl:1040](miniapps/seas/domain/elasticity_operator_assembly.inl))
  rebuild the RHS each step. These are face/boundary loops similar in
  shape to the wave operator's face flux; we'd need device kernels for
  full GPU coverage.
- **`ComputeTraction`** uses
  `DGElasticityIPCombinedIntegrator::ComputeTractionAtQuadPoints`
  ([integrator/dg_elasticity_ip_combined_integrator.hpp:434-554](miniapps/seas/integrator/dg_elasticity_ip_combined_integrator.hpp))
  followed by `FaceQuadrature::GalerkinProject` ([fault/face_quadrature.hpp:216](miniapps/seas/fault/face_quadrature.hpp)).
  Both are per-fault-face local — good `forall(f, nf_fault, ...)` targets.
- **The sparse solve** can move to GPU **without touching SEAS code** as long
  as the build links a `HYPRE ≥ 2.31` compiled with CUDA/HIP, runs
  `Device::Configure("cuda")` at startup, and the user selects
  `solver_type = "cg"` (which dispatches to `HypreBoomerAMG + CGSolver`,
  the GPU-capable path) instead of MUMPS / SuperLU / STRUMPACK
  (CPU-only). The `Hypre::Init`/`Hypre::UseGPU` call has to come from a
  driver-side change.

---

## 8. Conventions and constraints carried over from existing code

- **Component-major ordering** of `Q` (`Q[c * ndof_total_ + e * ndof_per_el_ + i]`)
  is asserted by the wave operator. Any kernel must preserve this — MFEM's
  `Ordering::byNODES` matches it; `byVDIM` does not.
- **Tetrahedral, GaussLobatto L2** discretisation throughout
  ([wave_operator.inl:36](miniapps/seas/dynamic/wave_operator.inl)). PA's hex restriction
  blocks naive use.
- **Static mesh assumption (R-1510, R-1510 documented at
  [wave_operator.hpp:435-444](miniapps/seas/dynamic/wave_operator.hpp)):** ghost
  exchange topology is built once. GPU port can mirror this — allocate
  device topology buffers once at construction.
- **R-1501** caches ADER scratch buffers as `mutable Vector`s to avoid
  per-step malloc. Same `Vector`s should be marked `UseDevice(true)` once
  and reused on device.
- **Shared-fault canonicalisation** (`shared_fault_elem1_on_plus_`,
  `interior_fault_elem1_on_plus_`): per-face flag arrays already exist
  precisely to remove per-QP branching that hurt MPI consistency. They
  also help GPU because the per-QP branch becomes a per-face flag read.
- **MPI consistency** rules in the wave operator
  ([wave_operator.hpp:385-398](miniapps/seas/dynamic/wave_operator.hpp), R-1600
  pairwise-collective contract): all ranks must call collective routines
  even on empty fault sets. GPU port must not bypass these regardless of
  device residency.
- **No silent fallbacks** policy from [CLAUDE.md](CLAUDE.md): if PA does
  not work for our integrators, we must error out, not silently revert
  to assembled. Same on GPU: a missing kernel must abort, not run a
  host fallback that masquerades as device execution.

---

## 9. Gotchas and non-obvious behaviour

1. **`MFEM_FORALL` body must compile with both host and device compilers.**
   That means no `std::vector`, no `std::sort`, no `dynamic_cast`, no
   exceptions, no virtual function calls inside the kernel. The existing
   `GodunovFlux::Evaluate` and friction solvers contain `if`/`for` only
   and inline math, so they are mostly portable — but each must be
   tagged `MFEM_HOST_DEVICE` and avoid `Vector`/`DenseMatrix` member
   methods that touch the memory manager.
2. **`Vector::GetData()` always returns the host pointer** and silently
   syncs the memory manager. Using it inside a `forall` kernel reads
   stale (possibly-not-on-device) data. Every hot-path `.GetData()` site
   in `wave_operator.inl` is therefore a porting hazard — there are about
   30 of them.
3. **Hypre's GPU path requires HYPRE ≥ 2.31.0**, and `Hypre::Init()` must
   happen *before* any `HypreParMatrix` is constructed. The current
   build script ([`build_frontera.sh`](build_frontera.sh)) is a CPU build;
   a GPU build needs `MFEM_USE_CUDA=YES` (or `MFEM_USE_HIP=YES`),
   `HYPRE_DIR` pointing at a GPU build, and `CUDA_ARCH` (e.g. `sm_80`).
4. **MUMPS, SuperLU_DIST, STRUMPACK have no GPU path.** Any benchmark that
   currently configures `solver_type = "mumps"` (default for BP5,
   see `bp5_500m.toml`) cannot be GPU-accelerated on the solve until the
   driver switches to `solver_type = "cg"` (which already builds a
   `HypreBoomerAMG` preconditioner — fine on GPU).
5. **CUDA-aware MPI** is needed to keep `ExchangeFaceNbrData` device-side.
   Without it, ghost exchanges trigger device→host→MPI→host→device round
   trips inside the time-step hot loop, often eating the entire GPU gain.
6. **Brent's method has variable iteration count → warp divergence on GPU.**
   The plan v4 (Phase 6.6) already calls this out; the dual-solver API
   in `friction_solver.hpp` (Brent for CPU, NewtonRaphson / NRStable for
   GPU) is the planned mitigation. NR-Stable already exists and is the
   default for TPV104.
7. **MFEM PA is hex-only for the vector elasticity integrator we'd want**
   (per [https://mfem.org/performance/](https://mfem.org/performance/)).
   We cannot just write `cached_a_->SetAssemblyLevel(AssemblyLevel::PARTIAL)`
   and expect BP5 to fly — we'd have to (a) remesh BP5 with hex,
   (b) write PA kernels for `DGElasticityIPCombinedIntegrator`, or
   (c) accept assembly stays on CPU. Option (c) is the realistic short
   path for BP5.
8. **`dFEM` (MFEM 4.9+, see CHANGELOG)** is the newer matrix-free / AD
   framework and *does* target tetrahedra + GPU. It may be worth a
   side-investigation as an alternative to hand-rolled `forall` kernels,
   but it's experimental (`mfem::future` namespace, API can change).

---

## 10. Recommended task plan

Phased to minimise breaking the existing working CPU code. Each phase
keeps the byte-for-byte CPU path intact and adds a device path that's
selected at runtime via `Device::Configure(...)`.

### Phase 0 — Build & smoke (1–2 weeks)

- Add `MFEM_USE_CUDA` (and/or `MFEM_USE_HIP`) flags to the SEAS Makefile.
  Provide a `build_gpu.sh` mirroring `build_frontera.sh`.
- Rebuild HYPRE with CUDA/HIP (>= 2.31.0).
- Add `--device <cuda|hip|cpu>` CLI flag to **all four drivers**
  (`seas_driver.cpp`, `tpv102_driver.cpp`, `tpv104_driver.cpp`,
  `tpv205_driver.cpp`) and call `Device::Configure(device_config)` +
  `Device::Print()` at startup (pattern from `examples/ex9p.cpp:330`).
- Add `Mpi::Init()` ordering check (must come before `Device::Configure`).
- Run smoke test: confirm CPU runs are byte-identical when `--device cpu`
  is selected explicitly.
- Deliverable: every benchmark builds with CUDA, runs with `--device cpu`
  with no regression.

### Phase 1 — BP5 sparse solve on GPU (1 week, biggest near-term win)

- In `seas_driver.cpp` and `elasticity_operator_assembly.inl`, gate
  `solver_type = "mumps"` to CPU and select `cg+amg` automatically when
  `Device::Allows(Backend::CUDA_MASK | Backend::HIP_MASK)`.
- Add `Hypre::Init()` + (where applicable) `Hypre::SetUseGPU(true)` before
  `HypreParMatrix` construction.
- Mark RHS / displacement / Lagrange vectors `UseDevice(true)`.
- Acceptance: BP5 500m TOML runs with `--device cuda`, traction time-series
  agrees with CPU reference to `1e-10` relative (HYPRE AMG is not bit-exact
  with MUMPS — use the existing `scripts/regression_check.py` 5% tolerance
  band as outer envelope, but also compare the *change* across CPU↔GPU on
  the same `cg+amg` solver to verify <1e-10).
- Deliverable: BP5 produces a measurable speedup on the linear solve
  fraction.

### Phase 2 — Per-DOF friction kernels on GPU (1 week)

- `RateStateFaultOperator::ComputeRHS` ([rate_state_fault.hpp:354](miniapps/seas/fault/rate_state_fault.hpp))
  is the single tightest per-step CPU loop after the linear solve. Lift
  `for (int i = 0; i < num_nodes_; i++)` into `mfem::forall(i, n, ...)`.
- Switch to `FrictionSolver::Method::NewtonRaphsonStable` automatically on
  device (fixed-iter Newton).
- Mark `tau_pre_`, `slip_rate_`, `a_values`, `eta_values`, `Dc_values` all
  `UseDevice(true)` once at construction.
- Same change applies to the per-fault-QP loop in
  `FaultFaceFlux::EvaluateADER` (dynamic rupture).
- Acceptance: friction RHS on GPU matches Brent CPU to `1e-8` per the
  existing Newton-vs-Brent test (`tests/unit/test_friction_law.cpp`).

### Phase 3 — Dynamic rupture `Mult` volume + mass-inverse on GPU (2 weeks)

- Port `ComputeVolumeRHS`, `ApplyMassInverse`, `AssembleElementMassInverse`
  to `forall(e, ne, ...)`. Pack `elem_mass_inv_` as a flat `DenseTensor`
  with `UseDevice(true)`.
- Pre-evaluate per-element `B`, `G`, `detJ` tables via
  `mesh.GetGeometricFactors(...)` so they live on device.
- Wire a `--volume-on-device` toggle so CPU/GPU paths can be A/B-tested.
- Acceptance: TPV102 / TPV104 / TPV205 volume RHS on GPU bit-matches CPU
  to `1e-12`; full `Mult` (with face flux still on CPU via a copy-back)
  matches CPU energy conservation in the existing reflecting-box test.

### Phase 4 — Dynamic rupture face flux on GPU (3–4 weeks, hardest)

- Pre-bucket faces by (`Interior`, `Boundary-Absorbing`, `Boundary-FreeSurface`,
  `Fault-Interior`, `Fault-Shared`, `Central-Flux-MF`) at construction.
  Each bucket gets its own `Array<int>` of mesh face indices + flat
  geometric-factor buffers.
- Write one `forall(f, n_bucket, ...)` kernel per bucket. Fault branches
  call the device-friendly NR friction solver.
- Replace the per-face `std::map<int,int>` lookups by flat `Array<int>`
  indexed lookups before kernel launch.
- For shared faces, require CUDA-aware MPI **or** explicitly
  `Memory::HostRead()` + `MPI_Isend` from host (slower; needs a build-time
  toggle).
- Acceptance: TPV104 production run (60s wallclock target) on GPU agrees
  with the existing chunhui-benchmark baseline at the same tolerances as
  the CPU run.

### Phase 5 — ADER on GPU (2 weeks, depends on 3 + 4)

- `ComputeADERTimeIntegrated` and `ApplySpatialDerivative` are element-local
  — same kernel pattern as Phase 3.
- `ComputeADERVolumeUpdate` reuses Phase 3 volume kernel on `I` instead of
  `Q`.
- `ComputeADERFaceFluxRHS` / `ComputeADERSharedFaceFluxRHS` reuse Phase 4
  face kernels with `I` and `dt` arguments and the ADER imposed-state
  side-channel ([wave_operator.hpp:350-352](miniapps/seas/dynamic/wave_operator.hpp)).
- Acceptance: substep-iterator TPV104 (round-11 R-602/R-603 path) matches
  CPU within the same tolerances; per-substep ghost-exchange MPI contract
  preserved.

### Phase 6 — BP5 DG assembly / traction recovery on GPU (open-ended)

- Two options:
  - **A. Stay assembled.** Accept that `AssembleStiffness` and RHS rebuild
    run on host; only the sparse solve, the per-DOF friction loop, and
    the traction recovery move to device. This is the **realistic
    near-term path** for BP5.
  - **B. Port custom DG integrators to PA.** Write `AssemblePA` /
    `AddMultPA` for `DGElasticityIPCombinedIntegrator` and
    `DGElasticityBR2Integrator`. Massive effort, blocked on
    tet-PA in MFEM for vector elasticity.
- For option A: port `DGElasticityIPCombinedIntegrator::ComputeTractionAtQuadPoints`
  and `FaceQuadrature::GalerkinProject` to `forall(f, nf_fault, ...)`.

### Phase 7 — Performance & profiling (continuous)

- Use Nsight Compute (CUDA) or Omniperf / ROCprof (HIP) to measure
  kernel occupancy, memory bandwidth, register pressure.
- Set baseline targets from plan v4 Phase 6: **≥ 5× speedup vs single
  core CPU** for TPV102. Compare against published SeisSol Fused-GEMMs
  performance (60% improvement over their previous GPU baseline).

---

## 11. Open questions

1. **Hex meshing.** Are we willing to re-mesh BP5 / TPV102 / TPV104 / TPV205
   with hexahedra to unlock MFEM's PA path? Likely not — the SCEC BP5
   geometry assumes simplex grading, and SeisSol / Tandem both run tets.
   If we stick with tets, all GPU kernels are custom `forall`.
2. **MFEM's `dFEM`** (4.9+) targets tets+GPU; should it be evaluated as an
   alternative to hand-rolled `forall` for the volume term? Worth a 1-week
   prototype before committing to Phase 3's manual kernels.
3. **Target hardware.** NVIDIA A100/H100 (CUDA) vs AMD MI250X (HIP)?
   Affects `CUDA_ARCH` / `HIP_ARCH` build flags and influences kernel
   tile-size tuning. Frontera (CPU-only) is irrelevant for this; we need
   a different test bed (Perlmutter, Frontier, Aurora).
4. **PETSc TS GPU integration.** The BP5 driver uses `PetscODESolver`
   (`drivers/seas_driver.cpp:229`). PETSc supports CUDA via Hypre/PETSc
   stack but our `seas_op.Mult()` must keep state on device. Confirm
   PETSc's `Vec` type matches MFEM's device memory.
5. **CUDA-aware MPI availability** on target clusters. If not available,
   per-step ghost-exchange transfer cost may dominate; we'd need to
   batch-amortise (R-1601 already does this for the substep ghost
   exchange — extend to RK stages).
6. **Branch naming.** The user said "create new branch with
   seas-mfem-paraview"; suggest something like
   `feature/gpu-port` or `feature/cuda-hip-baseline`. Phase 0–2 work
   could land on one branch; Phase 3–5 might want a sub-branch per kernel.

---

## 12. Sources

- [MFEM GPU support docs](https://mfem.org/gpu-support/) — Vector::Read/Write,
  `mfem::forall`, `Device::Configure` contract.
- [MFEM partial assembly](https://mfem.org/performance/) — PA hex-only
  limitation, DG decomposition with `L2FaceRestriction`.
- [MFEM GPU Tips & Tricks (LLNL)](https://software.llnl.gov/news/2021/02/17/mfemgpu/).
- [SeisSol GPU docs](https://seissol.readthedocs.io/en/latest/gpus.html) —
  CUDA 11.8 / ROCm 6.0 baseline, MPI+CUDA/HIP/SYCL.
- [Dorozhinskii et al. 2024, *Fused GEMMs for ADER-DG in SeisSol*](https://onlinelibrary.wiley.com/doi/full/10.1002/cpe.8037)
  — 60% GPU performance improvement, validated on Northridge 1994. Useful
  reference for ADER-DG GPU kernel design.
- [Tandem GPU enhancement (KONWIHR / ChEESE-2P)](https://www.konwihr.de/konwihr-projects/gpu-performance-and-feature-enhancement-of-the-earthquake-cycle-simulation-software-tandem/)
  — quasi-dynamic SEAS GPU porting in flight on SuperMUC-NG.
- [Modave et al. 2016, *GPU performance analysis of nodal DG for acoustic/elastic*](https://www.sciencedirect.com/science/article/abs/pii/S0098300416300668)
  — Bernstein–Bézier kernels, 28× single-precision speedup on multi-GPU.
- [MFEM CHANGELOG (4.9 release)](CHANGELOG) — `dFEM` introduction, particle
  methods, hyperbolic boundary integrators (relevant for future face PA).
- [Existing internal plan: dynamic_rupture_plan_v4.md §Phase 6](miniapps/seas/document/system_dev/dynamic_rupture_plan_v4.md)
  — pre-existing Phase 6 GPU plan, supersede with this analysis.
