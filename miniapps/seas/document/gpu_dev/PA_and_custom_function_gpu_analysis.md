# Partial Assembly (PA) & Custom-Function GPU-Consistency Analysis

**Scope of this document.** It answers two specific questions for the
SEAS-MFEM-SAFS codebase:

1. **To incorporate MFEM partial assembly (PA), which parts of the code
   must change?**
2. **Which custom functions we wrote must change to be consistent with
   MFEM's GPU support?**

It covers **both** solve paths:

- **Dynamic rupture** — `drivers/spatial_dyn_driver.cpp` →
  `dynamic/wave_operator.{hpp,inl}` (and TPV102/104/205).
- **Quasi-dynamic rupture (BP5)** — `drivers/seas_driver.cpp` →
  `solver/seas_operator.hpp` → `domain/elasticity_operator*.inl`.

**Relationship to the existing plans.** The `gpu_dev/` directory already
contains `gpu_implementation_plan_v1/v2/v3.md`. Those plans deliberately
choose a **hand-written `mfem::forall` matrix-free** strategy and do
**not** use MFEM's `AssemblyLevel::PARTIAL` machinery. This document
explains **why** literal MFEM PA is mostly *inapplicable* here, identifies
the few places where it *is* applicable, and inventories the custom
functions that block GPU execution regardless of which strategy is
chosen. Read this as the "PA-feasibility + custom-kernel readiness"
companion to the v3 phase plan, not a replacement for it.

> **Terminology guard.** Throughout, "PA" = MFEM's
> `BilinearForm::SetAssemblyLevel(AssemblyLevel::PARTIAL)` — operator
> action via cached quadrature-point data and `AddMultPA`, *no global
> sparse matrix*. "Matrix-free `forall`" = our own
> `mfem::forall(...)` device kernels (what v3 calls Phases 3–6). They
> are different things; the codebase today uses **neither**.

---

## TL;DR (the four findings that drive everything)

1. **The code uses zero PA and zero device today.** A tree-wide grep for
   `SetAssemblyLevel`, `AssemblyLevel::PARTIAL`, `AssemblePA`,
   `AddMultPA`, `mfem::forall`, `MFEM_FORALL`, `MFEM_HOST_DEVICE`,
   `UseDevice`, `Device::` returns **no hits** anywhere under
   `miniapps/seas`. Everything is full assembly + host loops.

2. **Both meshes are tetrahedral (simplex) DG/L2 spaces.** BP5 uses
   `DG_FECollection(order, 3, GaussLobatto)` with `vdim=3, byNODES`
   (`domain/elasticity_operator_setup.inl:11,13`); the wave operator uses
   `L2_FECollection(order, dim, GaussLobatto)` on a *scalar* space with a
   hand-rolled 9-component layout (`dynamic/wave_operator.inl:36,41`).
   MFEM's *fast* PA kernels are tensor-product (hex/quad) sum-factorized
   kernels; on simplices PA falls back to dense per-element
   (`DofToQuad::FULL`) evaluation, which is far less of a win.

3. **BP5's DG face terms have no PA path — not in our integrators and
   not in MFEM's.** Our `DGElasticityBR2Integrator`,
   `DGElasticityIPCombinedIntegrator`, `BR2InteriorFaceIntegrator`,
   `DGElasticityIPPenaltyIntegrator` implement only `AssembleFaceMatrix`
   (full assembly). MFEM's own `DGElasticityIntegrator` *also* implements
   only `AssembleFaceMatrix` — it has **no** `AssemblePA`. Only the
   *volume* `ElasticityIntegrator` is PA-capable. So the BP5 stiffness
   operator cannot be put into PA mode as-is.

4. **BP5's solver is incompatible with PA by construction.** Default is
   `SolverType::MUMPS_BLR` (`domain/elasticity_operator.hpp:111`);
   production TOMLs set `solver_type = "mumps"`
   (`config/bp5_*.toml`). MUMPS / SuperLU / STRUMPACK are *sparse-direct*
   factorizations, and BoomerAMG needs the assembled CSR. **PA produces
   no matrix.** Switching BP5 to PA forces a switch to a *matrix-free
   iterative* solver + a *matrix-free preconditioner* — a much larger
   change than "flip a flag."

**Net consequence.** "Incorporate PA" means different things on the two
paths:

| Path | What "PA" literally means here | Verdict |
|------|--------------------------------|---------|
| **Dynamic (WaveOperator)** | It is *already* matrix-free hand-rolled — there is no `BilinearForm` to set `AssemblyLevel::PARTIAL` on. "PA" reduces to *device-porting the existing element/face/QP loops* (PA-*style*), plus optionally using MFEM `MassIntegrator` PA for the mass solve. | **Feasible & is the v3 plan.** No literal PA. |
| **BP5 (elasticity)** | Switch the cached `ParBilinearForm` to `AssemblyLevel::PARTIAL`. | **Blocked** by (a) no DG-face PA integrator and (b) sparse-direct solver. Only the *volume* term is PA-able; the *coupling* (slip RHS, traction recovery, Dirichlet) is hand-rolled and must be device-ported anyway. |

---

## 0. Background you need before reading the rest

### 0.1 What MFEM PA actually requires

To put an operator in PA mode you must satisfy *all* of:

- The operator is a `BilinearForm`/`ParBilinearForm` and you call
  `a.SetAssemblyLevel(AssemblyLevel::PARTIAL)` *before* `a.Assemble()`.
- **Every** integrator added to it implements `AssemblePA` +
  `AddMultPA` (and `AssembleDiagonalPA` if you want a diagonal
  preconditioner). If *one* integrator lacks PA, `Assemble()` aborts.
- You consume the operator through `a.Mult(x,y)` (operator action), and
  you solve with a **matrix-free iterative** method (CG/GMRES). There is
  **no `SpMat()` / `HypreParMatrix`** to hand to MUMPS or BoomerAMG.
- Preconditioning is via `OperatorJacobiSmoother` (from
  `AssembleDiagonalPA`), Chebyshev, or a **Low-Order-Refined (LOR)**
  AMG built from a sparse low-order proxy.
- Memory lives in MFEM `Memory<T>` and is accessed with
  `Read()/Write()/ReadWrite()`; kernels are `mfem::forall` with
  `MFEM_HOST_DEVICE` bodies. `Device::Configure("cuda"|"hip")` is called
  once at startup.

### 0.2 MFEM's GPU device-memory contract (applies to *all* custom code)

A function/kernel is "GPU-consistent" in MFEM only if:

- Data is an MFEM `Vector`/`Memory`/`DenseTensor` and the kernel takes
  the result of `.Read()/.Write()/.ReadWrite()` (which returns a pointer
  valid in the *current* memory space), **not** `.GetData()` /
  `.HostRead()` (host-only).
- The loop is `mfem::forall(N, [=] MFEM_HOST_DEVICE (int i){...})`, with a
  body that calls only `MFEM_HOST_DEVICE`-qualified functions.
- The body uses **no** heap allocation, **no** `std::vector` /
  `std::function` / `std::map`, **no** virtual dispatch, **no**
  `mfem::DenseMatrix`/`mfem::Vector` *temporaries* (those allocate);
  small fixed-size work uses C arrays / `real_t a[9]`.
- Writes that can collide across threads use `mfem::AtomicAdd<real_t>`.

Measured against that contract, essentially every hot-path function in
both solve paths is currently host-only. The inventory in §3 lists each
one and what specifically violates the contract.

---

## 1. Part A — Quasi-dynamic (BP5): where PA could go and why it mostly cannot

### 1.1 The assembly objects (PA-relevant inventory)

`ElasticityOperator::AssembleStiffness()`
(`domain/elasticity_operator_assembly.inl:9–441`) builds **one** cached
form and converts it to a matrix:

```cpp
cached_a_ = std::make_unique<BilinFormType>(fes_.get());            // :23  (Par)BilinearForm

// volume term — PA-CAPABLE in MFEM:
cached_a_->AddDomainIntegrator(new ElasticityIntegrator(λ, μ));     // :37-38

// DG face terms — NO PA anywhere:
cached_a_->AddInteriorFaceIntegrator(new DGElasticityBR2Integrator(...));      // :47-49 (BR2)
cached_a_->AddBdrFaceIntegrator(new DGElasticityBR2BoundaryIntegrator(...));   // :54-57
//   or, in IP mode:
cached_a_->AddInteriorFaceIntegrator(new DGElasticityIPCombinedIntegrator(...));// :67-69
cached_a_->AddBdrFaceIntegrator(new DGElasticityIPCombinedIntegrator(...));     // :73-76

cached_a_->Assemble(0);  cached_a_->Finalize();                    // :85-86  → SparseMatrix
cached_a_->ParallelAssemble(cached_Ah_);                           // :216    → HypreParMatrix
```

The matrix is assembled **once** and cached
(`stiffness_assembled_`, `cached_a_`, `cached_Ah_` —
`domain/elasticity_operator.hpp:468–470`), then reused across every RK
stage and time step. `SetAssemblyLevel` is **never** called.

### 1.2 The integrator PA-support matrix (the hard blocker)

| Integrator | Where | Base | Methods implemented | PA? |
|---|---|---|---|---|
| `ElasticityIntegrator` (volume) | MFEM `fem/integ/bilininteg_elasticity_pa.cpp:40` | `BilinearFormIntegrator` | `AssemblePA`, `AddMultPA`, `AddMultTransposePA`, `AssembleDiagonalPA` | **YES** |
| `DGElasticityIntegrator` (MFEM's own DG face) | MFEM `fem/bilininteg.hpp:3691` | `BilinearFormIntegrator` | `AssembleFaceMatrix` only | **NO** |
| `DGElasticityBR2Integrator` | `integrator/dg_elasticity_br2_integrator.hpp:47` | `BilinearFormIntegrator` | `AssembleFaceMatrix` only | **NO** |
| `DGElasticityIPCombinedIntegrator` | `integrator/dg_elasticity_ip_combined_integrator.hpp:30` | `BilinearFormIntegrator` | `AssembleFaceMatrix` + custom RHS/traction helpers | **NO** |
| `BR2InteriorFaceIntegrator` / `BR2BoundaryFaceIntegrator` | `integrator/dg_br2_integrator.hpp:41,81` | `BilinearFormIntegrator` | `AssembleFaceMatrix` only | **NO** |
| `DGElasticityIPPenaltyIntegrator` | `integrator/dg_elasticity_ip_penalty_integrator.hpp:50` | `BilinearFormIntegrator` | `AssembleFaceMatrix` only | **NO** |

> **Contrast:** MFEM's *scalar* `DGDiffusionIntegrator` *does* have a DG
> face PA path (`AssemblePAInteriorFaces` / `AssemblePABoundaryFaces`,
> MFEM `fem/bilininteg.hpp:3531–3538`). So DG-face PA is *possible in
> principle*, but **no elasticity DG-face PA exists** — neither ours nor
> MFEM's. Writing one is a research-grade task (the BR2 lifting operator
> and the IP elasticity-tensor penalty are both substantial).

**Conclusion for the stiffness operator:** because *one* of the
integrators (the DG face term) lacks `AssemblePA`, `Assemble()` in PA
mode would abort. PA cannot be enabled on the full BP5 stiffness without
first authoring an elasticity-DG-face PA integrator.

### 1.3 The solver blocker (independent of §1.2)

Even if every integrator had PA, the BP5 solve consumes an **assembled
matrix**:

- Default `SolverType::MUMPS_BLR` (`domain/elasticity_operator.hpp:111`),
  production TOMLs `solver_type = "mumps"`.
- `AssembleStiffness` wires MUMPS / SuperLU / STRUMPACK / GMRES+BlockILU /
  (CG|GMRES)+BoomerAMG onto the CSR/HypreParMatrix
  (`domain/elasticity_operator_assembly.inl:315–429`).

All of those need `SpMat()` or `HypreParMatrix`. **PA gives neither.** To
go PA you must:

1. Switch to a matrix-free Krylov solver (`CGSolver`/`GMRESSolver` whose
   `SetOperator` is the PA `BilinearForm` itself).
2. Replace the preconditioner. The DG elasticity stiffness is
   ill-conditioned (penalty scales like `p²/h`), so Jacobi/Chebyshev
   alone will be slow; the realistic option is **LOR-AMG**
   (`ParLORDiscretization` + `HypreBoomerAMG` on the low-order proxy) —
   but LOR for *DG elasticity* is itself non-trivial and unverified here.

This is why v3 §"Out of scope" keeps `AssembleStiffness` on the CPU and
v3 Map-1 Box-3 simply links a HIP-built HYPRE and flips `mumps`→`cg` so
**BoomerAMG runs on-device on the still-assembled matrix** — i.e. v3
*keeps full assembly* and offloads the *solver*, sidestepping PA entirely.
That is the pragmatic BP5 GPU story; PA is not it.

### 1.4 What in BP5 *is* PA-able / device-able

- **Volume term only.** `ElasticityIntegrator::AssemblePA` requires
  `Ordering::byNODES` and `vdim == dim` — BP5 satisfies both
  (`elasticity_operator_setup.inl:13`,
  `bilininteg_elasticity_pa.cpp:42,47`). On tets it uses `DofToQuad::FULL`
  (dense, non-tensor) — correct but not sum-factorized-fast. Useful only
  if you also build a matrix-free solver (see §1.3); on its own it buys
  nothing because the matrix is assembled once and cached.
- **The per-step coupling work is hand-rolled, not a `BilinearForm`**,
  so "PA" doesn't apply to it — it must be **device-ported as `forall`
  kernels** (this is v3 Phase 6, *not* PA):
  - Slip-contribution RHS — `AssembleSlipContributionIP/BR2(+Shared)`
    (`elasticity_operator_assembly.inl:447–1019`): per-fault-face loops,
    fresh `DGElasticityIPCombinedIntegrator` allocated *per face*
    (`:454`), BR2 path does raw tensor contractions over QPs with
    `elem_mass_inv_[e]` (`std::vector<DenseMatrix>`).
  - Dirichlet loading — `AssembleDirichletLoading`
    (`:1067–1596`): same structure over boundary + interior-skeleton
    faces.
  - Traction recovery — `ComputeTractionImpl`
    (`elasticity_operator_traction.inl:653–1600`): per-face loop calling
    `DGElasticityIPCombinedIntegrator::ComputeTractionAtQuadPoints(Decomposed)`
    then `ProjectTractionToFaultDOFs` (which does a per-face
    `DenseMatrixInverse`, `dg_elasticity_ip_combined_integrator.hpp:740`).
  - Friction RHS — the per-DOF rate-and-state solve (see §3.3), called
    from `SEASQuasiDynamicOperator::Mult` (`solver/seas_operator.hpp`).

### 1.5 BP5 recommendation (ordered)

1. **Do not pursue literal PA for BP5 first.** Adopt the v3 path: keep
   full assembly, link GPU-HYPRE, flip `mumps`→`cg`, and let BoomerAMG
   run on-device. This is the 40%-of-step win with the least risk.
2. **Device-port the coupling kernels** (slip RHS, Dirichlet, traction
   recovery, friction RHS) with `mfem::forall` — these are the §3
   custom functions. This is the bulk of the per-step CPU time that
   *isn't* the solve.
3. **Only if** the assembled-matrix solve becomes the bottleneck *and*
   memory pressure from the CSR matters, invest in (a) an
   elasticity-DG-face `AssemblePA` integrator and (b) an LOR-AMG
   matrix-free preconditioner. Treat as research, not porting.

---

## 2. Part B — Dynamic rupture (WaveOperator): "PA" = device-porting hand-rolled kernels

### 2.1 There is no `BilinearForm` to make partial

`WaveOperator::Mult(Q, dQdt)` (`dynamic/wave_operator.inl:854–936`) and
`AdvanceADER` are **entirely hand-rolled element/face/QP loops**. There
is no `BilinearForm`/`ParBilinearForm` anywhere in the wave operator, so
`SetAssemblyLevel(PARTIAL)` has nothing to attach to. The operator is
*already* matrix-free; "incorporate PA" here means **convert the existing
host loops into MFEM device kernels** (the PA-*style* matrix-free pattern
v3 calls Phases 3–5), and optionally route the mass solve through MFEM's
PA `MassIntegrator`.

The RK4 step decomposes as `dQ/dt = M⁻¹ · (ComputeVolumeRHS −
ComputeFaceFluxRHS)`:

| Sub-function | File:line | Loop | PA-style port |
|---|---|---|---|
| `ComputeVolumeRHS` | `wave_operator.inl:942–1005` | `for(e) for(q)`; reads `Q.GetData()`, stack `Vector shape`/`DenseMatrix dshape`, member `Ax_,Ay_,Az_` | `forall(ne,...)` one thread/element; this is the canonical PA volume kernel (∇φ·A·Q) |
| `ComputeFaceFluxRHS` | `wave_operator.inl:2254–3370` | one `for(f)` over all faces, branch on regular/fault/boundary/mixed | bucket faces by type once, then 4 branch-free `forall` kernels; scatter via `AtomicAdd` |
| `ComputeSharedFaceFluxRHS` | `wave_operator.inl:2923–3370` | `for(shared f)` + `ExchangeFaceNbrData()` | `forall` + GPU-aware MPI |
| `ApplyMassInverse` | `wave_operator.inl:5096–5121` | `for(e) for(c=0..8)` small `Minv.Mult` | `forall(ne,...)` register matvecs, **or** MFEM PA `MassIntegrator` (see §2.2) |
| ADER CK recursion | `wave_operator.inl` (`ApplySpatialDerivative`, `ComputeADERTimeIntegrated`) | `for(order) for(dir) for(e)` | reuse the volume `forall`; ping-pong device scratch |

### 2.2 The mass matrix — the *one* place literal MFEM PA fits cleanly

The element mass inverse is a host `std::vector<DenseMatrix>
elem_mass_inv_` (`wave_operator.hpp:802`), built once via
`DenseMatrixInverse` (`wave_operator.inl:5127–5157`) and applied as 9
small `DenseMatrix::Mult` per element (`:5096–5121`).

Two device options:

- **(a) Keep custom, device-port:** pack the inverses into one flat
  device `Vector` and apply with `forall(ne, ...)` (v3 Phase 3). Lowest
  risk; preserves exact numerics.
- **(b) Use MFEM PA `MassIntegrator`:** for an explicit DG scheme the
  mass matrix is block-diagonal; MFEM's `MassIntegrator` has full PA
  (`AssemblePA/AddMultPA/AssembleDiagonalPA`, MFEM
  `fem/bilininteg.hpp:2416–2433`) and `VectorMassIntegrator` likewise
  (`:2643–2646`). But MFEM PA gives `M·x` (apply), **not** `M⁻¹·x`. You'd
  need a CG-with-Jacobi per step or an explicit inverse — for DG you
  almost always want the explicit element inverse, so **(a) is
  preferred**. Note option (b) also requires the state to live in a
  `vdim=9` (or 9 scalar) FE space; today the wave op uses a *scalar* L2
  space with a manual component-major layout (`Q[c*ndof_total + dof]`),
  so (b) would also force a layout change.

### 2.3 Data-layout fact that dominates GPU performance

The state vector is **component-major**: `Q[c*ndof_total_ + dof_offset +
i]` (see `ApplyMassInverse`, `wave_operator.inl:5096–5121`, and the
`QIndex` enum in `dynamic/wave_state.hpp:27`). For coalesced GPU loads in
a one-thread-per-element kernel this is acceptable (each thread strides
by `ndof_total_` across the 9 components), but it is the layout every
ported kernel must respect; do **not** silently re-order to
interleaved `[dof*9 + c]` without auditing every kernel and the MPI
ghost-exchange packing (`wave_operator.inl:2025,2977,4534`).

---

## 3. Question [2] — Custom functions that must change for GPU consistency

This is the union over **both** paths. Each row states the file, what
breaks the device contract (§0.2), and the required change. None of these
are PA per se — they are the prerequisites for *any* device execution,
PA or matrix-free.

### 3.1 Dynamic-rupture custom functions

| Function / class | File:line | Why it's host-only today | Required change |
|---|---|---|---|
| `WaveOperator::ComputeVolumeRHS` | `wave_operator.inl:942` | `Q.GetData()` raw host ptr; stack `Vector`/`DenseMatrix` temporaries per QP | `forall(ne)`; `Q.Read()`; fixed-size `real_t` scratch; member `Ax_/Ay_/Az_` copied to device POD arrays |
| `WaveOperator::ComputeFaceFluxRHS` | `wave_operator.inl:2254` | single branching face loop; `GetFaceElementTransformations` per face; scatter races | precompute per-type face buckets + per-face geometry/normals once; 4 branch-free `forall`; `AtomicAdd` scatter |
| `WaveOperator::ComputeSharedFaceFluxRHS` | `wave_operator.inl:2923` | `ExchangeFaceNbrData()` host bounce; `std::set shared_mesh_face_set_` lookups; `std::memcpy` packing | GPU-aware MPI (`MPICH_GPU_SUPPORT_ENABLED`); flat `Array<int>` instead of `std::set`; device pack kernels |
| `WaveOperator::ApplyMassInverse` | `wave_operator.inl:5096` | `std::vector<DenseMatrix> elem_mass_inv_` | flat device `Vector` of inverses + `forall(ne)` (§2.2) |
| `WaveOperator::ApplySpatialDerivative` / ADER CK | `wave_operator.inl` (1900+/2200+/3300+) | `mutable Vector` scratch (`ader_*_buf_`, `ck_*_buf_`, `wave_operator.hpp:846–853`) lazy host alloc | `UseDevice(true)` on the buffers; reuse volume `forall`; ping-pong on device |
| `GodunovFlux::Interior/Central/Absorbing*/FreeSurface*` | `dynamic/godunov_flux.hpp:53–161`, `.cpp` | **member** `DenseMatrix Ax_,Ax_plus_,Ax_minus_,ref_star_[3]` (heap); methods are non-`HOST_DEVICE`; rotation builds `DenseMatrix T(9)` per call | mark `MFEM_HOST_DEVICE`; store the 9×9 matrices as `real_t[81]` POD passed by value/const-ref into the kernel; rotation into fixed `real_t[81]` stack arrays. *Signatures already use `const real_t*` — that part is GPU-friendly.* |
| `GodunovFlux::BuildRotation/Inverse/Frame`, `BuildJacobian` | `godunov_flux.hpp:190–227` | `static` but write into `DenseMatrix&`; called per face per QP | `MFEM_HOST_DEVICE` static writing into `real_t[81]`; or precompute per-face frames once on host and upload |
| `FaultFaceFlux::Evaluate / CompleteFromTrial / CompleteFromTheta / EvaluateADER` | `dynamic/fault_face_flux.{hpp:138–400,cpp}` | owns `FrictionSolver solver_` (`fault_face_flux.hpp:423`) called per fault QP; uses MFEM containers | `MFEM_HOST_DEVICE`; embed the friction solve (§3.3) inline; POD scratch; no virtual dispatch |
| `FrictionSolver::Solve/SolveBrent/SolveNR*` | `dynamic/friction_solver.{hpp:71–94,cpp}` | member `DieterichRuinaFriction qd_friction_`; `Brent` path uses `std::function`-style delegation; Brent has data-dependent iteration count → **warp divergence** | use the **fixed-iteration Newton** (`SolveNR`/`SolveNRStable`, already present) on device; mark `MFEM_HOST_DEVICE`; pure-scalar `real_t` args (already!) — but the *callees* in `friction/dieterich_ruina.hpp` & `friction/friction_coeff_stable.hpp` must also be `MFEM_HOST_DEVICE` and allocation-free |
| `FaultBasis` frame/embed (`EmbedSlip`, `GetBasis`, `qp_data`) | `fault/fault_basis.hpp`, used at `wave_operator.inl:2502–2603` | per-QP `DenseMatrix` rotations; `std::vector`-backed per-face basis | precompute per-fault-QP canonical frame + sign flag on host (already partly cached: `interior_fault_elem1_on_plus_`, `fault_interior_face_to_basis_idx_`, `wave_operator.hpp:949,975`) and upload as flat device arrays |
| `spatial_nucleation` `gradual_overstress` accumulator | `dynamic/spatial_nucleation.{hpp,cpp}` | per-DOF host loop writing `DOFData[i].tau{1,2}_nuc` each sub-step | `forall(num_fault_dofs)`; `MFEM_HOST_DEVICE` smoothStep+Gaussian; device `DOFData` SoA |

### 3.2 BP5 custom functions

| Function / class | File:line | Why it's host-only today | Required change |
|---|---|---|---|
| `ElasticityOperator::AssembleSlipContribution{IP,BR2}{,Shared}` | `elasticity_operator_assembly.inl:447–1019` | per-face loop; **fresh integrator allocated per face** (`:454`); BR2 raw QP tensor contractions with `elem_mass_inv_[e]` (`std::vector<DenseMatrix>`) | `forall` over fault faces; hoist integrator constants to device POD; flat device mass-inverse; `AtomicAdd` into RHS |
| `ElasticityOperator::AssembleDirichletLoading` | `elasticity_operator_assembly.inl:1067–1596` | same per-face hand-rolled lifting; time-dependent `g_D(x,t)` evaluated in-loop | pre-evaluate `g_D` on host once/step → flat array; `forall` |
| `ElasticityOperator::ComputeTractionImpl` | `elasticity_operator_traction.inl:653–1600` | per-face `DGElasticityIPCombinedIntegrator::ComputeTractionAtQuadPoints(Decomposed)`; `GetSubVector`; per-QP `DenseMatrix` | `MFEM_HOST_DEVICE` traction-at-QP kernel; `forall` over fault faces |
| `DGElasticityIPCombinedIntegrator::ComputeTractionAtQuadPoints(Decomposed)` | `dg_elasticity_ip_combined_integrator.hpp:434,565` | `DenseMatrix grad1/grad2`, `Vector shape1/shape2` temporaries per call | hoist into a free `MFEM_HOST_DEVICE` function with fixed-size scratch |
| `DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs` | `dg_elasticity_ip_combined_integrator.hpp:714` | static, but per-face `DenseMatrixInverse` (`:740`) | precompute the small face-projection inverse once per face on host; device apply |
| `DGElasticity*` / `BR2*` `AssembleFaceMatrix` | `integrator/*.hpp` | full-assembly only; many `DenseMatrix`/`Vector` temporaries; raw `GetColumn()/GetData()` | only needed on device if you pursue the §1.2 DG-face PA integrator; otherwise these stay CPU (stiffness assembled once at startup) |
| `ElasticityOperator::PrecomputeMassInverse` | `elasticity_operator_setup.inl:1202–1297` | builds `std::vector<DenseMatrix> elem_mass_inv_` (BR2 lifting) | startup-only; pack to device once if the slip/traction kernels run on device |
| `SEASQuasiDynamicOperator::Mult` friction RHS | `solver/seas_operator.hpp:318` | per-fault-DOF Brent loop (`fault_->ComputeRHS`) | fixed-iteration Newton in `forall(num_fault_dofs)` (same kernel family as §3.3) |

### 3.3 The friction solve — shared cross-cutting kernel

Both paths solve the same scalar rate-and-state balance per fault
point/DOF: `Θ = |σ_n|·f(V,ψ) + η·V`. The codebase **already** ships a
device-friendly variant:

- `FrictionSolver::SolveNR` / `SolveNRStable`
  (`dynamic/friction_solver.hpp:81,93`) — fixed-form Newton with an
  analytic derivative (`ResidualDerivative`, `:103`). All arguments are
  scalar `real_t` — no containers. This is the right device kernel body.
- `FrictionSolver::SolveBrent` (`:76`) is the CPU production path
  (data-dependent iteration count → warp divergence; delegates to
  `DieterichRuinaFriction`). **Keep Brent on host; route to fixed-iter
  Newton on device.** This is exactly the v3 Phase-2 `EvolutionKind`/
  method-dispatch decision.

> **Numerical guard (do not lose this):** project CLAUDE.md mandates
> Brent for the CPU friction solver (Newton fails for large ψ/a, debug
> v1). The device path must therefore be **validated** to the
> CPU↔GPU equivalence tolerance, *not* assumed equal. The fixed-iter
> Newton must reproduce the Brent root within tolerance across the full
> ψ/a range, including the degenerate `Fb≥0` frictionless-limit bracket
> guard. This is a verification requirement, not a free swap.

---

## 4. Call graphs (the PA / device decision points)

```
BP5 quasi-dynamic (per RK stage)
SEASQuasiDynamicOperator::Mult                      solver/seas_operator.hpp:318
  ├─ ElasticityOperator::Solve                      elasticity_operator_traction.inl:394
  │    ├─ AssembleSlipContribution{IP,BR2}          [hand-rolled face loop]  → forall  (3.2)
  │    ├─ AssembleDirichletLoading                  [hand-rolled face loop]  → forall  (3.2)
  │    └─ solver_->Mult (MUMPS/CG+AMG)              [needs ASSEMBLED matrix] ✗ PA-incompatible (1.3)
  ├─ ElasticityOperator::ComputeTraction            elasticity_operator_traction.inl:653 → forall (3.2)
  │    └─ DGElasticityIPCombined::ComputeTractionAtQuadPoints / ProjectTractionToFaultDOFs
  └─ fault_->ComputeRHS                             [per-DOF friction]       → forall (3.3)

   Stiffness K = ElasticityIntegrator(volume, PA-able)  ⊕  DG-face integrators(NO PA)
                 → PA blocked unless a DG-face PA integrator is authored (1.2)


Dynamic rupture (per RK4 stage)            ── no BilinearForm anywhere ──
WaveOperator::Mult                                  wave_operator.inl:854
  ├─ ComputeVolumeRHS        :942   → forall(ne)                         (2.1)
  ├─ ComputeFaceFluxRHS      :2254  → bucket + 4× forall + AtomicAdd     (2.1)
  │     └─ fault branch :2472 → FaultFaceFlux::Evaluate → FrictionSolver (3.1/3.3)
  │     └─ interior/bdry      → GodunovFlux::Interior/Absorbing/FreeSurface (3.1)
  ├─ ComputeSharedFaceFluxRHS:2923  → forall + GPU-aware MPI             (2.1)
  └─ ApplyMassInverse        :5096  → forall(ne) (or MFEM Mass PA)       (2.2)
```

---

## 5. Conventions & gotchas (must survive any port)

- **Tet/simplex meshes ⇒ no sum-factorization.** Both paths are tets
  (`elasticity_operator_setup.inl:11`, `wave_operator.inl:36`). MFEM PA on
  simplices uses dense `DofToQuad::FULL`; the big tensor-PA speedups
  (hex sum-factorization) do **not** apply. A hex mesh would change this,
  but that is a meshing decision with its own fault-geometry cost.
- **DG-face PA is absent for elasticity** — ours and MFEM's
  (`DGElasticityIntegrator`, `bilininteg.hpp:3702`). Scalar
  `DGDiffusionIntegrator` has it (`:3531`), so it is *possible* but
  unwritten. This is the single biggest blocker to literal BP5 PA.
- **BP5 solver needs the assembled matrix** (MUMPS_BLR default). PA
  removes the matrix; you must co-deliver a matrix-free solver +
  preconditioner. Don't enable PA "halfway."
- **`.GetData()` / `.HostRead()` are host-only.** Every hot-path read in
  both operators uses `Q.GetData()` (e.g. `wave_operator.inl:944`) or
  MFEM `Vector::operator()`; on device these must become `.Read()` /
  `.Write()` results fed into `forall`.
- **`std::vector<DenseMatrix>` / `DenseMatrixInverse` / fresh per-face
  integrators / `std::set` / `std::function`** all violate the device
  contract. Replace with flat device `Vector`/`Array<int>` and fixed-size
  `real_t` scratch.
- **Face scatter races.** A face writes both its elements' slots; two
  faces share an element. Device scatter ⇒ `mfem::AtomicAdd<real_t>`.
  CPU loops don't need it because they're serial.
- **Component-major Q layout** (`Q[c*ndof_total+dof]`) and **`byNODES`
  ordering** are load-bearing for both correctness and coalescing — keep
  them.
- **FP64 only.** Rate-and-state friction and the stress balance are
  FP64-sensitive (debug history). Do not introduce mixed precision in the
  friction or flux kernels without explicit validation.
- **Friction: Brent on host, Newton on device** — and *verify*
  equivalence, per CLAUDE.md's "Brent required" rule. Treat the device
  Newton as a new numerical path, not a drop-in.
- **Sign conventions** (slip ∥ traction, dip = +z, DG face-normal sign,
  σ_n>0 compression) live in scalar kernel bodies; porting must not touch
  them. See `miniapps/seas/CLAUDE.md` "Sign Conventions."

---

## 6. Recommended sequencing (synthesizes with v3 phases)

1. **Phase 0 (both):** `Device::Configure`, `Vector::UseDevice(true)`
   plumbing; no math change.
2. **Dynamic first** (it's already matrix-free, no solver entanglement):
   port `ComputeVolumeRHS` → `ApplyMassInverse` → `ComputeFaceFluxRHS`
   (bucketed) → friction Newton → ADER. = v3 Phases 3–5. This is where
   "PA-style" actually lands.
3. **BP5 solver offload** (no PA): GPU-HYPRE + `mumps`→`cg` so BoomerAMG
   runs on-device on the assembled matrix. = v3 Phase 1.
4. **BP5 coupling kernels** (slip RHS, Dirichlet, traction, friction) →
   `forall`. = v3 Phase 6. This is the §3.2 list.
5. **(Research, optional, last):** elasticity DG-face `AssemblePA` +
   LOR-AMG matrix-free preconditioner — the *only* route to literal BP5
   PA. Justify with profiling first.

---

## 7. Open questions

- **Is a hex mesh acceptable for either benchmark?** Tensor-PA is only a
  large win on hex/quad. The fault-conforming SAFS/BP5 geometry is
  currently tet (Tandem-style). Quantify the simplex-PA penalty before
  committing.
- **Does an LOR-AMG preconditioner converge for DG elasticity at BP5
  penalty scales?** Unverified. Needed before any literal BP5 PA.
- **Will fixed-iter Newton match Brent across the full ψ/a range** —
  including the `Fb≥0` degenerate bracket and nucleation-zone extremes —
  to the equivalence tolerance? Must be measured, not assumed.
- **`MassIntegrator` PA vs. custom element-inverse for the wave mass
  solve** — does routing through MFEM's PA mass + CG-Jacobi beat the
  explicit per-element inverse, given DG's block-diagonal mass? Likely
  no; confirm.
- **GPU-aware MPI on the target machine** (Frontier MPICH GPU support):
  is `ExchangeFaceNbrData` device-buffer-clean, or does it still bounce
  through host? Affects whether `ComputeSharedFaceFluxRHS` keeps its
  kernel speedup at scale.

---

### Appendix — file map referenced by this analysis

| Concern | Files |
|---|---|
| BP5 stiffness assembly + solver | `domain/elasticity_operator_assembly.inl`, `domain/elasticity_operator.hpp` |
| BP5 RHS / traction (hand-rolled) | `domain/elasticity_operator_assembly.inl`, `domain/elasticity_operator_traction.inl` |
| BP5 coupling driver | `solver/seas_operator.hpp`, `drivers/seas_driver.cpp` |
| Custom DG integrators (no PA) | `integrator/dg_elasticity_br2_integrator.hpp`, `integrator/dg_elasticity_ip_combined_integrator.hpp`, `integrator/dg_br2_integrator.hpp`, `integrator/dg_elasticity_ip_penalty_integrator.hpp` |
| Dynamic operator (matrix-free) | `dynamic/wave_operator.{hpp,inl}`, `drivers/spatial_dyn_driver.cpp` |
| Flux / friction / fault custom fns | `dynamic/godunov_flux.{hpp,cpp}`, `dynamic/fault_face_flux.{hpp,cpp}`, `dynamic/friction_solver.{hpp,cpp}`, `dynamic/spatial_nucleation.{hpp,cpp}`, `fault/fault_basis.hpp`, `dynamic/wave_state.hpp` |
| MFEM PA reference | `fem/integ/bilininteg_elasticity_pa.cpp` (volume PA), `fem/bilininteg.hpp` (`DGElasticityIntegrator` :3691 no PA; `DGDiffusionIntegrator` :3531 has face PA; `MassIntegrator` :2416 PA) |
| Existing GPU plans | `document/gpu_dev/gpu_implementation_plan_v{1,2,3}.md` |
