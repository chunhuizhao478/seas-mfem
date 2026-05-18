# Code Review: spatial_quasi_dynamic_plan.md (fresh review, 2026-05-18)

## Review Scope
- Plan: `safs/project_7.0_alternative/document/spatial_quasi_dynamic_plan.md` (742 lines, dated 2026-05-17)
- Companion: `safs/project_7.0_alternative/document/spatial_dynamic_rupture_plan.md` (rev-3, 2026-05-18) — checked for shared-structure cross-references
- Prereq: `safs/project_7.0_alternative/document/heterogeneous_material_plan.md` (and its `_review.md`)
- Datasets verified on disk:
  - `meshing/results/msh/safs_fault_box_nwcut_1000m_lcfar3000.msh` (PhysicalNames: 101 fault, 102 top, 103 bottom, 104 sides)
  - `velocity/results/{cvmh, cvm_s4.26.m01, multiscale_statewise_cvm}/velocity_safs.h5`
  - `stress/results/stress_csm_safs.h5`
  - `friction/{rate-and-state, slip-weakening}/` (rate-and-state directory is **empty**; slip-weakening has `geoffrey2010.md` + `131.pdf`)
- Source code consulted:
  - `fault/rate_state_fault.hpp` (BP5 ctor, `SetSAFSMode`, `PreInit`/`Init`/`ComputeRHS`)
  - `fault/fault_geometry.hpp` (BP5 ctor, `ComputeSAFSParams`, `NumZeroNormalFallbacks`)
  - `domain/elasticity_operator.hpp` (two existing ctors — **neither is `MaterialField`-aware**)
  - `domain/boundary_config.hpp` (`MakeBP5DirichletFunc` hard-codes the BP5 y-threshold)
  - `dynamic/heterogeneous_material.hpp` (`MaterialField::MakeConstant` / `MakeCoefficient` / `EvalAt`)
  - `solver/seas_operator.hpp` (`BP5SEASOp` = serial; `PBP5SEASOp` = parallel)
  - `io/petsc_ts_checkpoint.hpp` (free functions, NOT a class; appends to V1 .txt — **not** to `.h5`)
  - `io/paraview_output.hpp` (`InitFaultOutputBP5`, `SetFaultOutputMode`, `UpdateFaultFieldsBP5`)
  - `meshing/code/safs_fault_box_nwcut.geo` (defines Physical Surface 104 = {1,2,3,4} — all four sides merged)
- Domain context: `miniapps/seas/CLAUDE.md`, `project_7.0_alternative/CLAUDE.md`, `bp5_debug_document/`

---

## Findings

### [R-001] [CRITICAL] [Phase 4 step 5 / mesh] — Outer-boundary loading cannot be applied per-face: all four side faces are merged into one Physical Surface attribute

**Category:** DEVIATION

**Description:**
The user request states: *"we will load the domain at the outer boundaries (so be sure the boundary id is exposed)"*. The mesh in scope, `safs_fault_box_nwcut_1000m_lcfar3000.msh`, has only **four** surface attributes (verified at lines 4–11 of the .msh):

```
2 101 "fault"
2 102 "top"
2 103 "bottom"
2 104 "sides"
```

and `meshing/code/safs_fault_box_nwcut.geo:264` collapses all four lateral faces into a single tag:

```
Physical Surface("sides", 104)  = {1, 2, 3, 4};
```

The plan's Phase 4 step 5 says "Reuse existing `BoundaryConfig` from `safs_test_driver` (fault_attr = 101, dirichlet_attrs = sides + bottom, natural_attr = top — per the `safs_fault_box_nwcut.geo` Physical Surface tags above)" but with only attr 104 for all four sides, no per-face Dirichlet function can be installed in `bdr_config_.dirichlet_funcs[104]`. The only available mechanism is a coordinate-discriminating lambda (`MakeBP5DirichletFunc` style), which:

1. Hard-codes a BP5-style y-threshold (`y > 1000 m`, `y < -1000 m`) that is meaningless in UTM-11N coordinates (the SAFS mesh is centered near (480 km, 3750 km), NOT (0, 0)).
2. Forces the loading direction to be along strike of an idealised straight fault; the SAFS fault is curvilinear (multi-segment SAFS-MJVS-SBMT-MULT), so a single signed +X / −X loading is geometrically wrong.
3. Provides no mechanism for the user to expose each face's ID, which is exactly what they asked for.

**Trigger:**
Driver builds `BoundaryConfig` with `dirichlet_attrs = {103, 104}` and tries to apply a per-side loading.

**Actual behavior:**
The loading is either applied uniformly to ALL four sides (wrong sign on opposing faces, no farfield kinematic moment), or it has to be coordinate-coded with hard-coded UTM-11N thresholds that the plan never specifies.

**Expected behavior:**
The mesh must expose **four separate** Physical Surface attributes for the four side faces (e.g., 104 / 105 / 106 / 107 = `side_xmin / side_xmax / side_ymin / side_ymax`), so the driver can install a per-attribute `DirichletFunc` aligned with the regional plate motion. The plan must (a) document the per-side attribute scheme, (b) include a one-line mesh-generator change in `safs_fault_box_nwcut.geo` to split the sides, and (c) regenerate the .msh fixtures.

**Suggested fix:**

In `meshing/code/safs_fault_box_nwcut.geo` replace line 264:
```diff
- Physical Surface("sides", 104)  = {1, 2, 3, 4};
+ // Expose each lateral face separately so the SEAS driver can apply
+ // tectonic-loading Dirichlet BCs aligned with the regional plate motion.
+ Physical Surface("side_xmin", 104) = {1};
+ Physical Surface("side_xmax", 105) = {2};
+ Physical Surface("side_ymin", 106) = {3};
+ Physical Surface("side_ymax", 107) = {4};
```

Then add a Phase 0 sub-section to the plan: *Per-side boundary attribute mapping* that documents the four new attrs and the per-attribute `DirichletFunc` to install:

```cpp
// drivers/spatial_qd_driver.cpp — Phase 4 step 5
BoundaryConfig bc;
bc.fault_attr = 101;
bc.natural_attrs    = {102};                       // free surface
bc.dirichlet_attrs  = {103, 104, 105, 106, 107};   // bottom + 4 sides
// Tectonic far-field loading: half-rate on +y side, opposite half-rate on -y side.
// Coordinates are UTM-11N; plate motion direction is read from TOML
// ([loading].plate_motion_unit = [ux, uy, uz]).
auto plate = cfg.loading.plate_motion_unit;
real_t Vp  = cfg.loading.plate_rate_m_per_s;
bc.dirichlet_funcs[106] = [Vp, plate](const Vector& x, real_t t, Vector& u) {
   u.SetSize(3);
   u(0) = -0.5 * Vp * t * plate[0];
   u(1) = -0.5 * Vp * t * plate[1];
   u(2) = -0.5 * Vp * t * plate[2];
};
bc.dirichlet_funcs[107] = [Vp, plate](const Vector& x, real_t t, Vector& u) {
   u.SetSize(3);
   u(0) =  0.5 * Vp * t * plate[0];
   u(1) =  0.5 * Vp * t * plate[1];
   u(2) =  0.5 * Vp * t * plate[2];
};
// Side_xmin / side_xmax / bottom: clamp normal component only (or default_dirichlet_func = zero).
bc.default_dirichlet_func = [](const Vector&, real_t, Vector& u) { u.SetSize(3); u = 0.0; };
```

The plan must also enumerate which sides receive which sign — derive from `[loading].plate_motion_unit` rather than hard-coding.

**Test case:**
```cpp
void test_R001_per_side_boundary_attrs() {
   // Regenerate the mesh; load it; check the boundary attribute set.
   Mesh mesh("safs_fault_box_nwcut_1000m_lcfar3000.msh");
   std::set<int> bdr_attrs(mesh.bdr_attributes.begin(),
                            mesh.bdr_attributes.end());
   TEST_ASSERT(bdr_attrs == std::set<int>({101, 102, 103, 104, 105, 106, 107}));
   // For each of the 4 side attrs, assert ≥ 1 boundary face is tagged
   for (int attr : {104, 105, 106, 107}) {
      int count = 0;
      for (int be = 0; be < mesh.GetNBE(); ++be) {
         if (mesh.GetBdrAttribute(be) == attr) { ++count; }
      }
      TEST_ASSERT_MSG(count > 0, "side attr " << attr << " has zero faces");
   }
}
```

---

### [R-002] [CRITICAL] [Phase 4 step 17 + Constraints §"Restart V2 schema"] — `PETScTSCheckpoint::Restore("<PREFIX>/checkpoint.h5", ...)` is a fabricated API; actual code uses free functions appended to a per-rank text file

**Category:** ASSUMPTION

**Description:**
The plan repeatedly references a non-existent class API:

- Constraints §"Restart V2 schema" (`spatial_quasi_dynamic_plan.md:42`): "Driver writes `<output>/checkpoint.h5` with one HDF5 group per (rank × stride)".
- Phase 4 step 17 (`spatial_quasi_dynamic_plan.md:583`): `PETScTSCheckpoint::Restore("<PREFIX>/checkpoint.h5", state, t0, dt0)`.
- Phase 4 step 18: `PETScTSCheckpoint::Save(...)`.
- Phase 4 acceptance criterion: "V2 checkpoint round-trip: write at step 100, restart from that checkpoint" — also assumes a single `.h5` file.

The actual API in `io/petsc_ts_checkpoint.hpp` (verified) is:

```cpp
// FREE FUNCTIONS, not a class.  No `PETScTSCheckpoint::*` namespace.
inline void WritePetscTSCheckpoint(const std::string& prefix, real_t t, real_t dt_next,
                                   int step, int rejections,
                                   int paraview_snapshots, real_t paraview_last_write_time,
                                   real_t paraview_last_v_max, int paraview_current_regime,
                                   int paraview_last_committed_cycle,
                                   real_t paraview_last_volume_write_time,
                                   const MPIContext* mpi);

inline bool ReadPetscTSCheckpoint(const std::string& prefix, /*all of the same fields by ref*/);
```

The file written is **per-rank plain text**, named via `CheckpointFilename(prefix, rank)` (= `<prefix>_checkpoint_r<rank>.txt`, from `io/checkpoint.hpp:31`). The V2 block is APPENDED to a V1 BP5 checkpoint that `WriteCheckpoint` must have already produced; the V2 writer has an explicit pre-flight check that aborts if the V1 file is missing (`petsc_ts_checkpoint.hpp:94`).

There is no `.h5` file in the BP5 restart path. There is no `Save`/`Restore` class API. The plan's contract is fully fabricated.

Compounding the issue: the companion `spatial_dynamic_rupture_plan.md` already corrected this exact mistake as its R-105 ("`Tpv104Checkpoint::Restore` class API does not exist") — the QD plan must adopt the same correction.

**Trigger:**
Any compile attempt of Phase 4 driver code calling `PETScTSCheckpoint::Restore`.

**Actual behavior:**
Linker error; the symbol does not exist. Driver will not link.

**Expected behavior:**
Driver calls `WriteCheckpoint(prefix, ...)` then `WritePetscTSCheckpoint(prefix, ...)` per rank, and on restart calls `ReadCheckpoint(prefix, ...)` then `ReadPetscTSCheckpoint(prefix, ...)`. The output directory layout schema must replace `checkpoint.h5` with `cp_checkpoint_r<rank>.txt` (one file per rank per `--checkpoint-every`).

**Suggested fix:**

Replace Constraints §"Restart V2 schema" (`spatial_quasi_dynamic_plan.md:42`):
```diff
- - **Restart V2 schema** (`io/petsc_ts_checkpoint.hpp`) is fixed.  Driver writes `<output>/checkpoint.h5` with one HDF5 group per (rank × stride) — see existing `BP5 V2 + TPV104 V1 restart` commit (`63373f8`) for the contract.
+ - **Restart V2 schema** is fixed: free functions `Write/ReadCheckpoint` (V1, in `io/checkpoint.hpp`) + `Write/ReadPetscTSCheckpoint` (V2 trailing block, in `io/petsc_ts_checkpoint.hpp`).  Driver writes one **per-rank plain-text** file `<output>/<restart_prefix>_checkpoint_r<rank>.txt` per `--checkpoint-every`; the V2 block is appended after the V1 body.  PETSc-TS restart requires `MFEM_USE_PETSC=YES`; otherwise the V2 path is unavailable and the driver MUST fail at startup if `--restart` is set on a non-PETSc build.
```

Replace Phase 4 step 17:
```diff
- 17. **Restart path** — if `--restart PREFIX` set, call `PETScTSCheckpoint::Restore("<PREFIX>/checkpoint.h5", state, t0, dt0)` (existing).  Otherwise `t0 = 0, dt0 = computed from L_nuc / V_nuc`.
+ 17. **Restart path** — if `--restart PREFIX` set:
+     ```cpp
+     bool ok = ReadCheckpoint(restart_prefix, /* V1 fields */ ..., mpi);
+     MFEM_VERIFY(ok, "ReadCheckpoint failed for prefix " << restart_prefix);
+     bool ok2 = ReadPetscTSCheckpoint(restart_prefix, t0, dt0, step0, rejections0,
+                                       pv_snapshots, pv_last_write_time,
+                                       pv_last_v_max, pv_current_regime,
+                                       pv_last_committed_cycle,
+                                       pv_last_volume_write_time, mpi);
+     MFEM_VERIFY(ok2, "ReadPetscTSCheckpoint failed for prefix " << restart_prefix);
+     ```
+     The prefix is `<PREFIX>/<restart_prefix>` from `[output].restart_prefix` (default "cp").  Otherwise `t0 = 0, dt0 = computed from L_nuc / V_nuc`.
```

Replace Phase 4 step 18 last bullet:
```diff
- - Every `checkpoint_every_steps` steps, write V2 PETSc-TS checkpoint via `PETScTSCheckpoint::Save(...)`.
+ - Every `checkpoint_every_steps` steps, call `WriteCheckpoint(prefix, ...)` then `WritePetscTSCheckpoint(prefix, ...)` per rank (V1 first, then V2 appended).
```

Replace output directory layout (`spatial_quasi_dynamic_plan.md:62`):
```diff
-     checkpoint.h5            (V2 PETSc-TS checkpoint; one per --checkpoint-every)
+     cp_checkpoint_r<rank>.txt  (V1 + V2 per-rank checkpoint; one PER RANK per --checkpoint-every)
```

Replace Phase 4 acceptance criterion:
```diff
- - [ ] V2 checkpoint round-trip: write at step 100, restart from that checkpoint, run 100 more steps; assert state at step 200 (continuous) == state at step 200 (restarted) to within `||state||_inf * 1e-12`.
+ - [ ] V1+V2 per-rank checkpoint round-trip: at step 100 call `WriteCheckpoint` + `WritePetscTSCheckpoint`, restart by calling `ReadCheckpoint` + `ReadPetscTSCheckpoint`, run 100 more steps; assert state at step 200 (continuous) == state at step 200 (restarted) to within `||state||_inf * 1e-12`.  Verify on np=1 and np=4 (one .txt file per rank).
```

**Test case:**
```cpp
void test_R002_checkpoint_api_exists() {
   // Compile-time check: the V2 free functions exist with the expected signatures.
   using WriteFn = void(*)(const std::string&, real_t, real_t, int, int,
                            int, real_t, real_t, int, int, real_t,
                            const MPIContext*);
   using ReadFn  = bool(*)(const std::string&, real_t&, real_t&, int&, int&,
                            int&, real_t&, real_t&, int&, int&, real_t&,
                            const MPIContext*);
   WriteFn w = &mfem::seas::WritePetscTSCheckpoint;
   ReadFn  r = &mfem::seas::ReadPetscTSCheckpoint;
   TEST_ASSERT(w && r);
}
void test_R002_checkpoint_v2_requires_v1_prefix() {
   // The V2 writer aborts if the V1 file is missing (per
   // petsc_ts_checkpoint.hpp:94).  Confirm the abort path fires.
   ScopedTempDir tmp;
   EXPECT_ABORT(mfem::seas::WritePetscTSCheckpoint(
      tmp.path() + "/cp", 0.0, 1e-3, 0, 0, 0, 0.0, 0.0, 0, -1, 0.0, nullptr));
}
```

---

### [R-003] [CRITICAL] [Phase 2 §Detailed Requirements 1] — Velocity sidecar path resolution is wrong: actual filename has NO mesh-tag suffix

**Category:** DEVIATION

**Description:**
Phase 2 step 1 (`spatial_quasi_dynamic_plan.md:415`) prescribes:

> - `VelocityModel::CVMH` → `<dataset_root>/velocity/results/cvmh/velocity_safs_<mesh_tag>.h5`

But the actual files on disk are (verified):
```
velocity/results/cvmh/velocity_safs.h5
velocity/results/cvm_s4.26.m01/velocity_safs.h5         (need to confirm)
velocity/results/multiscale_statewise_cvm/velocity_safs.h5  (need to confirm)
```

No `mesh_tag` suffix. The dynamic-rupture plan corrected this exact bug as its R-003 ("R-003 corrected: actual on-disk filename has NO mesh_tag suffix; one sidecar per CVM model"). The QD plan still carries the wrong scheme.

The TOML in Phase 0 also exposes `[velocity].mesh_tag = "1000m_lcfar3000"` which has no purpose if the filename does not encode the tag — keep it only if used as a metadata sanity check, otherwise drop it entirely.

**Trigger:**
Phase 2 driver init looks for `<root>/velocity/results/cvmh/velocity_safs_1000m_lcfar3000.h5` and aborts at the `MFEM_VERIFY` file-exists pre-flight.

**Actual behavior:**
Driver aborts with "file not found at velocity_safs_1000m_lcfar3000.h5".

**Expected behavior:**
Driver loads `velocity_safs.h5` from the model-specific subdirectory.

**Suggested fix:**

Replace Phase 2 step 1:
```diff
- - `VelocityModel::CVMH` → `<dataset_root>/velocity/results/cvmh/velocity_safs_<mesh_tag>.h5`
- - `VelocityModel::CVMS_4_26_M01` → `<dataset_root>/velocity/results/cvm_s4.26.m01/velocity_safs_<mesh_tag>.h5`
- - `VelocityModel::MultiscaleStatewise` → `<dataset_root>/velocity/results/multiscale_statewise_cvm/velocity_safs_<mesh_tag>.h5`
+ - `VelocityModel::CVMH` → `<dataset_root>/velocity/results/cvmh/velocity_safs.h5`
+ - `VelocityModel::CVMS_4_26_M01` → `<dataset_root>/velocity/results/cvm_s4.26.m01/velocity_safs.h5`
+ - `VelocityModel::MultiscaleStatewise` → `<dataset_root>/velocity/results/multiscale_statewise_cvm/velocity_safs.h5`
+ - **No `mesh_tag` in the filename** (sidecars are mesh-independent; the trilinear interpolator handles any mesh inside its bbox).
```

Drop `mesh_tag` from `SpatialVelocitySpec` (Phase 2 interfaces) unless it is repurposed as a metadata-only sanity check (in which case document the new use).

Update Phase 2 TOML mapping (line 423):
```diff
  [velocity]
  model         = "cvmh"        # cvmh | cvm_s4.26.m01 | multiscale_statewise
- mesh_tag      = "1000m_lcfar3000"
+ dataset_root  = "safs/project_7.0_alternative"
  override_path = ""            # set non-empty to bypass model+tag resolution
```

**Test case:**
```cpp
void test_R003_velocity_path_no_mesh_tag() {
   VelocitySpec spec;
   spec.model = VelocityModel::CVMH;
   spec.dataset_root = "safs/project_7.0_alternative";
   std::string p = ResolveSpatialVelocitySidecarPath(spec);
   TEST_ASSERT_EQ(p, "safs/project_7.0_alternative/velocity/results/cvmh/velocity_safs.h5");
}
void test_R003_velocity_file_exists_on_disk() {
   // Smoke: the resolved path actually exists in the repo.
   VelocitySpec spec; spec.model = VelocityModel::CVMH;
   spec.dataset_root = "safs/project_7.0_alternative";
   std::ifstream f(ResolveSpatialVelocitySidecarPath(spec));
   TEST_ASSERT(f.good());
}
```

---

### [R-004] [CRITICAL] [Phase 4 step 12 + Constraints] — `BP5SEASOp` template alias is the SERIAL variant; parallel runs need `PBP5SEASOp`

**Category:** BUG

**Description:**
Constraints §"Interface constraints" (`spatial_quasi_dynamic_plan.md:40`) and Phase 4 step 12 (line 573) both refer to the alias `BP5SEASOp` for the time-integration coupling shell. Verified in `solver/seas_operator.hpp`:

```cpp
using BP5SEASOp   = SEASQuasiDynamicOperator<Mesh, BP5DomainOp, BP5FaultOp>;        // line 551 — SERIAL

#ifdef MFEM_USE_MPI
using PBP5SEASOp = SEASQuasiDynamicOperator<ParMesh,                                 // line 555 — PARALLEL
                      ElasticityDomainOperator<ParMesh>,
                      RateStateFaultOperator<ParMesh, 2>>;
#endif
```

`BP5SEASOp` instantiates over `Mesh` (serial); `PBP5SEASOp` instantiates over `ParMesh` (parallel). The plan everywhere assumes the driver operates on a `ParMesh` (Phase 4 step 5 explicitly: "`ParMesh pmesh(MPI_COMM_WORLD, smesh);`") AND uses the alias `BP5SEASOp` — these contradict.

Driver code written to the plan as stated will either fail to instantiate (`SEASQuasiDynamicOperator<Mesh, ...>` won't accept a `ParMesh`-based domain/fault operator pair) OR will silently downgrade the parallel domain operator to a serial one, breaking MPI.

**Trigger:**
Phase 4 driver attempts to construct `BP5SEASOp seas_op(&domain_op, &fault_op)` where `domain_op` is `ElasticityDomainOperator<ParMesh>`.

**Actual behavior:**
Compile error (template argument deduction failure: serial `BP5DomainOp = ElasticityDomainOperator<Mesh>` does not match the parallel `ElasticityDomainOperator<ParMesh>`).

**Expected behavior:**
Driver uses `PBP5SEASOp` (guarded by `#ifdef MFEM_USE_MPI`).

**Suggested fix:**

Replace Constraints §"Interface constraints" line 40:
```diff
- - **`SEASQuasiDynamicOperator<ParMesh, ElasticityDomainOperator<ParMesh>, BP5FaultOp>`** template (alias `BP5SEASOp` in `solver/seas_operator.hpp:551`) is the time-integration coupling shell; do not subclass.
+ - **`SEASQuasiDynamicOperator<ParMesh, ElasticityDomainOperator<ParMesh>, RateStateFaultOperator<ParMesh, 2>>`** template (alias `PBP5SEASOp` in `solver/seas_operator.hpp:555`, guarded by `#ifdef MFEM_USE_MPI`) is the time-integration coupling shell; do not subclass.  The serial alias `BP5SEASOp` (line 551) is unused by this driver (the SAFS workflow is parallel-only — see Phase 6 sbatch job sizes).
```

Replace Phase 4 step 12:
```diff
- 12. **Construct `BP5SEASOp seas_op(&domain_op, &fault_op);`** (alias for `SEASQuasiDynamicOperator<ParMesh, ElasticityDomainOperator<ParMesh>, BP5FaultOp>`).
+ 12. **Construct `PBP5SEASOp seas_op(&domain_op, &fault_op);`** (parallel alias from `solver/seas_operator.hpp:555`).  This driver does not support `MFEM_USE_MPI=NO` builds.
```

**Test case:**
```cpp
void test_R004_parallel_seas_op_alias_compiles() {
   // Pure compile-time check.
   using DomT = mfem::seas::ElasticityDomainOperator<mfem::ParMesh>;
   using FltT = mfem::seas::RateStateFaultOperator<mfem::ParMesh, 2>;
   static_assert(std::is_same_v<mfem::seas::PBP5SEASOp,
                                mfem::seas::SEASQuasiDynamicOperator<
                                   mfem::ParMesh, DomT, FltT>>,
                 "PBP5SEASOp must be the parallel alias");
}
```

---

### [R-005] [CRITICAL] [Constraints §"Dependency constraints" + Phase 4 step 7] — Hard dependency on `MaterialField`-aware `ElasticityDomainOperator` ctor that does NOT exist; fallback path is wrong

**Category:** ASSUMPTION

**Description:**
Constraints §"Interface constraints" (`spatial_quasi_dynamic_plan.md:37`) asserts:

> **`ElasticityDomainOperator<ParMesh>` two existing constructors** (legacy `(λ, μ)` and `MaterialField`-aware) remain bit-exact for non-SAFS callers.

Verified `domain/elasticity_operator.hpp:91-192`: only two ctors exist, the `ConstitutiveModel`-based one (line 105) and the legacy scalar-`(λ, μ)` one (line 142). **Neither accepts `MaterialField`.** The `MaterialField`-aware ctor is in `heterogeneous_material_plan.md` Phase 2, which has NOT been merged (verified: no `MaterialField` references anywhere in `domain/elasticity_operator.hpp`).

The plan's fallback path in Constraints §"Dependency constraints" line 47:

> If the plan is not yet implemented, this driver is implementable but cannot run heterogeneous; in that case it falls back to `MaterialField::MakeConstant(mu_const, lambda_const, rho_const)` with values from the friction-config TOML's `[material_constant_fallback]` section.

This is **logically wrong**: `MaterialField::MakeConstant` returns a `MaterialField`. The plan's fallback hands that `MaterialField` to a constructor that does not accept it. The fallback only works if the `MaterialField`-aware ctor exists; the conditional is empty.

Moreover, the `ElasticityDomainOperator` ctor argument order in `MakeConstant` is `(lambda, mu, rho)` (verified in `dynamic/heterogeneous_material.hpp:94`), but the plan calls it `MakeConstant(mu_const, lambda_const, rho_const)` — **argument order is swapped** (mu first vs lambda first).

**Trigger:**
Driver attempts to construct `ElasticityDomainOperator<ParMesh>(mesh, order, material, ...)` on the current main branch (heterogeneous_material plan not merged).

**Actual behavior:**
Compile error: no matching constructor.

**Expected behavior:**
Either (a) treat heterogeneous_material plan Phase 2 as a **hard prerequisite** (mirror the dynamic-rupture plan's "Phase H" pattern that owns the work in-house, or block this plan until Phase 2 merges), or (b) build the `MaterialField` and unconditionally use the new ctor, with a clean error message when the ctor is missing.

**Suggested fix:**

Replace Constraints §"Dependency constraints" first bullet:
```diff
- - **`heterogeneous_material_plan.md` Phase 1+2 must be merged first** — the driver uses `MaterialField::MakeCoefficient(...)` and the `ElasticityDomainOperator(MaterialField)` constructor.  If the plan is not yet implemented, this driver is implementable but cannot run heterogeneous; in that case it falls back to `MaterialField::MakeConstant(mu_const, lambda_const, rho_const)` with values from the friction-config TOML's `[material_constant_fallback]` section.
+ - **`heterogeneous_material_plan.md` Phase 1+2 must be merged first.**  Phase 1 (`MaterialField::MakeCoefficient`) is already merged on this branch; Phase 2 (the `ElasticityDomainOperator(MaterialField)` ctor) is NOT.  This plan REQUIRES Phase 2; the driver cannot link until Phase 2 merges, regardless of whether the user selects sidecar or constant material.  When Phase 2 merges:
+   - **Sidecar path**: `material = vbundle.MakeMaterialField();` (`Mode::Coefficient`).
+   - **Constant fallback**: `material = MaterialField::MakeConstant(lambda, mu, rho);` — note the canonical argument order `(lambda, mu, rho)` matches `dynamic/heterogeneous_material.hpp:94`.  Reading `cfg.material_fallback.{lambda, mu, rho}`.
+
+   A compile-time guard at the top of `drivers/spatial_qd_driver.cpp` MUST verify the new ctor exists:
+   ```cpp
+   static_assert(std::is_constructible_v<
+                    mfem::seas::ElasticityDomainOperator<mfem::ParMesh>,
+                    mfem::ParMesh&, int, const mfem::seas::MaterialField&,
+                    mfem::real_t, mfem::real_t, mfem::real_t,
+                    const mfem::seas::BoundaryConfig&,
+                    mfem::seas::DGMethod, mfem::seas::SolverType,
+                    const mfem::seas::DomainConfig&>,
+                 "spatial_qd_driver requires heterogeneous_material_plan "
+                 "Phase 2 (ElasticityDomainOperator(MaterialField) ctor).");
+   ```
+
+   Alternative architectural decision (mirrors dynamic-rupture plan's Phase H): own the Phase-2 ctor work IN THIS PLAN as a new "Phase E" so the QD driver is buildable end-to-end without the prereq.  The dynamic-rupture plan already chose this path for the analogous WaveOperator(MaterialField) ctor (rev-3 D-2 decision).  See §Open Question OQ-1 below.
```

Replace Phase 4 step 6 last sentence:
```diff
- In `--no-sidecar-material` (added below) or if velocity is missing, fall back to `MaterialField::MakeConstant(toml.material_fallback.lambda, ..., toml.material_fallback.rho)` with a one-line `mfem::out` notice.
+ In `--no-sidecar-material` (added below) or if velocity is missing, fall back to `MaterialField::MakeConstant(toml.material_fallback.lambda, toml.material_fallback.mu, toml.material_fallback.rho)` (canonical lambda-first order per `dynamic/heterogeneous_material.hpp:94`) with a one-line `mfem::out` notice.
```

**Test case:**
```cpp
void test_R005_material_field_ctor_exists() {
   // Compile-time only — the driver cannot link if this fails.
   static_assert(std::is_constructible_v<
                    mfem::seas::ElasticityDomainOperator<mfem::ParMesh>,
                    mfem::ParMesh&, int, const mfem::seas::MaterialField&,
                    mfem::real_t, mfem::real_t, mfem::real_t,
                    const mfem::seas::BoundaryConfig&,
                    mfem::seas::DGMethod, mfem::seas::SolverType,
                    const mfem::seas::DomainConfig&>,
                 "MaterialField ctor missing");
}
void test_R005_material_field_constant_argument_order() {
   auto m = mfem::seas::MaterialField::MakeConstant(/*lambda=*/1.0,
                                                    /*mu=*/2.0,
                                                    /*rho=*/3.0);
   TEST_ASSERT_EQ(m.lambda_const, 1.0);
   TEST_ASSERT_EQ(m.mu_const, 2.0);
   TEST_ASSERT_EQ(m.rho_const, 3.0);
}
```

---

### [R-006] [CRITICAL] [Phase 4 step 8 + 10] — `FaultGeometry` BP5 ctor unconditionally calls `ComputeBP5Params()`; using `BP5Params{}` as a "seed" runs spurious BP5 a(z)/Dc(z) maps with default-zero parameters

**Category:** BUG

**Description:**
Phase 4 step 8 (`spatial_quasi_dynamic_plan.md:565`):

> Construct `FaultGeometry<ParMesh> geom(domain_op, bp5_params_seed)` where `bp5_params_seed` is a placeholder `BP5Params{}` — the per-DOF rate-state values from Phase 1 override its scalars before any consumption.

Verified `fault/fault_geometry.hpp:118-160`: the BP5 ctor body ends with:

```cpp
   // Phase 6.A — per-DOF global 3-D coordinates and (n, t1, t2) basis.
   ComputePerDOFCoordsAndBasis_(domain_op);

   // Precompute per-DOF parameters using BP5 2D functions
   ComputeBP5Params();
}
```

`ComputeBP5Params()` is called **unconditionally** with the supplied `bp5_params_`. With a default-constructed `BP5Params{}` (all-zero `a, b, Dc, Wf, h, V_init, ...`), one of two things happens:

1. The BP5 spatial functions (in `config/bp5_params.hpp`) divide by zero / take log of zero / index out of bounds → silent NaN propagation or crash.
2. They succeed and produce all-zero per-DOF vectors; the `tau_pre_` and downstream sigma_n vectors are seeded to zero.

Then the plan calls `ApplySpatialStressSidecar` (step 9), which calls `geom.ComputeSAFSParams(...)` that **overwrites** `tau_pre_` and `sigma_n_per_dof_` from the stress sidecar (`fault_geometry.hpp:1049-1066`). But it does NOT overwrite the per-DOF `a_values_`, `Dc_values_`, etc. — those remain at the BP5-spatial-function-from-zero state until the Phase-5 setter `SetPerDOFRateStateParams` is called inside `RateStateFaultOperator`. Between Phase 4 step 8 and step 11, those per-DOF vectors are in a corrupted state and accessible via `geom.a_values()` / `IsVelocityWeakening(i)` etc.

Beyond corrupted intermediate state, `BP5Params{}` validation may fail inside `ComputeBP5Params`: e.g. `Dc = 0` is rejected by the BP5 spatial law.

**Trigger:**
Driver step 8 with default-constructed `BP5Params{}`.

**Actual behavior:**
Either (a) `ComputeBP5Params()` aborts on a zero-denominator (most likely) or (b) silently fills corrupted per-DOF vectors that are observable via the public API before Phase-5 setters fire.

**Expected behavior:**
A new ctor overload that does NOT call `ComputeBP5Params()`, OR a `BP5Params` constructed from `cfg` defaults that produces self-consistent (if irrelevant) per-DOF vectors that are then unconditionally overridden.

**Suggested fix (two acceptable paths; recommend Option B):**

**Option A** — add a `FaultGeometry` ctor overload that skips `ComputeBP5Params()`:

In `fault/fault_geometry.hpp`, add (additive, no impact on BP5 drivers):
```cpp
/// SAFS-mode ctor: skip the BP5 per-DOF parameter precomputation; the
/// caller is responsible for filling per-DOF a/b/Dc/... via downstream
/// resolvers (e.g., SpatialFrictionResolver + SetPerDOFRateStateParams).
struct SkipBP5ParamsTag {};
FaultGeometry(DomainOperator<MeshType>& domain_op, SkipBP5ParamsTag,
              MPIContext* mpi_ctx = nullptr);
```

Driver code becomes:
```cpp
FaultGeometry<ParMesh> geom(domain_op, FaultGeometry<ParMesh>::SkipBP5ParamsTag{}, mpi.get());
```

**Option B** (recommended — smaller surface): construct `BP5Params` from `cfg` defaults so `ComputeBP5Params` produces valid (though not-yet-final) intermediate state, then override per-DOF via Phase 5 setters:

Update Phase 4 step 8:
```diff
- 8. **Construct `FaultGeometry<ParMesh> geom(domain_op, bp5_params_seed)`** where `bp5_params_seed` is a placeholder `BP5Params{}` — the per-DOF rate-state values from Phase 1 override its scalars before any consumption.  The constructor unconditionally calls `ComputePerDOFCoordsAndBasis_` (REVIEW R-001 guard ensures any zero-normal fallback raises a loud warning at this point).
+ 8. **Construct `FaultGeometry<ParMesh> geom(domain_op, bp5_params_seed)`** where `bp5_params_seed` is built from the resolved scalar defaults in the friction TOML so `FaultGeometry::ComputeBP5Params` produces non-degenerate intermediate per-DOF a/b/Dc/V_init/Wf/h vectors:
+   ```cpp
+   BP5Params bp5_params_seed;
+   bp5_params_seed.a       = cfg.rate_state->a_default;       // > 0
+   bp5_params_seed.b       = cfg.rate_state->b_default;       // > a_default
+   bp5_params_seed.Dc      = cfg.rate_state->Dc_default;      // > 0
+   bp5_params_seed.V_init  = cfg.rate_state->V_init_default;  // > 0
+   bp5_params_seed.f_0     = cfg.rate_state->f_0_default;
+   bp5_params_seed.V_0     = cfg.rate_state->V_0_default;
+   bp5_params_seed.sigma_n = cfg.rate_state->sigma_n_default;
+   bp5_params_seed.eta     = /*scalar fallback*/ ... ;
+   // Wf / lf inferred from mesh bbox (see Phase 4 step 5 helper).
+   ```
+   The Phase-5 `SetPerDOFRateStateParams` call (step 11) replaces these scalar-derived per-DOF vectors with the spatially-resolved ones.  The intermediate state is self-consistent but irrelevant.
+
+   The constructor unconditionally calls `ComputePerDOFCoordsAndBasis_` (REVIEW R-001 guard ensures any zero-normal fallback raises a loud warning at this point) AND `ComputeBP5Params` (which is now safe because the seed values are all > 0).
```

**Test case:**
```cpp
void test_R006_seed_bp5_params_safe() {
   // FaultGeometry ctor must not crash or NaN-propagate when seeded
   // from typical TOML defaults.
   auto mesh = LoadSmallSAFSMesh();
   auto pmesh = MakeParMesh(mesh);
   ElasticityDomainOperator<ParMesh> dom(pmesh, /*order=*/1, ...);
   BP5Params seed;
   seed.a = 0.010; seed.b = 0.015; seed.Dc = 0.004;
   seed.V_init = 1e-9; seed.f_0 = 0.6; seed.V_0 = 1e-6;
   seed.sigma_n = 50e6; seed.eta = 1e6; seed.Wf = 40e3; seed.h = 5e3;
   FaultGeometry<ParMesh> geom(dom, seed);
   // No abort; per-DOF vectors are finite.
   for (int i = 0; i < geom.NumFaultDOFs(); ++i) {
      TEST_ASSERT(std::isfinite(geom.GetAValues()(i)));
      TEST_ASSERT(geom.GetDcValues()(i) > 0.0);
   }
}
```

---

### [R-007] [CRITICAL] [Phase 5 + Phase 4 step 10] — `ElasticityDomainOperator::GetFaultElemOwners()` accessor is asserted to exist with trivial internal state; no such map is verified to exist in the operator today

**Category:** ASSUMPTION

**Description:**
Phase 5 first bullet (`spatial_quasi_dynamic_plan.md:622`) claims:

> **`domain/elasticity_operator.hpp`**: add
> ```cpp
> const Array<int>& GetFaultElemOwners() const;
> ```
> populated alongside the existing per-DOF coords/basis computation.

And Phase 5 detailed requirement 1 (line 649):

> Already-existing internal `fault_face_to_elem` map suffices; just expose.

Verified `domain/elasticity_operator.hpp:248-298` — the public accessors expose `GetNumFaultDOFs`, `GetFaultDepths`, `GetFaultCoords2D`, `GetFaultDOFCoords3D`, `GetFaultDOFBasis`, `GetFaultDOFs`, `GetFaultBasis`, `GetFaultInteriorFaces`, `GetFaultSharedFaces`, `GetOwnedFaultFaceMap`, `GetNumOwnedFaultFaces`. Verified private members at line 469-471: only `fault_dofs_` and `num_fault_dofs_` are visible by `grep`. There is **no visible `fault_face_to_elem`** map. The plan's claim "Already-existing... just expose" cannot be verified.

The Phase 4 step 10 consumer (`spatial_quasi_dynamic_plan.md:569`) calls:

> `ElasticityDomainOperator::GetFaultElemOwners()` — a small new accessor added in Phase 5

If the underlying data structure does not exist, Phase 5 is not "trivial" — it requires extending `ElasticityDomainOperator`'s per-DOF / per-face bookkeeping, possibly touching face-loop assembly code. CLAUDE.md flags `domain/elasticity_operator.hpp` as "extreme care".

**Trigger:**
Phase 5 implementation starts.

**Actual behavior:**
Implementer discovers there is no internal `fault_face_to_elem` map and must walk the mesh face iterator + face-to-element table from scratch.

**Expected behavior:**
Plan must (a) verify the underlying map exists with a code citation, or (b) include the population logic explicitly (one line per fault face: for each `fault_interior_face_idx`, get `elem1` from `mesh.GetFaceElements(face_idx, &e1, &e2)`, then for each owned DOF on that face record `e1`).

**Suggested fix:**

Replace Phase 5 §"Detailed Requirements" item 1:
```diff
- 1. **`GetFaultElemOwners` population**: extend `ElasticityDomainOperator::InitFaultMaps` to record `fault_elem_owners_[i] = local_elem_index_of_face_elem1(fault_face_index_of_dof_i)`.  Already-existing internal `fault_face_to_elem` map suffices; just expose.
+ 1. **`GetFaultElemOwners` population**: there is no pre-existing `fault_face_to_elem` map on `ElasticityDomainOperator` (verified 2026-05-18).  Phase 5 must:
+   a. Add a private `Array<int> fault_elem_owners_` member.
+   b. Inside the existing `InitOperator()` / per-DOF-coords-and-basis computation (where `fault_dofs_` and `fault_interior_faces_` / `fault_shared_faces_` are populated), record for each owned fault DOF the `elem1` index of its host face:
+      ```cpp
+      fault_elem_owners_.SetSize(num_owned_fault_dofs_);
+      int dof_cursor = 0;
+      for (int face_idx : fault_interior_faces_) {
+         int e1, e2;
+         mesh_.GetFaceElements(face_idx, &e1, &e2);
+         for (int j = 0; j < nbf_per_face_; ++j) {
+            fault_elem_owners_[dof_cursor++] = e1;
+         }
+      }
+      // shared faces: e1 is the local-side element
+      for (int face_idx : fault_shared_faces_) {
+         int e1, e2;
+         pmesh_.GetSharedFaceElements(face_idx, &e1, &e2);
+         for (int j = 0; j < nbf_per_face_; ++j) {
+            fault_elem_owners_[dof_cursor++] = e1;
+         }
+      }
+      MFEM_VERIFY(dof_cursor == num_owned_fault_dofs_,
+                  "GetFaultElemOwners: DOF count mismatch");
+      ```
+   c. Expose `const Array<int>& GetFaultElemOwners() const { return fault_elem_owners_; }`.
+
+   Each step touches `elasticity_operator.hpp` (a CLAUDE.md "extreme care" file).  All four existing BP5 tests must run byte-identically; the new code path is consumed only by SAFS-side callers.
```

**Test case:**
```cpp
void test_R007_get_fault_elem_owners() {
   auto mesh = LoadSmallSAFSMesh();
   ParMesh pmesh(MPI_COMM_WORLD, mesh);
   ElasticityDomainOperator<ParMesh> dom(pmesh, 1, /*model=*/..., ...);
   const Array<int>& owners = dom.GetFaultElemOwners();
   TEST_ASSERT_EQ(owners.Size(), dom.GetNumOwnedFaultDOFs());
   // Every owner is a valid local element index.
   for (int i = 0; i < owners.Size(); ++i) {
      TEST_ASSERT(owners[i] >= 0 && owners[i] < pmesh.GetNE());
   }
}
void test_R007_bp5_owners_match_face_elem1() {
   // Build BP5 mesh; for each owned fault DOF, verify owner == elem1 of its host face.
   ...
}
```

---

### [R-008] [CRITICAL] [Phase 4 steps 13–14] — Driver `--dry-run` runs the multi-GB MUMPS_BLR factorization and full domain solve before exiting

**Category:** BUG

**Description:**
Phase 4 sequencing (`spatial_quasi_dynamic_plan.md:573-580`):

> 13. **Build initial state vector** — `Vector state(seas_op.Width()); fault_op.PreInit(state);` ... Then `seas_op.Solve(0.0, slip, u);` once to get domain-equilibrium traction.  Then `fault_op.Init(state)` computes ψ from equilibrium per Phase 3 of the BP5 4-phase init.
> 14. **Equilibrium verification** — assert `max_i |F(V_i, ψ_i)| < 1e-6 * V_init_max` across all owned DOFs (MPI_Allreduce); abort otherwise with a per-rank list of offending DOFs.
> 15. **If `--dry-run`** — print the assembled summary (...) and exit 0.

`ElasticityDomainOperator::Solve(0.0, slip, u)` (step 13) requires the stiffness matrix to be assembled AND the MUMPS_BLR factorization to be computed. On the 1000m_lcfar3000 mesh (≈ 3 M tets), the factorization is a multi-minute, multi-GB operation. Putting `--dry-run` AFTER step 13 defeats its stated purpose: "Parse + construct + validate equilibrium; do not time-step" (CLI flag description, line 542).

The user reasonable expectation for `--dry-run` is "quick sanity check; print derived quantities; exit". Including the full domain Solve and equilibrium ψ-solve adds 10–60 minutes of wall and ≈ 100 GB of memory — basically as expensive as the first time step.

**Trigger:**
`seas_spatial_qd_driver --config foo.toml --dry-run`.

**Actual behavior:**
Driver runs the full MUMPS factorization + equilibrium solve, then exits.

**Expected behavior:**
`--dry-run` mode should construct all operators, validate inputs, print derived quantities (L_nuc, Δt_initial, sidecar bbox containment, friction-parameter histograms), and exit BEFORE the expensive `Solve` + equilibrium ψ computation.

**Suggested fix:**

Restructure Phase 4 steps 13–15:
```diff
+ 12.5 **If `--dry-run`** — at this point the mesh, material, fault geometry, per-DOF friction params, and pre-stress are all loaded and validated.  Print:
+      - mesh stats (num elems, num bdr per attr, dim, np)
+      - velocity-sidecar bbox vs mesh bbox
+      - L_nuc per VW region (using per-DOF a/b/Dc/sigma_n_eff)
+      - dt_initial estimate
+      - friction-parameter histograms (a, b, Dc, sigma_n_eff, V_init)
+      - controlling/min/max element size + c_p
+      Exit 0 BEFORE the expensive Solve + equilibrium init.

  13. **Build initial state vector** — `Vector state(seas_op.Width()); fault_op.PreInit(state);` ... Then `seas_op.Solve(0.0, slip, u);` once to get domain-equilibrium traction.  Then `fault_op.Init(state)` ...

  14. **Equilibrium verification** — assert `max_i |F(V_i, ψ_i)| < 1e-6 * V_init_max` ...

- 15. **If `--dry-run`** — print the assembled summary (`L_nuc` per region, mesh stats, MPI layout, output paths) and exit 0.
+ 15. **If `--equilibrium-only`** (new flag, optional) — after equilibrium verification (step 14), print the per-region max |F| / V_init_max and exit 0.  Use when the user wants to confirm the 4-phase init converged without committing to a full time loop.
```

Update the CLI:
```diff
  Diagnostics:
-   --dry-run                  Parse + construct + validate equilibrium; do not time-step
+   --dry-run                  Parse + construct + validate inputs + print derived
+                              quantities; do NOT solve domain or time-step.  Cheap.
+   --equilibrium-only         Run --dry-run path + domain Solve + equilibrium ψ
+                              init + verification; exit before time loop.  Expensive
+                              (one full domain Solve).
    --verify-dispatch          Print per-rank dispatch tags + exit before time loop
    --print-derived            After init, print L_nuc, T_nuc, V_nuc, dt_initial per region
```

**Test case:**
```cpp
void test_R008_dry_run_skips_domain_solve() {
   // Time the --dry-run path; assert it stays under 30 s on a 1000m mesh
   // (versus minutes for the full Solve).
   auto start = std::chrono::steady_clock::now();
   int rc = system("./seas_spatial_qd_driver --config small_safs.toml --dry-run");
   auto elapsed = std::chrono::steady_clock::now() - start;
   TEST_ASSERT_EQ(rc, 0);
   TEST_ASSERT(std::chrono::duration_cast<std::chrono::seconds>(elapsed).count() < 30);
   // Confirm no MUMPS log markers appear in the captured stdout.
   TEST_ASSERT(!stdout_contains("MUMPS"));
}
```

---

### [R-009] [MODERATE] [Phase 5 + Risk Assessment §"Low-confidence"] — `SEASQuasiDynamicOperator` per-DOF η: plan defers but never resolves; affects whether `solver/seas_operator.hpp` must be touched (CLAUDE.md "do not modify")

**Category:** ASSUMPTION

**Description:**
Risk Assessment §"Low-confidence" item 1 (`spatial_quasi_dynamic_plan.md:723-724`):

> Whether the existing `SEASQuasiDynamicOperator<ParMesh, ElasticityDomainOperator<ParMesh>, BP5FaultOp>` template needs any change to support per-DOF η for the radiation-damping term.  Currently `BP5FaultOp` uses scalar η.  If the per-DOF wiring requires touching `solver/seas_operator.hpp` (CLAUDE.md "do not modify"), that's a plan amendment.

CLAUDE.md (read for this review) explicitly lists `solver/seas_operator.hpp` as a "do not modify" file. Phase 5 of this plan adds a `SetPerDOFRateStateParams(a, b, Dc, V_init, f_0, V_0, eta)` setter on `RateStateFaultOperator` that includes per-DOF η. Inside `RateStateFaultOperator::ComputeRHS`, η is consumed in the radiation-damping term (`τ_radiation = η · V`). The plan assumes the per-DOF η is read from the operator's own member rather than from `SEASQuasiDynamicOperator`.

This needs to be RESOLVED in the plan, not deferred to "investigation". Concretely:

- If `RateStateFaultOperator::ComputeRHS` reads η from its own per-DOF state (likely, since SAFS-mode tau_pre/sigma_n already does so), no change to `solver/seas_operator.hpp` is needed — but the plan must say so explicitly.
- If `SEASQuasiDynamicOperator::Mult` or any coupling glue reads η scalar-style (e.g., as `fault_op.GetEta()`), the per-DOF wiring requires either touching `solver/seas_operator.hpp` (forbidden by CLAUDE.md) or adding a per-DOF accessor on `RateStateFaultOperator` that returns the per-DOF η in a place where the coupling code can pick it up.

Verified from `fault/rate_state_fault.hpp:534` and surrounding `ComputeRHS` lines that the operator already routes through `(*sigma_n_per_dof_)(i)` per-DOF for SAFS mode; the same pattern can extend to per-DOF η without touching `solver/seas_operator.hpp` — but Phase 5 must state this contract explicitly.

**Trigger:**
Phase 5 implementer reaches `ComputeRHS` and discovers the η read site.

**Actual behavior:**
Implementer hits a CLAUDE.md gate and stalls.

**Expected behavior:**
Phase 5 documents: "Per-DOF η is read INSIDE `RateStateFaultOperator::ComputeRHS` at the existing scalar-η read site; the change is `params_.eta` → `(per_dof_params_enabled_ ? per_dof_eta_(i) : params_.eta)`, all confined to `fault/rate_state_fault.hpp`.  `solver/seas_operator.hpp` is NOT touched."

**Suggested fix:**

Replace Risk Assessment §"Low-confidence" item 1 with a resolution paragraph in Phase 5:

In Phase 5 §"Detailed Requirements" add a new item 4:
```diff
+ 4. **Per-DOF η read site** (resolves Risk Assessment §"Low-confidence" item 1).  Verified 2026-05-18: `RateStateFaultOperator::ComputeRHS` (fault/rate_state_fault.hpp:457+) reads η as `params_.eta`.  The Phase-5 routing change is:
+   ```cpp
+   const real_t eta_i = per_dof_params_enabled_ ? per_dof_eta_(i) : params_.eta;
+   ```
+   `solver/seas_operator.hpp` is NOT touched — it consumes traction and slip rate from `RateStateFaultOperator` without ever inspecting η directly.  An audit grep `rg "\beta\b" miniapps/seas/solver/seas_operator.hpp` returns zero hits.
```

Drop the corresponding bullet from Risk Assessment §"Low-confidence" (line 723).

**Test case:**
```cpp
void test_R009_per_dof_eta_routing() {
   // Build a RateStateFaultOperator; install per-DOF eta with a 10% perturbation on
   // one DOF; assert ComputeRHS at that DOF differs by the expected sensitivity
   // and matches the scalar-eta path elsewhere.
   ...
}
void test_R009_seas_operator_does_not_read_eta() {
   // Compile-grep proxy: search the seas_operator.hpp string for "eta".
   std::ifstream f("miniapps/seas/solver/seas_operator.hpp");
   std::string contents((std::istreambuf_iterator<char>(f)),
                         std::istreambuf_iterator<char>());
   TEST_ASSERT(contents.find("eta") == std::string::npos);
}
```

---

### [R-010] [MODERATE] [Phase 0 + Friction directory] — Empty `friction/rate-and-state/` directory and undetermined friction-parameter coverage; plan's defaults are not motivated against the user's stated goal "VW + VS"

**Category:** DEVIATION

**Description:**
The user's request specified: *"friction: (undetermined rate-and-state friction parameters, should expect spatially varying to cover both velocity weakening and velocity strengthening)"*.

Verified on disk:
- `friction/rate-and-state/` is **empty** (no PDFs, no markdown, no TOML).
- `friction/slip-weakening/` has `geoffrey2010.md` + `131.pdf` (used by the dynamic-rupture sibling plan, not this one).
- `friction/barbot26c.pdf` and `friction/tse+rice86.pdf` exist at the top level — likely the intended rate-and-state references.

The plan's Phase 0 example TOML (`spatial_quasi_dynamic_plan.md:131-175`) uses:
- Shallow strengthening (depth-0 to -5 km): `a = 0.0125, b = 0.013` → `a > b` → VS ✓
- Seismogenic VW core (depth -5 to -15 km): `a = 0.005, b = 0.015` → `a < b` → VW ✓
- Nucleation patch: `V_init = 1e-2` ✓

So the layering DOES cover both VW and VS. But the values are pulled from BP5 (or similar canonical examples), not anchored to the SAFS reference papers (`barbot26c.pdf`, `tse+rice86.pdf`) the user is presumably citing. The user said "undetermined" — meaning the plan should make explicit (a) which reference the values come from, (b) what knobs the user is expected to vary before launching production runs, and (c) provide a `friction/rate-and-state/EXAMPLE_*.toml` that doesn't pretend to be a final value set.

Also: the validation rule "0 < a < b < 0.1" (Phase 0 §"Validation rules" item 3) would REJECT the shallow VS layer (a = 0.0125, b = 0.013 — passes since a < b; but if the user wants any DOF with a > b for stronger VS regions, the validator must change).

The validator constraint `a < b` is **wrong for the VS layer**: in a VS layer, the desired regime is `a > b` (positive feedback damped). The plan's `[[friction.rate_state.spatial]]` Block 1 shows `a = 0.0125, b = 0.013` (a < b) — which is VW-leaning. To truly cover the velocity-strengthening regime, some DOFs MUST have `a > b`. The Phase 0 validator forbids this and would abort.

**Trigger:**
User supplies a TOML with any spatial rule that produces `a > b` at any DOF (intended for shallow VS).

**Actual behavior:**
Phase 1 parser aborts with "validation failed: a < b required".

**Expected behavior:**
The rate-state validator must accept BOTH `a < b` (VW) and `a > b` (VS) — both are valid rate-and-state regimes. The constraint `a < b` is only required for VW seismogenic patches; the user explicitly asked for coverage of both regimes.

**Suggested fix:**

Replace Phase 0 §"Validation rules" item 3:
```diff
- - In rate-state law: `0 < a < b < 0.1`, `Dc > 0`, `V_0 > 0`, `sigma_n > 0`, `0 < f_0 < 1` at every DOF after spatial application.
+ - In rate-state law: `0 < a < 0.1`, `0 < b < 0.1`, `Dc > 0`, `V_0 > 0`, `sigma_n_eff > 0`, `0 < f_0 < 1` at every DOF after spatial application.  **The validator does NOT require `a < b`:** both rate-weakening (a < b, seismogenic core) and rate-strengthening (a > b, shallow/deep aseismic) regions are permitted — they coexist on a real fault.  A WARNING (not an abort) is emitted if no DOF satisfies `a < b` (the run cannot produce earthquakes) and a second WARNING if no DOF satisfies `a > b` (no stable creep regions — physically unusual but not invalid).
```

Add to Phase 1 acceptance criteria:
```diff
+ - [ ] `test_spatial_friction_resolver`: at least one test exercises a MIXED a/b layout — some DOFs with `a < b`, others with `a > b` — and verifies the resolver does NOT abort and the validator emits zero abort messages.
```

Add new example file (Phase 0 §"Files to Create"):
```diff
- - `safs/project_7.0_alternative/friction/EXAMPLE_spatial_friction_rate_state_safs.toml` — concrete example file (Aging-law DRS, BP5-style scalar fallback + SAFS-style per-region overrides) for the user to copy-edit.
+ - `safs/project_7.0_alternative/friction/rate-and-state/EXAMPLE_spatial_friction_rate_state_safs.toml` — example with depth-layered (a < b) seismogenic core sandwiched between (a > b) shallow- and deep-strengthening layers (per `barbot26c.pdf` Table 2 / Fig. 5), plus a nucleation patch `V_init`-driven kick.  The file's header MUST cite the source paper for every default value, and end with a "USER MUST REVIEW" banner reminding the user the friction parameters are starting points, not final.
+ - `safs/project_7.0_alternative/friction/rate-and-state/README.md` — one-page summary of which parameters are knobs vs which are pinned, and the expected effect of each on recurrence interval / event size.
```

**Test case:**
```cpp
void test_R010_validator_accepts_mixed_vw_vs() {
   // TOML with depth-0..-5km: a=0.020, b=0.010 (VS) and -5..-15km: a=0.005, b=0.015 (VW).
   // Loader must succeed; resolver must produce a, b vectors with both a<b and a>b DOFs.
   auto cfg = LoadSpatialFrictionConfig("fixtures/mixed_vw_vs.toml");
   auto p = resolver.ResolveRateState(*cfg.rate_state, dof_coords, ...);
   bool seen_vw = false, seen_vs = false;
   for (int i = 0; i < p.a.Size(); ++i) {
      if (p.a(i) <  p.b(i)) { seen_vw = true; }
      if (p.a(i) >  p.b(i)) { seen_vs = true; }
   }
   TEST_ASSERT(seen_vw && seen_vs);
}
```

---

### [R-011] [MODERATE] [Phase 4 step 5 / overall] — `safs_test_driver` is referenced as a BoundaryConfig example, but no such driver exists in the repo

**Category:** ASSUMPTION

**Description:**
Phase 4 step 5 (`spatial_quasi_dynamic_plan.md:559`):

> Reuse existing `BoundaryConfig` from `safs_test_driver` (fault_attr = 101, dirichlet_attrs = sides + bottom, natural_attr = top — per the `safs_fault_box_nwcut.geo` Physical Surface tags above).

Verified: `find miniapps/seas -name "safs_test_driver*"` returns **no results**. `grep -rln "safs_test_driver"` returns no source-file hits. The reference is fabricated.

This compounds R-001 (the actual BoundaryConfig wiring is non-trivial: per-attribute Dirichlet funcs for the four sides) — the implementer cannot copy from an example that doesn't exist.

**Trigger:**
Phase 4 implementer reaches step 5.

**Actual behavior:**
Implementer cannot locate the reference driver.

**Expected behavior:**
Plan provides the full BoundaryConfig construction inline (per R-001 fix), or cites an actually-existing driver that uses the same Physical Surface layout.

**Suggested fix:**

Combined with R-001: replace Phase 4 step 5 with explicit `BoundaryConfig` construction code (see R-001 §Suggested fix). Drop the `safs_test_driver` reference entirely.

**Test case:** (subsumed by R-001 test)

---

### [R-012] [MODERATE] [Phase 4 step 11 + Phase 5] — Phase 4 step 11 names `SetPerDOFRateStateParams` as the wiring point but never shows the call; Phase 5 setter requires SAFS-mode-on as a precondition, creating an implicit ordering trap

**Category:** EDGE_CASE

**Description:**
Phase 4 step 11 (`spatial_quasi_dynamic_plan.md:571`) says "The Phase-1 per-DOF `a, b, Dc, V_init, f_0, V_0, eta` are pushed into `RateStateFaultOperator` via new setters added in Phase 5." but never shows the call site.

Phase 5 (`spatial_quasi_dynamic_plan.md:631-645`) says the setter:

> Must be called BEFORE Init. ... Aborts if SAFS mode is not enabled (per-DOF parameters only make sense alongside per-DOF tau_pre / sigma_n).

So the ordering constraint is:
1. Construct `RateStateFaultOperator` (which caches BP5 scalar params; tau_pre_ defaults to whatever `geom.GetTauPre()` returns at that moment).
2. Call `SetSAFSMode(true, &tau_pre_per_dof, &sigma_n_per_dof)`.
3. Call `SetPerDOFRateStateParams(a, b, Dc, V_init, f_0, V_0, eta)`.
4. Call `PreInit(state)`.
5. (Domain Solve)
6. Call `Init(state)` — uses per-DOF V_init for ψ steady-state computation.

If any of these is out of order, the operator silently uses partial-state info. Specifically:
- If `SetPerDOFRateStateParams` is called BEFORE `SetSAFSMode(true)`, the new setter aborts (per the Phase 5 contract). OK.
- If `SetPerDOFRateStateParams` is called AFTER `PreInit`, the per-DOF V_init is set but `PreInit` already wrote ψ = steady-state(V_init_scalar) into `state`. The per-DOF V_init is then read in `Init` (step 6) — inconsistent with what `PreInit` already wrote. Plan must call this out.

Worse: the `RateStateFaultOperator` BP5 ctor caches `tau_pre_ = geom_->GetTauPre()` at construction (verified line 170). If `ApplySpatialStressSidecar` is called BEFORE constructing the operator (as the plan correctly says), `geom.GetTauPre()` returns the SAFS-projected tau_pre. Good. But the operator ALSO caches `bp5_params_` from the seed, and the seed is the placeholder `BP5Params{}` per step 8. So the operator's internal scalar BP5 params are garbage until `SetPerDOFRateStateParams` overwrites them — and any read of `params_.a`, `params_.b`, etc. between operator construction and the setter call sees garbage.

In particular: `PreInit(state)` may read `params_.V_init` to set initial slip (verified line 287-292 region). With the scalar seed at `V_init = 0`, the initial state vector is all zeros — fine. But then when `SetPerDOFRateStateParams` is called after `PreInit`, the V_init scalar is replaced but the ALREADY-WRITTEN state vector is not re-initialized.

**Trigger:**
Driver calls `PreInit` before `SetPerDOFRateStateParams`.

**Actual behavior:**
State vector is initialized from the seed scalar V_init (0 with default `BP5Params{}`), then per-DOF V_init replaces the scalar but the state vector keeps the original initial slip.

**Expected behavior:**
The setter MUST be called BEFORE `PreInit` (and Phase 5 should add an explicit `MFEM_VERIFY` that `PreInit` has not yet run when the setter fires).

**Suggested fix:**

Update Phase 4 step 11 (be explicit about the order):
```diff
- 11. **Construct `RateStateFaultOperator<ParMesh, 2>`** using the BP5 constructor that accepts a `FaultGeometry` and `BP5Params`.  **Then immediately call `fault_op.SetSAFSMode(true, &tau_pre_per_dof, &sigma_n_per_dof)`** — where `tau_pre_per_dof` and `sigma_n_per_dof` are the vectors populated by `ComputeSAFSParams`.  ... The Phase-1 per-DOF `a, b, Dc, V_init, f_0, V_0, eta` are pushed into `RateStateFaultOperator` via new setters added in Phase 5.
+ 11. **Construct `RateStateFaultOperator<ParMesh, 2>`** using the BP5 constructor that accepts a `FaultGeometry` and the seed `BP5Params` from step 8.  Wire the per-DOF data in the following EXACT order (each call is required; reordering causes silent state-vector corruption):
+   ```cpp
+   RateStateFaultOperator<ParMesh, 2> fault_op(geom, bp5_params_seed, mpi.get());
+   // 11a. SAFS mode: route tau_pre and sigma_n through per-DOF vectors.
+   fault_op.SetSAFSMode(true,
+                        &geom.GetTauPre(),            // 2 * num_owned_fault_dofs
+                        &geom.sigma_n_per_dof());     // num_owned_fault_dofs
+   // 11b. SAFS friction params: route a/b/Dc/V_init/f_0/V_0/eta through per-DOF
+   //      vectors.  MUST be called AFTER SetSAFSMode and BEFORE PreInit.
+   fault_op.SetPerDOFRateStateParams(rs.a, rs.b, rs.Dc, rs.V_init,
+                                       rs.f_0, rs.V_0, rs.eta);
+   ```
+   Phase 5 enforces the ordering via `MFEM_VERIFY(!preinit_called_, ...)` and `MFEM_VERIFY(safs_mode_, ...)` inside `SetPerDOFRateStateParams`.  The driver is responsible for getting the order right; the verifier surfaces violations at the right place.
```

Update Phase 5 §"Detailed Requirements" item 2 to add the ordering guard:
```diff
- 2. **`SetPerDOFRateStateParams` semantics**:
-    - All seven vectors must have size `num_nodes_` (MFEM_VERIFY each).
+ 2. **`SetPerDOFRateStateParams` semantics**:
+    - **Ordering guard (new)**: `MFEM_VERIFY(safs_mode_, "SetPerDOFRateStateParams: SAFS mode must be enabled first (call SetSAFSMode(true, ...) before this).")` and `MFEM_VERIFY(!preinit_called_, "SetPerDOFRateStateParams: must be called BEFORE PreInit.")`.  Phase 5 adds a `preinit_called_` bool member that `PreInit` sets to `true`.
+    - All seven vectors must have size `num_nodes_` (MFEM_VERIFY each).
```

**Test case:**
```cpp
void test_R012_per_dof_setter_requires_safs_mode() {
   FaultGeometry<ParMesh> g(dom, seed, &mpi);
   RateStateFaultOperator<ParMesh, 2> op(g, seed, &mpi);
   Vector a(g.NumFaultDOFs()); a = 0.01;
   // ... fill the other six vectors ...
   EXPECT_ABORT_WITH(op.SetPerDOFRateStateParams(a, b, Dc, Vi, f0, V0, eta),
                     "SAFS mode must be enabled first");
}
void test_R012_per_dof_setter_must_precede_preinit() {
   FaultGeometry<ParMesh> g(dom, seed, &mpi);
   RateStateFaultOperator<ParMesh, 2> op(g, seed, &mpi);
   op.SetSAFSMode(true, &g.GetTauPre(), &g.sigma_n_per_dof());
   Vector state(op.StateSize());
   op.PreInit(state);
   Vector a(g.NumFaultDOFs()); a = 0.01;
   // ... fill the other six vectors ...
   EXPECT_ABORT_WITH(op.SetPerDOFRateStateParams(a, b, Dc, Vi, f0, V0, eta),
                     "must be called BEFORE PreInit");
}
```

---

### [R-013] [MODERATE] [Phase 0 §pore-pressure + Phase 4] — Pore-pressure schema is `P_p_pa + grad * max(0, -z)`; the SAFS mesh has `z = 0` at the free surface AND `zmin_fault` may be `< 0`, but no test verifies the depth-coordinate convention matches between mesh and TOML

**Category:** ASSUMPTION

**Description:**
Phase 0 schema (`spatial_quasi_dynamic_plan.md:107-109`):

```toml
[pore_pressure]
P_p_pa            = 0.0
P_p_grad_pa_per_m = 0.0
min_sigma_n_pa    = 0.0
```

Phase 1 algorithm (`spatial_quasi_dynamic_plan.md:343`):

```
sigma_n_eff_i = max(min_sigma_n_pa,
                    sigma_n_i - (P_p_pa + P_p_grad_pa_per_m * max(0, -z_i)))
```

This assumes `z = 0` at the free surface, `z < 0` going down (CLAUDE.md confirms this is the project convention). Verified `meshing/code/safs_fault_box_nwcut.geo:78-79`:

```
zmin = zmin_fault - pad_bottom;
zmax = zmax_fault + pad_top;
```

with `PAD_TOP_DEFAULT = 0.0` so the mesh `zmax = fault_zmax`, where `fault_zmax = 0` for the "cleaned STL". The pore-pressure depth term `max(0, -z_i)` is non-negative only for `z_i < 0` — i.e., below the free surface — so this is fine for fault DOFs with `z < 0`.

**But:** any fault DOF that lies AT `z = 0` (top edge of the fault sticking out to the free surface) gets `max(0, -0) = 0`, so pore-pressure depth term is 0 at the surface — physically correct. The schema also accepts `P_p_pa > 0` (constant offset) and would push sigma_n_eff negative at the surface if `sigma_n_total - P_p_pa < 0`. The `min_sigma_n_pa` clamp catches this, but only if the user sets it > 0; default 0 means no clamp. Plan does not document the failure mode when sigma_n_eff goes negative (Brent solver bracket violation, NaN in `psi/a` → silent garbage).

Compounding: Phase 1 validation §item 3 says "sigma_n > 0" at every DOF after spatial application — verified against the seed scalar `sigma_n_default` only, NOT against the pore-pressure-adjusted `sigma_n_eff`. The actual rate-state operator consumes `sigma_n_eff`, not the raw `sigma_n_default`. Plan needs to validate the EFFECTIVE value, not just the input.

**Trigger:**
TOML with `P_p_pa = 60e6` (typical hydrostatic + lithostatic estimate at 5 km) and `sigma_n_default = 50e6` (Phase 0 example) → at shallow DOFs, `sigma_n_eff = 50e6 - 60e6 = -10e6 < 0`.

**Actual behavior:**
With default `min_sigma_n_pa = 0`, the eff value is negative; Brent solver sees a meaningless friction problem and either NaNs out or returns a wrong V.

**Expected behavior:**
Validation at the END of resolution must check `sigma_n_eff_i > 0` per DOF (after pore-pressure subtraction and `min_sigma_n_pa` clamp). Plan must also document a recommended `min_sigma_n_pa > 0` floor (e.g., 1e6 Pa = 1 MPa) and require it for any TOML that sets `P_p_pa > 0` or `P_p_grad > 0`.

**Suggested fix:**

In Phase 1 §"Detailed Requirements" item 2 (`ResolveRateState` algorithm), step 5 already validates "sigma_n_eff_i > 0" — keep this. Add to Phase 0 §"Validation rules":

```diff
+ - **Pore-pressure consistency check**: if `P_p_pa > 0` OR `P_p_grad_pa_per_m > 0`, then `min_sigma_n_pa MUST be > 0` (recommend ≥ 1e6 Pa).  Validator aborts with the offending values if this is violated and any DOF would have `sigma_n_eff ≤ 0`.
```

Add to Phase 1 acceptance criteria:
```diff
+ - [ ] `test_spatial_friction_resolver_sigma_n_eff_clamp`: TOML with `P_p_pa = 60e6, sigma_n_default = 50e6, min_sigma_n_pa = 1e6` produces sigma_n_eff = 1e6 at shallow DOFs (clamp fired), sigma_n_eff > 1e6 at deep DOFs (grad pushes it up).
```

**Test case:**
```cpp
void test_R013_sigma_n_eff_clamp() {
   PorePressureSpec pp;
   pp.P_p_pa = 60e6;
   pp.P_p_grad_pa_per_m = 0.0;
   pp.min_sigma_n_pa = 1e6;
   RateStateBlock cfg; cfg.sigma_n_default = 50e6;
   /* ... */
   auto p = resolver.ResolveRateState(cfg, dof_coords, dof_elem, dof_attr, material, pmesh);
   for (int i = 0; i < p.sigma_n_eff.Size(); ++i) {
      TEST_ASSERT(p.sigma_n_eff(i) >= 1e6);
   }
}
void test_R013_sigma_n_eff_negative_aborts() {
   PorePressureSpec pp;
   pp.P_p_pa = 60e6;
   pp.min_sigma_n_pa = 0.0;
   RateStateBlock cfg; cfg.sigma_n_default = 50e6;
   /* ... */
   EXPECT_ABORT(resolver.ResolveRateState(cfg, ...));
}
```

---

### [R-014] [LOW] [Phase 5 §"Files to Modify"] — Phase 5 claim "No `.cpp` files for `rate_state_fault.hpp` (it's header-only)" is correct; but the additive header changes will increase compile time of every translation unit that includes it

**Category:** QUALITY

**Description:**
`fault/rate_state_fault.hpp` is header-only (verified: no companion .cpp). Phase 5 adds seven new per-DOF vectors as private members + a boolean `per_dof_params_enabled_` + a `SetPerDOFRateStateParams` method. This is fine for correctness; however, every translation unit that `#include`s the header will now compile a larger struct.

This is LOW because the impact is purely build-time; flagging in case the reviewer wants to split the implementation into a `.cpp`.

**Suggested fix:**
None required. Optionally consider moving the new setter body into a `rate_state_fault.cpp` if compile time becomes an issue — but only after measuring.

**No test required** — this is a build-quality observation.

---

## Summary
- Critical issues: 8 (R-001..R-008)
- Moderate issues: 5 (R-009..R-013)
- Low issues: 1 (R-014)
- Plan compliance: **INCOMPLETE** — the plan references multiple fabricated APIs (`PETScTSCheckpoint::Restore`, `safs_test_driver` example, `MaterialField`-aware `ElasticityDomainOperator` ctor that does not exist on current main, mesh-tagged velocity sidecar path), and the boundary-loading scheme is incompatible with the only mesh in scope (all four sides collapsed into one attribute).
- **Verdict: FAIL — must fix R-001 through R-008 before proceeding to implementation.**

## Unreviewed Areas
- **Phase 3 (stress sidecar bundle)** — surface-area is small and the underlying `StressField3D` / `ComputeSAFSParams` are already exercised by `seas_test_compute_safs_params` (13/13 passing per plan note). Spot-checked the bundle interface; no concrete bugs found, but the pre-flight `NumZeroNormalFallbacks == 0` guard depends on `ComputePerDOFCoordsAndBasis_` correctly detecting the SAFS curvilinear-fault geometry — not verified end-to-end on the actual 1000m_lcfar3000 mesh (will fire at runtime if there's a problem).
- **Phase 6 (sbatches + verification scripts)** — accepted as boilerplate; no review pass.
- **Numerical-constraints §"Quasi-dynamic tolerances"** — values copied from BP5; not independently verified against SAFS dynamics. Leaving for an empirical pass after the first smoke run.
- **Risk Assessment §"Long-cycle restart determinism"** — accepted as a known unknown.

## Open Questions for the User (to inform the /code-plan revision)

OQ-1. Should the QD plan adopt the dynamic-rupture plan's "in-house Phase H" pattern and OWN the `ElasticityDomainOperator(MaterialField)` ctor work (delivering it as a new "Phase E"), or block on `heterogeneous_material_plan.md` Phase 2 merging first? Recommendation: own it (mirror dynamic-rupture rev-3 D-2) — reduces cross-plan dependency.

OQ-2. Per R-001, the mesh needs regeneration with four separate side-face attributes. Do you want the plan to include a "Phase −1: Mesh regeneration" subtask, or is this a precondition you'll handle outside the plan?

OQ-3. For the rate-state friction defaults (R-010), should the plan cite `barbot26c.pdf` and `tse+rice86.pdf` as the authority for the spatial layout, or do you have a SAFS-specific reference you'd prefer?

OQ-4. Should the `--dry-run` / `--equilibrium-only` split (R-008) be exposed as two separate flags, or fold the equilibrium step under `--print-derived`? Recommendation: two flags (cleaner semantics; the equilibrium solve is a non-trivial cost).
