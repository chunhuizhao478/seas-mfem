# Code Review: Phase 2b Stages 1b / 2 / 3 — 2026-06-04

## Review Scope
- Plan: `document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md` §Phase 2b (heterogeneous material), Stages 1b (IP method-gate), 2 (driver wiring), 3 (analytic slab test).
- Files reviewed:
  - `domain/elasticity_operator.hpp` — heterogeneous-ctor IP method-gate + ctor doc.
  - `drivers/spatial_seas_driver.cpp` — §4.0/4.4/4.7 material build + ctor selection.
  - `Makefile` — QD driver material-object link set; the two new test targets.
  - `config/safs_qd/bp5_phase2b_depthprofile_smoke.toml` — heterogeneous smoke config.
  - `tests/unit/test_elasticity_heterogeneous_ctor.cpp`, `tests/unit/test_elasticity_heterogeneous_slab.cpp`.
- Domain context: CLAUDE.md (no silent fallback; extreme-care domain/), the QD plan, the archived heterogeneous_material_plan.md (R-002), `spatial_dyn_driver.cpp` material-build reference, `spatial_friction.hpp` (VelocitySpec/MaterialSpec defaults), `linear_elastic.hpp`.
- Method: three adversarial passes; build + run exercised (het_ctor 8/8, slab 27/27, constant + depth-profile driver dry-runs np1/np4, config 40/40); **R-501 empirically reproduced** (constant config without `use_sidecar=false` hard-aborts).

## Findings

### [R-501] MODERATE [DEVIATION] drivers/spatial_seas_driver.cpp:§4.4 — sidecar gating ignores `[material].kind="constant"`; a constant config without `use_sidecar=false` attempts a sidecar and hard-aborts

**Category:** DEVIATION (plan: "when Constant, keep the LinearElastic ctor") / EDGE_CASE

**Description:**
The Stage-2 gate is `sidecar_requested = cfg.velocity.use_sidecar && cfg.material.kind != DepthProfile1D`. But `VelocitySpec::use_sidecar` **defaults to `true`** (`spatial_friction.hpp:246`) and `MaterialSpec::kind` defaults to `Constant`. So a config that sets `[material].kind="constant"` (or omits `[material]`) but does NOT set `[velocity].use_sidecar=false` has `sidecar_requested == true` → the driver tries to load a velocity sidecar instead of using the constant material. The explicit `kind=Constant` is ignored. Worse: the failure is a hard `MFEM_ABORT` from inside `LoadSpatialVelocityBundle`'s spec validation ("[velocity] must set 'dataset_root' … when 'use_sidecar' is true"), which is **not** a `std::exception` and therefore **bypasses the `try/catch`** — so the code's own guiding message ("Set [velocity].use_sidecar=false …") and the `return 4` path never execute; the process `MPI_ABORT`s (exit 1). Pre-Stage-2 such a config aborted at the §4.0 gate with a clear "set use_sidecar=false" message; Stage 2 replaced that with a cryptic sidecar-spec abort.

**Trigger (reproduced):** `[material]` constant (default) + `[velocity]` present without `use_sidecar=false` (so `use_sidecar` defaults true) → `./seas_spatial_seas_driver --config … --dry-run` → `MPI_ABORT` exit 1, message "[velocity] must set either 'dataset_root' … when 'use_sidecar' is true (the default)".

**Actual behavior:** constant config attempts a sidecar load and hard-aborts.

**Expected behavior:** `kind=Constant` → constant material (the `LinearElastic` ctor), regardless of the `use_sidecar` default.

**Suggested fix:** honor an explicit `kind=Constant` in the gate (a sidecar then requires a non-Constant kind, matching the explicit user intent):
```diff
       const bool sidecar_requested =
          cfg.velocity.use_sidecar
-         && cfg.material.kind != spatial::MaterialKind::DepthProfile1D;
+         && cfg.material.kind != spatial::MaterialKind::DepthProfile1D
+         && cfg.material.kind != spatial::MaterialKind::Constant;
```
(The constant smoke config still works; the depth-profile config still works; a sidecar now requires `kind="sidecar_hdf5"` — or any non-Constant kind — with `use_sidecar=true`.)

**Test case:**
```cpp
// test_R501_constant_kind_uses_constant_material:
//   write a TOML with [material].kind="constant" (or omit [material]) and a
//   [velocity] block WITHOUT use_sidecar=false; run --dry-run; EXPECT exit 0
//   and "[spatial_seas] material: constant" (currently: MPI_ABORT exit 1).
```

---

### [R-502] LOW [POSSIBLE] drivers/spatial_seas_driver.cpp:§4.4 — `kind="sidecar_hdf5"` with `use_sidecar=false` silently falls back to constant

**Category:** EDGE_CASE / ASSUMPTION

**Description:**
`MaterialKind::SidecarHDF5` is never branched on directly; the sidecar is reached only via `use_sidecar`. So `[material].kind="sidecar_hdf5"` + `[velocity].use_sidecar=false` (a contradictory config) silently uses the constant fallback with no warning — the opposite silent-fallback direction from R-501. CLAUDE.md forbids silent fallback. (After the R-501 fix this remains: kind=SidecarHDF5 + use_sidecar=false → constant.)

**Suggested fix:** detect the contradiction and abort (or warn loudly):
```diff
+      MFEM_VERIFY(!(cfg.material.kind == spatial::MaterialKind::SidecarHDF5
+                    && !cfg.velocity.use_sidecar),
+                  "spatial_seas_driver: [material].kind=\"sidecar_hdf5\" "
+                  "requires [velocity].use_sidecar=true.");
```

**Test case:** config with `kind="sidecar_hdf5"` + `use_sidecar=false` → expect a clear abort, not a silent constant run.

---

### [R-503] LOW [QUALITY] tests/unit/test_elasticity_heterogeneous_{ctor,slab}.cpp — `CreateTestMesh3D`/`AddFaultBoundaryElements` triplicated

**Category:** QUALITY (duplication)

**Description:**
The small-fault-mesh helpers `AddFaultBoundaryElements` + `CreateTestMesh3D` are now copy-pasted in three files (`test_elasticity_operator.cpp`, `test_elasticity_heterogeneous_ctor.cpp`, `test_elasticity_heterogeneous_slab.cpp`). A fix to the fixture (e.g., a BC-tag convention change) must be replicated 3×, and the copies can silently drift.

**Suggested fix:** extract to a shared header, e.g. `tests/unit/fault_mesh_fixture.hpp`, and `#include` it from all three. (Mechanical; no behavior change.)

---

### [R-504] LOW [POSSIBLE] tests/unit/test_elasticity_heterogeneous_ctor.cpp:main — bit-for-bit Vector compare across two operators assumes identical DOF ordering

**Category:** ASSUMPTION

**Description:**
`MaxAbsDiff(rhs_slip_c, rhs_slip_h)` / `u_c` vs `u_h` / `trac_c` vs `trac_h` compare Vectors from two SEPARATE operators built on two SEPARATE `Mesh` objects, index-by-index. This is correct ONLY because `CreateTestMesh3D(2,2,2,...)` is deterministic and produces identical topology+numbering for both. The assumption is currently sound but implicit; if mesh construction ever became ordering-nondeterministic (e.g., a future reordering pass), the test would report spurious failures.

**Suggested fix:** add an explicit comment documenting the identical-construction assumption, and assert structural equality first:
```cpp
Check(op_const.GetFESpace().GetVSize() == op_het.GetFESpace().GetVSize(),
      "both operators have identical FE-space size (identical meshes)");
```
(Already partially covered by the `u_c.Size()==u_h.Size()` check; make the intent explicit.)

---

### [R-505] LOW [POSSIBLE] domain/elasticity_operator.hpp:ElasticityDomainOperator(Coefficient&,...) — heterogeneous ctor defaults `method=DGMethod::BR2`, which always aborts

**Category:** QUALITY / EDGE_CASE

**Description:**
The heterogeneous ctor signature defaults `DGMethod method = DGMethod::BR2`, but the in-body guard aborts for any method other than IP. So a caller that uses the heterogeneous ctor WITHOUT explicitly passing `DGMethod::IP` (relying on the default) hits the abort. IP is currently the only supported heterogeneous method, so defaulting to it would make the common case work without the explicit argument and avoid a surprising abort-on-default.

**Suggested fix:**
```diff
       ElasticityDomainOperator(MeshType &mesh, int order,
                                Coefficient &lambda_c, Coefficient &mu_c,
                                real_t Vp, real_t Wf, real_t lf,
                                const BoundaryConfig &bdr_config,
-                               DGMethod method = DGMethod::BR2,
+                               DGMethod method = DGMethod::IP,
                                SolverType solver_type = SolverType::MUMPS_BLR,
                                const DomainConfig &config = {})
```
(The QD driver and both tests already pass `DGMethod::IP` explicitly, so this is a safety/clarity change, not a behavior change for current callers.)

---

## Summary
- Critical issues: **0**.
- Moderate issues: **1** — R-501 (sidecar gating ignores `kind=Constant`; constant config without `use_sidecar=false` hard-aborts; empirically reproduced).
- Low issues: **4** — R-502 (SidecarHDF5 + use_sidecar=false silent constant), R-503 (test-helper triplication), R-504 (cross-operator Vector-compare assumption), R-505 (hetero ctor default method aborts).
- Plan compliance:
  - Stage 1b — **FULL**: IP method-gate (heterogeneous IP allowed, BR2 aborts); the het-ctor==constant-ctor bit-for-bit test (8/8) validates the per-qp coefficient plumbing through the integrators.
  - Stage 2 — **PARTIAL**: material build + ctor selection work for the depth-profile path (np1/np4 verified) and the constant path (verified), BUT the gating mishandles `kind=Constant` with the `use_sidecar=true` default (R-501).
  - Stage 3 — **FULL**: the slab test asserts `strike = μ(depth)·γ` exactly and the deep/shallow ratio `= μ₁/μ₀` (27/27); genuinely exercises per-qp heterogeneity.
- Verdict: **PASS WITH FIXES** — fix R-501 before Phase 8 (SAF) and before documenting the constant-config workflow (it changes the failure mode for the common minimal config). R-502/503/504/505 are hardening.

## Unreviewed Areas
- The velocity-sidecar **HDF5 load path** itself (`LoadSpatialVelocityBundle` → `DataField3D`/`StressField3D`) — wired and linked, but not exercised locally (needs a sidecar `.h5`); it reaches the same `Mode::Coefficient` → coefficient ctor as the verified depth-profile path. Recommend a Frontera/with-sidecar smoke before Phase 8.
- BR2 + heterogeneous — guarded abort (further work), not exercised.
- Phase 8 SAF curved-tet mesh behavior of the per-qp coefficient — out of scope (Phase 8).
```
