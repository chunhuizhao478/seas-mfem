# Code Review: 2026-05-19 — spatial_dyn_driver heterogeneous-Riemann wiring vs working TPV* drivers

**Branch context:** `feature/heterogeneous_riemann_solver`.  The entire
purpose of this branch is to route every TPV* / SAFS dynamic-rupture run
through the heterogeneous (`MaterialField`, `BoundaryConfig`)
`WaveOperator` ctor (`dynamic/wave_operator.inl`, Phase R.2 bi-material
exact linearised Riemann solver per Pelties et al. 2012 / SeisSol).  The
`(λ, μ, ρ, bc)` scalar ctor exists for byte-parity regression against
the native drivers ONLY.  Any TPV* TOML that re-enables the scalar ctor
defeats the branch's purpose and should be flagged.

## Review Scope

- Plan documents:
  - `miniapps/seas/safs/project_7.0_alternative/document/spatial_dynamic_rupture_plan.md` (rev-3).
  - `miniapps/seas/safs/project_7.0_alternative/document/PLAN_phase_R_exact_bimaterial_riemann_rev3.md`.
  - `miniapps/seas/safs/project_7.0_alternative/debug_document/spatial_paraview_compaction_parity_plan_2026-05-18.md`.
- Files reviewed:
  - `miniapps/seas/drivers/spatial_dyn_driver.cpp` (2373 lines).
  - `miniapps/seas/tpv102/configs/tpv102.toml`.
  - `miniapps/seas/tpv104/configs/tpv104.toml`.
  - `miniapps/seas/tpv205/configs/tpv205.toml`.
  - `miniapps/seas/tpv31/configs/tpv31.toml`.
- Reference (native drivers — known good):
  - `miniapps/seas/drivers/tpv102_driver.cpp`.
  - `miniapps/seas/drivers/tpv104_driver.cpp`.
  - `miniapps/seas/drivers/tpv205_driver.cpp`.
  - `miniapps/seas/config/tpv102_params.hpp`, `tpv104_params.hpp`, `tpv205_params.hpp`.
  - `miniapps/seas/dynamic/tpv102_setup.hpp`, `tpv205_setup.hpp`, `tpv205_friction.hpp`.
- Supporting code:
  - `miniapps/seas/spatial/code/spatial_friction.{hpp,cpp}` (parser + resolver).
  - `miniapps/seas/spatial/code/spatial_stress.{hpp,cpp}` (stress sources).
  - `miniapps/seas/dynamic/spatial_setup.hpp` (fault-DOF init).
  - `miniapps/seas/dynamic/spatial_nucleation.{hpp,cpp}` (nucleation resolvers).
  - `miniapps/seas/fault/fault_geometry_safs_templated.inl` (Cauchy → fault-local projection).
- Domain context: `miniapps/seas/CLAUDE.md` (Canonical Coordinate System, Bi-material Riemann (Phase R) section).

---

## Findings

### [R-001] [CRITICAL] [tpv205.toml / tpv31.toml] — `interior_flux = "scalar"` defeats the entire branch purpose (heterogeneous Riemann solver)

**Category:** DEVIATION

**Description:**
This branch (`feature/heterogeneous_riemann_solver`) exists to route every interior face through `BimaterialFlux::ApplyPerFaceFlux` via the new heterogeneous ctor `WaveOperator<MeshT>(mesh, order, MaterialField, BoundaryConfig)` introduced by Phase R.2 (`drivers/spatial_dyn_driver.cpp:941-965` dispatches between scalar and heterogeneous ctors; the scalar ctor exists for byte-parity regression ONLY).

The TOML parser defaults `cfg.numerics.interior_flux` to `"bimaterial"` (`spatial/code/spatial_friction.hpp:225`, parsed at `spatial_friction.cpp:1147`).  TPV102 and TPV104 do NOT set this key and therefore correctly route through the heterogeneous ctor — that is what this branch is for.

But `tpv205.toml:170` and `tpv31.toml:224` both explicitly set:

```toml
interior_flux = "scalar"    # scalar WaveOperator ctor — byte-identical interior flux
```

The TPV205 TOML's accompanying comments (`tpv205.toml:167-170`) frame this as a temporary "REVIEW R-004 / R-005 / R-006: opt in to native TPV205 dispatch" — a vestige from a prior round of byte-parity validation against `drivers/tpv205_driver.cpp` that has outlived its purpose on this branch.  TPV31 mirrors the comment.  Both TOMLs ship in a state where they:

1. Compile against the heterogeneous WaveOperator ctor.
2. Run the new spatial driver.
3. But execute the **same interior-face flux as the native scalar driver** — i.e., the heterogeneous Riemann path is never exercised.

A TPV205 / TPV31 production run with the canonical TOML, on this branch, therefore proves nothing about the bi-material Riemann solver and produces output byte-equivalent to the native scalar driver (within FP rounding).  The CLAUDE.md "Bi-material Riemann (Phase R)" section's claim that "Verified on TPV205 (homogeneous regression — relative parity ≤ 1.4e-12...)" was the byte-parity verification gate; with the gate passed, the TOMLs should now route through the heterogeneous ctor by default, which means **removing** the `interior_flux = "scalar"` override.

**Trigger:**
```
seas_spatial_dyn_driver --config tpv205/configs/tpv205.toml
seas_spatial_dyn_driver --config tpv31/configs/tpv31.toml
```

**Actual behavior:**
The driver's startup banner prints
```
[wave] interior_flux=scalar: using scalar WaveOperator(λ, μ, ρ) ctor for byte parity ...
```
and constructs the scalar ctor (`spatial_dyn_driver.cpp:958-963`).  The heterogeneous (`MaterialField`) ctor is NEVER reached.  `wave.GetGodunovFluxPool()` is null; `BimaterialFlux::ApplyPerFaceFlux` is never called.

**Expected behavior:**
The driver constructs the heterogeneous ctor (`std::make_unique<WaveOperator<ParMesh>>(pmesh, order, material, bc)`); the bi-material flux dispatch routes every interior face through `BimaterialFlux::ApplyPerFaceFlux`.  For TPV205 with a homogeneous `MaterialField::MakeConstant(λ, μ, ρ)`, this is mathematically equivalent to the scalar ctor (parity ≤ 1e-12), so the existing TPV205 trace-comparison tests against SCEC reference data continue to pass.

**Suggested fix:**
Remove the override from both TOMLs.  Optionally retain it as a debug knob for future byte-parity regressions, but it must NOT be the canonical configuration.

```diff
 [numerics]
 ader_order     = 2
 mixed_flux     = "none"      # REQUIRED for bi-material runs (R.2 R-002)
 cfl            = 0.5
 use_pml        = false
-# REVIEW R-004 / R-005 / R-006: opt in to native TPV205 dispatch.
 cfl_safety     = "dg"        # cfl /= (3*(2p+1)) — matches tpv205_driver
 fault_iterator = "one-shot"  # O = 1 sub-step — matches tpv205_driver default
-interior_flux  = "scalar"    # scalar WaveOperator ctor — byte-identical interior flux
+# interior_flux omitted → defaults to "bimaterial" (heterogeneous Riemann
+# solver, the purpose of branch feature/heterogeneous_riemann_solver).
+# Set explicitly to "scalar" only for byte-parity regression against the
+# native TPV205 driver.
```

(Apply identical edit to `tpv31.toml`.)

**Test case:**
```c++
// tests/integration/test_canonical_tpvs_use_heterogeneous_ctor.cpp
TEST_CASE("R-001: canonical TPV205/TPV31 TOMLs route through heterogeneous Riemann ctor") {
    for (const auto& cfg_path : { "tpv102/configs/tpv102.toml",
                                   "tpv104/configs/tpv104.toml",
                                   "tpv205/configs/tpv205.toml",
                                   "tpv31/configs/tpv31.toml" }) {
        auto cfg = spatial::LoadSpatialFrictionConfig(cfg_path);
        REQUIRE(cfg.numerics.interior_flux == "bimaterial");
    }
}

// And the runtime dispatch check (per-rank assertion fires inside ctor):
TEST_CASE("R-001: spatial_dyn instantiates BimaterialFlux on canonical TPV205") {
    auto wave = run_spatial_dyn_construction_only("tpv205/configs/tpv205.toml");
    REQUIRE(wave->GetGodunovFluxPool() != nullptr);     // heterogeneous ctor populates this
}
```

---

### [R-002] [CRITICAL] [tpv31.toml] — `interior_flux = "scalar"` is mutually exclusive with `material.kind = "depth_profile_1d"` — driver hard-aborts at startup

**Category:** BUG

**Description:**
This is a direct consequence of R-001 for TPV31 but warrants a separate ID because it's a hard-stop crash, not a silent semantic drift.  The driver's scalar-ctor guard at `spatial_dyn_driver.cpp:944-950` reads:

```c++
MFEM_VERIFY(material.mode == MaterialField::Mode::Constant,
            "spatial_dyn_driver: [numerics].interior_flux = \"scalar\" "
            "requires homogeneous material (Mode::Constant); current "
            "material kind cannot be reduced to scalar (λ, μ, ρ).  "
            "Either set material.kind = \"constant\" or remove the "
            "interior_flux opt-in to use the bimaterial path.  "
            "(REVIEW R-006)");
```

For `tpv31.toml`:
- `[material] kind = "depth_profile_1d"` → `MakeDepthProfile1DMaterial(...)` returns a `Mode::Coefficient` material (`spatial_dyn_driver.cpp:872-893`).
- `[numerics] interior_flux = "scalar"` → `use_scalar_ctor = true` (line 941).

The verify at line 944 fires immediately.  **TPV31 cannot run through this driver at all** with the canonical TOML, even before reaching the stress-source dispatch.

**Trigger:**
`seas_spatial_dyn_driver --config tpv31/configs/tpv31.toml`.

**Actual behavior:**
`MFEM_ABORT` at `spatial_dyn_driver.cpp:944` before any time step.

**Expected behavior:**
TPV31 runs through the heterogeneous ctor with its depth-varying `(λ(z), μ(z), ρ(z))` material — exactly the case the heterogeneous Riemann solver was built for.

**Suggested fix:**
Remove `interior_flux = "scalar"` from `tpv31.toml`.  See R-001 diff.  After the fix, the next-blocker R-003 (missing depth-proportional stress branch) will surface.

**Test case:**
```c++
// tests/integration/test_tpv31_loads.cpp
TEST_CASE("R-002: TPV31 canonical TOML does not abort during WaveOperator construction") {
    bool aborted = false;
    try {
        run_spatial_dyn_dry_run("tpv31/configs/tpv31.toml");
    } catch (const mfem_abort_exception&) { aborted = true; }
    REQUIRE_FALSE(aborted);
}
```

---

### [R-003] [CRITICAL] [spatial_dyn_driver.cpp:1167-1200 / tpv31.toml] — No driver branch for `StressSourceKind::DepthProportionalToShearModulus`; TPV31 hard-aborts at stress projection

**Category:** BUG

**Description:**
`tpv31.toml:204` declares `[stress] kind = "depth_proportional"`, which parses to `StressSourceKind::DepthProportionalToShearModulus` (`spatial_friction.cpp:284`).  The driver's stress-source dispatch at `drivers/spatial_dyn_driver.cpp:1167-1200` has only three branches:

```c++
if      (kind == ConstantTensor)                  { ... }              // line 1167
else if (kind == ConstantTensorWithPatches)       { ... }              // line 1180
else                                              { ApplyCsmStressSidecar(...) ; }  // line 1197
```

The fallback `ApplyCsmStressSidecar` (`spatial/code/spatial_stress.cpp:33-34`) hard-aborts because

```c++
MFEM_VERIFY(spec.kind == StressSourceKind::SidecarHDF5,
            "ApplyCsmStressSidecar: StressSpec.kind must be SidecarHDF5");
```

The required `spatial::DepthProportionalToShearModulusStressSource` (declared at `spatial_stress.hpp:193-227`, implemented at `spatial_stress.cpp:159-202`) is **never instantiated by this driver**.

Even if R-002 is fixed and TPV31 reaches stress projection, this is the next abort.  TPV31 (Phase R.5) is the canonical test case for `μ(depth)`-scaled pre-stress through the heterogeneous Riemann path; without this branch the whole TPV31 pipeline is dead code.

**Trigger:**
`seas_spatial_dyn_driver --config tpv31/configs/tpv31.toml` (after R-002 is fixed).

**Actual behavior:**
`MFEM_ABORT` from `ApplyCsmStressSidecar` with message "StressSpec.kind must be SidecarHDF5".

**Expected behavior:**
The driver constructs `DepthProportionalToShearModulusStressSource` with the six `sigma_*_per_mu` components and a `mu_at_xyz(x,y,z)` callback that evaluates the current `MaterialField`, then calls `geom.ComputeParams(src, ...)`.

**Suggested fix:**
Insert a new branch between `ConstantTensorWithPatches` and the `else` fallthrough.  The callback can be implemented either through a coordinate-based point locator over `pmesh` or, more cheaply for TPV31, by routing through `depth_profile_wrapper->mu_at_depth(...)` when `depth_profile_wrapper != nullptr`.

```diff
    else if (cfg.stress.kind ==
             spatial::StressSourceKind::ConstantTensorWithPatches)
    {
       ...
       geom.ComputeParams(src, ...);
    }
+   else if (cfg.stress.kind ==
+            spatial::StressSourceKind::DepthProportionalToShearModulus)
+   {
+      // TPV31-style:  σ(x,y,z) = sigma_*_per_mu * μ(x,y,z) / μ_ref.
+      auto mu_at_xyz =
+         [&material, &pmesh](real_t x, real_t y, real_t z) -> real_t {
+            return spatial::EvalMuAtPoint(material, pmesh, x, y, z);
+         };
+      const auto& dp = cfg.stress.depth_proportional;
+      spatial::DepthProportionalToShearModulusStressSource src(
+         dp.sigma_xx_per_mu, dp.sigma_yy_per_mu, dp.sigma_zz_per_mu,
+         dp.sigma_xy_per_mu, dp.sigma_yz_per_mu, dp.sigma_xz_per_mu,
+         dp.mu_ref_pa, std::move(mu_at_xyz));
+      geom.ComputeParams(src,
+                         cfg.stress.pore_pressure.P_p_pa,
+                         cfg.stress.pore_pressure.P_p_grad_pa_per_m,
+                         cfg.stress.pore_pressure.min_sigma_n_pa);
+   }
    else
    {
+      MFEM_VERIFY(cfg.stress.kind == spatial::StressSourceKind::SidecarHDF5,
+                  "spatial_dyn_driver: unhandled stress kind " <<
+                  static_cast<int>(cfg.stress.kind));
       spatial::ApplyCsmStressSidecar(cfg.stress, geom);
    }
```

`spatial::EvalMuAtPoint(material, pmesh, x, y, z)` does not yet exist; the simplest correct implementation is a point-locator over `pmesh.GetNE()` followed by `material.EvalAt(elem, T, ip, λ, μ, ρ)`.  For TPV31 production (depth-only material), a fast path routed through `depth_profile_wrapper->EvalAtDepth(depth)` is preferable.

**Test case:**
```c++
// tests/unit/test_tpv31_stress_dispatch.cpp
TEST_CASE("R-003: spatial_dyn instantiates DepthProportionalStressSource for kind=depth_proportional") {
    spatial::SpatialFrictionConfig cfg = make_minimal_tpv31_cfg();
    cfg.stress.kind = spatial::StressSourceKind::DepthProportionalToShearModulus;
    cfg.stress.depth_proportional.sigma_xy_per_mu = 30.0e6;     // 30 MPa
    cfg.stress.depth_proportional.mu_ref_pa       = 32.03812032e9;
    auto geom = run_spatial_dyn_stress_only(cfg);
    REQUIRE(geom.HasParams());
    REQUIRE(geom.GetTauPre().Size() == 2 * geom.NumFaultDOFs());
    for (int i = 0; i < geom.NumFaultDOFs(); ++i) {
        REQUIRE(std::abs(geom.GetTauPre()(2*i + 1) - 30.0e6) < 1.0);
    }
}
```

---

### [R-004] [MODERATE] [tpv*.toml × 4] — `[output].paraview_enabled` missing; intended ParaView fault output silently dropped

**Category:** BUG

**Description:**
All four TOMLs declare per-collection ParaView modes (`paraview_volume = "off"`, `paraview_bulk = "off"`, `paraview_fault = "hdf5"`) and per-collection dt values (`paraview_fault_dt = "0.001s"` etc.).  None of them set `paraview_enabled = true`.

The parser defaults `paraview_enabled = false` (`spatial_friction.cpp:1216`).  The driver enforces a hard master gate at `spatial_dyn_driver.cpp:1597-1602`:

```c++
if (!cfg.output.paraview_enabled) {
   volume_pv_enabled = false;
   fault_pv_enabled  = false;
   bulk_pv_enabled   = false;
}
```

Consequently `any_pv_requested == false`, `pv_out` is never constructed, and rank 0 prints `"ParaView output: OFF"`.  The user reads `paraview_fault = "hdf5"` in the TOML, expects an HDF5 fault file, and gets nothing.  This was explicitly designed as a hard master gate (see `spatial_paraview_compaction_parity_plan_2026-05-18.md:169-173`), but every canonical TOML on this branch fails to flip the gate, making the per-collection settings dead code unless the user passes `--paraview` on the command line.

For the SCEC TPV* benchmarks this is particularly damaging: the SCEC-trace comparison deliverable is built from the fault.vtkhdf output.  Producing zero fault files in the canonical run is a silent regression vs the native drivers, which write traces unconditionally.

**Trigger:**
`seas_spatial_dyn_driver --config <any of the four configs>` with no `--paraview*` CLI flag.

**Actual behavior:**
`ParaView output: OFF` on rank 0; the output directory contains no `ParaView/` or `ParaView_bulk/` subdirectories at the end of the run.

**Expected behavior:**
Each TOML emits the HDF5 fault output it advertises.

**Suggested fix:**
Add `paraview_enabled = true` to every TOML's `[output]` block.

```diff
 [output]
 output_dir              = "tpv205/out"
 restart_prefix          = "cp"
+paraview_enabled        = true
 paraview_volume         = "off"
 paraview_bulk           = "off"
 paraview_fault          = "hdf5"
```

(Apply identical edit to `tpv102.toml`, `tpv104.toml`, `tpv31.toml`.)

**Test case:**
```python
# tests/integration/test_toml_paraview_enabled.py
def test_R004_paraview_enabled_present(toml_path):
    cfg = parse_toml(toml_path)
    out = cfg["output"]
    has_any_mode_on = any(out.get(k, "off") != "off"
                          for k in ("paraview_volume", "paraview_bulk", "paraview_fault"))
    if has_any_mode_on:
        assert out.get("paraview_enabled") is True, (
            f"{toml_path} configures a non-off ParaView mode but does NOT "
            f"set paraview_enabled=true; the master gate will drop the output.")
```

---

### [R-005] [MODERATE] [spatial_dyn_driver.cpp::main] — Zero-fault-faces case is not aborted; native drivers do

**Category:** QUALITY / DEVIATION

**Description:**
The native drivers (e.g., `drivers/tpv205_driver.cpp:1216-1218`) hard-abort when the mesh has no faces matching `bc.fault_attr`:

```c++
MFEM_VERIFY(n_fault_g > 0,
            "No fault faces with attr=" << bc.fault_attr
            << " found.  Check --bc-fault or mesh Physical Surface tags.");
```

The spatial driver only verifies `bc.fault_attr > 0` (line 830) and reports `num_fault_global` as an informational printout (line 1136-1142).  If a TOML's `[boundary].fault_attr` does not match the mesh, the run continues silently with zero fault DOFs, executes 1200+ macro-steps doing useless wave-only work, writes a final checkpoint, and exits 0.  The user discovers the silent failure only after consuming compute budget.

**Trigger:**
TOML with `[boundary].fault_attr = 99` (any value not in the mesh's Physical Surface tags).

**Actual behavior:**
Run completes with no error.  No fault DOFs, no ParaView fault output, no station traces — but exit code 0.

**Expected behavior:**
Hard-abort with the same message the native drivers emit.

**Suggested fix:**
```diff
    if (rank == 0)
    {
       std::cout << "[fault] QPs per face = " << nbf_per_face
                 << ", num_fault_global = " << num_fault_global
                 << " (local = " << num_fault_local
                 << ", shared = " << num_shared_fault << ")\n";
    }
+   MFEM_VERIFY(num_fault_global > 0,
+               "spatial_dyn_driver: no fault faces with attr="
+               << bc.fault_attr << " found in mesh '" << cfg.mesh.path
+               << "'.  Check [boundary].fault_attr in the TOML and the "
+               "mesh's Physical Surface tags.");
```

**Test case:**
```c++
TEST_CASE("R-005: spatial_dyn aborts when bc.fault_attr matches zero mesh faces") {
    auto cfg = make_minimal_tpv205_cfg();
    cfg.boundary.fault_attr = 999;        // mesh has no Physical Surface 999
    bool aborted = false;
    try { run_spatial_dyn_minimal(cfg); }
    catch (const mfem_abort_exception &e) {
        aborted = true;
        REQUIRE(std::string(e.what()).find("no fault faces") != std::string::npos);
    }
    REQUIRE(aborted);
}
```

---

### [R-006] [MODERATE] [tpv205.toml / tpv31.toml] — Mesh files do not exist on disk; driver fails at first `Mesh smesh(...)` call

**Category:** BUG

**Description:**
- `tpv205.toml:115` references `tpv205/mesh/tpv2053d_200m.msh`; the directory contains only `tpv2053d_200m.geo` and `tpv2053d_100m.geo`.
- `tpv31.toml:176` references `tpv31/mesh/tpv31_50m.msh`; the directory contains only `tpv31_50m.geo`.

TPV102 and TPV104 ship their `.msh` files (`tpv102_1000m.msh`, `tpv104_1000m.msh`, both 10.6 MB).  The asymmetry is a user-experience bug — a developer running the canonical TPV205/TPV31 configs gets an unhelpful MFEM file-open error rather than a clear "regenerate the mesh" message.

Both TOMLs document the regeneration command in their `[mesh]` block comment, but the driver itself surfaces no preflight check.  Additionally, per `miniapps/seas/CLAUDE.md` "Known limitation — Gmsh `.msh` format", MFEM requires Gmsh v2.2 format (`-format msh22`); a developer who forgets the flag and generates v4 sees a misleading "vertices indices are not unique" error from `mfem/mesh/mesh_readers.cpp:1628`.

**Trigger:**
Fresh checkout (or any environment without the regenerated mesh):
```
seas_spatial_dyn_driver --config tpv205/configs/tpv205.toml
seas_spatial_dyn_driver --config tpv31/configs/tpv31.toml
```

**Actual behavior:**
`Mesh smesh(cfg.mesh.path.c_str(), 1, 1);` (line 801) aborts with an MFEM file-open or parse error.

**Expected behavior:**
Either (a) commit the regenerated `.msh` artifacts (matches the TPV102/TPV104 convention), or (b) emit a preflight error message naming the regeneration command.

**Suggested fix:**
Preferred: commit `tpv2053d_200m.msh` and `tpv31_50m.msh` (v2.2 format).  Failing that, add a preflight check before `Mesh smesh(...)`:

```diff
+   if (rank == 0 && !std::filesystem::exists(cfg.mesh.path)) {
+      std::cerr << "ERROR: mesh file '" << cfg.mesh.path
+                << "' does not exist.  Generate via:\n"
+                << "    gmsh -format msh22 -3 <input>.geo -o <output>.msh\n"
+                << "  (Gmsh v2.2 — required by MFEM, see CLAUDE.md "
+                << "\"Known limitation — Gmsh .msh format\".)\n";
+   }
+   MPI_Barrier(comm);
    Mesh smesh(cfg.mesh.path.c_str(), 1, 1);
```

And add a Makefile rule per benchmark:

```makefile
tpv205/mesh/tpv2053d_200m.msh: tpv205/mesh/tpv2053d_200m.geo
	gmsh -format msh22 -3 $< -o $@
tpv31/mesh/tpv31_50m.msh: tpv31/mesh/tpv31_50m.geo
	gmsh -format msh22 -3 $< -o $@
```

**Test case:**
```python
# tests/integration/test_mesh_artifacts_committed.py
def test_R006_referenced_mesh_files_exist():
    for toml_path in ["tpv102/configs/tpv102.toml",
                      "tpv104/configs/tpv104.toml",
                      "tpv205/configs/tpv205.toml",
                      "tpv31/configs/tpv31.toml"]:
        cfg = parse_toml(toml_path)
        mesh_path = cfg["mesh"]["path"]
        assert Path(mesh_path).is_file(), (
            f"{toml_path} references mesh '{mesh_path}' which does not exist.")
```

---

### [R-007] [MODERATE] [spatial_dyn_driver.cpp::stations] — TPV102 / TPV104 SCEC station traces are not written (only TPV205 is wired)

**Category:** DEVIATION / FEATURE GAP

**Description:**
The spatial driver wires `TPV205StationWriter` (lines 2158-2182) only when `cfg.problem.tag == "tpv205"`.  The native TPV102 (`drivers/tpv102_driver.cpp:2564-2574`) and TPV104 drivers each write their own SCEC trace files, which downstream comparison scripts consume directly.  Running TPV102 or TPV104 through `seas_spatial_dyn_driver` produces **no station trace files**, breaking SCEC cross-code comparison parity.

`cfg.problem.tag = "tpv102"` and `"tpv104"` ARE set in the respective TOMLs, but the driver has no `if (cfg.problem.tag == "tpv102")` / `"tpv104"` branch.

This is a moderate issue because it removes the SCEC-deliverable output for these two benchmarks; the run otherwise produces correct physics (via the heterogeneous Riemann path, per R-001's intent).

**Trigger:**
`seas_spatial_dyn_driver --config tpv102/configs/tpv102.toml`.

**Actual behavior:**
The `tpv102/out` directory contains no `tpv102_traces_*.dat` files.  The native TPV102 driver writes one per station (typically ~16-64 stations).

**Expected behavior:**
Per-station SCEC trace files written analogously to the TPV205 path.

**Suggested fix:**
Add tpv102/tpv104 station-writer branches mirroring lines 2158-2182.  `DefaultStations_TPV102` / `Tpv102StationWriter` (and TPV104 equivalents) already exist in `dynamic/tpv102_setup.hpp` / `dynamic/tpv104_setup.hpp`.

```diff
+   Tpv102StationWriter tpv102_station_writer;
+   const bool tpv102_stations_active = (cfg.problem.tag == "tpv102");
+   if (tpv102_stations_active) {
+      auto stations = DefaultStations_TPV102();
+      tpv102_station_writer.Open(cfg.output.output_dir, "tpv102",
+                                 stations, fault_coords,
+                                 num_fault_local, comm);
+      if (restart_prefix.empty())
+         tpv102_station_writer.WriteStep(cfg.time.t_initial, dof_data);
+   }
+   // analogous block for TPV104
```
And hoist per-step writes alongside `tpv205_station_writer.WriteStep`.

**Test case:**
```c++
TEST_CASE("R-007: spatial_dyn produces SCEC traces for TPV102") {
    run_spatial_dyn(config="tpv102/configs/tpv102.toml",
                    tfinal=0.1, output_dir="/tmp/tpv102_test");
    REQUIRE(std::filesystem::exists("/tmp/tpv102_test/tpv102_traces_x2_0_x3_7.5.dat"));
}
```

---

### [R-008] [LOW] [spatial_dyn_driver.cpp:1363 / spatial_setup.hpp::InitializeFaultDOFs_Spatial_RS] — RS init uses element-centroid material evaluation, not fault-QP IP

**Category:** ASSUMPTION

**Description:**
For the LSW path (line 1355), the driver passes the per-DOF `dof_ips` cache so `InitializeFaultDOFs_Spatial` evaluates `MaterialField::EvalAt` at the actual fault QP.  For the RS path (line 1363), the driver calls the non-IP-aware overload — `seed_static_dof_fields` then falls back to the bulk element CENTROID (`spatial_setup.hpp:67-69`).

For TPV102 / TPV104 (rate-state + constant material) this is harmless.  But it is a latent bug for any future RS config combining `material.kind = "depth_profile_1d"` (or `"sidecar_hdf5"`) with rate-state friction: the per-DOF impedances `d.Zp_*, d.Zs_*, d.eta_p, d.eta_s` will be biased by the element-scale heterogeneity.  This is documented as a "best per-element stop-gap" in the inline comment at lines 64-66 of `spatial_setup.hpp`, but the driver should already be using the IP-aware overload because that information is available.

**Trigger:**
`material.kind = "depth_profile_1d"` AND `meta.law = "rate_state"` in the same TOML.

**Actual behavior:**
Per-DOF impedances evaluated at the bulk centroid.

**Expected behavior:**
Per-DOF impedances evaluated at the fault QP IP, matching the LSW path.

**Suggested fix:**
Add an IP-aware overload of `InitializeFaultDOFs_Spatial_RS` (analogous to the LSW one at `spatial_setup.hpp:358-394`) and call it from the driver:

```diff
- spatial::InitializeFaultDOFs_Spatial_RS<ParMesh>(
-    dof_data, num_fault_total, dof_to_elem, material, pmesh,
-    rs, geom.GetTauPre(), geom.sigma_n_per_dof());
+ spatial::InitializeFaultDOFs_Spatial_RS<ParMesh>(
+    dof_data, num_fault_total, dof_to_elem, material, pmesh,
+    rs, geom.GetTauPre(), geom.sigma_n_per_dof(),
+    dof_ips);
```

**Test case:**
```c++
TEST_CASE("R-008: RS init evaluates material at fault QP IP, not centroid") {
    // Build a tet whose centroid is at depth 1500 m but whose fault QP is at depth 750 m,
    // using a depth_profile_1d material that differs between those depths.
    auto cfg = make_tpv_rs_depth_profile_cfg();
    auto dof_data = run_spatial_dyn_rs_init(cfg);
    REQUIRE(dof_data[0].Zs_plus == Approx(expected_Zs_at_depth_750).margin(1e-6));
}
```

---

### [R-009] [LOW] [tpv205.toml — barrier rules] — Boundary DOFs at x=±15 km / z=−15 km are marked as barrier; native driver puts them inside the rupture area

**Category:** EDGE_CASE / DEVIATION

**Description:**
The native `InRuptureArea_TPV205` (`config/tpv205_params.hpp:158-163`) uses `<=` on the boundary:
```c++
return std::abs(along_strike) <= 15.0e3 && down_dip <= 15.0e3 && down_dip >= 0.0;
```
The TOML barrier rules use `x_max_m = -15000.0` etc., which the resolver interprets as `x ≤ -15000` (`SpatialRule::matches`).  At DOFs exactly on the boundary (`x = ±15000` or `z = -15000`), the spatial driver flags them as **barrier** (μ_s = 1.0e6 sentinel), but the native driver puts them inside the rupture area (μ_s = 0.677).

A 200 m structured mesh places nodes exactly at x = ±15 km and z = −15 km, so the affected set is non-trivial (on the order of dozens of QPs along each barrier edge).  The downstream effect is that the rupture-area edge is shifted inward by half an element, slightly altering the rupture-front kinematics at the boundary.

**Trigger:**
A DOF at exactly `(x, y, z) = (-15000, 0, -7500)`.

**Actual behavior:**
`mu_s = 1.0e6` (barrier sentinel).

**Expected behavior:**
`mu_s = 0.677` (interior rupture area).

**Suggested fix:**
Nudge the bounds inward by a small epsilon so the boundary DOFs fall into the rupture area:

```diff
 [[friction.slip_weakening.spatial]]
 kind     = "barrier"
-x_max_m  = -15000.0   # left  barrier: x < -15 km
+x_max_m  = -15000.001 # left  barrier (epsilon << mesh spacing)

 [[friction.slip_weakening.spatial]]
 kind     = "barrier"
-x_min_m  =  15000.0   # right barrier: x > +15 km
+x_min_m  =  15000.001

 [[friction.slip_weakening.spatial]]
 kind     = "barrier"
-z_max_m  = -15000.0
+z_max_m  = -15000.001  # depth > 15 km
```

Alternative: extend `SpatialRule::matches` to use strict `<` for `Kind::Barrier` only.

**Test case:**
```c++
TEST_CASE("R-009: TPV205 boundary DOFs at x=±15km / z=−15km are inside the rupture area") {
    auto lsw = make_tpv205_lsw_block();
    Vector coords(3); coords(0) = -15000.0; coords(1) = 0.0; coords(2) = -7500.0;
    Array<int> attrs(1); attrs[0] = 103;
    auto p = SpatialFrictionResolver{}.ResolveSlipWeakening(lsw, coords, attrs);
    REQUIRE(p.mu_s(0) == Approx(0.677));   // NOT 1.0e6 sentinel
}
```

---

## Summary

- **Critical issues:** 3
  - R-001 (TPV205 & TPV31 opt out of heterogeneous Riemann via `interior_flux = "scalar"`).
  - R-002 (TPV31 hard-aborts at startup — `scalar` ctor + `depth_profile_1d` material are mutually exclusive).
  - R-003 (no driver branch for `DepthProportionalToShearModulus` — TPV31 hard-aborts at stress projection).
- **Moderate issues:** 4
  - R-004 (`paraview_enabled` missing from all four TOMLs; ParaView output silently dropped).
  - R-005 (zero-fault-faces case not aborted; native drivers do).
  - R-006 (TPV205 / TPV31 mesh artifacts not committed).
  - R-007 (TPV102 / TPV104 SCEC station traces not produced by spatial driver).
- **Low issues:** 2
  - R-008 (RS init uses element-centroid material eval, not fault-QP IP).
  - R-009 (TPV205 barrier-rule boundary off-by-one vs native driver).
- **Plan compliance:** PARTIAL — the heterogeneous Riemann interior-flux path is implemented (Phase R.2/R.4) and TPV102/TPV104 reach it via the parser default, but the TPV205 and TPV31 TOMLs explicitly route around it (R-001/R-002).  TPV31's depth-proportional stress source is not wired (R-003).  Output deliverables (ParaView fault, SCEC station traces) are missing for the canonical configs (R-004/R-007).
- **Verdict:** FAIL — TPV205 and TPV31 do not exercise the heterogeneous Riemann solver despite that being the branch's stated purpose (R-001), and TPV31 cannot run at all (R-002 + R-003).  Must fix R-001, R-002, R-003, and R-004 before any canonical TPV* run on this branch is meaningful.

## Unreviewed Areas

- `WaveOperator<MeshT>(MaterialField, BoundaryConfig)` interior-face dispatch internals — assumed correct per the existing Phase R parity tests (`T-PHASEH-SCALAR-PARITY`, `test_phaser_dispatch_smoke`); not re-verified line-by-line.
- The `BimaterialFlux::ApplyPerFaceFlux` precomputation across MPI shared faces — relied on the inline comment at `wave_operator.inl:349-350` and CLAUDE.md "Bi-material Riemann (Phase R)" section.
- `ApplyCsmStressSidecar` end-to-end behavior — only the kind guard was inspected; the actual sidecar I/O path is unchanged from prior work and assumed correct.
- Checkpoint / restart correctness on this branch — not exercised in this review (the schema R-105 fix is referenced in the code comments but the actual restart flow was not run).
- ZFP / deflate compression flag wiring — visually inspected for parity with the native drivers; no test executed.
- `Tpv104SubStepIterator::AdvanceWithSubStepStates` callback overload — assumed to mirror the `Tpv102` and `Tpv205` variants based on call-site signature, not re-verified line-by-line.
