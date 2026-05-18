# Code Review: spatial_dynamic_rupture_plan rev-3 (2026-05-18)

## Review Scope

- **Plan:** `safs/project_7.0_alternative/document/spatial_dynamic_rupture_plan.{md,pdf}` (rev-3, 2026-05-18, 152 KB / 3020 LOC)
- **Files reviewed:**
  - `miniapps/seas/spatial/code/spatial_friction.{hpp,cpp}` (Phase 0/1)
  - `miniapps/seas/spatial/code/spatial_velocity.{hpp,cpp}` (Phase 2)
  - `miniapps/seas/spatial/code/spatial_stress.{hpp,cpp}` (Phase 3/3b)
  - `miniapps/seas/fault/fault_geometry_safs_templated.inl` (Phase 3b)
  - `miniapps/seas/tests/unit/test_spatial_friction_config.cpp`
  - `miniapps/seas/tests/unit/test_spatial_friction_resolver.cpp`
  - `miniapps/seas/tests/unit/test_spatial_velocity_bundle.cpp`
  - `miniapps/seas/tests/unit/test_spatial_stress_bundle.cpp`
  - `miniapps/seas/tests/unit/test_spatial_constant_stress_source.cpp`
  - `miniapps/seas/safs/project_7.0_alternative/document/spatial_friction_config_schema.md`
  - `miniapps/seas/safs/project_7.0_alternative/friction/EXAMPLE_*.toml` (×3)
- **Domain context consulted:**
  - `miniapps/seas/CLAUDE.md` (sign conventions, fault basis convention, CFL rules)
  - `miniapps/seas/io/field_coefficient.{hpp,cpp}` (existing `FieldProjector::ProjectFaultPreStress`)
  - `miniapps/seas/io/material_coefficients.hpp` (Lambda/Mu/Rho FromSidecar)
  - `miniapps/seas/io/data_field_3d.hpp` (DataField3D field naming)
  - `miniapps/seas/fault/fault_geometry.hpp` + `fault_geometry_safs.inl`
  - Tandem reference: `/Users/chunhuizhao/projects/tandem/`
  - SAFS sidecar builder: `velocity/code/build_velocity_cvmh.py` (field-name conformance)

## Findings

---

### [R-001] [CRITICAL] [spatial_dynamic_rupture_plan.md:1706+ / Phases H, 4, 5, 6] — Major plan scope NOT IMPLEMENTED

**Category:** DEVIATION

**Description:**
Phases 0, 1, 2, 3, and 3b have implementations on disk, but the plan's downstream phases — **Phase H** (heterogeneous `WaveOperator` ctor + `GodunovFluxPool` + `LSW_ForcedRupture` dispatch), **Phase 4** (the `seas_spatial_dyn_driver` executable), **Phase 5** (`dynamic/spatial_setup.hpp` free functions + new `FaultGeometry` ctor overload + `DRIVER_TAG_V1` checkpoint extension), and **Phase 6** (sbatches + verify/compare scripts) — are absent.  Specifically, none of the following exists:

| Expected file (plan-required) | Status |
|-------------------------------|--------|
| `dynamic/godunov_flux_pool.{hpp,cpp}`              | MISSING |
| New `WaveOperator(MeshType&, int, const MaterialField&, const BoundaryConfig&)` ctor | NOT in `dynamic/wave_operator.hpp` |
| `FaultFrictionLaw::LSW_ForcedRupture` enum value   | NOT in `dynamic/wave_operator.hpp` |
| `FaultFaceFlux::EvaluateADER_LSW_ForcedRupture(...)` | NOT in `dynamic/fault_face_flux.{hpp,cpp}` |
| `DOFData::T_forced_rupture` / `t0_decay_forced`    | NOT in `dynamic/fault_face_flux.hpp` |
| `drivers/spatial_dyn_driver.cpp`                   | MISSING |
| `spatial/code/spatial_dyn_driver_init.{hpp,cpp}`   | MISSING |
| `dynamic/spatial_setup.hpp`                        | MISSING |
| New `FaultGeometry` ctor taking pre-built per-DOF arrays | NOT in `fault/fault_geometry.hpp` |
| `fault_dof_ip(i)` accessor + `fault_dof_ip_` storage | MISSING |
| `WriteTpv104Checkpoint(..., driver_tag)` extension | NOT in `io/tpv104_checkpoint.hpp` |
| `jobs/safs/spatial_dyn_*.sbatch` (×3)              | MISSING |
| `spatial/code/scripts/verify_spatial_dyn_smoke_safs.py` | MISSING |
| `spatial/code/scripts/compare_spatial_dyn_velocity_models.py` | MISSING |
| `tests/unit/test_phaseh_*.cpp` (×4)                | MISSING |
| `tests/unit/test_spatial_setup.cpp`                | MISSING |

Only `spatial/code/scripts/verify_constant_tensor_projection.py` was created in the scripts directory.

The Phase-1 helper `LSWFrictionCoefficient_ForcedRupture(...)` is implemented in `spatial_friction.hpp`, but it currently has **no consumer**: Phase H's `EvaluateADER_LSW_ForcedRupture` (which would call the helper) does not exist.

**Trigger:**
Running `make seas_spatial_dyn_driver` or `make test-spatial-dyn-driver` will fail — the targets are not defined in `miniapps/seas/Makefile`.

**Actual behavior:**
The work-in-progress committed on this branch covers the data-preparation layers (TOML parser, per-DOF friction resolver, velocity bundle, CSM stress wrapper, constant-tensor stress source + templated `ComputeSAFSParams<StressSource>` overload).  Nothing downstream of those layers is built; the dynamic-rupture driver does not exist.

**Expected behavior:**
Per `spatial_dynamic_rupture_plan.md` §Preconditions (rev-3), Phase H, Phase 5, and Phase 4 are the in-plan deliverables — rev-3 explicitly removed the prior `static_assert` Precondition gates so the plan is supposed to be self-contained.  Phase 4 ships the executable; Phase 6 ships the verification artifacts.

**Suggested fix:**
This finding is informational — it cannot be "fixed" by editing one file; it is a status report that the implementation is at Phase 3b of 6.  Either (a) implement Phases H, 5, 4, 6 as scheduled, or (b) explicitly mark the in-tree implementation as "Phases 0–3b only" and remove Phase H/4/5/6 acceptance criteria from the "must-pass before merge" list in any companion review/QA doc.  No diff is provided.

```diff
- (no actionable diff for this finding; this is a scope-completeness report)
```

**Test case:**
```python
def test_R001_plan_phases_present():
    # Phase H sentinel
    assert (REPO / "miniapps/seas/dynamic/godunov_flux_pool.hpp").exists()
    # Phase 4 sentinel
    assert (REPO / "miniapps/seas/drivers/spatial_dyn_driver.cpp").exists()
    # Phase 5 sentinel
    assert (REPO / "miniapps/seas/dynamic/spatial_setup.hpp").exists()
    # Phase 6 sentinel
    assert (REPO / "miniapps/seas/spatial/code/scripts/verify_spatial_dyn_smoke_safs.py").exists()
```

---

### [R-002] [CRITICAL] [spatial_friction.cpp:resolve_forced_impl] — `ip.Init(0)` evaluates material at reference-element CORNER, not centroid

**Category:** BUG

**Description:**
`ResolveForcedRupture` per-DOF code evaluates the heterogeneous material with an `IntegrationPoint` initialized via `ip.Init(0)`:

```cpp
IntegrationPoint ip;
ip.Init(0);   // reference centroid placeholder (R-111: when
              // FaultGeometry::fault_dof_ip(i) is wired, use it).
real_t lam_e, mu_e, rho_e;
material.EvalAt(e, *T, ip, lam_e, mu_e, rho_e);
```

`IntegrationPoint::Init(int i)` (`mfem/fem/intrules.hpp`) is defined as:

```cpp
void Init(int const i) { x = y = z = weight = 0.0; index = i; }
```

So `ip.Init(0)` sets the reference-space point to **`(0, 0, 0)` — the origin of the reference element, which is a CORNER for both reference tets and reference hexes**.  The reference centroids are `(1/4, 1/4, 1/4)` and `(1/2, 1/2, 1/2)` respectively.

For `MaterialField::MakeConstant(...)` the value is constant inside the element and the bug is silent (every IP gives the same `lam, mu, rho`).  For `MaterialField::MakeCoefficient(...)` — the actual SAFS path that wraps `LambdaFromSidecar`/`MuFromSidecar`/`RhoFromSidecar` over `DataField3D` trilinear-interpolated CVM voxels — `Coefficient::Eval(T, ip)` internally calls `T.Transform(ip, phys)` to get the physical (x, y, z), so evaluating at the corner returns the **corner's** Vp/Vs/ρ, not the DOF's.  The resulting `Vs` is wrong by the heterogeneity scale of one element.

**Trigger:**
Any CVMH/CVM-S sidecar run with `[nucleation]` enabled and `MaterialField::MakeCoefficient(...)` (the default driver path per plan §Step 8).

**Actual behavior:**
`T_forced(r) = r / (0.7 · Vs_corner) + …`, where `Vs_corner` is the Vs at the bulk element's reference-origin corner — could be a different layer than where the fault DOF actually sits.  For a typical sediment basin with Vs ≈ 1 km/s near the surface and Vs ≈ 3.5 km/s deeper, the relative error in `T_forced` is up to ~250%.

**Expected behavior:**
Evaluate `material.EvalAt` at the fault DOF's actual reference IP, i.e. `FaultGeometry::fault_dof_ip(i)` (which Phase 5 R-111 would have provided).  Until Phase 5 wires that, the closest correct stop-gap is the element centroid (still imperfect because the DOF is on a face, but at least inside the element):

**Suggested fix:**
```diff
- IntegrationPoint ip;
- ip.Init(0);   // reference centroid placeholder (R-111: when
-               // FaultGeometry::fault_dof_ip(i) is wired, use it).
+ IntegrationPoint ip;
+ // True reference centroid for the bulk element's geometry.  For
+ // heterogeneous material the centroid is the best per-element stop-
+ // gap until Phase 5 wires FaultGeometry::fault_dof_ip(i) (R-111).
+ const Geometry::Type gtype = mesh.GetElementBaseGeometry(e);
+ ip = *Geometries.GetCenter(gtype);
  real_t lam_e, mu_e, rho_e;
  material.EvalAt(e, *T, ip, lam_e, mu_e, rho_e);
```
(Apply the same change to the eta_auto path in `resolve_rs_impl` — see [R-003].)

Long-term fix: extend the resolver signature to accept a `const std::vector<IntegrationPoint>* fault_dof_ip` parameter and, when non-null, use it instead of the centroid.  Drivers can then pass `geom.fault_dof_ip()` once Phase 5 ships the accessor.

**Test case:**
```cpp
// tests/unit/test_spatial_friction_resolver.cpp
static void F_3_layered_material_uses_correct_Vs_at_DOF()
{
   // Build a 2-layer material with Vs_top = 1000 m/s above z=0 and
   // Vs_bot = 3500 m/s below z=0, mapped through a CoefficientFromExpr.
   // Place a fault DOF at z=-100m (well inside the bottom layer).  Build
   // a single-element bulk mesh whose corner sits at z=-200 (bottom
   // layer) and centroid at z=-100 (bottom layer).
   // With Vs=3500 m/s, T_forced(r=2km) should be ~2000/(0.7*3500) =
   // 0.82s.  With the buggy corner evaluation it could be similar if the
   // corner is also bottom-layer; the test must place the corner in a
   // DIFFERENT layer than the DOF (e.g., a hex straddling z=0 with the
   // base corner at z=-200 sees ip.Init(0) -> corner at (x,y,z=-200),
   // which the test exercises).
   ...
   const real_t T_correct = r / (0.7 * 3500.0)
                           + 0.081 * r_crit / (0.7 * 3500.0)
                             * (1.0 / (1.0 - (r/r_crit)*(r/r_crit)) - 1.0);
   TEST_NEAR(p.T_forced_s(0), T_correct, 1e-9,
             "T_forced uses Vs at the DOF, not at the corner");
}
```

---

### [R-003] [CRITICAL] [spatial_friction.cpp:resolve_rs_impl] — same `ip.Init(0)` bug in `eta_auto` path

**Category:** BUG

**Description:**
`resolve_rs_impl` re-uses the same buggy idiom for the `eta_auto = true` path:

```cpp
if (cfg.eta_auto)
{
   ...
   IntegrationPoint ip;
   ip.Init(0);   // reference centroid placeholder; for true
                 // per-DOF reference IP, FaultGeometry::fault_dof_ip(i)
                 // (Phase 5) will be wired in.
   real_t lam_e, mu_e, rho_e;
   material.EvalAt(e, *T, ip, lam_e, mu_e, rho_e);
   ...
   eta_i = 0.5 * std::sqrt(mu_e * rho_e);
}
```

Same root cause as [R-002]: `ip.Init(0)` lands at the reference corner, not the centroid, so `mu_e` and `rho_e` come from the wrong physical location for `Mode::Coefficient`.  `eta_i = 0.5 · √(μ·ρ)` therefore picks up the corner's impedance instead of the DOF's.

**Trigger:**
Rate-and-state runs with `[friction.rate_state] eta = "auto"` and CVMH/CVM-S sidecar material (the default per the rate-state EXAMPLE TOML).

**Actual behavior:**
Per-DOF eta is computed from corner material, biasing the radiation-damping coefficient at every DOF.  Effect on slip rates depends on the layer contrast across the corner→DOF edge, but is the same order as the heterogeneity scale.

**Expected behavior:**
Evaluate at the DOF location (or, as stop-gap, the element centroid — see [R-002]).

**Suggested fix:**
```diff
- IntegrationPoint ip;
- ip.Init(0);   // reference centroid placeholder; for true
-               // per-DOF reference IP, FaultGeometry::fault_dof_ip(i)
-               // (Phase 5) will be wired in.
+ IntegrationPoint ip;
+ const Geometry::Type gtype = mesh.GetElementBaseGeometry(e);
+ ip = *Geometries.GetCenter(gtype);
```

**Test case:**
```cpp
static void R_9_eta_auto_uses_correct_material_at_DOF()
{
   // Mirror F_3 with eta_auto; assert that eta_i ≈ 0.5 * sqrt(mu_DOF *
   // rho_DOF), NOT 0.5 * sqrt(mu_corner * rho_corner), when the two
   // differ by ≥ 10x.
   ...
}
```

---

### [R-004] [MODERATE] [spatial_friction.cpp:parse_root] — Mandatory blocks silently accept struct defaults instead of aborting

**Category:** ASSUMPTION / DEVIATION

**Description:**
The schema (`spatial_friction_config_schema.md` §"Top-level layout") declares **`[stress]`, `[time]`, `[material_constant_fallback]`, `[numerics]`** as **mandatory** blocks.  But `parse_root` wraps each block in `if (root.contains("<block>")) { ... }` and only validates fields **inside** the conditional.  When the block is omitted entirely:

| Block | Result of omission |
|-------|--------------------|
| `[stress]`                       | `kind = ConstantTensor`, all 6 σ = 0 (silent zero pre-stress) |
| `[time]`                         | `tfinal = 12.0`, `dt_max = 0.1` (silently use defaults) |
| `[material_constant_fallback]`   | `lambda = mu = 32e9`, `rho = 2670` (silently use crust defaults) |
| `[numerics]`                     | `ader_order = 2`, `cfl = 0.5`, `mixed_flux = "none"` (silently use defaults) |

For `[mesh]`, `[velocity]`, and `[output]` the post-block `MFEM_VERIFY` does fire on empty-string defaults, so they are effectively required.  The four blocks above are not.

The `[stress]` case is the worst: a user who forgets the entire block gets a constant-tensor stress source with zero everywhere — the simulation will run, write outputs, and produce zero traction with no warning.  The intra-block validators (e.g. "constant_tensor requires all six sigma_*_pa keys") are skipped because we never enter the conditional.

**Trigger:**
Any TOML file missing one of the four mandatory blocks.

**Actual behavior:**
Parse succeeds; silent defaults are used.

**Expected behavior:**
`MFEM_ABORT("Top-level [<block>] block missing")` mirroring the existing `MFEM_VERIFY(root.contains("meta"), ...)` pattern.

**Suggested fix:**
Add a per-block presence guard before each `if (root.contains("...")) {` body in `parse_root`:

```diff
+ MFEM_VERIFY(root.contains("material_constant_fallback"),
+             "Top-level [material_constant_fallback] block missing");
  if (root.contains("material_constant_fallback"))
  {
     const auto& m = root.at("material_constant_fallback");
     ...
  }
+ MFEM_VERIFY(root.contains("mesh"),    "Top-level [mesh] block missing");
+ MFEM_VERIFY(root.contains("velocity"),"Top-level [velocity] block missing");
+ MFEM_VERIFY(root.contains("stress"),  "Top-level [stress] block missing");
+ MFEM_VERIFY(root.contains("numerics"),"Top-level [numerics] block missing");
+ MFEM_VERIFY(root.contains("time"),    "Top-level [time] block missing");
+ MFEM_VERIFY(root.contains("output"),  "Top-level [output] block missing");
```

(Once each guard is in place, the post-block `MFEM_VERIFY`s on empty-string defaults can stay as defense-in-depth.)

**Test case:**
```cpp
// tests/unit/test_spatial_friction_config.cpp — add four sub-tests:
static void T_13_missing_stress_block_aborts()
{
   std::string toml = MinimalLSWHeader();
   // Strip the [stress] block from the header
   const auto stress_pos = toml.find("[stress]");
   const auto next_block = toml.find("[", stress_pos + 1);
   toml.erase(stress_pos, next_block - stress_pos);
   toml += MinimalLSWBlock();
   TEST_ASSERT(ParseAbortsInChild(toml),
               "missing [stress] block must abort");
}
// (similar for [time], [material_constant_fallback], [numerics])
```

---

### [R-005] [MODERATE] [test_spatial_stress_bundle.cpp:S_2,S_3] — Tests do NOT exercise `ApplyCsmStressSidecar`; they call `MFEM_VERIFY` directly

**Category:** EDGE_CASE / QUALITY (test coverage gap)

**Description:**
`S_2_empty_path_aborts` and `S_3_wrong_kind_aborts` claim to test `ApplyCsmStressSidecar`'s pre-flight checks but do not actually call the function:

```cpp
// S-2 body:
const bool aborted = RunInChild([]()
{
   MFEM_VERIFY(!std::string("").empty(),  // tautology that runs the
                                          // same MFEM_VERIFY macro path
                                          // ApplyCsmStressSidecar uses
                                          // when spec.sidecar_path == ""
               "ApplyCsmStressSidecar: spec.sidecar_path is empty");
});
TEST_ASSERT(aborted, "empty path MFEM_VERIFY must abort");

// S-3 body:
const bool aborted = RunInChild([]()
{
   StressSpec spec;
   spec.kind = StressSourceKind::ConstantTensor;
   MFEM_VERIFY(spec.kind == StressSourceKind::SidecarHDF5,
               "ApplyCsmStressSidecar: StressSpec.kind must be SidecarHDF5");
});
```

Both tests trigger `MFEM_VERIFY(false, ...)` unconditionally and assert it aborts.  Neither calls `ApplyCsmStressSidecar(spec, geom)`.  Result: the actual pre-flight guards inside `apply_csm_impl` (lines 31–45 of `spatial_stress.cpp`) are completely untested.  The plan §Phase 3 §Acceptance demands **3 tests — sidecar loads, ComputeSAFSParams populates correctly, zero-normal pre-flight aborts**; only the first criterion (S-1) is genuinely exercised.

**Trigger:**
Any silent breakage of the real pre-flight logic inside `apply_csm_impl` would go undetected (e.g., if someone refactors the abort message text or comment-outs the `MFEM_VERIFY(geom.NumZeroNormalFallbacks() == 0, ...)` check).

**Actual behavior:**
S-2 and S-3 always pass because `MFEM_VERIFY(false, ...)` always aborts.  They are tautological.

**Expected behavior:**
Tests should construct a `StressSpec` and a `FaultGeometry<Mesh>` (the serial overload is available), then call `ApplyCsmStressSidecar(spec, geom)` inside the forked child.

**Suggested fix:**
```diff
 static void S_2_empty_path_aborts()
 {
    std::cout << "\n[S-2] empty sidecar_path triggers MFEM_VERIFY abort\n";
    const bool aborted = RunInChild([]()
    {
-      MFEM_VERIFY(!std::string("").empty(),  // tautology ...
-                  "ApplyCsmStressSidecar: spec.sidecar_path is empty");
+      StressSpec spec;
+      spec.kind = StressSourceKind::SidecarHDF5;
+      spec.sidecar_path = "";                      // EMPTY → must abort
+      // Build a 1-element serial mesh + a BP5-style FaultGeometry so
+      // ApplyCsmStressSidecar has something to consume.
+      mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
+         1, 1, 1, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0);
+      // (Construct the FaultGeometry however the existing
+      //  test_compute_safs_params test does.)
+      FaultGeometry<mfem::Mesh> geom = MakeDummyFaultGeometry(mesh);
+      ApplyCsmStressSidecar(spec, geom);            // ← THE ACTUAL CALL
    });
    TEST_ASSERT(aborted, "empty path MFEM_VERIFY must abort");
 }
```
(Same shape for S-3: set `spec.kind = ConstantTensor` and call the function.)

If constructing a FaultGeometry in this unit test is infeasible without invasive scaffolding, gate S-2/S-3 behind `#ifdef MFEM_USE_MPI` and add a `FaultGeometry`-construction helper in a shared test fixture file.

**Test case:**
(See the suggested fix — the corrected S-2/S-3 ARE the test cases.)

---

### [R-006] [MODERATE] [spatial_friction.cpp:parse_root, time block] — Missing `dt_initial` validator

**Category:** EDGE_CASE

**Description:**
The schema (`spatial_friction_config_schema.md` line 154) validates `dt_initial`: `> 0` OR string `"auto"`.  The parser parses `dt_initial` via `toml_time_seconds(t, "dt_initial", -1.0)` (line 577) which returns `-1.0` for `"auto"` and any parsed numeric value otherwise.  But **no validator is applied** to the parsed result:

```cpp
cfg.time.dt_initial = toml_time_seconds(t, "dt_initial", -1.0);
// (no MFEM_VERIFY on dt_initial)
```

So `dt_initial = "0"`, `dt_initial = "-5s"`, or `dt_initial = 0.0` all pass.  Only `dt_initial = -1.0` should be the "auto" sentinel; anything else `<= 0` (other than `-1`) is meaningless.

**Trigger:**
A user typo like `dt_initial = "0s"` or `dt_initial = -0.001`.

**Actual behavior:**
Parser accepts the value; whether the (eventual) driver crashes or runs depends on Phase 4 wiring which is not implemented yet.

**Expected behavior:**
Reject any non-`-1` non-positive value:

```cpp
MFEM_VERIFY(cfg.time.dt_initial > 0.0 || cfg.time.dt_initial == -1.0,
            "[time].dt_initial must be > 0 or the string \"auto\"; "
            "got " << cfg.time.dt_initial);
```

**Suggested fix:**
```diff
   if (root.contains("time"))
   {
      const auto& t = root.at("time");
      cfg.time.tfinal     = toml_time_seconds(t, "tfinal",     12.0);
      cfg.time.t_initial  = toml_real        (t, "t_initial",  0.0);
      cfg.time.dt_initial = toml_time_seconds(t, "dt_initial", -1.0);
      cfg.time.dt_max     = toml_time_seconds(t, "dt_max",     0.1);
   }
   MFEM_VERIFY(cfg.time.tfinal > 0.0,
               "[time].tfinal must be > 0; got " << cfg.time.tfinal);
   MFEM_VERIFY(cfg.time.dt_max > 0.0,
               "[time].dt_max must be > 0; got " << cfg.time.dt_max);
   MFEM_VERIFY(cfg.time.t_initial >= 0.0,
               "[time].t_initial must be >= 0; got " << cfg.time.t_initial);
+  MFEM_VERIFY(cfg.time.dt_initial > 0.0
+              || cfg.time.dt_initial == static_cast<real_t>(-1.0),
+              "[time].dt_initial must be > 0 or the string \"auto\" "
+              "(-1 sentinel); got " << cfg.time.dt_initial);
+  MFEM_VERIFY(cfg.time.dt_initial <= cfg.time.dt_max
+              || cfg.time.dt_initial == static_cast<real_t>(-1.0),
+              "[time].dt_initial (" << cfg.time.dt_initial
+              << ") must be <= dt_max (" << cfg.time.dt_max << ")");
```

**Test case:**
```cpp
static void T_13_dt_initial_zero_aborts()
{
   std::string toml = MinimalLSWHeader();
   // Replace `dt_initial = "auto"` with `dt_initial = 0.0`.
   const auto pos = toml.find("dt_initial=\"auto\"");
   toml.replace(pos, std::string("dt_initial=\"auto\"").size(),
                "dt_initial=0.0");
   toml += MinimalLSWBlock();
   TEST_ASSERT(ParseAbortsInChild(toml),
               "dt_initial = 0 must abort (must be > 0 or \"auto\")");
}
```

---

### [R-007] [MODERATE] [spatial_friction.cpp:toml_bool] — `toml_bool` does not verify the underlying type; silent misuse risk

**Category:** ASSUMPTION

**Description:**
Other `toml_*` helpers check the toml value type and abort with a clean message on mismatch:

```cpp
real_t toml_real(const toml::value& tbl, const std::string& key, real_t default_val)
{
   if (!tbl.contains(key)) { return default_val; }
   const auto& v = tbl.at(key);
   if (v.is_floating()) { return static_cast<real_t>(v.as_floating()); }
   if (v.is_integer()) { return static_cast<real_t>(v.as_integer()); }
   MFEM_ABORT("TOML key '" << key << "' must be a number (float or int)");
   return default_val;
}
```

But `toml_bool` calls `as_boolean()` directly with no type guard:

```cpp
bool toml_bool(const toml::value& tbl, const std::string& key, bool default_val)
{
   if (!tbl.contains(key)) { return default_val; }
   return tbl.at(key).as_boolean();
}
```

If a user writes `use_pml = "false"` (a string), toml11's `as_boolean()` throws/aborts with a less-helpful library message instead of the project's standard `MFEM_ABORT("TOML key '...' must be a boolean")`.

**Trigger:**
Any TOML with a quoted boolean (`use_pml = "true"`) or numeric (`use_pml = 1`) value where a bool was expected.

**Actual behavior:**
toml11's terse internal abort surfaces instead of the project's diagnostic message.

**Expected behavior:**
Same shape as `toml_real`/`toml_int`: explicit type check + `MFEM_ABORT`.

**Suggested fix:**
```diff
 bool toml_bool(const toml::value& tbl, const std::string& key, bool default_val)
 {
    if (!tbl.contains(key)) { return default_val; }
-   return tbl.at(key).as_boolean();
+   const auto& v = tbl.at(key);
+   if (v.is_boolean()) { return v.as_boolean(); }
+   MFEM_ABORT("TOML key '" << key << "' must be a boolean (true/false)");
+   return default_val;
 }
```

**Test case:**
```cpp
static void T_14_use_pml_string_aborts()
{
   std::string toml = MinimalLSWHeader();
   // Replace use_pml=false with use_pml="false"
   const auto pos = toml.find("use_pml=false");
   toml.replace(pos, std::string("use_pml=false").size(),
                "use_pml=\"false\"");
   toml += MinimalLSWBlock();
   TEST_ASSERT(ParseAbortsInChild(toml),
               "use_pml=\"false\" (quoted) must abort with a clean message");
}
```

---

### [R-008] [MODERATE] [spatial_friction.cpp:apply_csm_impl] — Stress sidecar load is not protected against the BP2 / non-3-D ctor

**Category:** EDGE_CASE

**Description:**
`ApplyCsmStressSidecar` checks `geom.NumZeroNormalFallbacks() == 0` and calls `geom.ComputeSAFSParams(field, ...)`.  The non-templated `ComputeSAFSParams(const StressField3D&, ...)` overload itself has `MFEM_VERIFY(is_bp5_, ...)` and `FieldProjector::ProjectFaultPreStress` has its own `MFEM_ABORT(... "BP2 / antiplane ctor cannot be used here")`, so a BP2-built FaultGeometry will abort eventually.  However, the order of guards is awkward: `NumZeroNormalFallbacks()` is invoked **before** `is_bp5_` is verified, and that accessor's defined behaviour on a BP2 FaultGeometry is unclear from the header — if it dereferences uninitialised per-DOF arrays, the abort message may be confusing.

This is on the same risk surface as Phase 3's "no silent fallback" rule.  A clearer error chain would be to check `is_bp5_` (or equivalent) first.

**Trigger:**
Pass a BP2 FaultGeometry into `ApplyCsmStressSidecar`.

**Actual behavior:**
Possibly the `NumZeroNormalFallbacks() == 0` check passes (because no per-DOF normals were ever computed), then `ComputeSAFSParams` aborts with its own message.  Or `NumZeroNormalFallbacks()` faults.  Unverified.

**Expected behavior:**
A single, precise abort: "ApplyCsmStressSidecar requires a BP5 (3-D) FaultGeometry — the BP2 ctor cannot be used here."

**Suggested fix:**
```diff
@@ apply_csm_impl @@
   MFEM_VERIFY(spec.kind == StressSourceKind::SidecarHDF5,
               "ApplyCsmStressSidecar: StressSpec.kind must be SidecarHDF5");
   MFEM_VERIFY(!spec.sidecar_path.empty(),
               "ApplyCsmStressSidecar: spec.sidecar_path is empty");
+  MFEM_VERIFY(geom.IsBP5(),
+              "ApplyCsmStressSidecar: requires the 3-D / BP5 FaultGeometry "
+              "ctor.  The BP2 / antiplane ctor cannot be used here.");
   ...
```
(Verify the `IsBP5()` accessor name — `fault/fault_geometry.hpp` should already expose `is_bp5_` via a public getter; if not, add one alongside this fix.)

**Test case:**
```cpp
static void S_4_bp2_geom_aborts()
{
   const bool aborted = RunInChild([]() {
      StressSpec spec;
      spec.kind = StressSourceKind::SidecarHDF5;
      spec.sidecar_path = "<...committed CSM sidecar...>";
      FaultGeometry<mfem::Mesh> geom = MakeBP2FaultGeometry();  // antiplane
      ApplyCsmStressSidecar(spec, geom);  // must abort cleanly
   });
   TEST_ASSERT(aborted, "BP2 geom must abort BEFORE the StressField3D ctor");
}
```

---

### [R-009] [LOW] [test_spatial_friction_resolver.cpp:R_1_rs_defaults] — Test name claims `eta='auto'` but explicitly overrides to `eta_auto=false`

**Category:** QUALITY (test coverage gap; misleading test name)

**Description:**
`R_1_rs_defaults` is titled "[R-1] RS defaults + eta='auto' per-DOF computation", but the body explicitly disables `eta_auto`:

```cpp
RateStateBlock cfg;
...
cfg.eta_auto = true;
...
// Simpler: skip RS eta=auto in this serial test using a small
// overload — we instead exercise eta_auto=false explicitly.
cfg.eta_auto = false;
cfg.eta_default = 5.0e6;
```

Result: the `eta='auto'` per-DOF computation path (`resolve_rs_impl` lines 897–916) is not exercised by any unit test.  The Plan §Phase 1 §Acceptance counted "8 RS tests" assuming this case was covered.

This compounds with [R-003] — the `eta_auto` path silently has the wrong IP, but no test fires.

**Trigger:**
A regression in the `eta_auto` formula would go undetected.

**Suggested fix:**
Either (a) rename the test to reflect what it actually does (`R_1_rs_defaults_eta_explicit`), or (b) implement the eta_auto path on the serial-mesh overload (which already exists in the implementation):

```diff
 static void R_1_rs_defaults()
 {
-   std::cout << "\n[R-1] RS defaults + eta='auto' per-DOF computation\n";
+   std::cout << "\n[R-1a] RS defaults + eta = explicit positive number\n";
    ...
    cfg.eta_auto = false;
    cfg.eta_default = 5.0e6;
    ...
 }
+
+static void R_1b_eta_auto_per_dof()
+{
+   std::cout << "\n[R-1b] RS eta='auto' computes 0.5 * sqrt(mu*rho) per DOF\n";
+   const int N = 2;
+   Vector dofs; Array<int> attr, elem;
+   make_synthetic_dofs(N, 2000.0, dofs, attr, elem);
+   Vector sn_total(N); sn_total = 50e6;
+   RateStateBlock cfg;
+   cfg.a_default = 0.010; cfg.b_default = 0.015;
+   cfg.Dc_default = 0.004; cfg.V_init_default = 1e-9;
+   cfg.f_0_default = 0.6;  cfg.V_0_default = 1e-6;
+   cfg.sigma_n_default = 50e6;
+   cfg.eta_auto = true; cfg.eta_default = 0.0;
+   auto mat = MaterialField::MakeConstant(32e9, 32e9, 2670.0);
+   TinyMeshHolder mh;
+   PorePressureSpec pp;
+   SpatialFrictionResolver R;
+   mfem::Mesh& srl = mh.mesh();
+   auto p = R.ResolveRateState(cfg, dofs, elem, attr, mat, srl, pp, sn_total);
+   const real_t eta_expect = 0.5 * std::sqrt(32e9 * 2670.0);
+   for (int i = 0; i < N; ++i)
+      TEST_NEAR(p.eta(i), eta_expect, 1e-6, "eta = 0.5*sqrt(mu*rho)");
+}
```

**Test case:**
(See the proposed `R_1b_eta_auto_per_dof` above.)

---

### [R-010] [LOW] [fault_geometry_safs_templated.inl:50–53] — Templated overload's empty-N early-out sets `tau_pre_` to size 0; non-templated does NOT — minor behavioural drift

**Category:** QUALITY (consistency)

**Description:**
The new templated `ComputeSAFSParams<StressSource>` early-outs as:

```cpp
if (num_fault_dofs_ == 0)
{
   sigma_n_per_dof_.SetSize(0);
   tau_pre_.SetSize(0);            // ← templated zeroes tau_pre_
   safs_params_computed_ = true;
   return;
}
```

The existing non-templated `ComputeSAFSParams(const StressField3D&, ...)` early-out (`fault_geometry_safs.inl:37–42`) does NOT clear `tau_pre_`:

```cpp
if (num_fault_dofs_ == 0)
{
   sigma_n_per_dof_.SetSize(0);    // ← only sigma_n_per_dof_; tau_pre_ untouched
   safs_params_computed_ = true;
   return;
}
```

This is benign for the SAFS smoke path (`num_fault_dofs_ == 0` is unusual), but it is a behavioural difference between the two overloads on the same input — a violation of the plan's "templated body mirrors the non-templated body" promise (§Phase 3b Detailed Req. 2).  If any caller is sensitive to `tau_pre_.Size()` after a zero-DOF init, the two overloads will produce different results.

**Suggested fix:**
Match the non-templated early-out exactly:
```diff
   if (num_fault_dofs_ == 0)
   {
      sigma_n_per_dof_.SetSize(0);
-     tau_pre_.SetSize(0);
      safs_params_computed_ = true;
      return;
   }
```

(Or, alternatively, change the non-templated path to also `tau_pre_.SetSize(0)` so both overloads pin the same post-condition.  The latter is preferable for predictability but touches the byte-exact gated file.)

---

### [R-011] [LOW] [test_spatial_friction_resolver.cpp:L_3 comment] — In-line comment claims a rule "violates" the inequality but actually satisfies it

**Category:** QUALITY (misleading comment)

**Description:**
In `L_3_last_match_wins`:

```cpp
// (NOTE: rule B's mu_s = 0.9 violates the per-DOF inequality mu_d
//  (0.5) < mu_s (0.9), so this is valid; pick distinct from 1.3
//  so we can check last-match wins.)
```

"0.5 < 0.9" SATISFIES the inequality `mu_d < mu_s`.  Saying it "violates" is wrong; the rest of the sentence ("so this is valid") suggests the author meant the opposite.  Confusing for a future reader debugging the validator.

**Suggested fix:**
```diff
-   // (NOTE: rule B's mu_s = 0.9 violates the per-DOF inequality mu_d
-   //  (0.5) < mu_s (0.9), so this is valid; pick distinct from 1.3
-   //  so we can check last-match wins.)
+   // (NOTE: rule B's mu_s = 0.9 SATISFIES the per-DOF inequality
+   //  mu_d (0.5) < mu_s (0.9), so this is valid; pick distinct from
+   //  rule A's 1.3 so we can check last-match wins.)
```

---

### [R-012] [LOW] [spatial_velocity.cpp:load_impl] — Documented deviation from plan §Phase 2 Detailed Req. 3 ("reuse `FieldProjector::AbortContainmentFailure`") — `AbortContainmentFailure` is private; a duplicate is rolled in spatial_velocity.cpp

**Category:** QUALITY (documented deviation; arguably should bubble the helper to `public:`)

**Description:**
Plan §Phase 2 Detailed Req. 3 says:

> All three fields must pass `ContainsBBox` against the pmesh; reuse `FieldProjector::AbortContainmentFailure` (do not roll a new abort path; matches `heterogeneous_material_plan.md` R-009 round-3).

The implementation rolls a new abort path with the comment:

```cpp
// The shared helper
// `FieldProjector::AbortContainmentFailure` is private; to keep
// io/field_coefficient.hpp byte-exact for TPV/BP5 we inline an
// equivalent precise-message abort here.
```

`AbortContainmentFailure` IS private (`io/field_coefficient.hpp:271 private:`; declaration on `:293`).  But making it `public:` is a single-line change that does NOT affect byte-exactness of any compiled output (it only changes whether external code can name the symbol).  The current duplicate-message approach risks the two error messages drifting apart over time.

**Suggested fix:**
Move `AbortContainmentFailure` (and its sibling `AbortRangeFailure`) from `private:` to `public:` in `io/field_coefficient.hpp`, then call it from `spatial_velocity.cpp`:

```diff
@@ io/field_coefficient.hpp @@
  public:
     ...
-private:
     /// Throw the "mesh not contained" abort with a formatted bbox table.
     static void AbortContainmentFailure(
        const DataField3D& field,
        real_t mxmin, real_t mxmax,
        real_t mymin, real_t mymax,
        real_t mzmin, real_t mzmax);
     /// Throw the "post-projection out of declared bounds" abort.
     static void AbortRangeFailure(
        const std::string& field_name,
        real_t observed_lo, real_t observed_hi,
        real_t declared_lo, real_t declared_hi);
+private:
     ...

@@ spatial/code/spatial_velocity.cpp @@
   auto check = [&](const DataField3D& field)
   {
      const bool inside = field.ContainsBBox(mxmin, mxmax,
                                             mymin, mymax,
                                             mzmin, mzmax);
      if (!inside)
      {
-        const auto& bb = field.BBox();
-        MFEM_ABORT("LoadSpatialVelocityBundle: field '"
-                   << field.FieldName()
-                   << "' bbox does not contain the mesh.\n"
-                   << ...
-                   << "  Either widen the velocity sidecar or shrink "
-                   << "the mesh.");
+        FieldProjector::AbortContainmentFailure(
+           field, mxmin, mxmax, mymin, mymax, mzmin, mzmax);
      }
   };
```

This is LOW because the current code does not produce wrong results — it merely duplicates a message that can drift.

---

## Summary

- Critical issues: **3**  (R-001 missing Phases H/4/5/6, R-002 IP corner-not-centroid in `ResolveForcedRupture`, R-003 same bug in `eta_auto`)
- Moderate issues: **5**  (R-004 missing-block parser silence, R-005 fake stress-bundle tests, R-006 missing dt_initial validator, R-007 toml_bool no type guard, R-008 ApplyCsmStressSidecar BP2-geom guard)
- Low issues: **4**  (R-009 misleading test name, R-010 templated empty-N drift, R-011 misleading comment, R-012 duplicate abort message)
- Plan compliance: **PARTIAL** — Phases 0–3b have implementations; Phases H, 4, 5, 6 do not.  Phase-1 forced-rupture helper exists but no consumer.
- Verdict: **FAIL — must fix R-001/R-002/R-003 before any single-event SAFS dynamic-rupture run.**  R-002 and R-003 will silently bias `T_forced(r)` and `eta` whenever the material is non-constant (i.e., the entire intended SAFS production configuration).  R-001 is informational and acknowledges that the implementation is incomplete relative to the rev-3 plan.

## Unreviewed Areas

- **Phase H heterogeneous WaveOperator code** — not implemented yet, nothing to review.
- **Phase 4 driver source (`drivers/spatial_dyn_driver.cpp`)** — does not exist.
- **Phase 5 free functions (`dynamic/spatial_setup.hpp`)** — does not exist.
- **Phase 5 checkpoint `DRIVER_TAG_V1` extension on `io/tpv104_checkpoint.hpp`** — header is at its pre-extension form (no `driver_tag` default-argument parameters).
- **Phase 5 new `FaultGeometry` ctor overload** — `fault/fault_geometry.hpp` has the templated `ComputeSAFSParams<StressSource>` but no new BP5 ctor taking pre-built per-DOF arrays.
- **Phase 6 sbatches and verification scripts** — none of `jobs/safs/spatial_dyn_*.sbatch`, `verify_spatial_dyn_smoke_safs.py`, or `compare_spatial_dyn_velocity_models.py` exist.
- **`verify_constant_tensor_projection.py`** — the script exists at `spatial/code/scripts/verify_constant_tensor_projection.py` but was not exercised in this review (Python; pillar-level pytest path).
- **TPV/BP5 byte-exact regression gates** — `T-FAULTGEO-LEGACY-CTOR-BYTE-IDENTICAL`, `T-CHECKPOINT-V1-BYTE-IDENTICAL`, `T-CHECKPOINT-DRIVER-TAG-BACK-COMPAT`, `T-PHASEH-SCALAR-PARITY` — none can yet run because Phases 5/H are not implemented.
- **Templated overload's behaviour vs `seas_test_compute_safs_params` (13/13)** — verifying byte-identical output of the existing 13 tests requires `make test-tpv-bp5-byte-exact`; not run as part of this static review.
