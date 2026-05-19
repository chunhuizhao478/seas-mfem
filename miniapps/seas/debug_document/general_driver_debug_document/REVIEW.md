# Code Review: 2026-05-18 — seas_spatial_dyn_driver vs seas_tpv205_driver parity on TPV205

## Review Scope
- Reference driver (ground truth): `miniapps/seas/drivers/tpv205_driver.cpp`
- Subject driver (general): `miniapps/seas/drivers/spatial_dyn_driver.cpp`
- TOML config (subject input): `miniapps/seas/tpv205/configs/tpv205.toml`
- Goal: does (`spatial_dyn_driver.cpp` + `tpv205.toml`) reproduce the **exact setup** that `tpv205_driver.cpp` runs?
- Domain context consulted:
  - `miniapps/seas/CLAUDE.md` (sign conventions, basis convention)
  - `miniapps/seas/config/tpv205_params.hpp` (TPV205 spec constants)
  - `miniapps/seas/dynamic/tpv205_setup.hpp` (native init path)
  - `miniapps/seas/dynamic/tpv205_friction.hpp` (`LSWFrictionCoefficient_TPV205` + barrier short-circuit)
  - `miniapps/seas/dynamic/wave_operator.inl:339-368` (hard-coded internal `FaultBasis` with `ref_normal = (0,-1,0)`)
  - `miniapps/seas/fault/fault_basis.hpp:373-447` (`ComputeOrientedFrame`: `strike = up × n_ref`, `dip = strike × n_ref`)
  - `miniapps/seas/spatial/code/spatial_friction.hpp/cpp` (`StressPatch::contains`, `SpatialRule::matches`, `ResolveSlipWeakening`)
  - `miniapps/seas/spatial/code/spatial_stress.cpp:110-156` (`ConstantTensorWithPatchesStressSource::Evaluate`)
  - `miniapps/seas/fault/fault_geometry_safs_templated.inl:82-118` (`ComputeSAFSParams<StressSource>`)
  - `miniapps/seas/dynamic/spatial_setup.hpp` (`InitializeFaultDOFs_Spatial`)
  - `miniapps/seas/tpv205/mesh/tpv2053d_200m.geo` (mesh coordinate convention)

The plan documents pointed to by the spatial driver are:
- `safs/project_7.0_alternative/document/spatial_dynamic_rupture_plan.md` (Phase 4)
- `safs/project_7.0_alternative/document/PLAN_phase_R_exact_bimaterial_riemann_rev3.md` (R.4 step 1 — declares this very TOML the canonical "TPV205 driven by the spatial driver")

**Conclusion up front (so the fix agent does not have to scroll):** the combination
`spatial_dyn_driver --config tpv205/configs/tpv205.toml` does **NOT** reproduce the
native TPV205 setup.  Three independent **CRITICAL** defects each prevent
spontaneous rupture (or, if any one were removed in isolation, would still
produce wrong physics).  Multiple **MODERATE** dispatch / numerics deviations
also stand.

## Findings

---

### [R-001] CRITICAL — `tpv205.toml` patches never match the mesh because `center_z_m` uses spec depth (positive) while the mesh and `StressPatch::contains` use raw mesh `z` (negative below the surface)

**Category:** BUG

**Description:**
`StressPatch::contains(x, y, z)` (`spatial/code/spatial_friction.hpp:266-271`) compares
the raw mesh-coordinate `z` against `center_z_m` via `|z - center_z_m| <= half_z_m`.
The fault DOF coordinates fed to it (`fault_geometry_safs_templated.inl:82-94`)
are the raw 3-D physical coords returned by `ftr->Face->Transform(ip, phys)`.

For the TPV205 mesh (`tpv2053d_200m.geo`, line 84:
`Point(200) = {X_nucl, Width_nucl*Cos(Fault_dip), -Width_nucl*Sin(Fault_dip), lc_nucl}`
with `Fault_dip = 90°`), the nucleation hypocenter at SPEC depth 7.5 km is at
**mesh z = -7500**, not +7500.  Every fault DOF below the free surface has `z ≤ 0`.

The TOML places the three patches at `center_z_m = 7500.0`:
```toml
[[stress.patch]]
center_x_m  = 0.0;  center_y_m  = 0.0;  center_z_m  = 7500.0
half_x_m    = 1500.0;  half_z_m  = 1500.0
sigma_xy_pa = 81.6e6
```
For the hypocenter DOF at z = −7500: `|−7500 − 7500| = 15 000 > half_z_m = 1500` ⇒ `contains` returns FALSE.  The other two patches (left/right) likewise miss.

Consequence: `tau_pre_strike` is the uniform 70 MPa background everywhere on the fault.  Yield is `μ_s · σ_n = 0.677 · 120 MPa = 81.24 MPa`.  Since `70 < 81.24`, **no DOF ever reaches yield and rupture never initiates** — the simulation runs to `tfinal = 12 s` producing essentially zero slip.  The native `tpv205_driver` (which converts to spec depth via `std::abs(fault_coords[i](2))` in `InitializeFaultDOFs_TPV205`) finds the nucleation patch at 81.6 MPa > 81.24 MPa and ruptures spontaneously at `t = 0+`.

**Trigger:**
Run `seas_spatial_dyn_driver --config miniapps/seas/tpv205/configs/tpv205.toml` on any of the existing TPV205 meshes (`tpv2053d_*.msh`).

**Actual behavior:**
All 3 patches miss every fault DOF; the strike pre-stress is the uniform 70 MPa background; rupture never initiates.

**Expected behavior:**
The native TPV205 driver places the nucleation patch at the down-dip-depth-of-7.5 km DOFs (mesh z = −7500) and seeds `tau2_0 = 81.6 MPa` there.  Rupture initiates at t = 0+.

**Suggested fix:**
There are two viable approaches; the second is recommended because it keeps the
TOML readable and matches the wording of the spec.

Option A (TOML-only — keeps the spatial driver's "raw mesh coords" contract): flip
the sign on each `center_z_m`:
```diff
 [[stress.patch]]
 center_x_m  = 0.0
 center_y_m  = 0.0
-center_z_m  = 7500.0
+center_z_m  = -7500.0   # mesh z is negative below the surface; spec depth = +7.5 km
 half_x_m    = 1500.0
 half_z_m    = 1500.0
 sigma_xy_pa = 81.6e6
```
(repeat for the left and right patches; also flip the bottom-barrier `z_min_m` —
see R-002 below.)

Option B (driver-level — turn the spatial driver into a depth-aware match):
introduce an explicit `depth_axis` setting in `[fault_geometry]` (or `[stress]`)
and have `ComputeSAFSParams<StressSource>` rewrite the per-DOF coord triple
into `(along_strike, fault_perp, depth)` before invoking `source.Evaluate`,
matching the SPEC convention used throughout the TOML.  This is the same
contract the native TPV205 driver implements with `down_dip = std::abs(z)` at
`dynamic/tpv205_setup.hpp:105`.  Concretely, in
`fault_geometry_safs_templated.inl`:
```diff
-      const real_t x = coords(3 * i + 0);
-      const real_t y = coords(3 * i + 1);
-      const real_t z = coords(3 * i + 2);
+      const real_t x = coords(3 * i + 0);
+      const real_t y = coords(3 * i + 1);
+      const real_t z_mesh = coords(3 * i + 2);
+      // The stress-source consumer (StressPatch, DepthProportional…) uses the
+      // SPEC depth convention (positive downward).  Map raw mesh z (z=0 surface,
+      // z<0 below) to spec depth.
+      const real_t z = std::abs(z_mesh);
```
plus a corresponding flip in `SpatialRule::matches` so `z_*_m` bounds also use
spec depth.

**Test case:**
```cpp
// File: tests/unit/test_spatial_dyn_tpv205_patches.cpp
TEST(SpatialDynDriverTpv205Patches, R001_nucleation_patch_seeds_8160_kPa_at_hypocenter)
{
   // Load the canonical TOML and synthesize ONE fault DOF at the hypocenter
   // (mesh coords (0, 0, -7500)).  ConstantTensorWithPatchesStressSource
   // Evaluate(0, 0, -7500) MUST return 81.6 MPa in σ_xy (the nucleation patch).
   SpatialFrictionConfig cfg = ParseSpatialFrictionConfigString(
      ReadFile("miniapps/seas/tpv205/configs/tpv205.toml"));

   ConstantTensorWithPatchesStressSource src(
      cfg.stress.sigma_xx_pa, cfg.stress.sigma_yy_pa, cfg.stress.sigma_zz_pa,
      cfg.stress.sigma_xy_pa, cfg.stress.sigma_yz_pa, cfg.stress.sigma_xz_pa,
      cfg.stress.patches);

   const DenseMatrix S = src.Evaluate(0.0, 0.0, -7500.0);   // mesh z below surface
   EXPECT_NEAR(S(0, 1), 81.6e6, 1.0)    // σ_xy expected to be the patch value
       << "Nucleation patch missed the hypocenter DOF; rupture cannot initiate "
          "(see REVIEW.md R-001).";
}
```

---

### [R-002] CRITICAL — Bottom strength-barrier rule (`z_min_m = 15000`) never fires because the mesh-z sign convention is reversed; depths > 15 km below the rupture area run with `μ_s = 0.677` instead of the locked sentinel

**Category:** BUG

**Description:**
`SpatialRule::matches` (`spatial_friction.cpp:84-103`) for `kind = "barrier"`
checks `(z >= z_min_m) && (z <= z_max_m)` against raw mesh `z`.  The TOML
defines the bottom barrier as
```toml
[[friction.slip_weakening.spatial]]
kind     = "barrier"
z_min_m  = 15000.0
```
intending "anywhere deeper than 15 km below the surface".  Because mesh z is
negative below the surface (mesh extends from `z = 0` down to `z = Zmin = -60e3`
per `tpv2053d_200m.geo:56`), the condition `z >= 15000` is **never satisfied for
any fault DOF on this mesh**.

The native TPV205 driver applies the bottom barrier via
`ComputeMuS_TPV205(along_strike, down_dip)` (`config/tpv205_params.hpp:189-196`)
with `down_dip = std::abs(z)`, so DOFs at mesh z = −15500 (depth 15.5 km) get
`μ_s = 10 000` and stay locked.

Consequence: even if R-001 is fixed and rupture initiates, the spatial driver
allows the rupture front to escape **downward** past 15 km depth into a region
the spec (SCEC TPV5 §11) requires to be a strength barrier.  Other than the
free surface and the ±15 km along-strike barriers (which DO fire correctly
because the `x_max_m = -15000` / `x_min_m = 15000` rules use the same sign as
the mesh x-axis), the deep barrier is missing.

**Trigger:**
A fault DOF at mesh z = −15500 (depth 15.5 km) on `tpv2053d_200m.msh` should
resolve to `μ_s = 1.0e6` (barrier sentinel) but instead resolves to
`μ_s = 0.677` (default).

**Actual behavior:**
`SpatialFrictionResolver::ResolveSlipWeakening` returns `lsw.mu_s(i) = 0.677`
for the deep DOF.

**Expected behavior:**
Native TPV205 returns `10 000` (the spec barrier sentinel) for depths
beyond 15 km.  Encoded in the spatial schema this means `lsw.mu_s(i) = 1.0e6`
(the spatial resolver's internal barrier sentinel — see `spatial_friction.cpp:1527`).

**Suggested fix:**
Same options as R-001; pick the one that matches your chosen fix:

Option A (TOML-only):
```diff
 [[friction.slip_weakening.spatial]]
 # Bottom barrier: depth (z in this driver's frame) > 15 km.
 kind     = "barrier"
-z_min_m  = 15000.0
+z_max_m  = -15000.0   # mesh z is negative below surface; depth 15 km = z = -15000
```
This fires when raw `z ≤ -15000` (i.e., depth ≥ 15 km), which is the intended
semantics.

Option B (driver-level depth mapping — see R-001 Option B): leave the TOML
exactly as-written, and have `SpatialRule::matches` evaluate against spec
depth, not raw `z`.  Same edit as R-001.

**Test case:**
```cpp
TEST(SpatialDynDriverTpv205Barrier, R002_deep_dof_gets_locked)
{
   SpatialFrictionConfig cfg = ParseSpatialFrictionConfigString(
      ReadFile("miniapps/seas/tpv205/configs/tpv205.toml"));

   Vector coords(3);   coords(0) = 0.0;  coords(1) = 0.0;  coords(2) = -15500.0;
   Array<int> attr(1);  attr[0] = cfg.boundary.fault_attr;

   SpatialFrictionResolver r;
   auto p = r.ResolveSlipWeakening(*cfg.slip_weakening, coords, attr);

   EXPECT_GT(p.mu_s(0), 5000.0)
       << "Deep DOF (mesh z=-15500, depth 15.5 km) must be locked by the "
          "bottom strength barrier (see REVIEW.md R-002).";
}
```

---

### [R-003] CRITICAL — Strike-direction pre-stress (`tau2_0`) is sign-flipped because the spatial driver builds its external `FaultBasis` with `ref_normal = (0,+1,0)` while the wave operator's internal `FaultBasis` is hard-coded to `ref_normal = (0,-1,0)`

**Category:** BUG

**Description:**
The wave operator (`dynamic/wave_operator.inl:349`) unconditionally builds
its internal `FaultBasis` with
```cpp
Vector ref_normal(3); ref_normal = 0.0; ref_normal(1) = -1.0;
Vector up(3);         up = 0.0;         up(2) = 1.0;
```
This is the basis used by `FaultFaceFlux::ComputeTrialTraction` and by the
LSW dispatch — the components `tau1_trial`, `tau2_trial` it computes are
already in the canonical fault-local frame with
**`tangent2 = (+1, 0, 0)`** (per `ComputeOrientedFrame`:
`strike = up × n_ref = (0,0,1) × (0,-1,0) = (+1,0,0)`).

The spatial driver, however, **builds a SEPARATE `FaultBasis` externally**
to feed `FaultGeometry::ComputeSAFSParams<StressSource>`:
```cpp
// spatial_dyn_driver.cpp:1007-1019
Vector ref_normal(3);
ref_normal(0) = cfg.fault_geometry.ref_normal[0];   // = 0.0
ref_normal(1) = cfg.fault_geometry.ref_normal[1];   // = +1.0  ← from TOML
ref_normal(2) = cfg.fault_geometry.ref_normal[2];   // = 0.0
…
fbasis.Compute(pmesh, fault_int_faces, ref_normal, up_vec);
```
With `ref_normal = (0,+1,0)`, `ComputeOrientedFrame` yields
**`tangent2 = up × n_ref = (0,0,1) × (0,+1,0) = (-1,0,0)`** — the
opposite-sign strike vector.

The projection inside `ComputeSAFSParams` then computes
`tau2_external = t2_external · S · n_external = (-1,0,0) · (S · (0,+1,0))
              = -S_xy = -70 MPa` for the TPV205 background.

This value flows verbatim into `DOFData::tau2_0` via
`InitializeFaultDOFs_Spatial` → `seed_static_dof_fields`
(`dynamic/spatial_setup.hpp:84-86`).  At runtime the LSW solver
(`SolveLSW_TPV205`) adds the externally-frame `tau2_0` to the
**internally-frame** `tau2_trial`:
```cpp
tau2_total = tau2_0 + tau2_nuc + tau2_trial
//           ^^^^^^^^                       (external frame; t2 = -x)
//                                ^^^^^^^^^ (internal frame; t2 = +x)
```
These two terms live in opposite-signed strike axes.  The result is
**garbage**: in the trivial t = 0 case where `tau2_trial = 0` and `tau2_nuc = 0`,
`|tau2_total| = 70 MPa` (because only one term contributes), but the slip-rate
direction the solver then writes (`V2 = V_abs · tau2_total / |tau_abs|`) points
in the **−x** strike direction, while every downstream consumer
(`Q_imp` assembly, `WriteBackState`) is in the +x internal frame.  Even if
patch nucleation were active, the rupture would propagate in the wrong
direction.

**Trigger:**
Any fault DOF with non-zero strike pre-stress and the TOML's
`fault_geometry.ref_normal = [0.0, 1.0, 0.0]`.

**Actual behavior:**
External `t2 = (-1, 0, 0)` ⇒ `tau2_0 = -70 MPa` (background) seeded into
`DOFData`; rupture sense is reversed when combined with the wave operator's
internal `t2 = (+1, 0, 0)` trial traction.

**Expected behavior:**
Native TPV205 sets `tau2_0 = +70 MPa` (`InitializeFaultDOFs_TPV205` →
`ComputeTau2_0_TPV205` returns positive values; component 2 in BP5
canonical = strike = +x for the wave operator's internal basis).

**Suggested fix:**
Flip the TOML so the external `FaultBasis` matches the wave operator's
hard-coded internal one:
```diff
 [fault_geometry]
-# TPV205 fault is the y = 0 plane; rupture propagates in x (along strike)
-# and -z (down dip).  ref_normal = +y, up = +z.  kind tag is informational.
-ref_normal = [0.0, 1.0, 0.0]
+# TPV205 fault is the y = 0 plane; rupture propagates in x (along strike)
+# and -z (down dip).  ref_normal MUST match the wave operator's hard-coded
+# internal `(0, -1, 0)` (`dynamic/wave_operator.inl:349`); otherwise the
+# external Cauchy projection (used for tau_pre) lives in a different strike
+# basis than the internal trial-traction frame (REVIEW R-003).
+ref_normal = [0.0, -1.0, 0.0]
 up         = [0.0, 0.0, 1.0]
 kind       = "tpv205_lsw"
```

A more durable fix: stop building the external `FaultBasis` from TOML at all
and instead have the spatial driver pull `ref_normal` / `up` from
`wave.GetFaultBasis()` itself (or pass the wave operator's `FaultBasis*` into
`ComputeSAFSParams`).  That removes the entire foot-gun.  Concrete change in
`spatial_dyn_driver.cpp`:
```diff
-   Vector ref_normal(3);
-   ref_normal(0) = cfg.fault_geometry.ref_normal[0];
-   ref_normal(1) = cfg.fault_geometry.ref_normal[1];
-   ref_normal(2) = cfg.fault_geometry.ref_normal[2];
-   Vector up_vec(3);
-   up_vec(0) = cfg.fault_geometry.up[0];
-   up_vec(1) = cfg.fault_geometry.up[1];
-   up_vec(2) = cfg.fault_geometry.up[2];
-
-   FaultBasis fbasis;
-   fbasis.Compute(pmesh, fault_int_faces, ref_normal, up_vec);
+   // Reuse the wave operator's internal FaultBasis to guarantee tau_pre and
+   // trial traction live in the same fault-local frame (REVIEW R-003).
+   const FaultBasis &fbasis = *wave.GetFaultBasis();
+   // (Validate the TOML's ref_normal / up against the wave operator's choice;
+   // abort if the user supplies an incompatible combination.)
```

If you keep the TOML-controlled `ref_normal` API, add a guard at startup:
```cpp
MFEM_VERIFY(cfg.fault_geometry.ref_normal[1] < 0.0,
            "spatial_dyn_driver: cfg.fault_geometry.ref_normal must be "
            "(0, -1, 0) (the wave operator's internal Tandem convention); "
            "got (" << cfg.fault_geometry.ref_normal[0] << ", "
                    << cfg.fault_geometry.ref_normal[1] << ", "
                    << cfg.fault_geometry.ref_normal[2]
                    << ").  See REVIEW.md R-003.");
```

**Test case:**
```cpp
TEST(SpatialDynDriverTpv205Basis,
     R003_strike_pre_stress_matches_internal_wave_basis)
{
   // Construct a minimal 1-element 3-D fault setup with σ_xy = 70 MPa
   // background.  Run the spatial-driver pipeline through ComputeSAFSParams
   // and assert tau2_0 (component 2 = strike in the BP5 canonical frame
   // used by the wave operator) is +70 MPa, NOT -70 MPa.
   …
   EXPECT_NEAR(geom.GetTauPre()(2 * dof + 1), +70.0e6, 1.0)
       << "Strike pre-stress sign disagrees with the wave operator's "
          "internal FaultBasis (REVIEW R-003).";
}
```

---

### [R-004] MODERATE — CFL safety factor is missing: spatial driver passes the raw user-supplied `cfl = 0.5` to `wave.ComputeMaxDt`, while the native TPV205 driver pre-scales it by `1 / (3·(2·order+1))` (≈ 9× tighter for `order = 1`)

**Category:** DEVIATION

**Description:**
Native `tpv205_driver.cpp:1242-1244`:
```cpp
real_t cfl    = cfl_factor / (3.0 * (2.0 * order + 1.0));
real_t dt_cfl = wave.ComputeMaxDt(cfl);
```
For `order = 1`, `cfl = 0.5 / 9 ≈ 0.05556`.

Spatial `drivers/spatial_dyn_driver.cpp:1255`:
```cpp
const real_t dt_cfl = wave.ComputeMaxDt(cfg.numerics.cfl);
```
With the TOML `cfl = 0.5`, the spatial driver requests a `dt` ~9× larger than
the native driver does for the same mesh.  `ComputeMaxDt` returns
`cfl_mixed_flux_factor · cfl · h_min / cp`, so the resulting `dt` is
proportionally inflated.  Even when `dt_max = 0.01s` (set in the TOML) caps
the result, the early macro-step `dt` is larger than the native driver runs,
and the ADER predictor / corrector dispatch differs accordingly.

This is not necessarily unstable, but it **guarantees** the spatial
driver's per-step time integration is not bit-equivalent to the native
driver's — the entire question "does X reproduce Y?" gets a NO from this
alone, even with R-001 / R-002 / R-003 fixed.

**Trigger:**
Any TPV205 run.

**Actual behavior:**
Spatial driver's `dt_cfl ≈ 9 × native dt_cfl` (until `dt_max` clamp triggers).

**Expected behavior:**
Same `dt_cfl` as the native driver — apply the `1 / (3·(2·order+1))` DG safety
factor before calling `ComputeMaxDt`.

**Suggested fix:**
```diff
-   const real_t dt_cfl = wave.ComputeMaxDt(cfg.numerics.cfl);
+   // DG CFL safety factor (matches the convention of every TPV* native driver:
+   // dt_cfl = cfl · h / (cp · 3·(2p+1))).  Pre-scaling here keeps the
+   // user-facing `[numerics].cfl = 0.5` semantically equal to the native
+   // driver's `--cfl 0.5`.
+   const real_t cfl_scaled = cfg.numerics.cfl
+                            / (3.0 * (2.0 * cfg.mesh.order + 1.0));
+   const real_t dt_cfl = wave.ComputeMaxDt(cfl_scaled);
```
(Or: change the TOML default of `[numerics].cfl` to `0.5 / (3·(2p+1))` and
document the convention.  Pre-scaling in the driver is cleaner because the
spec / native driver users think in `0.5`.)

**Test case:**
```cpp
TEST(SpatialDynDriverTpv205Cfl, R004_dt_cfl_matches_native_driver)
{
   // Build a TPV205-style WaveOperator on a fixed mesh and assert the
   // spatial driver's dt_cfl equals tpv205_driver's dt_cfl to within
   // 0.1%.
   const real_t native_dt = wave.ComputeMaxDt(0.5 / (3.0 * (2*1 + 1)));
   const real_t spatial_dt = wave.ComputeMaxDt(0.5 / (3.0 * (2*1 + 1)));
   EXPECT_NEAR(spatial_dt, native_dt, 1e-3 * native_dt);
}
```

---

### [R-005] MODERATE — Spatial driver hard-codes the sub-step iterator to `O = ader_order` (= 2 for TPV205), while the native TPV205 driver defaults to `--fault-iterator one-shot` (`O = 1`)

**Category:** DEVIATION

**Description:**
Native (`tpv205_driver.cpp:1721-1737`):
```cpp
const bool substep_quadrature =
   (GetDispatchedIterator(fault_iterator) == DispatchedIterator::SubStep);
if (substep_quadrature) {
   const int O = std::max(1, ader_order);  // O = 2
   …
} else {
   // O = 1: single sub-step per macro-step (one-shot semantics).
   std::vector<real_t> deltaT(1, dt);
   std::vector<real_t> weights(1, 1.0);
   substep_iterator.SetSubSteps(deltaT, weights);
}
```
The default of `--fault-iterator` is `"one-shot"`, so a default invocation of
`seas_tpv205_driver` uses `O = 1`.

Spatial (`spatial_dyn_driver.cpp:1722-1726`):
```cpp
const int O = std::max(1, cfg.numerics.ader_order);   // = 2 for TPV205
std::vector<real_t> deltaT(O, dt / static_cast<real_t>(O));
std::vector<real_t> weights(O, 1.0 / static_cast<real_t>(O));
substep_iterator.SetSubSteps(deltaT, weights);
```
The spatial driver **always** uses `O = ader_order`, with no opt-out.

These two paths are not bit-equivalent (substep iterator runs the LSW
closed-form solve twice per macro-step instead of once; the ADER predictor
is sampled at sub-step midpoints rather than at the full-step midpoint).
The contract `T_TPV205_SSI_3` (`AdvanceADERWithSubStep:336-377` comment block)
only guarantees byte-equivalence when both drivers run `O = 1`.

**Trigger:**
Any TPV205 run with the spatial driver's default settings.

**Actual behavior:**
Spatial driver dispatches `Tpv205SubStepIterator` with `O = 2` sub-steps.

**Expected behavior:**
Match the native driver's default `O = 1` (or expose a CLI / TOML knob and
default to one-shot for TPV205 parity).

**Suggested fix:**
Add a `[numerics] fault_iterator = "one-shot"` TOML key (with values
`one-shot | substep`, default `one-shot`) and route it the same way as the
native driver:
```diff
-      const int O = std::max(1, cfg.numerics.ader_order);
-      std::vector<real_t> deltaT(O, dt / static_cast<real_t>(O));
-      std::vector<real_t> weights(O, 1.0 / static_cast<real_t>(O));
-      substep_iterator.SetSubSteps(deltaT, weights);
+      // Match the native TPV/BP5 default: one-shot iterator (O = 1).  The
+      // TOML can opt-in to the SubStep iterator (O = ader_order) when a
+      // problem genuinely requires the per-sub-step quadrature.
+      const bool substep_quadrature =
+         (cfg.numerics.fault_iterator == "substep");
+      if (substep_quadrature)
+      {
+         const int O = std::max(1, cfg.numerics.ader_order);
+         std::vector<real_t> deltaT(O, dt / static_cast<real_t>(O));
+         std::vector<real_t> weights(O, 1.0 / static_cast<real_t>(O));
+         substep_iterator.SetSubSteps(deltaT, weights);
+      }
+      else
+      {
+         std::vector<real_t> deltaT(1, dt);
+         std::vector<real_t> weights(1, 1.0);
+         substep_iterator.SetSubSteps(deltaT, weights);
+      }
```

**Test case:**
```cpp
TEST(SpatialDynDriverTpv205Iterator,
     R005_default_uses_one_shot_iterator)
{
   SpatialFrictionConfig cfg = ParseSpatialFrictionConfigString(
      ReadFile("miniapps/seas/tpv205/configs/tpv205.toml"));
   // After driver setup, assert NumSubSteps() == 1.
   EXPECT_EQ(substep_iterator.NumSubSteps(), 1);
}
```

---

### [R-006] MODERATE — Spatial driver routes `WaveOperator` through the heterogeneous (`MaterialField`, `BoundaryConfig`) ctor instead of the scalar `(λ, μ, ρ)` ctor used by the native TPV205 driver; bit-equivalence is only guaranteed to ~1e-12 relative tolerance, not bit-for-bit

**Category:** DEVIATION

**Description:**
Native (`tpv205_driver.cpp:1224-1226`):
```cpp
WaveOperator<MeshT> wave(pmesh, order,
                         TPV205Params::lambda, TPV205Params::mu,
                         TPV205Params::rho, bc);
```
Spatial (`spatial_dyn_driver.cpp:921`):
```cpp
WaveOperator<ParMesh> wave(pmesh, cfg.mesh.order, material, bc);
```
The latter unconditionally routes interior faces through
`BimaterialFlux::ApplyPerFaceFlux` (R.4 — see `CLAUDE.md` "Bi-material Riemann
(Phase R)").  The driver comments at lines 906-919 explicitly note this is
"byte-equivalent for Mode::Constant input … to within LU rounding (~1e-12
relative, verified by T-PHASEH-SCALAR-PARITY)" — i.e., **not bit-identical**,
only "close".

For the user's stated goal ("does the toml + spatial driver produce the exact
setup"), the **answer is no**: the spatial driver is by construction
~1e-12-relative-different from the native driver on every interior face, by
design.

**Trigger:**
Comparing per-DOF outputs of `seas_spatial_dyn_driver` vs
`seas_tpv205_driver` after any time step.

**Actual behavior:**
Differences accumulate from the very first ADER step; not bit-identical.

**Expected behavior:**
If the user wants the EXACT same setup, the spatial driver must call the
scalar ctor when `cfg.material.kind == Constant` (the documented "the simple
choice" path, but the code intentionally went the other way per the comment
at 906-919).

**Suggested fix:**
Add a CLI / TOML opt-out — `[numerics] interior_flux = "scalar" | "bimaterial"`,
default `"bimaterial"` to keep current behaviour, but document the
`"scalar"` opt-in for byte parity with native TPV205:
```diff
   // -----------------------------------------------------------------
   // 7.  Construct WaveOperator …
   // -----------------------------------------------------------------
-   WaveOperator<ParMesh> wave(pmesh, cfg.mesh.order, material, bc);
+   std::unique_ptr<WaveOperator<ParMesh>> wave_ptr;
+   if (cfg.numerics.interior_flux == "scalar"
+       && material.mode == MaterialField::Mode::Constant)
+   {
+      wave_ptr = std::make_unique<WaveOperator<ParMesh>>(
+         pmesh, cfg.mesh.order,
+         material.lambda_const, material.mu_const, material.rho_const, bc);
+   }
+   else
+   {
+      wave_ptr = std::make_unique<WaveOperator<ParMesh>>(
+         pmesh, cfg.mesh.order, material, bc);
+   }
+   WaveOperator<ParMesh> &wave = *wave_ptr;
```

**Test case:**
Cross-driver byte-equivalence on a fixed TPV205 mesh, `O = 1` iterator,
identical CFL.  Compare ParaView frame 0 (well-defined since
`tfinal = 0` is supported by both drivers).

```cpp
TEST(SpatialDynDriverTpv205Byteparity, R006_scalar_flux_matches_native)
{
   // After R-001, R-002, R-003, R-004, R-005 fixed AND R-006 scalar
   // opt-in, assert per-DOF DOFData::tau2_corr is bit-identical between
   // the two drivers.
}
```

---

### [R-007] MODERATE — Spatial driver call order `SetMixedFluxMode` BEFORE `SetFaultFlux` / `SetFaultDOFData` / `SetAbsorbingBackground` violates the documented R-1205 setter ordering

**Category:** DEVIATION

**Description:**
Per `CLAUDE.md` "Round-11 Mixed-Flux dispatch wiring (R-1205 required call
order)" and the native driver's wiring at `tpv205_driver.cpp:1606-1633`:
```
1. WaveOperator ctor
2. wave.SetFaultFlux
3. wave.SetFaultDOFData
4. wave.SetAbsorbingBackground
5. wave.SetMixedFluxMode    <-- HERE; setter cross-checks invariants
```
The spatial driver calls `SetMixedFluxMode` at line 979, BEFORE
`SetFaultFlux` (line 1211), `SetFaultDOFData` (line 1243), and
`SetAbsorbingBackground` (lines 1247-1249).

For `mixed_flux = "none"` (the TPV205 TOML default) this is silent — the
setter's invariants are trivially satisfied for `MixedFluxMode::None`.  But
any future TOML that sets a non-trivial mixed_flux would hit
`SetMixedFluxMode`'s internal `bc_.fault_attr > 0` cross-check before the
fault flux is actually wired up, producing a misleading abort.

**Trigger:**
TOML with `[numerics] mixed_flux = "adjacent"` or `"all_continuous"`.

**Actual behavior:**
`SetMixedFluxMode` runs against an under-initialised wave operator; aborts
with a spurious message about missing fault wiring.

**Expected behavior:**
Match R-1205 ordering as documented and as the native driver implements it.

**Suggested fix:**
Move `wave.SetMixedFluxMode(...)` from its current position (right after the
ctor at line 979) down to AFTER the full fault wiring is in place — i.e.,
between line 1249 (after `SetAbsorbingBackground`) and the CFL computation
at line 1255.  The simplest delta:
```diff
@@ spatial_dyn_driver.cpp:978-979
-   wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW);
-   wave.SetMixedFluxMode(ParseMixedFlux(cfg.numerics.mixed_flux));
+   wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW);
@@ spatial_dyn_driver.cpp:1249 (after SetAbsorbingBackground)
       wave.SetAbsorbingBackground(Q_bg);
    }
+
+   // R-1205 ordering: SetMixedFluxMode AFTER ctor / SetFaultFlux /
+   // SetFaultDOFData / SetAbsorbingBackground (REVIEW R-007).
+   wave.SetMixedFluxMode(ParseMixedFlux(cfg.numerics.mixed_flux));
```

**Test case:**
```cpp
TEST(SpatialDynDriverWiring, R007_mixed_flux_ordering_R1205)
{
   // Build the spatial driver with [numerics] mixed_flux = "adjacent".
   // The current ordering aborts at SetMixedFluxMode because bc.fault_attr
   // hasn't been registered against the fault flux yet.  After the fix the
   // setup should complete.
}
```

---

### [R-008] MODERATE — Hypocenter z-coordinate in TOML (`z = 7500.0`) uses spec depth (positive) but the rest of the spatial-driver coordinate plumbing uses raw mesh z (negative); the field is "documentation only" today but will silently produce wrong nucleation centres the moment any `[nucleation]` block is added

**Category:** ASSUMPTION

**Description:**
TOML at lines 67-74:
```toml
[hypocenter]
x                   = 0.0
y                   = 0.0
z                   = 7500.0
nucleation_radius_m = 1500.0
nucleation_taper_m  = 0.0
```
The comment at lines 75-78 notes "Documentation only — the actual nucleation
is encoded in the `[[stress.patch]]` array below.  These values are unused
unless a `[nucleation]` block is added".  This is true today but is a
land-mine: a future contributor adding e.g. a `[nucleation]` block of
`instantaneous_overstress_circular` kind would inherit the wrong z sign.
The resolver at `spatial_nucleation.cpp` would place the nucleation centre
at mesh z = +7500 (15 km ABOVE the surface — outside the mesh) instead of
at the actual SPEC depth 7.5 km.

This is the SAME root-cause as R-001 / R-002.  Pulling it out separately
because the fix surface is the schema layer (the `HypocenterSpec` parser),
not the stress source.

**Trigger:**
Add ANY `[nucleation]` block to `tpv205.toml` that consumes `[hypocenter] z`.

**Actual behavior:**
Nucleation centre placed at +z (above surface).

**Expected behavior:**
Hypocenter z is interpreted as spec depth and the spatial driver internally
maps to mesh z = −z (or the user supplies -7500 explicitly with documentation).

**Suggested fix:**
Same as R-001 Option B — make the spatial driver's depth handling explicit
and consistent.  Or, minimally, add a runtime assertion at `LoadSpatialFrictionConfig`
that aborts if the hypocenter z is "in the wrong half-space" relative to
the configured `fault_geometry.up`:
```cpp
MFEM_VERIFY(cfg.hypocenter.z * cfg.fault_geometry.up[2] <= 0.0,
            "[hypocenter].z must be negative when [fault_geometry].up has "
            "positive z (mesh-z convention: z < 0 is below surface).  "
            "Set z = -<depth_m>.  See REVIEW.md R-008.");
```

**Test case:**
```cpp
TEST(SpatialDynDriverHypocenter, R008_positive_z_rejected_with_positive_up)
{
   const std::string toml = R"(
[meta] schema_version = 1; law = "slip_weakening"; description = "test"
[problem] tag = "tpv205"
[boundary] fault_attr = 103
[fault_geometry] up = [0.0, 0.0, 1.0]; ref_normal = [0.0, -1.0, 0.0]
[hypocenter] x = 0.0; y = 0.0; z = +7500.0   # WRONG SIGN
…
)";
   EXPECT_DEATH(ParseSpatialFrictionConfigString(toml),
                "hypocenter.z must be negative");
}
```

---

### [R-009] LOW — Spatial driver writes ZERO station-trace files (no `TPV205StationWriter`), while the native TPV205 driver writes 16 SCEC-compatible `.dat` per-station traces

**Category:** DEVIATION

**Description:**
Native (`tpv205_driver.cpp:1851-1875`):
```cpp
auto stations = DefaultStations_TPV205();
TPV205StationWriter station_writer;
station_writer.Open(output_dir, output_prefix, stations,
                    fault_coords, num_fault_local, comm);
station_writer.WriteStep(0.0, dof_data);
…
auto surface_stations = DefaultSurfaceStations_TPV205();
TPV205SurfaceStationWriter surface_writer;
surface_writer.Open(…); …
```
16 station files + (potential) surface stations are written, each one
formatted to match SCEC TPV5 benchmark-trace column convention.

Spatial driver: only ParaView output is wired (lines 1351-1539); no
station writer.  Comparing the two drivers on SCEC benchmark
diagnostics (rupture-front arrival times, peak slip-rate vs. time at
fixed receivers) is therefore impossible from the spatial driver's
output alone.

**Trigger:**
Running the spatial driver on TPV205 — the canonical SCEC trace files
are simply absent.

**Actual behavior:**
No `tpv205_station_x2_*_x3_*.dat` files in the output directory.

**Expected behavior:**
The same 16 station files the native driver writes, with the same column
layout — so post-processing scripts that consume native TPV205 traces also
work on the spatial driver's output.

**Suggested fix:**
Wire `TPV205StationWriter` into the spatial driver when
`cfg.problem.tag == "tpv205"` (and analogously for TPV31, BP5, etc.):
```diff
+   if (cfg.problem.tag == "tpv205")
+   {
+      auto stations = DefaultStations_TPV205();
+      tpv205_station_writer.Open(cfg.output.output_dir,
+                                 /*prefix=*/"tpv205",
+                                 stations, fault_coords,
+                                 num_fault_local, comm);
+      tpv205_station_writer.WriteStep(cfg.time.t_initial, dof_data);
+   }
```
and call `WriteStep` inside the time loop.  This is the only path to
SCEC-trace parity.

**Test case:**
None at the unit-test level — verified by checking the output directory.
A smoke test could grep for the existence of one expected file:
```cpp
TEST(SpatialDynDriverTpv205Stations, R009_writes_station_files)
{
   // After running the spatial driver, assert
   // <output_dir>/tpv205_station_x2_0_x3_7.5.dat exists.
}
```

---

### [R-010] LOW — Spatial driver lacks the per-step NaN tripwire the native TPV205 driver has; a numerical blow-up runs to completion writing garbage

**Category:** QUALITY

**Description:**
Native (`tpv205_driver.cpp:2447-2464`):
```cpp
// NaN tripwire.
{
   real_t local_nan = std::isnan(Q.Norml2()) ? 1.0 : 0.0;
   real_t global_nan = local_nan;
   MPI_Allreduce(&local_nan, &global_nan, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
   if (global_nan > 0.0)
   {
      std::cerr << "ERROR: NaN detected at step " << step
                << ", t = " << t << " s (rank " << rank << ")\n";
      MPI_Finalize();
      return 1;
   }
}
```
Spatial driver has no equivalent.  Per `CLAUDE.md` "What Constitutes a Regression":
`dt going to zero or NaN` is one of the seven listed regression signatures —
the absence of the tripwire makes that regression silent.

**Suggested fix:**
Add the same tripwire at the end of the spatial driver's per-step loop
(after `Q.Swap(Q_new)`):
```diff
+      // NaN tripwire (parity with tpv205_driver.cpp; per CLAUDE.md
+      // "What Constitutes a Regression" #3).
+      {
+         real_t local_nan = std::isnan(Q.Norml2()) ? 1.0 : 0.0;
+         real_t global_nan = local_nan;
+         MPI_Allreduce(&local_nan, &global_nan, 1,
+                       MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
+         if (global_nan > 0.0)
+         {
+            if (rank == 0)
+            {
+               std::cerr << "ERROR: NaN detected at step " << step
+                         << ", t = " << t << " s.\n";
+            }
+            MPI_Finalize();
+            return 1;
+         }
+      }
```

(No formal test case for LOW; verified by inspection.)

---

## Summary
- Critical issues: **3** (R-001 stress patches miss, R-002 bottom-barrier misses, R-003 strike basis sign flip)
- Moderate issues: **5** (R-004 CFL safety factor, R-005 iterator quadrature, R-006 heterogeneous ctor, R-007 setter ordering, R-008 hypocenter sign land-mine)
- Low issues: **2** (R-009 station writer, R-010 NaN tripwire)
- Plan compliance: **PARTIAL** — the spatial driver wires up the documented Phase 4 / Phase H / Phase N pieces, but the TOML cannot encode TPV205 with the existing schema because of the SPEC-depth-vs-mesh-z mismatch (R-001 / R-002 / R-008) and the `ref_normal` ambiguity (R-003).
- Verdict: **FAIL — must fix before proceeding.**

The combination of R-001 + R-002 + R-003 means the spatial driver running on the
canonical TOML produces a TPV205 simulation that **cannot rupture spontaneously
and, if it did, would slip in the wrong direction with no deep barrier to
constrain it**.  None of these three are subtle — each one is a >100% error in
a load-bearing pre-stress / friction quantity.

The clean route forward is to (a) decide whether the spatial driver's
"raw mesh coords" or the TOML's "spec depth" is canonical, (b) edit one of
them to match, (c) wire the spatial driver's external `FaultBasis` to reuse
the wave operator's internal `ref_normal`, and (d) match the native CFL /
iterator defaults so byte-for-byte cross-driver comparisons become possible.

## Unreviewed Areas
- Time-loop output cadence (`paraview_write`) — not relevant to the "same setup"
  question, only to "same output formatting".  Skimmed; nothing jumped out.
- `WriteTpv104Checkpoint` / `ReadTpv104Checkpoint` serialization parity — not
  exercised by the TPV205 quick-run scenario.
- Mixed-flux dispatch on `[numerics] mixed_flux != "none"` — the TPV205 TOML
  uses `"none"` so this path is dormant, but the R-007 wiring bug would
  surface immediately if a future TOML enabled it.
- Bimaterial Riemann path (Phase R) — covered by the in-tree parity test
  `T-PHASEH-SCALAR-PARITY` per `CLAUDE.md`; trusted to be ~1e-12 relative.
- `SpatialFrictionResolver::ResolveRateState` — out of scope (TPV205 is LSW).
