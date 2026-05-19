# Code Review: `seas_spatial_dyn_driver` vs native TPV102 / TPV104 drivers — 2026-05-18

## Review Scope

- **Goal of the audit (user objective)**: Verify that `tpv102.toml` (or
  `tpv104.toml`) **combined with** `drivers/spatial_dyn_driver.cpp` produces
  the **EXACT SAME** end-to-end setup that `drivers/tpv102_driver.cpp` (or
  `tpv104_driver.cpp`) produces when run with default CLI defaults.  The
  TPV102 / TPV104 / TPV205 sign convention is **the canonical one** —
  no special-case relaxations.

- **Plan / spec documents consulted**
  - `miniapps/seas/CLAUDE.md` (canonical coordinate system + sign
    conventions, project-wide).
  - `safs/project_7.0_alternative/document/PLAN_phase_R_exact_bimaterial_riemann_rev3.md`
    (Phase R goal: single driver for TPV205 / TPV31 / SAFS via TOML).
  - SCEC TPV101/102 spec (https://strike.scec.org/cvws/tpv101_102docs.html);
    SCEC TPV104 spec.

- **Files reviewed**
  - `drivers/spatial_dyn_driver.cpp`
  - `drivers/tpv102_driver.cpp`
  - `drivers/tpv104_driver.cpp`
  - `tpv102/configs/tpv102.toml`
  - `tpv104/configs/tpv104.toml`
  - `config/tpv102_params.hpp` (canonical TPV102 spec constants)
  - `config/tpv104_params.hpp` (canonical TPV104 spec constants)
  - `dynamic/tpv102_setup.hpp` (native `InitializeFaultDOFs`)
  - `dynamic/tpv104_setup.hpp` (native `InitializeFaultDOFs_TPV104` +
    `PopulateVwSideChannel_TPV104`)
  - `dynamic/tpv102_nucleation.hpp` (native
    `ApplyNucleationIncremental_TPV102`)
  - `dynamic/tpv104_nucleation.hpp`
  - `dynamic/spatial_nucleation.hpp` / `.cpp`
    (`ResolveGradualOverstressCompactCircular` +
    `ApplyGradualOverstressCompactCircularIncrement` + `SmoothStep`)
  - `dynamic/spatial_setup.hpp` (`InitializeFaultDOFs_Spatial_RS`)
  - `spatial/code/spatial_friction.hpp` / `.cpp`
    (`SpatialFrictionResolver::ResolveRateState`, `SCECBoxcar`,
    `BoxcarTaperFactor`).
  - `spatial/code/spatial_stress.cpp`
    (`ConstantTensorStressSource::Evaluate`).
  - `fault/fault_geometry_safs_templated.inl`
    (`FaultGeometry::ComputeSAFSParams<StressSource>`).
  - `fault/fault_basis.hpp` (`FaultBasis::ComputeOrientedFrame`).
  - `dynamic/fault_face_flux.cpp`
    (`FaultFaceFlux::ComputeTrialTraction`, `CompleteFromTrial`,
    `CompleteFromVabs`).
  - `dynamic/wave_operator.inl`
    (`ComputeMaxDt`, `EvaluateBulkAtFaultQPsCanonical`, `FaultBasis`
    construction).
  - `dynamic/godunov_flux.cpp` (`BuildRotation`, `BuildRotationInverse`).

- **Domain context**: TPV102 = regularised RS aging-law; TPV104 = slip-law
  with strong rate weakening (FL=103).  Both use a vertical y=0 fault
  with `ref_normal = (0,-1,0)`, `up = (0,0,1)`, strike `t2 = (+1,0,0)`,
  dip `t1 = (0,0,-1)`.  SEAS uses **compression-positive on σ_xx, σ_yy,
  σ_zz** AND **shear signs kept (right-lateral σ_xy = +75 MPa is
  positive in both compression-positive and tension-positive
  conventions)** — quoted from `miniapps/seas/CLAUDE.md`.  The native
  TPV102/104 drivers are the byte-exact regression contract for the
  spatial driver (per `CLAUDE.md` "the native TPV104 (rate-state) and
  TPV205 drivers ... continue to use their own nucleation paths
  verbatim — the spatial-driver schema does NOT touch the TPV*/BP5
  byte-exact regression contract").

---

## Findings

### [R-001] [CRITICAL] [fault/fault_geometry_safs_templated.inl:102-117] — Pre-stress projection produces the WRONG SIGN on `tau_pre.strike` for the canonical TPV convention; spatial driver runs every TPV102/104 case backwards (left-lateral instead of right-lateral)

**Category:** BUG (sign convention)

**Description:**
The native TPV102 driver encodes its strike pre-stress as
```cpp
// dynamic/tpv102_setup.hpp:76-78
d.tau1_0   = 0.0;                      // no dip pre-stress
d.tau2_0   = TPV102Params::tau_ini;    // +75e6  ← POSITIVE
```
(and identically `TPV104Params::tau_ini = +40e6` in
`dynamic/tpv104_setup.hpp`).  The TPV102/TPV104 toml files encode the same
quantity as the Cauchy tensor component
```toml
sigma_xy_pa = 75.0e6        # along-strike shear τ_ini (POSITIVE)
```
(line 100 of `tpv102/configs/tpv102.toml`; line 78 of
`tpv104/configs/tpv104.toml`).

The spatial driver projects this tensor onto the canonical fault frame
via `FaultGeometry::ComputeSAFSParams<StressSource>`:
```cpp
// fault/fault_geometry_safs_templated.inl:97-103
real_t Sn[3];
for (int r = 0; r < 3; ++r)
   Sn[r] = S(r,0)*n[0] + S(r,1)*n[1] + S(r,2)*n[2];
const real_t sigma_n_total = n[0]*Sn[0]  + n[1]*Sn[1]  + n[2]*Sn[2];
const real_t tau1          = t1[0]*Sn[0] + t1[1]*Sn[1] + t1[2]*Sn[2];
const real_t tau2          = t2[0]*Sn[0] + t2[1]*Sn[1] + t2[2]*Sn[2];
```
With the canonical basis (per `CLAUDE.md`):
`n = (0,-1,0)`, `t1 = (0,0,-1)`, `t2 = (+1,0,0)`, and σ tensor
{σ_yy = +120e6, σ_xy = +75e6, others 0}:

```
Sn  = S · n
    = ( S(0,1)·(-1) , S(1,1)·(-1) , S(2,1)·(-1) )
    = ( -75e6      , -120e6     , 0           )

sigma_n_total = n · Sn = (-1)·(-120e6)  = +120e6   ✓ matches native
tau1          = t1 · Sn = (-1)·0        =    0      ✓ matches native
tau2          = t2 · Sn = (+1)·(-75e6) = -75e6      ✗ WRONG SIGN
                                                       (native = +75e6)
```

`InitializeFaultDOFs_Spatial_RS` then plumbs this directly into
`DOFData::tau2_0`:
```cpp
// dynamic/spatial_setup.hpp:131
d.tau2_0   = tau_pre_strike;          // = -75e6   (vs native +75e6)
```

`FaultFaceFlux::CompleteFromTrial` then sums everything:
```cpp
// dynamic/fault_face_flux.cpp:166-167
s.tau2_total = data.tau2_0 + data.tau2_nuc + s.tau2_trial;
```
and `CompleteFromVabs` decomposes the slip vector **parallel** to
`tau2_total`:
```cpp
// dynamic/fault_face_flux.cpp:234-235
s.V1 = s.V_abs * (s.tau1_total) / (strength + data.eta_s * s.V_abs);
s.V2 = s.V_abs * (s.tau2_total) / (strength + data.eta_s * s.V_abs);
```

Consequence at the hypocenter for TPV102/TPV104 with the standard
gradual-overstress nucleation:

| quantity                              | native TPV102/104 | spatial via toml |
|---------------------------------------|-------------------|------------------|
| `tau2_0` (background strike pre-stress) | +75 MPa          | **−75 MPa**     |
| `tau2_nuc` after full ramp (TPV102: Δτ=+25 MPa, strike) | +25 MPa | +25 MPa |
| `tau2_total` at t = T_nuc, ignoring `tau2_trial` | **+100 MPa** | **−50 MPa** |
| Friction strength `|σ_n|·f` (≈ 81 MPa)  | ≈ 81 MPa         | ≈ 81 MPa        |
| Slip-rate sign (`V2 ∝ tau2_total`)     | **+** (right-lateral) | **−** (left-lateral) |
| Rupture broken?                        | yes (100 > 81)   | **no** (|−50| < 81) — fails to nucleate |

For TPV102 the nucleation `delta_tau_strike_pa = +25e6` (the same sign
as the background) is **inert to the projected pre-stress**: instead of
pushing the total past the static threshold it *reduces* `|tau2_total|`
from 75 MPa to 50 MPa.  Even if rupture did nucleate (e.g., via
sufficiently large stress drop), it would propagate as **left-lateral**
slip on the +y block, **opposite to the SCEC spec**.

This bug is masked in the TPV205 verification of Phase R because TPV205
exclusively exercises the BULK wave operator's homogeneous-limit
collapse (the `T-PHASEH-SCALAR-PARITY` gate measures *bulk* parity, not
fault-rupture-direction parity against the native driver's
`tpv205_driver.cpp::Tau2_0_TPV205`).  TPV102 / TPV104 are the **first
TPV configs whose rupture dynamics — not just the bulk waves — must
match the native driver**, and that is exactly where this bug surfaces.

**Trigger:**
Any TOML pre-stress block with `kind = "constant_tensor"` and a non-zero
shear (`sigma_xy_pa`, `sigma_xz_pa`, or `sigma_yz_pa`), under the
canonical ref_normal/up.  Specifically:
- `tpv102/configs/tpv102.toml`: `sigma_xy_pa = 75.0e6` and
  `delta_tau_strike_pa = 25.0e6` → no rupture.
- `tpv104/configs/tpv104.toml`: `sigma_xy_pa = 40.0e6` and
  `delta_tau_strike_pa = 45.0e6` → rupture nucleates but in the
  opposite direction relative to the native driver.

**Actual behavior:**
`tau_pre.strike = -σ_xy` (negative).  `V2` initial direction is along
`-t2 = -x` (left-lateral).  Nucleation Δτ adds same-sign perturbation
that **subtracts** from `|tau2_pre|` instead of adding to it (TPV102
case) or causes rupture in the opposite direction (TPV104 case).

**Expected behavior:**
`tau_pre.strike = +σ_xy` for σ_xy > 0 right-lateral pre-stress
(matches `tpv102_setup.hpp::InitializeFaultDOFs::d.tau2_0 =
TPV102Params::tau_ini` and `tpv205_params.hpp::tau_back = +70e6`).

**Suggested fix:**
The native drivers encode `tau2_0 = +σ_xy` (i.e., the **traction the
+y block exerts on the −y block in the +strike direction**, sign
convention = "right-lateral driving stress is positive").  The
correct projection that recovers this is

  `tau2 = − t2 · ( S · n ) = t2 · ( S · (−n) )`

i.e., project against the inward (+y-pointing) normal instead of the
outward (−y-pointing) `ref_normal`.

Apply the fix in `fault/fault_geometry_safs_templated.inl`:

```diff
       real_t Sn[3];
       for (int r = 0; r < 3; ++r)
       {
          Sn[r] = S(r, 0) * n[0] + S(r, 1) * n[1] + S(r, 2) * n[2];
       }
       const real_t sigma_n_total = n[0]*Sn[0]  + n[1]*Sn[1]  + n[2]*Sn[2];
-      const real_t tau1          = t1[0]*Sn[0] + t1[1]*Sn[1] + t1[2]*Sn[2];
-      const real_t tau2          = t2[0]*Sn[0] + t2[1]*Sn[1] + t2[2]*Sn[2];
+      // CLAUDE.md project-wide sign convention: positive tau2_0 represents
+      // right-lateral driving stress in the +strike direction, matching
+      // the native TPV102/104/205 driver's
+      //   d.tau2_0 = TPV{102,104,205}Params::tau_ini   (positive).
+      // The Cauchy projection T = σ·n with n = (0,-1,0) returns the
+      // traction the +y side exerts on the −y side (left-lateral by
+      // Newton's 3rd-law convention), giving tau2 = −σ_xy.  Flip the
+      // sign so the spatial driver matches the native sign convention
+      // unambiguously.
+      const real_t tau1          = -(t1[0]*Sn[0] + t1[1]*Sn[1] + t1[2]*Sn[2]);
+      const real_t tau2          = -(t2[0]*Sn[0] + t2[1]*Sn[1] + t2[2]*Sn[2]);
       // ... sigma_n stays as-is (n·S·n is sign-invariant under n → −n).
```

The same flip must be inspected at every SAFS source path that builds
`tau_pre_`:
- `ApplyCsmStressSidecar` (CSM sidecar projection).
- The non-templated `ComputeSAFSParams(const StressField3D&, …)`
  overload (`fault/fault_geometry_safs.inl`).
- The `[stress.patch]` overlay path (`ConstantTensorWithPatchesStressSource`)
  — patches go through the same projection so a single fix at the
  projection site is sufficient.

`sigma_n` is **unchanged** by this fix (it is `n·S·n`, invariant under
n→−n) so the normal-stress branch (R-002 below for nucleation) is not
disturbed.

**Test case (proposed unit test for `tests/unit/test_compute_safs_params.cpp`):**
```cpp
// T-R001 — canonical TPV102 stress source projects to native sign.
TEST(ComputeSAFSParamsSign, TPV102CanonicalGivesPositiveStrike)
{
   // Single fault DOF on the y=0 plane at the hypocenter (0,0,-7500).
   // Basis: canonical TPV (ref_normal=(0,-1,0), up=(0,0,1)).
   Vector dof_coords_3d(3);  dof_coords_3d(0)=0; dof_coords_3d(1)=0;
                             dof_coords_3d(2)=-7500.0;
   DenseMatrix dof_basis(9, 1);
   // n  = (0,-1,0)
   dof_basis(0,0)=0; dof_basis(1,0)=-1; dof_basis(2,0)=0;
   // t1 = dip = (0,0,-1)
   dof_basis(3,0)=0; dof_basis(4,0)=0;  dof_basis(5,0)=-1;
   // t2 = strike = (+1,0,0)
   dof_basis(6,0)=1; dof_basis(7,0)=0;  dof_basis(8,0)=0;

   FaultGeometry<Mesh> geom(...);  // ndof=1, BP5 ctor
   ConstantTensorStressSource src(
      /*sxx=*/0.0, /*syy=*/120.0e6, /*szz=*/0.0,
      /*sxy=*/75.0e6, /*syz=*/0.0, /*sxz=*/0.0);
   geom.ComputeSAFSParams(src, /*P_p=*/0.0, /*grad=*/0.0, /*min_sn=*/0.0);

   const Vector& tp = geom.GetTauPre();
   const Vector& sn = geom.sigma_n_per_dof();

   EXPECT_NEAR(sn(0),       120.0e6, 1e-6);  // sigma_n_total: invariant
   EXPECT_NEAR(tp(2*0+0),   0.0,     1e-6);  // tau_dip = 0
   EXPECT_NEAR(tp(2*0+1), +75.0e6,   1e-6);  // tau_strike = +75e6 (R-001)
   //                ^^^^^ MUST be POSITIVE — matches native
   //                drivers/tpv102_driver.cpp `d.tau2_0 = +tau_ini`.
}
```

A second test for the negative-σ_xy mirror case should also be added so
the fix is anchored on both signs (catches an accidental |·| or
abs-then-sign-restore implementation).

---

### [R-002] [CRITICAL] [drivers/spatial_dyn_driver.cpp:1370-1383] — Initial slip-rate direction follows the wrong-signed `tau_pre`; left-lateral creep initial condition instead of right-lateral

**Category:** BUG (downstream of R-001)

**Description:**
After R-001's mis-signed projection, the spatial driver's init loop
seeds the initial slip-rate vector parallel to `tau_pre`:
```cpp
// drivers/spatial_dyn_driver.cpp:1367-1380
// Slip rate direction PARALLEL to tau_pre (R-801 BP5 frame …)
d.slip_rate = V_init_i;
if (tau_abs > 0.0)
{
   d.V1 = (V_init_i / tau_abs) * tau1_pre;
   d.V2 = (V_init_i / tau_abs) * tau2_pre;
}
```
This is the CORRECT rule (`CLAUDE.md` — "Slip rate direction: PARALLEL
to traction tau"), but `tau2_pre` carries R-001's wrong sign, so the
init produces `V2 = (V_init_i / 75e6) · (−75e6) = −V_init_i`.

Native `InitializeFaultDOFs` hard-codes
```cpp
// dynamic/tpv102_setup.hpp:104-106
d.slip_rate = TPV102Params::V_ini;
d.V1 = 0.0;                            // no dip slip rate
d.V2 = TPV102Params::V_ini;            // along-strike initial slip rate (+x)
```
and identically TPV104's `InitializeFaultDOFs_TPV104`.

**Trigger:** Identical to R-001 (any TPV-canonical config).

**Actual behavior:**
`V2 = -V_init_i` at every fault DOF.

**Expected behavior:**
`V2 = +V_init_i` (same as native).

**Suggested fix:**
None at this site — the rule `V_vec = (V_abs/|τ|) τ_vec` is correct.
The fix lives in R-001 (sign of `tau_pre.strike`).  After R-001 lands,
this site automatically produces `+V_init_i`.

**Test case (regression for `tests/unit/test_spatial_setup.cpp`):**
```cpp
// T-R002 — after the spatial RS init, initial V_vec aligns with
// native TPV102 (V2 = +V_init, V1 = 0) for the canonical config.
TEST(SpatialRSInit, TPV102InitialSlipRateMatchesNative)
{
   // Build a 1-DOF fault on the canonical frame; ConstantTensorStressSource
   // with σ_xy = +75e6; rs params from tpv102.toml's defaults.
   ...
   InitializeFaultDOFs_Spatial_RS<Mesh>(dof_data, /*ndof=*/1, …);

   // After R-001 fix:
   EXPECT_NEAR(dof_data[0].V2, +1.0e-12, 1e-24);   // = +V_init
   EXPECT_NEAR(dof_data[0].V1,  0.0,     1e-24);
   EXPECT_NEAR(dof_data[0].slip_rate, 1.0e-12, 1e-24);
}
```

---

### [R-003] [CRITICAL] [tpv102/configs/tpv102.toml + tpv104/configs/tpv104.toml] — Missing `cfl_safety = "dg"`; spatial driver applies a 9× larger dt than the native driver at p=1

**Category:** DEVIATION (TOML completeness)

**Description:**
The native TPV102 / TPV104 drivers compute the time step as
```cpp
// drivers/tpv102_driver.cpp:1265
real_t cfl    = cfl_factor / (3.0 * (2.0 * order + 1.0));
real_t dt_cfl = wave.ComputeMaxDt(cfl);
```
i.e., the `1/(3·(2p+1))` DG order-stability factor is applied to `cfl`
*before* feeding it to the wave operator (`cfl_factor = 0.5` is the CLI
default; at `order = 1` this yields effective cfl = 0.5/9 ≈ 0.0556).

The spatial driver reads `cfl_safety` from `[numerics]` and only
applies the DG factor when it is `"dg"`:
```cpp
// drivers/spatial_dyn_driver.cpp:1463-1475
real_t cfl_for_max_dt = cfg.numerics.cfl;
if (cfg.numerics.cfl_safety == "dg")
{
   cfl_for_max_dt /= (3.0 * (2.0 * cfg.mesh.order + 1.0));
   …
}
const real_t dt_cfl = wave.ComputeMaxDt(cfl_for_max_dt);
```
Default `cfl_safety = "raw"` (`spatial/code/spatial_friction.hpp:202`).

`tpv102.toml` and `tpv104.toml` set
```toml
[numerics]
cfl        = 0.5
```
**without** `cfl_safety = "dg"`, so the driver feeds `cfl = 0.5`
directly to `ComputeMaxDt`.  With h_min ≈ 1000 m and cp = 6000 m/s,
this gives `dt_cfl = 0.5 · 1000/6000 ≈ 0.0833 s` (vs. native's
≈ 0.00926 s).

The `dt_max = "0.01s"` cap in both TOMLs hides this from a casual
inspection — it traps `dt` at 0.01 s — but this is still ~8 % larger
than native's `dt = 0.00926 s` and the *number of substeps over the
whole 12 s run differs by ~8 %*.  Byte parity is therefore impossible
even before considering R-001.

By contrast, `tpv205.toml:168` correctly declares
```toml
cfl_safety     = "dg"        # cfl /= (3*(2p+1)) — matches tpv205_driver
```

**Trigger:** Any run of `seas_spatial_dyn_driver --config
tpv102.toml` or `--config tpv104.toml` with no `--cfl` CLI override.

**Actual behavior:** `dt = 0.01 s` (capped by `dt_max`), but the
underlying `dt_cfl` is computed without the DG safety factor.  Removing
`dt_max` would put `dt` at ~0.083 s — well past the CFL stability limit
of the DG operator at p = 1.

**Expected behavior:** `dt = 0.00926 s` matching native, i.e.,
`cfl_safety = "dg"` opted-in.

**Suggested fix:**

Edit `miniapps/seas/tpv102/configs/tpv102.toml`:
```diff
 [numerics]
 ader_order = 2
 mixed_flux = "none"
 cfl        = 0.5
+cfl_safety = "dg"           # cfl /= (3*(2p+1)) — matches tpv102_driver.cpp:1265
 use_pml    = false
```

Edit `miniapps/seas/tpv104/configs/tpv104.toml`:
```diff
 [numerics]
 ader_order = 2
 mixed_flux = "none"
 cfl        = 0.5
+cfl_safety = "dg"           # cfl /= (3*(2p+1)) — matches tpv104_driver.cpp:1325
 use_pml    = false
```

**Test case:**
```cpp
// T-R003 — TPV102/104 tomls request the DG safety factor.
TEST(SpatialTomlByteParity, TPV102CflSafetyIsDG)
{
   const auto cfg = spatial::LoadSpatialFrictionConfig(
      "tpv102/configs/tpv102.toml");
   EXPECT_EQ(cfg.numerics.cfl_safety, "dg");
   EXPECT_DOUBLE_EQ(cfg.numerics.cfl, 0.5);
}
TEST(SpatialTomlByteParity, TPV104CflSafetyIsDG)
{
   const auto cfg = spatial::LoadSpatialFrictionConfig(
      "tpv104/configs/tpv104.toml");
   EXPECT_EQ(cfg.numerics.cfl_safety, "dg");
   EXPECT_DOUBLE_EQ(cfg.numerics.cfl, 0.5);
}
```

---

### [R-004] [MODERATE] [drivers/spatial_dyn_driver.cpp:2128-2138] — Spatial driver UNCONDITIONALLY runs the per-substep iterator path; native default is one-shot, so a default-vs-default comparison cannot byte-match

**Category:** DEVIATION (dispatch path)

**Description:**
The native TPV102 / TPV104 drivers default to the one-shot
`wave.AdvanceADER` path:
```cpp
// drivers/tpv102_driver.cpp:747-748
std::string fault_iterator =
   GetStringArg(argc, argv, "--fault-iterator", "one-shot");
```
and only enter `AdvanceADERWithSubStep` if the user opts in with
`--fault-iterator substep`.

The spatial driver always composes a per-sub-step iterator into the
production time loop:
```cpp
// drivers/spatial_dyn_driver.cpp:2134-2138
AdvanceADERWithSubStep_Spatial(wave, set_substeps, do_iterate,
                               substep_deltaT, substep_weights,
                               Q, dt_step,
                               cfg.numerics.ader_order, t, Q_new,
                               nuc_cb);
```
At `ader_order = 2` (the TPV102/104 toml setting), the substep path
differs from the one-shot path at O(dt²) — verified by the
`T_TPV104_SSI_3` contract: bit-equivalence holds **only at O = 1**.

This means even if R-001 and R-003 are fixed, the spatial-driver TPV102
run will **never byte-match** a default-CLI `tpv102_driver.cpp` run.
The valid comparison is `tpv102_driver.cpp --fault-iterator substep`
vs. `spatial_dyn_driver.cpp --config tpv102.toml` (and analogously for
TPV104).  This is a real comparison contract issue — anyone trying to
"check" the spatial driver against the native driver's default-output
file will see disagreement growing with `dt²`.

**Trigger:** Any run with `ader_order >= 2`.

**Actual behavior:** Always substep dispatch.

**Expected behavior (one of):**
(a) Surface the dispatch choice in the TOML/CLI of the spatial driver —
default `one-shot` to match native byte-by-byte; `substep` opt-in.
(b) Keep substep-only but document the comparison contract change
prominently in `tpv102.toml` / `tpv104.toml` headers.  The byte-exact
reference run must be regenerated by re-running the native driver with
`--fault-iterator substep`.

**Suggested fix (recommended: option (a))**:
Add a `[numerics] fault_iterator = "one-shot"` / `"substep"` knob to
`SpatialFrictionConfig::NumericsBlock`, parse it in
`spatial/code/spatial_friction.cpp` with default `"one-shot"`, then
gate the dispatch in the driver:

```diff
+   const bool use_substep =
+      (cfg.numerics.fault_iterator == "substep");
    ...
    for (int step = step0; step < nsteps; ++step)
    {
       …
-      AdvanceADERWithSubStep_Spatial(wave, set_substeps, do_iterate,
-                                     substep_deltaT, substep_weights,
-                                     Q, dt_step,
-                                     cfg.numerics.ader_order, t, Q_new,
-                                     nuc_cb);
+      if (use_substep)
+      {
+         AdvanceADERWithSubStep_Spatial(wave, set_substeps, do_iterate,
+                                        substep_deltaT, substep_weights,
+                                        Q, dt_step,
+                                        cfg.numerics.ader_order, t, Q_new,
+                                        nuc_cb);
+      }
+      else
+      {
+         // One-shot must still install nuc_cb's accumulation BEFORE the
+         // call: tau{1,2}_nuc are read inside wave.AdvanceADER -> EvaluateADER.
+         // Apply ONE smoothStep increment for the whole macrostep.
+         nuc_cb(t + dt_step, dt_step);
+         wave.AdvanceADER(Q, dt_step, cfg.numerics.ader_order, Q_new);
+      }
    }
```

Then `tpv102.toml` adds `fault_iterator = "one-shot"` to byte-match
the native default.

**Test case (regression for `tests/unit/test_spatial_dyn_driver_toml.cpp`):**
```cpp
// T-R004 — fault_iterator default is "one-shot"; tomls can override.
TEST(SpatialDynNumerics, DefaultFaultIteratorIsOneShot)
{
   const auto cfg = spatial::LoadSpatialFrictionConfig(/* tpv102 toml */);
   EXPECT_EQ(cfg.numerics.fault_iterator, "one-shot");
}
```

---

### [R-005] [MODERATE] [drivers/spatial_dyn_driver.cpp:1037-1067 + tpv104/configs/tpv104.toml + tpv102/configs/tpv102.toml] — Redundant `[fault_geometry]` block in tomls duplicates a project-wide invariant and risks the hard-fail guard tripping under future TOML edits

**Category:** QUALITY / DEVIATION (TOML completeness)

**Description:**
Both tomls explicitly set
```toml
[fault_geometry]
ref_normal = [0.0, -1.0, 0.0]
up         = [0.0,  0.0, 1.0]
kind       = "tpv102_rs"  / "tpv104_srw"
```
The driver enforces (CLAUDE.md "canonical coordinate system, project-wide")
that these values must be exactly `(0,-1,0)` and `(0,0,1)` and aborts
otherwise:
```cpp
// drivers/spatial_dyn_driver.cpp:1046-1067
MFEM_VERIFY(n_matches && up_matches,
            "spatial_dyn_driver: [fault_geometry].ref_normal / up = ("
            << … << ") does NOT match the wave operator's internal "
            "Tandem convention (0,-1,0) / (0,0,1) at "
            "dynamic/wave_operator.inl:349-350. …");
```
The values in the tomls are not consumed for anything else (the driver
ignores the user-provided ref_normal/up and reads `wave.GetFaultBasis()`
directly per R-003 of the spatial driver's own internal review).  Per
`CLAUDE.md`:
> "TOML configs (TPV102/104/205, BP5/SAFS) should NOT re-state these
>  values; they inherit them from the schema default."

Listing them in the toml has two negative effects:
1. A future copy-paste / typo (e.g., `ref_normal = [0.0, +1.0, 0.0]`,
   matching the SAFS legacy schema default) will silently hard-fail the
   driver at startup, rather than the driver falling through to the
   schema default.
2. Anyone reading the toml learns the wrong lesson — that these are
   per-config knobs.

The `kind` string (`"tpv102_rs"` / `"tpv104_srw"`) is also
non-functional — search for `tpv102_rs` / `tpv104_srw` in the spatial
codebase produces zero hits (it is documentation-only).

**Trigger:** TOML edit that perturbs `ref_normal` / `up` from the
schema canonical values.

**Actual behavior:** Driver aborts at startup with the R-003 message.

**Expected behavior:** TOML omits these fields; driver/schema supplies
canonical defaults.

**Suggested fix:**

```diff
-[fault_geometry]
-# TPV102 fault: vertical y=0 plane, +z = up (towards free surface),
-# ref_normal points into the - side (BP5 / TPV102 convention).
-# tangent1 = dip = (0, 0, -1)  (down into earth)
-# tangent2 = strike = (+1, 0, 0)
-ref_normal = [0.0, -1.0, 0.0]
-up         = [0.0,  0.0, 1.0]
-kind       = "tpv102_rs"
+# [fault_geometry] intentionally omitted — uses the project-wide
+# canonical frame (CLAUDE.md "Canonical Coordinate System"):
+#   ref_normal = (0, -1, 0), up = (0, 0, 1).
+# The spatial driver hard-fails if a toml overrides these.
```

**Test case:** existing test at
`tests/unit/test_spatial_dyn_driver_toml.cpp` should be extended:
```cpp
// T-R005 — TPV102/104 tomls do NOT carry a [fault_geometry] block, so
// they inherit the canonical project-wide default.
TEST(SpatialTomlConvention, TPV102OmitsFaultGeometry)
{
   const auto cfg = spatial::LoadSpatialFrictionConfig(
      "tpv102/configs/tpv102.toml");
   EXPECT_EQ(cfg.fault_geometry.ref_normal[0],  0.0);
   EXPECT_EQ(cfg.fault_geometry.ref_normal[1], -1.0);  // schema default
   EXPECT_EQ(cfg.fault_geometry.ref_normal[2],  0.0);
}
```

---

### [R-006] [MODERATE] [POSSIBLE] [drivers/spatial_dyn_driver.cpp:1330-1383] — Spatial RS init re-derives ψ from local arrays, BUT does NOT trip-check the per-DOF re-derivation against `dof_data[i].a` / `Dc` set by the underlying `InitializeFaultDOFs_Spatial_RS` — a silent decoupling can cause ψ-vs-a inconsistency

**Category:** ASSUMPTION (silent-coupling)

**Description:**
After the underlying `InitializeFaultDOFs_Spatial_RS` returns, the
spatial driver overrides `d.psi` per-DOF using `rs.a(i)`, `rs.V_0(i)`,
`rs.V_init(i)`:
```cpp
// drivers/spatial_dyn_driver.cpp:1343-1365
const real_t a_i      = rs.a(i);
const real_t V0_i     = rs.V_0(i);
const real_t V_init_i = std::max(rs.V_init(i), 1.0e-300);
const real_t sn_i     = sn_eff(i);
const real_t tau1_pre = tau_pre(2 * i + 0);
const real_t tau2_pre = tau_pre(2 * i + 1);
const real_t tau_abs  = std::sqrt(tau1_pre * tau1_pre
                                  + tau2_pre * tau2_pre);
…
d.psi = a_i * (absC + std::log(x / 2.0 * -sign_c
                               * std::expm1(-2.0 * absC)));
```
Critically, this ψ inversion uses the **per-DOF resolver values**, not
the values written into `dof_data[i].a / Dc / sigma_n0`.  If the
underlying `InitializeFaultDOFs_Spatial_RS` ever writes a different
`d.a` (e.g., a future refactor that adds a per-DOF cap, smoothing, or
material-coupled `a`), ψ will silently drift from the friction-law
fixed point.

The native TPV102 driver computes ψ from the same source-of-truth that
populates `d.a` (both via `ComputeA(along_strike, down_dip)` and
`ComputeInitialPsi(d.a)` inside `InitializeFaultDOFs`):
```cpp
// dynamic/tpv102_setup.hpp:96-100
d.a  = ComputeA(along_strike, down_dip);
d.Dc = TPV102Params::Dc;
d.psi = ComputeInitialPsi(d.a);
```

The native code is self-consistent by construction: any change to `a`
forces a matching change to ψ.  The spatial driver's split (resolver
→ underlying init → driver-level ψ override) breaks that invariant.

**Trigger:**
- Today: spatial drivers' `InitializeFaultDOFs_Spatial_RS` writes
  `d.a = rs.a(i)`, so the override is consistent.  However the contract
  is `d.psi` uses the SAME `a_i` the friction solver reads.
- Future: anything that introduces per-DOF post-init mutation to
  `d.a`, `d.Dc`, or `sigma_n0` (e.g., friction-clipping, sub-grid
  averaging, an LSW/RS coupling fix) before the time loop starts.

**Actual behavior:** Currently consistent (no observed regression),
but the safety invariant is missing.

**Expected behavior:** Either compute ψ from `d.a, d.Dc, d.sigma_n0,
d.tau{1,2}_0` (the friction solver's own inputs), OR add an explicit
post-init check.

**Suggested fix:** Compute ψ from `d`'s own fields:

```diff
-      const real_t a_i      = rs.a(i);
-      const real_t V0_i     = rs.V_0(i);
-      const real_t V_init_i = std::max(rs.V_init(i),
-                                       static_cast<real_t>(1.0e-300));
-      const real_t sn_i     = sn_eff(i);
+      const real_t a_i      = d.a;
+      const real_t V0_i     = rs.V_0(i);  // V_0 stays on resolver side; not
+                                          // duplicated into DOFData.
+      const real_t V_init_i = std::max(rs.V_init(i),
+                                       static_cast<real_t>(1.0e-300));
+      const real_t sn_i     = std::abs(d.sigma_n0);
        const real_t tau1_pre = tau_pre(2 * i + 0);
        const real_t tau2_pre = tau_pre(2 * i + 1);
```
and add a `verify_equilibrium_error < 1e-6` post-condition (mirrors BP5
debug v13).

**Test case:**
```cpp
// T-R006 — after the spatial RS init, the friction equation
//   tau = sigma_n * a * asinh((V/(2 V_0)) exp(psi/a))
// is satisfied at every fault DOF to rel-error < 1e-12.
TEST(SpatialRSInit, EquilibriumResidualBelowTol)
{
   ...
   InitializeFaultDOFs_Spatial_RS<Mesh>(dof_data, ndof, …);
   driver_psi_override(...);  // the L1343-1365 block
   for (int i = 0; i < ndof; ++i)
   {
      const auto& d = dof_data[i];
      const real_t V = d.slip_rate;
      const real_t f = d.a * std::asinh(V / (2*V0)
                                         * std::exp(d.psi / d.a));
      const real_t tau_pred = std::abs(d.sigma_n0) * f;
      const real_t tau_obs  = std::sqrt(d.tau1_0*d.tau1_0
                                         + d.tau2_0*d.tau2_0);
      EXPECT_NEAR(tau_pred, tau_obs, 1e-6 * tau_obs);
   }
}
```

---

### [R-007] [MODERATE] [tpv102/configs/tpv102.toml + tpv104/configs/tpv104.toml] — Tomls leave `f_w_default` set to 0.2 unconditionally, but TPV102 (aging law) does not use `f_w` — silent default may mask a future regression where the wrong friction law reads `f_w`

**Category:** QUALITY (defensive defaults)

**Description:**
The spatial RS resolver validates `f_w_default in (0,1)` ONLY when
`state_evolution == slip_law_srw`:
```cpp
// spatial/code/spatial_friction.cpp:628-630
if (state_evolution == StateEvolutionKind::SlipLawStrongRateWeakening)
{
   MFEM_VERIFY(out.f_w_default > 0.0 && out.f_w_default < 1.0, …);
}
```
For aging law (TPV102), `f_w_default` is read by `toml_real` with
default 0.2 and stored but never consumed.  This is a `tpv102.toml`
hygiene issue: a future refactor that wires the aging law iterator
through `SlipLawSRWPsi` (because they share a base class) would silently
pick up `f_w = 0.2` instead of failing loudly.

TPV104 toml DOES set `f_w_default = 0.2` correctly (spec value).

Note: this is independent of R-001..R-005.

**Trigger:** A future change that exposes `f_w` to the aging-law path.

**Actual behavior:** TPV102 toml inherits `f_w_default = 0.2` silently.

**Expected behavior:** TPV102 toml omits `f_w_default` (aging law does
not use it), OR explicitly sets it to a sentinel that the aging-law
init treats as "must not be consulted."  Tpv104 toml correctly sets
`f_w_default = 0.2`.

**Suggested fix (TPV102 only):** No change required in TPV102 toml
today (the field is silently ignored under the aging-law branch).  Add
a comment in the schema docstring so the next reader does not delete
the validator under the assumption it is dead code:

(`spatial/code/spatial_friction.hpp` doc comment):
```
/// f_w_default and V_w_default are CONSUMED ONLY when state_evolution
/// == slip_law_srw.  The aging-law iterator ignores both fields; do
/// not assume they are honoured at the aging-law DOF.
```

**Test case:** N/A — pure documentation hardening.

---

### [R-008] [LOW] [tpv102/configs/tpv102.toml:128-129 + tpv104/configs/tpv104.toml:107] — `checkpoint_every_steps = 10000` is larger than the entire run length (TPV102 nsteps ≈ 1200 at the configured dt); no checkpoint is ever written

**Category:** QUALITY (config hygiene)

**Description:**
At `tfinal = 12 s` and `dt = 0.01 s` (capped by `dt_max`),
`nsteps ≈ 1200`.  With `checkpoint_every_steps = 10000`, the spatial
driver's per-step checkpoint trigger at L2157
```cpp
if (cfg.output.checkpoint_every_steps > 0
    && (step + 1) % cfg.output.checkpoint_every_steps == 0)
```
never fires for either TPV102 or TPV104.  The end-of-run checkpoint at
L2181 does fire (it gates on `> 0`, not on the cadence) so a single
final checkpoint is produced.

This is intentional in `tpv104.toml`'s sibling but worth flagging
explicitly so a future user enabling restart does not wonder why no
intermediate checkpoints exist.

**Actual behavior:** Only the final checkpoint is written.

**Expected behavior:** Either lower the cadence (e.g., `every_steps =
200`) or annotate the toml.

**Suggested fix:** purely cosmetic; add a comment to both tomls:
```toml
# Run length is ~1200 steps; this large value means no INTERMEDIATE
# checkpoint is written.  Only the end-of-run checkpoint fires.
# Lower to e.g. 200 if you want intermediate restart points.
checkpoint_every_steps  = 10000
```
No test required.

---

### [R-009] [MODERATE] [tpv102/configs/tpv102.toml:114 + tpv104/configs/tpv104.toml:91-92] — `dt_max = "0.01s"` is the wrong knob to enforce CFL; with `cfl_safety = "raw"` it silently masks a CFL violation, and after fixing R-003 the cap may force a smaller-than-cfl dt that hurts performance

**Category:** ASSUMPTION (knob coupling)

**Description:**
The `[time] dt_max` knob caps the auto-selected `dt_cfl`.  In the
current toml,
- `cfl_safety` defaults to `"raw"` (R-003) → `dt_cfl ≈ 0.083 s`.
- `dt_max = 0.01 s` → spatial driver runs at dt = 0.01 s.

This is **8× larger than native** (0.00926 s).  Worse: without the
`dt_max` cap, the driver would run at the unstable 0.083 s.  So
`dt_max` is silently absorbing what should be the CFL stability factor.

After R-003 lands (`cfl_safety = "dg"`), `dt_cfl ≈ 0.00926 s` and the
0.01 s cap is then **inactive** (since 0.00926 < 0.01).  So the
`dt_max` value becomes vestigial.

**Trigger:** Removing `cfl_safety = "dg"` after a future toml edit.

**Actual behavior:** Today: `dt = 0.01 s` (close to native but not
identical).  After R-003 fix: `dt = 0.00926 s` (matches native).

**Expected behavior:** The CFL stability should come from the CFL knob
(`cfl` + `cfl_safety`), not the run-length cap `dt_max`.

**Suggested fix:** After R-003 lands, simplify either by removing
`dt_max` or by setting it to a value that ONLY caps for long-run
output-cadence purposes (e.g., `0.05 s` for the ParaView fault-dt
hint):
```diff
 [time]
 tfinal     = "12s"
 t_initial  = 0.0
 dt_initial = "auto"
-dt_max     = "0.01s"
+# dt is selected by [numerics].cfl / cfl_safety; dt_max is intentionally
+# left at the schema default (infinity) so the CFL machinery is the
+# single source of truth for stability.
```
Re-test: with R-003 applied + this change, `dt` must equal the native
driver's `dt_cfl` to within FP rounding.

**Test case:** added under R-003 (the same assertion of `dt == dt_cfl_native`).

---

### [R-010] [LOW] [drivers/spatial_dyn_driver.cpp:736-746 (banner)] — Spatial driver's startup banner advertises `nucleation: gradual_overstress (enabled)` regardless of the actual nucleation kind chosen, so the operator cannot verify from logs which nucleation path was taken

**Category:** QUALITY (log ambiguity)

**Description:**
```cpp
// drivers/spatial_dyn_driver.cpp:736-740
<< "nucleation:       "
<< (cfg.nucleation.enabled
    ? "gradual_overstress (enabled)"
    : "DISABLED")
```
This hardcodes `gradual_overstress` in the banner string even when the
TOML selected `gradual_overstress_compact_circular` (the TPV102/104
spec kind).  An operator scanning the log for "compact_circular" would
think the spatial driver fell back to the Gaussian flavour and stop
investigating R-001/R-002.

**Trigger:** Always under TPV102/104 tomls (which use
`gradual_overstress_compact_circular`).

**Actual behavior:** Banner always says `gradual_overstress`.

**Expected behavior:** Banner echoes `cfg.nucleation.kind`.

**Suggested fix:**

```diff
-      << "nucleation:       "
-      << (cfg.nucleation.enabled
-          ? "gradual_overstress (enabled)"
-          : "DISABLED")
+      << "nucleation:       "
+      << (cfg.nucleation.enabled
+          ? std::string{spatial::ToString(cfg.nucleation.kind)} + " (enabled)"
+          : std::string{"DISABLED"})
```
with `spatial::ToString(NucleationKind)` added next to the enum.

**Test case:** N/A (purely diagnostic).

---

### [R-011] [LOW] [drivers/spatial_dyn_driver.cpp:1442-1469] — `print_derived` post-CFL diagnostics are gated by `print_derived && is_lsw`; the TPV102/104 rate-state branches get no equivalent sanity check on dt vs nucleation rise time

**Category:** QUALITY (missing diagnostic)

**Description:**
The `PrintDerivedAndCheck` block at L1437-1470 only runs for the LSW
path.  TPV102/104 run through the RS branch, so the operator gets no
warning when, e.g.,
- `dt > T_nuc_s / 20` (under-resolving the nucleation ramp)
- `t_reflect < tfinal` and `--pml` not set (R-107 reflection warning
  IS issued earlier at L977 — this works for RS too).

This is a missed opportunity to flag the `dt = 0.01 s` vs `T_nuc = 1
s` ratio (= 100 sub-divisions of the nucleation ramp at native CFL;
under R-003-uncorrected toml, this drops to ~120).

**Trigger:** Any RS run.

**Actual behavior:** No RS-specific diagnostic.

**Expected behavior:** Either add an RS-side `PrintDerivedAndCheck`,
or note in the driver doc that the LSW-only gate is intentional and
TPV102/104 must be inspected via separate logs.

**Suggested fix:** Out of scope for the current review; flagged for
future hardening.

**Test case:** N/A.

---

## Summary

| | Count |
|--|--|
| Critical issues   | 3 (R-001, R-002, R-003) |
| Moderate issues   | 5 (R-004, R-005, R-006, R-007, R-009) |
| Low issues        | 3 (R-008, R-010, R-011) |
| Plan compliance   | **INCOMPLETE** — three critical sign / CFL / dispatch deviations make `tpv102.toml + spatial_dyn_driver.cpp` produce a setup that is **not** the exact-same setup as `tpv102_driver.cpp`. |
| Verdict           | **FAIL — must fix R-001, R-002, R-003 before treating the spatial driver as a TPV102/104 substitute**.  R-001 is the root cause of R-002 and of the nucleation sign mismatch.  R-003 + R-004 are the bare-minimum extra fixes required for byte parity with the native driver. |

**Byte-equivalence path forward (post-fix):**
1. Apply R-001 (sign of `tau_pre` projection).  R-002 then resolves
   automatically.
2. Apply R-003 (`cfl_safety = "dg"` in both tomls).
3. Apply R-004 (default `fault_iterator = "one-shot"` matching native).
4. Re-run `seas_spatial_dyn_driver --config tpv102.toml` AND
   `seas_tpv102_driver --mesh tpv102/mesh/tpv102_1000m.msh --tfinal 12`
   with **identical** CFL/output-dt/mesh, then run a per-station
   sigma_xy, V2, slip2, psi trace diff between the two runs' .dat
   files.  Acceptance: RMS error < 1e-12 relative per station per field
   (the LU-rounding tolerance for ADER homogeneous matches the
   T-PHASEH-SCALAR-PARITY gate).
5. Repeat for TPV104 (`tpv104.toml`).
6. Once both pass, the spatial driver is a drop-in replacement for the
   native TPV102/104 drivers, and the same toml plumbing extends
   directly to SAFS / TPV31 once the canonical-frame outliers are
   resolved.

## Unreviewed Areas

- **`paraview_write` writers** (`drivers/spatial_dyn_driver.cpp:1980-2060`):
  spot-checked the field publishing but did not verify field-name parity
  with `tpv102_driver.cpp` / `tpv104_driver.cpp`'s ParaView paths.  In
  scope only if the byte-parity workflow demands ParaView-frame
  comparison (the `.dat` station files are the SCEC-mandated
  comparison target).
- **`ApplyCsmStressSidecar`** (used when `stress.kind = "sidecar_hdf5"`):
  not exercised by TPV102/104 tomls (`kind = "constant_tensor"`) but
  shares R-001's sign-projection bug if the SidecarHDF5 path projects
  the same way; left out of this review because the user's objective
  is TPV102/104 alone.
- **TPV31 outlier (CLAUDE.md §Outliers)**: trips R-003 guard already;
  not part of this review.
- **Substep iterator's per-sub-step ψ cadence** (R7-007 in
  `tpv104_driver.cpp:2545-2554`): the macro-step ψ vs per-sub-step ψ
  diff is O(dt²) and contributes to the R-004 dispatch-path divergence
  budget but is independent of R-001/R-002/R-003.
- **MaterialField::EvalAt at the fault QP**: spatial driver uses
  `dof_ips[i]` (the per-DOF reference IntegrationPoint) — correct per
  CLAUDE.md R-111; native TPV102/104 driver uses the homogeneous
  constants directly via `TPV{102,104}Params`.  For TPV102/104's
  `material = "constant"`, `MaterialField::Mode::Constant` returns the
  same constants — byte-equivalent.
