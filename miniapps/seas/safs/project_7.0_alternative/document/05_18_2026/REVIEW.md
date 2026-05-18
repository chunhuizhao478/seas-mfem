# Code Review: PLAN_first_safs_run.md (2026-05-18)

## Review Scope

- Plan: `safs/project_7.0_alternative/document/05_18_2026/PLAN_first_safs_run.md`
- Files verified against:
  - `drivers/spatial_dyn_driver.cpp:228-285` — actual `dof_basis` construction
  - `dynamic/spatial_setup.hpp:283-284, 329-330, 386-387` — `tau_pre` interleaved layout
  - `dynamic/wave_state.hpp:23-37` — `NUM_STATE` enum ordering
  - `dynamic/tpv205_substep_iterator.cpp:241-327` — `AdvanceWithSubStepStates` body + sub-step loop
  - `spatial/code/spatial_friction.hpp/.cpp` — current nucleation enum + parser
  - `io/paraview_output.hpp` — `SetFaultParamsBP5`, `SetTotalRunTime`, `AdaptiveSchedule`
- Domain context: `miniapps/seas/CLAUDE.md`, `seas-mfem-safs/CLAUDE.md`

## Findings

### [R-001] CRITICAL `PLAN §Phase N Interfaces, Detailed Requirements #3-#4` — `dof_basis` storage layout is wrong by a factor of 3 × transposition

**Category:** BUG

**Description:**
The plan documents `dof_basis` as `Height() == 3 * N`, `Width() == 3` ("3x3 sub-block per DOF: row 0 = normal, row 1 = dip, row 2 = strike"). The **actual** layout in `drivers/spatial_dyn_driver.cpp:228, 257, 275-277` is `dof_basis(9, N)` — 9 ROWS × N COLUMNS, with per-DOF i:
- `dof_basis(0..2, i)` = normal (3 components)
- `dof_basis(3..5, i)` = tangent1 = dip (3 components)
- `dof_basis(6..8, i)` = tangent2 = strike (3 components)

The plan's resolver pseudocode + `GaussianFactorFaceLocal` callers would index the wrong memory and produce garbage F(r) values — every nucleation amplitude would be wrong, but the dry-run would not detect it because the resolver itself would not crash.

**Trigger:**
Implementer follows Phase N Detailed Req. #4 literally: "for each [DOF] build the basis row from `dof_basis` (3x3 sub-block per DOF: row 0 = normal, row 1 = dip, row 2 = strike)".

**Actual behavior:**
Resolver reads from `dof_basis.GetRow(3*i + 1)` (or similar) which is row 1, NOT dip. Garbage Gaussian decay.

**Expected behavior:**
Resolver reads from column `i` of `dof_basis`, offset by 3 for dip and 6 for strike. MFEM `DenseMatrix` is column-major, so `&dof_basis(0, i)` is a contiguous 9-element pointer for column i; `&dof_basis(3, i)` is contiguous dip; `&dof_basis(6, i)` is contiguous strike.

**Suggested fix:**

In the plan's Phase N Interfaces block, replace the `ResolveGradualOverstress` signature comment + Detailed Requirement #4 with the correct layout:

```diff
 /// `dof_coords_3d.Size() == 3 * N`; `dof_basis.Height() == 3 * N`,
 /// `dof_basis.Width() == 3` (basis is rows of `[n; t1; t2]` per DOF;
 /// matches `BuildPerDOFFaultTables` layout).
+/// `dof_coords_3d.Size() == 3 * N`; `dof_basis.Height() == 9`,
+/// `dof_basis.Width() == N`.  Per DOF i (matching driver:228, 257,
+/// 275-277): `dof_basis(0..2, i)` = normal, `dof_basis(3..5, i)` =
+/// tangent1 = dip, `dof_basis(6..8, i)` = tangent2 = strike.  MFEM
+/// DenseMatrix is column-major, so `&dof_basis(0, i)`, `&dof_basis(3, i)`,
+/// `&dof_basis(6, i)` are contiguous 3-vector pointers usable by
+/// GaussianFactorFaceLocal.
```

And in Detailed Requirement #4, replace the basis-row extraction language with:

```diff
-4. **`ResolveGradualOverstress` body** — iterate the N DOFs; for each
-   build the basis row from `dof_basis` (3x3 sub-block per DOF: row 0 =
-   normal, row 1 = dip, row 2 = strike — confirm layout from
-   `BuildPerDOFFaultTables` in `dynamic/spatial_setup.hpp`); call
-   `GaussianFactorFaceLocal`; ...
+4. **`ResolveGradualOverstress` body** — iterate the N DOFs; for each
+   pass `&dof_basis(3, i)` (dip) and `&dof_basis(6, i)` (strike) to
+   `GaussianFactorFaceLocal`.  The layout is `(9, N)` column-major
+   per `drivers/spatial_dyn_driver.cpp:228, 257, 275-277`: rows 0..2
+   = normal, rows 3..5 = dip, rows 6..8 = strike.  ...
```

**Test case:**
```cpp
TEST(spatial_nucleation, R001_basis_layout_centre_dof) {
   // 1-DOF fixture: centre at (10, 0, 0); basis: dip = +y, strike = +z.
   Vector dof_coords_3d(3); dof_coords_3d(0) = 10.0; dof_coords_3d(1) = 0.0; dof_coords_3d(2) = 0.0;
   DenseMatrix dof_basis(9, 1);
   dof_basis = 0.0;
   dof_basis(0, 0) = 1.0;   // normal +x
   dof_basis(4, 0) = 1.0;   // dip   +y
   dof_basis(8, 0) = 1.0;   // strike +z
   GradualOverstressSpec spec;
   spec.center_x_m = 10.0; spec.center_y_m = 0.0; spec.center_z_m = 0.0;
   spec.radius_dip_m = spec.radius_strike_m = 1000.0;
   spec.delta_tau_strike_pa = 25.0e6;
   spec.T_nuc_s = 1.0; spec.t0_smooth_s = 0.5;
   auto p = ResolveGradualOverstress(spec, /*enabled=*/true, dof_coords_3d, dof_basis);
   ASSERT_NEAR(p.radial(0), 1.0, 1e-12);                   // F(0) = 1
   ASSERT_NEAR(p.amplitude_strike(0), 25.0e6, 1e-6);       // F(0)·Δτ_strike
   ASSERT_NEAR(p.amplitude_dip(0), 0.0, 1e-12);             // Δτ_dip = 0
}
```

---

### [R-002] CRITICAL `PLAN §Phase D Acceptance T-D02` — Wrong field used to trigger outside-asperity equilibrium failure

**Category:** BUG (TEST INVALID)

**Description:**
T-D02 says: "a TOML with `delta_tau_strike_pa = 1.0e10` (absurdly large) AND `mu_s_default = 0.6` aborts at `--dry-run` with a `max OUTSIDE asperity = ...` message that includes the ratio."

But `delta_tau_strike_pa` is the **nucleation** amplitude (target value for the `tau1_nuc / tau2_nuc` accumulator). The outside-asperity equilibrium ratio uses `tau_pre` from the `[stress]` block, not `tau_nuc`. Inflating `delta_tau_strike_pa` would inflate the nucleation patch's overshoot (testing the `INSIDE` budget check from Phase D §3, which is T-D03) but would NOT change the outside-asperity ratio.

The test as written cannot trigger the gate it claims to test. The check at Phase D §2 reads `tau_pre_per_dof(i)` from `geom.GetTauPre()`, not from `nuc_params.amplitude_*(i)`.

**Trigger:**
Implementer writes T-D02 verbatim; the test fails to abort because the equilibrium-gate code path is never exercised.

**Actual behavior:**
With `delta_tau_strike_pa = 1.0e10` and a normal `[stress]` block, `tau_pre` magnitudes are normal (a few MPa). `mu_s_default = 0.6 · σ_n_eff` of order 30 MPa easily dominates → outside ratio < 1 → no abort.

**Expected behavior:**
T-D02 should use a `[stress] kind = "constant_tensor"` with extreme shear values so the projected per-DOF `|τ_pre|` exceeds `μ_s · σ_n_eff` outside the Gaussian's `3·radius_*` support.

**Suggested fix:**

```diff
-- [ ] **T-D02 — Pre-flight ABORTS on supercritical setup**: a TOML with
-  `delta_tau_strike_pa = 1.0e10` (absurdly large) AND
-  `mu_s_default = 0.6` aborts at `--dry-run` with a
-  `max OUTSIDE asperity = ...` message that includes the ratio.
+- [ ] **T-D02 — Pre-flight ABORTS on supercritical setup**: a TOML with
+  `[stress] kind = "constant_tensor", sigma_xy_pa = 5.0e8` (500 MPa pure
+  shear — physically unreasonable, will produce per-DOF `|τ_pre|` >
+  100 MPa everywhere) AND `mu_s_default = 0.4`, `mu_d_default = 0.3`
+  aborts at `--dry-run` with a `max OUTSIDE asperity = X.XX` message
+  where the ratio exceeds 1.  (Inflating `delta_tau_strike_pa` does
+  NOT trigger this gate — that flag controls the nucleation
+  amplitude, which is the INSIDE-asperity check T-D03 below.)
```

**Test case:**
```python
def test_R002_outside_asperity_supercritical_aborts():
    # Pre-stress 500 MPa shear, friction strength ~30 MPa => ratio ~17
    toml = """
        [meta] schema_version = 1; law = "slip_weakening"
        [stress] kind = "constant_tensor"
                 sigma_xx_pa = 0.0; sigma_yy_pa = 0.0; sigma_zz_pa = 0.0
                 sigma_xy_pa = 5.0e8; sigma_yz_pa = 0.0; sigma_xz_pa = 0.0
        [friction.slip_weakening]
                 mu_s_default = 0.4; mu_d_default = 0.3; d_c_default = 0.5
        # ... mandatory mesh/material/velocity/numerics/time/output ...
    """
    rc = subprocess.run(["./seas_spatial_dyn_driver", "--config", toml_path,
                         "--no-sidecar-material", "--dry-run", "--print-derived"])
    assert rc.returncode != 0
    assert "max OUTSIDE asperity" in rc.stderr
```

---

### [R-003] MODERATE `PLAN §Phase N Files to Modify (bullet on spatial_setup.hpp)` — `InitializeFaultDOFs_Spatial` arg-list change is ambiguous; impl path not chosen

**Category:** ASSUMPTION

**Description:**
The existing `InitializeFaultDOFs_Spatial<MeshT>` signature in `dynamic/spatial_setup.hpp:259-289, 358-388` requires `const Vector& T_forced_s, const Vector& t0_decay_s` as mandatory args (validated via `MFEM_VERIFY(T_forced_s.Size() == ndof, ...)`). The plan says:

> pass empty `T_forced_s` / `t0_decay_s` to `InitializeFaultDOFs_Spatial` (the LSW path consumes only LSW params when the friction-law is `LSW` not `LSW_ForcedRupture`), or remove those args from the call entirely if the signature changes (see below)

But the signature change is NOT specified anywhere "below". The validator at `spatial_setup.hpp:194-200` will abort on empty vectors (`MFEM_VERIFY(T_forced_s.Size() == ndof)`).

Two viable paths:
- **(a)** Add `InitializeFaultDOFs_Spatial` overload without `T_forced_s/t0_decay_s` (additive; preserves existing TPV path if any).
- **(b)** Pass dummy vectors initialized to `1.0e9` (the existing "never forced" sentinel) and `0.0`. The LSW dispatch (now always `FaultFrictionLaw::LSW`) never reads them, so the values are dormant.

**Trigger:**
Implementer reads Phase N "removed" list, deletes `ResolveForcedRupture` call, then has no `fr.T_forced_s / fr.t0_decay_s` to pass.

**Actual behavior:**
Compile error (`fr` undefined) or runtime abort from the size check.

**Expected behavior:**
Plan should pick one path. Recommend **(b)**: it's a 4-line driver change, no new overload, no API surface bloat. The dummy vectors are populated alongside the LSW resolver call.

**Suggested fix:**

Add to Phase N Detailed Requirement #7 (driver-side wiring):

```diff
+   // After step 11 (LSW resolver) and before step 14 (InitializeFaultDOFs_Spatial):
+   // Populate dummy "never forced" vectors so InitializeFaultDOFs_Spatial's
+   // size validator passes.  These vectors are inert under the LSW
+   // dispatch (FaultFrictionLaw::LSW, not LSW_ForcedRupture) — the
+   // EvaluateADER_LSW kernel never reads tau1_nuc / tau2_nuc as forced-
+   // rupture sentinels; they are only the per-sub-step accumulator
+   // channel that gradual_overstress writes into.
+   Vector dummy_T_forced(num_fault_total);  dummy_T_forced = 1.0e9;
+   Vector dummy_t0_decay(num_fault_total);  dummy_t0_decay = 0.0;
+
+   // Step 14 (InitializeFaultDOFs_Spatial) — use dummies in place of fr.*:
+   spatial::InitializeFaultDOFs_Spatial<ParMesh>(
+      dof_data, num_fault_total, dof_to_elem, material, pmesh,
+      lsw, geom.GetTauPre(), geom.sigma_n_per_dof(),
+      dummy_T_forced, dummy_t0_decay,
+      dof_ips);
```

And in Phase N Files to Modify, replace the `spatial_setup.hpp` bullet:

```diff
-- `miniapps/seas/dynamic/spatial_setup.hpp` (if any) — confirm
-  `InitializeFaultDOFs_Spatial`'s LSW overload does NOT require
-  `T_forced_s / t0_decay_s` arguments; if it does, add a third overload
-  without them OR pass dummy zero-sized vectors (driver convention).  ...
+- `miniapps/seas/dynamic/spatial_setup.hpp` — NO CHANGE.  The existing
+  signature requires `T_forced_s / t0_decay_s`; the driver supplies
+  dummy `(1.0e9, 0.0)` vectors (per Detailed Req. #7).  These are inert
+  under the LSW dispatch path that Phase N enforces.
```

**Test case:**
```cpp
TEST(spatial_dyn_driver, R003_dummy_forced_rupture_vectors_inert) {
   // Driver with gradual_overstress nucleation must run end-to-end with
   // dummy T_forced = 1e9 / t0_decay = 0 passed to InitializeFaultDOFs_Spatial.
   // Assert DOFData::tau{1,2}_nuc are not corrupted by the dummies.
   auto exit_code = run_driver_dry_run({"--config", EXAMPLE_TOML,
                                        "--no-sidecar-material", "--dry-run"});
   ASSERT_EQ(exit_code, 0);
}
```

---

### [R-004] MODERATE `PLAN §Phase D Detailed Req. #4` — σ_1 eigensolver vague; MFEM API not cited; azimuth convention test depends on unspecified atan2 ordering

**Category:** QUALITY (BLOCKS IMPLEMENTATION)

**Description:**
Plan says "use the stable `SymmetricEigenvalueDecomposition` helper if MFEM provides one, else hand-roll using the Smith 1961 trigonometric formula". MFEM's actual API is `mfem::DenseMatrix::Eigensystem(Vector& ev, DenseMatrix& evec)` (or `CalcEigenvalues`). The implementer must not "hand-roll Smith 1961" if MFEM has it.

Additionally, the azimuth convention is "radians clockwise from north in the horizontal plane" — this requires `atan2(eigvec_x, eigvec_y)` (note: x then y, NOT y then x). T-D06 expects 45° for σ_1 along (+x, +y, 0); a careless `atan2(eigvec_y, eigvec_x)` would also return 45° for this specific case but would be WRONG for σ_1 along (-x, +y, 0) (which should be azimuth 315°, not 135°).

**Trigger:**
Implementer hand-rolls Smith 1961 (slow + bug-prone) OR uses wrong atan2 ordering for azimuth.

**Actual behavior:**
Untested code path for non-NE σ_1; possible silent miscomputation.

**Expected behavior:**
Plan cites the actual MFEM eigensolver API and shows the exact atan2 formula.

**Suggested fix:**

```diff
-4. **σ_1 azimuth + max-shear direction** for the input stress tensor.
-   For `kind = constant_tensor`, compute the principal stresses
-   `σ_1 >= σ_2 >= σ_3` from the 6-component tensor via the cubic
-   characteristic equation (closed form for symmetric 3×3 — use the
-   stable `SymmetricEigenvalueDecomposition` helper if MFEM provides one,
-   else hand-roll using the Smith 1961 trigonometric formula). ...
+4. **σ_1 azimuth + max-shear direction** for the input stress tensor.
+   For `kind = constant_tensor`, build the symmetric 3x3 from the six
+   `[stress].sigma_*_pa` keys and call MFEM's eigensolver:
+   ```cpp
+   mfem::DenseMatrix S(3,3);
+   S(0,0)=σ_xx; S(1,1)=σ_yy; S(2,2)=σ_zz;
+   S(0,1)=S(1,0)=σ_xy; S(1,2)=S(2,1)=σ_yz; S(0,2)=S(2,0)=σ_xz;
+   mfem::Vector eigvals(3); mfem::DenseMatrix eigvecs(3,3);
+   S.Eigensystem(eigvals, eigvecs);   // ascending eigenvalues
+   const int i1 = 2;                  // σ_1 = max (largest = most compressive)
+   const real_t v1x = eigvecs(0, i1);
+   const real_t v1y = eigvecs(1, i1);
+   const real_t v1z = eigvecs(2, i1);
+   ```
+   Azimuth (clockwise from north in EAST-NORTH frame; x=E, y=N):
+   ```cpp
+   const real_t azimuth_rad = std::atan2(v1x, v1y);   // NB: (x, y) ORDER
+   real_t azimuth_deg = azimuth_rad * 180.0 / M_PI;
+   if (azimuth_deg < 0.0) { azimuth_deg += 360.0; }   // wrap to [0, 360)
+   ```
+   Plunge (degrees from horizontal):
+   ```cpp
+   const real_t v1h = std::sqrt(v1x*v1x + v1y*v1y);
+   const real_t plunge_rad = std::atan2(-v1z, v1h);   // +ve plunge = down
+   ```
+   For `kind = sidecar_hdf5`, ...
```

**Test case:**
```cpp
TEST(spatial_print_derived, R004_azimuth_convention_negative_x) {
   // sigma_1 along (-1, +1, 0)/sqrt(2) => azimuth 315 deg (NW), NOT 135 deg.
   real_t azimuth_deg = ComputeSigma1AzimuthDeg(
      /*sxx=*/ 50e6, /*syy=*/ 50e6, /*szz=*/ 0,
      /*sxy=*/ -25e6, /*syz=*/ 0,    /*sxz=*/ 0);
   EXPECT_NEAR(azimuth_deg, 315.0, 0.1);
}
```

---

### [R-005] MODERATE `PLAN §Phase Z Files to Create` — `verify_spatial_dyn_smoke_safs.py` "stdlib + h5py only" is misleading; h5py is third-party

**Category:** ASSUMPTION

**Description:**
Plan says "stdlib + h5py only". h5py is NOT in the Python stdlib. The `pythonenv` conda env (per `CLAUDE.md` Environment Setup) has h5py for the mesh / velocity / stress pipelines, but the runbook recommends `mfem-dev` for building/running which may or may not have h5py.

**Trigger:**
User runs `python verify_spatial_dyn_smoke_safs.py ...` in `mfem-dev` env; ImportError.

**Actual behavior:**
Verifier crashes with `ModuleNotFoundError: No module named 'h5py'`.

**Expected behavior:**
Plan explicitly states env requirement and documents fallback.

**Suggested fix:**

```diff
-3. **Verifier `verify_spatial_dyn_smoke_safs.py`** — stdlib + h5py only.
+3. **Verifier `verify_spatial_dyn_smoke_safs.py`** — Python 3 stdlib +
+   **third-party h5py** (h5py is NOT in stdlib).  Run under the
+   `pythonenv` conda env (per CLAUDE.md Environment Setup) which has
+   h5py for the mesh/velocity/stress pipelines.  At the top of the
+   script, fail-fast with a helpful message if the import fails:
+   ```python
+   try:
+       import h5py
+   except ImportError:
+       sys.exit("verify_spatial_dyn_smoke_safs.py requires h5py.  "
+                "Run under `conda activate pythonenv` (per "
+                "miniapps/seas/CLAUDE.md Environment Setup).")
+   ```
```

**Test case:**
```bash
# Verify graceful fail when h5py absent
$ python -c "import sys; sys.modules['h5py'] = None; \
             exec(open('verify_spatial_dyn_smoke_safs.py').read())"
# Should exit with the helpful message, not a stack trace.
```

---

### [R-006] LOW `PLAN §Phase Z Files to Create (production sbatch)` — node-count "derivation" language implies impossible runtime call from shell

**Category:** QUALITY

**Description:**
Plan says: "do NOT hardcode if you can derive it from `pmesh.GetGlobalNE()` divided by a target elements-per-rank (≈ 2500). Document the calculation in the sbatch header comment." But `pmesh.GetGlobalNE()` is a C++/MFEM call — shell scripts cannot invoke it.

The intent is clearly **offline computation**: the implementer runs the driver once (or inspects `[derived]` log output) to learn the global element count, divides by 2500, and bakes the node count into the sbatch with a header comment documenting the source.

**Trigger:**
Implementer reads the sbatch instruction and tries to find a way to compute the count in shell (impossible).

**Actual behavior:**
Confusion; possible attempt to invoke `seas_spatial_dyn_driver --dry-run` from inside the sbatch (which conflates the dev-queue submission with a pre-flight).

**Expected behavior:**
Plan says explicitly "offline".

**Suggested fix:**

```diff
-- `miniapps/seas/jobs/safs/spatial_dyn_production_normal_24hr_safs.sbatch`
-  — Frontera normal queue, 32 nodes × 50 ranks = 1600 ranks, 24 hr wall,
-  `--tfinal 12s`, V1 checkpoint every 10000 steps.  The 32-node count is
-  derived from the SAFS 1000 m mesh element count (≈ 4 M tets) and the
-  per-rank memory budget (≈ 2.5 GB/rank for the wave state under
-  constant material; documented in the sbatch header).  Do NOT hardcode
-  if you can derive it from `pmesh.GetGlobalNE()` divided by a target
-  elements-per-rank (≈ 2500).  Document the calculation in the sbatch
-  header comment.
+- `miniapps/seas/jobs/safs/spatial_dyn_production_normal_24hr_safs.sbatch`
+  — Frontera normal queue, 32 nodes × 50 ranks = 1600 ranks, 24 hr wall,
+  `--tfinal 12s`, V1 checkpoint every 10000 steps.  Node count derived
+  OFFLINE: implementer runs `seas_spatial_dyn_driver --dry-run
+  --print-derived` on a single Frontera dev node FIRST, reads
+  `num_global_elements = X` from the [derived] log, computes `nodes =
+  ceil(X / (50 ranks/node · 2500 elem/rank))`, bakes the result into the
+  sbatch header with a comment citing the dry-run log line.  The
+  2500-elem/rank target is empirical (from BP5/TPV104 scaling); revisit
+  after the first production run if step time is dominated by MPI
+  communication.
```

---

### [R-007] LOW `PLAN §Constraints "Numerical constraints" — "Accumulator zero-bound clamp"` — Threshold direction is unstated for `t > T_nuc_s` early-return

**Category:** EDGE_CASE

**Description:**
Phase N Detailed Req. #5 says the accumulator should "early-return when `t_substep_end >= T_nuc_s + dt_substep`". But the smoothStep is at exactly 1 for any `t >= T_nuc_s` (since `t0_smooth_s <= T_nuc_s` — wait, actually `t0_smooth_s` IS the smoothStep `t0` parameter, NOT `T_nuc_s`).

There's a subtle ambiguity: the `T_nuc_s` field's role is documented as "ramp completes at t = T_nuc_s" but `SmoothStep(t, t0_smooth_s)` returns 1 for `t >= t0_smooth_s`, not for `t >= T_nuc_s`. If the schema is `t0_smooth_s = T_nuc_s / 2` (the example), then the ramp actually completes at `t = t0_smooth_s = 0.5s`, not at `T_nuc_s = 1.0s`. The "T_nuc_s" field then becomes a no-op duration label.

The plan needs to clarify: which of `T_nuc_s` and `t0_smooth_s` is the smoothStep midpoint vs the ramp duration?

Convention from `tpv104_nucleation.hpp:115`: `SmoothStep_TPV104(t_substep_end, TPV104Params::nuc_T)` — passes `nuc_T` as `t0`. So `nuc_T` IS the smoothStep argument; the ramp completes at `t = nuc_T`.

In the plan's TOML schema, the user supplies `T_nuc_s = 1.0` and `t0_smooth_s = 0.5`. But the accumulator must pass ONE of these to `SmoothStepIncrement`. The plan's `ApplyGradualOverstressIncrement` signature takes BOTH (`T_nuc_s` and `t0_smooth_s`) — but the smoothStep only consumes one of them.

**Trigger:**
Implementer writes `SmoothStepIncrement(t, dt, T_nuc_s)` (using T_nuc_s) but EXAMPLE TOML has `T_nuc_s = 1.0, t0_smooth_s = 0.5`; ramp actually completes at t = 1.0s, not at t = 0.5s. OR vice versa.

**Actual behavior:**
Ambiguous; depends on which field the implementer passes.

**Expected behavior:**
Plan picks one. Recommend: drop `t0_smooth_s` from the schema entirely (it's redundant), let `T_nuc_s` be the smoothStep `t0` parameter (matches TPV104 convention).

**Suggested fix:**

In the schema doc + EXAMPLE TOML + parser + plan:
```diff
-`T_nuc_s`             | s    | —       | required; `> 0`; recommended `≤ tfinal/10`|
-`t0_smooth_s`         | s    | —       | required; `> 0`; recommended `≈ T_nuc_s/2`|
+`T_nuc_s`             | s    | —       | required; `> 0`; recommended `≤ tfinal/10`.  This is the smoothStep `t0` parameter — the ramp `SmoothStep(t, T_nuc_s)` goes from 0 at t=0 to 1 at t=T_nuc_s.
```

And in `ApplyGradualOverstressIncrement` signature, drop the redundant `t0_smooth_s` parameter (pass only `T_nuc_s`):
```diff
 void ApplyGradualOverstressIncrement(
    std::vector<DOFData>&                  dof_data,
    const GradualOverstressPerDOFParams&   params,
    real_t                                 T_nuc_s,
-   real_t                                 t0_smooth_s,
    real_t                                 t_substep_end,
    real_t                                 dt_substep);
```

In the smoothStep call inside the accumulator body: `SmoothStepIncrement(t_substep_end, dt_substep, T_nuc_s)`.

This is a SIMPLIFICATION; the plan's two-field schema was over-parameterized and invites the bug above.

---

## Summary

- Critical issues: 2 (R-001, R-002)
- Moderate issues: 3 (R-003, R-004, R-005)
- Low issues: 2 (R-006, R-007)
- Plan compliance: PARTIAL — plan is implementable in principle but ships with at least one ambient-true layout bug (R-001) and one invalid test (R-002) that would mislead the implementer.
- Verdict: **PASS WITH FIXES** — apply R-001, R-002, R-003, R-007 BEFORE handing the plan to `/code-implement`; R-004, R-005, R-006 can be applied in the same fix round but are not blocking.

## Unreviewed Areas

- The parity plan referenced by Stream C was NOT re-reviewed in this round (only the deltas this plan introduces against it). The original parity plan's review history (round-by-round documents in `debug_document/`) is the authoritative source for parity plan internals.
- The `verify_spatial_dyn_smoke_safs.py` fixture (`.vtkhdf` file) referenced by T-Z03 is not yet committed; the test cannot be exercised until the fixture exists.
- σ_1 azimuth convention for `kind = "sidecar_hdf5"` is deferred to a `[derived] TBD` line per Phase D §4; no review of that deferred path here.
