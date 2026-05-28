# Code Review: Phase 6 (config schema + D3.2 fault-local prestress) — increments 1+2 — 2026-05-27

> NOTE: prior unrelated reviews (`REVIEW.md` TPV205, `REVIEW_phase4_drename_*`,
> `REVIEW_phase5_*`) are left intact. Point `/code-fix` at THIS file.

## Review Scope
- **Plan:** `miniapps/seas/document/fullelasticity_dev/PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.md`, §Phase 6 (req 2, D3.2 fault-local prestress).
- **Worktree:** `/Users/chunhuizhao/projects/seas-mfem-safs` (branch `safs`). `/code-fix` operates here.
- **Files reviewed (Phase-6 increments 1+2):**
  - `fault/fault_geometry.hpp` — `ComputeParamsFaultLocal` (committed `2f2961d`)
  - `spatial/code/spatial_friction.hpp` — `StressSourceKind::FaultLocalPrestress` + `StressSpec` fields (committed `2f2961d`)
  - `spatial/code/spatial_friction.cpp` — `parse_stress_kind` + `parse_stress` branch (uncommitted)
  - `tests/unit/test_compute_safs_params.cpp` — T-66 (uncommitted)
  - `tests/unit/test_spatial_friction_config.cpp` — FLP-1/FLP-2 (uncommitted)
- **Domain context:** root + miniapps/seas `CLAUDE.md` (sign conventions: σ_n>0 compression, tau1=dip/tau2=strike); the `ProjectFaultPreStress` projection path (`io/field_coefficient.cpp`); `InitializeFaultDOFs_Spatial` (the downstream `tau_pre_`/`sigma_n_per_dof_` consumer).

## Independent verification performed
- `tau_pre_` layout `(2i)=dip→tau1_0`, `(2i+1)=strike→tau2_0`, `sigma_n_per_dof_(i)` confirmed against `field_coefficient.cpp:504-505` and `spatial_setup.hpp:90-91`; T-66 verifies it at runtime.
- `parse_stress` restructure (`if/else-if/else`) leaves the `ConstantTensor`/`SidecarHDF5` branches byte-unchanged — `spatial_friction_config` **66/66** (no regression to existing config parsing).
- Tests this session: `compute_safs_params` **20/20** (T-66 D3.2 sign), `spatial_friction_config` **66/66** (FLP-1 parse + FLP-2 conflict-abort).

**No CRITICAL production bug found.** The seeding is bit-correct and `ComputeParamsFaultLocal` is not yet wired into the driver (it is tested in isolation; production wiring is a later Phase-6/8 increment). Findings are a defensive-validation gap, a documented partial implementation, and two parse-robustness nits.

## Findings

### [R-001] MODERATE [fault/fault_geometry.hpp:ComputeParamsFaultLocal] — effective normal stress σ_n−P_p is not guarded > 0 (asymmetric with the projection path's floor)

**Category:** EDGE_CASE / ASSUMPTION

**Description:**
`sigma_n_eff = sigma_n - P_p` is written into `sigma_n_per_dof_(i)` with **no check that it is positive**. The friction solve requires a compressive (σ_n > 0) effective normal stress. The Cauchy-projection path (`FieldProjector::ProjectFaultPreStress`) clamps with `min_sigma_n_pa`; the fault-local path has no equivalent guard. The parse branch checks `sigma_n_pa > 0` but NOT `sigma_n_pa − P_p > 0` (P_p is parsed independently into `cfg.stress.pore_pressure`), so a config with `P_p ≥ sigma_n_pa` would silently seed a tensile/zero effective normal stress once the driver wires `ComputeParamsFaultLocal(..., P_p)`.

**Trigger:** `geom.ComputeParamsFaultLocal(tau_strike, tau_dip, sigma_n, P_p)` with `P_p >= sigma_n` (e.g., `sigma_n=120e6, P_p=120e6` → `sigma_n_eff=0`; `P_p=130e6` → negative).

**Actual behavior:** `sigma_n_per_dof_(i) = sigma_n_eff ≤ 0` is seeded silently; downstream the friction solve operates on a non-physical (tensile) normal stress.

**Expected behavior:** abort with a precise message (or clamp), mirroring the projection path's `min_sigma_n_pa` protection.

**Suggested fix:** add a guard after computing `sigma_n_eff`.
```diff
       const real_t sigma_n_eff = sigma_n - P_p;
+      MFEM_VERIFY(sigma_n_eff > 0.0,
+                  "FaultGeometry::ComputeParamsFaultLocal: effective normal "
+                  "stress sigma_n - P_p must be > 0 (compression positive); got "
+                  "sigma_n=" << sigma_n << ", P_p=" << P_p
+                  << " -> sigma_n_eff=" << sigma_n_eff);
       sigma_n_per_dof_.SetSize(num_fault_dofs_);
```

**Test case:**
```python
def test_R001_fault_local_tensile_effective_normal_aborts():
    # ComputeParamsFaultLocal with P_p >= sigma_n must abort (effective
    # normal stress would be non-positive), not silently seed tension.
    geom = make_bp5_fault_geometry()
    assert aborts(lambda: geom.ComputeParamsFaultLocal(
        tau_strike=75e6, tau_dip=0.0, sigma_n=120e6, P_p=120e6))
    # control: P_p < sigma_n seeds sigma_n_eff = 104e6 without aborting
    geom.ComputeParamsFaultLocal(75e6, 0.0, 120e6, 16e6)
    assert geom.sigma_n_per_dof()[0] == 104e6
```

---

### [R-002] MODERATE [DEVIATION] [spatial/code/spatial_stress.*] — req 2 only partially implemented: `FaultLocalPrestressSource` holder + rectangular `tau_strike` patches are missing

**Category:** DEVIATION

**Description:**
Plan §Phase 6 req 2 specifies "Add `FaultLocalPrestressSource` (`spatial_stress.{hpp,cpp}`) with constant background `{tau_strike_pa, tau_dip_pa, sigma_n_pa}` **+ optional rectangular-`tau_strike` patches** (`center_*`/`half_*`, last-match-wins)." Increments 1+2 implemented the **uniform-background** path (`ComputeParamsFaultLocal` + parse), but **not** the `FaultLocalPrestressSource` class nor the rectangular patches. The implementer documented this deferral, but per plan §2-D3.2 **TPV205 requires positive `tau_strike` patches** (background +70 MPa; central +81.6, left +78, right +62 MPa) — so TPV205 cannot yet express its pre-stress through this path, which blocks the Phase-8 TPV205 case.

**Trigger:** Attempting to configure a TPV205 fault-local pre-stress with spatially-varying `tau_strike` patches.

**Actual behavior:** Only a spatially-uniform fault-local background is supported.

**Expected behavior:** `FaultLocalPrestressSource` holds the background + a list of rectangular patches; a per-DOF override pass (last-match-wins) layers patch `tau_strike` onto the uniform seed.

**Suggested fix:** add the patch struct + source + a per-DOF override pass (next Phase-6 increment). Sketch:
```diff
+ // spatial_friction.hpp StressSpec (or a dedicated FaultLocalPrestressSpec):
+ struct FaultLocalPatch {
+    real_t center_x_m = NaN, center_y_m = NaN, center_z_m = NaN;
+    real_t half_x_m = +inf, half_y_m = +inf, half_z_m = +inf;
+    real_t tau_strike_pa = 0.0;   // overrides the background inside the box
+ };
+ std::vector<FaultLocalPatch> fault_local_patches;   // last-match-wins
```
then a `FaultGeometry` overload (or post-pass) that, for each fault DOF, starts from the uniform `tau_strike` and applies the last matching patch.

**Test case:**
```python
def test_R002_tpv205_patch_overrides_background():
    # central patch (+81.6 MPa) overrides the +70 MPa background inside its box;
    # a DOF outside all patches keeps +70 MPa.
    spec = parse_fault_local_with_patches(...)
    seed = compute_fault_local_with_patches(geom, spec)
    assert strike_prestress_at(seed, inside_central_patch) == 81.6e6
    assert strike_prestress_at(seed, outside_all_patches)  == 70.0e6
```

---

### [R-003] LOW [spatial/code/spatial_friction.cpp:parse_stress] — a forgotten `tau_strike_pa` silently yields a locked fault

**Category:** EDGE_CASE

**Description:**
The `FaultLocalPrestress` parse branch requires only `sigma_n_pa`; `tau_strike_pa` and `tau_dip_pa` default to 0.0. A config that omits `tau_strike_pa` therefore parses successfully but seeds zero shear pre-stress → the fault never reaches frictional failure (locked). For a benchmark whose whole point is a shear-driven rupture, this is a silent foot-gun.

**Trigger:** A `[stress] kind="fault_local_prestress"` block with `sigma_n_pa` but no `tau_strike_pa`.

**Actual behavior:** `tau_strike_pa = 0` silently; no rupture.

**Expected behavior:** require `tau_strike_pa` (it is the essential driver), or emit a rank-0 warning when it is absent/zero.

**Suggested fix:**
```diff
       MFEM_VERIFY(s.contains("sigma_n_pa"),
                   "[stress] kind=\"fault_local_prestress\" requires "
                   "sigma_n_pa (compression POSITIVE)");
+      MFEM_VERIFY(s.contains("tau_strike_pa"),
+                  "[stress] kind=\"fault_local_prestress\" requires "
+                  "tau_strike_pa (the shear driver; right-lateral POSITIVE)");
```
(If a zero-shear fault-local case is ever legitimate, downgrade to a rank-0 warning instead of an abort.)

---

### [R-004] LOW [spatial/code/spatial_friction.cpp:parse_stress] — asymmetric mutual-exclusion: ConstantTensor/SidecarHDF5 branches do not reject stray fault-local keys

**Category:** QUALITY

**Description:**
The new `FaultLocalPrestress` branch rejects stray Cauchy `sigma_*_pa` and `sidecar_path`, but the `ConstantTensor` and `SidecarHDF5` branches do **not** reject stray `tau_strike_pa`/`tau_dip_pa`/`sigma_n_pa`. A `constant_tensor` config that mistakenly sets `tau_strike_pa` would silently ignore it — inconsistent with the carefully-guarded fault-local branch.

**Trigger:** `kind="constant_tensor"` (or `"sidecar_hdf5"`) with a stray `tau_strike_pa`.

**Actual behavior:** the fault-local key is silently ignored.

**Expected behavior:** reject the unknown-for-this-kind keys, symmetric with the fault-local branch.

**Suggested fix:** in the `ConstantTensor` and `SidecarHDF5` branches, add
```diff
+      MFEM_VERIFY(!s.contains("tau_strike_pa") && !s.contains("tau_dip_pa")
+                  && !s.contains("sigma_n_pa"),
+                  "[stress] kind=\"constant_tensor\"/\"sidecar_hdf5\" must NOT "
+                  "set tau_strike_pa/tau_dip_pa/sigma_n_pa (those are for "
+                  "kind=\"fault_local_prestress\")");
```

---

## Summary
- Critical issues: 0
- Moderate issues: 2 (R-001 σ_n−P_p guard; R-002 partial req-2 / missing patches)
- Low issues: 2 (R-003 missing `tau_strike_pa` requirement; R-004 asymmetric key rejection)
- **Plan compliance:** PARTIAL — the uniform fault-local path (`ComputeParamsFaultLocal` + parse + sign tests) is complete and correct; the `FaultLocalPrestressSource` holder + rectangular patches half of req 2 is deferred (R-002), so TPV205's patched pre-stress is not yet expressible.
- **Verdict:** PASS WITH FIXES — no correctness bug on the tested path; apply R-001 (defensive guard) before the driver wires `ComputeParamsFaultLocal`, and complete R-002 (patches) before the Phase-8 TPV205 case. R-003/R-004 are cheap robustness hardening.

## Unreviewed Areas
- **Driver wiring of `ComputeParamsFaultLocal`** — not yet implemented (no production caller passes `cfg.stress.{tau_strike_pa,tau_dip_pa,sigma_n_pa}` + `P_p`); R-001 becomes live when it is.
- **Remaining Phase-6 requirements** (1, 3, 4, 5, 6, 7) — not yet implemented; out of scope for this review.
- **The 6 sibling TPV/SAFS TOMLs** — `spatial_friction_config` 66/66 covers the parser, but no end-to-end parse of an actual TPV TOML (those need reqs 1/3/4/5/6 to exist first).
