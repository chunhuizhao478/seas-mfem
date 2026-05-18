# Code Review: spatial_dynamic_rupture_plan.md — Adversarial Pre-Implementation Audit (2026-05-18)

## Review Scope

- **Plan reviewed:** `safs/project_7.0_alternative/document/spatial_dynamic_rupture_plan.md`
  (the PDF is the rendered twin of the .md; the .md is authoritative).
- **User objective (drives severity):**
  > "Dynamic-rupture simulation similar to **tpv104, tpv102, tpv205**, with
  > `safs_fault_box_nwcut_1000m_lcfar3000.msh`; spatial varying velocity (compare
  > `cvm_s4.26.m01` vs `multiscale_statewise_cvm`); for testing assume **constant
  > background stress + rotation onto fault**, later CSM; **focus on linear
  > slip-weakening** (TPV205-style, parameters in
  > `friction/slip-weakening/geoffrey2010.md`); later also rate-and-state. Need
  > functions general enough to take different datasets for quick switch."
- **User decisions accepted before this review** (will be baked into the
  `/code-plan` revision after this review lands):
  - **D-1** Constant-tensor stress is a **first-class TOML mode**
    `[stress] kind = "constant_tensor"`.
  - **D-2** Heterogeneous-material prereq is handled with a **hard Preconditions
    gate** (`static_assert` build-time guard); Phase 4 is blocked until
    `heterogeneous_material_plan.md` Phases 1+3 merge.
  - **D-3** Friction validator is **relaxed** (drop `mu_s < 1`; keep
    `0 < mu_d < mu_s`, `d_c > 0`); accept `d_o` as a deprecation-tagged alias
    for `d_c` so `geoffrey2010.md` loads verbatim.
- **Companion / sibling plans referenced for shared phases:**
  `spatial_quasi_dynamic_plan.md` (Phases 0/1/2/3 are claimed shared);
  `heterogeneous_material_plan.md` (prereq).
- **Code surface audited for API/contract claims the plan makes:**
  - `dynamic/wave_operator.hpp/.cpp/.inl` (`WaveOperator` template, ctor,
    `SetFaultFlux`, `SetFaultFrictionLaw`, `FaultFrictionLaw` enum,
    `MixedFluxMode`, `ComputeMaxDt`, `GodunovFlux` member)
  - `dynamic/godunov_flux.hpp` (class shape; verify absence of
    `GodunovFluxPool`)
  - `dynamic/fault_face_flux.hpp/.cpp` (`FaultFaceFlux` class, `DOFData`
    struct, `EvaluateADER`, `EvaluateADER_LSW`, impedance fields)
  - `dynamic/tpv205_setup.hpp`, `tpv104_setup.hpp` (the free-function
    `InitializeFaultDOFs_TPV*` per-driver init pattern that actually populates
    `DOFData`)
  - `dynamic/tpv102_substep_iterator.hpp`, `tpv104_substep_iterator.cpp`
  - `io/tpv104_checkpoint.hpp` (`Tpv104CheckpointFilename`, per-rank `.txt`
    layout via `CheckpointFilename`)
  - `io/checkpoint.hpp` (V1 free-function `WriteCheckpoint` / `ReadCheckpoint`,
    per-rank `_checkpoint_r<rank>.txt`)
  - `io/data_field_3d.hpp`, `io/stress_field_3d.hpp`
  - `io/material_coefficients.hpp/.cpp` (`{Lambda,Mu,Rho}FromSidecar`)
  - `dynamic/heterogeneous_material.hpp` (`MaterialField`, `Mode::Coefficient`,
    `MakeConstant`, `MakeCoefficient`, `EvalAt`)
  - `domain/elasticity_operator.hpp`, `domain/boundary_config.hpp`
  - `fault/fault_geometry.hpp`, `fault/fault_geometry_safs.inl`
    (BP5 ctor + `ComputeSAFSParams`, `HasSAFSParams`, `sigma_n_per_dof`,
    `GetTauPre`, `fault_dof_coords_3d`, `fault_dof_basis`, `NumZeroNormalFallbacks`)
  - `drivers/tpv104_driver.cpp` (CLI shape, `weakly_canonical` safety check,
    `FaultFaceFlux fault_flux(rho, cp, cs)` construction pattern,
    `wave.SetFaultFlux(&fault_flux)`)
  - `drivers/tpv205_driver.cpp` (LSW pattern, `SetFaultFrictionLaw(LSW)`)
  - `meshing/results/msh/safs_fault_box_nwcut_1000m_lcfar3000.msh`
    (Physical Surface tags 101/102/103/104 + UTM bbox)
  - `velocity/results/{cvmh,cvm_s4.26.m01,multiscale_statewise_cvm}/`
    (actual filenames are `velocity_safs.h5`, no mesh-tag suffix)
  - `friction/slip-weakening/geoffrey2010.md` (actual `mu_s=1.1, mu_d=0.5,
    d_o=0.5m`)
- **Domain context:** `seas-mfem-safs/CLAUDE.md`, `miniapps/seas/CLAUDE.md`,
  `velocity/README.md`, `stress/README.md`, project README.
- **Companion review:** `spatial_quasi_dynamic_review.md` (discarded per user
  decision; only its R-003, R-004, R-006, R-008 carry over because Phases
  0/1/2/3 are shared text).

This is a **pre-implementation review**: no `spatial_*` source has been
written. Pass 1 = "does the plan match the user's stated objective?", Pass 2 =
"if implemented as written, will the code work against the actual data and the
actual codebase?", Pass 3 = "does the plan document itself stand up to
inheritance / cross-reference?"

---

## Findings

### [R-101] CRITICAL `Phase 0/1/2 reuse of `spatial_quasi_dynamic_plan.md`` — Shared phases inherit critical defects

**Category:** DEVIATION (transitive)

**Description:**
This plan delegates Phases 0/1/2/3 (TOML schema, friction parser/resolver,
velocity bundle, stress bundle) to `spatial_quasi_dynamic_plan.md`. Several
CRITICAL findings against the QD plan apply unchanged here because the
delegated text is what the dynamic-rupture plan will consume:

| QD-review ID | Title | This-plan impact |
|---|---|---|
| QD R-002 | `geoffrey2010.md` `mu_s=1.1` fails `mu_s < 1` validator | LSW IS the production default for this driver → bites on first run. **Mitigated by user decision D-3 (relax validator + accept `d_o`)**; ensure the dynamic plan's `/code-plan` revision actually pulls the relaxed validator into its own §Phase 0. |
| QD R-003 | Velocity sidecar `velocity_safs.h5` (no mesh-tag) vs plan's `velocity_safs_<mesh_tag>.h5` | Identical for this plan; both plans share the same `LoadSpatialVelocityBundle` path resolver. |
| QD R-004 | `SidecarMaterialBundle` / `MakeMaterialField` facade not in code | Identical; this plan's Phase 3 says "no new code, reuse QD Phases 2/3" but the QD Phase 2 wraps a non-existent facade. |
| QD R-008 | `BoundaryConfig::fault_attr` default mismatch + `safs_test_driver` doesn't exist + single-sides Dirichlet sign | Identical; this plan reads SAFS tags 101/102/103/104 but never reconciles the BP5-style `BoundaryConfig` defaults. |

**Trigger:** Implementer reaches Phases 0/1/2/3 with only this plan in hand and
takes the cross-references at face value.

**Actual behavior:** All four QD defects surface in the dynamic-rupture
implementation too.

**Expected behavior:** The `/code-plan` revision must either (a) inline the
fixed shared text into this plan, or (b) keep the cross-reference but bump the
QD plan to the fixed version FIRST and pin the reference (e.g. "Phase 0 as
specified in `spatial_quasi_dynamic_plan.md` Phase 0 **rev-2** dated
2026-05-18 or later").

**Suggested fix:** Pick option (a) for this revision — inline the four shared
phases so the dynamic-rupture plan is self-contained against drift. The
duplicated text is small (≈ 250 LOC of plan markdown) compared to the
synchronisation risk.

**Test case (consistency CI):**
```python
def test_R101_shared_phase_consistency():
    qd = open("spatial_quasi_dynamic_plan.md").read()
    dyn = open("spatial_dynamic_rupture_plan.md").read()
    # Any Phase-0 TOML block that appears in dyn must appear verbatim in qd
    # (or vice-versa, post-inlining):
    for block in extract_toml_blocks(dyn):
        assert block in qd, f"Drift detected: {block[:80]}..."
```

---

### [R-102] CRITICAL `Overview §35, Constraints §47, Phase 4 §6,11` — `WaveOperator(MaterialField)` constructor and `GodunovFluxPool` do not exist

**Category:** ASSUMPTION (unverified prereq, blocks Phase 4 directly)

**Description:**
The plan repeatedly leans on two artifacts that come from
`heterogeneous_material_plan.md` Phase 3:
- Line 35: "Wiring `MaterialField::Mode::Coefficient` into `WaveOperator`
  (delegated entirely to `heterogeneous_material_plan.md` Phase 3 —
  `GodunovFluxPool`)."
- Line 47: "`WaveOperator<ParMesh>` `(MaterialField)` constructor from
  `heterogeneous_material_plan.md` Phase 3 is the only one used here."
- Phase 4 step 6 (line 460): "Construct `WaveOperator<ParMesh>` with the
  `MaterialField` constructor (heterogeneous_material_plan Phase 3)."

Actual code (`dynamic/wave_operator.hpp:88`):
```cpp
WaveOperator(MeshType &mesh, int order,
             real_t lambda, real_t mu, real_t rho,
             const BoundaryConfig &bc);
```
— scalar `(lambda, mu, rho)`, no `MaterialField` overload.

Actual `grep` for `GodunovFluxPool` in `miniapps/seas/`:
- Zero hits in `dynamic/`, `solver/`, `io/`, `domain/`, `fault/`, `drivers/`.
- The string only appears inside `safs/project_7.0_alternative/document/`
  (this plan, the QD sibling, and the unit-test plan all reference it) and
  two old `system_dev/dynamic_rupture_plan_v{1,3}.md` design docs.
- Only `class GodunovFlux` exists (`dynamic/godunov_flux.hpp:36`), and it is a
  single per-`WaveOperator` instance member (`wave_operator.hpp:611:
  GodunovFlux flux_;`), NOT a pool keyed by `(λ, μ, ρ)` triple.

So both prereqs are still vapor. Per user decision **D-2**, this plan must
declare a **hard Preconditions gate** before Phase 4.

**Trigger:** Implementer reaches Phase 4 step 6 OR Phase 5 step that asks
WaveOperator to expose `GodunovFluxPool`-keyed behavior.

**Actual behavior:** Compile/link error; if guarded by a runtime check, abort
at startup.

**Expected behavior:** Phase 4 is blocked at compile time until the prereq
exists. Phases 0/1/2/3 (TOML, friction, velocity bundle skeleton, stress
bundle) can proceed in parallel.

**Suggested fix (in the plan):** Insert a §Preconditions section between
§Constraints and §Phase 0:

```diff
+ ## Preconditions (must be true before Phase 4 implementation starts)
+
+ Each precondition has a build-time guard.  CI fails when a precondition
+ regresses.  No implementer should land Phase 4 changes that bypass these.
+
+ 1. `heterogeneous_material_plan.md` Phase 3 is **merged to main** — the
+    new ctor
+    `WaveOperator(MeshType&, int order, const MaterialField&,
+                  const BoundaryConfig&)`
+    exists, and the `GodunovFluxPool` exists.  Build-time guard,
+    inserted at the top of `drivers/spatial_dyn_driver.cpp`:
+    ```cpp
+    static_assert(std::is_constructible_v<WaveOperator<ParMesh>,
+                  ParMesh&, int, const MaterialField&, const BoundaryConfig&>,
+                  "Phase 4 of spatial_dynamic_rupture_plan requires "
+                  "heterogeneous_material_plan Phase 3 to be merged.");
+    ```
+ 2. `heterogeneous_material_plan.md` Phase 1 is merged — i.e.
+    `MaterialField`, `Mode::Coefficient`, `MakeCoefficient`, `EvalAt(elem,
+    T, ip, ...)` exist as documented.  Verified by:
+    `make test-heterogeneous-material` exits 0 with all 34/34 tests passing.
+ 3. `FaultGeometry<ParMesh>` BP5 constructor + `ComputeSAFSParams` are
+    available (already merged on this branch; verified by
+    `seas_test_compute_safs_params`).
+
+ Phases 0/1/2/3 do not depend on Preconditions 1-2; they MAY proceed in
+ parallel with the prereq merge.  Phase 4 MUST NOT be merged before all
+ three preconditions hold.
```

**Test case:**
```bash
# Pre-Phase-4 guard:
g++ -fsyntax-only drivers/spatial_dyn_driver.cpp 2>&1 | \
  grep -q "static_assert.*WaveOperator.*MaterialField" && exit 1 || exit 0
```

---

### [R-103] CRITICAL `Constraints §51, Phase 2 §Goal/Files-to-Modify` — `FaultFaceFlux::InitializeImpedancesFromMaterial` does NOT exist; replaceable scalar API is the ctor itself

**Category:** ASSUMPTION (API misnamed) + DEVIATION

**Description:**
Plan line 51: "`FaultFaceFlux` impedance computation — currently uses constant
`(λ, μ, ρ)` from the operator (`InitializeImpedancesFromMaterial`). For SAFS
with `Mode::Coefficient`, impedances must be **per-fault-QP** evaluated."

Plan Phase 2 §Files to Modify, line 315: "Replace the existing
`FaultFaceFlux::InitializeImpedancesFromMaterial(scalar_lambda, scalar_mu,
scalar_rho)` with a per-fault-QP variant…"

Actual code:
- `FaultFaceFlux::FaultFaceFlux(real_t rho, real_t cp, real_t cs)` is the
  constructor (`dynamic/fault_face_flux.hpp:124`). Impedances are set in the
  ctor body (or implicitly through `DOFData` — see below).
- `DOFData` has **per-DOF** fields `Zp_plus, Zp_minus, Zs_plus, Zs_minus, eta_p,
  eta_s` (`dynamic/fault_face_flux.hpp:27-29`). These are written **per DOF
  per driver** by `InitializeFaultDOFs_TPV205` (e.g. `dynamic/tpv205_setup.hpp:96-101`):
  ```cpp
  d.Zp_plus  = TPV205Params::Zp;
  d.Zp_minus = TPV205Params::Zp;
  d.Zs_plus  = TPV205Params::Zs;
  d.Zs_minus = TPV205Params::Zs;
  d.eta_p    = TPV205Params::eta_p;
  d.eta_s    = TPV205Params::eta_s;
  ```
- There is NO `InitializeImpedancesFromMaterial` method anywhere
  (`grep -rn "InitializeImpedances" miniapps/seas/` returns zero matches).

The plan's diagnosis ("currently uses constant (λ, μ, ρ) from the operator
**via `InitializeImpedancesFromMaterial`**") is fictional, and the proposed
fix (replace that method with `InitializeImpedancesPerQP`) is therefore
chasing a nonexistent baseline.

What the plan **should** propose:
1. Either teach `InitializeFaultDOFs_*` to take a `MaterialField` and call
   `material.EvalAt(elem, T, ip, ...)` per DOF (this is the SAFS-equivalent
   of the existing per-driver setup function), OR
2. Add a new `FaultFaceFlux::InitializeImpedancesPerDOF(...)` setter that
   overwrites the impedance fields in `dof_data_` AFTER
   `InitializeFaultDOFs_*` runs, AND audit the call order so the new setter
   runs second.

Either way the change is **per-DOF**, not "per-QP" — the existing storage is
already per-DOF (one `DOFData` per fault QP). The plan's "per-QP" naming
suggests an additional storage layer that does not exist.

**Trigger:** Phase 2 implementation tries to delete a method that doesn't
exist and add a "per-QP" overlay on top of per-DOF storage.

**Actual behavior:** Compile error (deleted method not found) or silent
double-init race if the new setter writes the same `DOFData` fields the
TPV-specific `InitializeFaultDOFs_*` already wrote.

**Expected behavior:** Plan acknowledges the existing per-DOF storage and
defines exactly one of (a) extend `InitializeFaultDOFs_*` (SAFS-specific
free function), (b) add `FaultFaceFlux::InitializeImpedancesPerDOF` (member
setter, asserts it runs after `InitializeFaultDOFs_*`).

**Suggested fix (in the plan):**
```diff
- ## Phase 2: Heterogeneous fault-QP impedance (delivered, not deferred)
+ ## Phase 2: Heterogeneous fault-DOF impedance (delivered, not deferred)

  …

- - Replace the existing `FaultFaceFlux::InitializeImpedancesFromMaterial`
-   (scalar_lambda, scalar_mu, scalar_rho) with a per-fault-QP variant…
+ - There is NO existing `InitializeImpedancesFromMaterial` to replace.
+   Impedances today are written into per-DOF `DOFData[i].Z{p,s}_{plus,minus}`
+   and `eta_{p,s}` by the per-driver `InitializeFaultDOFs_TPV*` free
+   functions (e.g. `dynamic/tpv205_setup.hpp:96`).  The SAFS path adds an
+   analogous per-driver helper `InitializeFaultDOFs_Spatial` that takes a
+   `MaterialField` and calls `material.EvalAt(elem, T, ip, ...)` to fill
+   per-DOF impedances.  No member function is added to `FaultFaceFlux`.
+
+ ```cpp
+ // dynamic/spatial_setup.hpp (NEW)
+ namespace mfem::seas::spatial {
+ void InitializeFaultDOFs_Spatial(
+    std::vector<DOFData>&            dof_data,
+    int                              ndof,
+    const std::vector<mfem::Vector>& fault_coords,
+    const Array<int>&                dof_to_elem,
+    const MaterialField&             material,
+    ParMesh&                         pmesh,
+    const SlipWeakeningPerDOFParams& lsw,           // R-104
+    const Vector&                    tau_pre_per_dof, // [2*ndof]
+    const Vector&                    sigma_n_eff);    // [ndof]
+ }
+ ```
+
+ The function:
+   1. For each DOF i:
+      - look up bulk element e = dof_to_elem[i];
+      - get ElementTransformation T from pmesh.GetElementTransformation(e);
+      - get reference centroid IntegrationPoint ip;
+      - call material.EvalAt(e, T, ip, lam, mu, rho);
+      - compute cp = sqrt((lam + 2 mu) / rho), cs = sqrt(mu / rho);
+      - set d.Zp_plus = d.Zp_minus = rho * cp;  d.Zs_plus = d.Zs_minus = rho * cs;
+      - set d.eta_p = 0.5 * rho * cp;  d.eta_s = 0.5 * rho * cs;
+   2. Write per-DOF LSW fields (lsw_mu_s, lsw_mu_d, lsw_d_c) from `lsw`.
+   3. Write per-DOF tau1_0, tau2_0, sigma_n0 from `tau_pre_per_dof`, `sigma_n_eff`.
+   4. Zero rate-state fields (a, psi, Dc) per the TPV205 defensive pattern.
+
+ This subsumes Phase 1's `InitializeFromSpatialStress` (R-105 makes it
+ redundant; merge the two setters).
```

**Test case:**
```cpp
TEST(R103_PerDOFImpedance, MatchesScalarOnConstantMaterial) {
    // Build a MaterialField with MakeConstant(lam, mu, rho).
    // Call InitializeFaultDOFs_Spatial.
    // Assert dof_data[i].eta_s == 0.5 * rho * sqrt(mu/rho) for every i.
    // Assert dof_data[i].eta_p == 0.5 * rho * sqrt((lam+2*mu)/rho) for every i.
}
TEST(R103_PerDOFImpedance, VariesAcrossLayer) {
    // MaterialField with two-layer Vp(z): layer A at z>-5km, layer B at z<-5km.
    // Assert dof_data at z=-2km has eta_p of layer A;
    //        dof_data at z=-10km has eta_p of layer B (differs by > 1e-3 rel).
}
```

---

### [R-104] CRITICAL `Phase 1 §Files to Modify` (lines 224-247) — Two-stage init (`InitializeFromSpatialStress` then `InitializeFromSpatialStressRateState`) on the same `dof_data_` racing the TPV-specific `InitializeFaultDOFs_*`

**Category:** EDGE_CASE (silent double-init / partial overwrite)

**Description:**
Phase 1 proposes two new public setters on `FaultFaceFlux`:
- `InitializeFromSpatialStress(tau_pre, sigma_n_eff, lsw_mu_s, lsw_mu_d,
  lsw_d_c, lsw_cohesion)` — LSW branch
- `InitializeFromSpatialStressRateState(...)` — RS branch

Each "Aborts if called twice" (line 254-255). But:
1. The existing TPV205 path calls **the free function**
   `InitializeFaultDOFs_TPV205(dof_data, num_fault_total, fault_coords)` to
   populate `dof_data_`. If the SAFS driver also calls a TPV-style setup
   function (it must, to fill the impedances per R-103) AND
   `InitializeFromSpatialStress`, both write to overlapping fields
   (`tau1_0`, `tau2_0`, `sigma_n0`, `lsw_mu_s`, etc.). Last-writer-wins is
   not guaranteed.
2. The "Aborts if called twice" guard only catches second calls **to itself**,
   not the cross-route conflict with the free-function setup.

This makes the call order critical and unspecified.

**Trigger:** Driver does `InitializeFaultDOFs_TPV205(...)` (or any analog) then
calls `InitializeFromSpatialStress(...)`; OR vice-versa; OR forgets one.

**Actual behavior:** Silent partial init; some fields hold TPV205 defaults
(e.g. `lsw_mu_s = ComputeMuS_TPV205(...)` patch-based) while others hold SAFS
sidecar values. The simulation runs but produces nonsense.

**Expected behavior:** Single init pipeline that writes every field exactly
once.

**Suggested fix:** Per R-103, MERGE the two `Initialize*FromSpatialStress`
setters into a single SAFS-specific free function `InitializeFaultDOFs_Spatial`
that writes every `DOFData` field in one pass. Drop the
`FaultFaceFlux::InitializeFromSpatialStress` and
`FaultFaceFlux::InitializeFromSpatialStressRateState` member setters
entirely. The "aborts if called twice" guard becomes "aborts if any
`DOFData[i].lsw_mu_s != NaN-sentinel` before the call".

```diff
- void InitializeFromSpatialStress(const Vector& tau_pre_per_dof,
-                               const Vector& sigma_n_eff_per_dof,
-                               const Vector& lsw_mu_s_per_dof,
-                               const Vector& lsw_mu_d_per_dof,
-                               const Vector& lsw_d_c_per_dof,
-                               const Vector& lsw_cohesion_per_dof);
- void InitializeFromSpatialStressRateState(...);
+ /// REMOVED.  Replaced by free function `InitializeFaultDOFs_Spatial(...)`
+ /// in `dynamic/spatial_setup.hpp` (R-103, R-104).  That function is the
+ /// SAFS analog of `InitializeFaultDOFs_TPV205` / `InitializeFaultDOFs_TPV104`.
```

**Test case:**
```cpp
TEST(R104_NoDoubleInit, AbortsIfMixedWithTPV205) {
    std::vector<DOFData> dof_data(N);
    InitializeFaultDOFs_TPV205(dof_data, N, coords);
    EXPECT_DEATH(InitializeFaultDOFs_Spatial(dof_data, N, ...),
                 "dof_data already initialised");
}
```

---

### [R-105] CRITICAL `Phase 4 §Detailed Req. 14, Edge Cases` (lines 468, 480) — Single-file `tpv104_checkpoint.h5` and `Tpv104Checkpoint::Restore` API do not exist; checkpoints are per-rank `.txt`

**Category:** ASSUMPTION (API misnamed) — identical to QD R-006 but rephrased for V1 schema

**Description:**
Plan line 468: "if `--restart PREFIX`, restore via
`Tpv104Checkpoint::Restore("<PREFIX>/tpv104_checkpoint.h5", state, t0, dt0,
schedule_state)`."

Plan line 466: "the checkpoint file is named `tpv104_checkpoint.h5` (V1
schema — preserved name for tool compatibility)."

Actual API (`io/tpv104_checkpoint.hpp`):
- Free functions `WriteTpv104CheckpointImpl(prefix, ..., mpi)` and
  `ReadTpv104CheckpointImpl(prefix, ..., mpi)`.
- Path layout (`io/tpv104_checkpoint.hpp:60-63`):
  ```cpp
  inline std::string Tpv104CheckpointFilename(const std::string &prefix,
                                              int rank)
  {
     return CheckpointFilename(prefix, rank);   // {prefix}_checkpoint_r{rank}.txt
  }
  ```
- `io/checkpoint.hpp:31-36`:
  ```cpp
  // Format: {prefix}_checkpoint_r{rank}.txt
  inline std::string CheckpointFilename(const std::string &prefix, int rank)
  ```

There is **no class `Tpv104Checkpoint`** and **no `.h5` file** in the V1
schema; files are per-rank `.txt`. The plan's `--restart PREFIX` semantics
(treating PREFIX as a directory containing `tpv104_checkpoint.h5`) is wrong.

The user-facing rule established by the existing `tpv104_driver.cpp` (line 501)
is `--restart PREFIX` where PREFIX is the prefix string passed to
`Tpv104CheckpointFilename(prefix, rank)`, NOT a directory.

**Trigger:** Implementer reads §Output Directory Layout (line 66) or §Phase 4
step 14.

**Actual behavior:** Driver looks for a nonexistent `.h5` file; restart fails.

**Expected behavior:** Plan reflects per-rank `.txt` path layout AND the
existing free-function API.

**Suggested fix:**
```diff
- the checkpoint file is named `tpv104_checkpoint.h5`
+ Checkpoint files are per-rank `.txt` produced by
+ `Tpv104CheckpointFilename(prefix, rank) = "{prefix}_checkpoint_r{rank}.txt"`.
+ No HDF5 file is involved in the V1 schema; the suffix is `.txt` for
+ historical reasons.

...

- 14. **Restart path** — if `--restart PREFIX`, restore via
-     `Tpv104Checkpoint::Restore("<PREFIX>/tpv104_checkpoint.h5", state, t0,
-     dt0, schedule_state)`.
+ 14. **Restart path** — if `--restart PREFIX`, restore via the existing
+     free functions (each rank reads its own file):
+     ```cpp
+     ReadTpv104Checkpoint(restart_prefix, t0, dt0, step, state,
+                          paraview_snapshots, paraview_last_write_time,
+                          paraview_last_v_max, paraview_current_regime,
+                          paraview_last_committed_cycle,
+                          paraview_last_volume_write_time, mpi_ctx);
+     ```
+     PREFIX is a string (not a directory).  The driver's
+     `weakly_canonical` safety check on `--output-dir` vs PREFIX's parent
+     mirrors `tpv104_driver.cpp:521-545`.
```

**Test case:**
```python
def test_R105_checkpoint_layout():
    # After 4-rank write at prefix="run/cp":
    for r in range(4):
        assert os.path.exists(f"run/cp_checkpoint_r{r}.txt")
    assert not os.path.exists("run/tpv104_checkpoint.h5")
```

---

### [R-106] CRITICAL `Constraints §50, Phase 4 §Edge Case driver_tag` (lines 50, 480) — `driver_tag` schema extension to V1 checkpoint is risky and not specified

**Category:** EDGE_CASE / DEVIATION

**Description:**
Plan line 50: "`Tpv104Checkpoint` (V1 schema) is the restart mechanism.
Reused verbatim with a **new `driver_tag = "spatial_dyn"` field in the
schema's metadata** so a checkpoint cannot be restored into the wrong driver."

Plan line 480: "`--restart PREFIX` schema mismatch — checkpoint contains
`driver_tag` metadata; if `driver_tag != "spatial_dyn"`, abort. This
requires adding a 32-byte `driver_tag` string to the V1 checkpoint schema
header — additive, length-prefixed, defaults to "tpv104" for back-compat
with existing TPV checkpoints."

Issues:
1. This is a **schema mutation to a critical I/O contract** used by TPV104
   AND TPV205 AND TPV102 AND the petsc_ts_restart_plan_2026-05-16 V2 BP5
   path. Any bug in the back-compat default ("tpv104" for missing field)
   silently corrupts TPV restart for all existing checkpoint files in flight.
2. The plan calls this "additive, length-prefixed" but the V1 schema is a
   text file (`*.txt`). Length-prefixing a string in a text file is
   unconventional and contradicts the file format.
3. The existing PETSc-TS V2 extension (`WritePetscTSCheckpoint`) appends a
   `PETSC_TS_V2` magic followed by named fields (one per line). The clean
   way to add `driver_tag` is to add it as a new named field in that
   extension, not to "extend the V1 header".
4. The plan's "32-byte `driver_tag` string" implies a binary fixed-width
   format, but everything around it is text.

**Trigger:** Implementer modifies the V1 schema; old TPV104 checkpoint files
fail to round-trip.

**Actual behavior:** Either back-compat default works and tag is useless, OR
back-compat breaks all existing checkpoints.

**Expected behavior:** Tag is implemented in a way that strictly preserves V1
back-compat:
- Tag lives in a **new optional V1-extension block** (e.g. `DRIVER_TAG_V1\n<tag>\n`)
  appended at the end of the per-rank `.txt` file.
- Readers that don't find the magic default to `driver_tag = ""` (no
  enforcement) — NOT to `"tpv104"`, because mis-tagging is silent-corruption-
  risk.
- Driver-side check is `if (file_tag != "" && file_tag != my_tag) abort();`
  — i.e. allow restoring untagged checkpoints, but refuse cross-driver
  tagged restores.

**Suggested fix (in the plan):**
```diff
- This requires adding a 32-byte `driver_tag` string to the V1 checkpoint
- schema header — additive, length-prefixed, defaults to "tpv104" for back-
- compat with existing TPV checkpoints.
+ Implementation:
+   1. The V1 per-rank `.txt` schema is preserved bit-identical.  An OPTIONAL
+      extension block `DRIVER_TAG_V1\n<tag string>\n` is appended at the end
+      of the file by `WriteTpv104Checkpoint` after the existing trailing
+      fields.  No format change to the V1 body.
+   2. `ReadTpv104Checkpoint` returns `driver_tag = ""` (empty) when the
+      block is absent — i.e. it does NOT fabricate a default tag like
+      "tpv104".  Mis-tagging is silent-corruption-risk.
+   3. Driver-side check: aborting only when both file_tag and runtime_tag are
+      non-empty AND differ.  Untagged old TPV104 files restore into the
+      tpv104 runtime without issue; spatial_dyn explicitly tags new files
+      with "spatial_dyn"; spatial_dyn refuses to restore a file tagged
+      "tpv104" or "tpv205".
+   4. Unit-test gate added before merge: T-CHECKPOINT-DRIVER-TAG-BACK-COMPAT.
+      Reads a stored pre-extension TPV104 checkpoint and verifies the
+      tpv104_driver round-trips without aborting.
```

**Test case:**
```python
def test_R106_driver_tag_backcompat():
    # Pre-extension TPV104 checkpoint (no DRIVER_TAG_V1 block):
    write_tpv104_checkpoint_pre_extension("legacy")
    tag = read_driver_tag("legacy")
    assert tag == ""  # not "tpv104"
    # tpv104_driver allows restore from "":
    restore_with_runtime_tag(legacy, runtime_tag="tpv104")  # no abort
    # spatial_dyn driver also allows restore from "" (back-compat):
    restore_with_runtime_tag(legacy, runtime_tag="spatial_dyn")  # no abort

def test_R106_driver_tag_cross_driver_refused():
    write_tpv104_checkpoint_with_tag("new", tag="tpv104")
    with pytest.raises(AssertionError, match="tag mismatch"):
        restore_with_runtime_tag("new", runtime_tag="spatial_dyn")
```

---

### [R-107] CRITICAL `Phase 4 §Detailed Req. 4` (line 458) — SAFS `BoundaryConfig` not specified; existing `BoundaryConfig::fault_attr = 3` defaults break the dynamic-rupture init

**Category:** ASSUMPTION (carries over from QD R-008; impact is harder here)

**Description:**
Plan line 458: "Load mesh → `ParMesh pmesh`. Apply the `BoundaryConfig` for
the SAFS .geo (fault_attr=101, top=102, bottom=103, sides=104)."

Two CRITICAL issues, both already noted for the QD plan (QD R-008) but with
LARGER impact in the dynamic-rupture context:

1. `BoundaryConfig::fault_attr` default (`domain/boundary_config.hpp:47`) is
   `3`. The plan instructs to "apply" — but does not say how. If the driver
   default-constructs and forgets to set `fault_attr = 101`, the WaveOperator
   will identify **zero fault faces**. In the dynamic-rupture context this
   means **no Riemann fault solver fires at all**; the simulation runs the
   wave equation in pure homogeneous-bulk mode and "succeeds" silently.
2. The SAFS mesh uses a SINGLE attribute (104) for all four vertical sides.
   For dynamic-rupture, the BCs are typically **absorbing** (or PML — see
   `dynamic/pml_layer.hpp`) on the side and bottom, **free-surface** on top.
   The plan does not say whether sides=104 should be Absorbing, PML, or
   FreeSurface. TPV104/205 use absorbing BCs on sides (with the optional
   `--pml` flag adding PML).

The plan is silent on PML for SAFS. Without absorbing or PML on the lateral
sides, reflected waves contaminate the rupture for `t > L_box / cp ≈
50 km / 6 km·s⁻¹ = 8 s` — close to the typical simulation window.

**Trigger:** First Phase 4 dry-run.

**Actual behavior (case 1):** Zero fault faces; WaveOperator returns "no
spontaneous rupture detected" silently. (case 2): Reflected waves from
sides break event nucleation for `tfinal > 8s`.

**Expected behavior:** `BoundaryConfig` is constructed explicitly with the
SAFS tags, with absorbing/PML on attrs 103+104.

**Suggested fix:**
```diff
- 4. **Load mesh** → `ParMesh pmesh`.  Apply the `BoundaryConfig` for the
-    SAFS .geo (fault_attr=101, top=102, bottom=103, sides=104).
+ 4. **Load mesh** → `ParMesh pmesh`.  Construct `BoundaryConfig` explicitly:
+    ```cpp
+    BoundaryConfig bc;
+    bc.fault_attr      = 101;
+    bc.natural_attrs   = {102};         // top: free surface
+    bc.absorbing_attrs = {103, 104};    // bottom + four sides: absorbing
+    // (PML available as an opt-in via --pml; constructed separately.)
+    ```
+    If the SAFS box dimensions (≈ 358 km E-W x 247 km N-S x 42 km depth) are
+    smaller than `cp * tfinal`, the user MUST add `--pml` to avoid reflection
+    contamination.  Driver prints a one-line warning at startup if
+    `min_box_dim / cp_max < tfinal` AND `--pml` is not set.
```
Cross-check with `dynamic/pml_layer.hpp` for the actual `BoundaryConfig`
field name (`absorbing_attrs` is what I expect from convention; verify
during implementation).

**Test case:**
```cpp
TEST(R107_BoundaryConfig, FaultAttrIsSet) {
    auto bc = ConstructSAFSBoundaryConfig(pmesh);
    EXPECT_EQ(bc.fault_attr, 101);
    EXPECT_TRUE(bc.absorbing_attrs.count(103));
    EXPECT_TRUE(bc.absorbing_attrs.count(104));
    EXPECT_EQ(bc.natural_attrs.count(102), 1);
}
TEST(R107_WaveOp, NonzeroFaultFacesOnSafsMesh) {
    auto wave = ConstructWaveOpFromSafsMesh();
    EXPECT_GT(wave.NumFaultFaces(), 0) << "WaveOp got zero fault faces; "
        "BoundaryConfig.fault_attr is likely wrong";
}
```

---

### [R-108] CRITICAL `User Objective vs Phase 3` — User's constant-tensor + rotation testing path (user decision D-1) not represented in the plan

**Category:** DEVIATION (user objective)

**Description:**
Per user decision **D-1**: constant-background-tensor stress with rotation
onto the per-DOF fault normal is **first-class** for testing. The plan as
written (`Phase 3 §Goal`, line 378: "No new code. Re-document that this plan
reuses Phases 2 (LoadSpatialVelocityBundle) and 3 (ApplySpatialStressSidecar)
from the QD plan.") only wires the CSM HDF5 sidecar.

The constant-tensor path must be addable without rebuilding `FaultGeometry::
ComputeSAFSParams` — that function is templated on a "field that exposes
`Evaluate(x, y, z) -> DenseMatrix`". Either (a) make
`ComputeSAFSParams` truly template on a `Source` concept, OR (b) build a
constant-tensor `StressField3D`-compatible adapter (a tiny class with the
same `Evaluate` signature). The QD plan's R-009 already proposes (b).

**Trigger:** User wants to validate the per-DOF projection by running
SAFS with a known constant tensor (analytic resolved-traction comparison)
BEFORE trusting the CSM sidecar.

**Actual behavior:** Plan offers no such test path.

**Expected behavior:** Phase 3 of the dynamic-rupture plan ADDS a constant-
tensor source and gates Phase 4 on its acceptance:
```toml
[stress]
kind = "constant_tensor"   # "constant_tensor" | "sidecar_hdf5"
# Cauchy tensor in Pa, EAST-NORTH-UP frame, compression POSITIVE:
sigma_xx_pa = ...
sigma_yy_pa = ...
sigma_zz_pa = ...
sigma_xy_pa = ...
sigma_yz_pa = ...
sigma_xz_pa = ...
# Pore-pressure block as before.
```
With the per-DOF projection going through the same `ComputeSAFSParams`
machinery so that the projection correctness test runs on both paths.

**Suggested fix:** Insert a new sub-phase Phase 3b in this plan (do NOT
delegate it to the QD plan — the QD plan needs it too, but the dynamic-
rupture plan owns the first-class user-facing entry point per D-1):

```diff
+ ### Phase 3b: Constant-tensor stress source (testing primitive — owned by this plan)
+
+ Per user decision D-1, the dynamic-rupture driver is the FIRST consumer of
+ the constant-tensor stress source.  See the QD plan's R-009 for the
+ algorithm; this phase delivers it.
+
+ Files to Create:
+   - `spatial/code/spatial_constant_stress_source.hpp/.cpp` — a thin class
+     `ConstantTensorStressSource` whose `Evaluate(x, y, z)` returns the
+     same constant 3x3 DenseMatrix at every (x, y, z).  Implements the
+     `StressSource3D` concept (see below).
+
+ Files to Modify:
+   - `fault/fault_geometry.hpp`: introduce a `StressSource3D` concept (or
+     a small `IStressSource3D` virtual base) such that
+     `FaultGeometry::ComputeSAFSParams` accepts EITHER a `StressField3D`
+     OR a `ConstantTensorStressSource`.  Easiest path: template
+     `ComputeSAFSParams` on the source type and rely on duck-typing of
+     `Evaluate(x, y, z) -> DenseMatrix`.  Pin the existing call sites
+     against bit-exact regression by keeping a non-templated overload
+     that forwards to the templated body.
+
+ Acceptance:
+   - [ ] On a planar fault DOF subset (e.g. all DOFs with |y - y0| < 5m),
+     per-DOF normal traction matches the analytic `n · σ · n` to L∞ ≤ 1e-10
+     relative.
+   - [ ] On the curvilinear SAFS fault, per-DOF projection agrees with a
+     Python reference (mirrors `stress/code/hickman_and_zoback_regional_
+     stress_projection.py::resolve_traction`) at 64 random DOFs (L∞ ≤ 1e-8
+     relative).
+   - [ ] `seas_spatial_dyn_driver --dry-run` with
+     `[stress] kind = "constant_tensor"` exits 0 and prints a non-degenerate
+     histogram of `geom.sigma_n_per_dof()`.
+
+ Dependencies: Phase 0 (the [stress] TOML block).
+ Required by:  Phase 4.
```

**Test case (cross-language):**
```python
def test_R108_constant_tensor_projection_matches_python_ref():
    sigma = np.array([[100e6, 10e6, 0], [10e6, 60e6, 0], [0, 0, 80e6]])
    # Pick 64 random DOFs from the SAFS fault.
    for (n, dof_idx) in sample_random_dofs(64):
        py_sn = float(n @ sigma @ n)
        cpp_sn = read_dry_run_sigma_n(dof_idx)
        assert abs(py_sn - cpp_sn) / abs(py_sn) < 1e-8
```

---

### [R-109] CRITICAL `User Objective vs Plan §Acceptance` — Two-CVM comparison workflow (`cvm_s4.26.m01` vs `multiscale_statewise_cvm`) not delivered

**Category:** DEVIATION (user objective, identical to QD R-010)

**Description:**
The user explicitly asked for COMPARISON between two CVMs. The plan's
Phase 4 §Acceptance Criteria asks only for `--dry-run` exit-0 on cvmh; no
acceptance criterion runs both CVMs and diffs their stiffness / wave-speed /
slip-rate / event-onset-time.

**Trigger:** User wants the deliverable.

**Actual behavior:** User scripts the comparison by hand.

**Expected behavior:** First-class comparison phase.

**Suggested fix:** Add to Phase 6 (sbatch + verification):
```diff
+ - `jobs/safs/spatial_dyn_smoke_cvm_compare_dev_2hr_safs.sbatch` —
+   sequentially runs BOTH `cvm_s4.26.m01` and `multiscale_statewise_cvm`
+   into `out_a/` and `out_b/`, then runs
+   `spatial/code/scripts/compare_spatial_dyn_velocity_models.py` to:
+   - load `out_a/fault.vtkhdf` and `out_b/fault.vtkhdf`,
+   - diff slip-rate over space-time, emit L∞, L2, and 95th-percentile metrics,
+   - emit a one-page PNG summary (rupture-front contour overlay).
+   - Exit code 0 if the two runs differ by ≥ 1% L2 (proves the velocity
+     model actually entered the result).
```

**Test case:**
```python
def test_R109_two_cvm_runs_produce_distinct_results():
    run_a = LoadFault("out_a/fault.vtkhdf")
    run_b = LoadFault("out_b/fault.vtkhdf")
    rel_l2 = compute_relative_l2(run_a.slip_rate, run_b.slip_rate)
    assert rel_l2 > 0.01, ("velocity model has no effect on result — "
                           "check that MaterialField is actually wired "
                           "into WaveOperator")
```

---

### [R-110] CRITICAL `Phase 4 §Detailed Req. 7` (line 461) — "throw-away ElasticityDomainOperator for fault-DOF bookkeeping" duplicates the WaveOperator's bulk-side fault traversal; not actually wired together

**Category:** DEVIATION (architectural fragility)

**Description:**
Plan line 461: "dynamic-rupture still uses `ElasticityDomainOperator` for
fault-DOF bookkeeping (matching TPV104's pattern), even though the bulk
physics goes through `WaveOperator`. Construct a throw-away
ElasticityDomainOperator with the same MaterialField purely to provide the
fault-DOF coordinate / basis arrays."

This is **inaccurate** about the existing pattern:
- `tpv104_driver.cpp` does NOT construct an `ElasticityDomainOperator`. It
  builds `WaveOperator wave(...)` directly, then constructs
  `FaultFaceFlux fault_flux(TPV104Params::rho, TPV104Params::cp, TPV104Params::cs)`
  (line 1652) and calls `wave.SetFaultFlux(&fault_flux)` (line 1654). Fault
  coordinates come from `wave`'s own fault-face traversal, NOT from an
  elasticity operator.
- `WaveOperator` has its own `int GetNumFaultDOFs() const` (line 110) and
  exposes the bulk fault face set via `GetFaultFaceSet()` (line 199).

So the "throw-away ElasticityDomainOperator" is unjustified, and constructing
one carries non-trivial overhead: it assembles the stiffness matrix
(`InitOperator()` on line 134/191), which for a SAFS 1000m mesh is **several
GB of factored LU/BLR**. For a dynamic-rupture run that does NOT solve a
static problem, this is pure waste.

The actual reason the plan reaches for `ElasticityDomainOperator` is that
`FaultGeometry`'s BP5 constructor takes a `DomainOperator<MeshType>&`
(`fault_geometry.hpp:118-119`) and `WaveOperator` is NOT a `DomainOperator`.
`FaultGeometry::ComputeSAFSParams` and `geom.fault_dof_coords_3d()` depend on
`FaultGeometry`, which depends on `DomainOperator`'s
`GetFaultDOFCoords3D` / `GetFaultDOFBasis` / `GetNumOwnedFaultDOFs`.

There are two real fixes:
1. **Wide:** add an `IFaultDOFProvider` mixin or have `WaveOperator` derive
   from `DomainOperator`. Large blast radius.
2. **Narrow (recommended):** extract the per-DOF coordinate/basis computation
   from `ElasticityDomainOperator::InitFaultMaps` into a free function
   `ComputePerDOFCoordsAndBasisFromWave(WaveOperator&, ...)` that
   `FaultGeometry` can consume directly. No throw-away operator.

The plan must commit to one approach BEFORE Phase 4, because the choice
affects whether `MUMPS_BLR` is even linked in for this driver.

**Trigger:** Phase 4 implementation builds the SAFS dyn driver.

**Actual behavior (as written):** A multi-GB factored stiffness matrix is
built and thrown away every run; SAFS dyn driver is unnecessarily slow and
memory-hungry.

**Expected behavior:** Either no `ElasticityDomainOperator` (option 2), or a
documented "we DO factor the stiffness — yes, it costs N GB" choice.

**Suggested fix (in the plan):**
```diff
- 7. **Construct `FaultGeometry<ParMesh>`** (BP5 ctor with the SAFS-aware
-    `ElasticityDomainOperator`).  Note: dynamic-rupture still uses
-    `ElasticityDomainOperator` for fault-DOF bookkeeping (matching TPV104's
-    pattern), even though the bulk physics goes through `WaveOperator`.
-    Construct a throw-away ElasticityDomainOperator with the same
-    MaterialField purely to provide the fault-DOF coordinate / basis arrays;
-    this matches TPV104's existing two-operator pattern.
+ 7. **Construct fault-DOF bookkeeping** (does NOT match any existing
+    pattern; TPV104 does not construct an ElasticityDomainOperator).  Two
+    options, decide before merging:
+
+    Option A (recommended).  Add a thin free function
+    `BuildFaultDOFTablesFromWave(WaveOperator&, Array<int>& dof_to_elem,
+                                 Vector& dof_coords_3d, DenseMatrix& dof_basis)`
+    in `dynamic/spatial_setup.hpp`.  Mirror the logic of
+    `ElasticityDomainOperator::InitFaultMaps` (extract the fault-face
+    walk into a helper that both operators can call).  Then construct
+    `FaultGeometry` via a **new** BP5 ctor overload that takes the three
+    arrays directly:
+    ```cpp
+    FaultGeometry<ParMesh> geom(BP5Params{}, dof_coords_3d, dof_basis,
+                                dof_to_elem, mpi_ctx);
+    ```
+    Cost: ~50 LOC in `dynamic/spatial_setup.hpp` + a new FaultGeometry
+    ctor overload (additive).  Pinned by a bit-exact test
+    `T-FAULTGEO-NEW-CTOR-MATCHES-OLD-PATH` against the
+    `ElasticityDomainOperator`-driven construction on a small fixture.
+
+    Option B (fallback if Option A is too invasive).  Construct a
+    real `ElasticityDomainOperator` with the SAFS MaterialField and
+    accept the stiffness factorisation cost (several GB on SAFS 1000m).
+    Document the memory budget in the sbatch.  Use only the fault-DOF
+    accessors (`GetFaultDOFCoords3D`, `GetFaultDOFBasis`); never call
+    Solve / ComputeTraction.
+
+    Plan freezes on Option A. If A is unimplementable within the prereq
+    window, switch to B with explicit user sign-off.
```

**Test case:**
```python
def test_R110_no_unnecessary_factorisation():
    # mprof or /usr/bin/time -v on a dry-run:
    rss = run_dry_run_peak_rss(driver="spatial_dyn", mesh="1000m_lcfar3000")
    # Without ElasticityDomainOperator, peak RSS should be < N GB
    # (mesh + WaveOp arrays only).  With it, peak RSS is dominated by MUMPS.
    assert rss < 16e9, "Phase 4 step 7 chose Option B; verify this is intended"
```

---

### [R-111] MODERATE `Phase 2 §Detailed Req. body` (line 340) — `T->TransformBack(phys, ip)` API is not the MFEM standard

**Category:** BUG (API misnamed)

**Description:**
Plan body:
```cpp
mfem::Vector phys(3);
phys(0) = dof_coords_3d(3*i + 0);
phys(1) = dof_coords_3d(3*i + 1);
phys(2) = dof_coords_3d(3*i + 2);
T->TransformBack(phys, ip);          // physical → reference
```

`mfem::ElementTransformation::TransformBack(const Vector&, IntegrationPoint&)`
exists but returns an `int` status code and is documented as approximate. The
standard, reliable way to evaluate a `Coefficient` at a known DOF is to use the
DOF's reference coordinate directly (which `FaultGeometry` already knows from
the fault face's basis — stored in `dof_basis_`). Building a physical point
just to back-transform it is brittle and slow.

Better: capture the reference `IntegrationPoint` per DOF at the time the
fault-face walk happens (the same walk that builds `dof_coords_3d_`), store it
in `FaultGeometry`, and use it directly. This is also the approach
`MaterialField::EvalAt(elem, T, ip, ...)` expects — `ip` should be the
reference point, not a back-transformed estimate.

**Trigger:** Phase 2 implementation; MaterialField evaluation is off by the
back-transform tolerance.

**Actual behavior:** Per-DOF impedances have a small error (~ 1e-8 relative,
not always bit-exact); for higher-order elements with curvilinear meshes
(SAFS fault is curvilinear), error can be larger.

**Expected behavior:** Use the cached reference IP, not back-transform.

**Suggested fix:**
```diff
- mfem::Vector phys(3);
- phys(0) = dof_coords_3d(3*i + 0);
- phys(1) = dof_coords_3d(3*i + 1);
- phys(2) = dof_coords_3d(3*i + 2);
- T->TransformBack(phys, ip);          // physical → reference
+ // The reference IntegrationPoint per fault DOF is cached by FaultGeometry
+ // at construction time (`fault_dof_ip()`).  Use it directly — back-
+ // transform is approximate and unnecessary.
+ const IntegrationPoint& ip = geom.fault_dof_ip(i);
  real_t lam, mu, rho;
- material.EvalAt(e, *T, ip, lam, mu, rho);
+ material.EvalAt(dof_to_elem[i], *T, ip, lam, mu, rho);
```
And add `fault_dof_ip_` storage + `fault_dof_ip(i)` accessor to
`FaultGeometry` (additive; populated alongside `dof_coords_3d_` in
`ComputePerDOFCoordsAndBasis_`).

**Test case:**
```cpp
TEST(R111_PerDOFImpedance, ExactOnConstantMaterial) {
    auto mat = MaterialField::MakeConstant(lam, mu, rho);
    auto impedances = InitializeFaultDOFs_Spatial(...);
    for (int i = 0; i < N; ++i) {
        EXPECT_EQ(impedances[i].eta_p, 0.5 * rho * sqrt((lam + 2*mu) / rho));
        // (EXACT equality; back-transform would fail this.)
    }
}
```

---

### [R-112] MODERATE `Phase 0 §[time] §[output]` (sibling-shared) — Dynamic-rupture-specific defaults are inadequate; CFL handling is silent

**Category:** ASSUMPTION / DEVIATION

**Description:**
Sibling plan's `[time]` block (used by this plan):
```toml
[time]
tfinal      = "12s"
t_initial   = 0.0
dt_initial  = "auto"        # "auto" => 0.5 * CFL_local_min
dt_max      = "0.1s"        # safety cap
cfl         = 0.5
```

For the SAFS 1000m mesh + cvmh sidecar, `cp_max ≈ 8000 m/s` (basement-rock
P-wave); `h_min ≈ 200 m` (post-tet-quality); rough CFL Δt ≈
`0.5 * 200 / 8000 ≈ 1.25e-2 s`. For tfinal = 12 s, that's **960 steps** —
within budget.

But the sidecars (especially basin / sediment) can produce **cp_min ≈
1500 m/s in sediments**, AND `dt_max = "0.1s"` is then **8x the CFL Δt**.
The plan does not clarify whether `dt_max` clamps `dt_auto` from above (which
would yield instability if `dt_max > CFL Δt`) or `cfl_dt` from below.

Worse: ADER is explicit. `Tpv104SubStepIterator` (existing) uses a fixed Δt
across sub-steps. The plan must specify:
- Is Δt recomputed each macro step? (Required for heterogeneous CFL to
  matter at all.)
- What happens if a sediment sub-region produces Δt < `MFEM_MIN_DT`?

The sibling plan's "Low-confidence" risk list (line 574-576) notes
"`Tpv104SubStepIterator` interaction with heterogeneous CFL" but does not
resolve it; this plan inherits the unresolved risk.

**Trigger:** SAFS dyn driver loads cvmh sidecar with deep sediments.

**Actual behavior:** Either silent instability (dt > CFL) or extreme slowdown
(dt = sediment CFL applied globally).

**Expected behavior:** Plan specifies and tests.

**Suggested fix:**
```diff
+ ### Phase 4a: Heterogeneous-CFL handling (resolve sibling plan's
+ §Risk-Assessment Low-confidence item)
+
+ Goal: clarify how `WaveOperator::ComputeMaxDt(cfl)` interacts with
+ `Tpv104SubStepIterator` on heterogeneous media.  The substep iterator
+ today uses a fixed Δt across substeps (TPV104 production).  For SAFS,
+ either:
+
+ Option A (recommended).  Take Δt_global = wave.ComputeMaxDt(cfl) ONCE
+ at startup using the per-element CFL min reduced via MPI_Allreduce(MIN).
+ If `MaterialField` is static throughout the run (it is — see §Out of
+ Scope), this is sufficient.  Add a one-line audit log:
+   "Δt_global = 1.23e-3 s; controlling element: e=12345 at (x, y, z) =
+    (350e3, 3.7e6, -1200), local cp = 1500 m/s, h = 250 m."
+ This lets the user spot if a single sediment element is dominating Δt.
+
+ Option B.  Sub-stepping with local Δt (out of scope for this plan; a
+ future-work item shared with `Tpv104SubStepIterator` design).
+
+ Plan freezes on Option A.  Acceptance:
+ - [ ] On a constant-material mesh, Δt_global matches the homogeneous
+   formula `cfl * h_min / cp` to 1e-12 relative.
+ - [ ] On the SAFS 1000m mesh + cvmh sidecar, Δt_global ≥ 1e-6 s; if
+   the controlling-element log shows `cp < 1000 m/s`, abort with the
+   suggestion "the velocity sidecar contains a sub-sound-speed cell;
+   verify the sidecar build".
+ - [ ] `cfg.time.dt_max` is treated as an upper bound on Δt_global; if
+   `dt_max < Δt_global`, the plan uses `dt_max` and prints a one-line
+   "tightened by user dt_max" note.
```

---

### [R-113] MODERATE `Phase 0 §[friction.slip_weakening]` (sibling §Phase 0 line 130) — `mu_s_default = 0.677` does NOT match user's `geoffrey2010.md` (mu_s=1.1)

**Category:** DEVIATION (user objective)

**Description:**
Sibling plan's example TOML uses `mu_s_default = 0.677`, `mu_d_default =
0.525`, `d_c_default = 0.40`. These are TPV205 production values, not the
user's `geoffrey2010.md` values (`mu_s=1.1, mu_d=0.5, d_o=0.5m`).

Per user decision **D-3**, the validator must accept the user's values
verbatim (`d_o` alias for `d_c`, drop `mu_s < 1`). The plan's EXAMPLE TOML
should also use the user's values to give a concrete starting point. If the
user copy-edits the example and forgets to swap in `geoffrey2010.md` values,
they run TPV205 by accident.

**Trigger:** User copies `EXAMPLE_spatial_friction_slip_weakening_safs.toml`
and runs.

**Actual behavior:** Runs TPV205-style friction, not geoffrey2010.

**Expected behavior:** The user's named reference is the default.

**Suggested fix:**
```diff
  [friction.slip_weakening]
  # Scalar defaults (applied everywhere unless a spatial rule overrides).
- mu_s_default       = 0.677
- mu_d_default       = 0.525
- d_c_default        = 0.40       # m
+ mu_s_default       = 1.1   # geoffrey2010.md (D-3: validator allows mu_s >= 1)
+ mu_d_default       = 0.5   # geoffrey2010.md
+ d_c_default        = 0.5   # m, geoffrey2010.md (alias `d_o` accepted)
  cohesion_default   = 0.0        # Pa
```
And document in `friction/slip-weakening/geoffrey2010.md` itself that the
file is the source of these defaults; add a backlink.

**Test case:**
```python
def test_R113_example_toml_matches_geoffrey2010():
    toml = load("EXAMPLE_spatial_friction_slip_weakening_safs.toml")
    ref  = parse_md("friction/slip-weakening/geoffrey2010.md")
    assert toml["friction.slip_weakening"]["mu_s_default"] == ref["mu_s"]
    assert toml["friction.slip_weakening"]["mu_d_default"] == ref["mu_d"]
    assert toml["friction.slip_weakening"]["d_c_default"] == ref["d_o"]
```

---

### [R-114] MODERATE `Phase 0 §Validation Rule 2` (line 160) — Barrier sentinel `mu_s == 1.0e6` collides with relaxed validator (D-3 allows `mu_s > 1`)

**Category:** EDGE_CASE

**Description:**
Sibling plan §Phase 0, line 160:
> "`mu_s == 1.0e6` (the barrier sentinel from `dynamic/fault_face_flux.hpp:73`
> comment) marks a barrier DOF; the LSW closed-form solve treats this DOF as
> locked."

Per user decision **D-3**, the validator accepts `mu_s = 1.1` (and any
value ≥ 1). The barrier sentinel `1.0e6` is a magic number; the boundary
between "user wants `mu_s = 1.1`" (valid LSW DOF) and "user signals a barrier
DOF" (the `1.0e6` sentinel) is now blurry. A user who types `mu_s = 1e6` by
accident silently turns a region into a barrier.

**Trigger:** User types `mu_s = 1e6` (intending "very strong patch").

**Actual behavior:** That region is locked.

**Expected behavior:** Explicit `kind = "barrier"` spatial rule that doesn't
rely on a magic mu_s value. The numeric sentinel stays as the back-end
representation but is never user-facing.

**Suggested fix:**
```diff
- `mu_s == 1.0e6` (the barrier sentinel from `dynamic/fault_face_flux.hpp:73`
- comment) marks a barrier DOF; the LSW closed-form solve treats this DOF as
- locked.
+ Barrier DOFs are marked via the spatial-rule schema
+   ```toml
+   [[friction.slip_weakening.spatial]]
+   kind = "barrier"     # locks the affected DOFs; the solver treats them
+                        # as immovable regardless of stress.
+   z_min_m = -20000.0
+   z_max_m = -15000.0
+   ```
+ The resolver translates `kind = "barrier"` to `mu_s = 1.0e6` internally;
+ users never type that number.  The validator rejects user-supplied
+ `mu_s > 1.0e5` (one order below the sentinel) with the message
+ "to mark a barrier region, use `kind = "barrier"` instead of a large mu_s".
```

**Test case:**
```python
def test_R114_user_mu_s_1e6_is_rejected():
    with pytest.raises(AssertionError, match="kind = .barrier."):
        LoadSpatialFrictionConfig(toml_with_mu_s_1e6)

def test_R114_barrier_kind_resolves_to_sentinel():
    cfg = LoadSpatialFrictionConfig(toml_with_barrier_kind)
    p   = ResolveSlipWeakening(cfg.slip_weakening, ...)
    # DOFs inside the barrier box must read mu_s == 1.0e6:
    for i in barrier_dof_indices:
        assert p.mu_s[i] == 1.0e6
```

---

### [R-115] MODERATE `Phase 1 §Detailed Req. 2` (line 295) — `apply_nuc_override` semantics duplicate Phase 0's `[friction.slip_weakening] nuc_*_override_pa`

**Category:** QUALITY (two ways to do the same thing)

**Description:**
The plan defines TWO ways to override the nucleation patch:
1. `[friction.slip_weakening] nuc_tau_override_pa = 81.6e6, nuc_sigma_n_override_pa = 120.0e6`
   (top-level scalars, line 137-138)
2. `[[friction.slip_weakening.spatial]] kind = "box", apply_nuc_override = true`
   (per-rule flag, line 155)

Algorithm Phase 1 step 3 (line 296): "Nucleation override (per Phase 0's
`apply_nuc_override`) — the driver applies the override to `sigma_n_eff` and
`tau_pre` BEFORE calling `InitFaultFaceFluxFromSpatial`; the setter is unaware
of the nuc patch."

This is brittle. Two failure modes:
- User sets `apply_nuc_override = true` in a rule but forgets to define the
  `nuc_*_override_pa` scalars → silent NaN.
- User sets `nuc_*_override_pa` but defines two boxes with
  `apply_nuc_override = true` → which one wins?

**Suggested fix:** Collapse to one mechanism. Either (a) drop the top-level
`nuc_*_override_pa` scalars and require the user to set `tau_pre_x/y/z_pa`
and `sigma_n_pa` directly in the spatial rule, OR (b) drop
`apply_nuc_override` and have one global override box at the top level.

Option (a) is more general (multiple nucleation patches with different
overrides) and aligns with the existing `[[spatial]]` rule pattern. Pick (a):

```diff
- nuc_tau_override_pa = 81.6e6
- nuc_sigma_n_override_pa = 120.0e6
- ...
- [[friction.slip_weakening.spatial]]
- kind        = "box"        # nucleation patch
- ...
- apply_nuc_override = true  # uses the [friction.slip_weakening] nuc_*_override_pa
+ # Nucleation patch: override tau_pre and sigma_n_eff directly via the
+ # spatial rule (no separate top-level scalars).
+ [[friction.slip_weakening.spatial]]
+ kind            = "box"
+ x_min_m         = ...
+ ...
+ # Stress overrides (applied AFTER the stress-source projection):
+ tau_pre_dip_pa    = 0.0       # interpreted in the local fault basis (dip)
+ tau_pre_strike_pa = 81.6e6    # interpreted in the local fault basis (strike)
+ sigma_n_pa        = 120.0e6
```

---

### [R-116] MODERATE `Phase 4 §Edge Cases` (line 481) — `T-WAVEOP-BIMATERIAL-SHARED-FACE-MPI` is referenced but its scope is unspecified

**Category:** ASSUMPTION

**Description:**
Plan line 481: "Bi-material face neighbour MPI exchange (`GodunovFluxPool`
Phase 3 Detailed Req. 7a from heterogeneous_material_plan) — required at SAFS
partition boundaries; passes the new heterogeneous_material_plan test
`T-WAVEOP-BIMATERIAL-SHARED-FACE-MPI` which this driver implicitly
exercises."

The cited test is in another plan (heterogeneous_material_plan.md) and is
itself unimplemented. "Implicitly exercises" is not a falsifiable test
criterion.

What this plan must commit to:
- A specific bit-exact check that bi-material partition boundaries do not
  produce a discontinuity in slip-rate at the rank boundary.
- A regression that fails when the MPI exchange is broken.

**Suggested fix:** Add to Phase 4 acceptance:
```diff
+ - [ ] Bi-material shared-face regression: split the SAFS 1000m mesh
+   across 4 ranks; place the basement/basin material interface so it
+   crosses at least one partition boundary; assert that the slip-rate
+   field at t = 0.5s is C^0 across the partition boundary to within
+   1e-8 relative (no rank-edge discontinuity).  This is a STRONGER
+   criterion than the heterogeneous_material_plan's
+   `T-WAVEOP-BIMATERIAL-SHARED-FACE-MPI`; if the prereq test passes but
+   this one fails, the integration is broken regardless of the prereq.
```

---

### [R-117] MODERATE `Phase 5 §Goal` (line 502) — `wave.GetFaultFaceFlux()` mismatch with actual API `wave.SetFaultFlux(...) / GetFaultFlux()`

**Category:** BUG (API misnamed)

**Description:**
Plan Phase 4 step 11: "`wave.GetFaultFaceFlux().InitializeImpedancesPerQP(...)`".
Plan Phase 5 line 504: "`FaultFaceFlux& GetFaultFaceFlux(); // already present?
if not, add`".

Actual API (`dynamic/wave_operator.hpp:117-119`):
```cpp
void SetFaultFlux(FaultFaceFlux *ff) { fault_flux_ = ff; }
FaultFaceFlux *GetFaultFlux() { return fault_flux_; }
```
— named `GetFaultFlux`, returns a **pointer**, not a reference.

**Trigger:** Phase 4 step 11 code as written.

**Actual behavior:** Compile error: no member `GetFaultFaceFlux`.

**Expected behavior:** Use `GetFaultFlux()` and dereference, OR add the
`GetFaultFaceFlux()` alias (low priority — confusing).

**Suggested fix:**
```diff
- 11. **Initialise per-QP impedances** (Phase 2 here):
-     `wave.GetFaultFaceFlux().InitializeImpedancesPerQP(material, pmesh, ...)`
+ 11. **Initialise per-DOF impedances** (Phase 2 here, see R-103): the
+     FaultFaceFlux pointer is already plugged in via wave.SetFaultFlux at
+     ctor time; the per-DOF impedance + LSW + pre-stress write happens
+     through the SAFS-specific free function
+     `InitializeFaultDOFs_Spatial(*wave.GetFaultFlux()->dof_data(),
+                                  num_fault_total, ...)`.
+     There is no `GetFaultFaceFlux` member; the existing accessor is
+     `GetFaultFlux()` returning a pointer.
```
(And drop the Phase 5 §Goal line that proposes adding `GetFaultFaceFlux`.)

---

### [R-118] MODERATE `Phase 0/4 §Files-to-Modify` — `spatial/code/` directory and `tests/fixtures/safs_dyn/` directory do not exist; Makefile / .gitignore changes ill-specified

**Category:** QUALITY

**Description:**
Plan Phase 4 §Files to Modify says:
- "`miniapps/seas/Makefile`: Add `SPATIAL_DYN_DRIVER_SRC = drivers/spatial_dyn_driver.cpp`."
- "`miniapps/seas/.gitignore`: add `seas_spatial_dyn_driver` to the binary list."
- "Files to Create: `drivers/spatial_dyn_driver.cpp`,
  `spatial/code/spatial_dyn_driver_init.hpp/.cpp`"
- "`spatial/code/scripts/verify_spatial_dyn_smoke_safs.py`"

Verification:
- `ls miniapps/seas/spatial/` → does not exist.
- `ls miniapps/seas/tests/fixtures/` → does not exist either (the tests
  directory is `tests/unit/`).
- The plan's Phase 6 "committed reference under `tests/fixtures/safs_dyn/`"
  refers to a non-existent directory.

This is a low-blast-radius issue (`mkdir -p` is trivial) but it leaves
implementers wondering whether they should create the dirs, follow some
existing convention, or wait for sign-off.

**Suggested fix:** Add a §File-system layout section that explicitly creates:
- `miniapps/seas/spatial/code/`
- `miniapps/seas/spatial/code/scripts/`
- `miniapps/seas/tests/fixtures/safs_dyn/`
And add `.gitkeep` files to each so the layout is committed before any source.
Document this in Phase 0 (since both the QD and dynamic plans share it).

---

### [R-119] LOW `Constraints §Numerical 73` (line 73) — Mixed-flux mode default `none` is silently `upwind`; verify dispatch label

**Category:** QUALITY

**Description:**
Plan line 73: "**Mixed-flux mode**: default `none` (upwind everywhere, TPV104
production default)".

`wave_operator.hpp:144-148` documents `MixedFluxMode` enum with `None`,
`Adjacent`, `AllContinuous`. `None` is "byte-identical to pre-Mixed-Flux
behavior" — which is upwind-everywhere via `GodunovFlux`. So the plan label
"upwind everywhere" is right, but a future reader who only reads the plan may
think "none" means "no flux" (i.e. broken). Tighten:

```diff
- Mixed-flux mode: default `none` (upwind everywhere, TPV104 production default);
+ Mixed-flux mode: default `none` (i.e. MixedFluxMode::None: pure upwind
+ Godunov on every face, TPV104 production default);
```

---

### [R-120] LOW `Phase 4 §Detailed Req. 14` (line 468) — "suppress initial t=0 station write on restart" relies on a fix not cited verifiably

**Category:** QUALITY

**Description:**
Plan line 468: "Suppress initial `t=0` station write on restart (paraview-
compaction's existing fix, commit `6bd963f`)."

Citing a 7-character commit hash in a plan is fragile (the hash is short, the
plan may outlive the branch). Either cite the commit by FULL hash + short
title, or refer to the unit test that exercises the behavior.

**Suggested fix:**
```diff
- Suppress initial `t=0` station write on restart (paraview-compaction's
- existing fix, commit `6bd963f`).
+ Suppress initial `t=0` station write on restart (paraview-compaction
+ R-006, exercised by `test_tpv104_restart_station_write_idempotent` in
+ `tests/unit/test_paraview_output.cpp`).
```

---

## Summary

- **Critical issues:** 10 (R-101..R-110)
- **Moderate issues:** 8 (R-111..R-118)
- **Low issues:** 2 (R-119..R-120)
- **Plan compliance against user objective:** **PARTIAL** — the plan correctly
  targets dynamic-rupture / LSW (so R-001 from the QD review is resolved), and
  with user decisions D-1, D-2, D-3 in hand, R-002/R-005/R-009/R-010 from the
  QD review are also resolvable. What remains is the dynamic-rupture-specific
  surface: misnamed `WaveOperator(MaterialField)` ctor, fictional
  `GodunovFluxPool`, fictional `InitializeImpedancesFromMaterial`,
  inadequate boundary-condition handling, throw-away
  `ElasticityDomainOperator` waste, two-stage init race in `FaultFaceFlux`,
  and the V1-checkpoint `driver_tag` schema mutation risk.
- **Plan compliance against actual codebase APIs:** **INCOMPLETE** — three
  CRITICAL findings (R-102, R-103, R-105) name APIs that do not exist; one
  (R-117) misnames an existing API; one (R-106) proposes a schema mutation
  to a critical I/O contract without a back-compat test.
- **Verdict:** **FAIL — must fix before proceeding.**

The `/code-plan` revision must (in priority order):
1. Inline the shared phases from the QD plan and apply R-101's fixes
   (D-3 relaxed validator with `d_o` alias; QD R-003/R-004/R-008 fixes).
2. Add the §Preconditions section (D-2 hard gate) per R-102.
3. Add Phase 3b for the constant-tensor stress source (D-1) per R-108.
4. Resolve R-103 + R-104 + R-117 by collapsing the proposed
   `InitializeFromSpatialStress*` setters into a single SAFS-specific
   free function `InitializeFaultDOFs_Spatial` that mirrors the existing
   TPV205 pattern.
5. Resolve R-110 by adding Option A (free function + new FaultGeometry ctor
   overload) so the dynamic driver does NOT construct a throw-away
   ElasticityDomainOperator with MUMPS_BLR factorization.
6. Specify the SAFS `BoundaryConfig` explicitly + add reflection-time warning
   per R-107.
7. Fix R-105 + R-106 (checkpoint layout + back-compat-safe driver_tag).
8. Add the two-CVM comparison phase (R-109) and the bi-material shared-face
   regression (R-116).
9. Apply the moderate findings (R-111..R-118).

## Unreviewed Areas

- **`heterogeneous_material_plan.md` Phases 1/3** — only spot-checked via
  `heterogeneous_material_review.md` (Phase 1 implemented; Phase 3 forward-
  looking, not implemented). A full audit of that plan is OUT of scope here
  but is gated by D-2's Preconditions block.
- **`tpv102_driver.cpp`, `tpv102_substep_iterator.hpp`** — TPV102 cycles the
  same `FaultFaceFlux` machinery but has a different friction law setup; not
  exhaustively read.
- **`dynamic/pml_layer.hpp`** — exists; not opened; the R-107 fix assumes the
  PML wiring on `WaveOperator::SetPML(...)` is functional (it appears to be,
  based on the existing API surface).
- **`io/petsc_ts_checkpoint.hpp` V2** — referenced but not used by this
  plan (this plan uses V1); fully audited under the QD review.
- **`spatial_unit_test_plan.md`** — referenced in passing; not audited here
  because it sits downstream of these two driver plans.
- **`stress/code/build_stress_csm_safs.py`** — produces the CSM HDF5 sidecar
  but is upstream of the C++ runtime; not audited (the C++ side only opens
  the resulting HDF5).
