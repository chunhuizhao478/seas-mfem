# Code Review: Phase 8 — TPV205/102/104 benchmarks (configs, tests, driver wiring) — 2026-05-28

## Review Scope
- Plan: `miniapps/seas/document/fullelasticity_dev/PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.md`
  (Phase 8 §1727; status note §1770).
- Files reviewed (Phase 8 deliverables, uncommitted in the working tree):
  - `tpv205/configs/tpv205_spatial.toml`, `tpv102/configs/tpv102_spatial.toml`,
    `tpv104/configs/tpv104_spatial.toml`
  - `tests/unit/test_tpv_config_parse.cpp`, `tests/unit/test_tpv_toml_stress_sign.cpp`,
    `tests/unit/test_planar_tpv_basis.cpp`
  - `drivers/spatial_dyn_driver.cpp` (Phase-8 `[boundary]` wiring + the CFL/selector paths)
  - `Makefile` (new test targets)
- Oracles cross-checked: `config/tpv{205,102,104}_params.hpp`,
  `dynamic/tpv{205,102,104}_setup.hpp`, `drivers/tpv205_driver.cpp`.
- Domain context: `CLAUDE.md` (z<0=depth, frame conventions, byte-exact contract), the
  mesh z-ranges (probed: both meshes z∈[−60000,0]), the Phase-8 status note's documented deviations.

## What was verified correct (so the fix agent does not "fix" it)
- **All physics constants match the byte-exact oracle.** TPV205: μ_s=0.677, μ_d=0.525, d_c=0.40,
  τ_back=70/τ_nuc=81.6/τ_left=78@−7.5km/τ_right=62@+7.5km MPa, half=1.5km, depth 7.5km, rupture
  region |x|<15 ∧ depth<15. TPV102: a_vw=0.008/a_vs=0.016, b=0.012, Dc=0.02, τ_ini=75, V_ini=1e-12,
  nuc Δτ=25 MPa, R=3km. TPV104: a_in=0.01/a_out=0.02, V_w_in=0.1/V_w_out=1.0, b=0.014, L=0.4, f_w=0.2,
  τ_ini=40, V_ini=1e-16, nuc Δτ=45 MPa. Material μ/λ/ρ identical across all three.
- **Sign / depth convention correct.** Both meshes are z∈[−60000,0] (z<0 downward); the native uses
  `down_dip=|z|`; the TOMLs' `center_z=−7500` and box `z_max=−15000` are consistent.
- **VW hard-box geometry matches the native core** (|x|<15km ∧ depth<15km); box predicate
  `z≥z_min ∧ z≤z_max` with ±inf defaults is correct.
- **Stress dispatch (R-001 from the prior review) is intact**: `FaultLocalPrestress` →
  `ComputeParamsFaultLocal` + patch override (last-match-wins). The three new tests build and pass
  (config-parse 60/60, stress-sign 32/32, planar-basis 77/77); Makefile targets registered.

## Findings

### [R-001] MODERATE — [spatial_dyn_driver.cpp (stress→WaveOperator ctor / interior-flux dispatch)] — `interior_flux` selector is never consumed and there is no guard against `matrix`

**Category:** DEVIATION / BUG (latent)

**Description:**
`cfg.numerics.interior_flux` (Scalar/Matrix) is parsed (`spatial_friction.cpp`), guarded at parse
time (matrix+mixed_flux aborts; matrix requires non-Constant material), set in all three TPV TOMLs
(`interior_flux="scalar"`), and asserted by `test_tpv_config_parse` — but `grep interior_flux drivers/spatial_dyn_driver.cpp`
returns **zero** uses. The driver unconditionally constructs the scalar `WaveOperator`. A config
that requests `interior_flux="matrix"` with a non-Constant material (which *passes* the parser
guards) would then **silently run the scalar path** with no error.

**Trigger:** any future config with `interior_flux="matrix"` + `material.kind="DepthProfile1D"`
(or any non-Constant) and no `mixed_flux`.

**Actual behavior:** driver silently builds the scalar `WaveOperator`; the matrix request is ignored.

**Expected behavior:** until Phase 9 ports the matrix path, the driver must **abort** on
`interior_flux==Matrix` (fail loud, not silent-wrong-physics). All three Phase-8 TPV configs are
`scalar`, so this does not affect them today — it is a robustness hole that turns a Phase-9-not-done
state into a silent wrong result.

**Suggested fix:** add an explicit guard right before the `WaveOperator` is constructed:
```diff
+  // Phase 9 (matrix/bimaterial Riemann) is not yet wired into this driver.
+  // Fail loud rather than silently running the scalar path on a matrix request.
+  MFEM_VERIFY(cfg.numerics.interior_flux == spatial::InteriorFlux::Scalar,
+              "spatial_dyn_driver: interior_flux=\"matrix\" is not yet supported "
+              "(Phase 9 port pending); use interior_flux=\"scalar\".");
   // ... existing scalar WaveOperator construction ...
```

**Test case:**
```cpp
// tests/unit/test_tpv_config_parse.cpp (or a driver-guard test)
void test_R001_matrix_flux_is_rejected_until_phase9() {
   // A config that passes the parser's matrix guards (matrix + non-Constant material,
   // no mixed_flux) must still be rejected by the driver until Phase 9 lands.
   SpatialFrictionConfig cfg = MakeMinimalCfg();
   cfg.numerics.interior_flux = InteriorFlux::Matrix;
   cfg.material.kind = MaterialKind::DepthProfile1D;   // passes parse guard
   // The driver's interior-flux decision (extract into a helper SelectWaveOperatorKind)
   // must throw/abort for Matrix.
   ASSERT_ABORTS(SelectWaveOperatorKind(cfg));
}
```

---

### [R-002] MODERATE — [spatial_dyn_driver.cpp:1419-1420] — `cfl_safety` selector is never consumed; the DG factor is applied unconditionally

**Category:** DEVIATION

**Description:**
Plan D2 (§186) specifies `cfl_safety ∈ {raw, dg}` as an **opt-in**: `"dg"` applies the
`1/(3·(2p+1))` factor, `"raw"` does not. The driver hardcodes the factor regardless of the selector:
```cpp
const real_t dt_cfl = wave.ComputeMaxDt(
   cfg.numerics.cfl / (3.0 * (2.0 * cfg.mesh.order + 1.0)));   // always "dg"
```
`grep cfl_safety drivers/spatial_dyn_driver.cpp` returns only comments. So `cfl_safety="raw"`
silently produces the dg-factored dt. (Same root cause as R-001/R-003: the `[numerics]` selectors
are parsed + tested but not wired into the driver. Note: at least two SAFS configs also set
`cfl_safety`, so the no-op is not TPV-only.)

**Trigger:** any config with `cfl_safety="raw"`.

**Actual behavior:** dt is dg-factored regardless; `"raw"` is a silent no-op.

**Expected behavior:** branch the factor on the selector. The `"dg"` path must remain
byte-identical to the current hardcode (it matches `tpv205_driver.cpp:1242`, the gold), so all
Phase-8 TOMLs (`cfl_safety="dg"`) are unaffected — only the experimental `"raw"` escape hatch is
restored.

**Suggested fix:**
```diff
-  const real_t dt_cfl = wave.ComputeMaxDt(
-     cfg.numerics.cfl / (3.0 * (2.0 * cfg.mesh.order + 1.0)));
+  const real_t dg_factor =
+     (cfg.numerics.cfl_safety == spatial::CflSafety::Dg)
+        ? (1.0 / (3.0 * (2.0 * cfg.mesh.order + 1.0)))
+        : 1.0;   // "raw": no DG safety factor (D2)
+  const real_t dt_cfl = wave.ComputeMaxDt(cfg.numerics.cfl * dg_factor);
```

**Test case:**
```cpp
void test_R002_cfl_safety_raw_vs_dg_differ() {
   // Extract the dt-factor decision into a pure helper CflFactor(cfg) for testability.
   SpatialFrictionConfig cfg; cfg.mesh.order = 1; cfg.numerics.cfl = 0.25;
   cfg.numerics.cfl_safety = CflSafety::Dg;
   const real_t dg  = CflFactor(cfg);          // == 1/(3*3) = 1/9
   cfg.numerics.cfl_safety = CflSafety::Raw;
   const real_t raw = CflFactor(cfg);          // == 1.0
   ASSERT_NE(dg, raw);                          // currently EQUAL (bug)
   ASSERT_NEAR(dg, 1.0/9.0, 1e-15);
   ASSERT_NEAR(raw, 1.0, 1e-15);
}
```

---

### [R-003] MODERATE — [spatial_dyn_driver.cpp:2072] — `fault_iterator` selector is never consumed; TPV205 `"one-shot"` is silently ignored

**Category:** DEVIATION

**Description:**
`cfg.numerics.fault_iterator` (OneShot/Substep) is parsed, set in all three TOMLs (TPV205
`"one-shot"` per plan req 1a; TPV102/104 `"substep"`), and asserted by `test_tpv_config_parse` — but
the driver never reads it. Substep count is fixed at `const int O = std::max(1, cfg.numerics.ader_order)`
(=2 for all three), and the friction iterator is always the substep `LinearSlipWeakeningIterator`/RS
iterator via `RunSubSteps_`. So `fault_iterator` is a silent no-op. (Same root cause as R-001/R-002.)

There is also a **plan/gold inconsistency**: req 1a asks TPV205 to be `"one-shot"`, but the Phase-8
status note states the TPV205 gold is `*_mfadj_p1_O2` (= O2 **substep**). The driver doing O2-substep
matches the gold; the TOML's `"one-shot"` label is therefore misleading and unenforced.

**Trigger:** any config toggling `fault_iterator`; specifically TPV205 set to `"one-shot"`.

**Actual behavior:** always O=`ader_order` substep, regardless of the selector.

**Expected behavior:** the selector must either take effect or be rejected — a parsed,
documented, test-asserted knob that does nothing is a trap. Because honoring `"one-shot"` (O=1)
would *break* the O2 gold match, the safe resolution is to (a) make the selector honest and
(b) align the TPV205 TOML with its gold.

**Suggested fix (two parts):**
```diff
   // (a) driver: reject the unsupported value rather than silently ignoring it
+  MFEM_VERIFY(cfg.numerics.fault_iterator == spatial::FaultIteratorKind::Substep,
+              "spatial_dyn_driver: fault_iterator=\"one-shot\" is not implemented "
+              "in the spatial driver (it always sub-steps with O=ader_order); "
+              "use fault_iterator=\"substep\".");
```
```diff
   # (b) tpv205/configs/tpv205_spatial.toml — match the O2-substep gold
-  fault_iterator = "one-shot" # Phase 8 req 1a: TPV205 one-shot LSW
+  fault_iterator = "substep"  # gold is *_p1_O2 (O2 substep); the driver
+                              # always sub-steps. (Plan req 1a's "one-shot"
+                              # is not supported by the spatial driver.)
```
And update `test_tpv_config_parse.cpp:96-97` to expect `Substep`.
(If instead the team wants real one-shot support, that is a driver feature, not a mechanical fix —
flag back to the planner; the gold would also need regenerating.)

**Test case:**
```cpp
void test_R003_fault_iterator_one_shot_rejected() {
   SpatialFrictionConfig cfg = MakeMinimalLswCfg();
   cfg.numerics.fault_iterator = FaultIteratorKind::OneShot;
   ASSERT_ABORTS(ValidateFaultIteratorSupported(cfg));  // extract a helper
   cfg.numerics.fault_iterator = FaultIteratorKind::Substep;
   ValidateFaultIteratorSupported(cfg);                 // no throw
}
```

---

### [R-004] MODERATE — [tpv102/tpv104 configs `[friction.rate_state.spatial]`] — hard `box` VW/V_w transition breaks the byte-exact contract for TPV102/104

**Category:** DEVIATION (documented, but the byte-exact Phase-8 AC is consequently unmet)

**Description:**
The native TPV102/104 use SCEC Eq.(4)/(5): `a = a_vs + (a_vw−a_vs)·B_strike·B_dip` with a **C∞ tanh
boxcar** `B` of 3 km transition width (`config/tpv102_params.hpp:73-94`, `tpv104_params.hpp:119-160`),
and the same for `V_w`. The TOMLs approximate this with **hard `box` rules** (step at exactly 15 km),
because the resolver rejects `boxcar_taper` rules (Phase 6 made the kind config-only). The VW *core*
is correct, but the 15–18 km transition annulus is `a_vs`/`V_w_out` in the TOML where the native
ramps smoothly. This is documented in both TOML headers and the status note ("the smooth boxcar
V_w/a taper should be wired before claiming byte-exact TPV102/104 gold parity").

This is not a hidden bug, but it means **Phase 8 AC "TPV station traces match gold within tolerance"
cannot be met for TPV102/104** as configured — the rupture-edge tapering differs. It is recorded here
so the fix agent does NOT attempt a mechanical fix (it requires implementing `boxcar_taper` resolver
consumption, a deferred feature) and so the deferral is tracked against the AC.

**Trigger:** TPV102/104 gold-trace comparison near the VW/VS border (|x|≈15–18 km, depth≈15–18 km).

**Actual behavior:** `a`/`V_w` step at 15 km; the 3 km smooth transition is absent.

**Expected behavior:** smooth tanh transition (matches `Boxcar`/`Boxcar_TPV104`).

**Suggested fix:** DO NOT auto-fix. This requires wiring `boxcar_taper` consumption into
`SpatialFrictionResolver::ResolveRateState` (currently rejected at `spatial_friction.cpp:1686`) and
re-authoring the TPV102/104 `[friction.rate_state.spatial]` blocks as `boxcar_taper` rules. Track as
a Phase-8 follow-up / Phase-(8.5) feature, not a review fix.

**Test case:**
```cpp
// Demonstrates the gap (not a regression guard for a mechanical fix):
void test_R004_hardbox_vs_native_transition() {
   // Native a at x=16.5 km (mid-transition), depth 7.5 km: strictly between a_vw and a_vs.
   const real_t a_native = ComputeA(/*along_strike=*/16.5e3, /*down_dip=*/7.5e3); // ~0.012
   ASSERT_GT(a_native, 0.008); ASSERT_LT(a_native, 0.016);
   // Hard-box resolver result at the same point: jumps straight to a_vs.
   // (resolve the TPV102 cfg over a 1-DOF fault at (16.5e3, 0, -7.5e3))
   const real_t a_box = ResolveSingleDofA(tpv102_cfg, 16.5e3, 0.0, -7.5e3);  // == 0.016
   ASSERT_NEAR(a_box, 0.016, 1e-12);
   ASSERT_NE(a_box, a_native);   // documents the byte-exactness gap
}
```

---

### [R-005] LOW — [tpv102/tpv104 configs] — `σ_n` double-source: `[friction.rate_state].sigma_n_default` and `[stress].sigma_n_pa` both set independently

**Category:** ASSUMPTION

**Description:**
Both TPV102 and TPV104 set `sigma_n_pa = 120e6` in `[stress]` (used by `ComputeParamsFaultLocal`
to seed `geom.sigma_n_per_dof`) **and** `sigma_n_default = 120e6` in `[friction.rate_state]`. These
are two independent inputs for the same physical quantity. They agree here (and P_p=0), so there is
no current error, but a future edit that changes one and not the other would silently desynchronize
the seeded effective normal stress from the value the RS resolver/initial-ψ seed uses.

**Trigger:** editing one `σ_n` field but not the other.

**Actual behavior:** two sources of truth for σ_n; no cross-check.

**Expected behavior:** single source, or a parse-time consistency assert.

**Suggested fix:** add a consistency guard in the parser when both are present:
```diff
   // After parsing [stress] and [friction.rate_state]:
+  if (cfg.stress.kind == StressSourceKind::FaultLocalPrestress &&
+      cfg.rate_state.has_value())
+  {
+     MFEM_VERIFY(std::abs(cfg.rate_state->sigma_n_default
+                          - (cfg.stress.sigma_n_pa - cfg.stress.pore_pressure.P_p_pa))
+                 <= 1e-6 * cfg.stress.sigma_n_pa,
+        "[friction.rate_state].sigma_n_default must equal "
+        "[stress].sigma_n_pa - P_p for a fault_local_prestress RS config "
+        "(double-source consistency).");
+  }
```

---

## Summary
- Critical issues: 0
- Moderate issues: 4 (R-001 interior_flux unguarded matrix fallthrough; R-002 cfl_safety no-op;
  R-003 fault_iterator no-op + TPV205 one-shot/gold inconsistency; R-004 hard-box byte-exact gap —
  documented deferral)
- Low issues: 1 (R-005 σ_n double-source)
- Plan compliance: PARTIAL. The deterministic Phase-8 deliverables are present and the **physics
  values, signs, depth convention, and VW-core geometry are all byte-correct against the oracles**;
  the three new tests pass and are wired into the Makefile. However: (a) the three `[numerics]`
  method selectors the TOMLs rely on are parsed/tested but **not consumed by the driver**
  (R-001/R-002/R-003) — for the current configs they coincide with the hardcoded behavior so no
  wrong result occurs *today*, but they are dead knobs and one (interior_flux) is a silent-wrong
  fallthrough; (b) the byte-exact TPV102/104 AC is **unmet** because of the documented hard-box
  approximation (R-004); (c) per the status note, the 8-rank smoke, gold/SCEC tolerance match, and
  np=2 parity are deferred to a run session.
- Verdict: PASS WITH FIXES. R-001 (matrix guard) and R-002/R-003 (wire or reject the dead selectors)
  are quick, mechanical, and should land before any non-default `[numerics]` config is run. R-004 is
  a tracked feature deferral (do not auto-fix). R-005 is a hardening nicety.

## Unreviewed Areas
- End-to-end behavior (8-rank `tfinal=0.2s` smoke writing SCEC `.dat`, gold/SCEC trace tolerance,
  np=2 cross-rank parity) — deferred to a run session per the status note; not executable in this
  review (heavy multi-rank compute + tolerance judgment).
- The RS initial-ψ equilibrium seed for TPV104's `V_init=1e-16` / `a=0.01` (numerical-stability of
  the logsinh form) — Phase 1–3 code, exercised only at runtime; not re-derived here.
- `mesh.order=1` + `ader_order=2` interaction with the actual TPV meshes at runtime (construction
  succeeds per the status note; stability/accuracy is a run-session concern).
