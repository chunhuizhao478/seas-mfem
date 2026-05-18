# Code Review (Round 2): spatial_dynamic_rupture_plan rev-3 (2026-05-18, post-fix)

## Review Scope

- **Plan:** `safs/project_7.0_alternative/document/spatial_dynamic_rupture_plan.md`
- **Round-1 review (consumed):** `safs/project_7.0_alternative/debug_document/spatial_dynamic_rupture_review_2026-05-18.md`
- **Fix report (consumed):** `safs/project_7.0_alternative/debug_document/spatial_dynamic_rupture_fix_2026-05-18.md`
- **Files reviewed (fresh pass on the post-fix state):**
  - `miniapps/seas/spatial/code/spatial_friction.{hpp,cpp}`
  - `miniapps/seas/spatial/code/spatial_velocity.{hpp,cpp}`
  - `miniapps/seas/spatial/code/spatial_stress.{hpp,cpp}`
  - `miniapps/seas/fault/fault_geometry_safs_templated.inl`
  - `miniapps/seas/io/field_coefficient.hpp` (the public-visibility change)
  - `miniapps/seas/tests/unit/test_spatial_friction_config.cpp`
  - `miniapps/seas/tests/unit/test_spatial_friction_resolver.cpp`
  - `miniapps/seas/tests/unit/test_spatial_stress_bundle.cpp`
  - `miniapps/seas/tests/unit/test_spatial_velocity_bundle.cpp`
  - `miniapps/seas/tests/unit/test_spatial_constant_stress_source.cpp`
  - `miniapps/seas/safs/project_7.0_alternative/document/spatial_friction_config_schema.md`
- **Domain context consulted:**
  - `miniapps/seas/CLAUDE.md` (fault basis convention, sign conventions)
  - `mfem/fem/intrules.hpp` (`Geometry`/`Geometries.GetCenter` API)
  - `mfem/fem/geom.hpp` (`Geometry::Type`, `GeomCenter[]` storage)
  - `miniapps/seas/fault/fault_geometry.hpp` (`IsBP5()` accessor location)
  - `miniapps/seas/tests/unit/test_compute_safs_params.cpp` (BuildFixture pattern)

## Round-1 finding verification

All round-1 findings have been addressed:

| Finding | Status |
|---------|--------|
| R-001 (Phases H/4/5/6 missing) | Informational; still missing — out of scope for fix pass |
| R-002 (ResolveForcedRupture IP corner→centroid) | **FIXED** — `Geometries.GetCenter(...)` substituted; F-3 regression added |
| R-003 (resolve_rs_impl same bug) | **FIXED** — same substitution; R-9 regression added |
| R-004 (missing-block silence) | **FIXED partially** — 7 of 8 mandatory blocks now guarded; see R-201 below |
| R-005 (fake stress-bundle tests) | **FIXED** — S-2/S-3 now call `ApplyCsmStressSidecar` against a real geom; S-4 added |
| R-006 (dt_initial validator) | **FIXED** — but introduces R-202 below |
| R-007 (toml_bool type guard) | **FIXED** — T-16 test added |
| R-008 (IsBP5 pre-flight) | **FIXED** — added before `NumZeroNormalFallbacks()` |
| R-009 (test rename + R-1b) | **FIXED** — `R_1a` and `R_1b` both wired in main() |
| R-010 (templated empty-N drift) | **FIXED** — `tau_pre_.SetSize(0)` removed |
| R-011 (misleading L_3 comment) | **FIXED** |
| R-012 (`AbortContainmentFailure` public) | **FIXED** — used in spatial_velocity.cpp |

Test sweep ran clean after the fixes:

- `seas_test_spatial_friction_config`        — **24 / 24**
- `seas_test_spatial_friction_resolver`      — **110 / 110**
- `seas_test_spatial_velocity_bundle`        — **6 / 6**
- `seas_test_spatial_stress_bundle`          — **6 / 6**
- `seas_test_spatial_constant_stress_source` — **45 / 45**
- `seas_test_compute_safs_params`            — **13 / 13** (byte-exact regression intact)

Below are the **new bugs / regressions / missed-edges** I found in this round 2 pass.

---

## Findings

### [R-201] [MODERATE] [spatial_friction.cpp:parse_root] — `[pore_pressure]` block silently accepts defaults despite being schema-mandatory

**Category:** DEVIATION / EDGE_CASE

**Description:**
The R-004 fix added presence guards for seven top-level blocks (`material_constant_fallback`, `mesh`, `velocity`, `stress`, `numerics`, `time`, `output`), but **`[pore_pressure]` was not added**.  The schema (`spatial_friction_config_schema.md` §"Top-level layout") lists `[pore_pressure]` as **mandatory**:

```toml
[pore_pressure]                # mandatory
```

The parser still wraps the block in `if (root.contains("pore_pressure")) { ... }` (lines 492–498), so a TOML config missing `[pore_pressure]` silently parses with `P_p_pa = 0`, `P_p_grad_pa_per_m = 0`, `min_sigma_n_pa = 0` — the struct defaults.  Whether the user intended hydrostatic pressure or no pore pressure is invisible; the simulation runs with zero pore pressure.  This is the same shape of bug that round-1 R-004 fixed for the other six blocks.

**Trigger:**
A TOML config that omits the `[pore_pressure]` block entirely.  The Phase-1 unit tests do not cover this (`MinimalLSWHeader` always emits the block), so the parser will accept it silently in production use.

**Actual behavior:**
Parse succeeds; effective normal stress = total normal stress (no pore pressure subtracted).

**Expected behavior:**
`MFEM_ABORT("Top-level [pore_pressure] block missing")`, matching the other six mandatory-block guards.

**Suggested fix:**
Extend the round-1 R-004 fix to include `pore_pressure`:

```diff
   MFEM_VERIFY(root.contains("material_constant_fallback"),
               "Top-level [material_constant_fallback] block missing");
+  MFEM_VERIFY(root.contains("pore_pressure"),
+              "Top-level [pore_pressure] block missing");
   MFEM_VERIFY(root.contains("mesh"),
               "Top-level [mesh] block missing");
```

**Test case:**
```cpp
// tests/unit/test_spatial_friction_config.cpp — add T-17
static void T_17_missing_pore_pressure_block_aborts()
{
   std::cout << "\n[T-17] missing [pore_pressure] block aborts (R-201)\n";
   std::ostringstream oss;
   oss << "[meta]\nschema_version = 1\nlaw = \"slip_weakening\"\n"
       << "[material_constant_fallback]\nlambda=32e9\nmu=32e9\nrho=2670\n"
       // (no [pore_pressure] block)
       << "[mesh]\npath=\"/dev/null\"\norder=1\n"
       << "[velocity]\nmodel=\"cvmh\"\ndataset_root=\"/tmp/x\"\noverride_path=\"\"\n"
       << "[stress]\nkind=\"constant_tensor\"\n"
       << "sigma_xx_pa=0\nsigma_yy_pa=0\nsigma_zz_pa=0\n"
       << "sigma_xy_pa=0\nsigma_yz_pa=0\nsigma_xz_pa=0\n"
       << "[numerics]\nader_order=2\nmixed_flux=\"none\"\ncfl=0.5\n"
       << "[time]\ntfinal=\"12s\"\n"
       << "[output]\noutput_dir=\"out\"\n"
       << "[friction.slip_weakening]\n"
       << "mu_s_default=1.1\nmu_d_default=0.5\nd_c_default=0.5\n"
       << "cohesion_default=0\n";
   TEST_ASSERT(ParseAbortsInChild(oss.str()),
               "missing [pore_pressure] block must abort");
}
```

---

### [R-202] [MODERATE] [spatial_friction.cpp:parse_root, time block] — `dt_initial` validator collides "auto" sentinel with a literal `-1.0` user input

**Category:** BUG (introduced by the R-006 fix)

**Description:**
The R-006 fix added:

```cpp
MFEM_VERIFY(cfg.time.dt_initial > 0.0
            || cfg.time.dt_initial == static_cast<real_t>(-1.0),
            "[time].dt_initial must be > 0 or the string \"auto\" "
            "(-1 sentinel); got " << cfg.time.dt_initial);
```

The validator accepts `dt_initial == -1.0` as the "auto" sentinel.  But the user can also produce `-1.0` by writing it literally:

- `dt_initial = "auto"` → `SpatialTimeParseSeconds("auto")` → `-1.0` (intended)
- `dt_initial = -1.0` (TOML number) → `toml_time_seconds` → `-1.0` (UNINTENDED — user typed a nonsensical negative value)
- `dt_initial = "-1.0s"` (TOML string) → `SpatialTimeParseSeconds("-1.0s")` → `-1.0` (UNINTENDED)

All three flow through the validator identically; the second and third are silently accepted as "auto" instead of being rejected as nonsensical negative time-steps.  This is the kind of typo a user would never expect to silently succeed.

**Trigger:**
TOML `dt_initial = -1.0` or `dt_initial = "-1.0s"`.

**Actual behavior:**
Parser accepts; downstream Phase 4 driver treats it as "compute from CFL".

**Expected behavior:**
Numeric `-1.0` and the string `"-1.0s"` should be rejected with the same error as `dt_initial = 0` or `dt_initial = -2`.  Only the explicit string `"auto"` should produce the sentinel.

**Suggested fix:**
Track whether the sentinel came from the "auto" string at parse time.  Cleanest approach: introduce a `parse_dt_initial(...)` helper that returns a `std::optional<real_t>` (empty = "auto" sentinel, populated = positive value) and validate accordingly:

```diff
+  // Local helper: parse [time].dt_initial; returns -1.0 ONLY when the
+  // value was the literal string "auto".  Any numeric value (including
+  // -1.0) must be > 0 or the validator aborts.
+  real_t parse_dt_initial(const toml::value& t)
+  {
+     if (!t.contains("dt_initial")) { return -1.0; }    // missing => auto
+     const auto& v = t.at("dt_initial");
+     if (v.is_string())
+     {
+        const std::string s = v.as_string();
+        if (s == "auto") { return -1.0; }
+        const real_t parsed = SpatialTimeParseSeconds(s);
+        MFEM_VERIFY(parsed > 0.0,
+                    "[time].dt_initial = \"" << s << "\" must be > 0 or "
+                    "the literal string \"auto\"");
+        return parsed;
+     }
+     real_t val = 0.0;
+     if (v.is_floating())    { val = static_cast<real_t>(v.as_floating()); }
+     else if (v.is_integer()) { val = static_cast<real_t>(v.as_integer()); }
+     else
+     {
+        MFEM_ABORT("[time].dt_initial must be a number or a string");
+     }
+     MFEM_VERIFY(val > 0.0,
+                 "[time].dt_initial = " << val << " must be > 0 or the "
+                 "literal string \"auto\"");
+     return val;
+  }
@@ parse_root, [time] block @@
-     cfg.time.dt_initial = toml_time_seconds(t, "dt_initial", -1.0);
+     cfg.time.dt_initial = parse_dt_initial(t);
   }
@@ post-block validators @@
-  // dt_initial: positive seconds, OR the -1 sentinel produced by the
-  // string "auto".  Anything else (0, negative non-(-1), …) is a typo.
-  MFEM_VERIFY(cfg.time.dt_initial > 0.0
-              || cfg.time.dt_initial == static_cast<real_t>(-1.0),
-              "[time].dt_initial must be > 0 or the string \"auto\" "
-              "(-1 sentinel); got " << cfg.time.dt_initial);
+  // (dt_initial is fully validated inside parse_dt_initial.)
```

**Test case:**
```cpp
static void T_18_dt_initial_literal_negative_aborts()
{
   std::cout << "\n[T-18] dt_initial = -1.0 (literal) aborts (R-202)\n";
   std::string toml = MinimalLSWHeader();
   const std::string before = "dt_initial=\"auto\"";
   const std::string after  = "dt_initial=-1.0";
   const auto pos = toml.find(before);
   MFEM_VERIFY(pos != std::string::npos,
               "MinimalLSWHeader changed shape — re-grep dt_initial");
   toml.replace(pos, before.size(), after);
   toml += MinimalLSWBlock();
   TEST_ASSERT(ParseAbortsInChild(toml),
               "literal dt_initial = -1.0 must abort (not silently \"auto\")");
}
```

---

### [R-203] [MODERATE] [spatial_friction.cpp:parse_root, output block] — `paraview_*_dt` lacks the `> 0` validator the schema mandates

**Category:** EDGE_CASE / ASSUMPTION

**Description:**
The schema (`spatial_friction_config_schema.md` lines 168–170) validates `paraview_volume_dt`, `paraview_bulk_dt`, `paraview_fault_dt` as `> 0`.  The parser reads them via `toml_time_seconds` and then **applies no validator**:

```cpp
cfg.output.paraview_volume_dt = toml_time_seconds(o, "paraview_volume_dt", 0.05);
cfg.output.paraview_bulk_dt   = toml_time_seconds(o, "paraview_bulk_dt",   0.05);
cfg.output.paraview_fault_dt  = toml_time_seconds(o, "paraview_fault_dt",  0.001);
```

A user-typo `paraview_fault_dt = "0s"` or `paraview_fault_dt = -0.001` is silently accepted; the (eventual) driver gets a non-positive cadence and will misbehave (infinite snapshots, divide-by-zero in the cadence scheduler, etc.).  `toml_time_seconds` also accepts the `"auto"` sentinel (-1), which makes no sense for an output cadence but is silently accepted too.

**Trigger:**
TOML `paraview_*_dt = 0`, `paraview_*_dt = -1`, `paraview_*_dt = "0s"`, or `paraview_*_dt = "auto"`.

**Actual behavior:**
Parser accepts; downstream code receives 0 / -1 / -1.

**Expected behavior:**
`> 0` validator per the schema; reject "auto" for cadences.

**Suggested fix:**
```diff
   MFEM_VERIFY(cfg.output.paraview_volume_zfp_tol >= 0.0
               && cfg.output.paraview_bulk_zfp_tol >= 0.0
               && cfg.output.paraview_fault_zfp_tol >= 0.0,
               "[output] paraview_*_zfp_tol must all be >= 0");
+  MFEM_VERIFY(cfg.output.paraview_volume_dt > 0.0,
+              "[output].paraview_volume_dt must be > 0; got "
+              << cfg.output.paraview_volume_dt);
+  MFEM_VERIFY(cfg.output.paraview_bulk_dt > 0.0,
+              "[output].paraview_bulk_dt must be > 0; got "
+              << cfg.output.paraview_bulk_dt);
+  MFEM_VERIFY(cfg.output.paraview_fault_dt > 0.0,
+              "[output].paraview_fault_dt must be > 0; got "
+              << cfg.output.paraview_fault_dt);
   MFEM_VERIFY(cfg.output.max_snapshots >= 1,
```

**Test case:**
```cpp
static void T_19_paraview_fault_dt_zero_aborts()
{
   std::cout << "\n[T-19] paraview_fault_dt = 0 aborts (R-203)\n";
   std::string toml = MinimalLSWHeader();
   const std::string before = "paraview_fault_dt=\"0.001s\"";
   const std::string after  = "paraview_fault_dt=0.0";
   const auto pos = toml.find(before);
   MFEM_VERIFY(pos != std::string::npos,
               "MinimalLSWHeader changed shape — re-grep paraview_fault_dt");
   toml.replace(pos, before.size(), after);
   toml += MinimalLSWBlock();
   TEST_ASSERT(ParseAbortsInChild(toml),
               "paraview_fault_dt = 0 must abort");
}
```

---

### [R-204] [MODERATE] [test_spatial_friction_resolver.cpp:R_1b_rs_eta_auto_per_dof] — R-1b uses `MakeConstant`, which does NOT exercise the R-003 corner-vs-centroid bug

**Category:** EDGE_CASE / QUALITY (test coverage gap)

**Description:**
The R-009 fix added `R_1b_rs_eta_auto_per_dof` to cover the `eta_auto` code path.  But the test uses:

```cpp
const real_t lam = 32.0e9, mu = 32.0e9, rho = 2670.0;
auto mat = MaterialField::MakeConstant(lam, mu, rho);
```

`MaterialField::MakeConstant` puts the material in `Mode::Constant`, where `EvalAt` returns the scalar triple **without ever calling the Coefficient's `Eval(T, ip)`**:

```cpp
inline void EvalAt(int /*elem*/, ElementTransformation& T,
                   const IntegrationPoint& ip,
                   real_t& lambda_out, ..., real_t& rho_out) const
{
   if (mode == Mode::Constant) { lambda_out = lambda_const; ...; return; }
   ...
}
```

So the IP value (corner vs. centroid) is **irrelevant** for `MakeConstant` material.  R-1b would still pass with the original buggy `ip.Init(0)` code.  This means R-1b does not regression-test R-003 — only `R_9_eta_auto_uses_centroid` (which correctly uses `MakeCoefficient` with a z-varying Vs) catches the bug.

R-1b is still useful as a smoke test for the eta_auto code path, but its test name and comment suggest it covers R-003 ("R-009 fix opens room for this") when it does not.  A clearer split would help.

**Trigger:**
Reintroduce the `ip.Init(0)` bug → R-1b still passes; only R-9 catches it.

**Actual behavior:**
R-1b is a partial smoke test only.

**Expected behavior:**
Either (a) rename R-1b to clarify it is a "Mode::Constant eta_auto smoke" and rely on R-9 for the corner-vs-centroid regression, or (b) reshape R-1b to use `MakeCoefficient` so the formula path is exercised end-to-end with an IP-dependent coefficient.

**Suggested fix:**
Take option (a) — the rename is minimal and R-9 already catches the bug:

```diff
-// R-1b  RS defaults round-trip with eta='auto': asserts the per-DOF
-// formula eta_i = 0.5 * sqrt(mu * rho) is applied via the serial-mesh
-// overload.  Uses MakeConstant so the expected value is closed-form.
-static void R_1b_rs_eta_auto_per_dof()
+// R-1b  Smoke test for the RS eta_auto code path: with MakeConstant
+// material the IP value does not affect mu/rho, so this only proves
+// the resolver dispatches into the auto branch and computes
+// 0.5 * sqrt(mu * rho).  The corner-vs-centroid regression (R-003)
+// is covered separately by R_9_eta_auto_uses_centroid which uses
+// MakeCoefficient with a z-varying Vs.
+static void R_1b_rs_eta_auto_smoke_const_material()
 {
-   std::cout << "\n[R-1b] RS eta='auto' computes 0.5 * sqrt(mu*rho) per DOF\n";
+   std::cout << "\n[R-1b] RS eta='auto' smoke (Mode::Constant): formula dispatches\n";
    ...
@@ main @@
-   R_1b_rs_eta_auto_per_dof();
+   R_1b_rs_eta_auto_smoke_const_material();
```

**Test case:**
The existing R-9 already serves as the genuine regression test; no new test is required if option (a) is taken.

---

### [R-205] [LOW] [spatial_friction.cpp:parse_root, stress block] — Silent default `kind = "constant_tensor"` when `[stress].kind` is omitted

**Category:** EDGE_CASE / QUALITY

**Description:**
Inside the `[stress]` block:

```cpp
cfg.stress.kind = parse_stress_kind(toml_str(s, "kind", "constant_tensor"));
```

If the user supplies `[stress]` (now required after R-004 / R-201) but forgets `kind`, the parser silently defaults to `"constant_tensor"`.  The user then has to supply six `sigma_*_pa` keys; if they instead supplied `sidecar_path`, the validator emits "constant_tensor must NOT set sidecar_path" — a misleading error, because the actual issue is the missing `kind`.

The schema does not give `kind` a default value (line 120: "—" in the Default column), implying it should be required.

**Trigger:**
TOML `[stress]` block without `kind` key.

**Actual behavior:**
Defaults to `kind = "constant_tensor"`, then validator may abort with a misleading message about the wrong key.

**Expected behavior:**
Require `kind` explicitly; abort with a clean "missing kind" message.

**Suggested fix:**
```diff
   if (root.contains("stress"))
   {
      const auto& s = root.at("stress");
+     MFEM_VERIFY(s.contains("kind"),
+                 "[stress].kind is required; must be \"constant_tensor\" "
+                 "or \"sidecar_hdf5\"");
      cfg.stress.kind = parse_stress_kind(toml_str(s, "kind", "constant_tensor"));
```

(Defense in depth — the default value passed to `toml_str` then becomes unreachable, but keep it so the call signature stays uniform.)

---

### [R-206] [LOW] [test_spatial_stress_bundle.cpp:BuildFixture] — `BuildFixture` is called three times per process; each builds a fresh `ElasticityDomainOperator<Mesh>` with MUMPS_BLR factorisation

**Category:** QUALITY (test wall-time)

**Description:**
S-2, S-3, and S-4 each call `BuildFixture()` (lines ~96–117), which constructs a `Mesh` + `ElasticityDomainOperator<Mesh>` with `SolverType::MUMPS_BLR` + `FaultGeometry<Mesh>`.  MUMPS_BLR factorisation of the BP5 1000m operator is expensive (10+ s on this machine).  The triplicate construction adds ~30 s to the test wall-time, with no functional benefit — the fixture state is identical each call.

This is purely a test-runtime concern; no correctness implication.

**Suggested fix:**
Cache the fixture as a function-static so it builds once and is shared across the four tests:

```diff
-std::unique_ptr<Fixture> BuildFixture()
+Fixture* GetSharedFixture()
 {
-   auto fix = std::make_unique<Fixture>();
+   static std::unique_ptr<Fixture> fix;
+   if (fix) { return fix.get(); }
+   fix = std::make_unique<Fixture>();
    ...
-   return fix;
+   return fix.get();
 }
```

Update S-2/S-3/S-4 to use `auto* fix = GetSharedFixture();` instead of `auto fix = BuildFixture();`.

Note: when forking, the child inherits the parent's already-built fixture (copy-on-write).  The MUMPS_BLR factorisation state is preserved across fork — that has been the behaviour of every other fork-based test in this project, so no regression.

---

### [R-207] [LOW] [test_spatial_friction_resolver.cpp:F_3_T_uses_centroid_not_corner] — Unused local variable `const int N = 1;` and `Array<int> attr(1)`

**Category:** QUALITY

**Description:**
Both `F_3_T_uses_centroid_not_corner` and `R_9_eta_auto_uses_centroid` declare:

```cpp
const int N = 1;
Vector dofs(3);  Array<int> attr(1), elem(1);
```

`N` is never referenced after declaration; for F-3, `attr` is also never used (the `ResolveForcedRupture` signature takes `dof_to_elem` only, not `dof_to_attr`).  Harmless but misleading — a reader will spend a moment wondering whether `attr` carries semantic content.

**Suggested fix:**
```diff
 static void F_3_T_uses_centroid_not_corner()
 {
    std::cout << "\n[F-3] R-002 regression: T_forced reads Vs at element "
                 "centroid (NOT reference corner)\n";
    ...
-   const int N = 1;
-   Vector dofs(3);  Array<int> attr(1), elem(1);
+   Vector dofs(3);  Array<int> elem(1);   // attr unused by ResolveForcedRupture
    dofs(0) = 0.5; dofs(1) = 0.5; dofs(2) = 0.5;
-   attr[0] = 101; elem[0] = 0;
+   elem[0] = 0;
```

(For R-9, keep `attr` because `ResolveRateState` does take `dof_to_attr` — but drop the unused `N`.)

---

### [R-208] [LOW] [field_coefficient.hpp] — Public visibility of `AbortRangeFailure` is documented as "for same reason" but has no out-of-tree consumer yet

**Category:** QUALITY

**Description:**
The R-012 fix moved both `AbortContainmentFailure` and `AbortRangeFailure` to `public:` with comments:

```cpp
/// Public so out-of-tree consumers (e.g. spatial_velocity.cpp) can
/// reuse the canonical error message without duplicating its
/// formatting (R-012).
static void AbortContainmentFailure(...);

/// Public for the same reason as AbortContainmentFailure (R-012).
static void AbortRangeFailure(...);
```

`AbortContainmentFailure` now has an out-of-tree consumer (`spatial_velocity.cpp`).  `AbortRangeFailure` does not — no caller outside `field_coefficient.cpp` uses it.  Exposing it pre-emptively widens the API surface for no current benefit.  Either revert `AbortRangeFailure` to `private:` or document the expected consumer.

**Suggested fix:**
Revert `AbortRangeFailure` to `private:` until a real consumer appears:

```diff
   /// Throw the "mesh not contained" abort with a formatted bbox table.
   /// Public so out-of-tree consumers (e.g. spatial_velocity.cpp) can
   /// reuse the canonical error message without duplicating its
   /// formatting (R-012).
   static void AbortContainmentFailure(
      const DataField3D& field,
      real_t mxmin, real_t mxmax,
      real_t mymin, real_t mymax,
      real_t mzmin, real_t mzmax);

-  /// Throw the "post-projection out of declared bounds" abort.
-  /// Public for the same reason as AbortContainmentFailure (R-012).
-  static void AbortRangeFailure(
-     const std::string& field_name,
-     real_t observed_lo, real_t observed_hi,
-     real_t declared_lo, real_t declared_hi);
-
 private:
+  /// Throw the "post-projection out of declared bounds" abort.
+  static void AbortRangeFailure(
+     const std::string& field_name,
+     real_t observed_lo, real_t observed_hi,
+     real_t declared_lo, real_t declared_hi);
```

---

## Summary

- Critical issues: **0**
- Moderate issues: **4**  (R-201 pore_pressure missing-block, R-202 dt_initial sentinel collision, R-203 paraview_*_dt missing validator, R-204 R-1b test does not cover R-003)
- Low issues: **4**  (R-205 stress.kind silent default, R-206 BuildFixture triplicated, R-207 unused locals, R-208 AbortRangeFailure unnecessarily public)
- Plan compliance: **PARTIAL** — Phases 0–3b complete with regressions in place; Phases H, 4, 5, 6 still not implemented (round-1 R-001 acknowledged informational; not re-flagged in round 2).
- Verdict: **PASS WITH FIXES** — no remaining CRITICAL bugs.  The MODERATE findings are config-parser hardening: they prevent silent acceptance of malformed configs, but none corrupt scientific results.  Safe to proceed to Phase H implementation; the fixes from R-201..R-204 can land in parallel.

## Unreviewed Areas

- **Phase H source code** (`dynamic/wave_operator.{hpp,cpp,inl}`, `dynamic/fault_face_flux.{hpp,cpp}`, `dynamic/godunov_flux_pool.{hpp,cpp}`) — still not implemented; nothing to review.
- **Phase 4 driver source** (`drivers/spatial_dyn_driver.cpp`) — still not implemented.
- **Phase 5 free functions** (`dynamic/spatial_setup.hpp`, new `FaultGeometry` ctor, checkpoint `DRIVER_TAG_V1` extension) — still not implemented.
- **Phase 6 sbatches + verification scripts** — only `verify_constant_tensor_projection.py` exists; `verify_spatial_dyn_smoke_safs.py` and `compare_spatial_dyn_velocity_models.py` still missing.
- **`verify_constant_tensor_projection.py`** — Python script; not exercised in this review.
- **TPV/BP5 byte-exact regression suite** (`T-PHASEH-SCALAR-PARITY`, `T-FAULTGEO-LEGACY-CTOR-BYTE-IDENTICAL`, etc.) — only `seas_test_compute_safs_params` (13/13) was directly run; the other gates need Phase 5/H implementation before they can fire.
- **Parallel (MPI) behaviour of the resolver paths** — all tests run serial; the `ResolveRateState` / `ResolveForcedRupture` ParMesh overloads were not exercised in this round.
