# Code Review: PLAN_spatial_seas_quasidynamic_driver — impact on `spatial_dyn_driver` (2026-06-02)

## Review Scope
- **Focused question (user):** "double check — I don't want the changes to affect the `spatial_dyn_driver`
  related code/functions."
- Plan: `miniapps/seas/document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.{md,pdf}`
  (post R-001…R-010 + OBS-1/OBS-2 fixes).
- Status: still **plan-stage** (no implementation). Findings are about whether the plan's proposed edits to
  **shared files** would, if implemented literally, change the behavior or break the build of the existing
  dynamic-rupture driver `drivers/spatial_dyn_driver.cpp` (and the shared code it depends on).
- Shared files the plan modifies (not just adds new files): `domain/elasticity_operator.hpp`,
  `fault/fault_geometry.hpp`, `spatial/code/spatial_friction.{hpp,cpp}`, `spatial/code/spatial_stress.hpp`.
- Domain context: `CLAUDE.md`, project memory [C2]/[1] (BP5 + dynamic-driver no-touch rules), plan §Constraints
  ("`spatial_dyn_driver.cpp` … must build and pass `make test` unchanged. New getters/ctors/config keys are
  **additive**; do not alter existing signatures or defaults.").

## What is CONFIRMED SAFE for `spatial_dyn_driver` (verified, no action needed)
- **`ElasticityDomainOperator` new getters** (`GetFaultDOFToElem/Attr/IntegrationPoints`): `spatial_dyn_driver.cpp`
  does **not** `#include` or construct `ElasticityDomainOperator` (grep: 0 matches). Adding const getters cannot
  affect it.
- **No `-Werror`** in `miniapps/seas/Makefile`, and **every** `StressSourceKind` dispatch is an `if/else-if`
  chain (not a `switch`) — `spatial_dyn_driver.cpp:1337-1459`, `spatial_friction.cpp:929-1057`. So adding an
  enumerator cannot cause a `-Wswitch` build break.
- **`FaultGeometry` 7-arg prebuilt-array ctor** (`fault_geometry.hpp:207`) is what the dynamic driver actually
  calls (`spatial_dyn_driver.cpp:1330-1332`). The plan's new 4-arg `(domain, seed, mpi, bool)` ctor does not
  intersect that call.
- **`SolverSpec solver;` + new `TimeSpec` fields**: purely additive struct members; the dynamic driver reads
  only the existing `tfinal/t_initial/dt_initial/dt_max` (e.g. `spatial_dyn_driver.cpp:1848-1872`), which are
  untouched.
- **Parsing the new `[solver]`/`[time]` keys** in the shared parser is additive (populates new fields the
  dynamic driver ignores) — safe *by itself*. The problem is specifically the dynamic-key **warning** (R-001).

## Findings

### [R-001] MODERATE [Phase 1 §"Files to Modify" — dynamic-key warning in the shared parser] — The "dynamic-only key → ignored" warning is routed through `ParseSpatialFrictionConfigString`, which the dynamic driver also calls, so it would fire on every legitimate dynamic run

**Category:** DEVIATION (violates the plan's own "must build and pass `make test` unchanged" contract for `spatial_dyn_driver`)

**Description:**
Plan Phase 1 says (line 229) "Parse `[solver]` and the new `[time]` keys in `ParseSpatialFrictionConfigString`"
and (line 231) "On load, if any dynamic-only key is present (`numerics.ader_order`, `mixed_flux`,
`interior_flux`, `use_pml`, `time_integrator`), emit **one** rank-0 warning 'ignored by spatial_seas
(quasi-dynamic) driver' and proceed." The Phase-1 unit test (line 238) "asserts a dynamic-only key triggers the
warning path" — i.e. the warning is expected **inside the parse path**.

But `ParseSpatialFrictionConfigString` → `parse_root` is the **single shared parser** used by *both* drivers,
and it **already reads exactly those keys into `cfg.numerics`** that the dynamic driver depends on:
```cpp
// spatial/code/spatial_friction.cpp
cfg.numerics.ader_order = toml_int(n, "ader_order", 2);   // :1063
cfg.numerics.mixed_flux = toml_str(n, "mixed_flux", "none"); // :1064
cfg.numerics.use_pml    = toml_bool(n, "use_pml", false);  // :1066
... interior_flux (:1091) ... time_integrator (:1101) ...
```
(`NumericsSpec` fields at `spatial_friction.hpp:130-139`.) The dynamic driver reads `cfg.numerics.*` to drive
ADER order, flux choice, PML, integrator. So if the dynamic-key warning is added to the shared parser, **every
`spatial_dyn_driver` run that sets any of these (the normal case) prints the spurious warning** "ignored by
spatial_seas (quasi-dynamic) driver" — a user-facing behavior change for a driver the plan promises to leave
unchanged. Worse failure mode: if a future implementer reads "ignored" literally and *skips populating*
`cfg.numerics` when those keys are seen, the dynamic driver silently loses its numerics config (CRITICAL).

**Trigger:** Running `spatial_dyn_driver` with any normal dynamic TOML (which sets `numerics.ader_order` etc.),
once Phase 1's warning logic lives in the shared parser.

**Actual behavior (as planned):** Shared parser emits the QD "ignored" warning for the dynamic driver's own
valid keys (and the Phase-1 test would pass while wiring this contamination in).

**Expected behavior:** The QD driver warns about dynamic-only keys; the dynamic driver never does. The shared
parser keeps silently parsing all keys.

**Suggested fix (plan Phase 1):** move the warning out of the shared parser into QD-driver-specific code (or
behind an explicit driver-context flag).
```diff
-  - Parse `[solver]` and the new `[time]` keys in `ParseSpatialFrictionConfigString`.
+  - Parse `[solver]` and the new `[time]` keys in `ParseSpatialFrictionConfigString` (additive — populates
+    new fields; the dynamic driver ignores them).  Do NOT change how the shared parser treats the existing
+    `numerics.*` keys.
-  - On load, if any dynamic-only key is present (`numerics.ader_order`, `mixed_flux`, `interior_flux`, `use_pml`, `time_integrator`), emit **one** rank-0 warning "ignored by spatial_seas (quasi-dynamic) driver" and proceed.
+  - The "dynamic-only keys ignored" warning MUST be emitted by **`spatial_seas_driver.cpp` only** (a QD-driver
+    helper that inspects the parsed `cfg.numerics`/raw TOML *after* `LoadSpatialFrictionConfig` returns), NOT
+    inside the shared `ParseSpatialFrictionConfigString`/`parse_root`.  The shared parser must keep reading
+    `numerics.*` unchanged so `spatial_dyn_driver` is byte-for-byte unaffected (no spurious warning, no dropped
+    keys).  (If a parser-side check is unavoidable, gate it on an added `Driver{QD,DYN}` context arg defaulting
+    to the dynamic/legacy behavior.)
```
And update the Phase-1 acceptance test (line 238) to assert the warning fires from the **QD driver path**, and
add a negative assertion that the **dynamic** parse path emits no such warning.

**Test case:**
```cpp
TEST(SpatialSeasConfig, R001_DynamicDriverParseEmitsNoQDWarning) {
  // The shared parser must NOT warn for the dynamic driver's own valid keys.
  const std::string dyn_toml = "[numerics]\nader_order=3\nuse_pml=true\n...";
  testing::internal::CaptureStderr();
  auto cfg = ParseSpatialFrictionConfigString(dyn_toml);   // shared parser
  std::string err = testing::internal::GetCapturedStderr();
  EXPECT_EQ(cfg.numerics.ader_order, 3);                   // keys still parsed (not dropped)
  EXPECT_TRUE(cfg.numerics.use_pml);
  EXPECT_EQ(err.find("ignored by spatial_seas"), std::string::npos);  // NO QD warning
}
```

---

### [R-002] MODERATE [Phase 3c — `Bp5AnalyticStressSource` enum extension] — Adding `StressSourceKind::Bp5Analytic` to the shared enum silently changes the meaning of `spatial_dyn_driver`'s catch-all `else` branch

**Category:** DEVIATION / BUG (shared-enum coupling — introduced by the R-007 fix in the prior round)

**Description:**
Phase 3c (added last round to fix R-007) says: "`spatial/code/spatial_friction.{hpp,cpp}`: add
`StressSourceKind::Bp5Analytic` and parse `[stress].kind = "bp5_analytic"`." `StressSourceKind` is a **shared**
enum (`spatial_friction.hpp:206-212`, currently 4 values). The dynamic driver's stress dispatch is an
`if/else-if` chain whose **final `else` is a catch-all** that assumes "anything not the first three kinds is
`SidecarHDF5`":
```cpp
// drivers/spatial_dyn_driver.cpp:1337-1459
if      (kind == ConstantTensor)               { ... }
else if (kind == FaultLocalPrestress)          { ... }
else if (kind == DepthProportionalToShearModulus) { ... }
else { spatial::ApplyCsmStressSidecar(cfg.stress, geom); }   // assumes == SidecarHDF5
```
Adding a 5th enumerator (`Bp5Analytic`) makes that `else` match **both** `SidecarHDF5` *and* `Bp5Analytic`. The
shared parser would accept `kind="bp5_analytic"` for *any* driver, so a config that reaches the dynamic driver
with `kind="bp5_analytic"` is **silently misrouted** to `ApplyCsmStressSidecar` (which then fails or
mis-projects) instead of being rejected. This is precisely the "change affecting `spatial_dyn_driver`" the user
wants to avoid: the dynamic driver's source is untouched but its `else`-branch *semantics* change because of a
shared-enum edit. (Fixing it by hardening the dynamic driver's `else` to reject unknown kinds would *also*
touch the dynamic driver — also undesired.)

**Trigger:** Any `StressSourceKind::Bp5Analytic` config loaded by `spatial_dyn_driver` (e.g. a copy-pasted
verification TOML, or a future SAF config that mistakenly selects it).

**Actual behavior (as planned):** Dynamic driver's catch-all `else` treats `Bp5Analytic` as `SidecarHDF5`.

**Expected behavior:** The BP5-analytic prestress is reachable from the **QD driver only**, with **zero** change
to the shared `StressSourceKind` enum or the dynamic driver's dispatch.

**Suggested fix (plan Phase 3c):** keep the new `Bp5AnalyticStressSource` *class* (additive in
`spatial_stress.hpp`; the dynamic driver never references it) but do **not** extend the shared enum.
```diff
-- `spatial/code/spatial_friction.{hpp,cpp}`: add `StressSourceKind::Bp5Analytic` and parse `[stress].kind =
--  "bp5_analytic"` (no extra fields — it pulls the BP5 functions from `bp5_params.hpp`).
+- Do NOT add a value to the shared `StressSourceKind` enum (it would change `spatial_dyn_driver`'s catch-all
+  `else == SidecarHDF5` invariant).  Two enum-free options, in preference order:
+  (ii-preferred) reuse the EXISTING `StressSourceKind::SidecarHDF5`: ship a tiny generator that evaluates
+        `bp5_params::tau0_vec`/σ_n on the owned fault DOFs into a `tau0` sidecar, referenced by the BP5
+        verification TOML.  Zero shared-code change; the dynamic driver is untouched.
+  (i-alt)  add the `Bp5AnalyticStressSource` CLASS to `spatial_stress.hpp` (additive; dyn driver never
+        references it) and have `spatial_seas_driver.cpp` call `geom.ComputeParams(Bp5AnalyticStressSource{
+        bp5_params}, ...)` behind a QD-only `[stress].bp5_analytic = true` boolean — NOT a shared enum value
+        and NOT a branch in the dynamic driver's dispatch.
```
(Then update the Phase-3c acceptance + the Phase-7 cross-reference, which currently name `StressSourceKind::Bp5Analytic`/`kind="bp5_analytic"`, to the chosen enum-free form. Phase-7 already offers the sidecar as option (ii).)

**Test case:**
```cpp
TEST(SpatialStress, R002_StressSourceKindEnumUnchanged) {
  // The shared enum must still have exactly its 4 original values so the dynamic
  // driver's `else == SidecarHDF5` catch-all stays correct.
  EXPECT_EQ(static_cast<int>(spatial::StressSourceKind::DepthProportionalToShearModulus), 3);
  // (No StressSourceKind::Bp5Analytic.)  BP5-analytic prestress is QD-driver-only:
  Bp5AnalyticStressSource src(bp5_params);          // new class is fine (additive)
  // ... QD driver fills tau_pre_ via geom.ComputeParams(src, ...) without any enum/dyn-dispatch change.
}
```

---

### [R-003] LOW [Phase 3 / Appendix — new `FaultGeometry` 4-arg ctor default argument] — Dynamic driver is unaffected, but the new ctor must declare `compute_bp5_params` with NO default to avoid overload ambiguity with the existing 3-arg BP5 ctor

**Category:** ASSUMPTION (regression-contract guard for the shared `FaultGeometry`)

**Description:**
The plan's new ctor `FaultGeometry(DomainOperator<MeshType>&, const BP5Params& seed, MPIContext*, bool
compute_bp5_params)` (Architecture §new-code 4; Appendix) sits next to the existing 3-arg
`FaultGeometry(DomainOperator<MeshType>&, const BP5Params&, MPIContext* = nullptr)` (`fault_geometry.hpp:118`).
**For `spatial_dyn_driver` this is a non-issue** — it constructs via the 7-arg prebuilt-array ctor
(`fault_geometry.hpp:207`; `spatial_dyn_driver.cpp:1330-1332`), which the new ctor cannot shadow. But if the
implementer gives `compute_bp5_params` a **default value**, a 3-arg call `FaultGeometry(domain, params, &mpi)`
becomes **ambiguous** between the 3-arg ctor and the 4-arg-with-default ctor — breaking `seas_driver.cpp`
(BP5)/BP2 build (those *do* use the 3-arg ctor), which the plan's §Constraints also forbids. The Appendix
already shows the bool without a default; make that **normative**.

**Trigger:** Implementer adds `bool compute_bp5_params = true` (or any default) → BP5/BP2 3-arg call sites fail
to compile (ambiguous overload).

**Actual behavior:** Plan is silent on the no-default requirement; relies on the implementer not adding one.

**Expected behavior:** Explicitly require no default argument on `compute_bp5_params`.

**Suggested fix (plan Architecture §new-code 4):**
```diff
-   A domain-op ctor variant that **skips analytic BP5 params**:
-     `FaultGeometry(DomainOperator<MeshType>& domain, const BP5Params& seed, MPIContext* mpi, bool compute_bp5_params)`; when `compute_bp5_params==false` ...
+   A domain-op ctor variant that **skips analytic BP5 params**:
+     `FaultGeometry(DomainOperator<MeshType>& domain, const BP5Params& seed, MPIContext* mpi, bool compute_bp5_params)`
+     — **`compute_bp5_params` takes NO default argument** (a default would make the existing 3-arg BP5 ctor
+     `(domain, params, mpi=nullptr)` an ambiguous overload and break `seas_driver`/BP2; `spatial_dyn_driver`
+     uses the 7-arg prebuilt ctor and is unaffected either way).  When `compute_bp5_params==false` ...
```

**Test case:**
```cpp
TEST(FaultGeometry, R003_ThreeArgCtorStillUnambiguous) {
  // Existing 3-arg BP5 call site must still compile unambiguously after the new ctor lands.
  FaultGeometry<ParMesh> g(domain, bp5_params, &mpi);   // resolves to the 3-arg ctor, not the 4-arg one
  SUCCEED();
}
```

---

### [R-004] LOW [Phase 1 — `TimeSpec` field duplication] — New QD time fields duplicate existing `dt_initial`/`dt_max` semantics in the shared `TimeSpec`

**Category:** QUALITY (no impact on `spatial_dyn_driver`, but a maintenance/confusion hazard in shared code)

**Description:**
`TimeSpec` already has `dt_initial` and `dt_max` (`spatial_friction.hpp:153-159`), read by the dynamic driver
(`spatial_dyn_driver.cpp:1848-1872`). Phase 1 adds `dt_init` and `dt_max_years` for the QD driver. These are
**additive and do not collide** (the dynamic driver keeps reading `dt_initial`/`dt_max` unchanged — confirmed
safe), but two near-identical names in one shared struct (`dt_init` vs `dt_initial`, `dt_max` (seconds) vs
`dt_max_years`) invite mis-wiring (e.g. a QD user setting `dt_max` expecting years). Not a correctness bug for
the dynamic driver; flagged so the shared struct does not accumulate confusing parallel fields.

**Suggested fix (plan Phase 1):** prefer reusing the existing fields, or rename for clarity.
```diff
-  - Extend `TimeSpec` ... with QD knobs: `double rk45_atol = 1e-7; double rk45_rtol = 1e-50; double dt_init = -1; ... double dt_max_years = 0.1; ...`
+  - Extend `TimeSpec` with QD knobs: `double rk45_atol = 1e-7; double rk45_rtol = 1e-50;` and REUSE the
+    existing `dt_initial` (seconds; <0 ⇒ derive 0.01·L_nuc/V_nuc) and `dt_max` (seconds) rather than adding
+    parallel `dt_init`/`dt_max_years` (the QD driver converts its 0.1-yr cap to seconds and stores it in
+    `dt_max`).  If distinct QD fields are truly needed, name them unambiguously (e.g. `qd_dt_max_seconds`) so
+    they cannot be confused with the dynamic driver's `dt_initial`/`dt_max`.
```

**Test case:** (covered by R-001's negative test that the dynamic driver's `dt_initial`/`dt_max` reads are unchanged.)

---

## Summary
- Critical issues: 0
- Moderate issues: 2 (R-001 shared-parser warning contaminates the dynamic driver; R-002 shared-enum extension changes the dynamic driver's catch-all `else` semantics)
- Low issues: 2 (R-003 ctor default-arg guard — dynamic driver unaffected; R-004 `TimeSpec` field duplication — dynamic driver unaffected)
- Plan compliance with its own "no impact on `spatial_dyn_driver`" constraint: **PARTIAL** — most edits are
  verified additive/safe (ElasticityDomainOperator getters, SolverSpec/TimeSpec members, the 7-arg ctor,
  `[solver]`/`[time]` parsing, no `-Werror`), but **R-001** and **R-002** would change the dynamic driver's
  runtime behavior, and **R-002 is a regression I introduced last round via the R-007 fix.**
- **Verdict: PASS WITH FIXES** — apply R-001 and R-002 (both keep the change entirely inside the QD driver and
  off the shared parser/enum); R-003/R-004 are cheap hardening of the shared `FaultGeometry`/`TimeSpec`.
  After these, the dynamic driver is provably untouched (source, build, and runtime behavior).

## Unreviewed Areas
- I did not exhaustively diff **every** `cfg.numerics.*` / `cfg.stress.*` read in `spatial_dyn_driver.cpp`
  against the full set of QD config additions; I verified the keys named in Phase 1's warning list and the
  stress-dispatch sites. A full field-by-field diff is recommended at implementation time (the R-001 negative
  test is the regression guard).
- The **`bp5_analytic` numerical correctness** (does it reproduce `seas_driver.cpp`'s `tau0`?) is a Phase-7
  parity concern, not a dynamic-driver-impact concern; out of scope here (covered by `REVIEW_2026-06-02.md`).
- `LoadSpatialFrictionConfig` (the Phase-0 file loader wrapping `ParseSpatialFrictionConfigString`): assumed to
  be the same shared entry both drivers use; if the QD driver gets its own loader wrapper, that is the natural
  place for the R-001 warning.
