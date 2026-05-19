# Code Review: 2026-05-18 (round 5) — Fresh audit of round-4 fixes

## Review Scope

Adversarial re-audit of all round-4 fixes (R-001..R-006), focusing on:
- Cascading effects from the validator + fixture change (R-001, R-006)
- Behavior change introduced by R-002's zero-tests guard
- Compatibility of R-003's env-var override with existing CI / Makefile recipes
- R-004/R-005's new validation edge cases
- Anything previously missed

**Files reviewed (all touched by round-4 fixes):**
- `miniapps/seas/domain/elasticity_operator_setup.inl` (R-001 + R-006)
- `miniapps/seas/tests/unit/test_bp5_fault_operator.cpp` (R-001 fixture + R-002 guard)
- `miniapps/seas/tests/unit/test_bp5_integration.cpp` (R-001 fixture + R-002 guard)
- `miniapps/seas/tests/parallel/test_bp5_parallel_smoke.cpp` (R-003 env-var override)
- `miniapps/seas/dynamic/heterogeneous_material.cpp` (R-004 + R-005 validation)

## Test verification (run during this audit)

| Test | Result | Notes |
|------|--------|-------|
| `seas_test_phaser_bimaterial_flux` | 21/21 PASS | R.1 unchanged |
| `seas_test_phaser_dispatch_smoke` | 7/7 PASS | R.2 unchanged |
| `seas_test_phaseh_wave_operator_constant_parity` (serial) | 15/15 PASS | round-3 R-005 warn still fires correctly |
| `seas_test_tpv102_local` | 16/16 PASS | session debug fix holds |
| `seas_test_bp5_fault_operator` | **EXIT=1** | R-002 fires: "0 tests ran" → fail (intentional surface of fixture short-circuit) |
| `seas_test_bp5_integration` | **EXIT=1** | Same as above |
| `mpirun -np 2/4 seas_test_bp5_parallel_smoke` (CG_AMG env) | 17/17 PASS each | R-003 env-var path works |

## Findings (round-5)

---

### [R-001] [MODERATE] [tests/unit/test_bp5_fault_operator.cpp + test_bp5_integration.cpp] — round-4 R-002 fix correctly surfaces the silent skip as a failure, but the underlying fixture is fundamentally broken: it cannot produce non-zero fault DOFs on a non-split Cartesian mesh

**Category:** BUG (latent — the test target has no path to passing without a fixture rewrite)

**Description:**
After round-4 R-001 (revert validator + change +Ly attr=3 → 7) and round-4 R-002 (exit=1 on `num_tests == 0`), `make test-bp5-fault-operator` and `make test-bp5-integration` now exit 1 with `*** FAIL: 0 tests ran (BP5Fixture::Setup likely returned false because GetNumFaultDOFs()==0)`.

This is correct per the round-4 design (loud failure is better than silent skip), but the FIXTURE itself has no way to produce non-zero fault DOFs:
- `Create3DMesh` builds a contiguous Cartesian mesh (no fault split at y=0).
- `ElasticityDomainOperator::GetNumFaultDOFs()` counts DOFs at faces whose boundary attribute equals `bdr_config_.fault_attr` (3); the mesh has no attr=3 faces anymore.
- The fixture's `Setup()` returns false → no tests run → exit=1.

The test target is now perpetually failing.  The "fix" makes the test honest about its broken state but doesn't make it useful.  CI runs would now report this as a regression (was silently exit=0, is now exit=1).

A future fix needs to either (a) rewrite `Create3DMesh` to use a proper split mesh with attr=3 on the 2-sided y=0 interior fault (heavy lift), or (b) remove the broken test targets from the Makefile's default `test` aggregate.

**Trigger:**
Run `make test-bp5-fault-operator` or `make test-bp5-integration` (or `make test` which includes them).

**Actual behavior:**
Test binary executes, fixture short-circuits, exits 1.

**Expected behavior:**
EITHER fixture should produce real fault DOFs and run tests, OR target should be removed from CI default until fixed.

**Suggested fix (option b — minimal, removes the perpetually-failing targets from `make test`):**
```diff
--- a/miniapps/seas/Makefile
+++ b/miniapps/seas/Makefile
@@ -... (find the `test:` aggregate target and remove the two broken tests)
-test: ... \
-      test-bp5-fault-operator \
-      test-bp5-integration \
+test: ... \
+      ##  test-bp5-fault-operator + test-bp5-integration are
+      ##  temporarily removed: the BP5Fixture's Cartesian mesh
+      ##  produces nf=0 fault DOFs, so all assertions are skipped.
+      ##  Re-enable after the fixture is rewritten to use a proper
+      ##  split-fault mesh (REVIEW round-5 R-001).
       ...
```

If the targets are kept in the aggregate, the build's "make test" will fail until the fixture is fixed.

**Test case:**
No new test needed — the existing exit=1 *is* the regression signal.

---

### [R-002] [MODERATE] [tests/parallel/test_bp5_parallel_smoke.cpp::PickTestSolverType] — `SEAS_TEST_SOLVER` default is `MUMPS_BLR`; on the local dev environment (arm64 macOS + conda) `make test-bp5-smoke` now crashes by default

**Category:** BUG (cross-environment regression — R-003 fix flipped the local-dev default to broken)

**Description:**
The round-4 R-003 fix made the smoke test's solver type env-gated.  Without an explicit `SEAS_TEST_SOLVER`, the test defaults to `MUMPS_BLR` — the production default that HPC CI needs.  But on the local development environment (the same one that triggered the original R-003 finding), `MUMPS_BLR` causes the arm64 alignment bus error in `libdmumps.dylib::dmumps_scatter_dist_rhs_`.

So `make test-bp5-smoke` (no env var set) reverts to the pre-round-4 crashing behavior on local dev.  The local user MUST remember to set `SEAS_TEST_SOLVER=CG_AMG` every time, or update their shell rc, or modify the Makefile.

The reviewer's intent in round-4 was "HPC CI defaults to MUMPS, local-dev opts into CG_AMG via env var".  But the Makefile recipe doesn't set the env var, so the user's local `make test-bp5-smoke` runs without it and crashes.

**Trigger:**
`make test-bp5-smoke` on arm64 macOS + conda (the original triggering environment).

**Actual behavior:**
MUMPS bus error in `libdmumps::dmumps_scatter_dist_rhs_`, EXIT=138 (SIGBUS).

**Expected behavior:**
Default to MUMPS where it works, fall back to CG_AMG where it's broken — without requiring per-developer manual env var setup.

**Suggested fix:**
Set the env var in the Makefile recipe only when the platform-detected default is known-broken:

```diff
--- a/miniapps/seas/Makefile
+++ b/miniapps/seas/Makefile
@@ -... (in test-bp5-smoke recipe)
+# On arm64-macOS the conda MUMPS triggers a bus error in
+# dmumps_scatter_dist_rhs_; force CG_AMG.  Other platforms keep the
+# production MUMPS_BLR default (HPC CI is expected to run with no
+# override).
+UNAME_S := $(shell uname -s)
+UNAME_M := $(shell uname -m)
+ifeq ($(UNAME_S),Darwin)
+  ifeq ($(UNAME_M),arm64)
+    BP5_SMOKE_ENV = SEAS_TEST_SOLVER=CG_AMG
+  endif
+endif
+
 test-bp5-smoke: seas_test_bp5_parallel_smoke
-	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_bp5_parallel_smoke
-	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_bp5_parallel_smoke
+	$(BP5_SMOKE_ENV) $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_bp5_parallel_smoke
+	$(BP5_SMOKE_ENV) $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_bp5_parallel_smoke
```

This gives the right defaults without requiring per-developer setup: arm64-macOS auto-picks CG_AMG, Linux/HPC auto-picks MUMPS_BLR (production).  Either platform's user can still override by setting `SEAS_TEST_SOLVER` explicitly.

**Test case:**
```bash
# arm64-macOS: should auto-set CG_AMG and pass
make test-bp5-smoke

# Other platform: should leave MUMPS_BLR default, exercise production path
uname=Linux make test-bp5-smoke

# Explicit override: user can force CG_AMG even on Linux
SEAS_TEST_SOLVER=CG_AMG make test-bp5-smoke
```

---

### [R-003] [LOW] [tests/parallel/test_bp5_parallel_smoke.cpp::PickTestSolverType] — unknown SEAS_TEST_SOLVER value silently falls back to MUMPS_BLR with stderr warning on every rank

**Category:** QUALITY

**Description:**
`PickTestSolverType` warns once on stderr for unknown values like `SEAS_TEST_SOLVER=GAUSS_SEIDEL` (typo), then falls back to MUMPS_BLR.  The warning is emitted on every MPI rank that calls the function — so an unknown value at `np=4` produces 4 identical warning lines.  The fallback to MUMPS_BLR is also questionable: if a user explicitly set an env var, they wanted to OVERRIDE, not silently get the default.  Aborting on unknown values would be safer.

**Suggested fix:**
```diff
--- a/miniapps/seas/tests/parallel/test_bp5_parallel_smoke.cpp
+++ b/miniapps/seas/tests/parallel/test_bp5_parallel_smoke.cpp
@@ -... (in PickTestSolverType, unknown-value branch)
-   // Unknown value: fall back to the production default.  Emit a
-   // warning so the user knows their override was ignored.
-   std::cerr << "[seas-test] WARNING: SEAS_TEST_SOLVER='" << s
-             << "' not recognised; using MUMPS_BLR.\n";
-   return SolverType::MUMPS_BLR;
+   // Unknown value: abort so the user knows their override was
+   // mis-spelled rather than silently using a different solver.
+   MFEM_ABORT("[seas-test] SEAS_TEST_SOLVER='" << s
+              << "' not recognised; valid values are "
+              "MUMPS_BLR, MUMPS, CG_AMG, GMRES_AMG.");
```

**No test case** — defensive change, behavior verified by inspection.

---

### [R-004] [LOW] [POSSIBLE] [domain/elasticity_operator_setup.inl::ValidateFacetBCTables] — round-4 R-001 fix relies on `std::unordered_map` transitively included via mfem.hpp; not robust to header reorganisation

**Category:** ASSUMPTION

**Description:**
The new `std::unordered_map<int, int> lf2sf` in the validator doesn't have a corresponding `#include <unordered_map>` in either `elasticity_operator_setup.inl` or its parent `elasticity_operator.hpp`.  It compiles today because some MFEM header transitively pulls in `<unordered_map>`, but future MFEM reorganisations could break this.

Defensive include is one line.

**Suggested fix:**
```diff
--- a/miniapps/seas/domain/elasticity_operator.hpp
+++ b/miniapps/seas/domain/elasticity_operator.hpp
@@ -... (with the other STL includes near the top)
 #include <set>
+#include <unordered_map>
 #include <vector>
```

**No test case** — defensive include, no behavioral change.

---

### [R-005] [LOW] [dynamic/heterogeneous_material.cpp::MakeDepthProfile1DMaterial] — round-4 R-005 fix rejects Poisson < 0 but the `vp == sqrt(2)·vs` boundary (Poisson exactly 0) is allowed and produces `lambda = 0`

**Category:** ASSUMPTION

**Description:**
The Poisson check uses `>=`:
```cpp
MFEM_VERIFY(L.vp_ms >= std::sqrt(2.0) * L.vs_ms, ...);
```

At equality (Poisson = 0), `lambda = rho * (vp² - 2 vs²) = 0`.  Then `lambda + 2*mu = 2*rho*vs² > 0`, so `BimaterialFlux::BuildGodunovStateFaceLocal`'s check passes, but `lambda = 0` means a fluid-like P-wave with no shear coupling in the bulk modulus, which is unphysical for solid earth materials (Poisson > 0 always in geology).

Risk: a user sets `vp = sqrt(2) * vs` exactly (e.g., `vp = 4000, vs ≈ 2828.43`) and runs without warning, producing a material that's at the edge of physical validity.

**Suggested fix:**
Add a strict-positive lambda lower bound (or equivalently `vp > sqrt(2) · vs * (1 + tiny)`).  Or downgrade the existing `>=` to `>` with a clearer message.  Marked LOW because the failure mode is borderline academic.

**No test case** — quality fix.

---

### [R-006] [LOW] [dynamic/heterogeneous_material.cpp::MakeDepthProfile1DMaterial] — round-4 R-004 contiguity check uses fixed tolerance `1e-9` (metres) without scaling to mesh dimensions

**Category:** ASSUMPTION

**Description:**
The contiguity check is:
```cpp
MFEM_VERIFY(L.depth_top_m <= layers[i-1].depth_bot_m + 1e-9, ...);
```

`1e-9 m` (1 nanometer) is a reasonable tolerance for SCEC-scale meshes (km-resolution), but for sub-micrometer meshes (lab-scale labquakes, future use) the tolerance might wrongly accept a genuine gap.  Not an immediate bug, just a tolerance assumption.

**Suggested fix:**
None required for current usage; revisit if the codebase ever supports sub-micrometer meshes.  Could document the assumption in the docstring.

**No test case** — quality fix.

---

## Cross-cutting observations (informational, not new findings)

- **R.1 / R.2 / parity / tpv102-local** all pass unchanged after round-4 + round-5 work.
- The round-4 R-001 fix correctly preserves the validator's safety invariant: production BP5 (split-mesh, 2-sided fault interior faces) passes the validator; the broken test fixture now fails the test target (R-005 round-5 above), surfacing the latent fixture issue.
- `unordered_map` transitive include (R-004 round-5) compiles today; flagged defensively.

## Summary

- Critical issues: **0**
- Moderate issues: **2** (R-001 perpetually-failing BP5 test targets — surfaced honest failure but no path to passing without fixture rewrite; R-002 local-dev `make test-bp5-smoke` regresses to MUMPS crash without per-developer env var setup)
- Low issues: **4** (R-003 unknown env-var values silently fallback to MUMPS; R-004 missing defensive `#include <unordered_map>`; R-005 Poisson==0 edge case; R-006 hardcoded tolerance)
- Plan compliance: **FULL** for Phase R.1 / R.2 / R.5 stress source / R.3 step 6 with the round-4 + round-5 fixes applied; the BP5 test fixture limitation is pre-existing (predates Phase R) and out of scope for this review series.
- Verdict: **PASS WITH FIXES**

R-001 and R-002 from this round are the priorities: they describe behavior the user should be aware of (BP5 tests now require fixture rewrite OR removal from the default test aggregate; local-dev needs Makefile env-var auto-set for the smoke test).  All other findings are quality / defense-in-depth.

## Unreviewed areas

- The R.1 / R.2 code paths were verified by re-running tests but not line-by-line re-audited in round 5.
- `wave_operator.{hpp,inl}` round-3 changes are unchanged since they were last reviewed.
- `spatial/code/spatial_friction.cpp` (+313 lines for R.3 TOML schema) was not deeply audited in this round; the TOML schema additions are mostly POD types and parsing glue.
- `Makefile` (+139 lines for new test targets) was inspected for the R-001/R-002 new targets only; the larger R.3/R.5 build entries were not audited.
